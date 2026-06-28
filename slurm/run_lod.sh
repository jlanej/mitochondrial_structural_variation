#!/usr/bin/env bash
###############################################################################
# run_lod.sh — batteries-included LOD (limit-of-detection) sweep for all callers.
#
# For each (deletion x heteroplasmy x depth x replicate) it simulates a chrM BAM
# carrying that deletion at that VAF/depth (MitoHPC's deterministic simulator),
# runs the callers under two input arms (pipeline bwa-mem + MitoHPC circular-
# aware), and scores detection.
#
# Work is split PER (caller, depth): one SLURM array per combination, named
# lod-<caller>-d<depth>, so no single caller (e.g. eKLIPse) can monopolise a
# submission. Consolidation cascades fire as results land:
#   * cons-<caller>-d<depth>  after each (caller,depth) array  -> realtime runtime
#   * cons-<caller>           after a caller finishes all depths
#   * lod-final               after everything
# Every consolidation re-renders the cumulative report and prints the scope's
# runtime + a rough remaining-wall projection (so you can cancel a slow caller).
#
#   ./slurm/run_lod.sh                       # full grid (default)
#   ./slurm/run_lod.sh --quick               # small grid
#   ./slurm/run_lod.sh --hets 0,0.05,0.1 --depths 1000,2000 --reps 10
#
# Jobs are uncapped unless the projected total exceeds --max-jobs (500): then
# cells-per-job is raised (iterations distributed) until it fits, keeping one
# caller per job.
###############################################################################
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log() { printf '[lod %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# --- deletion catalogue (1-based retained breakpoints) ---
del_bp() {  # <name> <5|3> -> echo breakpoint, or return 1 if unknown
    case "$1" in
        del4977) [[ "$2" == 5 ]] && echo 8469 || echo 13447;;
        del6000) [[ "$2" == 5 ]] && echo 5999 || echo 10999;;
        *) return 1;;
    esac
}
KNOWN_DELS="del4977 del6000"

# --- defaults (LOD grid: 10 hets x 4 depths x 10 reps x 2 deletions = 800 cells) ---
# Heteroplasmy is dense at the low end (where the LOD lives); depth tops out at
# 2000x to expose how real-world depth warps caller runtime (see the report's
# runtime-vs-depth panel) without the runaway cost of 4000x.
HETS="0,0.01,0.02,0.03,0.05,0.08,0.10,0.20,0.30,0.50"
DEPTHS="250,500,1000,2000"
REPS=10
DELS="del4977,del6000"
ARMS="pipeline,circular"
CALLERS="mitohpc,eklipse,mitosalt,splicebreak2,mitomut,mitoseek"   # one array per caller
THREADS="${MITO_SV_LOD_THREADS:-24}"   # cores per array task (--cpus-per-task)
TPC="${MITO_SV_LOD_TPC:-2}"            # threads per cell -> concurrency = THREADS/TPC = 12
CELLS_PER_JOB="${MITO_SV_LOD_CPJ:-12}" # cells/task (= concurrency: one full wave per job)
MAX_JOBS="${MITO_SV_LOD_MAX_JOBS:-500}" # raise cells-per-job if the projection exceeds this
OUTDIR="${MITO_SV_LOD_OUTDIR:-$PWD/lod_out}"
IMAGE="${MITO_SV_IMAGE:-ghcr.io/jlanej/mitochondrial_structural_variation:latest}"
SIF="${MITO_SV_SIF:-}"
PARTITION="${MITO_SV_PARTITION:-}"; ACCOUNT="${MITO_SV_ACCOUNT:-}"
TIME="${MITO_SV_LOD_TIME:-12:00:00}"; MEM="${MITO_SV_LOD_MEM:-64G}"
LOCAL=0; DRYRUN=0

usage() { sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//' >&2; exit "${1:-2}"; }
while [[ $# -gt 0 ]]; do case "$1" in
    --hets) HETS="$2"; shift 2;;
    --depths) DEPTHS="$2"; shift 2;;
    --reps) REPS="$2"; shift 2;;
    --deletions) DELS="$2"; shift 2;;
    --arms) ARMS="$2"; shift 2;;
    --callers) CALLERS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tpc) TPC="$2"; shift 2;;
    --cells-per-job) CELLS_PER_JOB="$2"; shift 2;;
    --max-jobs) MAX_JOBS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
    --sif) SIF="$2"; shift 2;;
    --partition) PARTITION="$2"; shift 2;;
    --account) ACCOUNT="$2"; shift 2;;
    --time) TIME="$2"; shift 2;;
    --mem) MEM="$2"; shift 2;;
    --quick) HETS="0,0.02,0.03,0.05,0.08,0.10,0.20,0.30"; DEPTHS="500,1000,2000"; REPS=8; shift;;
    --local) LOCAL=1; shift;;
    --dry-run) DRYRUN=1; shift;;
    -h|--help) usage 0;;
    *) die "unknown arg: $1";;
esac; done

mkdir -p "$OUTDIR"; OUTDIR="$(cd "$OUTDIR" && pwd)"
IFS=',' read -ra DELA <<< "$DELS"; IFS=',' read -ra HETA <<< "$HETS"
IFS=',' read -ra DEPA <<< "$DEPTHS"; IFS=',' read -ra CALA <<< "$CALLERS"
IFS=',' read -ra ARMA <<< "$ARMS"; n_arms="${#ARMA[@]}"
# Validate + cache breakpoints once.
declare -a DB5 DB3
for i in "${!DELA[@]}"; do
    DB5[$i]="$(del_bp "${DELA[$i]}" 5)" || die "unknown deletion '${DELA[$i]}' (known: $KNOWN_DELS)"
    DB3[$i]="$(del_bp "${DELA[$i]}" 3)"
done

# One manifest per depth (cells are caller-independent; every caller's array for a
# depth reads the same manifest and runs only its caller on those cells).
mkdir -p "$OUTDIR/manifests"
for dp in "${DEPA[@]}"; do
    m="$OUTDIR/manifests/d${dp}.tsv"; : > "$m"
    for ((r=0; r<REPS; r++)); do
        for i in "${!DELA[@]}"; do
            for h in "${HETA[@]}"; do
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${DELA[$i]}" "${DB5[$i]}" "${DB3[$i]}" "$h" "$dp" "$r" >> "$m"
            done
        done
    done
done

cells_per_combo=$(( ${#DELA[@]} * ${#HETA[@]} * REPS ))
[[ "$cells_per_combo" -gt 0 ]] || die "empty grid"
n_combos=$(( ${#CALA[@]} * ${#DEPA[@]} ))
n_cons=$(( n_combos + ${#CALA[@]} + 1 ))   # per-(caller,depth) + per-caller + final
conc=$(( THREADS / TPC )); [[ "$conc" -ge 1 ]] || conc=1

# Pick cells-per-job: keep the default (max parallelism) unless the projected job
# count exceeds --max-jobs, then distribute more iterations per task until it fits.
cpj="$CELLS_PER_JOB"
while :; do
    tasks_per_combo=$(( (cells_per_combo + cpj - 1) / cpj ))
    total_tasks=$(( n_combos * tasks_per_combo ))
    total_jobs=$(( total_tasks + n_cons ))
    [[ "$total_jobs" -le "$MAX_JOBS" || "$cpj" -ge "$cells_per_combo" ]] && break
    cpj=$(( cpj + 1 ))
done
CELLS_PER_JOB="$cpj"
[[ "$CELLS_PER_JOB" -gt "$cells_per_combo" ]] && CELLS_PER_JOB="$cells_per_combo"
tasks_per_combo=$(( (cells_per_combo + CELLS_PER_JOB - 1) / CELLS_PER_JOB ))
total_tasks=$(( n_combos * tasks_per_combo ))
total_jobs=$(( total_tasks + n_cons ))

total_cells=$(( cells_per_combo * ${#DEPA[@]} ))
log "grid: $total_cells cells (${#DELA[@]} dels x ${#HETA[@]} hets x ${#DEPA[@]} depths x $REPS reps), each run by ${#CALA[@]} callers x $n_arms arms; $cells_per_combo cells/(caller,depth) combo, $n_combos combos"
log "-> $total_tasks array task(s) [$tasks_per_combo/combo, $CELLS_PER_JOB cells each, ${conc}x${TPC} on $THREADS cores] + $n_cons consolidations = $total_jobs jobs (cap $MAX_JOBS)"

need() { command -v "$1" >/dev/null 2>&1 || die "$1 not on PATH"; }
resolve_sif() {
    [[ -n "$SIF" ]] && { [[ -f "$SIF" ]] || die "sif not found: $SIF"; SIF="$(cd "$(dirname "$SIF")"&&pwd)/$(basename "$SIF")"; return; }
    need apptainer; SIF="$OUTDIR/mito_sv.sif"
    [[ -f "$SIF" ]] && log "reusing $SIF" || { log "pulling $IMAGE"; [[ "$DRYRUN" == 1 ]] || apptainer pull "$SIF" "docker://$IMAGE"; }
}

export MITO_SV_LOD_OUTDIR="$OUTDIR" MITO_SV_LOD_CPJ="$CELLS_PER_JOB" \
       MITO_SV_LOD_TPC="$TPC" MITO_SV_LOD_CONC="$conc" MITO_SV_LOD_ARMS="$ARMS"

# emit the planned submissions (used by --dry-run and as a launch manifest)
plan() {
    for caller in "${CALA[@]}"; do
        for dp in "${DEPA[@]}"; do
            printf '  array  lod-%-12s-d%-5s  %d task(s) -> manifest d%s.tsv\n' "$caller" "$dp" "$tasks_per_combo" "$dp"
            printf '  cons   cons-%-12s-d%-5s after that array (scope runtime)\n' "$caller" "$dp"
        done
        printf '  cons   cons-%-12s         after all depths of %s\n' "$caller" "$caller"
    done
    printf '  cons   lod-final                  after everything\n'
}

if [[ "$DRYRUN" == 1 ]]; then
    log "DRY-RUN plan ($total_jobs jobs):"; plan >&2; exit 0
fi

if [[ "$LOCAL" == 1 ]]; then
    resolve_sif; export MITO_SV_SIF="$SIF"; mkdir -p "$OUTDIR/logs"
    for caller in "${CALA[@]}"; do
        for dp in "${DEPA[@]}"; do
            scope="${caller}_d${dp}"
            export MITO_SV_LOD_MANIFEST="$OUTDIR/manifests/d${dp}.tsv" \
                   MITO_SV_LOD_CALLER="$caller" MITO_SV_LOD_SCOPE="$scope"
            for ((task=1; task<=tasks_per_combo; task++)); do
                log "LOCAL $scope task $task/$tasks_per_combo"
                SLURM_ARRAY_TASK_ID="$task" bash "$HERE/lod_array.sbatch"
            done
            export MITO_SV_LOD_SCOPE_CALLER="$caller" MITO_SV_LOD_SCOPE_DEPTH="$dp" \
                   MITO_SV_LOD_SCOPE_EXPECT="$(( cells_per_combo * n_arms ))"
            bash "$HERE/lod_consolidate.sbatch"
        done
    done
    export MITO_SV_LOD_SCOPE_CALLER="" MITO_SV_LOD_SCOPE_DEPTH="" MITO_SV_LOD_SCOPE_EXPECT=0
    bash "$HERE/lod_consolidate.sbatch"
    log "done (local). Report: $OUTDIR/lod_report/index.html"; exit 0
fi

need sbatch; resolve_sif; export MITO_SV_SIF="$SIF"
common=(--parsable --time "$TIME"); [[ -n "$PARTITION" ]] && common+=(--partition "$PARTITION"); [[ -n "$ACCOUNT" ]] && common+=(--account "$ACCOUNT")
conscommon=(--cpus-per-task 2 --mem 8G --time 02:00:00)
mkdir -p "$OUTDIR/logs"

# Dependency model. SLURM's `afterany:a:b:c` releases the job after ALL of a,b,c
# terminate (the "any" is exit-STATE — succeeded or failed — unlike afterok); the
# colon list is ANDed. Every consolidation re-reads ALL shards present at its start
# (not other consolidations' output), so each one depends only on the data-producing
# ARRAYS, not on other consolidations: cons-<caller>-d<depth> on its one array,
# cons-<caller> on all its depth arrays, lod-final on every array. That lets each
# fire as early as its inputs exist, and lod-final (after all arrays => all shards)
# always renders the complete report.
all_arr=""
for caller in "${CALA[@]}"; do
    caller_arr=""
    for dp in "${DEPA[@]}"; do
        scope="${caller}_d${dp}"
        export MITO_SV_LOD_MANIFEST="$OUTDIR/manifests/d${dp}.tsv" \
               MITO_SV_LOD_CALLER="$caller" MITO_SV_LOD_SCOPE="$scope"
        arr="$(sbatch "${common[@]}" --job-name "lod-${caller}-d${dp}" --array "1-${tasks_per_combo}" \
            --cpus-per-task "$THREADS" --mem "$MEM" \
            --output "$OUTDIR/logs/lod_${scope}_%A_%a.out" "$HERE/lod_array.sbatch")"
        log "submitted lod-${caller}-d${dp}: $arr ($tasks_per_combo tasks)"
        export MITO_SV_LOD_SCOPE_CALLER="$caller" MITO_SV_LOD_SCOPE_DEPTH="$dp" \
               MITO_SV_LOD_SCOPE_EXPECT="$(( cells_per_combo * n_arms ))"
        cons="$(sbatch "${common[@]}" "${conscommon[@]}" --job-name "cons-${caller}-d${dp}" \
            --dependency "afterany:${arr}" \
            --output "$OUTDIR/logs/cons_${scope}_%j.out" "$HERE/lod_consolidate.sbatch")"
        log "  -> cons-${caller}-d${dp}: $cons (after $arr)"
        caller_arr="${caller_arr:+$caller_arr:}$arr"; all_arr="${all_arr:+$all_arr:}$arr"
    done
    export MITO_SV_LOD_SCOPE_CALLER="$caller" MITO_SV_LOD_SCOPE_DEPTH="" \
           MITO_SV_LOD_SCOPE_EXPECT="$(( cells_per_combo * ${#DEPA[@]} * n_arms ))"
    cons="$(sbatch "${common[@]}" "${conscommon[@]}" --job-name "cons-${caller}" \
        --dependency "afterany:${caller_arr}" \
        --output "$OUTDIR/logs/cons_${caller}_%j.out" "$HERE/lod_consolidate.sbatch")"
    log "submitted cons-${caller}: $cons (after all depths)"
done

export MITO_SV_LOD_SCOPE_CALLER="" MITO_SV_LOD_SCOPE_DEPTH="" MITO_SV_LOD_SCOPE_EXPECT=0
fin="$(sbatch "${common[@]}" "${conscommon[@]}" --job-name "lod-final" \
    --dependency "afterany:${all_arr}" \
    --output "$OUTDIR/logs/lod_final_%j.out" "$HERE/lod_consolidate.sbatch")"
log "submitted lod-final: $fin. Report -> $OUTDIR/lod_report/index.html"
