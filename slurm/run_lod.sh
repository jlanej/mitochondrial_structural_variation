#!/usr/bin/env bash
###############################################################################
# run_lod.sh — batteries-included LOD (limit-of-detection) sweep for all callers.
#
# For each (deletion x heteroplasmy x depth x replicate) it simulates a chrM BAM
# carrying that deletion at that VAF/depth (MitoHPC's deterministic simulator),
# runs ALL callers under two input arms (pipeline bwa-mem + MitoHPC circular-
# aware), and scores detection. A dependent job aggregates the sweep into LOD
# curves and an interactive HTML report (docs/lod/index.html).
#
#   ./slurm/run_lod.sh                       # full grid (default)
#   ./slurm/run_lod.sh --quick               # small grid
#   ./slurm/run_lod.sh --hets 0,0.05,0.1 --depths 1000,2000 --reps 10
#
# Grid is fully configurable; jobs are chunked so the scheduler is not flooded,
# and each array task uses --threads cores (default 24).
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
THREADS="${MITO_SV_LOD_THREADS:-24}"   # cores per array task (--cpus-per-task)
TPC="${MITO_SV_LOD_TPC:-2}"            # threads per cell -> concurrency = THREADS/TPC = 12
CELLS_PER_JOB="${MITO_SV_LOD_CPJ:-12}" # cells/task (= concurrency: one full wave per job)
OUTDIR="${MITO_SV_LOD_OUTDIR:-$PWD/lod_out}"
IMAGE="${MITO_SV_IMAGE:-ghcr.io/jlanej/mitochondrial_structural_variation:latest}"
SIF="${MITO_SV_SIF:-}"
PARTITION="${MITO_SV_PARTITION:-}"; ACCOUNT="${MITO_SV_ACCOUNT:-}"
TIME="${MITO_SV_LOD_TIME:-12:00:00}"; MEM="${MITO_SV_LOD_MEM:-64G}"
MAX_ARRAY="${MITO_SV_LOD_MAX_ARRAY:-300}"   # don't throttle the array (cluster handles it)
LOCAL=0; DRYRUN=0

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//' >&2; exit "${1:-2}"; }
while [[ $# -gt 0 ]]; do case "$1" in
    --hets) HETS="$2"; shift 2;;
    --depths) DEPTHS="$2"; shift 2;;
    --reps) REPS="$2"; shift 2;;
    --deletions) DELS="$2"; shift 2;;
    --arms) ARMS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tpc) TPC="$2"; shift 2;;
    --cells-per-job) CELLS_PER_JOB="$2"; shift 2;;
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
manifest="$OUTDIR/cells.manifest.tsv"; : > "$manifest"
IFS=',' read -ra DELA <<< "$DELS"; IFS=',' read -ra HETA <<< "$HETS"; IFS=',' read -ra DEPA <<< "$DEPTHS"
# Validate + cache breakpoints once.
declare -a DB5 DB3
for i in "${!DELA[@]}"; do
    DB5[$i]="$(del_bp "${DELA[$i]}" 5)" || die "unknown deletion '${DELA[$i]}' (known: $KNOWN_DELS)"
    DB3[$i]="$(del_bp "${DELA[$i]}" 3)"
done
# Emit cells with DEPTH innermost so each contiguous cells-per-job chunk spans the
# full depth range — every array task then carries a balanced mix of cheap (250x)
# and expensive (2000x) cells instead of one task hogging all the deep ones.
for ((r=0; r<REPS; r++)); do
    for i in "${!DELA[@]}"; do
        for h in "${HETA[@]}"; do for dp in "${DEPA[@]}"; do
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${DELA[$i]}" "${DB5[$i]}" "${DB3[$i]}" "$h" "$dp" "$r" >> "$manifest"
        done; done
    done
done
n="$(wc -l < "$manifest" | tr -d ' ')"
[[ "$n" -gt 0 ]] || die "empty grid"
njobs=$(( (n + CELLS_PER_JOB - 1) / CELLS_PER_JOB ))
conc=$(( THREADS / TPC )); [[ "$conc" -ge 1 ]] || conc=1
log "grid: $n cells ($DELS x ${#HETA[@]} hets x ${#DEPA[@]} depths x $REPS reps), arms=$ARMS"
log "-> $njobs array task(s), $CELLS_PER_JOB cells each, $THREADS threads/task (${conc}x${TPC}); out=$OUTDIR"

need() { command -v "$1" >/dev/null 2>&1 || die "$1 not on PATH"; }
resolve_sif() {
    [[ -n "$SIF" ]] && { [[ -f "$SIF" ]] || die "sif not found: $SIF"; SIF="$(cd "$(dirname "$SIF")"&&pwd)/$(basename "$SIF")"; return; }
    need apptainer; SIF="$OUTDIR/mito_sv.sif"
    [[ -f "$SIF" ]] && log "reusing $SIF" || { log "pulling $IMAGE"; [[ "$DRYRUN" == 1 ]] || apptainer pull "$SIF" "docker://$IMAGE"; }
}

export MITO_SV_LOD_OUTDIR="$OUTDIR" MITO_SV_LOD_MANIFEST="$manifest" \
       MITO_SV_LOD_CPJ="$CELLS_PER_JOB" MITO_SV_LOD_TPC="$TPC" MITO_SV_LOD_CONC="$conc" \
       MITO_SV_LOD_ARMS="$ARMS"

if [[ "$LOCAL" == 1 ]]; then
    resolve_sif
    export MITO_SV_SIF="$SIF"
    for ((task=1; task<=njobs; task++)); do
        log "LOCAL task $task/$njobs"
        [[ "$DRYRUN" == 1 ]] || SLURM_ARRAY_TASK_ID="$task" bash "$HERE/lod_array.sbatch"
    done
    [[ "$DRYRUN" == 1 ]] || bash "$HERE/lod_consolidate.sbatch"
    log "done (local). Report: $OUTDIR/lod_report/index.html"; exit 0
fi

need sbatch; resolve_sif; export MITO_SV_SIF="$SIF"
common=(--parsable --time "$TIME"); [[ -n "$PARTITION" ]] && common+=(--partition "$PARTITION"); [[ -n "$ACCOUNT" ]] && common+=(--account "$ACCOUNT")
mkdir -p "$OUTDIR/logs"
if [[ "$DRYRUN" == 1 ]]; then log "DRY-RUN: sbatch array 1-$njobs%$MAX_ARRAY ($THREADS cpus, $MEM) + consolidate"; exit 0; fi
arr="$(sbatch "${common[@]}" --job-name mito-lod --array "1-${njobs}%${MAX_ARRAY}" \
    --cpus-per-task "$THREADS" --mem "$MEM" --output "$OUTDIR/logs/lod_%A_%a.out" "$HERE/lod_array.sbatch")"
log "submitted LOD array: $arr"
cons="$(sbatch "${common[@]}" --job-name mito-lod-report --dependency "afterany:${arr}" \
    --cpus-per-task 2 --mem 8G --time 02:00:00 --output "$OUTDIR/logs/lod_report_%j.out" "$HERE/lod_consolidate.sbatch")"
log "submitted LOD report job: $cons (after $arr). Report -> $OUTDIR/lod_report/index.html"
