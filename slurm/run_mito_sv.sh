#!/usr/bin/env bash
###############################################################################
# run_mito_sv.sh — batteries-included HPC launcher for the mtDNA SV pipeline.
#
# The only required argument is a directory of CRAM (and/or BAM) files. The
# script discovers samples, then:
#   1. mito-prep      — one array task per sample: CRAM/BAM -> chrM BAM + FASTQ
#                       (decoded ONCE, reused by every caller)
#   2. mito-<caller>  — one array PER CALLER over all samples, aftercorr on prep
#                       so a sample starts as soon as ITS prep is done. Splitting
#                       per caller means a slow caller (eKLIPse) never blocks a
#                       fast one, and each caller's results land independently.
#   3. cons-<caller>  — fires after each caller finishes ALL its CRAMs: rebuilds
#                       cohort_*.tsv + the interactive cohort_sv_summary.html, so
#                       results appear as soon as the first caller completes.
#   4. mito-final     — after everything: authoritative summary + cleanup.
#
#   ./slurm/run_mito_sv.sh /path/to/crams
#
# All non-software prerequisites (reference indexes, BLAST/LAST databases) are
# baked into the image or generated/downloaded at runtime by the pipeline.
###############################################################################
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log() { printf '[mito-sv %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ---- defaults (override via flags or env) ---------------------------------
INPUT_DIR=""
OUTDIR="${MITO_SV_OUTDIR:-$PWD/mito_sv_out}"
IMAGE="${MITO_SV_IMAGE:-ghcr.io/jlanej/mitochondrial_structural_variation:latest}"
SIF="${MITO_SV_SIF:-}"                 # prebuilt .sif; else pulled from $IMAGE
EXTS="${MITO_SV_EXTS:-cram}"           # comma list: cram,bam
CALLERS="${MITO_SV_CALLERS:-all}"      # comma list or 'all' -> one array each
REFERENCE="${MITO_SV_REFERENCE:-}"     # explicit CRAM reference (optional)
STRICT=0
LOCAL=0
DRYRUN=0
# SLURM resource knobs — one generous profile for every job (extract, prep,
# callers, consolidation): 24 threads, 64 GB, 10 h walltime.
PARTITION="${MITO_SV_PARTITION:-}"
ACCOUNT="${MITO_SV_ACCOUNT:-}"
THREADS="${MITO_SV_THREADS:-24}"
MEM="${MITO_SV_MEM:-64G}"
TIME="${MITO_SV_TIME:-10:00:00}"
# Consolidation is single-threaded Python over TSVs — give it a lean profile so it
# doesn't hold 24 idle cores (override if a very large cohort needs more).
CONS_CPUS="${MITO_SV_CONS_CPUS:-2}"
CONS_MEM="${MITO_SV_CONS_MEM:-16G}"
CONS_TIME="${MITO_SV_CONS_TIME:-02:00:00}"
# Global ceiling on simultaneously-RUNNING array tasks. The six caller arrays
# fire concurrently, so this budget is split evenly across them (%N per array)
# rather than applied per-array — keeping total concurrency at ~MAX_CONCURRENT
# instead of MAX_CONCURRENT*ncallers, which is what tripped AssocMaxSubmitJob.
MAX_CONCURRENT="${MITO_SV_MAX_CONCURRENT:-${MITO_SV_MAX_ARRAY:-500}}"
# AssocMaxSubmitJobLimit ceiling. Every SUBMITTED array element (pending+running)
# counts toward this, NOT just the %N-throttled running ones — so a big cohort can
# trip it even with a modest concurrency budget. Preflight refuses to launch above
# it (default matches the marcotte/user assoc cap; override if yours differs).
SUBMIT_LIMIT="${MITO_SV_SUBMIT_LIMIT:-5000}"
# Samples are PACKED into a fixed number of array tasks per stage: each task
# strides over the manifest and processes its slice of the cohort sequentially.
# This keeps the submitted-job count constant (~8*CHUNKS+7) no matter how large
# the cohort grows, so big runs never trip AssocMaxSubmitJobLimit. The trade-off
# is walltime: a task now runs ceil(n/CHUNKS) samples back-to-back, so raise
# --time for very large cohorts (or raise --chunks to keep slices small).
CHUNKS="${MITO_SV_CHUNKS:-100}"

ALL_CALLERS="mitohpc eklipse mitosalt splicebreak2 mitomut mitoseek"

usage() {
    cat >&2 <<EOF
Usage: run_mito_sv.sh <cram_dir> [options]

  <cram_dir>            directory searched recursively for input files (required)

  --outdir DIR         output root (default: \$PWD/mito_sv_out)
  --image URI          container image to pull (default: $IMAGE)
  --sif PATH           use an existing .sif instead of pulling
  --exts LIST          input extensions to find (default: cram; e.g. cram,bam)
  --callers LIST       callers to run, one array each (default: all six)
  --threads N          threads per caller step (default: 4)
  --reference FASTA    explicit reference for CRAM decoding (optional)
  --max-concurrent N   global cap on running array tasks, split across the
                       caller arrays (default: 500; alias: --max-array)
  --submit-limit N     AssocMaxSubmitJobLimit; refuse to launch if the total
                       submitted array elements would exceed it (default: 5000)
  --chunks N           array tasks per stage; samples are packed across them so
                       each runs ceil(n/N) samples in sequence (default: 100)
  --strict             fail a sample if its caller fails
  --local              run sequentially on THIS node (no sbatch); for testing
  --dry-run            print what would be submitted/run, then exit
  --partition NAME / --account NAME / --time T / --mem M   SLURM knobs
EOF
    exit "${1:-2}"
}

[[ $# -ge 1 ]] || usage 2
INPUT_DIR="$1"; shift || true
while [[ $# -gt 0 ]]; do case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
    --sif) SIF="$2"; shift 2;;
    --exts) EXTS="$2"; shift 2;;
    --callers) CALLERS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --reference) REFERENCE="$2"; shift 2;;
    --max-concurrent|--max-array) MAX_CONCURRENT="$2"; shift 2;;
    --submit-limit) SUBMIT_LIMIT="$2"; shift 2;;
    --chunks) CHUNKS="$2"; shift 2;;
    --strict) STRICT=1; shift;;
    --local) LOCAL=1; shift;;
    --dry-run) DRYRUN=1; shift;;
    --partition) PARTITION="$2"; shift 2;;
    --account) ACCOUNT="$2"; shift 2;;
    --time) TIME="$2"; shift 2;;
    --mem) MEM="$2"; shift 2;;
    -h|--help) usage 0;;
    *) die "unknown arg: $1";;
esac; done

[[ -d "$INPUT_DIR" ]] || die "input dir not found: $INPUT_DIR"
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
mkdir -p "$OUTDIR"; OUTDIR="$(cd "$OUTDIR" && pwd)"

# ---- caller list ----------------------------------------------------------
[[ "$CALLERS" == "all" ]] && CALLERS="$(printf '%s' "$ALL_CALLERS" | tr ' ' ',')"
IFS=',' read -ra CALA <<< "$CALLERS"

# ---- discover samples -> manifest ----------------------------------------
manifest="$OUTDIR/samples.manifest.tsv"
: > "$manifest"
seen=" "   # space-delimited set of names (bash-3.2 safe; no associative array)
IFS=',' read -ra EXTA <<< "$EXTS"
for ext in "${EXTA[@]}"; do
    while IFS= read -r -d '' f; do
        base="$(basename "$f")"; name="${base%.*}"; name="${name%.chrM}"
        name="$(printf '%s' "$name" | tr ' /' '__')"
        case "$seen" in
            *" $name "*) log "warn: duplicate sample name '$name' ($f) — skipping dup"; continue;;
        esac
        seen="$seen$name "
        printf '%s\t%s\n' "$name" "$f" >> "$manifest"
    done < <(find "$INPUT_DIR" -type f -iname "*.${ext}" -print0 | sort -z)
done

n="$(wc -l < "$manifest" | tr -d ' ')"
[[ "$n" -gt 0 ]] || die "no *.{$EXTS} files found under $INPUT_DIR"

# ---- pack samples into a fixed number of array tasks per stage ------------
# TASKS = array size for every stage (extract/prep/each caller). Each task
# strides over the manifest, so cohort size no longer drives the job count.
TASKS=$(( CHUNKS < n ? CHUNKS : n ))
per_chunk=$(( (n + TASKS - 1) / TASKS ))   # ceil(n/TASKS): samples per task

# ---- per-array %N throttles from the global concurrency budget ------------
# extract/prep run early and (largely) alone -> hand them the whole budget.
# The ${#CALA[@]} caller arrays run concurrently -> divide the budget so their
# combined running tasks stay at ~MAX_CONCURRENT.  Never exceed TASKS (pointless)
# and never drop below 1 (a 0 throttle would stall the array).
clamp() { local v="$1" hi="$2"; (( v < 1 )) && v=1; (( v > hi )) && v="$hi"; printf '%s' "$v"; }
SOLO_THR="$(clamp "$MAX_CONCURRENT" "$TASKS")"
CALLER_THR="$(clamp "$(( MAX_CONCURRENT / ${#CALA[@]} ))" "$TASKS")"

log "discovered $n sample(s); callers: ${CALA[*]}"
log "-> packed into $TASKS task(s)/stage (~${per_chunk} sample(s) each): mito-extract[1-$TASKS%$SOLO_THR] + mito-prep[1-$TASKS%$SOLO_THR] + ${#CALA[@]} caller array(s) [1-$TASKS each, %$CALLER_THR] + $(( ${#CALA[@]} + 1 )) consolidations; concurrency budget ${MAX_CONCURRENT} (~$(( CALLER_THR * ${#CALA[@]} )) across callers); work ${THREADS}c/${MEM}/${TIME}, cons ${CONS_CPUS}c/${CONS_MEM}/${CONS_TIME}; out=$OUTDIR"

# ---- preflight: stay under AssocMaxSubmitJobLimit -------------------------
# Submitted = extract(TASKS) + prep(TASKS) + callers(ncallers*TASKS) + cons(ncallers+1).
# %N throttles RUNNING tasks only; all elements still count as submitted. Because
# TASKS is capped by --chunks, this is bounded regardless of cohort size, but we
# still guard (e.g. a very large --chunks) and fail fast rather than mid-loop.
submitted=$(( 2 * TASKS + ${#CALA[@]} * TASKS + ${#CALA[@]} + 1 ))
log "projected submitted jobs: $submitted / AssocMaxSubmitJobLimit $SUBMIT_LIMIT"
if (( submitted > SUBMIT_LIMIT )); then
    max_chunks=$(( (SUBMIT_LIMIT - ${#CALA[@]} - 1) / (2 + ${#CALA[@]}) ))
    die "would submit $submitted jobs, over AssocMaxSubmitJobLimit=$SUBMIT_LIMIT. With ${#CALA[@]} callers, lower --chunks to <= ${max_chunks} (currently $TASKS), drop callers, or raise --submit-limit if your assoc allows."
fi

# ---- resolve the container image -----------------------------------------
need() { command -v "$1" >/dev/null 2>&1 || die "$1 not found on PATH"; }
resolve_sif() {
    if [[ -n "$SIF" ]]; then
        [[ -f "$SIF" ]] || die "sif not found: $SIF"
        SIF="$(cd "$(dirname "$SIF")" && pwd)/$(basename "$SIF")"; return
    fi
    need apptainer
    SIF="$OUTDIR/mito_sv.sif"
    [[ -f "$SIF" ]] && log "reusing existing $SIF" \
        || { log "pulling $IMAGE -> $SIF"; [[ "$DRYRUN" == 1 ]] || apptainer pull "$SIF" "docker://${IMAGE}"; }
}

export MITO_SV_OUTDIR="$OUTDIR" MITO_SV_MANIFEST="$manifest" \
       MITO_SV_INPUT_DIR="$INPUT_DIR" MITO_SV_THREADS="$THREADS" \
       MITO_SV_REFERENCE="$REFERENCE" MITO_SV_STRICT="$STRICT" \
       MITO_SV_NTASKS="$TASKS"   # array size: each task strides the manifest by this step

# ---- DRY-RUN --------------------------------------------------------------
if [[ "$DRYRUN" == 1 ]]; then
    log "DRY-RUN plan (work jobs ${THREADS}c/${MEM}/${TIME}; consolidations ${CONS_CPUS}c/${CONS_MEM}/${CONS_TIME}):"
    log "  array  mito-extract[1-$TASKS%$SOLO_THR]  (~${per_chunk} CRAM(s)/task; chrM out, once, idempotent)"
    log "  array  mito-prep[1-$TASKS%$SOLO_THR]     aftercorr:mito-extract (realign cached slices)"
    for c in "${CALA[@]}"; do
        log "  array  mito-${c} [1-$TASKS%$CALLER_THR]  aftercorr:mito-prep (~${per_chunk} sample(s)/task)"
        log "  cons   cons-${c}  afterany:mito-${c}  (refresh cohort_sv_summary.html)"
    done
    log "  cons   mito-final  afterany:(all caller arrays + per-caller cons) + cleanup"
    exit 0
fi

# ---- LOCAL mode: run here, sequentially (no scheduler) --------------------
if [[ "$LOCAL" == 1 ]]; then
    resolve_sif; export MITO_SV_SIF="$SIF"; mkdir -p "$OUTDIR/logs"
    # Same packing as SLURM: TASKS tasks, each striding the manifest by TASKS.
    for ((t=1; t<=TASKS; t++)); do
        log "LOCAL extract chunk $t/$TASKS"; SLURM_ARRAY_TASK_ID="$t" bash "$HERE/extract_job.sbatch"
    done
    for ((t=1; t<=TASKS; t++)); do
        log "LOCAL prep chunk $t/$TASKS"; SLURM_ARRAY_TASK_ID="$t" bash "$HERE/prep_job.sbatch"
    done
    for c in "${CALA[@]}"; do
        export MITO_SV_CALLER="$c"
        for ((t=1; t<=TASKS; t++)); do
            log "LOCAL $c chunk $t/$TASKS"
            SLURM_ARRAY_TASK_ID="$t" bash "$HERE/sample_job.sbatch" || log "$c chunk $t returned non-zero"
        done
    done
    export MITO_SV_CONS_SCOPE="final" MITO_SV_CONS_FINAL=1
    bash "$HERE/consolidate.sbatch"
    log "done (local). Summary: $OUTDIR/cohort_sv_summary.html"; exit 0
fi

# ---- SLURM mode -----------------------------------------------------------
need sbatch; resolve_sif; export MITO_SV_SIF="$SIF"
mkdir -p "$OUTDIR/logs"
sched=(--parsable)
[[ -n "$PARTITION" ]] && sched+=(--partition "$PARTITION")
[[ -n "$ACCOUNT" ]]   && sched+=(--account "$ACCOUNT")
common=("${sched[@]}" --cpus-per-task "$THREADS" --mem "$MEM" --time "$TIME")          # extract/prep/callers
conscommon=("${sched[@]}" --cpus-per-task "$CONS_CPUS" --mem "$CONS_MEM" --time "$CONS_TIME")  # consolidations

# initial chrM extraction: decode each (WGS) CRAM ONCE, idempotently
extract_id="$(sbatch "${common[@]}" --job-name mito-extract --array "1-${TASKS}%${SOLO_THR}" \
    --output "$OUTDIR/logs/extract_%A_%a.out" "$HERE/extract_job.sbatch")"
log "submitted mito-extract: $extract_id"

# realign each cached chrM slice (off the small BAM, not the CRAM)
prep_id="$(sbatch "${common[@]}" --job-name mito-prep --array "1-${TASKS}%${SOLO_THR}" \
    --dependency "aftercorr:${extract_id}" \
    --output "$OUTDIR/logs/prep_%A_%a.out" "$HERE/prep_job.sbatch")"
log "submitted mito-prep: $prep_id (aftercorr:$extract_id)"

deps=""   # all caller arrays + per-caller cons, for the final job
for c in "${CALA[@]}"; do
    export MITO_SV_CALLER="$c"
    aid="$(sbatch "${common[@]}" --job-name "mito-$c" --array "1-${TASKS}%${CALLER_THR}" \
        --dependency "aftercorr:${prep_id}" \
        --output "$OUTDIR/logs/${c}_%A_%a.out" "$HERE/sample_job.sbatch")"
    log "submitted mito-$c: $aid (aftercorr:$prep_id)"
    export MITO_SV_CONS_SCOPE="$c" MITO_SV_CONS_FINAL=0
    cid="$(sbatch "${conscommon[@]}" --job-name "cons-$c" \
        --dependency "afterany:${aid}" \
        --output "$OUTDIR/logs/cons_${c}_%j.out" "$HERE/consolidate.sbatch")"
    log "  -> cons-$c: $cid (after mito-$c)"
    deps="${deps:+$deps:}$aid:$cid"
done

export MITO_SV_CONS_SCOPE="final" MITO_SV_CONS_FINAL=1
fin="$(sbatch "${conscommon[@]}" --job-name mito-final \
    --dependency "afterany:${deps}" \
    --output "$OUTDIR/logs/final_%j.out" "$HERE/consolidate.sbatch")"
log "submitted mito-final: $fin"
log "outputs will appear in: $OUTDIR  (cohort_sv_summary.html refreshed as each caller finishes)"
