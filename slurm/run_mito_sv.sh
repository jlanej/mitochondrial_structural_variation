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
THREADS="${MITO_SV_THREADS:-4}"
REFERENCE="${MITO_SV_REFERENCE:-}"     # explicit CRAM reference (optional)
STRICT=0
LOCAL=0
DRYRUN=0
# SLURM resource knobs
PARTITION="${MITO_SV_PARTITION:-}"
ACCOUNT="${MITO_SV_ACCOUNT:-}"
TIME="${MITO_SV_TIME:-12:00:00}"
MEM="${MITO_SV_MEM:-16G}"
PREP_TIME="${MITO_SV_PREP_TIME:-04:00:00}"
CONS_TIME="${MITO_SV_CONS_TIME:-01:00:00}"
CONS_MEM="${MITO_SV_CONS_MEM:-8G}"
MAX_ARRAY="${MITO_SV_MAX_ARRAY:-50}"   # max concurrent tasks PER array (%N)

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
  --max-array N        max concurrent tasks per array (%N; default: 50)
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
    --max-array) MAX_ARRAY="$2"; shift 2;;
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
log "discovered $n sample(s); callers: ${CALA[*]}"
log "-> mito-prep[1-$n] + ${#CALA[@]} caller array(s) [1-$n each, %$MAX_ARRAY] + $(( ${#CALA[@]} + 1 )) consolidations; out=$OUTDIR"

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
       MITO_SV_REFERENCE="$REFERENCE" MITO_SV_STRICT="$STRICT"

# ---- DRY-RUN --------------------------------------------------------------
if [[ "$DRYRUN" == 1 ]]; then
    log "DRY-RUN plan:"
    log "  array  mito-prep[1-$n%$MAX_ARRAY]  (preprocess each CRAM once)"
    for c in "${CALA[@]}"; do
        log "  array  mito-${c} [1-$n%$MAX_ARRAY]  aftercorr:mito-prep"
        log "  cons   cons-${c}  afterany:mito-${c}  (refresh cohort_sv_summary.html)"
    done
    log "  cons   mito-final  afterany:(all caller arrays + per-caller cons) + cleanup"
    exit 0
fi

# ---- LOCAL mode: run here, sequentially (no scheduler) --------------------
if [[ "$LOCAL" == 1 ]]; then
    resolve_sif; export MITO_SV_SIF="$SIF"; mkdir -p "$OUTDIR/logs"
    for ((t=1; t<=n; t++)); do
        log "LOCAL prep $t/$n"; SLURM_ARRAY_TASK_ID="$t" bash "$HERE/prep_job.sbatch"
    done
    for c in "${CALA[@]}"; do
        export MITO_SV_CALLER="$c"
        for ((t=1; t<=n; t++)); do
            log "LOCAL $c $t/$n"
            SLURM_ARRAY_TASK_ID="$t" bash "$HERE/sample_job.sbatch" || log "$c sample $t returned non-zero"
        done
    done
    export MITO_SV_CONS_SCOPE="final" MITO_SV_CONS_FINAL=1
    bash "$HERE/consolidate.sbatch"
    log "done (local). Summary: $OUTDIR/cohort_sv_summary.html"; exit 0
fi

# ---- SLURM mode -----------------------------------------------------------
need sbatch; resolve_sif; export MITO_SV_SIF="$SIF"
mkdir -p "$OUTDIR/logs"
common=(--parsable)
[[ -n "$PARTITION" ]] && common+=(--partition "$PARTITION")
[[ -n "$ACCOUNT" ]]   && common+=(--account "$ACCOUNT")
conscommon=(--cpus-per-task 2 --mem "$CONS_MEM" --time "$CONS_TIME")

prep_id="$(sbatch "${common[@]}" --job-name mito-prep --array "1-${n}%${MAX_ARRAY}" \
    --cpus-per-task "$THREADS" --mem "$MEM" --time "$PREP_TIME" \
    --output "$OUTDIR/logs/prep_%A_%a.out" "$HERE/prep_job.sbatch")"
log "submitted mito-prep: $prep_id"

deps=""   # all caller arrays + per-caller cons, for the final job
for c in "${CALA[@]}"; do
    export MITO_SV_CALLER="$c"
    aid="$(sbatch "${common[@]}" --job-name "mito-$c" --array "1-${n}%${MAX_ARRAY}" \
        --dependency "aftercorr:${prep_id}" \
        --cpus-per-task "$THREADS" --mem "$MEM" --time "$TIME" \
        --output "$OUTDIR/logs/${c}_%A_%a.out" "$HERE/sample_job.sbatch")"
    log "submitted mito-$c: $aid (aftercorr:$prep_id)"
    export MITO_SV_CONS_SCOPE="$c" MITO_SV_CONS_FINAL=0
    cid="$(sbatch "${common[@]}" "${conscommon[@]}" --job-name "cons-$c" \
        --dependency "afterany:${aid}" \
        --output "$OUTDIR/logs/cons_${c}_%j.out" "$HERE/consolidate.sbatch")"
    log "  -> cons-$c: $cid (after mito-$c)"
    deps="${deps:+$deps:}$aid:$cid"
done

export MITO_SV_CONS_SCOPE="final" MITO_SV_CONS_FINAL=1
fin="$(sbatch "${common[@]}" "${conscommon[@]}" --job-name mito-final \
    --dependency "afterany:${deps}" \
    --output "$OUTDIR/logs/final_%j.out" "$HERE/consolidate.sbatch")"
log "submitted mito-final: $fin"
log "outputs will appear in: $OUTDIR  (cohort_sv_summary.html refreshed as each caller finishes)"
