#!/usr/bin/env bash
###############################################################################
# run_mito_sv.sh — batteries-included HPC launcher for the mtDNA SV pipeline.
#
# The only required argument is a directory of CRAM (and/or BAM) files. The
# script discovers the samples, fans out one SLURM array task per sample (each
# runs all five callers inside the Apptainer image), then submits a dependent
# consolidation job that builds the cohort summary.
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
CALLERS="${MITO_SV_CALLERS:-all}"
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
CONS_TIME="${MITO_SV_CONS_TIME:-01:00:00}"
CONS_MEM="${MITO_SV_CONS_MEM:-8G}"
MAX_ARRAY="${MITO_SV_MAX_ARRAY:-50}"   # max concurrent array tasks (%N)

usage() {
    cat >&2 <<EOF
Usage: run_mito_sv.sh <cram_dir> [options]

  <cram_dir>            directory searched recursively for input files (required)

  --outdir DIR         output root (default: \$PWD/mito_sv_out)
  --image URI          container image to pull (default: $IMAGE)
  --sif PATH           use an existing .sif instead of pulling
  --exts LIST          input extensions to find (default: cram; e.g. cram,bam)
  --callers LIST       callers to run (default: all)
  --threads N          threads per caller step (default: 4)
  --reference FASTA    explicit reference for CRAM decoding (optional)
  --strict             fail a sample if any caller fails
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

# ---- discover samples -> manifest ----------------------------------------
manifest="$OUTDIR/samples.manifest.tsv"
: > "$manifest"
declare -A seen=()
IFS=',' read -ra EXTA <<< "$EXTS"
for ext in "${EXTA[@]}"; do
    while IFS= read -r -d '' f; do
        base="$(basename "$f")"; name="${base%.*}"; name="${name%.chrM}"
        name="$(printf '%s' "$name" | tr ' /' '__')"
        if [[ -n "${seen[$name]:-}" ]]; then
            log "warn: duplicate sample name '$name' ($f) — skipping dup"; continue
        fi
        seen[$name]=1
        printf '%s\t%s\n' "$name" "$f" >> "$manifest"
    done < <(find "$INPUT_DIR" -type f -iname "*.${ext}" -print0 | sort -z)
done

n="$(wc -l < "$manifest" | tr -d ' ')"
[[ "$n" -gt 0 ]] || die "no *.{$EXTS} files found under $INPUT_DIR"
log "discovered $n sample(s) -> $manifest"
log "output root: $OUTDIR"
log "image: ${SIF:-$IMAGE}  callers: $CALLERS  threads: $THREADS"

# ---- resolve the container image -----------------------------------------
need() { command -v "$1" >/dev/null 2>&1 || die "$1 not found on PATH"; }

resolve_sif() {
    if [[ -n "$SIF" ]]; then
        [[ -f "$SIF" ]] || die "sif not found: $SIF"
        SIF="$(cd "$(dirname "$SIF")" && pwd)/$(basename "$SIF")"
        return
    fi
    need apptainer
    SIF="$OUTDIR/mito_sv.sif"
    if [[ -f "$SIF" ]]; then
        log "reusing existing $SIF"
    else
        log "pulling $IMAGE -> $SIF"
        [[ "$DRYRUN" == 1 ]] || apptainer pull "$SIF" "docker://${IMAGE}"
    fi
}

REF_BIND=()
[[ -n "$REFERENCE" ]] && REF_BIND=(--bind "$(dirname "$(readlink -f "$REFERENCE")")")

# ---- LOCAL mode: run here, no scheduler ----------------------------------
if [[ "$LOCAL" == 1 ]]; then
    resolve_sif
    need apptainer
    while IFS=$'\t' read -r name path; do
        sdir="$OUTDIR/$name"; mkdir -p "$sdir"
        args=(--input "$path" --sample "$name" --outdir "$sdir" --threads "$THREADS" --callers "$CALLERS")
        [[ "$STRICT" == 1 ]] && args+=(--strict)
        [[ -n "$REFERENCE" ]] && args+=(--reference "$REFERENCE")
        cmd=(apptainer exec --bind "$INPUT_DIR" --bind "$OUTDIR" "${REF_BIND[@]}"
             "$SIF" /opt/pipeline/run_sample.sh "${args[@]}")
        log "RUN $name: ${cmd[*]}"
        [[ "$DRYRUN" == 1 ]] || "${cmd[@]}" || log "sample $name returned non-zero"
    done < "$manifest"
    cons=(apptainer exec --bind "$OUTDIR" "$SIF"
          /opt/conda/envs/mitosv/bin/python /opt/pipeline/postprocess.py --root "$OUTDIR")
    log "CONSOLIDATE: ${cons[*]}"
    [[ "$DRYRUN" == 1 ]] || "${cons[@]}"
    log "done (local). Cohort summary in $OUTDIR/cohort_*.tsv"
    exit 0
fi

# ---- SLURM mode -----------------------------------------------------------
need sbatch
resolve_sif

# Shared sbatch flags (per-job --time is set explicitly on each submission).
common_sbatch=(--parsable)
[[ -n "$PARTITION" ]] && common_sbatch+=(--partition "$PARTITION")
[[ -n "$ACCOUNT" ]]   && common_sbatch+=(--account "$ACCOUNT")

export MITO_SV_SIF="$SIF" MITO_SV_OUTDIR="$OUTDIR" MITO_SV_MANIFEST="$manifest" \
       MITO_SV_INPUT_DIR="$INPUT_DIR" MITO_SV_CALLERS="$CALLERS" \
       MITO_SV_THREADS="$THREADS" MITO_SV_REFERENCE="$REFERENCE" \
       MITO_SV_STRICT="$STRICT"

array_cmd=(sbatch "${common_sbatch[@]}"
    --job-name mito-sv --array "1-${n}%${MAX_ARRAY}"
    --cpus-per-task "$THREADS" --mem "$MEM" --time "$TIME"
    --output "$OUTDIR/logs/sample_%A_%a.out"
    "$HERE/sample_job.sbatch")

if [[ "$DRYRUN" == 1 ]]; then
    log "DRY-RUN array: ${array_cmd[*]}"
    log "DRY-RUN consolidate: sbatch --dependency=afterany:<arrayid> $HERE/consolidate.sbatch"
    exit 0
fi

mkdir -p "$OUTDIR/logs"
array_id="$("${array_cmd[@]}")"
log "submitted sample array job: $array_id"

cons_id="$(sbatch "${common_sbatch[@]}" \
    --job-name mito-sv-consolidate \
    --dependency "afterany:${array_id}" \
    --cpus-per-task 2 --mem "$CONS_MEM" --time "$CONS_TIME" \
    --output "$OUTDIR/logs/consolidate_%j.out" \
    "$HERE/consolidate.sbatch")"
log "submitted consolidation job: $cons_id (after array $array_id)"
log "outputs will appear in: $OUTDIR  (cohort_*.tsv after consolidation)"
