#!/usr/bin/env bash
#SBATCH --job-name=mito-extract
###############################################################################
# extract_job.sbatch — one array task = pull chrM out of ONE (whole-genome) CRAM,
# ONCE, idempotently.
#
# Real WGS CRAMs are large; this is the only stage that decodes them. The tiny
# chrM slice is cached at $OUT/chrM/<sample>.chrM.bam and reused by prep + every
# caller — and by re-runs: an already-valid slice is skipped, so the big CRAM is
# never re-touched. Runs as a parallel array up front; everything else depends on
# it (aftercorr), so a sample's downstream work starts as soon as ITS slice lands.
###############################################################################
set -euo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"; : "${MITO_SV_MANIFEST:?}"; : "${MITO_SV_INPUT_DIR:?}"
THREADS="${MITO_SV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
TASK="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job}"
OUT="$MITO_SV_OUTDIR"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi

line="$(sed -n "${TASK}p" "$MITO_SV_MANIFEST")"
[[ -n "$line" ]] || { echo "ERROR: no manifest line $TASK" >&2; exit 1; }
name="${line%%$'\t'*}"; path="${line#*$'\t'}"
mkdir -p "$OUT/chrM"
out="$OUT/chrM/${name}.chrM.bam"

# Isolate this task's micromamba cache (shared $HOME lock) — before any micromamba.
cache="$OUT/.cache/extract_${TASK}"
export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

# Idempotent: a valid cached slice (bam + index, passes quickcheck) -> skip; the
# CRAM is never re-decoded on a re-run.
if [[ -f "$out" && -f "$out.bai" ]] && "$RT" exec --bind "$OUT" "$MITO_SV_SIF" \
        micromamba run -n mitosv samtools quickcheck "$out" 2>/dev/null; then
    echo "[$(date)] extract task=$TASK sample=$name: cached $out (skip)"
    exit 0
fi

work="$OUT/.cache/extract_work_${TASK}"; rm -rf "$work"; mkdir -p "$work"
binds=(--bind "$MITO_SV_INPUT_DIR" --bind "$OUT")
args=(--input "$path" --sample "$name" --outdir "$work" --threads "$THREADS" --extract-only "$out")
if [[ -n "${MITO_SV_REFERENCE:-}" ]]; then
    binds+=(--bind "$(dirname "$MITO_SV_REFERENCE")"); args+=(--reference "$MITO_SV_REFERENCE")
fi

echo "[$(date)] extract task=$TASK sample=$name input=$path -> $out"
"$RT" exec "${binds[@]}" "$MITO_SV_SIF" /opt/pipeline/preprocess.sh "${args[@]}"
rm -rf "$work"
echo "[$(date)] extract task=$TASK sample=$name done"
