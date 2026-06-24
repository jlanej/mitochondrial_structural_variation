#!/usr/bin/env bash
#SBATCH --job-name=mito-sv
###############################################################################
# sample_job.sbatch — one SLURM array task = one sample.
#
# Submitted by run_mito_sv.sh, which exports the MITO_SV_* environment and sets
# --array=1-N. Runs all requested callers inside the Apptainer image and writes
# results to $MITO_SV_OUTDIR/<sample>/.
###############################################################################
set -euo pipefail

: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"; : "${MITO_SV_MANIFEST:?}"
: "${MITO_SV_INPUT_DIR:?}"
THREADS="${MITO_SV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
CALLERS="${MITO_SV_CALLERS:-all}"
TASK="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job}"

# Resolve the container runtime (apptainer or legacy singularity).
if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi

line="$(sed -n "${TASK}p" "$MITO_SV_MANIFEST")"
[[ -n "$line" ]] || { echo "ERROR: no manifest line $TASK" >&2; exit 1; }
name="${line%%$'\t'*}"
path="${line#*$'\t'}"
sdir="$MITO_SV_OUTDIR/$name"
mkdir -p "$sdir"

echo "[$(date)] task=$TASK sample=$name input=$path threads=$THREADS callers=$CALLERS"

binds=(--bind "$MITO_SV_INPUT_DIR" --bind "$MITO_SV_OUTDIR")
args=(--input "$path" --sample "$name" --outdir "$sdir" --threads "$THREADS" --callers "$CALLERS")
[[ "${MITO_SV_STRICT:-0}" == "1" ]] && args+=(--strict)
if [[ -n "${MITO_SV_REFERENCE:-}" ]]; then
    binds+=(--bind "$(dirname "$MITO_SV_REFERENCE")")
    args+=(--reference "$MITO_SV_REFERENCE")
fi

"$RT" exec "${binds[@]}" "$MITO_SV_SIF" /opt/pipeline/run_sample.sh "${args[@]}"
echo "[$(date)] task=$TASK sample=$name finished"
