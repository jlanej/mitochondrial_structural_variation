#!/usr/bin/env bash
#SBATCH --job-name=mito-prep
###############################################################################
# prep_job.sbatch — one array task = preprocess ONE sample, ONCE.
#
# CRAM/BAM -> normalised chrM BAM + mito FASTQ in $OUT/prepared/<sample>/. The
# per-caller arrays (sample_job.sbatch) then consume these via run_sample
# --prepared-bam, so a real CRAM is decoded once, not once per caller.
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
pdir="$OUT/prepared/$name"; mkdir -p "$pdir"

# Isolate this task's micromamba cache: Apptainer binds the shared host $HOME, so
# without isolation every concurrent micromamba serialises on one mamba/proc lock.
# APPTAINERENV_/SINGULARITYENV_ inject the var INTO the container.
cache="$OUT/.cache/prep_${TASK}"
export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

binds=(--bind "$MITO_SV_INPUT_DIR" --bind "$OUT")
ppargs=(--input "$path" --sample "$name" --outdir "$pdir" --threads "$THREADS")
if [[ -n "${MITO_SV_REFERENCE:-}" ]]; then
    binds+=(--bind "$(dirname "$MITO_SV_REFERENCE")"); ppargs+=(--reference "$MITO_SV_REFERENCE")
fi

echo "[$(date)] prep task=$TASK sample=$name input=$path threads=$THREADS"
"$RT" exec "${binds[@]}" "$MITO_SV_SIF" /opt/pipeline/preprocess.sh "${ppargs[@]}"
echo "[$(date)] prep task=$TASK sample=$name -> $pdir"
