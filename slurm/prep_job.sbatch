#!/usr/bin/env bash
#SBATCH --job-name=mito-prep
###############################################################################
# prep_job.sbatch — one array task = realign ONE sample's cached chrM slice to
# rCRS -> normalised chrM BAM + mito FASTQ in $OUT/prepared/<sample>/.
#
# Reads the small $OUT/chrM/<sample>.chrM.bam produced by the mito-extract stage
# (NOT the original WGS CRAM — that is decoded only once, in extract). The
# per-caller arrays then consume these prepared inputs via run_sample
# --prepared-bam, so the CRAM is never re-touched.
###############################################################################
set -euo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"; : "${MITO_SV_MANIFEST:?}"
THREADS="${MITO_SV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
TASK="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job}"
OUT="$MITO_SV_OUTDIR"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi

line="$(sed -n "${TASK}p" "$MITO_SV_MANIFEST")"
[[ -n "$line" ]] || { echo "ERROR: no manifest line $TASK" >&2; exit 1; }
name="${line%%$'\t'*}"
chrm="$OUT/chrM/${name}.chrM.bam"
[[ -f "$chrm" ]] || { echo "ERROR: extracted chrM missing for $name ($chrm; extract failed?)" >&2; exit 1; }
pdir="$OUT/prepared/$name"; mkdir -p "$pdir"

# Isolate this task's micromamba cache: Apptainer binds the shared host $HOME, so
# without isolation every concurrent micromamba serialises on one mamba/proc lock.
cache="$OUT/.cache/prep_${TASK}"
export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

echo "[$(date)] prep task=$TASK sample=$name input=$chrm threads=$THREADS"
"$RT" exec --bind "$OUT" "$MITO_SV_SIF" /opt/pipeline/preprocess.sh \
    --input "$chrm" --sample "$name" --outdir "$pdir" --threads "$THREADS"
echo "[$(date)] prep task=$TASK sample=$name -> $pdir"
