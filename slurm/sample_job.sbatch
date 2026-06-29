#!/usr/bin/env bash
#SBATCH --job-name=mito-sv
###############################################################################
# sample_job.sbatch — one array task = one sample x ONE caller.
#
# Submitted by run_mito_sv.sh once per caller (MITO_SV_CALLER), --array=1-N with
# aftercorr on the prep array so each sample starts as soon as ITS prep is done.
# Consumes the prepared inputs (run_sample --prepared-bam) and writes to the
# isolated $OUT/by_caller/<caller>/<sample>/, so a slow caller (eKLIPse) never
# blocks a fast one and per-caller results land independently.
###############################################################################
set -euo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"; : "${MITO_SV_MANIFEST:?}"
: "${MITO_SV_CALLER:?MITO_SV_CALLER (a single caller) required}"
THREADS="${MITO_SV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
TASK="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job}"
OUT="$MITO_SV_OUTDIR"; CALLER="$MITO_SV_CALLER"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi

line="$(sed -n "${TASK}p" "$MITO_SV_MANIFEST")"
[[ -n "$line" ]] || { echo "ERROR: no manifest line $TASK" >&2; exit 1; }
name="${line%%$'\t'*}"

pdir="$OUT/prepared/$name"
bam="$pdir/${name}.chrM.bam"; r1="$pdir/${name}.R1.fastq.gz"; r2="$pdir/${name}.R2.fastq.gz"
[[ -f "$bam" && -f "$r1" && -f "$r2" ]] \
    || { echo "ERROR: prepared inputs missing for $name (prep failed?)" >&2; exit 1; }

sdir="$OUT/by_caller/$CALLER/$name"; mkdir -p "$sdir"
cache="$OUT/.cache/${CALLER}_${TASK}"
export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

args=(--prepared-bam "$bam" --prepared-r1 "$r1" --prepared-r2 "$r2"
      --sample "$name" --outdir "$sdir" --threads "$THREADS" --callers "$CALLER")
[[ "${MITO_SV_STRICT:-0}" == "1" ]] && args+=(--strict)

echo "[$(date)] task=$TASK caller=$CALLER sample=$name threads=$THREADS"
"$RT" exec --bind "$OUT" "$MITO_SV_SIF" /opt/pipeline/run_sample.sh "${args[@]}"
echo "[$(date)] task=$TASK caller=$CALLER sample=$name finished"
