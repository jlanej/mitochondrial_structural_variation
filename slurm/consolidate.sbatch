#!/usr/bin/env bash
#SBATCH --job-name=mito-sv-consolidate
###############################################################################
# consolidate.sbatch — cohort post-processing after all sample jobs finish.
#
# Submitted by run_mito_sv.sh with --dependency=afterany on the sample array.
# Parses every sample's caller outputs and writes $MITO_SV_OUTDIR/cohort_*.tsv.
###############################################################################
set -euo pipefail

: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi

echo "[$(date)] consolidating cohort under $MITO_SV_OUTDIR"
"$RT" exec --bind "$MITO_SV_OUTDIR" "$MITO_SV_SIF" \
    /opt/conda/envs/mitosv/bin/python /opt/pipeline/postprocess.py --root "$MITO_SV_OUTDIR"
echo "[$(date)] cohort summary written:"
ls -1 "$MITO_SV_OUTDIR"/cohort_*.tsv "$MITO_SV_OUTDIR"/cohort_summary.txt 2>/dev/null || true
