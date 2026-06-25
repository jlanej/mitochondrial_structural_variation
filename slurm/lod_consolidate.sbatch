#!/usr/bin/env bash
#SBATCH --job-name=mito-lod-report
###############################################################################
# lod_consolidate.sbatch — aggregate the LOD shards into the analysis + report.
#
# Concatenates the per-cell shard TSVs, fits the LOD curves, and renders the
# interactive report to $MITO_SV_LOD_OUTDIR/lod_report/index.html.
###############################################################################
set -euo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_LOD_OUTDIR:?}"
OUT="$MITO_SV_LOD_OUTDIR"
ARMS="${MITO_SV_LOD_ARMS:-pipeline,circular}"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: no apptainer/singularity" >&2; exit 1; fi
PYEXE=/opt/conda/envs/mitosv/bin/python

shopt -s nullglob
shards=("$OUT"/shards/shard_*.tsv)
[[ "${#shards[@]}" -gt 0 ]] || { echo "ERROR: no shards in $OUT/shards" >&2; exit 1; }
echo "[lod_consolidate] merging ${#shards[@]} shards"
# keep the first header only
awk 'FNR==1 && seen++ {next} {print}' "${shards[@]}" > "$OUT/sweep.tsv"
ncells="$(($(wc -l < "$OUT/sweep.tsv") - 1))"
echo "[lod_consolidate] sweep.tsv: $ncells caller-rows"

mkdir -p "$OUT/lod_report"
"$RT" exec --bind "$OUT" "$MITO_SV_SIF" "$PYEXE" /opt/pipeline/lod/analyze_lod.py \
    --sweep "$OUT/sweep.tsv" --outdir "$OUT/lod_report"
"$RT" exec --bind "$OUT" "$MITO_SV_SIF" "$PYEXE" /opt/pipeline/lod/make_lod_report.py \
    --sweep "$OUT/sweep.tsv" --cells "$OUT/lod_report/lod_cells.tsv" \
    --fits "$OUT/lod_report/lod_fits.tsv" --arms "$ARMS" \
    --out "$OUT/lod_report/index.html"
echo "[lod_consolidate] report -> $OUT/lod_report/index.html"
