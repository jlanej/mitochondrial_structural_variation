#!/usr/bin/env bash
#SBATCH --job-name=mito-lod
###############################################################################
# lod_array.sbatch — one LOD array task = a chunk of grid cells.
#
# Submitted by run_lod.sh with --array=1-N. Processes its slice of the cell
# manifest, running run_cell.sh concurrently (CONC cells x TPC threads ~ 24
# cores) inside the Apptainer image. Each cell writes its OWN shard TSV to avoid
# concurrent-append races; lod_consolidate.sbatch concatenates them.
###############################################################################
set -uo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_LOD_OUTDIR:?}"; : "${MITO_SV_LOD_MANIFEST:?}"
CPJ="${MITO_SV_LOD_CPJ:-30}"; CONC="${MITO_SV_LOD_CONC:-6}"; TPC="${MITO_SV_LOD_TPC:-4}"
ARMS="${MITO_SV_LOD_ARMS:-pipeline,circular}"
TASK="${SLURM_ARRAY_TASK_ID:-1}"
OUT="$MITO_SV_LOD_OUTDIR"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: no apptainer/singularity" >&2; exit 1; fi

start=$(( (TASK - 1) * CPJ + 1 )); end=$(( TASK * CPJ ))
mkdir -p "$OUT/shards" "$OUT/work"
wl="$(mktemp)"; i=0
while IFS=$'\t' read -r variant bp5 bp3 vaf depth rep; do
    [[ -n "$variant" ]] || continue
    i=$((i + 1))
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$variant" "$bp5" "$bp3" "$vaf" "$depth" "$rep" \
        "$OUT/shards/shard_${TASK}_${i}.tsv" >> "$wl"
done < <(sed -n "${start},${end}p" "$MITO_SV_LOD_MANIFEST")
echo "[lod_array] task $TASK: lines $start-$end ($i cells), conc=$CONC tpc=$TPC arms=$ARMS"

xargs -a "$wl" -P "$CONC" -L1 bash -c '
  "'"$RT"'" exec --bind "'"$OUT"'" "'"$MITO_SV_SIF"'" /opt/pipeline/lod/run_cell.sh \
    --variant "$1" --bp5 "$2" --bp3 "$3" --vaf "$4" --depth "$5" --rep "$6" \
    --shard "$7" --outdir "'"$OUT"'/work" --threads "'"$TPC"'" --arms "'"$ARMS"'"
' _
rm -f "$wl"
echo "[lod_array] task $TASK done"
