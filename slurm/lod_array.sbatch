#!/usr/bin/env bash
#SBATCH --job-name=mito-lod
###############################################################################
# lod_array.sbatch — one LOD array task = a chunk of grid cells for ONE caller.
#
# Submitted by run_lod.sh once per (caller, depth) as --array=1-N. Each task
# processes its slice of that depth's cell manifest, running run_cell.sh
# concurrently (CONC cells x TPC threads ~ 24 cores) inside the Apptainer image
# with --callers "$CALLER" (a SINGLE caller), so a slow caller can never block a
# fast one. Everything is namespaced by SCOPE=<caller>_d<depth>: shards, the
# per-cell work tree, and the progress dir, so concurrent (caller,depth) arrays
# never collide. Each cell writes its OWN shard TSV (no concurrent-append races);
# lod_consolidate.sbatch concatenates all shards seen so far.
###############################################################################
set -uo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_LOD_OUTDIR:?}"; : "${MITO_SV_LOD_MANIFEST:?}"
CPJ="${MITO_SV_LOD_CPJ:-12}"; CONC="${MITO_SV_LOD_CONC:-12}"; TPC="${MITO_SV_LOD_TPC:-2}"
ARMS="${MITO_SV_LOD_ARMS:-pipeline,circular}"
CALLER="${MITO_SV_LOD_CALLER:-all}"
SCOPE="${MITO_SV_LOD_SCOPE:-all}"
HEARTBEAT_SEC="${MITO_SV_LOD_HEARTBEAT:-120}"
TASK="${SLURM_ARRAY_TASK_ID:-1}"
OUT="$MITO_SV_LOD_OUTDIR"
WORK="$OUT/work/$SCOPE"                 # per-scope cell work tree (no cross-caller clash)
PROG="$OUT/progress/$SCOPE/t$TASK"; mkdir -p "$PROG" "$WORK"; : > "$PROG/cells.done"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: no apptainer/singularity" >&2; exit 1; fi

start=$(( (TASK - 1) * CPJ + 1 )); end=$(( TASK * CPJ ))
mkdir -p "$OUT/shards"
wl="$(mktemp)"; i=0
while IFS=$'\t' read -r variant bp5 bp3 vaf depth rep; do
    [[ -n "$variant" ]] || continue
    i=$((i + 1))
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$variant" "$bp5" "$bp3" "$vaf" "$depth" "$rep" \
        "$OUT/shards/shard_${SCOPE}_t${TASK}_${i}.tsv" >> "$wl"
done < <(sed -n "${start},${end}p" "$MITO_SV_LOD_MANIFEST")
echo "[lod_array] $SCOPE task $TASK: lines $start-$end ($i cells), caller=$CALLER conc=$CONC tpc=$TPC arms=$ARMS"

# --- progress monitor: every HEARTBEAT_SEC, report done/total + each active cell's
#     elapsed time, current phase, and (within an arm) which callers have finished
#     so a stuck cell/caller is obvious. Runs in the background; killed at the end.
NCELLS="$i"
heartbeat() {
    while sleep "$HEARTBEAT_SEC"; do
        local now done; now="$(date +%s)"; done="$(wc -l < "$PROG/cells.done" 2>/dev/null || echo 0)"
        echo "[heartbeat $SCOPE task $TASK] $(date +%H:%M:%S) — ${done}/${NCELLS} cells done; active:"
        local any=0 f ts cell ph el st cur
        for f in "$PROG"/*.run; do
            [[ -e "$f" ]] || continue; any=1
            IFS=$'\t' read -r ts cell ph < "$f" 2>/dev/null || continue
            el=$(( now - ts )); cur=""
            if [[ "$ph" == pipeline || "$ph" == circular ]]; then
                st="$WORK/$cell/$ph/status.tsv"
                [[ -f "$st" ]] && cur="callers done: $(awk -F'\t' 'NR>1{printf "%s ",$1}' "$st" 2>/dev/null)"
            fi
            printf '    %6ds  %-30s phase=%-9s %s\n' "$el" "$cell" "$ph" "$cur"
        done
        [[ "$any" == 0 ]] && echo "    (none active)"
    done
}
heartbeat & HB=$!
trap 'kill "$HB" 2>/dev/null || true' EXIT

xargs -a "$wl" -P "$CONC" -L1 bash -c '
  "'"$RT"'" exec --bind "'"$OUT"'" "'"$MITO_SV_SIF"'" /opt/pipeline/lod/run_cell.sh \
    --variant "$1" --bp5 "$2" --bp3 "$3" --vaf "$4" --depth "$5" --rep "$6" \
    --shard "$7" --outdir "'"$WORK"'" --threads "'"$TPC"'" --arms "'"$ARMS"'" \
    --callers "'"$CALLER"'" --progress "'"$PROG"'"
' _
kill "$HB" 2>/dev/null || true
rm -f "$wl"

# Mean cell wall-time for THIS task (a first runtime read-out; the per-(caller,depth)
# consolidation prints the full picture). One depth per array, so this is the scope.
echo "[lod_array] $SCOPE task $TASK done: $(wc -l < "$PROG/cells.done") cells. Mean cell time:"
awk -F'\t' '{n++; s+=$5} END{if(n) printf "    %s: %6.0fs/cell mean (n=%d)\n", "'"$SCOPE"'", s/n, n}' \
    "$PROG/cells.done" 2>/dev/null
