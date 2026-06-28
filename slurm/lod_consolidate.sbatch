#!/usr/bin/env bash
#SBATCH --job-name=mito-lod-report
###############################################################################
# lod_consolidate.sbatch — aggregate the LOD shards into the analysis + report.
#
# Fires many times during a sweep: once after each (caller,depth) array, once
# after each caller's depths finish, and once at the very end. Each run
# concatenates ALL shards present so far, fits the LOD curves, and re-renders the
# cumulative report — so results appear as soon as the first caller completes a
# depth. When MITO_SV_LOD_SCOPE_CALLER is set, it also prints that scope's runtime
# distribution and a rough projection of the remaining wall-time, the realtime
# signal for deciding whether to cancel the rest of the sweep.
#
# Concurrent consolidations serialize the analyse+render via flock so they never
# race on lod_report/; the scope runtime read-out always runs (it only reads).
###############################################################################
set -uo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_LOD_OUTDIR:?}"
OUT="$MITO_SV_LOD_OUTDIR"
ARMS="${MITO_SV_LOD_ARMS:-pipeline,circular}"
SCOPE_CALLER="${MITO_SV_LOD_SCOPE_CALLER:-}"
SCOPE_DEPTH="${MITO_SV_LOD_SCOPE_DEPTH:-}"
SCOPE_EXPECT="${MITO_SV_LOD_SCOPE_EXPECT:-0}"   # expected datapoints (cells x arms) for the scope
CONC="${MITO_SV_LOD_CONC:-12}"
JOBID="${SLURM_JOB_ID:-$$}"

# Isolate this job's micromamba cache so overlapping consolidations don't contend
# on the shared $HOME 'mamba/proc' lock (same fix as the cells).
export XDG_CACHE_HOME="$OUT/.cache/cons_$JOBID"; mkdir -p "$XDG_CACHE_HOME"
trap 'rm -rf "$XDG_CACHE_HOME"' EXIT

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: no apptainer/singularity" >&2; exit 1; fi
PYEXE=/opt/conda/envs/mitosv/bin/python

scope_label="${SCOPE_CALLER:+$SCOPE_CALLER}${SCOPE_DEPTH:+ @ depth $SCOPE_DEPTH}"
echo "[lod_consolidate] ${scope_label:-FULL} — $(date '+%F %T')"

consolidate() {   # the racy part: concat shards -> analyse -> render (under flock)
    shopt -s nullglob
    local shards=("$OUT"/shards/shard_*.tsv)
    [[ "${#shards[@]}" -gt 0 ]] || { echo "[lod_consolidate] no shards yet" >&2; return 1; }
    echo "[lod_consolidate] merging ${#shards[@]} shards"
    # Keep the first header only; drop any partial trailing rows from shards that a
    # still-running array is mid-write (intermediate consolidations run concurrently
    # with active writers). The header's column count is the gate, so a truncated
    # line (fewer fields) is skipped rather than crashing analyze on int("").
    awk -F'\t' 'FNR==1 && !seen {ncol=NF; seen=1; print; next}
                FNR==1 {next}
                NF==ncol {print}' "${shards[@]}" > "$OUT/sweep.tsv.tmp"
    mv -f "$OUT/sweep.tsv.tmp" "$OUT/sweep.tsv"
    local ncells; ncells="$(($(wc -l < "$OUT/sweep.tsv") - 1))"
    echo "[lod_consolidate] sweep.tsv: $ncells caller-rows"
    mkdir -p "$OUT/lod_report"
    "$RT" exec --bind "$OUT" "$MITO_SV_SIF" "$PYEXE" /opt/pipeline/lod/analyze_lod.py \
        --sweep "$OUT/sweep.tsv" --outdir "$OUT/lod_report"
    "$RT" exec --bind "$OUT" "$MITO_SV_SIF" "$PYEXE" /opt/pipeline/lod/make_lod_report.py \
        --sweep "$OUT/sweep.tsv" --cells "$OUT/lod_report/lod_cells.tsv" \
        --fits "$OUT/lod_report/lod_fits.tsv" --arms "$ARMS" \
        --out "$OUT/lod_report/index.html"
    echo "[lod_consolidate] report -> $OUT/lod_report/index.html"
}

if command -v flock >/dev/null 2>&1; then
    exec 9>"$OUT/.consolidate.lock"
    if flock -w 1800 9; then
        consolidate || true
        flock -u 9 2>/dev/null || true
    else
        echo "[lod_consolidate] render lock busy 30m; skipping render, printing runtime from last report" >&2
    fi
else
    consolidate || true
fi

# --- scope runtime read-out + projection (always runs; read-only) ------------
rt="$OUT/lod_report/lod_runtime.tsv"
if [[ -n "$SCOPE_CALLER" && -f "$rt" ]]; then
    echo "[lod_consolidate] runtime so far — $scope_label (both arms combined):"
    # lod_runtime.tsv cols: arm caller depth n mean median p25 p75 min max total
    awk -F'\t' -v c="$SCOPE_CALLER" -v d="$SCOPE_DEPTH" '
        $1=="all" && $2==c && (d=="" || $3==d) && $3!="all" {
            printf "    depth %-5s  n=%-4s  mean=%6.0fs  median=%6.0fs  p25=%-5s p75=%-5s  max=%6.0fs\n",
                   $3, $4, $5, $6, $7, $8, $10
        }' "$rt"
    # overall scope mean + projection (mean across the scope rows, weighted by n).
    # NB: 'exp' is an awk builtin (exponential) -- use 'xpct' for the expected count.
    awk -F'\t' -v c="$SCOPE_CALLER" -v d="$SCOPE_DEPTH" -v xpct="$SCOPE_EXPECT" -v conc="$CONC" '
        $1=="all" && $2==c && (d=="" || $3==d) && $3!="all" { n+=$4; s+=$5*$4 }
        END {
            if (n>0) {
                m=s/n; rem=(xpct+0>n)?xpct-n:0;
                printf "    => %d/%s datapoints; weighted mean %.0fs/invocation\n", n, (xpct+0?xpct:"?"), m
                if (xpct+0>0) printf "    => rough remaining wall @ conc=%s: ~%.1f h (%.0f datapoints left)\n", conc, rem*m/conc/3600.0, rem
            }
        }' "$rt"
fi
