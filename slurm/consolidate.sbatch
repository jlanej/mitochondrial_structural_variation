#!/usr/bin/env bash
#SBATCH --job-name=mito-sv-consolidate
###############################################################################
# consolidate.sbatch — cohort post-processing + interactive HTML summary.
#
# Fires as a cascade: once after each caller finishes ALL its CRAMs (cons-<caller>)
# and once at the very end (mito-final, MITO_SV_CONS_FINAL=1). Each run rebuilds
# the cumulative cohort_*.tsv tables AND the global cohort_sv_summary.html, so
# results appear as soon as the first caller completes. Concurrent consolidations
# are flock-serialised so they never race on the shared outputs.
###############################################################################
set -uo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"
OUT="$MITO_SV_OUTDIR"
SCOPE="${MITO_SV_CONS_SCOPE:-cohort}"
FINAL="${MITO_SV_CONS_FINAL:-0}"
JOBID="${SLURM_JOB_ID:-$$}"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi
PY=/opt/conda/envs/mitosv/bin/python

cache="$OUT/.cache/cons_${JOBID}"
export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

run() {
    echo "[$(date)] consolidate scope=$SCOPE -> $OUT"
    "$RT" exec --bind "$OUT" "$MITO_SV_SIF" "$PY" /opt/pipeline/postprocess.py --root "$OUT" || true
    "$RT" exec --bind "$OUT" "$MITO_SV_SIF" "$PY" /opt/pipeline/cohort_report.py \
        --root "$OUT" --out "$OUT/cohort_sv_summary.html" \
        --generated "$(date +%F)" || true
    echo "[$(date)] refreshed cohort_*.tsv + cohort_sv_summary.html (scope=$SCOPE)"
}

if command -v flock >/dev/null 2>&1; then
    exec 9>"$OUT/.consolidate.lock"
    if flock -w 1800 9; then run; flock -u 9 2>/dev/null || true
    else echo "[consolidate] render lock busy 30m; skipping (a later caller will refresh)" >&2; fi
else
    run
fi

# the final consolidation drops the bulky prepared inputs + per-task caches
if [[ "$FINAL" == "1" ]]; then
    echo "[$(date)] final consolidation: removing prepared inputs + caches"
    rm -rf "$OUT/prepared" "$OUT/.cache" 2>/dev/null || true
    ls -1 "$OUT"/cohort_*.tsv "$OUT/cohort_sv_summary.html" 2>/dev/null || true
fi
