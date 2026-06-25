#!/usr/bin/env bash
###############################################################################
# lod_smoke.sh — single-iteration LOD gate for CI.
#
#   bash test/lod_smoke.sh [IMAGE]
#
# Proves the LOD machinery runs end to end on ONE tiny cell: simulate a del4977
# BAM, run all callers under both input arms (pipeline + circular), score, then
# aggregate + render the report. The FULL sweep is the HPC job (slurm/run_lod.sh).
#
# HARD gate: the per-cell shard is produced with the expected rows and the
# analysis + report render without error. Detection by >=1 caller is reported
# (not gated — heuristic callers on a tiny single cell).
###############################################################################
set -uo pipefail
IMAGE="${1:-mito-sv:ci}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT" 2>/dev/null || sudo -n rm -rf "$OUT" 2>/dev/null || true' EXIT
fail=0
note() { printf '\n=== %s ===\n' "$*"; }
err()  { printf 'FAIL: %s\n' "$*" >&2; fail=1; }

echo "image: $IMAGE"; echo "out:   $OUT"

note "running one LOD cell in the container (del4977, vaf=0.5, depth=200, both arms)"
docker run --rm -v "$OUT:/out" --entrypoint bash "$IMAGE" -lc '
set -e
trap "chmod -R a+rwX /out 2>/dev/null || true" EXIT
/opt/pipeline/lod/run_cell.sh --variant del4977 --bp5 8469 --bp3 13447 \
    --vaf 0.5 --depth 200 --rep 0 --outdir /out/work --shard /out/sweep.tsv \
    --threads 2 --arms pipeline,circular
' || err "container LOD cell returned non-zero"

note "shard"
if [[ -f "$OUT/sweep.tsv" ]]; then
    cat "$OUT/sweep.tsv" | cut -f1-11 | column -t || cat "$OUT/sweep.tsv"
    rows="$(($(wc -l < "$OUT/sweep.tsv") - 1))"
    # 6 callers x 2 arms = 12 rows
    [[ "$rows" -eq 12 ]] || err "expected 12 shard rows (6 callers x 2 arms), got $rows"
    det="$(awk -F'\t' 'NR>1 && $11==1' "$OUT/sweep.tsv" | wc -l | tr -d ' ')"
    echo "callers detecting the deletion: $det / $rows"
    [[ "$det" -ge 1 ]] || echo "WARNING: no caller detected del4977 in this single cell" >&2
else
    err "no sweep.tsv produced"
fi

note "analyze + report (host)"
if [[ -f "$OUT/sweep.tsv" ]]; then
    python3 "$REPO/pipeline/lod/analyze_lod.py" --sweep "$OUT/sweep.tsv" --outdir "$OUT/an" \
        || err "analyze_lod failed"
    python3 "$REPO/pipeline/lod/make_lod_report.py" --sweep "$OUT/sweep.tsv" \
        --cells "$OUT/an/lod_cells.tsv" --fits "$OUT/an/lod_fits.tsv" \
        --out "$OUT/an/index.html" --generated "ci-single-iteration" \
        || err "make_lod_report failed"
    [[ -s "$OUT/an/lod_cells.tsv" && -s "$OUT/an/index.html" ]] \
        || err "LOD analysis/report artifacts missing"
    echo "lod_cells:"; head -3 "$OUT/an/lod_cells.tsv"
fi

note "result"
if [[ "$fail" == 0 ]]; then echo "LOD SMOKE PASSED"; exit 0; else echo "LOD SMOKE FAILED"; exit 1; fi
