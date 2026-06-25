#!/usr/bin/env bash
###############################################################################
# run_cell.sh — one LOD grid cell, end to end (runs INSIDE the container).
#
#   1. generate the cell BAM (simulate WT+deletion reads -> circular-aware BAM)
#   2. run all callers under each requested input arm:
#        pipeline : run_sample.sh re-normalises with bwa mem (production path)
#        circular : callers run directly on MitoHPC's circSam BAM (--prepared)
#   3. score each (arm, caller) vs truth -> append rows to --shard
#
# Heavy intermediates are deleted after scoring (the sweep is large); only the
# tiny TSV shard survives.
###############################################################################
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPE="$(cd "$HERE/.." && pwd)"
log() { printf '[run_cell %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
PY() { micromamba run -n mitosv python "$@"; }

VARIANT="" BP5="" BP3="" VAF="" DEPTH="" REP="" OUTDIR="" SHARD="" THREADS=4
ARMS="pipeline,circular" CALLERS="all" KEEP=0
while [[ $# -gt 0 ]]; do case "$1" in
    --variant) VARIANT="$2"; shift 2;;
    --bp5) BP5="$2"; shift 2;;
    --bp3) BP3="$2"; shift 2;;
    --vaf) VAF="$2"; shift 2;;
    --depth) DEPTH="$2"; shift 2;;
    --rep) REP="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --shard) SHARD="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --arms) ARMS="$2"; shift 2;;
    --callers) CALLERS="$2"; shift 2;;
    --keep) KEEP=1; shift;;
    *) log "unknown arg: $1"; exit 2;;
esac; done
: "${VARIANT:?}"; : "${BP5:?}"; : "${BP3:?}"; : "${VAF:?}"; : "${DEPTH:?}"; : "${REP:?}"
: "${OUTDIR:?}"; : "${SHARD:?}"

vtag="$(printf '%s' "$VAF" | tr '.' 'p')"
cell="${VARIANT}_v${vtag}_d${DEPTH}_r${REP}"
work="$OUTDIR/$cell"
mkdir -p "$work" "$(dirname "$SHARD")"
log "cell=$cell arms=$ARMS threads=$THREADS"

# 1. generate
gen="$work/gen"
if ! PY "$PIPE/lod/gen_cell.py" --variant "$VARIANT" --bp5 "$BP5" --bp3 "$BP3" \
        --vaf "$VAF" --depth "$DEPTH" --rep "$REP" --out "$gen" --threads "$THREADS"; then
    log "ERROR: generation failed for $cell"; exit 1
fi
bam="$gen/cell.bam"; r1="$gen/cell.r1.fastq.gz"; r2="$gen/cell.r2.fastq.gz"
truth="$gen/truth.tsv"

# 2 + 3. arms
IFS=',' read -ra ARMA <<< "$ARMS"
for arm in "${ARMA[@]}"; do
    adir="$work/$arm"
    if [[ "$arm" == "pipeline" ]]; then
        bash "$PIPE/run_sample.sh" --input "$bam" --sample "$cell" \
            --outdir "$adir" --threads "$THREADS" --callers "$CALLERS" >"$work/${arm}.log" 2>&1 || true
    elif [[ "$arm" == "circular" ]]; then
        bash "$PIPE/run_sample.sh" --prepared-bam "$bam" --prepared-r1 "$r1" --prepared-r2 "$r2" \
            --sample "$cell" --outdir "$adir" --threads "$THREADS" --callers "$CALLERS" >"$work/${arm}.log" 2>&1 || true
    else
        log "unknown arm: $arm"; continue
    fi
    PY "$PIPE/lod/score_cell.py" --sample-dir "$adir" --truth "$truth" \
        --arm "$arm" --sample "$cell" --out "$SHARD" || log "scoring failed for $arm"
done

[[ "$KEEP" == 1 ]] || rm -rf "$work"
log "done $cell -> shard $SHARD"
