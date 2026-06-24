#!/usr/bin/env bash
###############################################################################
# run_mitomut.sh — MitoMut (Python 3 + pysam + UCSC BLAT) in the mitosv env.
# Input : normalised rCRS chrM BAM (sorted + indexed). Output: *_results.txt
###############################################################################
set -euo pipefail
log() { printf '[mitomut %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

SAMPLE="" BAM="" OUTDIR="" R1="" R2="" THREADS=4
while [[ $# -gt 0 ]]; do case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --r1) R1="$2"; shift 2;;
    --r2) R2="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) shift;;
esac; done
: "${SAMPLE:?}"; : "${BAM:?}"; : "${OUTDIR:?}"

mkdir -p "$OUTDIR"
bam_abs="$(readlink -f "$BAM")"
out_abs="$(readlink -f "$OUTDIR")"

log "running MitoMut on $SAMPLE"
# -c chrM matches our normalised contig; -f mt.fasta is bundled rCRS (header >chrM).
micromamba run -n mitosv bash -c "cd '$out_abs' && python /opt/MitoMut/MitoMut.py \
    -f /opt/MitoMut/mt.fasta -c chrM -l 16569 -d '$out_abs' '$bam_abs'" \
    || log "MitoMut returned non-zero"

res="$(find "$out_abs" -maxdepth 1 -name '*_results.txt' | head -1 || true)"
if [[ -n "$res" ]]; then
    cp "$res" "$out_abs/mitomut_results.txt"; log "calls -> mitomut_results.txt"; log "done"
else
    log "WARNING: no *_results.txt produced"; exit 1
fi
