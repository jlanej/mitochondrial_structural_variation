#!/usr/bin/env bash
###############################################################################
# run_eklipse.sh — eKLIPse (Python 2, BLAST+, circos) in the py2tools env.
# Input : normalised rCRS chrM BAM (indexed). Output: eKLIPse_deletions.csv
###############################################################################
set -euo pipefail
log() { printf '[eklipse %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

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

mkdir -p "$OUTDIR" "$OUTDIR/tmp"
bam_abs="$(readlink -f "$BAM")"
out_abs="$(readlink -f "$OUTDIR")"
list="$out_abs/bam_list.tsv"
# eKLIPse -in format: one alignment per line, "<bam-path>\t<title>"
printf '%s\t%s\n' "$bam_abs" "$SAMPLE" > "$list"

log "running eKLIPse on $SAMPLE"
# Run from the eKLIPse dir so its local modules (pybam/spinner/tabulate) resolve.
micromamba run -n py2tools bash -c "cd /opt/eKLIPse && python eKLIPse.py \
    -in '$list' \
    -ref /opt/eKLIPse/data/NC_012920.1.gb \
    -out '$out_abs' -tmp '$out_abs/tmp' -thread '$THREADS' --nocolor" \
    || log "eKLIPse returned non-zero"

# eKLIPse writes into a fresh eKLIPse_<uuid>/ subdir; surface stable names.
res="$(find "$out_abs" -name eKLIPse_deletions.csv | head -1 || true)"
gen="$(find "$out_abs" -name eKLIPse_genes.csv | head -1 || true)"
[[ -n "$gen" ]] && cp "$gen" "$out_abs/eKLIPse_genes.csv"
if [[ -n "$res" ]]; then
    cp "$res" "$out_abs/eKLIPse_deletions.csv"; log "deletions -> eKLIPse_deletions.csv"; log "done"
else
    log "WARNING: no eKLIPse_deletions.csv produced"; exit 1
fi
