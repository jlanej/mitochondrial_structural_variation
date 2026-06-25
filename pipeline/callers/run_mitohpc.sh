#!/usr/bin/env bash
###############################################################################
# run_mitohpc.sh — MitoHPC reference SV caller (Python3 + pysam) in the mitosv env.
# Input : normalised rCRS chrM BAM (indexed). Output: <sample>.sv.tab
###############################################################################
set -euo pipefail
log() { printf '[mitohpc %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

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
prefix="$out_abs/$SAMPLE"

log "running MitoHPC callSV on $SAMPLE"
# callSV.sh reads HP_MT/HP_MTLEN to pick the contig + its reference; pin to our
# normalised rCRS chrM (16569). HP_PYTHON=python -> the mitosv env's pysam.
# HP_SV_PLOT='' explicitly DISABLES the optional samplot visualization
# (callSV.sh only plots when HP_SV_PLOT is non-empty) so the measured runtime is
# pure SV calling — comparable to the other callers, never inflated by plotting.
micromamba run -n mitosv bash -c "HP_SDIR=/opt/MitoHPC/scripts HP_RDIR=/opt/MitoHPC/RefSeq \
    HP_PYTHON=python HP_MT=chrM HP_MTLEN=16569 HP_SV_PLOT= \
    bash /opt/MitoHPC/scripts/callSV.sh '$SAMPLE' '$bam_abs' '$prefix'" \
    || log "callSV.sh returned non-zero"

if [[ -f "$prefix.sv.tab" ]]; then
    cp "$prefix.sv.tab" "$out_abs/mitohpc.sv.tab"
    [[ -f "$prefix.sv.vcf" ]] && cp "$prefix.sv.vcf" "$out_abs/mitohpc.sv.vcf"
    log "calls -> mitohpc.sv.tab"; log "done"
else
    log "WARNING: no .sv.tab produced"; exit 1
fi
