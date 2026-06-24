#!/usr/bin/env bash
###############################################################################
# run_mitoseek.sh — MitoSeek (Perl, legacy samtools 0.1.x) in the mitoseek env.
# Input : normalised rCRS chrM BAM (indexed). Runs in mito-only mode (-t 4).
# Outputs: mito1_structure_discordant_mates.txt, mito1_structure_large_deletion.sam
###############################################################################
set -euo pipefail
log() { printf '[mitoseek %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

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

# MitoSeek is locked to samtools 0.1.x semantics; the install recorded the
# resolved legacy samtools (isolated conda env or the bundled 0.1.18 ELF).
sam018="$(cat /opt/MitoSeek/.samtools_path 2>/dev/null || true)"
[[ -n "$sam018" && -x "$sam018" ]] || sam018="/opt/MitoSeek/Resources/samtools/samtools"

log "running MitoSeek on $SAMPLE (samtools: $sam018)"
# -t 4 mitochondria-only; SV detection (discordant mates + large deletions) is
# always on; -noch/-nocs/-noQC skip circos + QC plots (avoid GD plotting paths).
# MitoSeek creates an output dir named after the BAM basename in the CWD.
micromamba run -n mitoseek bash -c "cd '$out_abs' && perl /opt/MitoSeek/mitoSeek.pl \
    -i '$bam_abs' -t 4 -r rCRS -R rCRS -str 4 -strf 500 -d 5 \
    -noch -nocs -noQC -samtools '$sam018'" || log "mitoSeek.pl returned non-zero"

disc="$(find "$out_abs" -name 'mito1_structure_discordant_mates.txt' | head -1 || true)"
delsam="$(find "$out_abs" -name 'mito1_structure_large_deletion.sam' | head -1 || true)"
[[ -n "$disc" ]]   && cp "$disc"   "$out_abs/mitoseek_discordant_mates.txt"
[[ -n "$delsam" ]] && cp "$delsam" "$out_abs/mitoseek_large_deletion.sam"
if [[ -n "$disc$delsam" ]]; then
    log "done"
else
    log "WARNING: no MitoSeek structure outputs produced"; exit 1
fi
