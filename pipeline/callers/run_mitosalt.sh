#!/usr/bin/env bash
###############################################################################
# run_mitosalt.sh — MitoSAlt (Perl + R, LAST) in the mitosalt env.
# Input : mito FASTQ pair. Runs in "enriched" mode (LAST MT index only).
# Output: indel/<sample>.tsv  (deletion/duplication calls)
#
# MitoSAlt is CWD-relative and writes into bam/ bw/ tab/ indel/ log/ plot/ next
# to its scripts, so we stage a writable working dir per sample with symlinks to
# the install's delplot.R / perl script / genome index.
###############################################################################
set -euo pipefail
log() { printf '[mitosalt %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

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
: "${SAMPLE:?}"; : "${OUTDIR:?}"; : "${R1:?}"; : "${R2:?}"

mkdir -p "$OUTDIR"
out_abs="$(readlink -f "$OUTDIR")"
r1_abs="$(readlink -f "$R1")"
r2_abs="$(readlink -f "$R2")"
tag="$SAMPLE"

work="$out_abs/work"
mkdir -p "$work"/{bam,bw,tab,indel,log,plot,bin}
ln -sf /opt/MitoSAlt/delplot.R        "$work/delplot.R"
ln -sf /opt/MitoSAlt/MitoSAlt1.1.1.pl "$work/MitoSAlt1.1.1.pl"
ln -sf /opt/MitoSAlt/genome           "$work/genome"
# Per-job thread count.
sed "s/^threads = .*/threads = ${THREADS}/" /opt/MitoSAlt/config_pipeline.txt > "$work/config.txt"

log "running MitoSAlt on $SAMPLE"
mitosalt_ok=1
micromamba run -n mitosalt bash -c "cd '$work' && MITOSALT_XMX='${MITOSALT_XMX:-4g}' \
    perl MitoSAlt1.1.1.pl config.txt '$r1_abs' '$r2_abs' '$tag'" \
    || { mitosalt_ok=0; log "MitoSAlt returned non-zero"; }

# delplot.R only writes indel/<tag>.tsv when deletions are actually called, so a
# missing .tsv on a *completed* run means "no calls", not failure. The log ends
# with ':Finished' on a clean completion.
finished=0
grep -q ':Finished' "$work/log/${tag}.log" 2>/dev/null && finished=1

if [[ -f "$work/indel/${tag}.tsv" ]]; then
    cp "$work/indel/${tag}.tsv" "$out_abs/${tag}.mitosalt.tsv"
    log "calls -> ${tag}.mitosalt.tsv"; log "done"
elif [[ "$finished" == 1 || "$mitosalt_ok" == 1 ]]; then
    # Ran to completion with zero deletion clusters: emit an empty (header-only)
    # result so this is recorded as "operated, no calls" rather than a failure.
    printf 'sample\tcluster.id\talt.reads\tref.reads\theteroplasmy\tdel.start.range\tdel.end.range\tdel.size\tfinal.event\tfinal.start\tfinal.end\tfinal.size\tseq1\tseq2\tseq\n' \
        > "$out_abs/${tag}.mitosalt.tsv"
    log "no deletion clusters; wrote empty ${tag}.mitosalt.tsv"; log "done"
else
    log "WARNING: MitoSAlt did not complete and produced no ${tag}.tsv"; exit 1
fi
