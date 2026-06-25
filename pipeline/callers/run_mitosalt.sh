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

# --- under-the-hood diagnostics --------------------------------------------
# Surface where the read evidence was lost, EVERY run (not only failures), and
# leave a machine-readable sidecar. split alignments ~0 => LAST found no split
# reads; clusters < cluster_threshold => evidence too sparse to call; a non-empty
# breakpoint/cluster but empty .tsv => delplot.R filtered/errored (see Rout).
_cnt() { [[ -f "$1" ]] && { wc -l < "$1" | tr -d ' '; } || echo 0; }
tabz="$work/tab/${tag}.tab.gz"
n_split=$( [[ -f "$tabz" ]] && { zcat -f "$tabz" 2>/dev/null | grep -cv '^#' || true; } || echo 0 )
n_bp=$(_cnt "$work/indel/${tag}.breakpoint")
n_clu=$(_cnt "$work/indel/${tag}.cluster")
n_call=$(_cnt "$work/indel/${tag}.tsv"); [[ "$n_call" -gt 0 ]] && n_call=$((n_call - 1))
# Forensic signals for the two suspected drop points (see review):
#  * paired-name fraction: build_hash needs query names ending /1 or /2 (col 7);
#    if reformat addslash failed, this is ~0 and mates mis-bin -> no breakpoints.
#  * low-score arms: each split arm must clear score_threshold(=80)/lastal -e80;
#    a short arm of a breakpoint-spanning 150bp read can score below it and be
#    dropped, so count how many alignment rows score < 80.
named="?"; lowscore="?"
if [[ -f "$tabz" ]]; then
    tot="$(zcat -f "$tabz" 2>/dev/null | grep -cv '^#' || echo 0)"
    sl="$(zcat -f "$tabz" 2>/dev/null | awk '$1!~/^#/ && $7~/\/[12]$/' | wc -l | tr -d ' ')"
    ls="$(zcat -f "$tabz" 2>/dev/null | awk '$1!~/^#/ && $1<80' | wc -l | tr -d ' ')"
    named="${sl}/${tot}"; lowscore="${ls}/${tot}"
fi
{
    echo "----- MitoSAlt under-the-hood diagnostics ($tag) -----"
    echo "native clean exit   : $mitosalt_ok      (0 => perl returned non-zero)"
    echo "log reached :Finished: $finished"
    echo "split alignments    : $n_split   (tab/${tag}.tab.gz; ~0 => LAST emitted no split reads)"
    echo "  paired-name arms  : $named     (col7 ends /1 or /2; ~0 => pairing broken, mates mis-bin)"
    echo "  arms score < 80   : $lowscore  (killed by score_threshold/lastal -e80)"
    echo "breakpoints         : $n_bp      (indel/${tag}.breakpoint; split reads with both arms kept)"
    echo "clusters            : $n_clu     (indel/${tag}.cluster; need >= cluster_threshold to call)"
    echo "final calls         : $n_call    (indel/${tag}.tsv rows)"
} >&2
[[ -f "$work/log/${tag}.log" ]] && { echo "----- log/${tag}.log (tail 40) -----" >&2; tail -n 40 "$work/log/${tag}.log" >&2; }
[[ -f "$work/delplot.Rout" ]] && { echo "----- delplot.Rout (tail 30) -----" >&2; tail -n 30 "$work/delplot.Rout" >&2; }
printf 'caller\tnative_rc\tcompleted\tsplit_aln\tpaired_name_arms\tlowscore_arms\tbreakpoints\tclusters\tcalls\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    mitosalt "$mitosalt_ok" "$finished" "$n_split" "$named" "$lowscore" "$n_bp" "$n_clu" "$n_call" > "$out_abs/mitosalt.runstatus"
# ---------------------------------------------------------------------------

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
