#!/usr/bin/env bash
###############################################################################
# run_splicebreak2.sh — Splice-Break2 (bash + Python2 + Java8, MapSplice2).
# Input : mito FASTQ pair (realigned internally to bundled rCRS).
# Output: <sample>_LargeMTDeletions_WGS-only_NoPositionFilter.txt
#
# Splice-Break2 needs a WRITABLE SB_Path (it writes/removes temp.log there) and
# reads FASTQs named <sample>.R1.fastq[.gz]/<sample>.R2.fastq[.gz] from an input
# dir, so we stage a per-sample writable copy of the install + input dir.
###############################################################################
set -euo pipefail
log() { printf '[splicebreak2 %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# NB: Splice-Break2 hardcodes MapSplice's thread count (-p 8) in its driver, so
# --threads is accepted for a uniform wrapper interface but is a no-op here.
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

sbsrc=/opt/Splice-Break2/Splice-Break2-v3.0.2_PAIRED-END
sbwork="$out_abs/sb_install"
rm -rf "$sbwork"
# The driver needs a WRITABLE SB_Path but only ever CREATES files in it (its
# fresh temp.log; everything else under SB_Path is read-only binaries/reference).
# So hard-link the 136 MB install (`cp -al`) instead of byte-copying it — instant
# and ~zero extra space, which matters hugely under the LOD sweep (Splice-Break2
# runs thousands of times). `cp -al` falls back to a real copy across filesystems;
# drop any hard-linked temp.log so an append can never reach the shared install.
cp -al "$sbsrc" "$sbwork" 2>/dev/null || cp -a "$sbsrc" "$sbwork"
rm -f "$sbwork/temp.log"
SB_PATH="$(readlink -f "$sbwork")"

indir="$out_abs/in"; sbout="$out_abs/out"; sblog="$out_abs/log"
mkdir -p "$indir" "$sbout" "$sblog"
# The driver expects <sample>.R1.fastq / <sample>.R2.fastq. Splice-Break2's
# dedupe->reformat->MapSplice pipeline (and the read-name cleanup seds in the
# driver, which expect the pattern "<name>/1_dd0 /1") requires reads to carry
# /1 and /2 mate suffixes; our preprocessing FASTQ has clean names (samtools
# fastq -n), which makes MapSplice reject the pair with "Base name of two ends
# not consistent". Re-add the mate suffixes for Splice-Break2's copy only
# (MitoSAlt still consumes the clean-named preprocessing FASTQ).
zcat -f "$r1_abs" | awk 'NR%4==1{print $1"/1"; next} {print}' > "$indir/${SAMPLE}.R1.fastq"
zcat -f "$r2_abs" | awk 'NR%4==1{print $1"/2"; next} {print}' > "$indir/${SAMPLE}.R2.fastq"

log "running Splice-Break2 on $SAMPLE"
# Runs in py2tools env: provides Java 8 + Python 2 (and /usr/bin/python -> py2).
# bbmap (dedupe/bbduk/reformat) autodetects a JVM heap and sets -Xms == -Xmx to
# ~85% of free RAM; on memory-limited runners the JVM fails to RESERVE that heap
# at startup and aborts in ~2-3s (its stderr is swallowed into $sblog). Cap it
# via RQCMEM (MB) — bbmap's calcmem.sh uses RQCMEM to bound -Xmx/-Xms for every
# bbmap tool. Override with SB_MEM_MB for big cohorts on roomy HPC nodes.
SB_MEM_MB="${SB_MEM_MB:-2000}"
sb_rc=0
micromamba run -n py2tools bash -c "export RQCMEM='$SB_MEM_MB'; bash '$SB_PATH/Splice-Break2_paired-end.sh' \
    '$indir' '$sbout' '$sblog' '$SB_PATH' \
    --align=yes --ref=rCRS --fastq_keep=no --skip_preAlign=no" \
    || { sb_rc=$?; log "driver returned non-zero (rc=$sb_rc)"; }

res="$(find "$sbout" -name '*_LargeMTDeletions_WGS-only_NoPositionFilter.txt' | head -1 || true)"

# --- under-the-hood diagnostics (every run, not only failures) --------------
# MapSplice may run cleanly yet the deletion never reaches the LargeMTDeletions
# table; capture the junction count and the deletion-row count so a CI run shows
# exactly where it dropped, plus the driver / MapSplice / bbmap logs.
junc="$(find "$sbout" -name 'junctions.txt' -o -name '*junction*.txt' 2>/dev/null | head -1 || true)"
n_junc=$( [[ -n "$junc" && -f "$junc" ]] && { wc -l < "$junc" | tr -d ' '; } || echo 0 )
# del4977 junction present in MapSplice's post-filter junctions.txt? (rCRS bp ~8470/13447)
n_del4977j=$( [[ -n "$junc" && -f "$junc" ]] && { awk '$2<=8600 && $3>=13300 && $3-$2>4000' "$junc" 2>/dev/null | wc -l | tr -d ' '; } || echo 0 )
res_bytes=$( [[ -n "$res" && -f "$res" ]] && wc -c < "$res" | tr -d ' ' || echo 0 )
n_del=$( [[ -n "$res" && -f "$res" && "$res_bytes" -gt 0 ]] && { c=$(wc -l < "$res"); echo $(( c > 1 ? c - 1 : 0 )); } || echo 0 )
# The inner Splice-Break2_0725.sh logs an ordered set of progress markers and
# exits at the first failing stage; the LAST marker reached localizes the drop
# (its pipeline logic is verified to work on Linux, so the fault is a stage that
# exits/empties on our input). Capture the last stage, the failed-step count, and
# the key intermediate sizes so a single CI run pinpoints it.
nohup="$(find "$sbout" -name '*_nohup.log' 2>/dev/null | head -1 || true)"
last_stage="?"; steps_failed="?"
if [[ -n "$nohup" && -f "$nohup" ]]; then
    last_stage="$(grep -oE 'Running Samtools|Running CountBases.py|Adding per base coverage|Coverage with (rCRS|Nsub) reference|Using (average|left)_benchmark|Annotating|Applying Filter [12]|Create LargeMTDeletions[^ ]*|Annotation with Impact of Gene|Filtering Top 30[^ ]*|Junction Filter Steps Success' "$nohup" 2>/dev/null | tail -1 | tr ' ' '_')"
    steps_failed="$(grep -oE '[0-9]+ steps failed' "$nohup" 2>/dev/null | tail -1 | grep -oE '^[0-9]+' || echo NA)"
fi
# inner intermediate sizes (in the inner outputDir under $sbout); NA if absent.
ilines() { local f; f="$(find "$sbout" -name "$1" 2>/dev/null | head -1)"; [[ -n "$f" && -f "$f" ]] && { wc -l < "$f" | tr -d ' '; } || echo NA; }
pileup_lines="$(ilines pileup.txt)"; basecounts="$(ilines BaseCounts.pileup.txt)"
coverage="$(ilines Coverage.txt)"; modjunc="$(ilines modifiedJunc.txt)"; largedel="$(ilines large_deletions.txt)"
{
    echo "----- Splice-Break2 under-the-hood diagnostics ($SAMPLE) -----"
    echo "driver clean exit    : $([[ "$sb_rc" == 0 ]] && echo 1 || echo "0 (rc=$sb_rc)")"
    echo "MapSplice junctions  : $n_junc   ($(basename "${junc:-none}"); 0 => alignment/MapSplice produced none)"
    echo "  del4977-like junc  : $n_del4977j   (junction spanning ~8470..13447)"
    echo "LAST inner stage     : ${last_stage:-none}   (where Splice-Break2_0725.sh stopped)"
    echo "inner steps_failed   : $steps_failed"
    echo "intermediates        : pileup=$pileup_lines BaseCounts=$basecounts Coverage=$coverage modifiedJunc=$modjunc large_deletions=$largedel"
    echo "result file          : ${res_bytes} bytes; LargeMTDeletions rows: $n_del"
} >&2
while IFS= read -r lf; do
    [[ -f "$lf" ]] || continue
    echo "----- ${lf#"$out_abs"/} (tail 60) -----" >&2; tail -n 60 "$lf" >&2
done < <(find "$sblog" "$sbout" -type f \( -name '*.log' -o -name '*_nohup.log' \) 2>/dev/null | head -20)
printf 'caller\tnative_rc\tjunctions\tdel4977_junc\tlast_stage\tsteps_failed\tpileup\tbasecounts\tcoverage\tmodjunc\tlargedel\tresult_bytes\tcalls\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    splicebreak2 "$sb_rc" "$n_junc" "$n_del4977j" "${last_stage:-none}" "$steps_failed" "$pileup_lines" "$basecounts" "$coverage" "$modjunc" "$largedel" "$res_bytes" "$n_del" > "$out_abs/splicebreak2.runstatus"
# ---------------------------------------------------------------------------

ok=0
if [[ -n "$res" ]]; then
    cp "$res" "$out_abs/${SAMPLE}_LargeMTDeletions_WGS-only_NoPositionFilter.txt"
    log "calls -> ${SAMPLE}_LargeMTDeletions_WGS-only_NoPositionFilter.txt (${n_del} deletion rows)"; ok=1
else
    log "WARNING: no LargeMTDeletions output produced (see diagnostics above)"
fi
# Free the bulky per-sample install copy unless asked to keep it.
[[ "${SB_KEEP_INSTALL:-0}" == "1" ]] || rm -rf "$sbwork"
[[ "$ok" == 1 ]] || exit 1
log "done"
