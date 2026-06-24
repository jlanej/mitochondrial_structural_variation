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
cp -a "$sbsrc" "$sbwork"
SB_PATH="$(readlink -f "$sbwork")"

indir="$out_abs/in"; sbout="$out_abs/out"; sblog="$out_abs/log"
mkdir -p "$indir" "$sbout" "$sblog"
# The driver expects <sample>.R1.fastq / <sample>.R2.fastq (auto-gunzips .gz).
cp "$r1_abs" "$indir/${SAMPLE}.R1.fastq.gz"
cp "$r2_abs" "$indir/${SAMPLE}.R2.fastq.gz"

log "running Splice-Break2 on $SAMPLE"
# Runs in py2tools env: provides Java 8 + Python 2 (and /usr/bin/python -> py2).
# bbmap (dedupe/bbduk/reformat) autodetects a JVM heap and sets -Xms == -Xmx to
# ~85% of free RAM; on memory-limited runners the JVM fails to RESERVE that heap
# at startup and aborts in ~2-3s (its stderr is swallowed into $sblog). Cap it
# via RQCMEM (MB) — bbmap's calcmem.sh uses RQCMEM to bound -Xmx/-Xms for every
# bbmap tool. Override with SB_MEM_MB for big cohorts on roomy HPC nodes.
SB_MEM_MB="${SB_MEM_MB:-2000}"
micromamba run -n py2tools bash -c "export RQCMEM='$SB_MEM_MB'; bash '$SB_PATH/Splice-Break2_paired-end.sh' \
    '$indir' '$sbout' '$sblog' '$SB_PATH' \
    --align=yes --ref=rCRS --fastq_keep=no --skip_preAlign=no" || log "driver returned non-zero"

res="$(find "$sbout" -name '*_LargeMTDeletions_WGS-only_NoPositionFilter.txt' | head -1 || true)"
ok=0
if [[ -n "$res" ]]; then
    cp "$res" "$out_abs/${SAMPLE}_LargeMTDeletions_WGS-only_NoPositionFilter.txt"
    log "calls -> ${SAMPLE}_LargeMTDeletions_WGS-only_NoPositionFilter.txt"; ok=1
else
    log "WARNING: no LargeMTDeletions output produced; Splice-Break2 driver logs follow:"
    # bbmap stderr goes to $sblog/<sample>.log; the MapSplice error log and the
    # driver's nohup log (which records '#ERROR: Mapsplice Fail') live inside the
    # per-sample output dir. Surface all of them so the real failure is captured
    # in this caller log (and thus the committed *.FAILED.log).
    while IFS= read -r lf; do
        [[ -f "$lf" ]] || continue
        log "----- ${lf#"$out_abs"/} (tail) -----"; tail -n 40 "$lf" >&2
    done < <(find "$sblog" "$sbout" -type f \( -name '*.log' -o -name '*_nohup.log' \) 2>/dev/null | head -20)
fi
# Free the bulky per-sample install copy unless asked to keep it.
[[ "${SB_KEEP_INSTALL:-0}" == "1" ]] || rm -rf "$sbwork"
[[ "$ok" == 1 ]] || exit 1
log "done"
