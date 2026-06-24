#!/usr/bin/env bash
###############################################################################
# run_sample.sh — end-to-end per-sample driver (runs INSIDE the container).
#
#   1. preprocess one CRAM/BAM -> normalised rCRS chrM BAM + mito FASTQ pair
#   2. run each requested caller on the appropriate input
#   3. record per-caller status
#
# Output layout (under --outdir, which should be a per-sample directory):
#   preprocess/<sample>.chrM.bam(.bai), <sample>.R1/R2.fastq.gz
#   eklipse/ mitosalt/ splicebreak2/ mitomut/ mitoseek/   (one per caller)
#   status.tsv
###############################################################################
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log() { printf '[run_sample %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

INPUT="" SAMPLE="" OUTDIR="" REFERENCE="" THREADS=4 MT_CONTIG=""
CALLERS="eklipse,mitosalt,splicebreak2,mitomut,mitoseek"
STRICT=0

usage() {
    cat >&2 <<EOF
Usage: run_sample.sh --input <cram|bam> --outdir <dir> [options]

  --input PATH        input CRAM or BAM (required)
  --outdir DIR        per-sample output directory (required)
  --sample NAME       sample name (default: input basename)
  --reference FASTA   reference for CRAM decoding (optional)
  --mt-contig NAME    override mito contig name detection (optional)
  --threads N         threads per step (default 4)
  --callers LIST      comma list or 'all' (default: all five)
  --strict            exit non-zero if any caller fails
EOF
    exit "${1:-2}"
}

while [[ $# -gt 0 ]]; do case "$1" in
    --input)     INPUT="$2"; shift 2;;
    --sample)    SAMPLE="$2"; shift 2;;
    --outdir)    OUTDIR="$2"; shift 2;;
    --reference) REFERENCE="$2"; shift 2;;
    --mt-contig) MT_CONTIG="$2"; shift 2;;
    --threads)   THREADS="$2"; shift 2;;
    --callers)   CALLERS="$2"; shift 2;;
    --strict)    STRICT=1; shift;;
    -h|--help)   usage 0;;
    *) log "unknown arg: $1"; usage 2;;
esac; done

[[ -n "$INPUT" && -n "$OUTDIR" ]] || usage 2
[[ -f "$INPUT" ]] || { log "ERROR: input not found: $INPUT"; exit 1; }
if [[ -z "$SAMPLE" ]]; then
    SAMPLE="$(basename "$INPUT")"; SAMPLE="${SAMPLE%.cram}"; SAMPLE="${SAMPLE%.bam}"
    SAMPLE="${SAMPLE%.chrM}"
fi
# Sanitise: callers use the sample name in filenames.
SAMPLE="$(printf '%s' "$SAMPLE" | tr ' /' '__')"
[[ "$CALLERS" == "all" ]] && CALLERS="eklipse,mitosalt,splicebreak2,mitomut,mitoseek"

mkdir -p "$OUTDIR"
OUTDIR="$(readlink -f "$OUTDIR")"
PRE="$OUTDIR/preprocess"
log "sample=$SAMPLE input=$INPUT outdir=$OUTDIR callers=$CALLERS threads=$THREADS"

###############################################################################
# 1. preprocess
###############################################################################
pp_args=(--input "$INPUT" --sample "$SAMPLE" --outdir "$PRE" --threads "$THREADS")
[[ -n "$REFERENCE" ]] && pp_args+=(--reference "$REFERENCE")
[[ -n "$MT_CONTIG" ]] && pp_args+=(--mt-contig "$MT_CONTIG")
if ! bash "$HERE/preprocess.sh" "${pp_args[@]}"; then
    log "ERROR: preprocessing failed for $SAMPLE"
    exit 1
fi
BAM="$PRE/${SAMPLE}.chrM.bam"
R1="$PRE/${SAMPLE}.R1.fastq.gz"
R2="$PRE/${SAMPLE}.R2.fastq.gz"
[[ -f "$BAM" && -f "$R1" && -f "$R2" ]] || { log "ERROR: preprocessing outputs missing"; exit 1; }

###############################################################################
# 2. run callers
###############################################################################
declare -A WRAP=(
    [eklipse]="$HERE/callers/run_eklipse.sh"
    [mitosalt]="$HERE/callers/run_mitosalt.sh"
    [splicebreak2]="$HERE/callers/run_splicebreak2.sh"
    [mitomut]="$HERE/callers/run_mitomut.sh"
    [mitoseek]="$HERE/callers/run_mitoseek.sh"
)

status_tsv="$OUTDIR/status.tsv"
printf 'caller\tstatus\tseconds\n' > "$status_tsv"
overall_fail=0

IFS=',' read -ra REQ <<< "$CALLERS"
for caller in "${REQ[@]}"; do
    caller="$(printf '%s' "$caller" | tr -d '[:space:]')"
    wrap="${WRAP[$caller]:-}"
    if [[ -z "$wrap" ]]; then log "skip unknown caller: $caller"; continue; fi
    cdir="$OUTDIR/$caller"; mkdir -p "$cdir"
    log "=== caller: $caller ==="
    t0=$SECONDS
    if bash "$wrap" --sample "$SAMPLE" --bam "$BAM" --r1 "$R1" --r2 "$R2" \
            --outdir "$cdir" --threads "$THREADS" > "$cdir/${caller}.log" 2>&1; then
        st=ok
    else
        st=failed; overall_fail=1
        log "caller $caller FAILED (see $cdir/${caller}.log)"
    fi
    printf '%s\t%s\t%s\n' "$caller" "$st" "$((SECONDS - t0))" >> "$status_tsv"
done

log "status:"; cat "$status_tsv" >&2

if [[ "$STRICT" == 1 && "$overall_fail" == 1 ]]; then
    log "ERROR: --strict and at least one caller failed"
    exit 1
fi
exit 0
