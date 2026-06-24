#!/usr/bin/env bash
###############################################################################
# preprocess.sh — normalise one CRAM/BAM into the inputs every caller needs.
#
# Output (in --outdir):
#   <sample>.chrM.bam(.bai)   mito reads realigned to rCRS (contig "chrM",
#                             16569 bp); consumed by eKLIPse / MitoMut / MitoSeek
#   <sample>.R1.fastq.gz      mito read 1   } consumed by MitoSAlt / Splice-Break2
#   <sample>.R2.fastq.gz      mito read 2   } (they realign from FASTQ themselves)
#
# Why realign instead of slice-and-keep? Inputs vary (GRCh38 chrM, b37 MT, hg19
# 16571 Yoruba, mito-only subsets, different contig names). Re-aligning the mito
# reads to a single canonical rCRS gives every caller identical, predictable
# input regardless of how the source CRAM was produced.
###############################################################################
set -euo pipefail

SELF="preprocess"
log() { printf '[%s %s] %s\n' "$SELF" "$(date +%H:%M:%S)" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

REF="${MITO_SV_REF:-/opt/assets/rCRS.chrM.fa}"
THREADS=4
REFERENCE=""          # optional explicit CRAM reference
MT_CONTIG=""          # optional override of mito contig name

usage() {
    cat >&2 <<EOF
Usage: preprocess.sh --input <cram|bam> --sample <name> --outdir <dir>
                     [--reference <fasta>] [--threads N] [--mt-contig NAME]
EOF
    exit 2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)     INPUT="$2"; shift 2;;
        --sample)    SAMPLE="$2"; shift 2;;
        --outdir)    OUTDIR="$2"; shift 2;;
        --reference) REFERENCE="$2"; shift 2;;
        --threads)   THREADS="$2"; shift 2;;
        --mt-contig) MT_CONTIG="$2"; shift 2;;
        -h|--help)   usage;;
        *) die "unknown arg: $1";;
    esac
done

: "${INPUT:?--input required}"; : "${SAMPLE:?--sample required}"; : "${OUTDIR:?--outdir required}"
[[ -f "$INPUT" ]] || die "input not found: $INPUT"
mkdir -p "$OUTDIR"

# Run everything in the shared "mitosv" env (samtools/bwa).
run() { micromamba run -n mitosv "$@"; }

###############################################################################
# CRAM reference resolution
###############################################################################
declare -a REFARGS=()
is_cram=0
case "$INPUT" in
    *.cram|*.CRAM) is_cram=1;;
esac

if [[ "$is_cram" == 1 ]]; then
    if [[ -n "$REFERENCE" ]]; then
        [[ -f "$REFERENCE" ]] || die "reference not found: $REFERENCE"
        REFARGS=(-T "$REFERENCE")
        log "CRAM reference: $REFERENCE"
    else
        # No explicit reference: use a writable REF_CACHE, seeded with rCRS so
        # any rCRS-based CRAM decodes offline; fall back to the EBI ENA MD5
        # service for non-rCRS contigs (requires network).
        export REF_CACHE="${REF_CACHE:-$OUTDIR/refcache}"
        export REF_PATH="${REF_PATH:-${REF_CACHE}/%2s/%2s/%s:https://www.ebi.ac.uk/ena/cram/md5/%s}"
        mkdir -p "$REF_CACHE"
        seed_cache="$(run bash -lc 'ls "$CONDA_PREFIX"/share/samtools/seq_cache_populate.pl 2>/dev/null | head -1' || true)"
        if [[ -n "$seed_cache" ]]; then
            run perl "$seed_cache" -root "$REF_CACHE" "$REF" >/dev/null 2>&1 || \
                log "warn: could not seed rCRS into REF_CACHE (will rely on REF_PATH/EBI)"
        fi
        log "CRAM reference: REF_CACHE=$REF_CACHE (seeded rCRS), REF_PATH fallback=EBI ENA"
    fi
fi

###############################################################################
# 1. detect mito contig
###############################################################################
header="$(run samtools view -H ${REFARGS[@]+"${REFARGS[@]}"} "$INPUT")"
if [[ -n "$MT_CONTIG" ]]; then
    mt="$MT_CONTIG"
else
    mt="$(awk -F'\t' '
        /^@SQ/ {
            name=""; for (i=1;i<=NF;i++) if ($i ~ /^SN:/) name=substr($i,4)
            if (name ~ /^(chrM|chrMT|MT|M|NC_012920\.1|NC_012920)$/) { print name; exit }
        }' <<<"$header")"
fi
[[ -n "$mt" ]] || die "could not find a mitochondrial contig (chrM/MT/M) in $INPUT header"
log "mito contig: $mt"

###############################################################################
# 2. extract mito reads (use existing index if present, else build one in a
#    writable temp via -X so we never need to write next to the input)
###############################################################################
raw="$OUTDIR/${SAMPLE}.mito.raw.bam"
have_index=0
for ext in .crai .bai .csi; do
    [[ -e "${INPUT}${ext}" || -e "${INPUT%.*}${ext}" ]] && have_index=1
done

if [[ "$have_index" == 1 ]]; then
    log "extracting mito reads via existing index"
    run samtools view -@ "$THREADS" -b ${REFARGS[@]+"${REFARGS[@]}"} "$INPUT" "$mt" > "$raw"
else
    log "no index found; building one in outdir"
    tmpidx="$OUTDIR/${SAMPLE}.input.idx"
    if run samtools index -@ "$THREADS" -o "$tmpidx" "$INPUT" 2>/dev/null; then
        run samtools view -@ "$THREADS" -b -X ${REFARGS[@]+"${REFARGS[@]}"} "$INPUT" "$tmpidx" "$mt" > "$raw"
        rm -f "$tmpidx"
    else
        log "indexing failed; full-scan filtering by contig (slower)"
        run samtools view -h ${REFARGS[@]+"${REFARGS[@]}"} "$INPUT" \
            | awk -v m="$mt" 'BEGIN{OFS="\t"} /^@/ || $3==m' \
            | run samtools view -@ "$THREADS" -b - > "$raw"
    fi
fi

nreads="$(run samtools view -c "$raw" || echo 0)"
log "mito reads extracted: $nreads"
[[ "$nreads" -gt 0 ]] || die "no reads on contig '$mt' — wrong contig or empty input"

###############################################################################
# 3. mito reads -> FASTQ (name-collated)
###############################################################################
r1="$OUTDIR/${SAMPLE}.R1.fastq.gz"
r2="$OUTDIR/${SAMPLE}.R2.fastq.gz"
sgl="$OUTDIR/${SAMPLE}.singletons.fastq.gz"   # mate-lost (-s)
othr="$OUTDIR/${SAMPLE}.other.fastq.gz"        # READ_OTHER, i.e. flag-0 single-end (-0)
# -0 captures unpaired/single-end reads instead of discarding them (sending -0 to
# /dev/null would drop every read of a genuinely single-end mito input).
run bash -c "samtools collate -u -@ $THREADS '$raw' -O \
    | samtools fastq -@ $THREADS -1 '$r1' -2 '$r2' -0 '$othr' -s '$sgl' -n -"

count_reads() { local f="$1"; echo $(( $(run bash -c "zcat -f '$f' 2>/dev/null | wc -l" || echo 0) / 4 )); }
r1n="$(count_reads "$r1")"
# Unpaired reads = singletons (-s) + READ_OTHER (-0); align them single-end.
se="$OUTDIR/${SAMPLE}.se.fastq.gz"
cat "$sgl" "$othr" > "$se" 2>/dev/null || true   # concatenated gzip members are valid gzip
sen="$(count_reads "$se")"
log "FASTQ: $r1n read pairs, $sen unpaired reads"

###############################################################################
# 4. realign to canonical rCRS (contig chrM, 16569)
###############################################################################
[[ -f "${REF}.bwt" ]] || die "bwa index missing for $REF (expected ${REF}.bwt); build it with 'bwa index $REF'"
norm="$OUTDIR/${SAMPLE}.chrM.bam"
rg="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}"
tmp_p="$OUTDIR/${SAMPLE}.paired.bam"
tmp_s="$OUTDIR/${SAMPLE}.se.bam"
declare -a merge=()

if [[ "$r1n" -gt 0 ]]; then
    run bash -c "bwa mem -t $THREADS -R '$rg' '$REF' '$r1' '$r2' \
        | samtools sort -@ $THREADS -o '$tmp_p' -"
    merge+=("$tmp_p")
fi
if [[ "$sen" -gt 0 ]]; then
    run bash -c "bwa mem -t $THREADS -R '$rg' '$REF' '$se' \
        | samtools sort -@ $THREADS -o '$tmp_s' -"
    merge+=("$tmp_s")
fi
[[ "${#merge[@]}" -gt 0 ]] || die "no reads survived FASTQ conversion"

if [[ "${#merge[@]}" -eq 1 ]]; then
    mv "${merge[0]}" "$norm"
else
    run samtools merge -f -@ "$THREADS" "$norm" "${merge[@]}"
    rm -f "${merge[@]}"
fi
run samtools index -@ "$THREADS" "$norm"     # -> <sample>.chrM.bam.bai
rm -f "$raw" "$sgl" "$othr" "$se"            # keep R1/R2 (feed MitoSAlt/Splice-Break2)

log "done: $norm + ${SAMPLE}.R1/R2.fastq.gz"
echo "$norm"
