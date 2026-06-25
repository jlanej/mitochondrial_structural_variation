#!/usr/bin/env bash
###############################################################################
# 10_mitosv.sh — shared pipeline environment + reference indexes
#
# The "mitosv" env hosts everything the pipeline itself needs: CRAM/BAM I/O,
# realignment, FASTQ extraction, post-processing, AND the MitoMut caller
# (Python 3 + pysam + UCSC BLAT).
###############################################################################
set -euxo pipefail

micromamba create -y -n mitosv \
    python=3.10 \
    'samtools>=1.19' \
    'bcftools>=1.19' \
    htslib \
    bwa \
    minimap2 \
    pysam \
    numpy \
    pandas \
    ucsc-blat \
    perl
# perl: the LOD generator's circular-aware realign pipes through MitoHPC's
# circSam.pl (Perl); minimap2 + samtools above complete that path.

# Build the canonical rCRS (NC_012920.1, contig "chrM", 16569 bp) indexes used
# by the preprocessing realignment step. Bundled in the image -> no runtime DL.
micromamba run -n mitosv samtools faidx /opt/assets/rCRS.chrM.fa
micromamba run -n mitosv bwa index /opt/assets/rCRS.chrM.fa

# Sanity check the toolchain.
micromamba run -n mitosv samtools --version | head -1
micromamba run -n mitosv bwa 2>&1 | grep -i version || true
micromamba run -n mitosv blat 2>&1 | head -1 || true
micromamba run -n mitosv python -c "import pysam, numpy, pandas; print('pysam', pysam.__version__)"

micromamba clean -a -y
