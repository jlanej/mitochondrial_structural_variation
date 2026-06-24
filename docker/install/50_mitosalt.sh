#!/usr/bin/env bash
###############################################################################
# 50_mitosalt.sh — MitoSAlt (Perl + R, LAST aligner)
#
# We run MitoSAlt in "enriched" mode (the input is already mitochondrially
# enriched after our preprocessing), which uses ONLY the LAST MT index and
# skips the giant nuclear-genome HISAT2 step (whose upstream download URLs are
# long dead). The rCRS LAST index is built here from the bundled reference.
#
# The MitoSAlt source is COPYed to /opt/MitoSAlt by the Dockerfile.
###############################################################################
set -euxo pipefail

micromamba create -y -n mitosalt \
    perl \
    last \
    bbmap \
    'samtools>=1.9' \
    bedtools \
    ucsc-bedgraphtobigwig \
    r-base \
    r-plotrix \
    r-rcolorbrewer \
    bioconductor-biostrings

cd /opt/MitoSAlt
mkdir -p genome bin bam bw tab indel log plot

# Build the rCRS (chrM, 16569 bp) reference + LAST index MitoSAlt aligns to.
cp /opt/assets/rCRS.chrM.fa genome/human_mt_rCRS.fasta
micromamba run -n mitosalt samtools faidx genome/human_mt_rCRS.fasta
micromamba run -n mitosalt lastdb -uNEAR genome/human_mt_rCRS genome/human_mt_rCRS.fasta

# Enriched-mode config. Tool paths are bare command names resolved on PATH
# when the wrapper runs `micromamba run -n mitosalt`. DB paths are relative to
# the per-sample working dir, which symlinks genome/ -> /opt/MitoSAlt/genome.
cat > /opt/MitoSAlt/config_pipeline.txt <<'CFG'
#TOOLS
# hisat2 + sambamba are unused in enriched mode (nu_mt=no, cn_mt=no) and are
# intentionally NOT installed in the mitosalt env. Do not flip cn_mt to "yes"
# without adding sambamba (and a nuclear genome) to the env.
hisat2 = hisat2
lastal = lastal
lastsp = last-split
mfcv = maf-convert
reformat = reformat.sh
samtools = samtools
sambamba = sambamba
b2fq = bamToFastq
gcov = genomeCoverageBed
intersectBed = intersectBed
sortBed = sortBed
clusterBed = clusterBed
randomBed = randomBed
groupBy = groupBy
bg2bw = bedGraphToBigWig

#DATABASES
hsindex = genome/human_mt_rCRS
faindex = genome/human_mt_rCRS.fasta.fai
lastindex = genome/human_mt_rCRS
mtfaindex = genome/human_mt_rCRS.fasta.fai
gsize = genome/human_mt_rCRS.fasta.fai
MT_fasta = genome/human_mt_rCRS.fasta

#COMPUTATION
threads = 4

#MITOCHONDRIA FEATURES
refchr = chrM
msize = 16569
exclude = 5
orihs = 16081
orihe = 407
orils = 5730
orile = 5763

#SCORING AND FILTERING FEATURES
score_threshold = 80
evalue_threshold = 0.00001
split_length = 15
paired_distance = 1000
deletion_threshold_min = 30
deletion_threshold_max = 30000
breakthreshold = -2
cluster_threshold = 5
breakspan = 15
sizelimit = 10000
hplimit = 0.01
flank = 15
split_distance_threshold = 5

#STEPS
dna = yes
enriched = yes
nu_mt = no
rmtmp = yes
o_mt = yes
i_del = yes
cn_mt = no
CFG

# Make sure the R plotting/scoring deps load from the conda lib (no runtime
# Biostrings network install, which delplot.R would otherwise attempt).
micromamba run -n mitosalt Rscript -e 'suppressMessages({library(plotrix);library(RColorBrewer);library(Biostrings)}); cat("MitoSAlt R deps OK\n")'

micromamba clean -a -y
