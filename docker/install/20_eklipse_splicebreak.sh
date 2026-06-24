#!/usr/bin/env bash
###############################################################################
# 20_eklipse_splicebreak.sh — Python 2 / Java 8 toolbox
#
# eKLIPse (Python 2.7, BLAST+ >=2.3, circos) and Splice-Break2 (bash driver
# that shells out to /usr/bin/python [must be py2] and Java 8, with MapSplice2
# and bbmap bundled as ELF binaries) share one Python-2 conda environment.
###############################################################################
set -euxo pipefail

# Pin the last Python-2.7-compatible biopython/numpy so the py2.7 solve is
# deterministic (modern biopython has no py2.7 build and would force a fragile
# legacy cascade). blast/circos/samtools/openjdk are python-independent.
micromamba create -y -n py2tools \
    python=2.7 \
    biopython=1.76 \
    numpy=1.16 \
    tqdm \
    'blast>=2.9' \
    circos \
    'samtools>=1.9,<1.16' \
    openjdk=8

# --- eKLIPse ---------------------------------------------------------------
git clone https://github.com/dooguypapua/eKLIPse /opt/eKLIPse
git -C /opt/eKLIPse checkout "${EKLIPSE_SHA}"
# The bundled rCRS GenBank annotation we rely on:
test -f /opt/eKLIPse/data/NC_012920.1.gb
# eKLIPse imports local modules (pybam/spinner/tabulate) -> run from its dir.
# Assert the Python deps AND the external binaries eKLIPse shells out to are
# present, so a missing tool fails the build rather than the first analysis run.
micromamba run -n py2tools python -c "from Bio import SeqIO; import tqdm; print('biopython ok')"
micromamba run -n py2tools bash -c 'command -v blastn makeblastdb circos samtools'

# --- Splice-Break2 ---------------------------------------------------------
git clone https://github.com/brookehjelm/Splice-Break2 /tmp/Splice-Break2-src
git -C /tmp/Splice-Break2-src checkout "${SPLICEBREAK_SHA}"
mkdir -p /opt/Splice-Break2
tar -xzf /tmp/Splice-Break2-src/Paired-End_Download/Splice-Break2-v3.0.2_PAIRED-END.tar.gz \
    -C /opt/Splice-Break2
rm -rf /tmp/Splice-Break2-src
# Verify the extracted layout the wrapper depends on.
test -f /opt/Splice-Break2/Splice-Break2-v3.0.2_PAIRED-END/Splice-Break2_paired-end.sh
test -f /opt/Splice-Break2/Splice-Break2-v3.0.2_PAIRED-END/NC_012920.1/NC.fa
chmod -R u+rwX /opt/Splice-Break2

micromamba run -n py2tools java -version 2>&1 | head -1
micromamba clean -a -y
