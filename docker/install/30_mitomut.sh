#!/usr/bin/env bash
###############################################################################
# 30_mitomut.sh — MitoMut (Python 2/3 compatible; we run it under py3)
#
# Needs only: a clone, pysam, samtools (>=1.9 for `samtools fasta`/`cat -o`)
# and UCSC BLAT — all already provided by the "mitosv" env. Ships its own
# rCRS reference (mt.fasta, header ">chrM").
###############################################################################
set -euxo pipefail

git clone https://github.com/shane-e945/MitoMut /opt/MitoMut
git -C /opt/MitoMut checkout "${MITOMUT_SHA}"

test -f /opt/MitoMut/MitoMut.py
test -f /opt/MitoMut/mt.fasta

# Confirm it imports under the shared env (catches gross dependency drift).
micromamba run -n mitosv python -c "import pysam"
micromamba run -n mitosv python /opt/MitoMut/MitoMut.py -h >/dev/null 2>&1 || true
echo "MitoMut installed at /opt/MitoMut"
