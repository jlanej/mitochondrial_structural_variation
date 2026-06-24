#!/usr/bin/env bash
###############################################################################
# 60_mitohpc.sh — MitoHPC SV caller (the reference method, for a fair head-to-head)
#
# MitoHPC's structural-variant caller (scripts/callSV.sh -> scripts/callsv.py) is
# pure Python3 + pysam — both already in the "mitosv" env — and reads the chrM
# BAM in-process (split-read SA:Z junctions + coverage drop). We vendor just its
# scripts/ + RefSeq/ (chrM reference + FP masks + gene annotation) and run it on
# the same normalised rCRS chrM BAM the other BAM-based callers consume, so all
# six methods get identical input.
###############################################################################
set -euxo pipefail

git clone https://github.com/jlanej/MitoHPC /opt/MitoHPC-src
git -C /opt/MitoHPC-src checkout "${MITOHPC_SHA}"
# Keep only what the SV caller needs (scripts/ + RefSeq/) to stay lean.
mkdir -p /opt/MitoHPC
cp -r /opt/MitoHPC-src/scripts /opt/MitoHPC/scripts
cp -r /opt/MitoHPC-src/RefSeq  /opt/MitoHPC/RefSeq
rm -rf /opt/MitoHPC-src

test -f /opt/MitoHPC/scripts/callSV.sh
test -f /opt/MitoHPC/scripts/callsv.py
test -s /opt/MitoHPC/RefSeq/chrM.fa
test -s /opt/MitoHPC/scripts/sv.vcf            # VCF header template

# Caller runs under the mitosv env (python3 + pysam); confirm it imports + parses.
micromamba run -n mitosv python -c "import pysam"
micromamba run -n mitosv python /opt/MitoHPC/scripts/callsv.py --help >/dev/null
echo "MitoHPC SV caller installed at /opt/MitoHPC (scripts/ + RefSeq/)"
