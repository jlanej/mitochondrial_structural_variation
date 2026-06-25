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

# MITOHPC_REF tracks the sv-calling branch (CI passes the resolved commit SHA so
# changes to our MitoHPC repo land on the next build); falls back to the branch
# name for local builds. Record the resolved commit so the benchmarked version is
# always discoverable at /opt/MitoHPC/GIT_SHA.
MITOHPC_REF="${MITOHPC_REF:-sv-calling}"
git clone https://github.com/jlanej/MitoHPC /opt/MitoHPC-src
git -C /opt/MitoHPC-src checkout "${MITOHPC_REF}"
# Keep only what the SV caller needs (scripts/ + RefSeq/) to stay lean.
mkdir -p /opt/MitoHPC
cp -r /opt/MitoHPC-src/scripts /opt/MitoHPC/scripts
cp -r /opt/MitoHPC-src/RefSeq  /opt/MitoHPC/RefSeq
git -C /opt/MitoHPC-src rev-parse HEAD > /opt/MitoHPC/GIT_SHA
echo "MitoHPC checked out $(cat /opt/MitoHPC/GIT_SHA) (ref=${MITOHPC_REF})"
rm -rf /opt/MitoHPC-src

test -f /opt/MitoHPC/scripts/callSV.sh
test -f /opt/MitoHPC/scripts/callsv.py
test -s /opt/MitoHPC/RefSeq/chrM.fa
test -s /opt/MitoHPC/scripts/sv.vcf            # VCF header template

# Caller runs under the mitosv env (python3 + pysam); confirm it imports + parses.
micromamba run -n mitosv python -c "import pysam"
micromamba run -n mitosv python /opt/MitoHPC/scripts/callsv.py --help >/dev/null
echo "MitoHPC SV caller installed at /opt/MitoHPC (scripts/ + RefSeq/)"
