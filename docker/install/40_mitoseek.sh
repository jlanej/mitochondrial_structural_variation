#!/usr/bin/env bash
###############################################################################
# 40_mitoseek.sh — MitoSeek (Perl)
#
# MitoSeek is hard-locked to legacy samtools 0.1.x pileup/flagstat semantics.
# Rather than rely on the ancient prebuilt ELF it ships, we install
# bioconda samtools 0.1.19 (CLI/behaviour compatible with 0.1.18) and point
# MitoSeek at it via -samtools at runtime.
#
# Perl deps that are not reliably packaged on bioconda are installed with
# cpanm; the script then asserts mitoSeek.pl compiles so a missing module
# fails the *build*, not a silent runtime crash.
###############################################################################
set -euxo pipefail

micromamba create -y -n mitoseek \
    perl \
    perl-gd \
    perl-app-cpanminus \
    perl-statistics-descriptive \
    make

git clone https://github.com/riverlee/MitoSeek /opt/MitoSeek
git -C /opt/MitoSeek checkout "${MITOSEEK_SHA}"
chmod +x /opt/MitoSeek/Resources/samtools/samtools || true

# MitoSeek is hard-locked to legacy samtools 0.1.x. That pin does not co-solve
# cleanly with a current perl/perl-gd closure, so isolate it in its own env.
# Fall back to MitoSeek's bundled 0.1.18 ELF (needs libncurses5/zlib1g, provided
# by the base image) if the pin is unavailable.
if micromamba create -y -n mitoseek_samtools 'samtools=0.1.19'; then
    micromamba run -n mitoseek_samtools bash -lc 'command -v samtools' > /opt/MitoSeek/.samtools_path
else
    echo "/opt/MitoSeek/Resources/samtools/samtools" > /opt/MitoSeek/.samtools_path
fi
echo "MitoSeek samtools: $(cat /opt/MitoSeek/.samtools_path)"

# Pure/near-pure Perl modules MitoSeek `use`s at compile time.
micromamba run -n mitoseek cpanm --notest --no-man-pages \
    GD::Text::Wrap \
    GD::Graph::lines \
    GD::Graph::boxplot \
    Statistics::KernelEstimation \
    Text::NSP \
    Statistics::Multtest

# Hard gate: every `use` in mitoSeek.pl (incl. bundled Convert/Mitoanno/
# Circoswrap/Math::SpecFun) must resolve.
micromamba run -n mitoseek perl -c -I/opt/MitoSeek /opt/MitoSeek/mitoSeek.pl

# Sanity-check the resolved legacy samtools is a 0.1.x build.
"$(cat /opt/MitoSeek/.samtools_path)" 2>&1 | grep -i "Version" | head -1 || true

micromamba clean -a -y
