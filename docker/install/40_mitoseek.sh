#!/usr/bin/env bash
###############################################################################
# 40_mitoseek.sh — MitoSeek (Perl)
#
# We run MitoSeek with -noch -nocs -noQC, so its plotting subs (_boxplot /
# _density / _histogram) are never called — they live only inside `if($qc){...}`
# blocks. The fragile part of the install is the GD::Graph::* / GD::Text /
# Statistics::KernelEstimation stack (especially GD::Graph::boxplot, a
# poorly-maintained CPAN dist). Those are used purely as runtime method calls
# (`new GD::Graph::boxplot(...)`, `->new()`), so we eval-guard their `use` lines
# and never install them. We KEEP `use GD;` (from conda perl-gd) because the
# plotting subs reference GD's exported font barewords (gdGiantFont,
# gdMediumBoldFont) which must be defined at COMPILE time under `use strict`.
# Heteroplasmy genuinely needs Text::NSP (Fisher) + Statistics::Multtest (BH);
# Math::SpecFun::Beta and Convert/Mitoanno/Circoswrap ship in the repo.
#
# MitoSeek is hard-locked to legacy samtools 0.1.x semantics; we provide it in
# an isolated env (so the old pin never co-solves with perl) and fall back to
# MitoSeek's bundled 0.1.18 ELF if that pin is unavailable.
###############################################################################
set -euxo pipefail

micromamba create -y -n mitoseek \
    perl \
    perl-gd \
    perl-app-cpanminus \
    make

git clone https://github.com/riverlee/MitoSeek /opt/MitoSeek
git -C /opt/MitoSeek checkout "${MITOSEEK_SHA}"
chmod +x /opt/MitoSeek/Resources/samtools/samtools || true

# Make the plotting-only modules optional (their subs never run under -noQC).
# NB: keep `use GD;` (line 11) so the gd*Font barewords stay defined at compile.
perl -0pi -e '
  s/^use GD::Text::Wrap;/eval { require GD::Text::Wrap; };/m;
  s/^use GD::Graph::boxplot;/eval { require GD::Graph::boxplot; };/m;
  s/^use Statistics::KernelEstimation;/eval { require Statistics::KernelEstimation; };/m;
  s/^use GD::Graph::lines;/eval { require GD::Graph::lines; };/m;
' /opt/MitoSeek/mitoSeek.pl

# Legacy samtools 0.1.x in its own env; bundled ELF fallback. Both the create
# and the path capture are in the condition so neither can abort the build.
if micromamba create -y -n mitoseek_samtools 'samtools=0.1.19' \
   && micromamba run -n mitoseek_samtools bash -lc 'command -v samtools' > /opt/MitoSeek/.samtools_path; then
    echo "MitoSeek samtools (isolated conda env): $(cat /opt/MitoSeek/.samtools_path)"
else
    echo "/opt/MitoSeek/Resources/samtools/samtools" > /opt/MitoSeek/.samtools_path
    echo "MitoSeek samtools (bundled 0.1.18 ELF): $(cat /opt/MitoSeek/.samtools_path)"
fi

# The two pure-Perl modules the heteroplasmy path genuinely needs.
micromamba run -n mitoseek cpanm --notest --no-man-pages \
    Text::NSP \
    Statistics::Multtest

# Hard gate: every remaining `use` in mitoSeek.pl (Text::NSP, Statistics::Multtest,
# bundled Math::SpecFun::Beta + Convert/Mitoanno/Circoswrap) must resolve.
micromamba run -n mitoseek perl -c -I/opt/MitoSeek /opt/MitoSeek/mitoSeek.pl

# Sanity-check the resolved legacy samtools is a 0.1.x build (non-fatal).
"$(cat /opt/MitoSeek/.samtools_path)" 2>&1 | grep -i "Version" | head -1 || true

micromamba clean -a -y
