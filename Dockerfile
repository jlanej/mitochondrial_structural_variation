# syntax=docker/dockerfile:1
###############################################################################
# Mitochondrial structural-variant caller suite
#
# A single, batteries-included image bundling five mtDNA SV/deletion callers:
#   - eKLIPse        (Python 2, BLAST+, circos)        -> conda env "py2tools"
#   - Splice-Break2  (bash + Python 2 + Java 8, MapSplice2 bundled)
#   - MitoMut        (Python 3, pysam, UCSC BLAT)      -> conda env "mitosv"
#   - MitoSeek       (Perl, bundled samtools 0.1.18)   -> conda env "mitoseek"
#   - MitoSAlt       (Perl + R, HISAT2/LAST)           -> conda env "mitosalt"
#
# Because the callers have hard, mutually incompatible dependency sets
# (Python 2 vs 3, samtools 0.1.18 vs 1.x, etc.) each lives in its own
# micromamba environment. A modern "mitosv" env hosts the shared pipeline
# (preprocessing, BAM<->FASTQ, post-processing).
###############################################################################
FROM mambaorg/micromamba:1.5.8

USER root
ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# --- system build/runtime deps (glibc base required by bundled ELF binaries) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
        git wget curl unzip ca-certificates procps bash gawk bc \
        bsdextrautils \
        build-essential make gcc g++ \
        libgd-dev zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libdb-dev \
        libncurses5 libgomp1 \
    && rm -rf /var/lib/apt/lists/*
# bsdextrautils provides `column`, which Splice-Break2's inner script uses to
# format its deletion tables (column -t at lines 110/142/170). Without it those
# pipes write empty files -> CompareMT throws IndexOutOfBounds -> the result is
# clobbered to 0 bytes on every sample (confirmed via the committed inner log).

# Best-effort libtinfo.so.5 for MitoSeek's bundled 0.1.18 samtools ELF fallback
# (the primary path uses conda samtools 0.1.19). Absent on some Debian releases,
# so tolerate failure.
RUN apt-get update \
    && (apt-get install -y --no-install-recommends libtinfo5 || true) \
    && rm -rf /var/lib/apt/lists/*

# Channel config (conda-forge first, then bioconda, per bioconda guidance)
RUN micromamba config append channels conda-forge \
    && micromamba config append channels bioconda \
    && micromamba config set channel_priority flexible

# --- copy assets, vendored MitoSAlt, install scripts, and pipeline code ---
COPY assets/                /opt/assets/
COPY vendor/MitoSAlt_1.1.1/ /opt/MitoSAlt/
COPY docker/install/        /opt/install/
COPY pipeline/              /opt/pipeline/

# Pin THIRD-PARTY caller revisions for reproducibility (we don't control these).
ENV EKLIPSE_SHA=3606cb2edac983d2623ddc667b49206c3d01373c \
    SPLICEBREAK_SHA=7b4ee7aed77586e67dc9fa2710288317e133f7cf \
    MITOMUT_SHA=ba56a65a5fc5728b2807d1253f6233db56e1c391 \
    MITOSEEK_SHA=624efc623832e3ca7f1095460ce4bc4e68bf8503

# --- shared pipeline env + reference indexes (samtools/bwa/minimap2/pysam) ---
RUN bash /opt/install/10_mitosv.sh

# --- eKLIPse + Splice-Break2 (Python 2 / Java 8 toolbox) ---
RUN bash /opt/install/20_eklipse_splicebreak.sh

# --- MitoMut (runs inside the mitosv env; only needs a clone + BLAT) ---
RUN bash /opt/install/30_mitomut.sh

# --- MitoSeek (Perl + GD + bundled samtools 0.1.18) ---
RUN bash /opt/install/40_mitoseek.sh

# --- MitoSAlt (Perl + R + LAST, rCRS index built from bundled reference) ---
RUN bash /opt/install/50_mitosalt.sh

# --- MitoHPC reference SV caller (Python3 + pysam; runs in the mitosv env) ---
# MitoHPC is OUR repo, so the image TRACKS its sv-calling branch rather than
# pinning: CI resolves the latest sv-calling commit and passes it as MITOHPC_REF,
# so changes to MitoHPC are reflected on the next build (and the exact commit is
# recorded at /opt/MitoHPC/GIT_SHA). Placed here, just before its install step, so
# a new ref only busts this layer — the expensive caller envs above stay cached.
ARG MITOHPC_REF=sv-calling
ENV MITOHPC_REF=${MITOHPC_REF}
RUN bash /opt/install/60_mitohpc.sh

# Splice-Break2 hardcodes /usr/bin/python and expects it to be Python 2.
RUN ln -sf /opt/conda/envs/py2tools/bin/python2.7 /usr/bin/python

# Make the wrapper scripts executable and expose them on PATH.
RUN chmod -R +x /opt/pipeline && \
    ln -s /opt/pipeline/run_sample.sh   /usr/local/bin/mito-sv-run-sample && \
    ln -s /opt/pipeline/postprocess.py  /usr/local/bin/mito-sv-postprocess

ENV MITO_SV_HOME=/opt/pipeline \
    MITO_SV_ASSETS=/opt/assets \
    MITO_SV_REF=/opt/assets/rCRS.chrM.fa

WORKDIR /work
ENTRYPOINT ["/opt/pipeline/run_sample.sh"]
CMD ["--help"]
