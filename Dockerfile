# ==============================================================
# K-CHOPORE Dockerfile
# Comprehensive ONT Direct RNA Sequencing Analysis Pipeline
# ==============================================================
# This Dockerfile builds the complete K-CHOPORE environment with
# all tools for basecalling, QC, alignment, isoform analysis,
# epitranscriptomic modification detection, and differential expression.
#
# Base image: Ubuntu 22.04
# Maintained by: Pelayo Gonzalez de Lena Rodriguez
# ==============================================================

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Ensure 'python' command exists (many tools use #!/usr/bin/env python)
RUN ln -sf /usr/bin/python3 /usr/bin/python

LABEL maintainer="pelayovic"
LABEL description="K-CHOPORE: ONT Direct RNA-seq Analysis Pipeline"

# ==============================================================
# 1. System dependencies and essential tools
# ==============================================================
RUN apt-get update && apt-get install -y \
    curl wget git \
    jq \
    build-essential \
    gcc g++ \
    libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev \
    liblzma-dev libssl-dev libffi-dev libcurl4-openssl-dev \
    python3 python3-pip python3-venv python3-dev \
    default-jdk default-jre \
    r-base \
    bedtools \
    autoconf automake \
    libhdf5-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Install yq for YAML parsing
RUN curl -L https://github.com/mikefarah/yq/releases/download/v4.9.8/yq_linux_amd64 \
    -o /usr/bin/yq && chmod +x /usr/bin/yq

# ==============================================================
# 2. Samtools (alignment stats, sort, index, flagstat)
# ==============================================================
RUN wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar -xjf samtools-1.19.tar.bz2 && \
    cd samtools-1.19 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf samtools-1.19*

# ==============================================================
# 3. Minimap2 (splice-aware long-read aligner)
# ==============================================================
# NOTE: installed via apt, which is NOT version-pinnable here (Ubuntu 22.04 ships
# ~2.24). The analysis server runs minimap2 2.30 (envs/frozen/kchopore.txt).
# TODO (Pelayo, P2, optional): for full determinism, replace this apt install with
# a fixed-version minimap2 2.30 release download to match the server.
RUN apt-get update && apt-get install -y minimap2 && rm -rf /var/lib/apt/lists/*

# ==============================================================
# 4. Picard (BAM validation and manipulation)
# ==============================================================
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.7/picard.jar \
    -P /usr/local/bin/

# ==============================================================
# 5. Dorado (modern ONT basecaller for direct RNA)
# ==============================================================
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.0-linux-x64.tar.gz -P /opt/ && \
    tar -xzvf /opt/dorado-0.8.0-linux-x64.tar.gz -C /opt/ && \
    rm /opt/dorado-0.8.0-linux-x64.tar.gz && \
    ln -s /opt/dorado-0.8.0-linux-x64/bin/dorado /usr/local/bin/dorado
ENV PATH="/opt/dorado-0.8.0-linux-x64/bin:$PATH"

# ==============================================================
# 6. Guppy (legacy ONT basecaller, CPU version)
# ==============================================================
RUN rm -rf /opt/ont-guppy && \
    wget --no-check-certificate \
    https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.1.5_linux64.tar.gz \
    -O /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz && \
    tar -xvzf /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz -C /opt/ && \
    rm /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz
ENV PATH="/opt/ont-guppy-cpu/bin:$PATH"

# ==============================================================
# 7. StringTie2 (long-read isoform assembly)
# ==============================================================
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz && \
    tar -xzvf stringtie-2.2.1.Linux_x86_64.tar.gz && \
    cp stringtie-2.2.1.Linux_x86_64/stringtie /usr/local/bin/ && \
    rm -rf stringtie-2.2.1*

# ==============================================================
# 8. Nanopolish (signal-level analysis for m6Anet/xPore)
# ==============================================================
# Nanopolish must be compiled from source (no pip package available)
# Use serial build (-j1) to avoid race conditions between HDF5/eigen
# downloads and source compilation
# Pinned to tag v0.14.0 (last stable release) for a reproducible clone.
# TODO (Pelayo): the analysis server has no nanopolish binary in any conda env
# to read a true version from; v0.14.0 is the documented stable tag, not an
# observed install. Confirm against the binary actually used (nanopolish --version)
# before publishing, and adjust the tag if it differs.
RUN git clone --recursive https://github.com/jts/nanopolish.git /opt/nanopolish && \
    cd /opt/nanopolish && \
    git checkout v0.14.0 && git submodule update --init --recursive && \
    make -j1 && \
    ln -s /opt/nanopolish/nanopolish /usr/local/bin/nanopolish && \
    ln -s /opt/nanopolish/scripts/nanopolish_makerange.py /usr/local/bin/nanopolish_makerange.py

# VBZ compression plugin for HDF5 (required to read VBZ-compressed FAST5 files)
# The plugin is already bundled with Guppy at /opt/ont-guppy-cpu/lib/
# Point HDF5_PLUGIN_PATH there so nanopolish (and its bundled HDF5) can find it
ENV HDF5_PLUGIN_PATH="/opt/ont-guppy-cpu/lib"

# ==============================================================
# 9. Python package manager and Snakemake
# ==============================================================
# pulp and snakemake pinned to versions observed in the server `kchopore` env
# (envs/frozen/kchopore.txt: pulp 2.8.0, snakemake-minimal 9.13.4).
# The list_solvers->listSolvers sed targeted snakemake 5.x; snakemake 9.x already
# uses listSolvers, so it is now a guarded no-op (kept harmless via `|| true`).
# TODO (Pelayo): once a Docker build with snakemake 9.13.4 is validated, drop the
# sed line entirely.
RUN pip install --upgrade pip && \
    pip install pulp==2.8.0 && \
    pip install snakemake==9.13.4 && \
    sed -i 's/list_solvers/listSolvers/' /usr/local/lib/python3.10/dist-packages/snakemake/__init__.py || true

# ==============================================================
# 10. Bonito (alternative ONT basecaller)
# ==============================================================
# ont-bonito pinned to the observed 1.1.0 (envs/frozen/kchopore.txt).
RUN pip install ont-bonito==1.1.0 || true
# NOTE: Bonito models are large (~5GB). Download at runtime instead:
#   docker run -v bonito_models:/root/.local/share/bonito k-chopore bonito download --models

ENV PATH="/usr/local/bin/bonito:$PATH"

# ==============================================================
# 11. Python requirements (core pipeline)
# ==============================================================
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir setuptools wheel && \
    pip install --no-cache-dir -r /tmp/requirements.txt

# ==============================================================
# 12. NanoPlot, NanoComp, NanoFilt (ONT-specific QC tools)
# ==============================================================
# Pinned to versions observed in the working `kchopore` env (envs/frozen/kchopore.txt).
RUN pip install --no-cache-dir \
    NanoPlot==1.46.1 \
    NanoComp==1.25.6 \
    NanoFilt==2.8.0 \
    nanoget==1.19.3 \
    nanomath==1.4.0

# pycoQC installed with --no-deps so its legacy plotly==4.1.0 dependency pin is
# NOT enforced: in the working server env, pycoQC 2.5.2 runs fine against the
# plotly 6.x that NanoPlot needs (verified — NanoPlot 1.46.1 + plotly 6.3.1 +
# pycoQC 2.5.2 coexist in one env). setuptools must be reinstalled after NanoPlot
# (which may remove it) because pycoQC imports pkg_resources at runtime.
# pycoQC pinned to the observed 2.5.2.
RUN pip install --no-cache-dir setuptools && \
    pip install --no-cache-dir --no-deps pycoQC==2.5.2

# ==============================================================
# 13. m6Anet (m6A modification detection from nanopore signal)
# ==============================================================
# Server runs m6anet 2.1.0 in its own python-3.8 env (envs/frozen/m6anet.txt),
# where it pins torch==1.6.0. That torch is incompatible with this image's
# Python 3.10, so we keep the --no-deps install and supply compatible deps.
# m6anet is pinned to the observed 2.1.0; torch/scikit-learn are left unpinned
# because the image cannot use the env's torch 1.6.0 and no Python-3.10-compatible
# torch version was observed installed (pin from a validated image build).
RUN pip install --no-cache-dir --no-deps m6anet==2.1.0 && \
    pip install --no-cache-dir torch scikit-learn

# ==============================================================
# 14. xPore (differential RNA modification from nanopore)
# ==============================================================
# Pinned to xpore 2.1, the version installed in the server xpore env
# (envs/frozen/xpore.txt).
RUN pip install --no-cache-dir xpore==2.1 || \
    pip install --no-cache-dir --no-deps xpore==2.1

# ==============================================================
# 15. ELIGOS2 (RNA modification detection from error signatures)
# ==============================================================
# Pinned to the exact commit used by the analysis (server /home/usuario2/eligos2,
# HEAD = b205a5c, i.e. tag v2.1.0 + 8 commits: `git describe` -> v2.1.0-8-gb205a5c,
# committed 2024-01-29). eligos2 has no reliable release tag for this state, so we
# pin by SHA for a reproducible clone.
RUN git clone https://gitlab.com/piroonj/eligos2.git /home/eligos2 && \
    cd /home/eligos2 && \
    git checkout b205a5cec34cd2974e5ddae78f7fcf7beb49d9c8 && \
    pip install --no-cache-dir -r requirements.txt 2>/dev/null || true && \
    chmod +x /home/eligos2/Scripts/*.py 2>/dev/null || true

# ==============================================================
# 16. MultiQC (aggregate QC report)
# ==============================================================
RUN pip install --no-cache-dir multiqc==1.32

# ==============================================================
# 16. POD5 tools (modern ONT file format support)
# ==============================================================
RUN pip install --no-cache-dir pod5==0.3.39

# ==============================================================
# 17. Additional Python packages
# ==============================================================
# pysam pinned to the observed 2.x->0.23.3; tabulate already pinned.
# rpy2 and ont-fast5-api were NOT present in any captured server env, so no real
# version was observed — left unpinned. TODO (Pelayo): pin both from a validated
# image build (`pip show rpy2 ont-fast5-api`).
RUN pip install --no-cache-dir \
    pysam==0.23.3 \
    rpy2 \
    tabulate==0.9.0 \
    ont-fast5-api

# ==============================================================
# 18. R packages for differential expression
# ==============================================================
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")' && \
    Rscript -e 'BiocManager::install(c("DESeq2", "tximport", "apeglm"), ask=FALSE)' && \
    Rscript -e 'install.packages(c("samplesizeCMH", "ggplot2", "pheatmap", "EnhancedVolcano"), repos="https://cloud.r-project.org")' || true

# ==============================================================
# 19. FLAIR (isoform analysis for long reads)
# ==============================================================
RUN pip install --no-cache-dir flair-brookslab==2.0.0

# ==============================================================
# 20. Set up workspace
# ==============================================================
WORKDIR /workspace

# Copy project files (data/ is mounted at runtime via docker -v)
COPY scripts /workspace/scripts
COPY config /workspace/config
COPY requirements.txt /workspace/requirements.txt
COPY Snakefile /workspace/Snakefile

# Create data directory structure (will be overridden by volume mount)
RUN mkdir -p /workspace/data/raw/fastq /workspace/data/raw/fast5 \
    /workspace/data/raw/pod5 /workspace/data/raw/summaries \
    /workspace/data/reference/genome /workspace/data/reference/annotations \
    /workspace/data/reference/transcriptome

# ==============================================================
# 21. Fix setuptools for pycoQC (pkg_resources removed in v71+)
# ==============================================================
# Must be the LAST pip install to prevent later packages upgrading it
RUN pip install --no-cache-dir "setuptools<71"

# ==============================================================
# 22. Environment configuration
# ==============================================================
ENV PYTHONPATH="/usr/local/lib/python3.10/dist-packages"
ENV PATH="/home/eligos2:/home/eligos2/Scripts:$PATH"
ENV CUDA_HOME=/usr/local/cuda

# ==============================================================
# Verify installation
# ==============================================================
RUN echo "=== K-CHOPORE Build Verification ===" && \
    echo "Minimap2: $(minimap2 --version 2>&1 || echo 'not found')" && \
    echo "Samtools: $(samtools --version 2>&1 | head -1 || echo 'not found')" && \
    echo "Dorado: $(dorado --version 2>&1 || echo 'not found')" && \
    echo "StringTie: $(stringtie --version 2>&1 || echo 'not found')" && \
    echo "Nanopolish: $(nanopolish --version 2>&1 || echo 'not found')" && \
    echo "NanoPlot: $(NanoPlot --version 2>&1 || echo 'not found')" && \
    echo "NanoFilt: $(NanoFilt --version 2>&1 || echo 'not found')" && \
    echo "m6anet: $(pip show m6anet 2>&1 | grep Version || echo 'not found')" && \
    echo "MultiQC: $(multiqc --version 2>&1 || echo 'not found')" && \
    echo "Snakemake: $(snakemake --version 2>&1 || echo 'not found')" && \
    echo "R: $(R --version 2>&1 | head -1 || echo 'not found')" && \
    echo "=== K-CHOPORE Docker image built successfully! ==="

# === K-CHOPORE extra deps (analisis/figuras anadidos 2026) ===
# Pinned to versions observed in the server `viz` figure env (envs/frozen/viz.txt).
RUN pip install --no-cache-dir logomaker==0.8.7 upsetplot==0.9.0 seaborn==0.13.2 "kaleido==0.2.1"
RUN apt-get update && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev && rm -rf /var/lib/apt/lists/* && \
    Rscript -e 'install.packages(c("gprofiler2","RColorBrewer"), repos="https://cloud.r-project.org")'
