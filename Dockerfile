# Stage 1: Build stage
FROM python:3.10-slim-bookworm AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    bzip2 \
    libcurl4-openssl-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Set samtools version
ENV SAMTOOLS_VERSION=1.23.1

# Download, compile, and install samtools + htslib utilities
RUN curl -L https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -o samtools.tar.bz2 \
    && tar -xf samtools.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/opt/samtools --enable-libcurl \
    && make all all-htslib \
    && make install install-htslib

# Stage 2: Runtime stage
FROM python:3.10-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH="/opt/conda/bin:/opt/conda/envs/fifa-r/bin:/opt/samtools/bin:${PATH}"

# Install runtime dependencies for samtools, gcloud, and other utilities
RUN apt-get update && apt-get install -y \
    libcurl4 \
    libssl3 \
    libbz2-1.0 \
    liblzma5 \
    libtinfo6 \
    zlib1g \
    curl \
    wget \
    gnupg \
    ca-certificates \
    apt-transport-https \
    lsb-release \
    bzip2 \
    git \
    jq \
    build-essential \
    cmake \
    gfortran \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libssl-dev \
    libuv1-dev \
    libxml2-dev \
    pkg-config \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba for the pinned R runtime.
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# Install the R packages required by run_mobster.R. Keep the stack aligned with
# conda/FIFA_environment.yml (R 4.4) and avoid pinning Bioconductor package
# versions that were built against older R ABIs.
RUN set -eux; \
    for attempt in 1 2 3 4 5; do \
        micromamba create -y -n fifa-r -c conda-forge -c bioconda \
            r-base=4.4 \
            r-remotes \
            r-dplyr=1.1.4 \
            r-tidyr=1.3.1 \
            r-vroom=1.6.5 \
            bioconductor-genomicranges \
            bioconductor-variantannotation \
            r-sads=0.6.3 \
            r-crayon=1.5.3 \
            r-ggplot2=4.0.3 \
            r-ggthemes=5.1.0 \
            r-tidyverse=2.0.0 \
            r-ggpubr=0.6.0 \
            r-ggrepel=0.9.5 \
            r-cli=3.6.3 \
            r-magrittr=2.0.3 \
            r-cowplot=1.2.0 \
            r-reshape2=1.4.4 \
            r-clisymbols=1.2.0 \
            r-pio=0.1.0 \
            r-easypar=1.0.0 \
            r-ctree=1.1.0 \
            r-dndscv=0.1.0 \
        && break; \
        echo "micromamba create failed on attempt ${attempt}; retrying" >&2; \
        if [ "${attempt}" = "5" ]; then \
            exit 1; \
        fi; \
        micromamba clean --index-cache --yes || true; \
        sleep $((attempt * 10)); \
    done; \
    micromamba clean --all --yes

RUN /opt/conda/envs/fifa-r/bin/Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'https://cloud.r-project.org'); need <- c('GenomeInfoDbData', 'GenomicRanges', 'VariantAnnotation'); miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]; if (length(miss)) BiocManager::install(miss, ask = FALSE, update = FALSE)" \
    && /opt/conda/envs/fifa-r/bin/Rscript -e "remotes::install_github('caravagnalab/mobster@85c898f087b46e15f79144bf57bdc019678e2481', dependencies = FALSE, upgrade = 'never')" \
    && /opt/conda/envs/fifa-r/bin/Rscript -e "pkgs <- c('dplyr', 'tidyr', 'vroom', 'GenomeInfoDbData', 'GenomicRanges', 'VariantAnnotation', 'mobster'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop(paste('Missing required R packages:', paste(missing, collapse = ', ')))"

# Install Google Cloud SDK (provides `gcloud` CLI)
RUN set -eux; \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg; \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" > /etc/apt/sources.list.d/google-cloud-sdk.list; \
    apt-get update && apt-get install -y google-cloud-sdk && rm -rf /var/lib/apt/lists/*

# Install Python requirements in an isolated virtual environment.
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:/opt/conda/bin:/opt/conda/envs/fifa-r/bin:/opt/samtools/bin:${PATH}"
RUN /usr/local/bin/python3 -m venv "$VIRTUAL_ENV"
RUN python -m pip install --no-cache-dir --upgrade \
    "pip<25" \
    "setuptools<82" \
    wheel \
    pybind11 \
    "Cython<3" \
    "numpy==1.24.4" \
    google-auth \
    requests
RUN python -m pip install --no-cache-dir "setuptools<58"
RUN python -m pip install --no-cache-dir --no-build-isolation \
    "pysam==0.15.4" \
    && python -m pip install --no-cache-dir --upgrade "setuptools<82" \
    && python -c "import pkg_resources, pybind11, numpy, pysam"
RUN set -eux; \
    mkdir -p /tmp/build; \
    wget https://raw.githubusercontent.com/JenniferShelton/fifa/refs/heads/main/requirements.txt -O /tmp/build/requirements.txt; \
    python -m pip install --no-build-isolation --no-cache-dir -r /tmp/build/requirements.txt

# Copy compiled binaries and libraries from the builder stage
COPY --from=builder /opt/samtools /opt/samtools


# Install FIFA into an empty directory; cloning into the copied build context
# fails once that context contains files.
WORKDIR /opt/fifa
RUN git clone --depth 1 https://github.com/JenniferShelton/fifa.git .
RUN echo "alias fifa='python3 /opt/fifa/src/cli.py'" >> /etc/bash.bashrc