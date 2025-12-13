# syntax=docker/dockerfile:1.7
# =============================================================================
# VNtyper 2.0 - Multi-Stage Production Dockerfile
# =============================================================================
# Build Strategy:
# - Stage 1: Base conda environment setup
# - Stage 2: Build conda environments and compile dependencies
# - Stage 3: Minimal runtime image with only necessary components
#
# Expected outcome: 60-70% smaller image size, improved security
# =============================================================================

# =============================================================================
# Stage 1: Base Dependencies (Shared)
# =============================================================================
FROM condaforge/miniforge3:24.11.3-0@sha256:bcfc04b3562faa8dfd48f3767e2bd7b73a8529f9611b508d42e5bd227db1342d AS base

# OCI Labels for metadata
LABEL org.opencontainers.image.title="VNtyper" \
      org.opencontainers.image.description="VNtyper 2.0 - MUC1 VNTR genotyping pipeline for ADTKD-MUC1" \
      org.opencontainers.image.vendor="VNtyper Project" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.url="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.documentation="https://github.com/hassansaei/VNtyper" \
      com.vntyper.type="bioinformatics" \
      com.vntyper.pipeline="genotyping" \
      com.vntyper.build="multistage"

# Install mamba for faster package management
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
    conda install -y mamba -n base -c conda-forge && \
    conda clean -afy

# =============================================================================
# Stage 2: Build Environment (Contains all build tools)
# =============================================================================
FROM base AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build-time system dependencies with BuildKit cache
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        unzip \
        zip \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Define build arguments
ARG REPO_URL=https://github.com/hassansaei/VNtyper.git
ARG REPO_DIR=/opt/vntyper

# Clone repository
RUN git clone ${REPO_URL} ${REPO_DIR}

WORKDIR ${REPO_DIR}

# Copy environment files to temp location for better layer caching
RUN cp conda/environment_vntyper.yml /tmp/environment_vntyper.yml && \
    cp conda/environment_envadvntr.yml /tmp/environment_envadvntr.yml && \
    cp conda/environment_shark.yml /tmp/environment_shark.yml

# Create conda environments with aggressive cleanup
# Using cache mounts to speed up repeated builds
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
    --mount=type=cache,target=/root/.conda/pkgs,sharing=locked \
    mamba env create -f /tmp/environment_vntyper.yml && \
    mamba env create -f /tmp/environment_envadvntr.yml && \
    mamba env create -f /tmp/environment_shark.yml && \
    conda clean -afy && \
    find /opt/conda/envs -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda/envs -type f -name '*.pyc' -delete 2>/dev/null || true && \
    find /opt/conda/envs -type f -name '*.pyo' -delete 2>/dev/null || true && \
    rm -f /tmp/environment_*.yml

# Install VNtyper Python package
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install .

# Install adVNTR
RUN chmod +x vntyper/dependencies/advntr/install_advntr.sh && \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate envadvntr && \
    bash vntyper/dependencies/advntr/install_advntr.sh -o"

# Install reference genomes
RUN conda run -n vntyper vntyper --config-path /opt/vntyper/vntyper/config.json install-references \
    --output-dir /opt/vntyper/reference

# Extract adVNTR database if it was installed as a zip file
# Do this in builder stage where we have proper Python environment
RUN if [ -f /opt/vntyper/reference/vntr_db_advntr.zip ] && [ ! -d /opt/vntyper/reference/vntr_db_advntr ]; then \
        echo "Extracting adVNTR database in builder stage..."; \
        conda run -n vntyper python3 -c "import zipfile; z = zipfile.ZipFile('/opt/vntyper/reference/vntr_db_advntr.zip'); z.extractall('/opt/vntyper/reference'); print('Successfully extracted adVNTR database')" && \
        echo "adVNTR database extracted successfully" || \
        echo "WARNING: Failed to extract adVNTR database"; \
    fi

# Install FastAPI and web dependencies
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install \
        "fastapi[standard]==0.115.3" \
        "uvicorn[standard]==0.32.0" \
        redis==5.2.0 \
        celery==5.4.0 \
        python-multipart==0.0.12 \
        fastapi-limiter==0.1.6 \
        email_validator==2.2.0 \
        'passlib[bcrypt]==1.7.4' \
        pydantic==2.10.0

# Copy local entrypoint.sh and app directory to override cloned versions
# This allows local development changes without pushing to GitHub
# Note: Paths are relative to build context (repo root), not Dockerfile location
COPY docker/entrypoint.sh /opt/vntyper/docker/entrypoint.sh
COPY docker/app /opt/vntyper/app

# =============================================================================
# Stage 3: Minimal Runtime Image
# =============================================================================
FROM condaforge/miniforge3:24.11.3-0@sha256:bcfc04b3562faa8dfd48f3767e2bd7b73a8529f9611b508d42e5bd227db1342d AS runtime

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    PATH=/opt/conda/envs/vntyper/bin:/opt/conda/envs/environment_vntyper/bin:$PATH \
    DEFAULT_INPUT_DIR=/opt/vntyper/input \
    DEFAULT_OUTPUT_DIR=/opt/vntyper/output \
    REFERENCE_DIR=/opt/vntyper/reference

# Install ONLY runtime dependencies (minimal set - no build tools!)
# curl needed for health checks, ca-certificates for SSL
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy conda environments from builder (this is where the magic happens!)
# Conda environments are self-contained, so copying the entire envs directory works
COPY --from=builder /opt/conda/envs /opt/conda/envs

# Copy VNtyper application and references
# Using --chown to set ownership during copy for better layer efficiency
COPY --from=builder /opt/vntyper /opt/vntyper

# Verify adVNTR database presence (should have been extracted in builder stage)
RUN if [ -d /opt/vntyper/reference/vntr_db_advntr ]; then \
        echo "adVNTR database found and ready"; \
    else \
        echo "INFO: adVNTR database not found - adVNTR module will not be available (this is optional)"; \
    fi

# Set working directory
WORKDIR /opt/vntyper

# Set up default directories
RUN mkdir -p $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR

# Create non-root user with configurable UID/GID
ARG USERNAME=appuser
ARG USER_UID=1001
ARG USER_GID=1001

RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID --shell /bin/bash --create-home $USERNAME

# Set ownership of directories to non-root user
RUN chown -R $USERNAME:$USERNAME \
    $DEFAULT_INPUT_DIR \
    $DEFAULT_OUTPUT_DIR \
    $REFERENCE_DIR \
    /opt/vntyper

# Copy and set up entrypoint
COPY --from=builder /opt/vntyper/docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Copy FastAPI app
COPY --from=builder /opt/vntyper/docker/app /opt/vntyper/app
RUN chown -R $USERNAME:$USERNAME /opt/vntyper/app

# Expose port
EXPOSE 8000

# Switch to non-root user
USER $USERNAME

# Health check to verify container is running correctly
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD ["/usr/local/bin/entrypoint.sh", "health"]

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
