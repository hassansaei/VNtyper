# syntax=docker/dockerfile:1.4
# Use a specific Miniforge3 version for reproducibility
FROM condaforge/miniforge3:24.11.3-0

# OCI Labels for metadata
LABEL org.opencontainers.image.title="VNtyper" \
      org.opencontainers.image.description="VNtyper 2.0 - MUC1 VNTR genotyping pipeline for ADTKD-MUC1" \
      org.opencontainers.image.vendor="VNtyper Project" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.url="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.documentation="https://github.com/hassansaei/VNtyper" \
      com.vntyper.type="bioinformatics" \
      com.vntyper.pipeline="genotyping"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    PATH=/opt/conda/envs/vntyper/bin:/opt/conda/envs/environment_vntyper/bin:$PATH \
    DEFAULT_INPUT_DIR=/opt/vntyper/input \
    DEFAULT_OUTPUT_DIR=/opt/vntyper/output \
    REFERENCE_DIR=/opt/vntyper/reference

# Update and install necessary system packages with BuildKit cache mount
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

# Install mamba for faster conda package installation with BuildKit cache
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
    conda install -y mamba -n base -c conda-forge && \
    conda clean -afy

# Define build arguments with default values
ARG REPO_URL=https://github.com/hassansaei/VNtyper.git
ARG REPO_DIR=/opt/vntyper

# Set default repository URL and directory
ENV REPO_URL=${REPO_URL}
ENV REPO_DIR=${REPO_DIR}

# Clone the VNtyper repository using environment variables
RUN git clone ${REPO_URL} ${REPO_DIR}

# Set the working directory
WORKDIR ${REPO_DIR}

# Copy VNtyper environment files (using RUN cp to work from cloned repo)
RUN cp conda/environment_vntyper.yml /tmp/environment_vntyper.yml && \
    cp conda/environment_envadvntr.yml /tmp/environment_envadvntr.yml && \
    cp conda/environment_shark.yml /tmp/environment_shark.yml

# Create and install conda environments with BuildKit cache mounts
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
    --mount=type=cache,target=/root/.conda/pkgs,sharing=locked \
    mamba env create -f /tmp/environment_vntyper.yml && \
    mamba env create -f /tmp/environment_envadvntr.yml && \
    mamba env create -f /tmp/environment_shark.yml && \
    conda clean -afy && \
    find /opt/conda/envs -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda/envs -type f -name '*.pyc' -delete 2>/dev/null || true && \
    find /opt/conda/envs -type f -name '*.pyo' -delete 2>/dev/null || true && \
    rm /tmp/environment_vntyper.yml /tmp/environment_envadvntr.yml /tmp/environment_shark.yml

# Install VNtyper using pip within the vntyper environment with BuildKit cache
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install --no-cache-dir .

# Make adVNTR installation script executable
RUN chmod +x vntyper/dependencies/advntr/install_advntr.sh

# Install adVNTR environment
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate envadvntr && \
    bash vntyper/dependencies/advntr/install_advntr.sh -o"

# Install VNtyper references during build (hg19 and hg38 UCSC by default)
RUN conda run -n vntyper vntyper --config-path /opt/vntyper/vntyper/config.json install-references \
    --output-dir $REFERENCE_DIR

# Install FastAPI, Uvicorn, Celery, Redis, and other dependencies with BuildKit cache
# Note: Using fastapi[standard] to include fastapi-cli for 'fastapi dev' command
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install --no-cache-dir \
        "fastapi[standard]==0.115.3" \
        "uvicorn[standard]==0.32.0" \
        redis==5.2.0 \
        celery==5.4.0 \
        python-multipart==0.0.12 \
        fastapi-limiter==0.1.6 \
        email_validator==2.2.0 \
        'passlib[bcrypt]==1.7.4' \
        pydantic==2.10.0

# Set up default input, output, and reference directories
RUN mkdir -p $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR

# Allow dynamic UID and GID as build arguments
ARG USERNAME=appuser
ARG USER_UID=1001
ARG USER_GID=1001

# Create a group and user with the specified UID and GID
RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID --shell /bin/bash --create-home $USERNAME

# Copy entrypoint and app files (using COPY for better practice, but files are in cloned repo)
# Note: These files come from the cloned repo, so we use RUN cp from within the repo
RUN cp docker/entrypoint.sh /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

# Copy FastAPI app code to the container
RUN cp -r docker/app /opt/vntyper/app

# Set ownership of directories and app to appuser
RUN chown -R $USERNAME:$USERNAME $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR /opt/vntyper/app /opt

# Expose FastAPI port
EXPOSE 8000

# Switch to a non-root user
USER $USERNAME

# Health check to verify container is running correctly
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD ["/usr/local/bin/entrypoint.sh", "health"]

# Define the entry point to either run VNtyper, Celery worker, Celery Beat, or start the FastAPI server
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
