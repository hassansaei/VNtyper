# Use a specific Miniforge3 version for reproducibility
FROM condaforge/miniforge3:4.12.0-0

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    PATH=/opt/conda/envs/vntyper/bin:/opt/conda/envs/environment_vntyper/bin:$PATH \
    DEFAULT_INPUT_DIR=/opt/vntyper/input \
    DEFAULT_OUTPUT_DIR=/opt/vntyper/output \
    REFERENCE_DIR=/opt/vntyper/references

# Update and install necessary system packages
RUN apt-get update && \
    apt-get install -y build-essential git jq wget unzip && \
    rm -rf /var/lib/apt/lists/*

# Install mamba for faster conda package installation
RUN conda install -y mamba -n base -c conda-forge

# Define build arguments with default values
ARG REPO_URL=https://github.com/berntpopp/VNtyper.git
ARG REPO_DIR=/opt/vntyper

# Set default repository URL and directory
ENV REPO_URL=${REPO_URL}
ENV REPO_DIR=${REPO_DIR}

# Clone the VNtyper repository using environment variables
RUN git clone ${REPO_URL} ${REPO_DIR}

# Set the working directory
WORKDIR ${REPO_DIR}

# Copy VNtyper core code and environment files
COPY conda/environment_vntyper.yml /tmp/environment_vntyper.yml
COPY conda/environment_envadvntr.yml /tmp/environment_envadvntr.yml

# Create and install conda environments
RUN mamba env create -f /tmp/environment_vntyper.yml && \
    mamba env create -f /tmp/environment_envadvntr.yml && \
    conda clean -afy && \
    rm /tmp/environment_vntyper.yml /tmp/environment_envadvntr.yml

# Install VNtyper using setup.py within the vntyper environment
RUN conda run -n vntyper pip install .

# Make adVNTR installation script executable
RUN chmod +x vntyper/dependencies/advntr/install_advntr.sh

# Install adVNTR environment
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate envadvntr && \
    bash vntyper/dependencies/advntr/install_advntr.sh"

# Install VNtyper references during build
RUN conda run -n vntyper vntyper install-references \
    --output-dir $REFERENCE_DIR \
    --config-path /opt/vntyper/config/config.json \
    --skip-indexing

# Install FastAPI and Uvicorn for web service
RUN conda run -n vntyper pip install fastapi uvicorn[standard] redis celery

# Set up default input and output directories
RUN mkdir -p /opt/vntyper/input /opt/vntyper/output

# Set up entrypoint
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Copy FastAPI app code to the container
COPY docker/app /opt/vntyper/app

# Expose FastAPI port
EXPOSE 8000

# Switch to a non-root user
RUN useradd -m appuser
USER appuser

# Define the entry point to either run VNtyper or start the FastAPI server
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

# Add HEALTHCHECK to ensure the container is running correctly
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
  CMD ["vntyper", "--version"] || exit 1
