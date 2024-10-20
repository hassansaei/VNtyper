# Use a specific Miniforge3 version for reproducibility
FROM condaforge/miniforge3:4.12.0-0

# Set environment variables to prevent interactive prompts during installations
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    PATH=/opt/conda/envs/vntyper/bin:/opt/conda/envs/environment_vntyper/bin:$PATH

# Update and install necessary system packages
RUN apt-get update && \
    apt-get install -y build-essential git jq wget unzip && \
    rm -rf /var/lib/apt/lists/*

# Install mamba for faster package management
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

# Create the vntyper environment from the environment.yml
RUN mamba env create -f conda/environment_vntyper.yml && \
    conda clean -afy

# Install vntyper using setup.py within the vntyper environment
RUN conda run -n vntyper pip install .

# Create the advntr environment from the environment.yml
RUN mamba env create -f conda/environment_envadvntr.yml && \
    conda clean -afy

# Make the adVNTR installation script executable
RUN chmod +x vntyper/dependencies/advntr/install_advntr.sh

# Install adVNTR in the envadvntr environment
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate envadvntr && \
    bash vntyper/dependencies/advntr/install_advntr.sh"

# Create a non-root user for running the application
RUN useradd -m appuser
USER appuser

# Define the entry point to use vntyper as the default command
ENTRYPOINT ["/opt/conda/bin/conda", "run", "-n", "vntyper", "vntyper"]

# Add HEALTHCHECK to ensure the container is running correctly
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
  CMD ["vntyper", "--version"] || exit 1
