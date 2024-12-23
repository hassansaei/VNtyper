# Use a specific Miniforge3 version for reproducibility
FROM condaforge/miniforge3:4.12.0-0

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    PATH=/opt/conda/envs/vntyper/bin:/opt/conda/envs/environment_vntyper/bin:$PATH \
    DEFAULT_INPUT_DIR=/opt/vntyper/input \
    DEFAULT_OUTPUT_DIR=/opt/vntyper/output \
    REFERENCE_DIR=/opt/vntyper/reference

# Update and install necessary system packages in a single RUN command to reduce layers
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        unzip \
        zip \
        curl \
        && rm -rf /var/lib/apt/lists/*

# Install mamba for faster conda package installation
RUN conda install -y mamba -n base -c conda-forge

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

# Copy VNtyper environment files using RUN cp
RUN cp conda/environment_vntyper.yml /tmp/environment_vntyper.yml && \
    cp conda/environment_envadvntr.yml /tmp/environment_envadvntr.yml && \
    cp conda/environment_shark.yml /tmp/environment_shark.yml
 
# Create and install conda environments
RUN mamba env create -f /tmp/environment_vntyper.yml && \
    mamba env create -f /tmp/environment_envadvntr.yml && \
    mamba env create -f /tmp/environment_shark.yml && \
    conda clean -afy && \
    rm /tmp/environment_vntyper.yml /tmp/environment_envadvntr.yml /tmp/environment_shark.yml

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
    --config-path /opt/vntyper/config/config.json

# Install FastAPI, Uvicorn, Celery, Redis, and other dependencies
RUN conda run -n vntyper pip install fastapi==0.115.3 uvicorn==0.32.0 redis==5.2.0 celery==5.4.0 python-multipart==0.0.12 fastapi-limiter==0.1.6 email_validator==2.2.0 passlib[bcrypt]==1.7.4 pydantic==2.10.0

# Set up default input, output, and reference directories
RUN mkdir -p $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR

# Create a non-root user for running the application
RUN useradd -m appuser

# Set up entrypoint using RUN cp commands instead of COPY
RUN cp docker/entrypoint.sh /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

# Copy FastAPI app code to the container using RUN cp
RUN cp -r docker/app /opt/vntyper/app

# Set ownership of directories and app to appuser
RUN chown -R appuser:appuser $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR /opt/vntyper/app

# Expose FastAPI port
EXPOSE 8000

# Switch to a non-root user
USER appuser

# Define the entry point to either run VNtyper, Celery worker, Celery Beat, or start the FastAPI server
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
