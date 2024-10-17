# Use Miniforge3 as the base image from conda-forge
FROM condaforge/miniforge3:latest

# Set environment variables to prevent interactive prompts during installations
ENV DEBIAN_FRONTEND=noninteractive

# Update and install necessary system packages
RUN apt-get update && \
    apt-get install -y build-essential git && \
    rm -rf /var/lib/apt/lists/*

# Install mamba for faster package management
RUN conda install -y mamba -n base -c conda-forge

# Create and activate the vntyper environment with Python 3.9 and dependencies
RUN mamba create -y -n vntyper python=3.9 \
    bioconda::bwa \
    bioconda::samtools \
    bioconda::fastp \
    pandas \
    numpy \
    regex \
    biopython \
    setuptools \
    pysam \
    jinja2 \
    matplotlib \
    seaborn \
    igv-reports \
    plotly

# Create and activate the advntr environment with Python 2.7 and dependencies
RUN mamba create -y -n advntr python=2.7 \
    setuptools \
    matplotlib \
    biopython \
    pysam \
    cython \
    networkx=1.11 \
    scipy \
    joblib \
    scikit-learn \
    numpy \
    tk \
    pip \
    zlib \
    samtools=1.11 \
    muscle=3.8.1551

# Clone the VNtyper repository from GitHub
RUN git clone https://github.com/berntpopp/VNtyper.git /opt/vntyper

# Set the working directory
WORKDIR /opt/vntyper

# Install vntyper using setup.py within the vntyper environment
RUN conda run -n vntyper pip install .

# (Optional) If advntr is a separate package, install it similarly
# RUN conda run -n advntr pip install path_to_advntr_package

# Create a non-root user for running the application
RUN useradd -m appuser
USER appuser

# Install advntr into the base environment if it's a standalone tool
# This assumes advntr can run within the base environment or needs separate handling

# Define the entry point to use vntyper as the default command
ENTRYPOINT ["/opt/conda/bin/conda", "run", "-n", "vntyper", "vntyper"]
