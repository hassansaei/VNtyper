# Installation

## Prerequisites

VNtyper 2 requires the following to be installed on your system:

| Dependency   | Minimum Version | Purpose                          |
|-------------|-----------------|----------------------------------|
| Python      | 3.9+            | Core runtime                     |
| Java        | 11+             | Kestrel genotyping engine        |
| BWA         | 0.7.18+         | Read alignment                   |
| samtools    | 1.20+           | BAM/CRAM manipulation            |
| fastp       | 0.23+           | FASTQ quality control            |

!!! warning "External tools required"
    External tools (BWA, samtools, fastp, Java 11) must be installed separately when using pip or source installation. The Docker image includes all dependencies.

## Install VNtyper 2

=== "pip"

    Install directly from GitHub:

    ```bash
    pip install git+https://github.com/hassansaei/VNtyper.git
    ```

=== "From Source"

    Clone the repository and install in editable mode:

    ```bash
    git clone https://github.com/hassansaei/vntyper.git
    cd vntyper
    pip install -e .
    ```

=== "Conda"

    Use the provided environment file to create a Conda environment with all dependencies (including external tools):

    ```bash
    git clone https://github.com/hassansaei/vntyper.git
    cd vntyper

    # Create the environment
    conda env create -f conda/environment_vntyper.yml
    conda activate vntyper

    # Install VNtyper 2 into the environment
    pip install -e .
    ```

    The `environment_vntyper.yml` includes Python 3.9, BWA, samtools, fastp, OpenJDK 11, and all Python dependencies.

    Additional environment files are available for optional modules:

    - `conda/environment_envadvntr.yml` --- adVNTR genotyping module
    - `conda/environment_shark.yml` --- SHARK read filtering module

=== "Docker"

    Pull the pre-built image from Docker Hub or GitHub Container Registry:

    ```bash
    # Docker Hub
    docker pull saei/vntyper:latest

    # GitHub Container Registry
    docker pull ghcr.io/hassansaei/vntyper:latest
    ```

    Or build from source:

    ```bash
    git clone https://github.com/hassansaei/VNtyper.git
    cd VNtyper
    DOCKER_BUILDKIT=1 docker build -f docker/Dockerfile -t vntyper:latest .
    ```

## Verify Installation

Confirm VNtyper 2 is installed correctly:

```bash
vntyper --version
```

## Development Install

To install with development dependencies (Ruff linter, pytest, mypy):

```bash
git clone https://github.com/hassansaei/vntyper.git
cd vntyper
pip install -e .[dev]
```

Run the test suite to verify everything works:

```bash
make test-unit
```

See the [Contributing Guide](https://github.com/hassansaei/VNtyper/blob/main/CONTRIBUTING.md) for full development setup instructions.
