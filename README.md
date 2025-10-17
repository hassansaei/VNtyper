# VNtyper 2.0 - A Pipeline to genotype the MUC1-VNTR

**VNtyper 2.0** is an advanced pipeline designed to genotype MUC1 coding Variable Number Tandem Repeats (VNTR) in Autosomal Dominant Tubulointerstitial Kidney Disease (ADTKD-MUC1) using Short-Read Sequencing (SRS) data. This version is a refactored version of VNtyper v1 integrates enhanced variant calling algorithms, robust logging mechanisms, and streamlined installation processes to provide researchers with a powerful tool for VNTR analysis.

- We have developed a web server to provide free access to VNtyper, which runs in the background for ease of use.  
  Access it through the following link: [vntyper-online](https://vntyper.org/)

---

## Table of Contents

1. [Features](#features)  
2. [Installation](#installation)  
3. [Usage](#usage)  
4. [Pipeline Overview](#pipeline-overview)  
5. [Dependencies](#dependencies)  
6. [Linting and Code Formatting](#linting-and-code-formatting)  
7. [Pipeline Logic Diagram](#pipeline-logic-diagram)  
8. [Results](#results)  
9. [Notes](#notes)  
10. [Citations](#citations)  
11. [Contributing](#contributing)  
12. [License](#license)  
13. [Contact](#contact)

---

## Features

- **Variant Calling Algorithms:**
  - **Kestrel:** Mapping-free genotyping using k-mer frequencies.
  - **code-adVNTR (optional):** Profile-HMM-based method for VNTR genotyping.
  - **SHARK (optional, FASTQ-only):** Rapid filtering and read extraction for MUC1 region in exome/whole-genome data.

- **Comprehensive Logging:**
  - Logs both to the console and a dedicated log file.
  - Generates MD5 checksums for all downloaded and processed files.

- **Modern Packaging:**
  - Uses modern Python packaging with `pyproject.toml` (PEP 517/518/621).
  - No `setup.py` required - fully compatible with pip ≥21.3.
  - Provides Conda environment setup for easy dependency management.

- **Subcommands:**
  - `install-references`
  - `pipeline`
  - `fastq`
  - `bam`
  - `kestrel`
  - `report`
  - `cohort`
  - `online`

---

## Installation

VNtyper 2.0 uses modern Python packaging with `pyproject.toml` and can be installed using `pip` (≥21.3) or via Conda environments for streamlined dependency management.

### Using `pip`

1. **Clone the Repository:**

   ```bash
   mkdir vntyper
   git clone https://github.com/hassansaei/vntyper.git
   cd vntyper
   pip install .
   ```

2. **For Development:**

   ```bash
   # Install in editable mode with development dependencies
   pip install -e .[dev]
   ```

   This installs VNtyper with additional tools for development:
   - **Ruff**: Fast linter and formatter (replaces flake8 + black)
   - **pytest**: Testing framework with coverage support

---

## Usage

VNtyper 2.0 offers multiple subcommands that can be used depending on your input data and requirements. Below are the main subcommands available:

### 1. Running the Full Pipeline

To run the entire pipeline using a BAM file:

```bash
vntyper --config-path /path/to/config.json pipeline \
    --bam /path/to/sample.bam \
    --output-dir /path/to/output/dir \
    --threads 4 --fast-mode
```

Alternatively, using paired-end FASTQ files:

```bash
vntyper --config-path /path/to/config.json pipeline \
    --fastq1 /path/to/sample_R1.fastq.gz \
    --fastq2 /path/to/sample_R2.fastq.gz \
    --output-dir /path/to/output/dir \
    --threads 4 --fast-mode
```

The adVNTR genotyping is optional and skipped by default. To enable adVNTR genotyping, use the `--extra-modules advntr` option.

**New**: To enable SHARK filtering on FASTQ reads *before* the usual QC and alignment (for improved MUC1 detection), add `shark` to the `--extra-modules` flag (e.g., `--extra-modules shark`). This can be done as:

```bash
vntyper --config-path /path/to/config.json pipeline \
    --fastq1 /path/to/sample_R1.fastq.gz \
    --fastq2 /path/to/sample_R2.fastq.gz \
    --extra-modules shark \
    --threads 4 \
    --output-dir /path/to/output/dir
```

- SHARK will run first on the raw FASTQ files to extract and filter reads covering the MUC1 VNTR region.  
- **Important**: SHARK is only supported in FASTQ mode. If you try to use `--extra-modules shark` together with `--bam` or `--cram`, VNtyper will exit gracefully with a warning.

### 2. Running VNtyper with Docker

Docker image for VNtyper 2.0 is provided and can be pulled and used as follows:

```bash
# pull the docker image
docker pull saei/vntyper:main

# run the pipeline using the docker image
docker run -w /opt/vntyper --rm \
    -v /local/input/folder/:/opt/vntyper/input \
    -v /local/output/folder/:/opt/vntyper/output \
    saei/vntyper:latest \
    vntyper pipeline \
    --bam /opt/vntyper/input/filename.bam \
    -o /opt/vntyper/output/filename/
```

> **Important Host Volume Permissions Note:**  
> When mounting host directories into the container (using the `-v` flag), please ensure that the host directories (e.g., `/local/input/folder/` and `/local/output/folder/`) have the appropriate permissions so that they are writable by the container's non-root user.
>
> **Why Non-Root?**  
> VNtyper runs as a non-root user for enhanced security and to avoid file ownership issues on your host. Running as root may create files owned by root, leading to permission problems later.
>
> There are two ways to ensure proper permissions:
>
> 1. **Adjust Host Directory Permissions:**  
>    Change the ownership/permissions on the host directories so that the UID and GID match those expected by VNtyper in the container.
>
> 2. **Use the `--user` Flag:**  
>    Run the container with the `--user` flag to specify your current user’s UID and GID. For example:
>
>    ```bash
>    docker run --user $(id -u):$(id -g) -w /opt/vntyper --rm \
>      -v /local/input/folder/:/opt/vntyper/input \
>      -v /local/output/folder/:/opt/vntyper/output \
>      saei/vntyper:latest \
>      vntyper pipeline \
>      --bam /opt/vntyper/input/filename.bam \
>      -o /opt/vntyper/output/filename/
>    ```
>
> Using either method ensures VNtyper can write its log files (e.g., `pipeline.log`) and other outputs without encountering permission errors.

An Apptainer image can be generated from the Docker image as follows:

```bash
# create the apptainer sif image
apptainer pull docker://saei/vntyper:main

# run the pipeline using the apptainer image
apptainer run --pwd /opt/vntyper \
    -B /local/input/folder/:/opt/vntyper/input \
    -B /local/output/folder/:/opt/vntyper/output \
    vntyper_latest.sif vntyper pipeline \
    --bam /opt/vntyper/input/filename.bam \
    -o /opt/vntyper/output/filename/
```

### 3. Installing References

```bash
vntyper --config-path /path/to/config.json install-references \
    --output-dir /path/to/reference/install \
    --skip-indexing  # Optional: skip BWA indexing if needed
```

### 4. Generating Reports

```bash
vntyper --config-path /path/to/config.json report \
    --output-dir /path/to/output/dir
```

---

## Pipeline Overview

VNtyper 2.0 integrates multiple steps into a streamlined pipeline. The following is an overview of the steps involved:

1. **FASTQ Quality Control**: Raw FASTQ files are checked for quality.  
2. **(Optional) SHARK Filtering**: If `shark` is specified in `--extra-modules`, raw FASTQ reads are first filtered to extract MUC1-specific reads (especially relevant for exome or large WGS datasets).  
3. **Alignment**: Reads are aligned using BWA (if FASTQ files are provided).  
4. **Kestrel Genotyping**: Mapping-free genotyping of VNTRs.  
5. **(Optional) adVNTR Genotyping**: Profile-HMM-based method for VNTR genotyping (requires additional setup).  
6. **Summary Report Generation**: A final HTML report is generated to summarize the results.

---

## Dependencies

VNtyper 2.0 relies on several tools and Python libraries. Ensure that the following dependencies are available in your environment:

- Python >= 3.9  
- BWA  
- Samtools  
- Fastp  
- Pandas  
- Numpy  
- Biopython  
- Pysam  
- Jinja2  
- Matplotlib  
- Seaborn  
- IGV-Reports  

You can easily set up these dependencies via the provided Conda environment file.

---

## Linting and Code Formatting

VNtyper adheres to PEP8 style guidelines to ensure clean, readable, and maintainable code. We use **Ruff**, a modern, extremely fast linter and formatter written in Rust that replaces both flake8 and black.

### Quick Start with Makefile

VNtyper provides a Makefile for common development tasks:

```bash
# Show all available commands
make help

# Install in development mode
make install-dev

# Run linter
make lint

# Format code (auto-fix issues)
make format

# Check formatting without making changes
make format-check

# Run tests
make test          # All tests
make test-unit     # Unit tests only
make test-cov      # With coverage report

# Clean build artifacts
make clean
```

### Modern Python Packaging

VNtyper uses `pyproject.toml` for all project configuration (PEP 517/518/621):
- **No `setup.py` required**: Modern pip (≥21.3) and setuptools (≥64) fully support `pyproject.toml`
- **Single source of truth**: All configuration (build, dependencies, tools) in one file
- **Editable installs**: Fully supported via PEP 660

---

## Pipeline Logic Diagram

Below is a logical overview of the VNtyper pipeline:

```mermaid
graph TD
  A[Input: FASTQ/BAM] -->|Quality Control| B[Alignment BWA]
  B -->|Genotyping| C[Kestrel]
  C --> D[Optional: adVNTR]
  D --> E[Generate Summary Report]
  E --> F[Output: VCF, Summary HTML]
```

---

## Results

Once the pipeline completes, you will have:

- **BAM or FASTQ** slices containing MUC1-specific reads.  
- **VCF files** or **TSV files** with genotyping results (for Kestrel and optional adVNTR).  
- **HTML summary report** including:
  - **VNTR Region Coverage Statistics**: Detailed coverage metrics specifically for the VNTR region, including mean, median, standard deviation, minimum, and maximum coverage, as well as the percentage of the VNTR region with zero coverage.
  - **Genotyping Calls**: Results from Kestrel and optional adVNTR analyses.
  - **Quality Metrics**: When available, includes duplication rate, Q20/Q30 rates, and other quality indicators.
  - **Pipeline Log**: Comprehensive logging information about the pipeline execution.

---

## Notes

1. This tool is for **research use only**.  
2. Ensure **high-coverage WES/WGS or targeted data** is used to genotype MUC1 VNTR accurately.  
3. For questions or issues, refer to the GitHub repository for support.

---

## Citations

If you use VNtyper 2.0 in your research, please cite the following:

1. Saei H, Morinière V, Heidet L, et al. VNtyper enables accurate alignment-free genotyping of MUC1 coding VNTR using short-read sequencing data. iScience. 2023.  
2. Audano PA, Ravishankar S, et al. Mapping-free variant calling using haplotype reconstruction from k-mer frequencies. Bioinformatics. 2018.  
3. Park J, Bakhtiari M, et al. Detecting tandem repeat variants in coding regions using code-adVNTR. iScience. 2022.

---

## Contributing

We welcome contributions to VNtyper. Please refer to the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines.

---

## License

VNtyper is licensed under the BSD 3-Clause License. See the LICENSE file for more details.

---

## Contact

For questions or issues, please open an [issue on GitHub](https://github.com/hassansaei/vntyper/issues) or email the corresponding authors listed in the manuscript.
