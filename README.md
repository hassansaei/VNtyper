# VNtyper 2 - A Pipeline to genotype the MUC1-VNTR

[![CI tests](https://github.com/vntyper-org/VNtyper/actions/workflows/ci-tests.yml/badge.svg?branch=main)](https://github.com/vntyper-org/VNtyper/actions/workflows/ci-tests.yml)
[![PyPI version](https://img.shields.io/pypi/v/vntyper)](https://pypi.org/project/vntyper/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/vntyper)](https://pypi.org/project/vntyper/)
[![Documentation](https://img.shields.io/badge/docs-vntyper--org.github.io-blue)](https://vntyper-org.github.io/VNtyper/)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-green)](https://github.com/vntyper-org/VNtyper/blob/main/LICENSE)
[![Zenodo DOI 10.5281/zenodo.19744166](https://zenodo.org/badge/484326398.svg)](https://doi.org/10.5281/zenodo.19744166)

VNtyper 2 genotypes the coding variable number tandem repeat (VNTR) of *MUC1* from
short-read sequencing data, to find the frameshift variants that cause autosomal
dominant tubulointerstitial kidney disease (ADTKD-MUC1). It reads BAM, CRAM or paired
FASTQ files, calls variants with Kestrel (mapping-free k-mer haplotype assembly), and
can confirm them with code-adVNTR (profile-HMM repeat-unit genotyping). Each run writes
an HTML report with an embedded IGV.js view, a JSON summary with full provenance, and
VCF and TSV call tables.

VNtyper is for research use only. To analyse a single sample without installing
anything, use the web service at [vntyper.org](https://vntyper.org/).

## Installation

The Conda environment brings the external tools the pipeline calls (`bwa`, `samtools`,
`fastp`, `bcftools` and Java 11 for Kestrel):

```bash
git clone https://github.com/vntyper-org/VNtyper.git
cd VNtyper
mamba env create -f conda/environment_vntyper.yml
conda activate vntyper
pip install -e .
vntyper install-references --output-dir reference
```

`install-references` downloads the reference bundle into `./reference/` and checks
every file against its SHA-256 checksum. Run `vntyper pipeline` from the repository
root as well: the paths in `vntyper/config.json` are relative to the working directory.

`pip install vntyper` installs the package from PyPI. The `report`, `cohort` and
`calibrate` subcommands then work as they are; `pipeline` additionally needs the
external tools above on your `PATH`. The
[installation guide](https://vntyper-org.github.io/VNtyper/getting-started/installation/)
covers every route, including optional modules.

## Quick start

Download the test dataset (about 1.1 GB, from Zenodo) and genotype one of its BAM files:

```bash
make download-test-data
vntyper pipeline \
    --bam tests/data/example_6449_hg19_subset.bam \
    --output-dir results/example_6449 \
    --threads 4
```

The run takes a few seconds and ends with `Pipeline finished successfully.`; open
`results/example_6449/summary_report.html` to see the call. Keep `--output-dir` outside
the directory that holds the input alignment. For CRAM input, add `--cram` and
`--reference-fasta`; for FASTQ input, use `--fastq1` and `--fastq2`. Enable code-adVNTR
with `--extra-modules advntr`.

## Containers

Released images are published to GHCR as `latest` (the newest release), the immutable
`vX.Y.Z` and `X.Y.Z` tags naming one exact release, and the moving `X` and `X.Y` series
tags. Pin `vX.Y.Z` for a reproducible run. The `main` tag is rolling and tracks the
default branch. The images include the reference bundle. Mount input and output as
separate host directories:

```bash
docker pull ghcr.io/vntyper-org/vntyper:latest
docker run --user $(id -u):$(id -g) -w /opt/vntyper --rm \
    -v /local/inputs/:/opt/vntyper/input \
    -v /local/results/:/opt/vntyper/output \
    ghcr.io/vntyper-org/vntyper:latest \
    vntyper pipeline --bam /opt/vntyper/input/sample.bam -o /opt/vntyper/output/sample/
```

For Apptainer, run `apptainer pull docker://ghcr.io/vntyper-org/vntyper:latest` and
see [Docker and containers](https://vntyper-org.github.io/VNtyper/user-guide/docker/).

## Documentation

The documentation at [vntyper-org.github.io/VNtyper](https://vntyper-org.github.io/VNtyper/)
covers:

- [running the pipeline](https://vntyper-org.github.io/VNtyper/user-guide/running-pipeline/), including `--resume` and the optional SHARK read filter
- [output files](https://vntyper-org.github.io/VNtyper/user-guide/output-files/) and how confidence grades are assigned
- [the `report`, `cohort`, `online`, `install-references` and `calibrate` subcommands](https://vntyper-org.github.io/VNtyper/cli/)
- [the pipeline stages and methods](https://vntyper-org.github.io/VNtyper/pipeline/)
- [troubleshooting](https://vntyper-org.github.io/VNtyper/user-guide/troubleshooting/)

Related repositories:
[vntyper-data](https://github.com/vntyper-org/vntyper-data) (reference bundles),
[vntyper-online-frontend](https://github.com/vntyper-org/vntyper-online-frontend) and
[vntyper-online-backend](https://github.com/vntyper-org/vntyper-online-backend) (the web service),
and [MucOneUp](https://github.com/vntyper-org/MucOneUp) (MUC1 haplotype and read simulator).

## Getting help

Report bugs and ask questions in the
[issue tracker](https://github.com/vntyper-org/VNtyper/issues).

## Citing

If you use VNtyper 2 in your research, please cite:

> Popp B, Saei H, Teltsh O, et al. VNtyper 2 enables open-access short-read genotyping
> of MUC1 VNTR variants in ADTKD at high-speed. *medRxiv* (2026).
> [doi:10.64898/2026.05.27.26352937](https://doi.org/10.64898/2026.05.27.26352937)

and the original method:

> Saei H, Morinière V, Heidet L, et al. VNtyper enables accurate alignment-free
> genotyping of MUC1 coding VNTR using short-read sequencing data. *iScience*
> 26(7):107171 (2023). [doi:10.1016/j.isci.2023.107171](https://doi.org/10.1016/j.isci.2023.107171)

GitHub's "Cite this repository" button reads [`CITATION.cff`](https://github.com/vntyper-org/VNtyper/blob/main/CITATION.cff).
BibTeX entries, and citations for Kestrel, code-adVNTR and SHARK, are on the
[citation page](https://vntyper-org.github.io/VNtyper/about/citation/).

## License

VNtyper is released under the [BSD 3-Clause License](https://github.com/vntyper-org/VNtyper/blob/main/LICENSE).
Contributions are welcome; see [CONTRIBUTING.md](https://github.com/vntyper-org/VNtyper/blob/main/CONTRIBUTING.md)
and [AGENTS.md](https://github.com/vntyper-org/VNtyper/blob/main/AGENTS.md).
