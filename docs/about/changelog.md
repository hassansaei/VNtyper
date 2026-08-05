# Changelog

All notable changes to VNtyper 2 are documented on this page.

## 2.0.5 (Current)

**Docker image now runs Python 3.12.** Python 3.9 reached end of life on 2025-10-31.
Every package on the numerical path is unchanged - bwa, samtools, fastp, bcftools,
openjdk (Kestrel), pysam, numpy, pandas and biopython all resolve to the same versions -
and the pipeline reproduces the expected genotype exactly. Only the reporting stack
moved (matplotlib, plotly, igv-reports, jinja2). The package now requires Python 3.10+
and is tested on 3.10-3.13.

**CI/CD pipeline rebuilt.** The Docker image is split into a rarely-rebuilt base
(conda environments, adVNTR, reference genomes) and a per-commit application image, so a
commit now rebuilds in about 3 minutes instead of 40-70. Images are published only after
their tests pass, and the build uses the repository's own source rather than cloning
GitHub, so a pull request now tests the code in that pull request.

- Fixed `pyproject.toml` declaring `numpy<2.0.0` while the conda environment installed
  numpy 2.0.2: `pip install .` inside the image downgraded conda's numpy, so the
  published image ran a different numerical stack than its own recipe declared.
- Fixed the TSV/CSV summary parsers silently truncating rows whose field count did not
  match the header; a malformed row is now logged and skipped without discarding the
  rest of the file.
- Recovered 30 unit tests that carried no `unit` marker and had therefore never run in
  CI, covering the Issue #136/#145 genotype tie-breaking logic.
- Added a fast image smoke tier that verifies every reference path `config.json`
  declares actually exists in the image.
- Added coverage reporting with a ratcheting floor and an 80% target.
- `make check-all` now runs in seconds and works on a fresh clone; it previously
  required a 1.1 GB download and a Docker daemon.

## 2.0.4

- Migrated all logging to per-module loggers (`logging.getLogger(__name__)`); log
  records now carry the emitting module name instead of `root`.
- Replaced deprecated `datetime.utcnow()` calls; emitted timestamps are unchanged.
- Ruff lint rules are now selected explicitly instead of extending ruff's defaults, so
  lint results no longer change when ruff widens its default rule set.
- mypy configuration moved from Makefile flags into `pyproject.toml`.
- Fixed `.gitignore` ignoring the tracked `docs/` directory, which silently hid newly
  added documentation pages from Git.
- Fixed `make test-quiet`, which used a `--log-cli=false` flag pytest does not accept.
- Added `AGENTS.md` and `CLAUDE.md` with repository instructions for coding agents.

## 2.0.3

- Raised Kestrel max align and haplotype states to 40, fixing GDP inflation
  ([#156](https://github.com/hassansaei/VNtyper/issues/156)).
- Handled `pd.NA` and `None` in flagging condition evaluation to prevent `TypeError`s.
- Added a GDP inflation guard test with an anonymized hg38 sample.
- Added `CITATION.cff` with the Zenodo concept DOI.

## 2.0.2

- Removed specific motifs from `exclude_motifs_right` to prevent unwanted motif filtering.
- Clarified depth score conditions in confidence assignment
  ([#154](https://github.com/hassansaei/VNtyper/issues/154)).
- Added a dynamic version variable to the docs via `mkdocs-macros-plugin`.
- Modernized docs deployment to use GitHub Actions Pages.
- Added a CFLAGS workaround for adVNTR compilation errors on GCC 14+.

## 2.0.1

- Disabled duplicate flagging in Kestrel configuration.
- Cleared `motifs_for_alt_gg` array to prevent unintended variant filtering.
- Fixed flagging to occur before variant selection, preventing selection of flagged variants ([#145](https://github.com/hassansaei/VNtyper/issues/145)).
- Fixed confidence assignment to filter sub-threshold variants at the root rather than via override.
- Raised `Low_Depth_Conserved_Motifs` threshold from 0.2 to 0.4.
- Updated variant selection to prioritize depth score before haplotype count.
- Extracted frameshift deduplication into a helper function for clarity.
- Added comprehensive unit tests for the flagging module.

## 2.0.0

A complete rewrite and modernization of the VNtyper 2 pipeline.

- **Modern Python packaging** using `pyproject.toml` (PEP 517/518/621).
- **Refactored Kestrel postprocessing** with configurable thresholds via `kestrel_config.json`.
- **HTML reports** with embedded IGV integration for interactive variant inspection.
- **Cohort analysis** command with built-in pseudonymization for multi-sample studies.
- **Docker multi-stage build** for reproducible, lightweight container images.
- **Comprehensive test suite** with unit and integration tests, including downloadable test data.
- **Multiple reference assemblies** supported: hg19, hg38, GRCh37, and GRCh38.
- **Modular architecture** separating variant parsing, motif processing, scoring, flagging, and confidence assignment into individual modules.
- **Snakemake workflow** for scalable batch processing.

## 1.x

The original VNtyper release as described in Saei et al., *iScience* 2023. This version provided the initial alignment-free genotyping approach for MUC1 VNTR using Kestrel and adVNTR.
