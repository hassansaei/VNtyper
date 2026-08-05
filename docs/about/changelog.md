# Changelog

All notable changes to VNtyper 2 are documented on this page.

## 2.0.4 (Current)

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
