# Changelog

All notable changes to VNtyper 2 are documented on this page.

## 2.0.1 (Current)

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

The original VNtyper 2 release as described in Saei et al., *iScience* 2023. This version provided the initial alignment-free genotyping approach for MUC1 VNTR using Kestrel and adVNTR.
