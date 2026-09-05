---
hide:
  - navigation
  - toc
---

# VNtyper 2

**Genotype MUC1 coding VNTRs for ADTKD-MUC1 diagnosis using short-read sequencing.**

VNtyper 2 is a bioinformatics pipeline that detects frameshift mutations in the MUC1
Variable Number Tandem Repeat (VNTR) region, the genetic cause of Autosomal Dominant
Tubulointerstitial Kidney Disease (ADTKD-MUC1). It combines mapping-free k-mer
genotyping (Kestrel) with optional Profile-HMM validation (adVNTR) to deliver
confidence-scored variant calls from BAM, CRAM, or FASTQ input.

<div class="grid cards" markdown>

-   :material-dna:{ .lg .middle } **Mapping-Free Genotyping**

    ---

    Kestrel k-mer assembly avoids reference bias in repetitive VNTR regions,
    delivering empirically validated confidence scoring.

    [:octicons-arrow-right-24: How it works](pipeline/kestrel.md)

-   :material-file-multiple:{ .lg .middle } **Flexible Input**

    ---

    Accepts BAM, CRAM, or paired-end FASTQ files with support for 8 reference
    assembly conventions across UCSC, NCBI RefSeq, and Ensembl.

    [:octicons-arrow-right-24: Input formats](user-guide/input-formats.md)

-   :material-chart-bar:{ .lg .middle } **Interactive Reports**

    ---

    HTML reports with embedded IGV genome browser tracks, coverage metrics,
    and cohort-level summaries with optional pseudonymization.

    [:octicons-arrow-right-24: Output guide](user-guide/output-files.md)

</div>

## Quick Install

```bash
# Recommended: install with all binary dependencies via Conda
git clone https://github.com/hassansaei/vntyper.git
cd vntyper
conda env create -f conda/environment_vntyper.yml
conda activate vntyper
pip install -e .

# Install reference assets
vntyper install-references -d ./references

# Run pipeline
vntyper pipeline --bam inputs/sample.bam -o results/sample/
```

[:octicons-rocket-24: Get started](getting-started/quickstart.md){ .md-button .md-button--primary }
[:octicons-book-24: CLI Reference](cli/index.md){ .md-button }
