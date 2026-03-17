# Frequently Asked Questions

## What input coverage do I need?

A minimum of approximately **100x coverage** over the MUC1 VNTR region is recommended for reliable variant detection. Whole-genome sequencing at 30x or higher typically provides sufficient local coverage, while whole-exome sequencing may fall short depending on the capture kit's coverage of the MUC1 VNTR.

## BAM vs FASTQ -- which is faster?

**BAM input is faster.** When you provide a BAM or CRAM file, VNtyper 2 extracts only the reads overlapping the MUC1 region before processing. Starting from FASTQ files requires an additional alignment step, which adds runtime.

## Do I need adVNTR?

No. adVNTR is an **optional** module that provides independent, alignment-based validation of Kestrel calls. Enabling it adds approximately 9 minutes to the runtime. It can be useful for confirmation in research settings but is not required for routine genotyping.

## What does Low_Precision mean?

A variant classified as **Low_Precision** (also called Low Confidence) was detected but has marginal depth support -- the depth score falls between the low and high thresholds defined in `kestrel_config.json`. These calls may benefit from manual review or independent validation. See [Scoring & Confidence](../pipeline/scoring-and-confidence.md) for detailed threshold definitions.

## Can I use GRCh38?

Yes. VNtyper 2 supports both UCSC and NCBI naming conventions. Use the `--reference-assembly` flag with any of the following values:

- `hg19` or `GRCh37`
- `hg38` or `GRCh38`

Example:

```bash
vntyper pipeline --bam input.bam --reference-assembly hg38 -o output/
```

## Docker vs local install?

**Docker** bundles all external dependencies (BWA, samtools, fastp, Java 11, Kestrel JAR) into a single image, so you do not need to install them yourself. This is the easiest way to get started.

A **local install** gives you more control and avoids container overhead, but you must ensure that BWA, samtools, fastp, and Java 11 are available on your PATH.

## How do I interpret the HTML report?

The HTML report includes an embedded IGV viewer for visual inspection of variants and a summary table of detected mutations with confidence levels. For a detailed walkthrough, see [Output Files](../user-guide/output-files.md).

## SHARK fails with BAM input

SHARK requires FASTQ input. If you want to use the SHARK module, provide reads with the `--fastq1` and `--fastq2` flags along with `--extra-modules shark`:

```bash
vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
    --extra-modules shark -o output/
```

## How do I run multiple samples?

There are two approaches:

1. **Cohort command** -- VNtyper 2's built-in `cohort` subcommand processes a directory of samples and produces an aggregated summary. See [Cohort Analysis](../user-guide/cohort-analysis.md).
2. **Snakemake workflow** -- For large-scale batch processing with parallelization and cluster support, use the provided Snakemake workflow. See [Snakemake](../user-guide/snakemake.md).

## Where can I get help?

Open an issue on [GitHub](https://github.com/hassansaei/VNtyper/issues). Please search existing issues first to avoid duplicates.
