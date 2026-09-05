# Frequently Asked Questions

## What input coverage do I need?

A minimum of approximately **100x local coverage** over the MUC1 VNTR region is recommended for reliable variant detection. Whole-genome sequencing at 30x average depth typically yields sufficient VNTR coverage. Whole-exome sequencing coverage varies according to capture probe efficiency across the GC-rich MUC1 VNTR.

## BAM versus FASTQ: which is faster?

**BAM or CRAM input is significantly faster.** With aligned input, VNtyper extracts only reads covering the MUC1 locus and unmapped read pools. FASTQ input requires whole-sample read alignment before downstream processing, substantially increasing runtime.

## Do I need adVNTR?

No. adVNTR is an **optional** module providing independent, alignment-based validation of Kestrel calls. Enabling it adds approximately 9 minutes per sample. It is valuable for orthogonal confirmation in research workflows, but is not required for primary genotyping.

## What does Low_Precision mean?

A call classified as **Low_Precision** indicates that Kestrel detected a variant, but its depth score falls between the configured low and high confidence thresholds in `kestrel_config.json`. These calls represent borderline read support and benefit from visual inspection in the HTML report or orthogonal validation. See [Scoring & Confidence](../pipeline/scoring-and-confidence.md) for threshold specifications.

## Can I use GRCh38?

Yes. VNtyper 2 supports standard UCSC, NCBI, and Ensembl naming conventions via the `--reference-assembly` argument:

- `hg19` or `GRCh37` (including `hg19_ensembl` and `hg19_ncbi`)
- `hg38` or `GRCh38` (including `hg38_ensembl` and `hg38_ncbi`)

Example:

```bash
vntyper pipeline --bam inputs/sample.bam --reference-assembly hg38 -o results/sample/
```

## Docker versus local installation?

**Docker** bundles all system dependencies (BWA, samtools, fastp, Java 11, and the Kestrel JAR) within a self-contained image. It provides the most reproducible execution environment.

A **local installation** avoids containerization overhead, but requires that BWA, samtools, fastp, and Java 11 reside on your `PATH`.

## How do I interpret the HTML report?

The interactive HTML report displays an executive screening summary, variant calls with confidence tiers, quality control metrics, and an embedded IGV browser for direct read inspection. Refer to [Output Files](../user-guide/output-files.md) for a comprehensive field-by-field breakdown.

## Does SHARK support BAM input?

No. SHARK operates directly on paired-end FASTQ reads to extract MUC1-matching reads prior to alignment. If invoked with `--bam` or `--cram`, the pipeline logs a warning and exits. When using `--extra-modules shark`, provide paired FASTQ files:

```bash
vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
    --extra-modules shark -o output/
```

## How do I process multi-sample cohorts?

Execute individual samples with `vntyper pipeline`, then aggregate the sample output directories with `vntyper cohort`. This generates cohort-wide summary reports, category breakdowns, and allele frequency tables. See [Cohort Analysis](../user-guide/cohort-analysis.md).

## Where can I get help?

Open an issue on [GitHub](https://github.com/hassansaei/VNtyper/issues). Search open and resolved issues before creating new tickets.
