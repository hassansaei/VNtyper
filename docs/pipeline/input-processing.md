# Input Processing

VNtyper 2 accepts three input formats: aligned BAM files, aligned CRAM files, or paired-end FASTQ files. Each format follows a distinct preprocessing path before entering the genotyping stage.

## BAM/CRAM Region Extraction

For aligned input, the pipeline extracts reads overlapping the MUC1 VNTR locus using `samtools view` with region coordinates appropriate for the detected reference assembly.

| Assembly | Style | BAM Extraction Region | VNTR Core Region |
|----------|-------|-----------------------|------------------|
| hg19     | UCSC  | `chr1:155158000-155163000` | `chr1:155160500-155162000` |
| hg38     | UCSC  | `chr1:155184000-155194000` | `chr1:155188000-155192500` |
| GRCh37   | NCBI  | `NC_000001.10:155158000-155163000` | `NC_000001.10:155160500-155162000` |
| GRCh38   | NCBI  | `NC_000001.11:155184000-155194000` | `NC_000001.11:155188000-155192500` |

The BAM extraction region is a wider window (~5-10 kb) to capture flanking reads. The VNTR core region is the narrower interval used for coverage calculation and variant analysis.

The extraction process:

1. **Region slicing** -- `samtools view -P -b` extracts reads from the MUC1 region (using a BED file or coordinate string)
2. **Unmapped read recovery** -- In normal mode (not fast mode), unmapped and partially mapped reads are extracted and merged with the sliced BAM. This captures reads that may carry VNTR variants but failed to align
3. **FASTQ conversion** -- The merged BAM is name-sorted and converted to paired FASTQ using `samtools fastq`

!!! info "Why recover unmapped reads?"
    MUC1 VNTR variants (especially large insertions) can cause reads to fail alignment against the reference. Recovering unmapped reads ensures these variant-carrying reads are not lost before k-mer-based genotyping.

## FASTQ Passthrough

When FASTQ files are provided directly, the pipeline optionally applies SHARK filtering first (if enabled), then runs fastp QC, aligns reads with BWA, and re-extracts the MUC1 region from the resulting BAM.

## fastp Quality Control

FASTQ quality control uses [fastp](https://github.com/OpenGene/fastp) with the following configurable parameters:

| Parameter | Config Key | Description |
|-----------|-----------|-------------|
| Adapter trimming | `disable_adapter_trimming` | Disabled by default for pre-trimmed data |
| Deduplication | `deduplication` | Optical/PCR duplicate removal |
| Quality threshold | `qualified_quality_phred` | Minimum base quality score |
| Length filter | `length_required` | Minimum read length after trimming |
| Compression | `compression_level` | Output gzip compression level |

fastp produces a JSON report (`output.json`) with metrics used in the final HTML report, including Q20/Q30 rates, duplication rates, and passed-filter read counts.

## Coverage Calculation

After alignment, the pipeline computes per-base coverage over the VNTR region using `samtools depth`. The summary statistics include:

- **Mean and median coverage** -- primary indicators of sequencing depth
- **Standard deviation, min, max** -- coverage uniformity
- **Region length** -- total VNTR span in base pairs
- **Uncovered bases and percent uncovered** -- fraction of the VNTR with zero coverage

These metrics are written to `coverage_summary.tsv` and are used in the final report for quality assessment. A mean VNTR coverage below the configured threshold (default: 100x) triggers a warning in the report.

## BAM Header Analysis

For BAM/CRAM input, the pipeline parses the header (`@SQ` and `@PG` lines) to extract:

- **Reference assembly** -- detected via both text matching (keywords like `hg19`, `GRCh38`) and contig-length matching against known assemblies
- **Alignment pipeline** -- identified as BWA, Dragen, CLC, or Unknown

!!! warning "Dragen and CLC aligners"
    The Dragen pipeline has known issues aligning reads in the MUC1 VNTR region. CLC has not been extensively tested. VNtyper 2 recommends BWA-aligned input for optimal results. When a non-BWA aligner is detected, the report includes a warning.
