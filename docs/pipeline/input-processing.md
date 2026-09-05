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

### Atomic Artifact Publication

All intermediate and final alignment files (`.bam`) and their corresponding index files (`.bai`) are created atomically. Commands write to temporary, deterministic partial files ending in `.partial` (e.g. `output_sliced.bam.partial`). Upon successful completion and size verification, files are atomically installed to their permanent names using POSIX `os.replace`.

A leftover `.partial` file in the output directory indicates that an upstream process was killed (such as via `SIGKILL` or an out-of-memory event) or encountered an unhandled crash before completing. Public BAM and BAI files are guaranteed to be complete and fully formed.

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

After alignment, the pipeline computes per-base coverage over the VNTR region using `samtools depth -a`. The summary statistics include:

- **Mean and median coverage** -- primary indicators of sequencing depth
- **Standard deviation, min, max** -- coverage uniformity
- **Region length** -- total VNTR span in base pairs
- **Uncovered bases and percent uncovered** -- fraction of the VNTR with zero coverage
- **Coverage QC** -- `PASS` or `FAIL` against the configured thresholds

**Every statistic is computed over the region, not over the positions that carry reads.** The `-a` flag makes `samtools depth` emit a row for every position in the region, zero-depth ones included, and a depth table written without it is padded back to the region with zeros before anything is computed. A sample covered at 30x across 10% of the VNTR therefore reports a mean of 3, and `min` is `0` wherever any position is uncovered. Before VNtyper 2.0.8 the mean divided by the number of covered positions and reported 30 -- the depth where there happened to be reads, which is not a property of the region and was systematically too high exactly where coverage was patchy.

### Comparing coverage between assemblies

**Mean coverage over the VNTR region is not comparable between GRCh37 and GRCh38.** Measured on the same sample aligned to both, the published window mean runs about **2.7x higher on GRCh37** (2.498-2.888 across seven paired fixtures). Over the repeat array alone the gap is larger still, about **4.3x** (4.257-4.339), and the window figure is smaller only because the window also carries flanking sequence whose depth does not scale with repeat count.

That is not a defect in either number. The two assemblies represent different amounts of the repeat array -- roughly 13.5 60-bp units in GRCh37 against roughly 58 in GRCh38 -- and a sample's reads pile onto whichever copies the reference happens to contain. Concentrating the same reads onto a quarter of the sequence quadruples the depth. Harmonising the configured windows would not change this, because the difference is in the reference, not in the window.

Three columns exist so a reader can compare something, and they do not carry equal weight:

- **`vntr_flank_mean_depth`** -- mean depth over unique sequence flanking the array, `vntr_flank_bases` per side. This is the figure to compare **between assemblies**: unique sequence has no repeat copies to distribute reads between, so it measures the sample rather than the reference, and the cross-build ratio is unmoved by mapping-quality or overlap filtering. It is not a policy-free absolute -- counting overlapping mates once changes the value itself by 24-40% -- so compare it across runs made the same way, which `depth_counting_policy` records.
- **`vntr_array_depth_sum`** and **`vntr_array_depth_sum_per_unit_length`** -- summed depth over the array, and that sum over the fixed `depth_sum_reference_length`. Comparable between assemblies **only** under the counting policy named in `depth_counting_policy`, and the qualification is not decorative: under `samtools depth` defaults the two builds agree to within 1%, and requiring `MAPQ >= 1` moves them to a ratio of 0.31-0.55. The equality holds because every primary alignment is counted whatever its mapping quality, and GRCh37's array is a collapsed near-consensus core whose reads are genuinely ambiguous.

All of these are recorded as `NA` when no array geometry applies to the run -- a configuration without `vntr_array_coords`, or a coverage region that is not the configured VNTR window. `NA` rather than `0`: zero would say the region was examined and no reads were found. Note that `--custom-regions` and `--bed-file` change the *extraction* target while coverage is still resolved to the configured VNTR window, so those runs do report the figures.

The array bounds are heuristic rather than an authoritative annotation, and both were derived by one procedure applied identically to both assemblies. That symmetry is what matters -- moving both outward by 100 bp changes the comparison by about 1%, while moving one build's bound alone by 140 bp changes it by 17.5%.

These metrics are written to `coverage_summary.tsv` and are used in the final report for quality assessment. A mean VNTR coverage below the configured threshold (default: 100x), or an uncovered fraction above its own (default: 50%), fails the report's coverage QC.

## BAM Header Analysis

For BAM/CRAM input, the pipeline parses the header (`@SQ` and `@PG` lines) to extract:

- **Reference assembly** -- detected via both text matching (keywords like `hg19`, `GRCh38`) and contig-length matching against known assemblies
- **Alignment pipeline** -- identified as BWA, Dragen, CLC, or Unknown

!!! warning "Dragen and CLC aligners"
    The Dragen pipeline has known issues aligning reads in the MUC1 VNTR region. CLC has not been extensively tested. VNtyper 2 recommends BWA-aligned input for optimal results. When a non-BWA aligner is detected, the report includes a warning.
