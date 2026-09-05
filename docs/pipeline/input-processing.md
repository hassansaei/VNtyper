# Input Processing

VNtyper accepts three input formats: aligned BAM files, aligned CRAM files, or paired-end FASTQ files. Each format follows a dedicated preprocessing path before genotyping.

## BAM/CRAM Region Extraction

For aligned inputs, the pipeline extracts reads overlapping the MUC1 VNTR locus using `samtools view` with coordinates matched to the detected reference assembly.

| Assembly | Style | BAM Extraction Region | VNTR Core Region |
|----------|-------|-----------------------|------------------|
| hg19     | UCSC  | `chr1:155158000-155163000` | `chr1:155160500-155162000` |
| hg38     | UCSC  | `chr1:155184000-155194000` | `chr1:155188000-155192500` |
| GRCh37   | NCBI  | `NC_000001.10:155158000-155163000` | `NC_000001.10:155160500-155162000` |
| GRCh38   | NCBI  | `NC_000001.11:155184000-155194000` | `NC_000001.11:155188000-155192500` |

The extraction region spans a 5 to 10 kb window to capture flanking reads. The core region defines the target interval for coverage calculations and variant evaluation.

Extraction proceeds in three steps:

1. **Region slicing**: `samtools view -P -b` extracts reads overlapping the MUC1 coordinates.
2. **Unmapped read recovery**: In standard mode (non-fast mode), unmapped and partially mapped reads are extracted and merged with sliced reads. This preserves reads whose variants prevented standard alignment.
3. **FASTQ conversion**: The merged BAM is name-sorted and converted to paired FASTQ files via `samtools fastq`.

!!! info "Why recover unmapped reads?"
    MUC1 VNTR insertions often disrupt read alignment against the reference. Recovering unmapped mates and orphaned reads prevents variant dropout before k-mer genotyping.

### Atomic Artifact Publication

All intermediate and final alignment files (`.bam`) and index files (`.bai`) are created atomically. Subprocesses write to deterministic temporary files ending in `.partial` (for example, `output_sliced.bam.partial`). Once verified for zero exit status and non-zero file size, files are installed to permanent destinations via atomic POSIX `os.replace`.

A leftover `.partial` file indicates that an upstream process was terminated prematurely (for example, via `SIGKILL` or out-of-memory limits). Public BAM and BAI files are guaranteed complete.

## FASTQ Passthrough

When supplied directly with paired FASTQ inputs, the pipeline optionally applies SHARK filtering, runs fastp QC, aligns reads via BWA-MEM, and extracts the target MUC1 region from the resulting BAM.

## fastp Quality Control

FASTQ quality control runs via [fastp](https://github.com/OpenGene/fastp) with configurable parameters:

| Parameter | Config Key | Description |
|-----------|-----------|-------------|
| Adapter trimming | `disable_adapter_trimming` | Disabled by default for pre-trimmed data |
| Deduplication | `deduplication` | Optical and PCR duplicate removal |
| Quality threshold | `qualified_quality_phred` | Minimum base quality score |
| Length filter | `length_required` | Minimum read length after trimming |
| Compression | `compression_level` | Output gzip compression level |

fastp generates `output.json`, providing metrics for downstream reporting including Q20/Q30 rates, duplication rates, and passed-filter read fractions.

## Coverage Calculation

The pipeline measures per-base depth across the VNTR region using `samtools depth -a`. Summary metrics include:

- **Mean and median coverage**: Primary indicators of sequencing depth.
- **Standard deviation, min, max**: Metrics of coverage uniformity.
- **Region length**: Total VNTR interval length in base pairs.
- **Uncovered bases and percent uncovered**: Count and fraction of positions with zero read depth.
- **Coverage QC**: Categorical `PASS` or `FAIL` evaluated against configured thresholds.

Every statistic is computed across the entire region, not merely over covered bases. The `-a` flag forces `samtools depth` to emit every genomic position, including zero-depth sites. Truncated depth files are zero-padded to the full interval length before computing metrics. A sample with 30x depth across 10% of the VNTR yields a mean coverage of 3x, and minimum coverage is 0 whenever any base lacks coverage.

### Comparing Coverage Between Assemblies

Mean VNTR coverage cannot be compared directly between GRCh37 and GRCh38. For the same sequencing library, published window mean depth is approximately 2.7-fold higher on GRCh37 (2.498 to 2.888 across paired fixtures). Across the repeat array alone, depth is roughly 4.3-fold higher on GRCh37 (4.257 to 4.339).

This difference stems from reference structure: GRCh37 represents roughly 13.5 copies of the 60 bp unit, whereas GRCh38 contains approximately 58 copies. Concentrating identical read pools into one-fourth the reference sequence quadruples apparent depth.

Three metrics facilitate cross-assembly comparisons:

- **`vntr_flank_mean_depth`**: Mean read depth over unique flanking sequence (`vntr_flank_bases` per side). This metric is directly comparable between assemblies because unique sequence lacks copy-number variation. Values depend on the overlap-counting policy recorded in `depth_counting_policy`.
- **`vntr_array_depth_sum`** and **`vntr_array_depth_sum_per_unit_length`**: Summed depth across the array, and that sum divided by `depth_sum_reference_length`. These values are comparable between assemblies only under the counting policy recorded in `depth_counting_policy`. Under default `samtools depth` settings (counting all primary alignments), the two builds agree within 1%. Requiring `MAPQ >= 1` shifts the GRCh37-to-GRCh38 ratio to 0.31-0.55 because reads mapping to GRCh37's collapsed repeat array have low mapping quality.

When no array geometry applies to a run, these metrics are recorded as `NA` rather than `0`. Zero indicates an evaluated region with zero reads, whereas `NA` indicates unmeasured array metrics.

Metrics are saved to `coverage_summary.tsv`. A mean coverage below the configured cutoff (default: 100x) or an uncovered fraction exceeding its cutoff (default: 50%) triggers a `FAIL` verdict in coverage QC.

## BAM Header Analysis

For BAM and CRAM inputs, the pipeline inspects `@SQ` and `@PG` header records:

- **Reference assembly**: Detected via keyword matching (`hg19`, `GRCh38`) and contig length verification.
- **Alignment pipeline**: Identified as BWA, Dragen, CLC, or Unknown.

!!! warning "Dragen and CLC aligners"
    The Dragen pipeline frequently misaligns or clips reads in the repetitive MUC1 VNTR region. VNtyper recommends BWA-aligned inputs. The report issues an explicit warning whenever a non-BWA aligner is detected.
