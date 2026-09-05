# Input Formats

VNtyper 2 accepts BAM, CRAM, or paired-end FASTQ files as input. Provide exactly one input type per run.

=== "BAM"

    ```bash
    vntyper pipeline --bam inputs/sample.bam -o results/sample/
    ```

    Requirements:

    - Sorted and indexed (`.bam.bai` or `.bai` must exist alongside the BAM)
    - Aligned to a [supported reference assembly](reference-assemblies.md)
    - VNtyper 2 validates file integrity with `samtools quickcheck` before processing

=== "CRAM"

    ```bash
    vntyper pipeline --cram inputs/sample.cram -o results/sample/
    ```

    Requirements:

    - Sorted and indexed (`.cram.crai` or `.crai` must exist alongside the CRAM)
    - Aligned to a supported reference assembly
    - Reference FASTA must be accessible (CRAM compression depends on reference sequence)
    - A local reference named in the CRAM header must resolve inside the directory containing
      the CRAM. For a reference elsewhere, specify `--reference-fasta` or configure it in `config.json`.
    - Validated with `samtools quickcheck` before processing

=== "FASTQ"

    ```bash
    vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz -o results/
    ```

    Requirements:

    - Paired-end reads: both `--fastq1` and `--fastq2` are required
    - Gzipped (`.fastq.gz`) or uncompressed (`.fastq`) accepted
    - Validates format by verifying first 4 lines follow standard FASTQ structure
    - Reads are quality-controlled via fastp, then aligned with BWA

!!! warning "SHARK requires FASTQ input"
    The SHARK read-filtering module (`--extra-modules shark`) requires paired-end FASTQ input.
    If specified with BAM or CRAM input, the pipeline logs an error and exits immediately (`sys.exit(1)`).

!!! important "Keep alignment inputs and outputs in separate trees"
    For BAM or CRAM input, the output root directory must stay outside the directory containing the
    alignment file and all descendants of that directory. A BAM at `sample.bam` with output at
    `results/` shares the current working directory as a common root, causing the run to be rejected.
    Use distinct directory paths such as `inputs/sample.bam` and `results/sample/`.

## Input Validation

VNtyper 2 performs automated validation before executing any pipeline stage:

| Input Type | Validation Method |
|------------|-------------------|
| BAM / CRAM | `samtools quickcheck`: verifies file integrity, BAM/CRAM headers, and EOF markers |
| FASTQ      | Format verification: validates record structure across first 4 lines |

If validation fails, the pipeline aborts immediately with an informative diagnostic error.
