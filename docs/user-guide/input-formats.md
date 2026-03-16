# Input Formats

VNtyper 2 accepts BAM, CRAM, or paired-end FASTQ files as input. Provide exactly one input type per run.

=== "BAM"

    ```bash
    vntyper pipeline --bam sample.bam -o results/
    ```

    Requirements:

    - Sorted and indexed (`.bam.bai` or `.bai` must exist alongside the BAM)
    - Aligned to a [supported reference assembly](reference-assemblies.md)
    - VNtyper 2 validates the file with `samtools quickcheck` before processing

=== "CRAM"

    ```bash
    vntyper pipeline --cram sample.cram -o results/
    ```

    Requirements:

    - Sorted and indexed (`.cram.crai` must exist)
    - Aligned to a supported reference assembly
    - The original reference FASTA must be accessible (CRAM files are reference-dependent)
    - Validated with `samtools quickcheck`

=== "FASTQ"

    ```bash
    vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz -o results/
    ```

    Requirements:

    - Paired-end reads: both `--fastq1` and `--fastq2` are required
    - Gzipped (`.fastq.gz`) or uncompressed (`.fastq`) accepted
    - VNtyper 2 validates FASTQ format (checks first 4 lines for correct structure)
    - Reads are processed through fastp for quality control, then aligned with BWA

!!! note "SHARK module requires FASTQ input"
    The SHARK filtering module (`--extra-modules shark`) only works with FASTQ input.
    Attempting to use SHARK with BAM or CRAM input will cause VNtyper 2 to exit with an error.

## Input Validation

VNtyper 2 performs automatic validation before starting the pipeline:

| Input Type | Validation Method |
|------------|-------------------|
| BAM / CRAM | `samtools quickcheck` -- verifies file integrity and EOF marker |
| FASTQ      | Format check -- verifies the first 4 lines follow FASTQ structure |

If validation fails, the pipeline exits immediately with a descriptive error message.
