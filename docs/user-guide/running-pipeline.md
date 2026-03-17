# Running the Pipeline

## Minimal Run

```bash
vntyper pipeline --bam sample.bam -o results/
```

This uses default settings: hg19 assembly, 4 threads, Kestrel genotyping only.

## Common Options

```bash
vntyper pipeline --bam sample.bam -o results/ \
    --threads 8 \
    --fast-mode \
    --reference-assembly hg38
```

| Option | Effect |
|--------|--------|
| `--threads N` | Number of threads (default: 4) |
| `--fast-mode` | Skip filtering for unmapped and partially mapped reads (faster) |
| `--reference-assembly` | Assembly used for alignment (default: hg19). See [Reference Assemblies](reference-assemblies.md) |
| `--sample-name` | Label for results. Defaults to input filename stem |
| `--output-name` | Base name for intermediate output files (default: `processed`) |

## With adVNTR Module

```bash
vntyper pipeline --bam sample.bam -o results/ \
    --extra-modules advntr
```

adVNTR provides independent VNTR genotyping. Results are cross-matched with Kestrel calls. Use `--advntr-max-coverage 300` for faster adVNTR runs on high-coverage samples.

## With SHARK Module (FASTQ Only)

```bash
vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz -o results/ \
    --extra-modules shark
```

!!! warning
    SHARK only works with FASTQ input. Using `--extra-modules shark` with `--bam` or `--cram` will cause an error.

## Custom Regions

Override the default MUC1 VNTR region coordinates:

=== "Inline regions"

    ```bash
    vntyper pipeline --bam sample.bam -o results/ \
        --custom-regions chr1:155160500-155162000
    ```

    Multiple regions are comma-separated: `chr1:1000-2000,chr1:3000-4000`

=== "BED file"

    ```bash
    vntyper pipeline --bam sample.bam -o results/ \
        --bed-file regions.bed
    ```

`--custom-regions` and `--bed-file` are mutually exclusive.

## Archiving Results

```bash
vntyper pipeline --bam sample.bam -o results/ \
    --archive-results --archive-format tar.gz
```

Creates a compressed archive of the output directory after the pipeline completes. Supported formats: `zip` (default) and `tar.gz`.

## Additional Summary Formats

```bash
vntyper pipeline --bam sample.bam -o results/ \
    --summary-formats csv,tsv
```

Generates `pipeline_summary.csv` and/or `pipeline_summary.tsv` alongside the default JSON summary.

## Intermediate Files

| Option | Effect |
|--------|--------|
| `--keep-intermediates` | Retain BAM slices, temporary FASTQ files, etc. |
| `--delete-intermediates` | Delete intermediate files after processing (overrides `--keep-intermediates`) |

## Custom Configuration

```bash
vntyper --config-path /path/to/custom/config.json pipeline --bam sample.bam -o results/
```

Note that `--config-path` is a global option and must appear before the subcommand. See [Configuration](configuration.md) for details.

## Logging

```bash
vntyper -l DEBUG pipeline --bam sample.bam -o results/
```

Log levels: `DEBUG`, `INFO` (default), `WARNING`, `ERROR`, `CRITICAL`. The pipeline log is automatically written to `<output-dir>/pipeline.log`. Override with `-f /path/to/logfile`.
