# Running the Pipeline

## Minimal Run

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/
```

This uses default settings: hg19 assembly, 4 threads, Kestrel genotyping engine.

!!! important "Keep alignment inputs and outputs in separate trees"
    For BAM or CRAM input, the output root cannot be the directory containing the
    alignment or any of its descendants. A BAM at `sample.bam` with `-o results/` is
    rejected because both resolve to the current directory root. Use separate trees
    such as `inputs/sample.bam` and `results/sample/`.

## Common Options

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --threads 8 \
    --fast-mode \
    --reference-assembly hg38
```

| Option | Effect |
|--------|--------|
| `--threads N` | Number of threads to use (default: 4) |
| `--fast-mode` | Skip filtering for unmapped and partially mapped reads (faster execution) |
| `--reference-assembly` | Assembly used for alignment (default: hg19). See [Reference Assemblies](reference-assemblies.md) |
| `--sample-name` | Label for results. Defaults to input filename stem |
| `--output-name` | Base name for intermediate files. Fixed at `output`; any other value is rejected. Use `--output-dir` to separate runs |

## With adVNTR Module

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --extra-modules advntr
```

adVNTR provides independent Profile-HMM VNTR genotyping. Results are cross-matched with Kestrel calls. Use `--advntr-max-coverage 300` for faster adVNTR processing on high-coverage sequencing runs.

## With SHARK Module (FASTQ Only)

```bash
vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz -o results/ \
    --extra-modules shark
```

!!! warning "SHARK requires FASTQ input"
    SHARK read filtering only functions with paired-end FASTQ input. Supplying `--extra-modules shark` with `--bam` or `--cram` causes the pipeline to abort with an error.

## Custom Regions

Override default MUC1 VNTR coordinates:

=== "Inline regions"

    ```bash
    vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
        --custom-regions chr1:155160500-155162000
    ```

    Multiple regions are comma-separated: `chr1:1000-2000,chr1:3000-4000`

=== "BED file"

    ```bash
    vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
        --bed-file regions.bed
    ```

`--custom-regions` and `--bed-file` are mutually exclusive.

## Archiving Results

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --archive-results --archive-format tar.gz
```

Creates a compressed archive of the output directory after pipeline completion. Supported formats: `zip` (default) and `tar.gz`.

## Additional Summary Formats

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --summary-formats csv,tsv
```

Each requested format writes two tabular files alongside `pipeline_summary.json`:
`pipeline_summary.csv` / `.tsv` (one row per pipeline step, with run provenance in
leading `run_*` columns), and `pipeline_summary_rows.csv` / `.tsv` (one row per step,
result row, and field, capturing full multi-row adVNTR and cross-match results).
Unsupported formats are ignored. See [Output Files](output-files.md) for column definitions.

## Intermediate Files

| Option | Effect |
|--------|--------|
| `--keep-intermediates` | Compatibility flag: intermediate files are kept by default, so this flag changes nothing. Use `--delete-intermediates` to remove them. |
| `--delete-intermediates` | Delete intermediate files after processing (takes precedence when `--keep-intermediates` is also supplied). |

## Custom Configuration

```bash
vntyper --config-path /path/to/custom/config.json pipeline \
    --bam inputs/sample.bam -o results/sample/
```

`--config-path` is a global option and must appear before the subcommand. See [Configuration](configuration.md) for details.

## Logging

```bash
vntyper -l DEBUG pipeline --bam inputs/sample.bam -o results/sample/
```

Global log levels: `DEBUG`, `INFO` (default), `WARNING`, `ERROR`, `CRITICAL`. Global options must precede the subcommand. Pipeline logs write automatically to `<output-dir>/pipeline.log`. Override with `-f /path/to/logfile`.

## Resuming Execution

When a pipeline run is interrupted or when re-running an analysis with new reporting options, pass `--resume` to reuse completed stage outputs:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ --resume
```

### Reusable Stages

Resumption evaluates checkpoint integrity from `pipeline_summary.json` and reuses completed outputs for:
- **Alignment and extraction**: BAM/CRAM to FASTQ extraction or FASTQ alignment.
- **Kestrel genotyping**: Reuses `kestrel_result.tsv`, `output.vcf`, `output.bam`, and `kestrel_pre_result.tsv`.
- **adVNTR genotyping**: Reuses `output_adVNTR_result.tsv`.

Lightweight stages (BAM header parsing, coverage calculation, cross-match comparison, nomenclature determination, and HTML report rendering) always recompute to ensure up-to-date reporting.

### Run Identity Invariants and Refusals

Resumption aborts with an error if run identity diverges from the prior summary checkpoint:
- Pipeline version changes
- Input file path changes
- Sample name changes
- Reference key or assembly changes
- Decision profile checksum changes

If an output artifact from a reusable stage is missing or fails MD5 verification, that stage and all dependent downstream stages are recomputed.
