# pipeline

Run the full VNtyper 2 genotyping pipeline.

## Synopsis

```
vntyper [global-options] pipeline
    (--bam <file> | --cram <file> | --fastq1 <file> --fastq2 <file>)
    [-o <dir>] [-n <name>] [-s <name>]
    [--reference-assembly <assembly>] [--custom-regions <regions> | --bed-file <file>]
    [--threads <int>] [--fast-mode] [--keep-intermediates] [--delete-intermediates]
    [--extra-modules <module> ...]
    [--advntr-max-coverage <int>]
    [--archive-results] [--archive-format <format>]
    [--summary-formats <formats>]
```

## Input Options

Provide exactly one input type: BAM, CRAM, or paired FASTQ files.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--bam` | path | — | Path to the BAM file |
| `--cram` | path | — | Path to the CRAM file |
| `--fastq1` | path | — | Path to the first FASTQ file (paired-end) |
| `--fastq2` | path | — | Path to the second FASTQ file (paired-end) |

When using FASTQ input, both `--fastq1` and `--fastq2` are required.

## Output Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output-dir` | path | `out` | Output directory for the results |
| `-n, --output-name` | string | `processed` | Base name for the output files |
| `-s, --sample-name` | string | (from input filename) | Sample name for labeling results. If not provided, defaults to the input BAM or FASTQ filename stem |

## Reference & Region Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--reference-assembly` | choice | `hg19` | Reference assembly for BAM/CRAM alignment. Options: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ncbi`, `hg38_ncbi`, `hg19_ensembl`, `hg38_ensembl` |
| `--custom-regions` | string | — | Custom regions for MUC1 analysis as comma-separated values (e.g., `chr1:1000-2000,chr2:3000-4000`) |
| `--bed-file` | path | — | Path to a BED file specifying regions for MUC1 analysis |

`--custom-regions` and `--bed-file` are mutually exclusive.

## Processing Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--threads` | int | `4` | Number of threads to use |
| `--fast-mode` | flag | off | Enable fast mode (skips filtering for unmapped and partially mapped reads) |
| `--keep-intermediates` | flag | off | Keep intermediate files (e.g., BAM slices, temporary files) |
| `--delete-intermediates` | flag | off | Delete intermediate files after processing. Overrides `--keep-intermediates` |

## Optional Modules

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--extra-modules` | string | `[]` | Optional extra modules to include (e.g., `advntr`, `shark`). Can be repeated multiple times |
| `--advntr-max-coverage` | int | — | Max coverage (e.g., 300) for quick adVNTR mode. Only applies when `advntr` is in `--extra-modules` |

The `shark` module is not supported in BAM/CRAM mode; use FASTQ mode or remove the shark flag.

## Archive & Summary Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--archive-results` | flag | off | Create an archive of the results folder after pipeline completion |
| `--archive-format` | choice | `zip` | Format of the archive: `zip` or `tar.gz` |
| `--summary-formats` | string | `""` | Comma-separated list of additional summary output formats to generate (supported: `csv`, `tsv`). JSON is always generated |

## Examples

Run the pipeline with a BAM file using default settings:

```bash
vntyper pipeline --bam sample.bam -o results/
```

Run with paired FASTQ files and hg38 reference:

```bash
vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
    --reference-assembly hg38 -o results/ -s my_sample
```

Enable fast mode with multiple threads and archive results:

```bash
vntyper pipeline --bam sample.bam -o results/ --threads 8 \
    --fast-mode --archive-results --archive-format tar.gz
```

Run with the adVNTR module and coverage cap:

```bash
vntyper pipeline --bam sample.bam -o results/ \
    --extra-modules advntr --advntr-max-coverage 300
```

Generate additional summary formats and clean up intermediate files:

```bash
vntyper pipeline --bam sample.bam -o results/ \
    --summary-formats csv,tsv --delete-intermediates
```
