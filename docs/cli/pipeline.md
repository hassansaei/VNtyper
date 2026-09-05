# pipeline

Run the full VNtyper 2 genotyping pipeline.

## Synopsis

```
vntyper [global-options] pipeline
    (--bam <file> | --cram <file> | --fastq1 <file> --fastq2 <file>)
    [--reference-fasta <file>]
    [-o <dir>] [-n <name>] [-s <name>]
    [--reference-assembly <assembly>] [--custom-regions <regions> | --bed-file <file>]
    [--threads <int>] [--fast-mode] [--keep-intermediates] [--delete-intermediates]
    [--decision-profile <complete-profile.json>]
    [--extra-modules <module> ...]
    [--advntr-max-coverage <int>]
    [--archive-results] [--archive-format <format>]
    [--summary-formats <formats>]
    [--report-igv {embedded,sidecar,off}]
```

## Input Options

Provide exactly one input format: BAM, CRAM, or paired FASTQ files.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--bam` | path | None | Path to input BAM file |
| `--cram` | path | None | Path to input CRAM file |
| `--reference-fasta` | path | None | Path to reference FASTA for CRAM decoding |
| `--fastq1` | path | None | Path to first FASTQ mate file (paired-end) |
| `--fastq2` | path | None | Path to second FASTQ mate file (paired-end) |

When supplying FASTQ input, both `--fastq1` and `--fastq2` are required.

## Output Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output-dir` | path | `out` | Destination directory for results. For BAM/CRAM input, this must reside outside the directory containing the alignment file |
| `-n, --output-name` | string | `output` | Base name for output files. Fixed at `output`: any other value is rejected. Stage components rely on this literal naming; use `--output-dir` to separate distinct analyses |
| `-s, --sample-name` | string | (from filename) | Sample identifier for reporting. Defaults to input BAM or FASTQ file stem |

## Reference & Region Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--reference-assembly` | choice | `hg19` | Reference genome assembly for alignment. Options: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ncbi`, `hg38_ncbi`, `hg19_ensembl`, `hg38_ensembl`. For BAM/CRAM input, this is validated against alignment headers (see [Reference Assemblies](../user-guide/reference-assemblies.md#the-declared-assembly-check)) |
| `--custom-regions` | string | None | Custom MUC1 target coordinates as comma-separated 1-based inclusive intervals (e.g. `chr1:1000-2000,chr2:3000-4000`) |
| `--bed-file` | path | None | Path to BED file specifying intervals using 0-based half-open coordinates |

`--custom-regions` and `--bed-file` are mutually exclusive. Configured extraction windows retain literal BED coordinates. User-supplied inline strings convert to 1-based coordinates.

## Processing Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--threads` | int | `4` | Execution thread count for BWA, samtools, fastp, and adVNTR |
| `--fast-mode` | flag | off | Enable fast processing mode (bypasses filtering for unmapped and partially mapped reads) |
| `--keep-intermediates` | flag | off | Compatibility flag: intermediate files are preserved by default |
| `--delete-intermediates` | flag | off | Remove temporary BAM slices and staging files upon completion (takes precedence over `--keep-intermediates`) |
| `--resume` | flag | off | Resume execution from a prior run in `--output-dir`, reusing completed stages |
| `--decision-profile` | path | packaged profile | Select one explicit-custom or generated decision profile. Overlays or modified fixed-safety cutoffs trigger immediate abort |

Omitting `--decision-profile` applies the verified packaged profile. Runtime paths, reference indices, and coverage parameters are controlled via configuration or code rather than decision profile leaves.

## Optional Modules

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--extra-modules` | string | `[]` | Optional modules to execute: `advntr`, `shark`. Repeat the flag or supply comma-separated values (`--extra-modules advntr,shark`) |
| `--advntr-max-coverage` | int | None | Coverage cap (e.g. 300) for adVNTR quick mode. Only valid when `advntr` is active |

The `shark` module requires FASTQ input and is not supported in BAM/CRAM mode.

## Archive & Summary Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--archive-results` | flag | off | Package output folder into a compressed archive upon completion |
| `--archive-format` | choice | `zip` | Archive compression format: `zip` or `tar.gz` |
| `--summary-formats` | string | `""` | Comma-separated summary formats (`csv`, `tsv`). Writes `pipeline_summary.<fmt>` and `pipeline_summary_rows.<fmt>`. JSON generates automatically |
| `--report-igv` | choice | `embedded` | Alignment viewer deployment mode: `embedded`, `sidecar`, or `off` |

## Report Options

| Value | Behavior | Size |
|-------|----------|------|
| `embedded` (default) | Single self-contained HTML file embedding gzipped igv.js. Operates offline without external network assets | ~560 KB with alignments |
| `sidecar` | Summary report omits viewer and references `igv_report.html` generated in the same directory | ~75 KB, plus sidecar |
| `off` | Omits alignment viewer generation entirely | ~75 KB |

Reports never query remote hosts. All styles, scripts, and fonts resolve locally.

## Examples

Run the pipeline using BAM input:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/
```

Run with paired FASTQ files and hg38 reference:

```bash
vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
    --reference-assembly hg38 -o results/ -s my_sample
```

Process CRAM input with explicit reference:

```bash
vntyper pipeline --cram inputs/sample.cram \
    --reference-fasta /data/ref/GRCh38.fa \
    --reference-assembly hg38 -o results/sample/
```

Enable fast mode with 8 threads and archive results:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ --threads 8 \
    --fast-mode --archive-results --archive-format tar.gz
```

Include adVNTR module with coverage limit:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --extra-modules advntr --advntr-max-coverage 300
```

Apply an explicit decision profile:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --decision-profile reviewed-complete-profile.json
```

Generate tabular summary formats and remove intermediate files:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --summary-formats csv,tsv --delete-intermediates
```
