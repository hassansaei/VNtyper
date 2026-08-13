# report

Generate an HTML summary report and visualizations from pipeline output data.

## Synopsis

```
vntyper [global-options] report
    -o <dir>
    [--input-dir <dir>]
    [--report-file <name>]
    [--bed-file <file>] [--bam-file <file>] [--reference-fasta <file>]
    [--flanking <int>]
    [-s <name>]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output-dir` | path | (required) | Directory containing pipeline results |
| `--input-dir` | path | — | If provided, search this directory (and its subdirectories) for standard pipeline output filenames |
| `--report-file` | string | `summary_report.html` | Name of the output report file |
| `--bed-file` | path | — | Path to the BED file for IGV reports |
| `--bam-file` | path | — | Path to the BAM file for IGV reports |
| `--reference-fasta` | path | (from config) | Path to the reference FASTA file for IGV reports. Falls back to `muc1_reference_vntr` in config if not provided |
| `--flanking` | int | `50` | Flanking region size for IGV reports |
| `-s, --sample-name` | string | (derived) | What the report calls its sample, in the title, the heading and the header block |

## Naming the sample

The report's `<title>`, `<h1>` and header block all carry the sample name, so two
reports are distinguishable in two browser tabs and a printed one says whose it is.

`--sample-name` wins whenever it is given. Otherwise the name is derived from the
input file names the run recorded in `pipeline_summary.json`: one recognised
compound extension (`.fastq.gz`, `.fq.gz`, `.bam`, `.cram`, …) is stripped, then a
single **trailing** `_R1` or `_R2`. So `example_b178_hg19_subset_R1.fastq.gz`
becomes `example_b178_hg19_subset` and `S1.lane3.L001_R1.fastq.gz` becomes
`S1.lane3.L001`.

Anything the rule does not recognise is printed verbatim rather than guessed at —
an unrecognised extension, a mate marker that is not at the end
(`S1_S1_L001_R1_001.fastq.gz` stays `S1_S1_L001_R1_001`), or two mates whose stems
disagree. A summary naming no input files at all is reported as `unnamed sample`.
Pass `--sample-name` when you want something else.

## Provenance

The report states the reference assembly the run asked for, the assembly its BAM
header declared, the region it actually analysed and the summary's schema version.
Each is read from `pipeline_summary.json` and **only** from there: a run that did
not record one has that row read `not recorded by this run`, never a value taken
from the configuration, which would mislabel any `--reference-assembly` override.

A report also carries two timestamps — when the pipeline ran, and when this file
was rendered. Re-running `vntyper report` over a finished run changes only the
second.

## Auto-Discovery

When `--input-dir` is provided and `--bam-file` or `--bed-file` are not specified, the report command will attempt to auto-discover standard pipeline output files:

- **BAM file:** `<input-dir>/kestrel/output.bam`
- **BED file:** `<input-dir>/kestrel/output.bed`

## Examples

Generate a report from pipeline output:

```bash
vntyper report -o results/ --input-dir results/
```

Generate a report with custom IGV settings:

```bash
vntyper report -o results/ \
    --bam-file results/kestrel/output.bam \
    --bed-file results/kestrel/output.bed \
    --reference-fasta ref/muc1.fa \
    --flanking 100
```

Generate a report with a custom filename:

```bash
vntyper report -o results/ --input-dir results/ --report-file my_report.html
```
