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
