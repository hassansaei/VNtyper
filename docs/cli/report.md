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
    [--report-igv {embedded,sidecar,off}]
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
| `--report-igv` | choice | `embedded` | How the report carries its alignment browser. See [The alignment browser](#the-alignment-browser) |

## The alignment browser

The report is a file people archive, forward and reopen years later, so it carries
everything it needs. igv.js 3.0.2 is vendored inside VNtyper, stored gzipped, and
written into the report base64-encoded; the browser expands it with
`DecompressionStream`. **No mode fetches anything from the network** — not the library,
not a stylesheet, not a font, and not igv.js's own genome registry, which is switched
off explicitly.

| Value | What the report is | Size |
|-------|--------------------|------|
| `embedded` (default) | One file, and a complete alignment browser. | ~560 KB with alignments |
| `sidecar` | No library in the report; `igv_report.html` beside it is the alignment browser, and it is self-contained too (`create_report --standalone`). | ~75 KB, plus the sidecar |
| `off` | No alignment browser at all; `create_report` is not run. | ~75 KB |

A report for a run with no alignment data carries no library in any mode and is ~75 KB.

The report's Provenance block states the version and the **SHA-256 of the decompressed
library**, in full, so it can be checked against upstream. That digest is verified in
Python at render time, over the exact bytes the document later evaluates. It is
deliberately not re-checked in the browser: `crypto.subtle` is unavailable on `file://`
in several engines, and an archived report is opened from `file://` more often than
not, so a runtime check would be one that silently does not run on the machines it was
written for.

If a reader's browser predates `DecompressionStream` (Chrome 80, Safari 16.4, Firefox
113), the alignment panel says so in words and the variant tables stand unchanged. No
state of this report is blank space.

## Naming the sample

The report's `<title>`, `<h1>` and header block all carry the sample name, so two
reports are distinguishable in two browser tabs and a printed one says whose it is.

`--sample-name` wins whenever it is given. Otherwise the report uses the name the
*run itself* recorded in `pipeline_summary.json` — the same string Kestrel embedded
in its output filenames and VCF header, so one run cannot carry two identities. The
summary records where that name came from, in `sample_name_is_explicit`:

- `true` — the operator gave `vntyper pipeline --sample-name`. It is printed
  verbatim, whatever it looks like.
- `false` — the pipeline derived it from an input path, and it is a `Path.stem`
  rather than a finished name: `S1_R1.fastq.gz` records `S1_R1.fastq`. The
  derivation rule below finishes it.

A summary that records no name at all — every run archived before the field existed
— falls back to deriving one from the input file names it recorded. The rule is the
same either way: one recognised compound extension (`.fastq.gz`, `.fq.gz`, `.bam`,
`.cram`, …) is stripped, then a single **trailing** `_R1` or `_R2`. So
`example_b178_hg19_subset_R1.fastq.gz` becomes `example_b178_hg19_subset` and
`S1.lane3.L001_R1.fastq.gz` becomes `S1.lane3.L001`.

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
