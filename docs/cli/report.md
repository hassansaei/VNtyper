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
| `embedded` (default) | One file, and a complete alignment browser. | 575,641 B with alignment data in the verification specimen |
| `sidecar` | No library in the summary; `igv_report.html` beside it is the alignment browser, built from VNtyper's verified vendored igv.js with no network request. | Summary remains near the no-alignment size, plus the sidecar |
| `off` | No alignment browser at all; the complete result tables remain. | 78,365 B in the no-alignment specimen |

A report for a run with no alignment data carries no library in any mode. The last
verification specimen measured **78,365 bytes without alignment data** and **575,641
bytes with embedded alignment data**, against **2,002,405 bytes** fetched over 11 CDN
tags by the old report. These numbers are measurements to be refreshed for a release,
not size guarantees; sample content also affects the artifact.

The report's Provenance block states the version and the **SHA-256 of the decompressed
library**, in full, so it can be checked against upstream. That digest is verified in
Python at render time, over the exact bytes the document later evaluates. It is
deliberately not re-checked in the browser: `crypto.subtle` is unavailable on `file://`
in several engines, and an archived report is opened from `file://` more often than
not, so a runtime check would be one that silently does not run on the machines it was
written for.

The embedded mode has a 2023 cross-browser floor because it requires
`DecompressionStream`: Chrome 80+, Safari 16.4+ or Firefox 113+. On an older browser,
the alignment panel says so in words and the variant tables stand unchanged. No state
of this report is blank space.

## Reading and printing the report

The masthead names the sample and presents labelled Kestrel, adVNTR, concordance and
Coverage QC chips before the configured interpretive text. It deliberately produces
**no verdict word**: the internal emphasis state changes styling only. A screening-
provenance line then spells out the Kestrel state, adVNTR state and Coverage QC using
pipeline vocabulary. Coordinate fields absent from an older summary read **not recorded
by this run** rather than being guessed from today's configuration.

The per-sample report never hides a variant row. Flag highlighting can change emphasis,
but it cannot filter the evidence. Flag filtering remains available in the cohort report
as a multi-sample triage affordance.

For print, result tables and the full motif sequence are expanded and interactive controls
are omitted. Chromium repeats identity through `@page` margin boxes; a first-page document
block supplies the identity in browsers that do not implement those boxes, so Chromium
shows both on page 1. The report opens print-relevant sections expanded and scripting
temporarily reopens any the reader collapsed. With scripting disabled, a section that the
reader collapsed prints without its contents; print the report as it opens. See
[Printing and archiving](../pipeline/reports.md#printing-and-archiving) for the full
limitations.

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
