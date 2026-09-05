# report

Generate an HTML summary report and visualizations from pipeline output data.

## Synopsis

```
vntyper [global-options] report
    -o <dir>
    [--input-dir <dir>]
    [--report-file <name>]
    [--bed-file <file>] [--bam-file <file>] [--vcf-file <file>]
    [--reference-fasta <file>]
    [--flanking <int>]
    [-s <name>]
    [--report-igv {embedded,sidecar,off}]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output-dir` | path | (required) | Target directory containing pipeline results |
| `--input-dir` | path | None | Source directory to search recursively for pipeline output artifacts |
| `--report-file` | string | `summary_report.html` | File name for generated HTML report (written inside `--output-dir`) |
| `--bed-file` | path | None | Path to BED file for IGV tracking |
| `--bam-file` | path | None | Path to Kestrel resolved haplotype-record BAM |
| `--vcf-file` | path | None | Path to VCF file for IGV variant track (auto-discovered if omitted) |
| `--reference-fasta` | path | (from config) | Reference FASTA file for IGV viewer. Defaults to configured `muc1_reference_vntr` |
| `--flanking` | int | `50` | Flanking base-pair margin for alignment viewing |
| `-s, --sample-name` | string | (derived) | Explicit sample name displayed in report titles and headings |
| `--report-igv` | choice | `embedded` | Alignment viewer packaging mode (`embedded`, `sidecar`, `off`) |

## The alignment browser

The HTML report operates offline without remote asset fetching. Vendored igv.js 3.0.2 is stored compressed and unpacks in the browser using `DecompressionStream`. The report makes zero network requests: scripts, stylesheets, fonts, and genome references load locally.

| Value | Implementation | Size |
|-------|----------------|------|
| `embedded` (default) | Self-contained report embedding compressed viewer code. No secondary files or network required | ~560 KB with alignments |
| `sidecar` | Core report references external `igv_report.html` generated in the same directory | ~75 KB, plus sidecar |
| `off` | Viewer omitted; tables and call summaries remain fully populated | ~75 KB |

Kestrel `output.bam` contains resolved haplotype records, not sequencing reads. Its record counts are haplotype-record support; `XD` is minimum k-mer depth and does not weight votes or alter names or tiers. The report prints this concise clarification in the reading key even in `off` mode, and the IGV panel identifies the track as resolved haplotype-record alignments.

The report Provenance block records the full SHA-256 digest of the decompressed viewer library, verified during Python rendering.

Embedded rendering requires browser support for `DecompressionStream` (Chrome 80+, Safari 16.4+, Firefox 113+). Older browsers display an informational note while displaying variant tables without modification.

## Reading and printing the report

The report header displays status chips for Kestrel, adVNTR, concordance, and Coverage QC alongside configured interpretive notes without generating verdict words. Screening provenance fields display pipeline execution states directly.

Individual sample reports preserve all variant rows without filtering. Flags provide visual emphasis without removing call evidence.

For printing, tables and motif sequences expand fully while interactive buttons hide. Chromium print styles format headers through `@page` margin rules. Disable browser print-preview collapse options to capture complete records.

## Naming the sample

Sample naming priority:

1. Explicit `--sample-name` option (printed verbatim).
2. Sample name recorded in `pipeline_summary.json`.
3. Derivation from input file paths: strips standard extensions (`.fastq.gz`, `.bam`, `.cram`) and trailing mate tags (`_R1`, `_R2`).
4. Fallback: `unnamed sample` when no inputs are recorded.

## Provenance

Reports capture reference assembly, header build identity, target coordinates, and summary schema version directly from `pipeline_summary.json`.

Render timestamps record generation time without altering original run timestamps.

Schema-3 executions verify the run-local decision profile snapshot (`decision_profile.json`) against recorded SHA-256 hashes. Missing or modified snapshots trigger fatal errors; older schema records report `decision profile not recorded by legacy run`.

Subcommand `report` does not accept `--decision-profile`: regenerating reports preserves the original run decisions.

## Auto-Discovery

When `--input-dir` is provided and file options are omitted, the command searches for standard artifacts:

- **BAM file:** `<input-dir>/kestrel/output.bam` (Kestrel resolved haplotype records)
- **BED file:** `<input-dir>/kestrel/output.bed`
- **VCF file:** `<input-dir>/kestrel/output.vcf`

## Examples

Generate a report from pipeline results:

```bash
vntyper report -o results/ --input-dir results/
```

Generate report specifying IGV alignment files:

```bash
vntyper report -o results/ \
    --bam-file results/kestrel/output.bam \
    --bed-file results/kestrel/output.bed \
    --reference-fasta ref/muc1.fa \
    --flanking 100
```

Specify a custom output file name:

```bash
vntyper report -o results/ --input-dir results/ --report-file my_report.html
```
