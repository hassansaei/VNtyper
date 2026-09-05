# pipeline

Run the full VNtyper 2 genotyping pipeline.

## Synopsis

```
vntyper [global-options] pipeline
    (--bam <file> | --cram <file> | --fastq1 <file> --fastq2 <file>)
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
| `-o, --output-dir` | path | `out` | Output directory for the results. For BAM/CRAM input it must stay outside the directory containing the alignment and all descendants of that directory |
| `-n, --output-name` | string | `output` | Base name for the output files. **Fixed at `output`**: any other value is rejected. The report generator, the `report` subcommand and the Kestrel stage each name their files from that literal and take no basename argument, so moving only the stages that accept one would leave the report reading files nothing wrote — which VNtyper reports as a negative genotype rather than an error. Use `--output-dir` to separate runs. |
| `-s, --sample-name` | string | (from input filename) | Sample name for labeling results. If not provided, defaults to the input BAM or FASTQ filename stem |

## Reference & Region Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--reference-assembly` | choice | `hg19` | Reference assembly for BAM/CRAM alignment. Options: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ncbi`, `hg38_ncbi`, `hg19_ensembl`, `hg38_ensembl`. For BAM/CRAM input this is **checked against the alignment header** and a disagreement stops the run — see [Reference Assemblies](../user-guide/reference-assemblies.md#the-declared-assembly-check) |
| `--custom-regions` | string | — | Custom regions for MUC1 analysis as comma-separated, 1-based inclusive values (e.g., `chr1:1000-2000,chr2:3000-4000`); VNtyper converts them to BED internally |
| `--bed-file` | path | — | Path to a BED file specifying regions with standard 0-based, half-open coordinates; VNtyper copies its rows verbatim |

`--custom-regions` and `--bed-file` are mutually exclusive.
For compatibility with established read-routing results, configured default extraction windows retain their
historical literal BED starts, including defaults supplied through `--config-path`. User-supplied inline regions
use the documented 1-based conversion. A custom boundary can intersect the VNTR, so changing it can change the
genotype result.

## Processing Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--threads` | int | `4` | Number of threads to use |
| `--fast-mode` | flag | off | Enable fast mode (skips filtering for unmapped and partially mapped reads) |
| `--keep-intermediates` | flag | off | Compatibility flag: intermediate files (BAM slices, temporary files) are already kept by default, so this flag changes nothing. Use `--delete-intermediates` to remove them. |
| `--delete-intermediates` | flag | off | Delete intermediate files after processing (wins when `--keep-intermediates` is also given). |
| `--decision-profile` | path | packaged profile | Select exactly one complete `explicit-custom` or `generated` decision profile. Partial overlays and a profile whose fixed-safety fields differ from the packaged profile are rejected before pipeline output is created. |

Omitting `--decision-profile` uses the verified packaged profile and preserves the
package's existing decisions. The option does not merge with `--config-path`: runtime
paths, tools, references, coverage presentation, and the adVNTR model/version guards
remain in their existing configuration or code and are not profile fields. A generated
profile is inert after it is written; it affects a run only when an operator selects its
complete file with `--decision-profile`.

## Optional Modules

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--extra-modules` | string | `[]` | Optional extra modules to include: `advntr`, `shark`. Can be repeated (`--extra-modules advntr --extra-modules shark`) or given as a comma-separated list (`--extra-modules advntr,shark`). Names are case-insensitive, and an unknown name is rejected rather than ignored |
| `--advntr-max-coverage` | int | — | Max coverage (e.g., 300) for quick adVNTR mode. Only applies when `advntr` is in `--extra-modules` |

The `shark` module is not supported in BAM/CRAM mode; use FASTQ mode or remove the shark flag.

## Archive & Summary Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--archive-results` | flag | off | Create an archive of the results folder after pipeline completion |
| `--archive-format` | choice | `zip` | Format of the archive: `zip` or `tar.gz` |
| `--summary-formats` | string | `""` | Comma-separated list of additional summary output formats (supported: `csv`, `tsv`). Each format writes `pipeline_summary.<fmt>` (one row per step, run provenance first) and `pipeline_summary_rows.<fmt>` (one row per result field). JSON is always generated; unknown names are ignored without a message |
| `--report-igv` | choice | `embedded` | How the report carries its alignment browser: `embedded`, `sidecar` or `off`. See [Report Options](#report-options) |

## Report Options

The HTML report is written by the same code whether it comes from `vntyper pipeline` or
from `vntyper report`, so `--report-igv` means the same thing on both and defaults to
the same value.

| Value | What the report is | Size |
|-------|--------------------|------|
| `embedded` (default) | One file. igv.js travels inside it, gzipped, and expands when the report is opened. No second file, no network. | ~560 KB with alignments |
| `sidecar` | The summary report carries no library and points at `igv_report.html`, written beside it. That file is self-contained too. | ~75 KB, plus the sidecar |
| `off` | No alignment browser is produced at all, and `create_report` is not run. | ~75 KB |

A report for a run with no alignment data is ~75 KB whatever the mode: the library is
written only when there is a session for it to show.

**Neither report fetches anything.** Every mode produces a document that resolves
nothing off-machine — no CDN, no fonts, no genome registry. That is measured in a real
browser rather than asserted, because the request that used to escape came from inside
the minified igv.js bundle rather than from any tag in the document.

## Examples

The BAM examples use separate `inputs/` and `results/` trees. An output directory beside
the BAM itself would still be inside the protected patient input tree and is rejected.

Run the pipeline with a BAM file using default settings:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/
```

Run with paired FASTQ files and hg38 reference:

```bash
vntyper pipeline --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
    --reference-assembly hg38 -o results/ -s my_sample
```

Enable fast mode with multiple threads and archive results:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ --threads 8 \
    --fast-mode --archive-results --archive-format tar.gz
```

Run with the adVNTR module and coverage cap:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --extra-modules advntr --advntr-max-coverage 300
```

Run with one explicitly selected complete decision profile:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --decision-profile reviewed-complete-profile.json
```

Generate additional summary formats and clean up intermediate files:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --summary-formats csv,tsv --delete-intermediates
```

This writes `pipeline_summary.csv`, `pipeline_summary_rows.csv`, `pipeline_summary.tsv` and
`pipeline_summary_rows.tsv` beside `pipeline_summary.json`; see
[Output Files](../user-guide/output-files.md) for their columns.
