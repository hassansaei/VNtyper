# Report Generation

VNtyper produces two report types: a per-sample HTML summary report with embedded IGV visualization, and an optional multi-sample cohort summary report with aggregated statistics.

## Sample Report

The sample report (`summary_report.html`) is generated at the end of each run via a Jinja2 template, compiling results from all pipeline stages into a single document.

### Report Contents

**Variant Summary Table**

Displays Kestrel results in a sortable table with columns for motif, variant type, position, REF/ALT alleles, depth metrics, depth score, confidence level, flag status, full nomenclature records, and motif sequences. The nomenclature record includes the reconciled MUC1 name, confidence tier, explanatory flags, caller-specific names, ambiguity interval, repeat form, and notes. See [MUC1 Nomenclature](nomenclature.md).

The reading key under the result tables includes this concise artifact clarification:
Kestrel `output.bam` contains resolved haplotype records, not sequencing reads. Its
record counts are haplotype-record support; `XD` is minimum k-mer depth and does not
weight votes or alter names or tiers. This clarification remains visible when the
alignment browser is off or no BAM is available.

Confidence is displayed as a labelled indicator:

- High_Precision / High_Precision*: High-confidence call.
- Low_Precision: Marginal-confidence call.
- Negative: No variant detected.

If adVNTR was run, its results appear in a separate table showing VID, variant state, supporting read count, mean coverage, p-value, repeat unit, REF/ALT, flag status and the same complete mutation-naming record.

**Screening Summary**

The masthead displays sample metadata, computed result chips, and interpretive text. Chips indicate Kestrel status, adVNTR status, concordance, confidence grade, VNTR depth, mean depth, and flank depth. Unmatched configurations display `Screening rule: Not configured`. Chip words come from the `chip_labels` block of `report_config.json`; a token it does not name is shown with its underscores and hyphens as spaces (`High_Precision` reads **High precision**). VNTR depth reads **Adequate**, **Reduced**, **Insufficient** or **Not evaluated**.

The **Confidence grade** chip conveys sample-level confidence from `confidence_grade_rules` in `report_config.json`:

- **Finding tone** (amber highlight): `Finding` (`finding`), `Finding, corroborated` (`finding-corroborated`: Kestrel and adVNTR agree).
- **Caution tone** (warning highlight): `Finding, with caveats` (`finding-limited`), `No finding, inconclusive` (`no-finding-limited`), `No result` (`not-established`).
- **Neutral tone**: `No finding` (`no-finding`).

When custom configurations omit `confidence_grade_rules`, the chip is omitted without error.

The masthead prints no arbitrary verdict text. Styling reflects internal states (`finding`, `no-finding`, `indeterminate`), while text is populated from configuration rules based on:

- Kestrel result category (High_Precision, Low_Precision, flagged variants, negative)
- adVNTR result category (positive, negative, not performed)

The message never states coverage. What the measured coverage means for the result is stated beneath it by the `coverage_notes` sentences of `report_config.json`, filled in with the measured figure and the threshold applied. Mean depth qualifies a result without a finding and never a Kestrel call (an adVNTR-only call is graded `Finding, with caveats` at failing depth): a reduced depth adds one sentence, a failing depth adds a sentence and a follow-up and grades the result `No finding, inconclusive`. A region that is mostly uncovered is stated for every result and limits the grade of a call too, as does a run with no coverage measurement. A Kestrel candidate below the reporting floor is described in the message and in the Kestrel section, and does not change the grade. Every call ends with one sentence asking for confirmation by an independent method.

Each message is tracked both verbatim (`message`) and as segmented units (`segments`).

Stages that fail to produce readable results are treated as indeterminate rather than negative, withholding interpretive text and drawing the masthead in neutral style.

A screening provenance line reports raw caller findings and Coverage QC. The **Provenance** block records requested and detected assemblies, analysed target coordinates, and schema versions. Missing legacy fields display `not recorded by this run`.

**Cross-Match Summary**

When both callers run, reports indicate whether concordant variants were detected.

### QC Metrics

Quality metrics are evaluated across tools:

| Metric | Source | Configured threshold | Evaluation rule |
|--------|--------|----------------------|-----------------|
| Mean VNTR coverage | samtools depth | `thresholds.mean_vntr_coverage_by_assembly`, `thresholds.mean_vntr_coverage` | per assembly: adequate at >= `adequate`, reduced at >= `low`, failing below it (GRCh38 100/50, GRCh37 290/145); an unknown assembly fails below `mean_vntr_coverage` |
| Percent VNTR uncovered | samtools depth | `thresholds.percent_vntr_uncovered` | measured <= configured value |
| Duplication rate | fastp | `thresholds.duplication_rate` | measured <= configured value |
| Q20 rate | fastp | `thresholds.q20_rate` | measured >= configured value |
| Q30 rate | fastp | `thresholds.q30_rate` | measured >= configured value |
| Passed filter rate | fastp | `thresholds.passed_filter_reads_rate` | measured >= configured value |

Each metric is displayed with a status indicator based on its configured threshold. The four configured fastp cutoffs are required finite numeric fractions from 0 through 1; a missing, malformed, non-finite, or out-of-range cutoff is logged and raises `ValueError` during report rendering. Each nonmissing measured fastp rate from `output.json` must also be a finite numeric fraction from 0 through 1; malformed, non-finite, and out-of-range rates are logged and raise `ValueError` with their metric key. A `None` rate remains missing and renders as `N/A` with no status icon. JSON fractions enter this decision boundary as exact decimal values, without first being rounded through a binary float. Each measured rate and configured cutoff is rounded half-up on its exact decimal value to two decimal places of percent before display and comparison, and equality passes.

The passed-filter rate divides `filtering_result.passed_filter_reads` by `summary.before_filtering.total_reads`. Total reads must be non-negative integers. An absent `output.json` marks fastp metrics as unavailable; corrupt files raise `ValueError`.

Coverage metrics determine the **Coverage QC** status (`PASS`, `REDUCED` or `FAIL`), written to `coverage_summary.tsv` as `coverage_qc` and exported to cohorts as `cov_coverage_qc`. `REDUCED` is a mean between the assembly's `low` and `adequate` lines; it passes. `FAIL` is a mean below `low` or an uncovered percentage above its own. The mean's status indicator marks only a failing mean. See [Coverage thresholds](../user-guide/configuration.md#coverage-thresholds).

**BAM Header Information**

Displays detected reference assembly, alignment software, and associated warnings.

### IGV Integration

Embeds interactive IGV genome browser views via [igv-reports](https://github.com/igvteam/igv-reports):

- **BAM track**: Kestrel's resolved haplotype records in `output.bam`, not sequencing reads.
- **VCF track**: Filtered INDEL variants (`output_indel.vcf.gz` or `.vcf`).
- **BED track**: Variant locus file (`output.bed`).
- **FASTA reference**: MUC1 VNTR reference sequence.

The flanking display window defaults to 50 bp (`default_values.flanking`).

Configured via `--report-igv`:

| Mode | Behavior |
|------|----------|
| `embedded` (default) | Compresses igv.js library and alignment session directly within `summary_report.html`. |
| `sidecar` | Writes an external self-contained `igv_report.html` file alongside the report. |
| `off` | Omits the genome browser while retaining tabular results. |

Self-contained reports without alignments measure ~78 KB; embedded alignments measure ~576 KB. Decompression uses native browser `DecompressionStream` (Chrome 80+, Safari 16.4+, Firefox 113+).

### Printing and Archiving

The print stylesheet includes running margin headers with sample ID, assay name, reference build, VNtyper version, and timestamp. Pages are numbered `Page N of M` so a print that lost a sheet is visible as one: a sheet separated from the rest still identifies what it belongs to. Result tables print in full width without scrollbars. Collapsible sections print open. Pipeline execution logs print as a single-line reference back to the digital file.

### Custom template context compatibility

Custom Jinja template directories configured via `paths.template_dir` interact with a public compatibility interface. VNtyper 2.x provides the deprecated keys below even though the shipped template no longer reads them. They remain available throughout VNtyper 2.x and may be removed in **VNtyper 3.0.0**, not before:

- Color keys: `percent_vntr_uncovered_color`, `mean_vntr_coverage_color`, `duplication_rate_color`, `q20_color`, `q30_color`, `passed_filter_color`.
- Raw fastp fractions and display replacements:

  | Deprecated raw key (available through 2.x) | Display replacement |
  | --- | --- |
  | `duplication_rate` | `duplication_rate_display` |
  | `q20_rate` | `q20_rate_display` |
  | `q30_rate` | `q30_rate_display` |
  | `passed_filter_rate` | `passed_filter_rate_display` |

- Alignment and state keys: `igv_content`, `screening_state.kestrel_result`, and `screening_state.advntr_result`.

## Cohort Report

The cohort module (`cohort_summary.py`) aggregates runs across directories:

- **Sample result table**: Aggregated calls, confidence tiers, and nomenclature.
- **Donut charts**: Plotly visualizations of call distributions.
- **Coverage statistics**: Multi-sample VNTR coverage metrics.
- **Runtime metrics**: Step execution times.
- **Call frequency table**: Variant frequencies across the cohort.

When a cohort contains a BAM-specific nomenclature flag, its reading key gives the same resolved haplotype records and `XD` minimum k-mer depth explanation as the sample report, noting that it does not weight votes or alter names or tiers.

### Cohort Call Frequency

Aggregates variants across the cohort:

- **Grouping key**: Grouped by canonical `Molecular_Identity` when resolved; unresolved variants fall back to caller representation `(Motifs, POS, REF, ALT)`.
- **Grouping key kind**: Tracked as `molecular-identity` or `caller-representation` in `Grouping_Key_Kind`.
- **Denominator**: Total sample roster, including negative and indeterminate runs.
- **Rare variant flag**: Calls with `Frequency <= rare_allele_max_frequency` set `Below_Cutoff = "yes"`.

### Pseudonymization

Anonymizes sample identifiers for presentation using the prefix and first 12 hex characters of their unsalted SHA-256 digest (for example, `sample_a1b2c3d4e5f6`). Configurable via `cohort.pseudonym`. The mapping is strictly injective: identifier collisions raise an error.

### Report Configuration

Behavior is controlled via `report_config.json`, defining algorithm logic mappings, screening summary text rules, and fallback messages.
