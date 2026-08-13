# Report Generation

VNtyper 2 produces two types of reports: a per-sample HTML summary report with embedded IGV visualization, and an optional multi-sample cohort summary report with aggregated statistics.

## Sample Report

The sample report (`summary_report.html`) is generated at the end of each pipeline run using a Jinja2 HTML template. It integrates results from all pipeline stages into a single document.

### Report Contents

**Variant Summary Table**

The Kestrel results are displayed in a sortable table with columns for motif, variant type, position, REF/ALT alleles, motif sequence, depth metrics, depth score, confidence level, and flag status. Confidence labels are color-coded:

- High_Precision / High_Precision* -- highlighted in red (positive finding)
- Low_Precision -- highlighted in orange (requires validation)
- Negative -- no color (no variant detected)

If adVNTR was run, its results appear in a separate table showing VID, variant state, supporting read count, mean coverage, p-value, repeat unit, REF/ALT, and flag status.

**Screening Summary**

The report opens with a masthead: who the report is about, the computed state as a row of labelled chips, and then the interpretive text. The text is generated from a rule-based system defined in `report_config.json` that considers:

- Kestrel result category (High_Precision, Low_Precision, flagged variants, negative)
- adVNTR result category (positive, negative, not performed)
- Quality metrics pass/fail status

The screening summary states what the combined evidence supports, including where orthogonal validation is recommended. Each configured message is stored both verbatim (`message`) and as the ordered parts it is rendered in (`segments`); the two are kept in step by a round-trip check, and a configuration supplying only `message` still renders.

A stage that ran without producing a readable result is reported as neither a positive nor a negative. When either genotyping stage is in that state, the configured message is withheld -- it was selected by matching a state the run never established -- and the masthead is drawn in its indeterminate state instead.

**Cross-Match Summary**

When both Kestrel and adVNTR results are available, the report indicates whether at least one concordant variant was found between the two methods.

### QC Metrics

The report includes quality metrics from multiple sources:

| Metric | Source | Threshold |
|--------|--------|-----------|
| Mean VNTR coverage | samtools depth | >= 100x |
| Percent VNTR uncovered | samtools depth | <= 50% |
| Duplication rate | fastp | <= 10% |
| Q20 rate | fastp | >= 80% |
| Q30 rate | fastp | >= 70% |
| Passed filter rate | fastp | >= 80% |

Each metric is displayed with a color-coded indicator (green check or red warning) based on its threshold.

The two coverage rows also decide a verdict. The report shows it as a **Coverage QC** row reading `PASS` or `FAIL`, and the same value is written to `coverage_summary.tsv` as a `coverage_qc` column and reaches the cohort exports as `cov_coverage_qc`. `FAIL` on either metric is what makes the screening summary report the sample's quality metrics as below threshold.

Before VNtyper 2.0.8 the uncovered row was displayed and never enforced: only the mean decided the quality gate, so a sample with acceptable mean coverage and half the VNTR uncovered passed. The table above now describes what the code does.

**BAM Header Information**

For BAM/CRAM input, the report displays the detected reference assembly (from both text and contig matching), alignment pipeline, and any associated warnings.

### IGV Integration

The report embeds an interactive IGV genome browser view using the [igv-reports](https://github.com/igvteam/igv-reports) library. The IGV view is generated from:

- **BAM track** -- Kestrel's haplotype alignment output (`output.bam`)
- **VCF track** -- filtered INDEL variants (`output_indel.vcf.gz` or `.vcf`)
- **BED track** -- variant position file (`output.bed`)
- **FASTA reference** -- MUC1 VNTR reference sequence

The flanking region parameter (default: 50 bp) controls how much sequence context is shown around each variant position. This is configurable in `config.json` under `default_values.flanking`.

!!! info "VCF compression"
    If bcftools is installed, the INDEL VCF is compressed and sorted for optimal IGV performance. If bcftools is unavailable, the uncompressed VCF is used. The report generation handles both cases gracefully.

### Printing and archiving

The report is written to be printed: the printed sheet is the artefact that gets filed, forwarded and read years later, so it carries its own identity. **Every** page repeats the sample, the assay, the reference assembly, the VNtyper version and the time the run started in the page margin, and is numbered `Page N of M` so a print that lost a sheet is visible as one — a sheet separated from the rest still says what it is a sheet of. The detailed coverage table prints in place of the basic view, results tables print in full rather than clipped to the screen's column widths, and the on-screen switches are omitted because paper has no switches.

The running header is a paged-CSS feature that only Chromium implements, so the first page additionally states the same identity as an ordinary line of the document. In Chromium that line is redundant with the margin; in Firefox and WebKit, which drop page margin boxes silently, it is the only identity the printed record has.

The pipeline log is the deliberate exception. It prints as a one-line pointer back to the HTML original rather than as pages of DEBUG output, whether or not the reader had it expanded.

!!! warning "One thing a collapsed section costs a reader with JavaScript disabled"
    Sections such as **Variants** are collapsible. The report is served with them **open**, so an ordinary print carries them; if a reader collapses one and then prints, a small script reopens it for the duration of the print and closes it again afterwards.

    That script is the only mechanism that works. Measured in Chromium 151: a section the reader collapsed prints its heading and **none of its contents** when JavaScript is disabled, and no stylesheet can change that — `open` is an HTML attribute rather than a style, and the browser hides a closed section's contents in a way an author stylesheet cannot reach. Firefox and WebKit ignore the paged-CSS features this report relies on altogether.

    So with scripting off, print the report **as it opens** — do not collapse a section first. If a print is missing a section, its heading will still be on the page where the contents should have been.

## Cohort Report

The cohort summary module (`cohort_summary.py`) aggregates results from multiple pipeline runs into a single report. It scans a directory structure for `pipeline_summary.json` files and constructs:

- **Sample result table** -- aggregated variant calls, confidence levels, and flags across all samples
- **Donut charts** -- interactive Plotly visualizations showing the distribution of results (positive/negative/low precision) across the cohort
- **Coverage statistics** -- per-sample VNTR coverage metrics
- **Runtime statistics** -- pipeline execution times
- **Version and assembly tracking** -- VNtyper 2 versions and detected reference assemblies

### Pseudonymization

The cohort report can replace sample identifiers with a digest, which keeps names out of the tables an aggregated result is read from. It is obfuscation for readability, **not** de-identification -- the digest is unsalted and unkeyed, so a shared report must still be treated as identifying; see [what it does not protect against](../user-guide/cohort-analysis.md#what-pseudonymization-does-not-protect-against). Each sample name is replaced with a prefix (default `sample_`) followed by the first 12 hex characters of its SHA-256 digest. Both the algorithm and the width are configurable under `cohort.pseudonym` in `config.json`.

The mapping is injective by construction: if two samples would share a reported name, the run raises an error naming both rather than merging their genotypes into one row. See [Cohort Analysis](../user-guide/cohort-analysis.md#pseudonymization) for the configuration block and what to do about a collision.

### Report Configuration

Report behavior is controlled by `report_config.json`, which defines:

- **Algorithm logic rules** -- how Kestrel and adVNTR results map to categorical outcomes (e.g., "High_Precision" + "Not flagged" = positive)
- **Screening summary rules** -- condition-to-message mappings for the interpretive summary text
- **Default messages** -- fallback text when no rule matches
