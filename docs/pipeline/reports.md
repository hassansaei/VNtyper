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

An interpretive text block summarizes the clinical significance of the combined results. The summary is generated from a rule-based system defined in `report_config.json` that considers:

- Kestrel result category (High_Precision, Low_Precision, flagged variants, negative)
- adVNTR result category (positive, negative, not performed)
- Quality metrics pass/fail status

The screening summary provides actionable guidance, including recommendations for orthogonal validation when appropriate.

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

## Cohort Report

The cohort summary module (`cohort_summary.py`) aggregates results from multiple pipeline runs into a single report. It scans a directory structure for `pipeline_summary.json` files and constructs:

- **Sample result table** -- aggregated variant calls, confidence levels, and flags across all samples
- **Donut charts** -- interactive Plotly visualizations showing the distribution of results (positive/negative/low precision) across the cohort
- **Coverage statistics** -- per-sample VNTR coverage metrics
- **Runtime statistics** -- pipeline execution times
- **Version and assembly tracking** -- VNtyper 2 versions and detected reference assemblies

### Pseudonymization

The cohort report supports sample pseudonymization by hashing sample identifiers, allowing sharing of aggregated results without exposing patient identifiers. Each sample name is replaced with a prefix (default `sample_`) followed by the first 5 characters of its MD5 hash.

### Report Configuration

Report behavior is controlled by `report_config.json`, which defines:

- **Algorithm logic rules** -- how Kestrel and adVNTR results map to categorical outcomes (e.g., "High_Precision" + "Not flagged" = positive)
- **Screening summary rules** -- condition-to-message mappings for the interpretive summary text
- **Default messages** -- fallback text when no rule matches
