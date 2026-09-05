# Output Files

## Directory Structure

```
results/
├── pipeline_summary.json        # Machine-readable pipeline summary
├── pipeline_summary.csv         # Optional (--summary-formats csv): one row per step
├── pipeline_summary_rows.csv    # Optional (--summary-formats csv): one row per result field
├── pipeline_summary.tsv         # Optional (--summary-formats tsv): one row per step
├── pipeline_summary_rows.tsv    # Optional (--summary-formats tsv): one row per result field
├── pipeline.log                 # Pipeline execution log
├── summary_report.html          # Self-contained HTML report (IGV mode is configurable)
├── igv_report.html              # Optional self-contained sidecar (--report-igv sidecar)
├── predefined_regions_<assembly>.bed  # Region BED file (e.g., hg19, hg38)
├── kestrel/
│   ├── kestrel_result.tsv       # Final genotyping result
│   ├── kestrel_pre_result.tsv   # Pre-filter variants (all candidates)
│   ├── output.vcf               # Raw Kestrel VCF
│   ├── output_indel.vcf         # Filtered INDEL VCF
│   ├── output_indel.vcf.gz      # Compressed INDEL VCF (if bcftools available)
│   ├── output.bam               # Kestrel resolved haplotype-record BAM
│   ├── output.bam.bai           # BAM index
│   └── output.bed               # BED file for coverage visualization
├── fastq_bam_processing/
│   ├── output_R1.fastq.gz       # Extracted/processed R1 reads
│   ├── output_R2.fastq.gz       # Extracted/processed R2 reads
│   ├── output_other.fastq.gz    # Extracted unpaired reads (if present)
│   ├── output_single.fastq.gz   # Extracted single-end reads (if present)
│   ├── output_sliced.bam        # Region-sliced BAM extracted from input
│   ├── output.json              # fastp QC summary metrics
│   ├── output.html              # fastp QC HTML report
│   └── pipeline_info.json       # BAM header metadata (BAM/CRAM input)
├── alignment_processing/
│   └── output_sorted.bam        # BWA-aligned BAM (FASTQ input only)
├── coverage/
│   └── coverage_summary.tsv     # VNTR region coverage statistics
├── advntr/                      # Only when --extra-modules advntr
│   ├── output_adVNTR.tsv        # Raw adVNTR output
│   ├── output_adVNTR.vcf        # adVNTR variant calls in VCF format
│   ├── output_adVNTR_result.tsv # Processed adVNTR result
│   ├── cross_match_results.tsv  # Kestrel vs adVNTR comparison
│   └── advntr_model.db          # Resolved model database
└── provenance/                  # Checkpoint and decision profile artifacts
    ├── advntr_artifact_evidence.json
    └── decision_profile.json
```

## kestrel_result.tsv Columns

This is the primary output file. Each row represents a genotyped variant.

| Column | Description |
|--------|-------------|
| `Motif` | MUC1 repeat motif identifier (e.g., `1`, `2`, `3`) |
| `Variant` | Variant type (`Insertion` or `Deletion`) |
| `POS` | Position within the MUC1 reference sequence |
| `REF` | Reference allele |
| `ALT` | Alternate allele |
| `Motif_sequence` | 60 bp sequence of the repeat-unit half named in `Motif` (one half of the 120 bp reference pair record) |
| `Estimated_Depth_AlternateVariant` | k-mer-path depth supporting the alternate allele |
| `Estimated_Depth_Variant_ActiveRegion` | total k-mer depth in the variant active region |
| `Depth_Score` | Ratio of alternate depth to active region depth |
| `Confidence` | Confidence classification (see below) |
| `Flag` | Quality flag (`Not flagged` or a flag reason) |
| `haplo_count` | Number of haplotype calls supporting the same variant |
| `Nomenclature` | Reconciled variant name, e.g. `59dupC` (see [MUC1 Nomenclature](../pipeline/nomenclature.md)) |
| `Nomenclature_Tier` | `A`, `B`, or `C`: how much of the name may be stated |
| `Nomenclature_Flags` | `;`-separated reasons the tier is what it is |
| `Ambiguity_Interval` | Span in which every anchor is the same allele, e.g. `53_59`. Empty when the variant cannot shift |
| `Repeat_Form` | Tract copy-number change, e.g. `53C[7]>53C[8]`. Empty outside a detectable tract |
| `Nomenclature_Kestrel` | What Kestrel named on its own |
| `Nomenclature_adVNTR` | What adVNTR named on its own. Empty when the optional adVNTR module did not run |
| `Molecular_Identity` | Stable molecular identity serialization. Empty when unresolved |
| `Molecular_Identity_Status` | `unique`, `legacy-selected-among-multiple`, or `unresolved` |
| `Equivalent_Representation_Count` | Number of caller representations equivalent to the selected identity; integer `0` when unresolved |
| `Identity_Hypothesis_Count` | Number of distinct resolved identities considered for this caller result |

`Nomenclature` is the reconciled verdict; the two per-caller columns beside it say
what each caller reported, ensuring disagreements remain explicit in either result file. Both
files carry identical values for all three fields, enabling merged cohort tables.

`Ambiguity_Interval`, `Repeat_Form`, and `Nomenclature_adVNTR` are nullable: they are
written only where applicable, and left empty rather than padded with placeholders.

A negative run writes a narrower schema (first column `Motif`, omitting depth and
flag columns) and carries none of these fields: without a detected variant, there is no name or
molecular identity row.

Positive adVNTR result rows append the same four molecular-identity columns in the same
order. Negative adVNTR output retains its narrower schema. The Kestrel result
additionally carries `__Identity_*` capture and selection cells, and positive rows in both
caller result files carry `__Reconciled_Molecular_Identity`: the canonical serialization of
the whole-locus identity selected by the reconciler, or an empty cell after abstention. A run
under an explicit profile with dominance enabled also writes
`__Dominance_Abstention_Reason`, recording the closed abstention token or an empty cell. The
double-underscore prefix marks internal identity-capture and calibration-replay
persistence rather than public fields: they are copied into `pipeline_summary.json` alongside
other result columns, excluded from policy projections, omitted from HTML
reports, and never written to negative rows. The sample HTML
report, `pipeline_summary.json`, and cohort HTML/TSV/CSV/JSON exports carry all four
recorded identity values directly without inferring them from `POS`, `REF`, `ALT`, `Variant`, or
`Nomenclature`. In a sample HTML report, use table column controls to display the
quartet. The printed report folds these wide fields into a labelled
per-row appendix.

For compatibility, a schema-1 or schema-2 summary row missing any of the four fields
renders `legacy identity not recorded` across all four downstream identity cells. A complete
quartet is copied exactly, including empty unresolved identities, integer `0`
representation counts, and nonzero hypothesis counts.

## Kestrel output.bam evidence

Kestrel `output.bam` contains resolved haplotype records, not sequencing reads. A count
of those records is haplotype-record support. Its optional `XD` tag is the minimum k-mer
depth of one resolved haplotype and does not weight votes or alter names or tiers. The
HTML report states the same distinction in its reading key and labels its IGV track as
resolved haplotype-record alignments.

These units are separate from both Kestrel VCF k-mer-path depth and adVNTR supporting
read counts. The stable `Nomenclature_Flags` field can therefore contain
`thin-haplotype-record-support`, `low-haplotype-record-support`,
`low-kmer-path-support`, `low-read-support`, or `low-evidence-support`; see the
[authoritative flag table](../pipeline/nomenclature.md#flags) for source-specific
definitions.

## Confidence Levels

Confidence is assigned based on depth score thresholds calibrated from Saei et al., iScience 26, 107171 (2023).

| Level | Meaning |
|-------|---------|
| **High_Precision\*** | Depth score >= 0.00515 and alternate depth >= 100 |
| **High_Precision** | Depth score >= 0.00515, alternate depth 21 to <100, and region depth > 200 |
| **Low_Precision** | Variant detected with marginal depth or depth score support |
| **Negative** | No variant passed filtering thresholds |

!!! tip
    A result with confidence `Negative` indicates no MUC1-VNTR frameshift variant passed filtering thresholds: it does not guarantee the sample is biologically negative.

## How summary_report.html Displays Numbers

Every numeric column in the per-sample HTML report is formatted server-side and
written into the file by VNtyper.
Until issue #242, the report executed a client-side script that reformatted numeric
cells with `toFixed(4)` in the reader's browser, making on-screen numbers
dependent on external CDN reachability: the archived
file presented different precision online than offline.

Removing that script standardized column formatting. **No underlying values changed; only
rendering was adjusted.** `kestrel_result.tsv`, `pipeline_summary.json`, and adVNTR TSV
outputs remain untouched as the authoritative sources to parse. For scrapers parsing the HTML report,
the formatted columns are:

| Table | Column | Displayed before | Displayed now | Formatting Rule |
|-------|--------|------------------|---------------|-----------------|
| Kestrel | `POS`, `Estimated_Depth_AlternateVariant`, `Estimated_Depth_Variant_ActiveRegion` | `67`, `120`, `12000` | unchanged | Whole numbers, no decimals |
| Kestrel | `Depth_Score` | `0.01` | `0.010012` | Six decimal places. Confidence thresholds are calibrated to five decimals (0.00469 and 0.00515), so four decimals obscured threshold boundaries |
| adVNTR | `VID`, `NumberOfSupportingReads`, `POS` | `25561`, `14`, `67` | unchanged | Whole numbers, no decimals |
| adVNTR | `MeanCoverage` | `98.5`, `40` | `98.50`, `40.00` | Exactly two decimal places across all mean values |
| adVNTR | `Pvalue` | `0`, `0.0001` | `1e-09`, `0.000123` | Three significant figures; prevents tiny probabilities displaying as 0 |

### Column Order Migrations

`Motif_sequence` is now the **last** column of the per-sample report's Kestrel table,
following confidence, flag, and nomenclature fields. It previously sat sixth, between `ALT`
and `Estimated_Depth_AlternateVariant`.

Before correction to the selected 60 bp half, it contained the full 120 bp pair record.
In the sixth position, that long sequence forced `Confidence` and `Flag` off
standard 1280px displays. Moving it to the final column ensures primary diagnostic fields
remain visible without horizontal scrolling.

**This migration affects display only.** In result TSVs and reports, `Motif_sequence`
reflects the selected 60 bp half; no other fields changed. When scraping the HTML report by column
position, key on header text instead. The `cohort_summary.html` table defines its own column
order and is not affected.

Additional interface refinements:

- **adVNTR table headings use plain text**, matching the Kestrel layout:
  `NumberOfSupportingReads` displays as `Supporting Reads`, `Pvalue` as `P-value`, `RU` as
  `Repeat Unit`, and naming fields mirror Kestrel headings (`Nomenclature` is `MUC1 Name`,
  `Nomenclature_Tier` is `Tier`, `Ambiguity_Interval` is `Ambiguity`). `advntr_result.tsv`
  retains raw database names.
- **Semicolon-separated flags are spaced.** `Nomenclature_Flags` renders as
  `known-variant; motif-context-diverges; position-ambiguous`, permitting clean line breaks.
- **Numbers are tabular-aligned**, and column headers include hover definitions. Coded values
  (tier letters, nomenclature flags) are defined in an explicit **reading key printed beneath the tables**.
- **The `Flag` column displays full text reasons** beside visual icons rather than relying on tooltip hovers.
- **Flagged rows remain visible.** The report switch *Highlight flagged values* adjusts visual emphasis
  without removing rows from view.

## summary_report.html Artifact Modes

`--report-igv embedded` (default) produces a self-contained report. `sidecar`
writes a compact summary alongside a standalone `igv_report.html`, and `off` omits the
browser track while preserving all result tables. Verification benchmarks measure
**78,486 bytes without alignment data** and **575,762 bytes with embedded alignment
data**, compared to **2,002,405 bytes** retrieved across external CDN tags.

Embedded alignment viewing utilizes browser `DecompressionStream` support (supported
across Chrome 80+, Safari 16.4+, Firefox 113+). On older browsers, the alignment panel
displays a compatibility notice while the tables remain accessible.

## Pipeline Summary JSON

`pipeline_summary.json` records execution timestamps, input parameters, output paths, and parsed step outputs. It serves as the input format for [cohort analysis](cohort-analysis.md).

Key fields:

```json
{
  "schema_version": 2,
  "decision_policy": "legacy-selection-v1",
  "version": "{{ version }}",
  "pipeline_start": "2024-01-15T10:30:00",
  "pipeline_end": "2024-01-15T10:35:00",
  "input_files": { "bam": "sample.bam" },
  "steps": [
    {
      "step": "Kestrel Genotyping",
      "output_file": "kestrel/kestrel_result.tsv",
      "parsed_result": { "data": [...] }
    }
  ]
}
```

Molecular-identity records are additive to summary schema 2 without altering schema versioning. Current summaries record the packaged selection policy as
`decision_policy: legacy-selection-v1`. Older summaries lacking this key remain
supported and are marked as unrecorded.

### Flattened summary tables

`--summary-formats csv,tsv` writes two flattened files per format alongside `pipeline_summary.json`.
These files facilitate external analysis in tabular tools without JSON parsing.

`pipeline_summary.csv` / `pipeline_summary.tsv` contains one row per recorded step:

- **Run provenance (`run_*`)**: Top-level JSON fields repeated across
  every row (`run_schema_version`, `run_decision_policy`, `run_advntr_evidence_digest`,
  six `run_decision_profile_*` fields, `run_pipeline_start`, `run_pipeline_end`,
  `run_version`, `run_input_files_<kind>`, `run_sample_name`, `run_sample_name_is_explicit`,
  four `run_reference_*` fields, `run_region_resolved`, `run_kestrel_counting_mode`,
  and `run_advntr_model_*` when adVNTR ran). Nested keys flatten with `_`; arrays join with `; `;
  nulls are blank cells.
- **Step record**: `step`, `start`, `end`, `command`, `result_file`, `file_type`,
  `md5sum`, `result_file_missing` (`True` or `False`).
- **Parsed result (`parsed_result_*`)**: Single-row outputs expand to
  `parsed_result_data_<column>` cells; multi-row counts are recorded in
  `parsed_result_n_rows`; `#` comment headers join with pipe delimiters in
  `parsed_result_comments`; JSON sub-objects flatten with `_`.

`pipeline_summary_rows.csv` / `pipeline_summary_rows.tsv` records long-form tabular
results: one row per step, result row, and field, with columns `step`, `row_index` (0-based),
`field`, and `value`. This table captures multi-row adVNTR and cross-match records.

## Cohort Output Files

When executing `vntyper cohort`, the following summary files are generated in the output directory:

| File | Format | Description |
|------|--------|-------------|
| `cohort_summary.html` | HTML | Standalone summary report with cohort distributions and interactive tables |
| `cohort_kestrel.<csv\|tsv\|json>` | CSV/TSV/JSON | Aggregated Kestrel variant calls across all cohort samples |
| `cohort_advntr.<csv\|tsv\|json>` | CSV/TSV/JSON | Aggregated adVNTR variant calls across all cohort samples (if adVNTR present) |
| `cohort_stats.<csv\|tsv\|json>` | CSV/TSV/JSON | Aggregated execution parameters, reference files, and coverage statistics |
| `cohort_call_frequency.<csv\|tsv\|json>` | CSV/TSV/JSON | Grouped variant call frequency table across the cohort |
| `pseudonymization_table.tsv` | TSV | Sample name to pseudonym mapping table (when `--pseudonymize-samples` is used) |

### cohort_call_frequency Columns

The `cohort_call_frequency.<csv|tsv|json>` export contains exactly 14 columns detailing variant call distribution:

| Column | Description |
|--------|-------------|
| `Grouping_Key` | Molecular identity or caller representation used to group the call |
| `Grouping_Key_Kind` | Either `molecular-identity` or `caller-representation` |
| `Molecular_Identity` | Canonical molecular identity serialization, or empty if unresolved |
| `Motifs` | Repeat motif identifier(s) from caller |
| `POS` | Locus position |
| `REF` | Reference sequence |
| `ALT` | Alternate sequence |
| `Variant` | Variant type (e.g., insertion, deletion, SNV) |
| `Sample_Count` | Number of distinct cohort samples carrying this variant call |
| `Frequency` | Call frequency computed as `Sample_Count / cohort_size` |
| `Below_Cutoff` | `yes` if `Frequency <= rare_allele_max_frequency`, `no` otherwise |
| `Samples` | Semicolon-separated list of sample names (or pseudonyms) carrying the call |
| `Min_Depth_Score` | Minimum numeric depth score observed among samples with this call |
| `Max_Depth_Score` | Maximum numeric depth score observed among samples with this call |
