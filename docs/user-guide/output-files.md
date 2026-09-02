# Output Files

## Directory Structure

```
results/
├── pipeline_summary.json        # Machine-readable pipeline summary
├── pipeline_summary.csv         # Optional (--summary-formats csv)
├── pipeline_summary.tsv         # Optional (--summary-formats tsv)
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
│   └── pipeline_info.json       # BAM header metadata (BAM/CRAM input)
├── alignment_processing/
│   └── output_sorted.bam        # BWA-aligned BAM (FASTQ input only)
├── coverage/
│   └── coverage_summary.tsv     # VNTR region coverage statistics
└── advntr/                      # Only when --extra-modules advntr
    ├── output_adVNTR.tsv         # Raw adVNTR output
    ├── output_adVNTR_result.tsv  # Processed adVNTR result
    └── cross_match_results.tsv   # Kestrel vs adVNTR comparison
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
| `Nomenclature_Tier` | `A`, `B` or `C` — how much of the name may be stated |
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
what each caller reported, so a disagreement stays legible in either result file. Both
files carry the same values for all three, which is what lets a cohort table merging
them be read as one.

`Ambiguity_Interval`, `Repeat_Form` and `Nomenclature_adVNTR` are nullable: they are
written only where they mean something, and are left empty rather than padded with a
placeholder.

A negative run writes a different, narrower schema (first column `Motif`, no depth or
flag columns) and carries none of these — there is no variant, so there is no name or
molecular identity row.

Positive adVNTR result rows append the same four molecular-identity columns in the same
order. Negative adVNTR output retains its existing narrower schema. The Kestrel result
additionally carries `__Identity_*` capture and selection cells, and positive rows in both
caller result files carry `__Reconciled_Molecular_Identity`: the canonical serialization of
the whole-locus identity the reconciler selected, or an empty cell after abstention. The
double-underscore prefix marks internal identity-capture and calibration-replay
persistence, not public fields: they are copied into `pipeline_summary.json` with the
other result columns, excluded from the policy projection, not rendered in the HTML
report, and never written to a negative row. The per-sample HTML
report, `pipeline_summary.json`, and cohort HTML/TSV/CSV/JSON exports carry all four
recorded values; they never infer them from `POS`, `REF`, `ALT`, `Variant` or
`Nomenclature`. In a sample HTML report, use the table's column control to show the
quartet on screen. The printed table folds these width-heavy fields and the labelled
per-row appendix prints them in the same order, including the full molecular identity.

For compatibility, a schema-1 or schema-2 summary row missing any one of the four fields
renders `legacy identity not recorded` in all four downstream identity cells. A complete
quartet is copied exactly, including an empty unresolved identity, integer `0`
representation count, and a nonzero hypothesis count.

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
[authoritative flag table](../pipeline/nomenclature.md#flags) for their source-specific
meanings.

## Confidence Levels

Confidence is assigned based on empirically validated depth score thresholds from Saei et al., iScience 26, 107171 (2023).

| Level | Meaning |
|-------|---------|
| **High_Precision\*** | Depth score >= 0.00515 and alternate depth >= 100 |
| **High_Precision** | Depth score >= 0.00515, alternate depth 21 to <100, and region depth > 200 |
| **Low_Precision** | Variant detected with marginal depth or depth score support |
| **Negative** | No variant passed filtering thresholds |

!!! tip
    A result with confidence `Negative` means no MUC1-VNTR frameshift variant was detected -- it does not necessarily mean the sample is truly negative.

## How summary_report.html Displays Numbers

Every numeric column in the per-sample HTML report is now formatted server-side and
written into the file by VNtyper.
Until the fix for issue #242 the report shipped a small script that rewrote every numeric
cell of every table with `toFixed(4)` in the reader's browser, so the number on
screen depended on whether three content delivery networks were reachable -- the archived
file said one thing to a reader online and another to a reader offline.

Removing that script changes the printed form of some columns. **No value changed; only
its rendering did.** The `kestrel_result.tsv`, `pipeline_summary.json` and adVNTR TSV
outputs are untouched, and remain the source to parse. If you scrape the HTML report,
these are the columns whose text differs.

| Table | Column | Displayed before | Displayed now | Why |
|-------|--------|------------------|---------------|-----|
| Kestrel | `POS`, `Estimated_Depth_AlternateVariant`, `Estimated_Depth_Variant_ActiveRegion` | `67`, `120`, `12000` | unchanged | Whole numbers, no decimals |
| Kestrel | `Depth_Score` | `0.01` | `0.010012` | Six decimal places. The confidence calibration is stated to five (0.00469 and 0.00515), so four was coarser than the thresholds the value is judged against, and a score of 0.00001234 printed as `0` |
| adVNTR | `VID`, `NumberOfSupportingReads`, `POS` | `25561`, `14`, `67` | unchanged | Whole numbers, no decimals |
| adVNTR | `MeanCoverage` | `98.5`, `40` | `98.50`, `40.00` | Always two decimal places, so every mean states the same precision |
| adVNTR | `Pvalue` | `0`, `0.0001` | `1e-09`, `0.000123` | Three significant figures. The old `toFixed(4)` displayed `1e-9` as `0` |

### The Kestrel table's column order changed

`Motif_sequence` is now the **last** column of the per-sample report's Kestrel table,
after the confidence, flag and nomenclature fields. It used to be sixth, between `ALT`
and `Estimated_Depth_AlternateVariant`.

Before the field was corrected to the selected 60 bp half, it held the 120 bp pair record.
In its old sixth position, that long unbroken value pushed `Confidence` and `Flag` off the
right edge of a 1280px screen, so the two columns a reader opens the report for were the
ones that scrolled out of sight while the widest and least-scanned column sat in the
middle. Last, it is the column that scrolls -- which is the correct one to lose.

**The column-order migration is display-only.** Separately, the value correction above
changes `Motif_sequence` from the pair record to its named half in result TSVs and reports;
no other field changes. If you scrape the HTML report by column *position*, key on the
heading text instead. The `cohort_summary.html` Kestrel table declares its own column
order and is **not** affected by the display-order migration.

Further display changes in the same release:

- **The adVNTR table's headings are English**, matching the Kestrel table above it:
  `NumberOfSupportingReads` is `Supporting Reads`, `Pvalue` is `P-value`, `RU` is
  `Repeat Unit`, and the eight naming fields take the same headings the Kestrel table
  already used -- `Nomenclature` is `MUC1 Name`, `Nomenclature_Tier` is `Tier`,
  `Ambiguity_Interval` is `Ambiguity`, and so on. The two tables previously named the
  same field two different ways in one document. `advntr_result.tsv` is unchanged and
  keeps the source names; if you scrape the HTML, key on the new heading text.
- **A semicolon-separated list is spaced.** `Nomenclature_Flags` renders as
  `known-variant; motif-context-diverges; position-ambiguous` rather than with bare
  semicolons, so a cell narrow enough to wrap breaks between flags instead of inside
  one. The separator is unchanged; splitting on `;` and stripping is unaffected.
- **Numbers are aligned on their own axis** and set with tabular figures, and each
  column heading carries a one-line explanation as hover text. Every coded value the
  tables print -- the tier letter and each nomenclature flag -- is also spelled out in
  words in a **reading key printed underneath the tables**, so the explanation is on
  paper and needs no pointer.
- **The `Flag` column shows the reason in words**, beside a tick or a cross, instead of a
  glyph whose reason appeared only in a hover tooltip. The reason is therefore in the
  printed page, readable by a screen reader, and present when scripts do not run.
- **Flagged variant rows are never hidden.** The per-sample report's switch is now
  *Highlight flagged values*: it changes emphasis and removes nothing. A row-count line
  above each table states how many rows are shown out of how many exist. The
  [cohort report](cohort-analysis.md) keeps its show/hide filter, where hiding flagged
  rows across many samples is a triage aid rather than a way to hide one sample's evidence.
  Note that the cohort report still rounds its numbers in the browser; that is tracked
  separately.

## summary_report.html Artifact Modes

`--report-igv embedded` (the default) produces one self-contained report. `sidecar`
writes a small summary plus a self-contained `igv_report.html`, and `off` omits the
alignment browser without removing any table row. The last verification specimen was
**78,486 bytes without alignment data** and **575,762 bytes with embedded alignment
data**, compared with **2,002,405 bytes** fetched by the previous report's 11 CDN tags.
These are measured specimen sizes, not fixed limits; the tables and other sample content
also contribute to the file.

Embedded alignment viewing uses the browser's `DecompressionStream` support, with a
2023 cross-browser floor of Chrome 80+, Safari 16.4+ or Firefox 113+. On an older
browser, the alignment panel explains the limitation and the complete variant tables
remain readable.

## Pipeline Summary JSON

The `pipeline_summary.json` file records each pipeline step with timestamps, output paths, and parsed results. It is used by the [cohort analysis](cohort-analysis.md) module to aggregate results across samples.

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

Molecular-identity publication is additive to summary schema 2; it does not introduce
schema 3. Current summaries record the packaged selection policy as
`decision_policy: legacy-selection-v1`. Older summaries without that provenance remain
readable and are labelled as not recorded rather than being assigned the current policy.
