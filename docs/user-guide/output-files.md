# Output Files

## Directory Structure

```
results/
├── pipeline_summary.json        # Machine-readable pipeline summary
├── pipeline_summary.csv         # Optional (--summary-formats csv)
├── pipeline_summary.tsv         # Optional (--summary-formats tsv)
├── pipeline.log                 # Pipeline execution log
├── summary_report.html          # HTML report with IGV visualization
├── predefined_regions_<assembly>.bed  # Region BED file (e.g., hg19, hg38)
├── kestrel/
│   ├── kestrel_result.tsv       # Final genotyping result
│   ├── kestrel_pre_result.tsv   # Pre-filter variants (all candidates)
│   ├── output.vcf               # Raw Kestrel VCF
│   ├── output_indel.vcf         # Filtered INDEL VCF
│   ├── output_indel.vcf.gz      # Compressed INDEL VCF (if bcftools available)
│   ├── output.bam               # Kestrel alignment BAM
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
| `Motif_sequence` | Nucleotide sequence of the motif |
| `Estimated_Depth_AlternateVariant` | Read depth supporting the alternate allele |
| `Estimated_Depth_Variant_ActiveRegion` | Total read depth in the variant active region |
| `Depth_Score` | Ratio of alternate depth to active region depth |
| `Confidence` | Confidence classification (see below) |
| `Flag` | Quality flag (`Not flagged` or a flag reason) |
| `haplo_count` | Number of haplotype calls supporting the same variant |

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

Every figure in the per-sample HTML report is now written into the file by VNtyper.
Until the fix for issue #242 the report shipped a small script that rewrote every numeric
cell of every table to four decimal places in the reader's browser, so the number on
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
| adVNTR | `Pvalue` | `0`, `0.0001` | `1e-09`, `0.000123` | Three significant figures. Rounding to four decimals displayed a highly significant p-value as zero |

### The Kestrel table's column order changed

`Motif_sequence` is now the **last** column of the per-sample report's Kestrel table,
after `Flag`. It used to be sixth, between `ALT` and `Estimated_Depth_AlternateVariant`.

The real motif sequence is 121 bp. Sixth, it pushed `Confidence` and `Flag` off the right
edge of a 1280px screen, so the two columns a reader opens the report for were the ones
that scrolled out of sight while the widest and least-scanned column sat in the middle.
Last, it is the column that scrolls -- which is the correct one to lose.

**This is a display order only.** `kestrel_result.tsv` is unchanged, and so is every
other output; if you parse anything, parse those. If you scrape the HTML report by column
*position*, this is a breaking change and you should key on the heading text instead. The
`cohort_summary.html` Kestrel table declares its own column order and is **not** affected.

Two further display changes in the same release:

- **The `Flag` column shows the reason in words**, beside a tick or a cross, instead of a
  glyph whose reason appeared only in a hover tooltip. The reason is therefore in the
  printed page, readable by a screen reader, and present when scripts do not run.
- **Flagged variant rows are never hidden.** The per-sample report's switch is now
  *Highlight flagged values*: it changes emphasis and removes nothing. A row-count line
  above each table states how many rows are shown out of how many exist. The
  [cohort report](cohort-analysis.md) keeps its show/hide filter, where hiding flagged
  rows across many samples is a triage aid rather than a way to lose a patient's result.
  Note that the cohort report still rounds its numbers in the browser; that is tracked
  separately.

## Pipeline Summary JSON

The `pipeline_summary.json` file records each pipeline step with timestamps, output paths, and parsed results. It is used by the [cohort analysis](cohort-analysis.md) module to aggregate results across samples.

Key fields:

```json
{
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
