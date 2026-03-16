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
| `Variant` | Variant type (`insertion` or `deletion`) |
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

## Pipeline Summary JSON

The `pipeline_summary.json` file records each pipeline step with timestamps, output paths, and parsed results. It is used by the [cohort analysis](cohort-analysis.md) module to aggregate results across samples.

Key fields:

```json
{
  "version": "2.0.1",
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
