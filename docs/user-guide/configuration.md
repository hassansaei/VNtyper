# Configuration

VNtyper 2 uses two configuration files. To override the default configuration:

```bash
vntyper --config-path /path/to/custom/config.json pipeline --bam sample.bam -o results/
```

!!! warning "Modifying Thresholds"
    The depth score and confidence thresholds in `kestrel_config.json` are empirically validated values from Saei et al., iScience 26, 107171 (2023). Changing them may affect genotyping accuracy.

## Main Configuration (config.json)

Controls tool paths, reference data, processing parameters, and quality thresholds.

```json
{
  "default_values": {
    "threads": 4,
    "output_dir": "out",
    "output_name": "processed",
    "archive_format": "zip",
    "reference_assembly": "hg19"
  },
  "reference_data": {
    "muc1_reference_vntr": "reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
    "bwa_reference_hg19": "reference/alignment/chr1.hg19.fa",
    "bwa_reference_hg38": "reference/alignment/chr1.hg38.fa"
  },
  "tools": {
    "fastp": "fastp",
    "samtools": "samtools",
    "bwa": "bwa",
    "kestrel": "vntyper/dependencies/kestrel/kestrel.jar",
    "java_path": "java"
  },
  "bam_processing": {
    "compression_level": 6,
    "deduplication": true,
    "length_required": 50,
    "qualified_quality_phred": 20
  },
  "thresholds": {
    "mean_vntr_coverage": 100,
    "percent_vntr_uncovered": 50.0,
    "duplication_rate": 0.1,
    "q20_rate": 0.8,
    "q30_rate": 0.7
  }
}
```

### Key Sections

| Section | Purpose |
|---------|---------|
| `default_values` | Fallback values when CLI arguments are not provided |
| `reference_data` | Paths to BWA indexes, MUC1 motif references, adVNTR databases |
| `tools` | Executable paths for external tools |
| `bam_processing` | fastp QC parameters and assembly-specific region coordinates |
| `thresholds` | Quality thresholds for coverage and read quality metrics |
| `api` | Base URL for the online mode API (`https://vntyper.org/api`) |

## Kestrel Configuration (kestrel_config.json)

Controls Kestrel execution and the entire postprocessing pipeline (scoring, confidence, motif filtering, flagging).

```json
{
  "kestrel_settings": {
    "java_memory": "12g",
    "kmer_sizes": [20],
    "max_align_states": 60,
    "max_hap_states": 60
  },
  "confidence_assignment": {
    "depth_score_thresholds": {
      "low": 0.00469,
      "high": 0.00515
    },
    "alt_depth_thresholds": {
      "low": 20,
      "mid_low": 21,
      "mid_high": 100
    }
  },
  "flagging_rules": {
    "False_Positive_4bp_Insertion": "(REF == 'C') and (ALT == 'CGGCA')",
    "Low_Depth_Conserved_Motifs": "(Depth_Score < 0.4) and (Motif in ['1', '2', '3', '4', '6', '7', '8', '9'])"
  },
  "motif_filtering": {
    "exclude_motifs_right": ["Q", "8", "9", "7", "6p", "6", "V", "J", "I", "G", "E", "A"],
    "exclude_alts_combined": ["CCGCC", "CGGCG", "CGGCC"]
  }
}
```

### Kestrel Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `java_memory` | `12g` | JVM heap size for Kestrel |
| `kmer_sizes` | `[20]` | K-mer sizes to try (pipeline stops after first success) |
| `max_align_states` | `60` | Maximum alignment states in Kestrel |
| `max_hap_states` | `60` | Maximum haplotype states in Kestrel |

### Confidence Thresholds

| Threshold | Value | Used For |
|-----------|-------|----------|
| `depth_score_thresholds.low` | 0.00469 | Minimum depth score for Low_Precision |
| `depth_score_thresholds.high` | 0.00515 | Minimum depth score for High_Precision |
| `alt_depth_thresholds.low` | 20 | Minimum alternate depth for any positive call |
| `alt_depth_thresholds.mid_high` | 100 | Threshold for High_Precision* designation |
