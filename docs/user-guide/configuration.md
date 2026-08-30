# Configuration

VNtyper 2 uses two configuration files. To override the default configuration:

```bash
vntyper --config-path /path/to/custom/config.json pipeline \
    --bam inputs/sample.bam -o results/sample/
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
| `reference_data` | Paths to BWA indexes, MUC1 motif references, adVNTR databases, and SHARK's MUC1 region FASTAs. The snippet above shows only three keys; `vntyper install-references` writes one `bwa_reference_*` per genome it installed plus the shared adVNTR and SHARK keys (e.g. `bwa_reference_GRCh38`, `advntr_reference_vntr_hg38`, `muc1_region_fasta_hg38`) -- see [Reference Assemblies](reference-assemblies.md) for the complete key list and how a missing key falls back |
| `tools` | Executable paths for external tools |
| `bam_processing` | fastp QC parameters and assembly-specific region coordinates |
| `thresholds` | Quality thresholds for coverage and read quality metrics |
| `api` | Base URL for the online mode API (`https://vntyper.org/api`) |

#### Coverage thresholds

`mean_vntr_coverage` and `percent_vntr_uncovered` both decide the report's coverage QC verdict, and through it the `quality_metrics_pass` axis of the screening summary. The mean fails strictly *below* its threshold; the uncovered fraction fails strictly *above* its own, so a sample at exactly 100x and exactly 50.0% uncovered passes both. A metric that was never measured -- a run with no Coverage Calculation step -- does not fail the gate.

The verdict is evaluated on the figures the report displays, which are rounded to two decimal places, so no report prints `FAIL` beside a mean of `100.00` and a threshold of 100.

Since VNtyper 2.0.8 both keys are enforced. Before it, `percent_vntr_uncovered` had a threshold, was computed on every run and drove only a color-coded icon.

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
    "False_Positive_4bp_Insertion": {
      "all": [
        {"left": {"column": "REF"}, "operator": "eq", "right": {"literal": "C"}},
        {"left": {"column": "ALT"}, "operator": "eq", "right": {"literal": "CGGCA"}}
      ]
    },
    "Low_Depth_Conserved_Motifs": {
      "all": [
        {"left": {"column": "Depth_Score"}, "operator": "lt", "right": {"literal": 0.4}},
        {
          "left": {"column": "Motif"},
          "operator": "in",
          "right": {"literal": ["1", "2", "3", "4", "6", "7", "8", "9"]}
        }
      ]
    }
  },
  "artifact_flags": ["False_Positive_4bp_Insertion"],
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

### Flagging and Artifact Exclusion

| Key | Default | Description |
|-----|---------|-------------|
| `flagging_rules` | two rules | Named conditions that annotate a call's `Flag` column. Advisory by default: the call is still reported. |
| `artifact_flags` | `["False_Positive_4bp_Insertion"]` | Which of those flag names identify a **technical artifact** rather than a call. A row carrying one is excluded from `kestrel_result.tsv`. |
| `duplicate_flagging.enabled` | `false` | Whether to flag lower-priority calls sharing a `REF`/`ALT` as `Potential_Duplicate`. |

A name listed in `artifact_flags` must also be raised by a `flagging_rules` entry, or it never matches anything. Excluded rows are not lost: they remain in `kestrel_pre_result.tsv` with `flag_filter_pass = False`.

Each flag rule starts with one explicit boolean node: `all` and `any` take non-empty child lists, while `not` takes one child. Children can be nested boolean nodes or predicates, with boolean nesting capped at 32 nodes. Predicates contain exactly `left`, `operator`, and `right`; operands contain one `column` or `literal`. Available operators are `eq`, `lt`, `in`, and `casefold_eq`. Columns are validated against the frame available to the consumer before any row is processed. Null values make a predicate false, and VNtyper never coerces strings, numbers, or booleans between types. Configured numbers and non-null numbers read from a row must be finite; NaN is null only when it comes from a row, while infinities always abort.

Calls, attributes, indexing, imports, comprehensions, lambdas, regular expressions, arithmetic, and executable boolean syntax are not part of the schema and are never executed. Boolean logic is represented only by the bounded `all`, `any`, and `not` JSON nodes. The byte-exact last-release string

```json
"(REF == 'C') and (ALT == 'CGGCA')"
```

is accepted only as a migration input for `False_Positive_4bp_Insertion` and maps to the structured object shown above. Modified whitespace, added clauses, custom expression strings, and another flag's historical expression are rejected. New configuration must use structured rules. See [Variant Flagging](../pipeline/flagging.md#rule-schema-and-adding-custom-rules) for the complete schema and examples.

!!! note "Emptying `artifact_flags` restores the previous behaviour"
    Setting `"artifact_flags": []` makes every flag advisory again --- exactly how
    VNtyper behaved before the Issue #174 fix --- with no code change. The artifact
    decision lives in this file, so narrowing or withdrawing it is a configuration edit.

See [Variant Flagging](../pipeline/flagging.md) for the full description of the two flag classes.
