# Configuration

VNtyper 2 separates runtime configuration from decision profiles. To override the main runtime configuration:

```bash
vntyper --config-path /path/to/custom/config.json pipeline \
    --bam inputs/sample.bam -o results/sample/
```

`--config-path` replaces the main `config.json`; it does not select or overlay a decision profile. The Kestrel and optional-module sidecars retain excluded runtime values such as memory, executable arguments, and reference paths.

## Decision Profiles

Every pipeline run resolves one complete profile. With no option, VNtyper verifies and uses its packaged `vntyper/profiles/decision_profile.json`; this keeps all package defaults unchanged. To select another profile:

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ \
    --decision-profile /path/to/complete-profile.json
```

An explicit profile is a complete standalone file, not a patch. It must contain every packaged decision leaf and use `profile_kind: "explicit-custom"` or `profile_kind: "generated"`. Missing or unknown keys, altered field types or comparison semantics, and any change to a `fixed-safety` field fail before stage artifacts are written. This includes the `0.00469` reporting floor and independent GG gate, `0.00515`, alternate k-mer-path boundaries `20`, `21`, and `100`, active-region boundary `200`, BAM flank `8`, thin haplotype-record support `3`, and both source-specific tier-A support values `5`.

Generated profiles do not auto-activate. Creating or placing one in an output directory changes no run; an operator must pass that exact complete file to `--decision-profile`. The adVNTR model identity, fetch-window compatibility, and minimum binary version from Issue 268 are deliberately outside the profile and cannot be bypassed by one. Runtime paths, references, coverage presentation, input routing, and report wording are excluded for the same reason.

The resolved canonical bytes are snapshotted under `provenance/`, and schema-3 `pipeline_summary.json` records their ID, revision, kind, source, and SHA-256. Reports and cohorts verify that run-local snapshot rather than inferring policy from the package installed later.

The optional [`vntyper calibrate`](../cli/calibrate.md) workflow can generate such a profile from role-bound retained evidence. It changes only six declared dominance and whole-locus-abstention leaves; it does not change any shipped cutoff or default. Its output remains inactive until explicitly supplied to `pipeline --decision-profile`. The available simulations are development evidence, not independent external validation, so Issue #295 remains blocked. Reporting an interval is not a clinical safety claim.

!!! warning "Fixed safety thresholds"
    The shipped depth score and confidence boundaries are empirically validated values from Saei et al., iScience 26, 107171 (2023). They are recorded in the complete profile for reproducibility but cannot be changed by an explicit or generated profile.

## Main Configuration (config.json)

Controls tool paths, reference data, processing parameters, and quality thresholds.

```json
{
  "default_values": {
    "threads": 4,
    "output_dir": "out",
    "output_name": "output",
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
    "q30_rate": 0.7,
    "passed_filter_reads_rate": 0.8
  }
}
```

### Key Sections

| Section | Purpose |
|---------|---------|
| `default_values` | Fallback values when CLI arguments are not provided |
| `reference_data` | Paths to BWA indexes, MUC1 motif references, adVNTR databases, and SHARK MUC1 region FASTAs. The snippet above shows three keys. `vntyper install-references` writes one `bwa_reference_*` per installed genome plus shared adVNTR and SHARK keys (such as `bwa_reference_GRCh38`, `advntr_reference_vntr_hg38`, `muc1_region_fasta_hg38`). See [Reference Assemblies](reference-assemblies.md) for full key mappings and fallback logic |
| `tools` | Executable paths for external binaries |
| `bam_processing` | fastp QC parameters and assembly-specific region coordinates |
| `thresholds` | Quality thresholds for coverage and read quality metrics |
| `api` | Base URL for online mode API (`https://vntyper.org/api`) |

#### Coverage thresholds

`mean_vntr_coverage` and `percent_vntr_uncovered` decide the report coverage QC verdict, and through it the `quality_metrics_pass` axis of the screening summary. The mean fails strictly below its threshold; the uncovered fraction fails strictly above its own, so a sample at exactly 100x and exactly 50.0% uncovered passes both. An unmeasured metric (such as a run omitting coverage calculation) does not fail the gate.

The report evaluates verdicts on rounded figures (two decimal places), preventing false failures beside matching boundaries.

Since VNtyper 2.0.8 both keys are enforced. Previously, `percent_vntr_uncovered` drove only a color-coded indicator.

#### fastp report thresholds

The fastp rows in the HTML report require all four `thresholds` keys:
`duplication_rate`, `q20_rate`, `q30_rate`, and `passed_filter_reads_rate`.
Each is a finite numeric fraction from 0 through 1, inclusive. The report derives
both the displayed cutoff label and the status decision from that same configured
value. Each measured rate and configured cutoff is rounded half-up on its exact
decimal value to two decimal places of percent before display and comparison.
JSON threshold fractions enter that decimal boundary directly, without an
intermediate binary-float round trip. The configuration object remains
dict-compatible and retains its existing numeric values for other consumers.

`duplication_rate` warns strictly above its cutoff; the three rate metrics warn
strictly below theirs. A measured value exactly at its configured cutoff is OK.
Missing keys, strings, booleans, non-finite values, or fractions outside 0 to 1 are
logged and raise `ValueError` during rendering rather than substituting a
default cutoff.

## Packaged Kestrel decision projection

The packaged profile projects the following Kestrel decisions. `kestrel_config.json`
retains excluded Kestrel runtime settings; pipeline decision consumers use the one
resolved profile instead of reading mutable module globals.

```json
{
  "kestrel_settings": {
    "java_memory": "12g",
    "kmer_sizes": [20],
    "max_align_states": 40,
    "max_hap_states": 40
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
    "exclude_motifs_right": ["8", "9", "7", "6p", "6"],
    "exclude_alts_combined": ["CCGCC", "CGGCG", "CGGCC"]
  }
}
```

### Kestrel Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `java_memory` | `12g` | JVM heap size for Kestrel |
| `kmer_sizes` | `[20]` | K-mer sizes to evaluate (pipeline stops after first success) |
| `max_align_states` | `40` | Maximum alignment states in Kestrel |
| `max_hap_states` | `40` | Maximum haplotype states in Kestrel |

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
| `flagging_rules` | two rules | Named conditions annotating a call `Flag` column. Advisory by default: the call remains reported. |
| `artifact_flags` | `["False_Positive_4bp_Insertion"]` | Flag names identifying a technical artifact rather than a true variant. A row carrying an artifact flag is excluded from `kestrel_result.tsv`. |
| `duplicate_flagging.enabled` | `false` | Flag lower-priority calls sharing `REF`/`ALT` as `Potential_Duplicate`. |

A name listed in `artifact_flags` must also be evaluated by a `flagging_rules` entry. Excluded rows remain recorded in `kestrel_pre_result.tsv` with `flag_filter_pass = False`.

Each flag rule starts with one boolean node: `all` and `any` require non-empty child lists; `not` takes one child. Children can be nested boolean nodes or predicates, capped at 32 boolean nodes. Predicates contain `left`, `operator`, and `right`. Operands specify one `column` or `literal`. Supported operators are `eq`, `lt`, `in`, and `casefold_eq`. Columns are validated against the declared result schema before execution, rejecting malformed rules across all samples. Configured numbers must be finite; NaN is treated as null only when originating from a row.

Executable syntax, lambdas, regexes, and imports are prohibited. Boolean logic is restricted to structured `all`, `any`, and `not` nodes. The legacy string:

```json
"(REF == 'C') and (ALT == 'CGGCA')"
```

is accepted only for backward compatibility with `False_Positive_4bp_Insertion` and maps directly to the structured object shown above. Custom expression strings are rejected. See [Variant Flagging](../pipeline/flagging.md#rule-schema-and-custom-rules) for the complete schema.

!!! note "Emptying `artifact_flags` restores earlier behavior"
    Setting `"artifact_flags": []` makes every flag advisory, matching VNtyper behavior prior to Issue #174 without code edits.

See [Variant Flagging](../pipeline/flagging.md) for full descriptions of both flag classes.
