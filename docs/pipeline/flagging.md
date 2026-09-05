# Variant Flagging

Flagging applies configurable, post-hoc empirical rules to variant calls. While confidence scoring estimates signal strength, flagging annotates calls with recurring error profiles observed in Kestrel and adVNTR outputs.

There are **two classes of flag**:

| Class | Declared in | Effect on reported call |
|-------|-------------|-------------------------|
| **Advisory** | `flagging_rules` only | The variant is reported. The flag lowers selection ranking and appears in the `Flag` column. |
| **Artifact** | `flagging_rules` **and** `artifact_flags` | The variant is **excluded** from `kestrel_result.tsv`. It is retained in `kestrel_pre_result.tsv` with `flag_filter_pass = False`. |

An advisory flag highlights a call that may be valid but warrants scrutiny. An artifact flag identifies a technical false positive; reporting it as a call would produce a false finding.

!!! warning "Artifact flags exclude rather than deprioritize (Issue #174)"
    Prior to Issue #174, all flags were advisory. A sample whose sole variant was the 4-bp insertion artifact produced a `High_Precision` row with a flag, which mapped to `High_Precision_flagged` in `report_config.json` and was interpreted by `is_finding` as a positive call. Artifact flags ensure false positives are excluded prior to selection.

## How Flagging Works

Flagging rules are defined as structured boolean trees in `kestrel_config.json`. VNtyper compiles and validates the entire rule set against the result schema prior to processing rows. If the tree evaluates to `true`, the flag name is added to the `Flag` column. Multiple flags are comma-separated in configuration order. Variants matching no rules receive `Flag = "Not flagged"`.

!!! info "Flagging precedes variant selection"
    Flagging executes before final variant selection (Issue #145). When multiple candidates pass all filters, unflagged variants take precedence over flagged variants.

## Default Flagging Rules

Configured in `kestrel_config.json`:

### False_Positive_4bp_Insertion

```json
{
  "all": [
    {"left": {"column": "REF"}, "operator": "eq", "right": {"literal": "C"}},
    {"left": {"column": "ALT"}, "operator": "eq", "right": {"literal": "CGGCA"}}
  ]
}
```

Flags a specific 4-bp insertion (C > CGGCA) that recurs in GC-rich VNTR regions due to k-mer graph ambiguity.

This is the only default flag declared in `artifact_flags`. A row carrying it is excluded from `kestrel_result.tsv`. A sample with only this call reports a negative result.

### Low_Depth_Conserved_Motifs

```json
{
  "all": [
    {"left": {"column": "Depth_Score"}, "operator": "lt", "right": {"literal": 0.4}},
    {
      "left": {"column": "Motif"},
      "operator": "in",
      "right": {"literal": ["1", "2", "3", "4", "6", "7", "8", "9"]}
    }
  ]
}
```

Flags variants in conserved repeat unit motifs (motifs 1 to 9) when the depth score is below 0.4. Because these motifs are conserved across alleles, low-depth calls in these positions represent probable artifacts.

This flag is **advisory**: it deprioritizes candidates during selection but does not exclude them.

## Artifact Flags

`artifact_flags` lists the flag names that disqualify a candidate:

```json
{
  "artifact_flags": ["False_Positive_4bp_Insertion"]
}
```

A flag must appear in both `flagging_rules` (to be evaluated) and `artifact_flags` (to exclude). Names in `artifact_flags` that cannot be triggered are inert.

### Exclusion mechanics

The pipeline derives `flag_filter_pass`, which is `False` whenever the `Flag` column contains an artifact flag:

- Membership is evaluated per token in comma-separated flag strings.
- A row carrying both an artifact flag and an advisory flag is excluded.
- Rows with no flag pass the gate.

`flag_filter_pass` is evaluated alongside the other five structural gates in `FILTER_COLUMNS`. Excluded rows are retained in `kestrel_pre_result.tsv` with `flag_filter_pass = False`.

To make all flags advisory, set `"artifact_flags": []`.

## Duplicate Flagging

Groups calls by `REF` and `ALT`, sorts by `Depth_Score` descending, and flags secondary calls as `Potential_Duplicate`:

```json
{
  "enabled": false,
  "flag_name": "Potential_Duplicate",
  "group_by": ["REF", "ALT"],
  "sort_by": [
    {"column": "Depth_Score", "ascending": false}
  ]
}
```

Sorting keys use `Depth_Score` only (Issue #197). Duplicate flagging is disabled by default in shipped configurations.

## Rule Schema and Custom Rules

To add a custom rule:

```json
{
  "flagging_rules": {
    "My_Custom_Flag": {
      "all": [
        {"left": {"column": "Depth_Score"}, "operator": "lt", "right": {"literal": 0.005}},
        {"left": {"column": "Variant"}, "operator": "eq", "right": {"literal": "Insertion"}}
      ]
    }
  }
}
```

Rules start with a boolean operator (`all`, `any`, or `not`). Nesting is limited to 32 nodes. Predicates specify `left`, `operator`, and `right`. Operands reference a declared `column` or a `literal`. Supported operators:

| Operator | Definition |
|----------|------------|
| `eq` | Strict equality. Booleans compare only with booleans. |
| `lt` | Numeric less-than comparison. |
| `in` | Scalar membership within a homogeneous literal list. |
| `casefold_eq` | Case-insensitive string equality. |

Rules are inert JSON data; arbitrary expressions, lambdas, and arithmetic are prohibited.

## adVNTR Flagging Rules

Defined in `vntyper/modules/advntr/advntr_config.json`:

| Rule | Condition | Definition |
|------|-----------|------------|
| `Low_Coverage` | `NumberOfSupportingReads < 10` | Low supporting read count. |
| `Repeat_Unit_7` | `RU == '7'` | Call in repeat unit 7, a known recurring artifact. |
| `Polymorphic_Call` | `Variant in [...24 states...]` | adVNTR state matches recurrent non-specific background. |

adVNTR flags are advisory: they annotate rows and map to `Positive (Flagged)` without withdrawing the call.

### The `Polymorphic_Call` list

Derived from historical runs across the `renome` cohort, these 24 states represent recurrent events that do not establish molecular identity. Canonical state signatures are tracked in `vntyper/modules/advntr/advntr_artifact_evidence.json` with an accompanying SHA-256 digest.

Positive adVNTR calls record `Evidence_Disposition`:
- Matching an active state marks the call `identity-insufficient`. The call remains visible but cannot establish cross-caller concordance or promote a tier-A nomenclature call.
- Non-matching calls are marked `admissible`.

### Run provenance

Runs record the canonical evidence snapshot at `provenance/advntr_artifact_evidence.json` and log its SHA-256 digest in `advntr_evidence_digest`.
