# Variant Flagging

Flagging applies configurable, post-hoc empirical filters to variant calls. Unlike the scoring and confidence system (which determines whether a variant is real), flagging identifies calls that may be technically valid but warrant additional scrutiny due to known artifact patterns.

## How Flagging Works

Each flagging rule is defined as a named condition in `kestrel_config.json`. Rules are Python logical expressions evaluated against each row of the variant DataFrame. If a rule's condition evaluates to `True`, the corresponding flag name is appended to the variant's `Flag` column. Multiple flags are comma-separated. Variants matching no rules receive `Flag = "Not flagged"`.

!!! info "Flagging occurs before variant selection"
    As of VNtyper 2 (Issue #145 fix), flagging is applied **before** the final variant selection step. This ensures that when multiple candidate variants pass all filters, unflagged variants are preferred over flagged ones. Previously, a flagged variant could be selected as the best call because flags were added after selection.

## Current Flagging Rules

The default rules in `kestrel_config.json`:

### False_Positive_4bp_Insertion

```
(REF == 'C') and (ALT == 'CGGCA')
```

Flags a specific 4-bp insertion (C > CGGCA) that has been empirically observed as a recurrent false positive in the Kestrel output. This artifact likely arises from k-mer graph ambiguity in GC-rich regions of the VNTR.

### Low_Depth_Conserved_Motifs

```
(Depth_Score < 0.4) and (Motif in ['1', '2', '3', '4', '6', '7', '8', '9'])
```

Flags variants occurring in conserved repeat unit motifs (numbered motifs 1--9) when the depth score is below 0.4. These motifs are highly conserved across MUC1 VNTR alleles, making true pathogenic variants in these positions unlikely unless strongly supported by sequencing depth.

## Duplicate Flagging

A separate mechanism identifies potential duplicate variant calls. When enabled, variants are grouped by `REF` and `ALT` alleles, sorted by depth score (descending), and all but the first (highest-scoring) entry in each group are flagged as `Potential_Duplicate`.

The duplicate flagging configuration:

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

`sort_by` is `Depth_Score` descending only, per [@hassansaei on #197](https://github.com/hassansaei/VNtyper/issues/197):
"Fall back to the 1.3 Depth_Score-only rule [...] Do not use `Motifs` or `Motif` as a sort key." An earlier
revision of this config sorted on three keys, including the plural `Motifs`, which does not exist as a column
by the time flagging runs and raised `KeyError` if the toggle were ever enabled.

!!! note
    Duplicate flagging is **disabled** by default in the current configuration. Per the same decision on #197,
    it stays disabled in the shipped config ("We have already tested with this setup. I do not know what will
    happen if we turn it on!"). Enable it by setting `"enabled": true` in `kestrel_config.json`.

## The `regex_match` Helper

Flagging rules can use a built-in `regex_match(pattern, value)` function for pattern-based matching. For example:

```json
{
  "Motif_X_Pattern": "regex_match('^X', Motif) and Depth_Score < 0.01"
}
```

This flags variants where the motif name starts with "X" and the depth score is below 0.01. The function uses Python's `re.search` internally.

## Adding Custom Rules

To add a new flagging rule:

1. Open `vntyper/scripts/kestrel_config.json`
2. Add a new key-value pair under `"flagging_rules"`:

```json
{
  "flagging_rules": {
    "My_Custom_Flag": "(Depth_Score < 0.005) and (Variant == 'Insertion')"
  }
}
```

The condition string has access to all columns in the variant DataFrame at the time of evaluation, including: `REF`, `ALT`, `POS`, `Motif`, `Variant`, `Depth_Score`, `Confidence`, `Estimated_Depth_AlternateVariant`, `Estimated_Depth_Variant_ActiveRegion`, and `is_valid_frameshift`.

## Impact on Variant Selection

During the final variant selection step, the selection priority is:

1. Highest confidence level
2. **Unflagged preferred over flagged**
3. Highest depth score
4. Highest haplo_count
5. Lowest genomic position

This means a High_Precision unflagged variant will always be selected over a High_Precision flagged variant, even if the flagged variant has a higher depth score. This behavior ensures that known artifact patterns do not take priority over cleaner calls.
