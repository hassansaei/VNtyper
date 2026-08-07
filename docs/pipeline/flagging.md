# Variant Flagging

Flagging applies configurable, post-hoc empirical filters to variant calls. Unlike the scoring and confidence system (which determines whether a variant is real), flagging annotates calls with named patterns observed in Kestrel output.

There are **two classes of flag**, and they have different consequences:

| Class | Declared in | Effect on the reported call |
|-------|-------------|-----------------------------|
| **Advisory** | `flagging_rules` only | The call is reported. The flag lowers its priority during variant selection, and appears in the `Flag` column of the report. |
| **Artifact** | `flagging_rules` **and** `artifact_flags` | The call is **excluded** from `kestrel_result.tsv`. It survives in `kestrel_pre_result.tsv` with `flag_filter_pass = False`. |

An advisory flag identifies a call that may be technically valid but warrants additional scrutiny. An artifact flag identifies a call that is not a candidate variant at all --- a known technical artifact of the k-mer genotyping --- so reporting it as a call would be wrong rather than merely cautious.

!!! warning "Artifact flags exclude, they do not merely deprioritise (Issue #174)"
    Before the Issue #174 fix, *every* flag was advisory. A sample whose only call was
    the 4-bp insertion artifact produced one `High_Precision` row with a non-empty `Flag`,
    which `report_config.json` maps to the `High_Precision_flagged` screening state ---
    and `is_finding` treats that as a **positive** MUC1 finding. The HTML report
    therefore styled a known technical artifact as a positive result.

    Which flags are artifacts is configuration, not code: see
    [Artifact Flags](#artifact-flags) below.

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

This is the one flag the shipped configuration declares as an **artifact** (see [Artifact Flags](#artifact-flags)), so a row carrying it is excluded from `kestrel_result.tsv`. A sample whose only call is this insertion is reported as **negative**.

### Low_Depth_Conserved_Motifs

```
(Depth_Score < 0.4) and (Motif in ['1', '2', '3', '4', '6', '7', '8', '9'])
```

Flags variants occurring in conserved repeat unit motifs (numbered motifs 1--9) when the depth score is below 0.4. These motifs are highly conserved across MUC1 VNTR alleles, making true pathogenic variants in these positions unlikely unless strongly supported by sequencing depth.

This flag is **advisory**. A low-depth call in a conserved motif is still a call: it is reported, and the flag only lowers its priority if another candidate competes with it. Excluding these rows would delete real low-depth calls.

## Artifact Flags

`artifact_flags` lists the flag names that disqualify a row from being reported:

```json
{
  "artifact_flags": ["False_Positive_4bp_Insertion"]
}
```

A flag name must appear in **both** `flagging_rules` (which is what raises it) and `artifact_flags` (which is what makes it disqualifying). A name in `artifact_flags` that no rule can raise is inert.

### How the exclusion works

After flagging, VNtyper derives a boolean column `flag_filter_pass`, which is `False` exactly when the row's `Flag` value contains one of the declared artifact names:

- The `Flag` value is comma-separated, so membership is tested **per element**, not as a substring. An artifact named `X` does not match a flag named `XY`.
- A row carrying both an artifact flag and an advisory flag is excluded --- one artifact is enough.
- Rows with no `Flag` value, or a missing `Flag` column (a negative run carries none), pass the gate. Absence of evidence is not evidence of an artifact.

`flag_filter_pass` is then the sixth of the boolean gates that the final filter requires to all be `True`, alongside `is_frameshift`, `is_valid_frameshift`, `depth_confidence_pass`, `alt_filter_pass` and `motif_filter_pass`.

!!! note "Nothing is destroyed"
    Consistent with the pipeline's "stages mark, they do not filter" contract, the
    artifact gate **marks**; only the final filter drops rows. Every excluded row is
    written to `kestrel_pre_result.tsv` with `flag_filter_pass = False` beside the
    `Flag` value that explains it, so the evidence for a sample is always recoverable.

### Reverting or narrowing the decision

Emptying the list restores the previous behaviour, where every flag was advisory, with no code change:

```json
{
  "artifact_flags": []
}
```

The flag name is never written into the Python; `add_artifact_gate` reads the list from `kestrel_config.json`. Narrowing or withdrawing the artifact rule is therefore a configuration edit, made by whoever owns the domain judgement.

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

A new rule is **advisory** by default. It becomes an artifact rule only if you also add its name to `artifact_flags`, and you should only do so for a pattern that is not a candidate variant at all.

## Impact on the Reported Call

The two flag classes act at different points in the pipeline, in this order.

### Artifact flags: exclusion, at the final filter

An artifact-flagged row has `flag_filter_pass = False` and is dropped by the final filter, **before** selection runs. It cannot be selected, cannot appear in `kestrel_result.tsv`, and cannot make a sample positive. If it was the sample's only candidate, the sample is reported as negative.

### Advisory flags: priority, during selection

Rows that survive the filter are ranked, and the selection priority is:

1. Highest confidence level
2. **Unflagged preferred over flagged**
3. Highest depth score
4. Highest haplo_count
5. Lowest genomic position

This means a High_Precision unflagged variant will always be selected over a High_Precision flagged variant, even if the flagged variant has a higher depth score. This behavior ensures that flagged calls do not take priority over cleaner ones.

Step 2 still applies to advisory flags only, since artifact-flagged rows are already gone by the time selection runs.
