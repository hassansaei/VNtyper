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

Each flagging rule is a named, structured conjunction in `kestrel_config.json`. VNtyper validates the complete rule set against the DataFrame's available columns before processing any row. If every predicate in a rule evaluates to `true`, the corresponding flag name is appended to the variant's `Flag` column. Multiple flags are comma-separated in configuration order. Variants matching no rules receive `Flag = "Not flagged"`.

!!! info "Flagging occurs before variant selection"
    As of VNtyper 2 (Issue #145 fix), flagging is applied **before** the final variant selection step. This ensures that when multiple candidate variants pass all filters, unflagged variants are preferred over flagged ones. Previously, a flagged variant could be selected as the best call because flags were added after selection.

## Current Flagging Rules

The default rules in `kestrel_config.json`:

### False_Positive_4bp_Insertion

```json
{
  "all": [
    {"left": {"column": "REF"}, "operator": "eq", "right": {"literal": "C"}},
    {"left": {"column": "ALT"}, "operator": "eq", "right": {"literal": "CGGCA"}}
  ]
}
```

Flags a specific 4-bp insertion (C > CGGCA) that has been empirically observed as a recurrent false positive in the Kestrel output. This artifact likely arises from k-mer graph ambiguity in GC-rich regions of the VNTR.

This is the one flag the shipped configuration declares as an **artifact** (see [Artifact Flags](#artifact-flags)), so a row carrying it is excluded from `kestrel_result.tsv`. A sample whose only call is this insertion is reported as **negative**.

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

Flags variants occurring in conserved repeat unit motifs (numbered motifs 1--9) when the depth score is below 0.4. These motifs are highly conserved across MUC1 VNTR alleles, making true pathogenic variants in these positions unlikely unless strongly supported by sequencing depth.

This flag is **advisory**: by itself it does not exclude a row, and it only lowers the row's priority if another candidate competes with it. The separate position-based motif gate still excludes right-half calls in the motifs configured by `motif_filtering.exclude_motifs_right`; see [Kestrel Step 6](kestrel.md#step-6-motif-correction-and-annotation).

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

The artifact decision is never written into Python; `add_artifact_gate` reads the list from `kestrel_config.json`. Narrowing or withdrawing the artifact rule is therefore a configuration edit, made by whoever owns the domain judgement. A separate, flag-name-scoped compatibility map contains the byte-exact previous rule solely to migrate last-release configuration.

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

## Rule Schema and Adding Custom Rules

To add a new flagging rule:

1. Open `vntyper/scripts/kestrel_config.json`
2. Add a new key-value pair under `"flagging_rules"`:

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

Every rule has exactly one key, `all`, whose value is a non-empty list. Every predicate has exactly `left`, `operator`, and `right`; each operand contains exactly one `column` or `literal`. A column must be in the explicit allowlist formed from the DataFrame columns available when flagging begins. Common Kestrel columns include `REF`, `ALT`, `POS`, `Motif`, `Variant`, `Depth_Score`, `Confidence`, `Estimated_Depth_AlternateVariant`, `Estimated_Depth_Variant_ActiveRegion`, and `is_valid_frameshift`.

The only operators are:

| Operator | Meaning |
|----------|---------|
| `eq` | Strict same-family equality. Booleans compare only with booleans. |
| `lt` | Numeric less-than; booleans are not numbers. |
| `in` | Left scalar membership in a non-empty homogeneous right literal list. |
| `casefold_eq` | String equality after Unicode case-folding. |

`None`, `pd.NA`, and floating NaN make a predicate false. Values are never coerced: for example, `"7"` is not the number `7`. Missing columns, malformed rules, and incompatible non-null row values abort flagging instead of silently disabling a rule. Flag names must be non-empty strings, cannot contain commas, and cannot use the reserved result values `Not flagged` or `Not applicable`.

Rules are JSON data, never executable source. Calls, attributes, indexing, imports, comprehensions, lambdas, regular expressions, arithmetic, `or`, `not`, and nested boolean forms are unsupported. Code-shaped text inside a `literal` remains inert data; a rule supplied as such a string is rejected.

### Migrating a rule from the immediately preceding release

The last-release Kestrel string

```json
"(Depth_Score < 0.4) and (Motif in ['1', '2', '3', '4', '6', '7', '8', '9'])"
```

migrates to the structured `Low_Depth_Conserved_Motifs` object shown above. For upgrade compatibility, VNtyper accepts only each flag name's byte-exact string from the release immediately before Issue #286. Whitespace edits, renamed columns, added clauses, custom strings, and another flag's historical expression are rejected. New and edited configurations must use the structured form.

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

## adVNTR Flagging Rules

Everything above describes Kestrel's flags, declared in `kestrel_config.json`. The optional
adVNTR module has its own three rules in `vntyper/modules/advntr/advntr_config.json`,
evaluated by the same `add_flags` machinery against adVNTR's result frame.

None of them is an artifact flag in the `artifact_flags` sense: adVNTR has no equivalent
gate, so a flagged adVNTR call is still reported. `report_config.json` maps it to
`positive flagged`, which `is_finding` still counts as a finding — the flag downgrades a
call visibly, it does not withdraw it.

| Rule | Condition | What it asserts |
|------|-----------|-----------------|
| `Low_Coverage` | `NumberOfSupportingReads < 10` | Too few supporting reads to be confident. A threshold, not a claim about the state. |
| `Repeat_Unit_7` | `RU == '7'` | A call in repeat unit 7, an established recurrent artifact. |
| `Polymorphic_Call` | `Variant in [...24 states...]` | The state string matches one of a list of recurrent events treated as artifacts. |

### What `Polymorphic_Call` asserts, and where the list came from

@hassansaei on [#267](https://github.com/hassansaei/VNtyper/issues/267): a polymorphic or
artifact call is *"a recurrent event observed in positive and/or negative samples, likely
originating from motif differences at that locus rather than a true pathogenic event"*.
The entries were derived by running an older adVNTR — pre-2.0.4 — over the `renome`
cohort and recording the states that recurred across many samples. **No per-entry
observation count was retained.**

That provenance now lives in `vntyper/modules/advntr/advntr_calibration.json`, one record
per entry with a `status` of `confirmed_artifact` or `pending_renome_revalidation`.
Production code never reads that file; `tests/unit/test_advntr_polymorphic_calls.py`
asserts it and the live rule agree. It is the same pattern
`vntyper/scripts/calibration.json` uses for the Kestrel constants.

!!! warning "The key cannot separate an artifact from a pathogenic variant"
    adVNTR's State string records the *shape* of an event — the repeat unit, the length,
    and the **first inserted base only** (`advntr/vntr_finder.py`: *"If there are run of
    insertions, the sequence might differ, but we just take the first base"*). It does not
    record the allele.

    So a recurrent benign single-base insertion and a pathogenic single-base duplication in
    the same unit produce the **same string**, and no denylist keyed on that string can
    tell them apart, whatever cohort it was derived from. The MUC1 array carries repeat
    units with both `TTGGGGGG` and `TGGGGGGG`, so a G inserted into that run is placeable
    either way — @hassansaei's note on #267, and the same ambiguity seen from the sequence
    side.

    Measured on the 400-run simulated benchmark (VNtyper 2.0.22, GRCh38, refs-v2): of the
    172 carriers adVNTR detects, **8 carry a `Polymorphic_Call` row** — every one of the
    5 `dupA` carriers it detects, and 3 of 7 `delGCCCA` carriers. No control sample
    produced any adVNTR call at all.

    Keying the flag on something richer than the state string is
    [#267](https://github.com/hassansaei/VNtyper/issues/267)'s suggestion 4 and is not
    decided.

### Reachability: why the list is 24 entries and not 32

Flagging runs **after** the pathogenic-frame filter (`advntr_processing_ins` /
`advntr_processing_del`), which keeps only rows whose signed net indel change Δ satisfies
`Δ % 3 == 1`. An entry that does not is removed before `add_flags` is ever reached and can
never fire.

Seven of the 32 shipped entries were in that state, and one was listed twice. They were
removed in the #267 cleanup with the owner's agreement; replaying both lists over the
adVNTR output of all 400 simulated samples changes **no** `Flag` value.

`tests/unit/test_advntr_polymorphic_calls.py` runs the production filter arms over the
live list, so a future unreachable entry fails the build. That matters here because the
the historical failure mode was silent: a missing expression name could disable a rule,
which is how `Polymorphic_Call` shipped misspelled as `Poylmorhic_Call` until `742b872`,
and loose expression typing allowed `Repeat_Unit_7` to ship as `RU == 7` — comparing a
string column against an integer — until `52f822e`. Structured validation now rejects
both defects before processing rows.

### What is still open

23 of the 24 live entries are recorded as `pending_renome_revalidation`. @hassansaei asked
on #267 for them to be re-measured against the re-analysed renome cohort and decided case
by case; only `D58_2&D59_2` (and the separate `Repeat_Unit_7` rule) are confirmed. Until
then they remain flagged exactly as shipped. Acting on that decision is a data edit to
`advntr_config.json` and `advntr_calibration.json` — no code change.
