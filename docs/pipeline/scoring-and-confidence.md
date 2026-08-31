# Scoring and Confidence Assignment

VNtyper 2 uses two complementary scoring systems to evaluate variant calls: the **frame score** (biological validity) and the **depth score** (signal strength). Together, these determine the confidence label assigned to each variant.

## Frame Score

The frame score quantifies whether a variant introduces a frameshift in the MUC1 coding sequence:

$$
\text{Frame Score} = \frac{\text{len(ALT)} - \text{len(REF)}}{3}
$$

A variant is classified as a **frameshift** when `(len(ALT) - len(REF)) % 3 != 0`. Only frameshift variants are pathologically relevant for ADTKD-MUC1, because in-frame insertions or deletions preserve the reading frame and do not produce the truncated mucin-1 protein.

!!! info "Frameshift patterns in MUC1"
    The pipeline further classifies frameshifts by their modular arithmetic pattern:

    - **Insertion frameshifts**: direction > 0 and `frameshift_amount == 1` (i.e., 3n+1 bp inserted)
    - **Deletion frameshifts**: direction < 0 and `frameshift_amount == 2` (i.e., 3n+2 bp deleted)

    These specific patterns correspond to the known pathogenic mutation signatures in the MUC1 VNTR. Variants matching these patterns are marked as `is_valid_frameshift = True`.

!!! note "Why insertions and deletions use different remainders ([#181](https://github.com/hassansaei/VNtyper/issues/181))"
    An insertion of 3n+1 bases and a deletion of 3n+2 bases are frame-equivalent: both
    shift the reading frame by Delta = +1 (mod 3). That is the pathogenic ADTKD-MUC1
    reading frame, the one that produces the toxic MUC1-fs neo-protein (classically
    exemplified by dupC). The opposite pair -- an insertion of 3n+2 bases or a deletion
    of 3n+1 bases -- shifts into the other frame, which has not been established as
    pathogenic in patients and is treated as unknown / not clinically identified for
    ADTKD-MUC1. Rejecting a (3n+1)-bp deletion is therefore a considered, MUC1-specific
    choice, not accidental asymmetry -- and not a lost call.

### Haplo Count

After frame scoring, a **haplo_count** is computed for each variant: the number of times the exact same variant (POS, REF, ALT) appears across different haplotype calls from Kestrel. A higher haplo_count indicates more supporting evidence and is used as a tie-breaker during variant selection.

## Depth Score

The depth score represents the signal-to-noise ratio for a variant call:

$$
\text{Depth Score} = \frac{\text{Estimated Depth (Alternate Variant)}}{\text{Estimated Depth (Variant Active Region)}}
$$

where:

- **Estimated Depth (Alternate Variant)** is the number of k-mer paths supporting the alternate allele
- **Estimated Depth (Variant Active Region)** is the total k-mer depth across the active region where the variant was called

Infinite values (from division by zero) are replaced with NaN. Higher depth scores indicate stronger variant signals relative to background.

## Confidence Assignment

Confidence labels are assigned based on a combination of depth score, alternate allele depth, and active region depth. All thresholds are empirically derived from Saei et al., *iScience* 26, 107171 (2023).

### Threshold Configuration

The thresholds are defined in `kestrel_config.json`:

| Parameter | Config Key | Value |
|-----------|-----------|-------|
| Depth score (low) | `depth_score_thresholds.low` | 0.00469 |
| Depth score (high) | `depth_score_thresholds.high` | 0.00515 |
| Alt depth (low) | `alt_depth_thresholds.low` | 20 |
| Alt depth (mid_low) | `alt_depth_thresholds.mid_low` | 21 |
| Alt depth (mid_high) | `alt_depth_thresholds.mid_high` | 100 |
| Region depth threshold | `var_active_region_threshold` | 200 |

### Confidence Levels

| Level | Criteria | What the evidence supports |
|-------|----------|------------------------|
| **High_Precision*** | Alt depth >= 100 **and** Depth Score **> 0.00515** | Very high confidence call; strong supporting evidence |
| **High_Precision** | Alt depth >= 21 and < 100, **and** Depth Score **> 0.00515** | High confidence call: both calibrated thresholds are met |
| **Low_Precision** | Depth Score anywhere in the **closed** band `0.00469 <= Depth Score <= 0.00515`, at **any** alt depth; **or** alt depth <= 20 with Depth Score > 0.00515 | Variant detected but with marginal evidence; requires independent validation |
| **Negative** | Depth Score < 0.00469, or undefined (active-region depth of 0 -> NaN) | Signal below noise threshold; variant is likely an artifact |

The two High tiers use a **strict** `>` against the high threshold. The rules that
assign them are written with `>=` (see the ordered list below), but the mid-band
demotion runs **last** and its band is closed at the top, so a Depth Score of exactly
0.00515 ends up Low_Precision. This is the behaviour @hassansaei specified on
[#184](https://github.com/hassansaei/VNtyper/issues/184); see
"[What #184 changed](#what-184-changed)".

### Decision table

Depth Score is `alt depth / region depth`, so the two depths are **not** independent
axes — pick the Depth Score band first, then read across to the alternate depth.

| Depth Score | Alt depth <= 20 | Alt depth 21--99 | Alt depth >= 100 |
|-------------|-----------------|------------------|------------------|
| undefined (region depth 0) | Negative | Negative | Negative |
| < 0.00469 | Negative | Negative | Negative |
| **= 0.00469** exactly | Low_Precision | Low_Precision | Low_Precision |
| 0.00469 -- 0.00515 (interior) | Low_Precision | Low_Precision | Low_Precision |
| **= 0.00515** exactly | Low_Precision | Low_Precision | Low_Precision |
| > 0.00515 | Low_Precision | High_Precision | High_Precision* |

The region-depth threshold does not appear in this table because it never changes the
outcome for any of these cells — see the warning below for why.

!!! note "The alt-depth thresholds leave a gap between 20 and 21"
    `alt_depth_thresholds.low` is 20 and `mid_low` is 21, so no rule covers an alt
    depth strictly between them. A variant with, say, an alt depth of 20.5 and a Depth
    Score above 0.00515 matches no condition at all unless its region depth is at most
    200, so it is reported **Negative** on a deep active region and Low_Precision on a
    shallow one — the only place in the whole table where the region-depth rule decides
    anything. Alternate depths come from Kestrel's `Sample` field as k-mer-path depths
    and are always whole numbers in practice, so this gap is unreachable in production —
    but `confidence_assignment.py` performs no integer cast, so nothing in the module
    itself enforces it.

!!! warning "The region-depth threshold does not cap the confidence label"
    A region depth at or below 200 demotes a variant to Low_Precision, but that
    demotion is applied **first** and is overwritten by either High_Precision tier
    whenever one applies. A variant with an alt depth of 50 and a total active-region
    k-mer depth of 150 is reported as High_Precision today, not Low_Precision — the same
    label it would receive at a total active-region k-mer depth of 5000.

    It is in fact overwritten in **every** case where a High tier is in play, as
    arithmetic rather than as a coincidence: a non-zero region depth of at most 200
    with an alt depth of at least 21 forces `Depth Score >= 21/200 = 0.105`, twenty
    times the high threshold, so rule 2 or rule 5 always fires. Where the region
    demotion is *not* overwritten, the label would have been Low_Precision anyway
    (alt depth <= 20, or a Depth Score inside the mid-band). A region depth of 0 makes
    the Depth Score undefined, which keeps the variant Negative regardless.

    Tier precedence is **specified**: @hassansaei decided on
    [#183](https://github.com/hassansaei/VNtyper/issues/183) (2026-08-06) that the
    2.x last-wins sequential assignment shown above is the intended behaviour, and
    that the 1.3 absolute region-depth <=200 cap must not be restored:

    > "Keep the current 2.x last-wins logic — do not restore the absolute
    > region-depth <=200 cap from 1.3. In practice it is very unlikely that a
    > variant with region depth <200 is later promoted to High_Precision by a
    > subsequent rule. Where that pattern can appear, it is mostly for early
    > (beginning) and late conserved motifs; we already have a flagging rule
    > when Depth_Score is far from ~0.5 (50%) [...] the intentional behaviour
    > going forward is the 2.x sequential assignment as implemented today."

    The behaviour is pinned by `tests/unit/test_confidence_boundaries.py`; see
    issue #179.

!!! warning "Empirically derived thresholds"
    These thresholds were calibrated on a cohort of known-positive and known-negative samples as described in Saei et al. (2023). They are specific to the MUC1 VNTR assay and should not be applied to other genomic regions or variant types without re-calibration.

### Assignment Logic

The confidence assignment follows a layered rule system where conditions are applied in sequence and **later conditions overwrite earlier assignments**. All variants start as "Negative", and a variant whose Depth Score is below the low threshold (0.00469) — or whose Depth Score is undefined, because the active-region depth is 0 — stays Negative no matter which conditions match. The six conditions are applied in this order, and each is transcribed here exactly as `confidence_assignment.py` writes it:

1. `cond1` — Depth Score **exactly** 0.00469, **or** region depth <= 200 -> Low_Precision
2. `cond2` — Alt depth >= 100 **and** Depth Score >= 0.00515 -> High_Precision*
3. `cond3` — Alt depth 21--99 **and** Depth Score between 0.00469 and 0.00515 inclusive -> Low_Precision
4. `cond4` — Alt depth <= 20 -> Low_Precision
5. `cond5` — Alt depth 21--99 **and** Depth Score >= 0.00515 -> High_Precision
6. `cond_midband` — Depth Score between 0.00469 and 0.00515 **inclusive at both ends** -> Low_Precision

Three consequences of the ordering, none of them a property of any single rule:

- **Rules 2 and 5 are written with `>=` but behave like `>`.** Rule 6 runs last and its
  band is closed at the top, so it takes back every row that rules 2 and 5 promoted at
  exactly 0.00515. That is why the criteria table above states the High tiers with a
  strict `>`.
- **Rule 1's region-depth demotion is overwritten by rules 2 and 5** whenever either
  applies.
- **Rule 3 never decides a row.** Rule 6 covers the same closed Depth Score band at
  *every* alternate depth, writes the same label, and runs after it, so rule 3 is
  wholly subsumed. It is retained deliberately because it names the intent — see the
  comment on it in `confidence_assignment.py` and
  [#184](https://github.com/hassansaei/VNtyper/issues/184).

A boolean column `depth_confidence_pass` is set to `True` for all non-Negative variants, enabling downstream filtering.

### What #184 changed

Before [#184](https://github.com/hassansaei/VNtyper/issues/184), rule 6 covered the
**open** interval `0.00469 < Depth Score < 0.00515`. It now covers the **closed**
interval. The interior of the band and the lower endpoint were already Low_Precision
(the lower endpoint via rule 1), so the only Depth Score whose label moved is
**exactly 0.00515**, and at that point the change is a demotion:

| Alt depth at Depth Score = 0.00515 | Before #184 | After #184 |
|------------------------------------|-------------|------------|
| <= 20 | Low_Precision | Low_Precision (unchanged) |
| 21--99 | High_Precision | **Low_Precision** |
| >= 100 | High_Precision* | **Low_Precision** |

With whole-number depths, `Depth Score == 0.00515` requires `alt = 103k` and
`region = 20000k` (0.00515 is exactly 103/20000, and that division is exact in IEEE 754),
so the smallest reachable alternate depth at the boundary is 103. That is above the
`mid_high` threshold of 100, which means **the only tier that moves on real data is
`High_Precision*` -> `Low_Precision`**; the 21--99 row above is reachable only with a
fractional depth.

Because `select_single_best_variant` ranks candidates on confidence before anything
else, a demotion at this boundary can change which variant is reported for a sample
that carries more than one. `tests/unit/test_confidence_boundaries.py` pins both the
label and that selection consequence.

## Interpreting the confidence levels

The confidence level describes how much evidence supports the call:

- **High_Precision / High_Precision***: Both calibrated thresholds are met. Orthogonal validation (e.g., SNaPshot for dupC, long-read sequencing for complex variants) is recommended before the call is relied on.
- **Low_Precision**: The depth score falls inside the closed calibration band, or the alternate depth is at or below 20. Independent validation is essential before the call is relied on.
- **Negative**: No evidence of a pathogenic VNTR variant above the noise floor.

## Reference

Saei H. et al., *iScience* 26, 107171 (2023).
