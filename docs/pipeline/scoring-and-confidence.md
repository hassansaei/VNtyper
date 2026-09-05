# Scoring and Confidence Assignment

VNtyper evaluates variant calls using two complementary scoring mechanisms: the **frame score** (biological validity) and the **depth score** (signal-to-noise ratio). Together, these determine the confidence tier assigned to each call.

## Frame Score

The frame score evaluates whether a variant shifts the MUC1 coding reading frame:

$$
\text{Frame Score} = \frac{\text{len(ALT)} - \text{len(REF)}}{3}
$$

A variant is classified as a **frameshift** when `(len(ALT) - len(REF)) % 3 != 0`. Only frameshift variants cause ADTKD-MUC1; in-frame indels preserve the reading frame and do not produce the pathogenic truncated neo-protein.

!!! info "Frameshift patterns in MUC1"
    The pipeline classifies frameshifts by modular direction:

    - **Insertion frameshifts**: direction > 0 and `frameshift_amount == 1` (3n+1 bp insertion)
    - **Deletion frameshifts**: direction < 0 and `frameshift_amount == 2` (3n+2 bp deletion)

    These patterns match the established pathogenic mutation signatures in the MUC1 VNTR. Matching variants receive `is_valid_frameshift = True`.

!!! note "Why insertions and deletions use different remainders (Issue #181)"
    An insertion of 3n+1 bases and a deletion of 3n+2 bases are frame-equivalent: both shift the reading frame by +1 (mod 3). This specific shift creates the toxic MUC1-fs neo-protein (exemplified by canonical dupC). The opposite pair (an insertion of 3n+2 bases or a deletion of 3n+1 bases) shifts into the alternate frame, which has not been identified as pathogenic in patients. Rejecting a (3n+1)-bp deletion is a deliberate biological filter, not accidental asymmetry.

### Haplo Count

After frame scoring, VNtyper computes `haplo_count`: the frequency with which the identical variant call (POS, REF, ALT) occurs across Kestrel's reconstructed haplotype assemblies. Higher counts indicate stronger assembly support and serve as a tie-breaker during variant selection.

## Depth Score

The depth score measures the signal strength of a variant relative to its locus:

$$
\text{Depth Score} = \frac{\text{Estimated Depth (Alternate Variant)}}{\text{Estimated Depth (Variant Active Region)}}
$$

Definitions:

- **Estimated Depth (Alternate Variant)**: The number of k-mer paths supporting the alternate allele.
- **Estimated Depth (Variant Active Region)**: Total k-mer depth across the active region containing the call.

Division-by-zero occurrences evaluate to NaN. Higher depth scores indicate higher relative variant abundance.

## Confidence Assignment

Confidence tiers are assigned based on depth score, alternate allele depth, and active region depth. Thresholds are derived from Saei et al., *iScience* 26, 107171 (2023).

### Threshold Configuration

Defined in `kestrel_config.json`:

| Parameter | Config Key | Value |
|-----------|-----------|-------|
| Reporting floor | `reporting_floor` | 0.00469 |
| Depth score (low) | `depth_score_thresholds.low` | 0.00469 |
| Depth score (high) | `depth_score_thresholds.high` | 0.00515 |
| Alt depth (low) | `alt_depth_thresholds.low` | 20 |
| Alt depth (mid_low) | `alt_depth_thresholds.mid_low` | 21 |
| Alt depth (mid_high) | `alt_depth_thresholds.mid_high` | 100 |
| Region depth threshold | `var_active_region_threshold` | 200 |

Three calibrated parameters share the base value `0.00469`:

1. `confidence_assignment.reporting_floor` (0.00469): Minimum noise floor; variants below this are labeled `Negative`.
2. `confidence_assignment.depth_score_thresholds.low` (0.00469): Lower boundary of the `Low_Precision` band (`[low, high]`).
3. `alt_filtering.gg_depth_score_threshold` (0.00469): Minimum depth score required for `ALT == "GG"` variants.

### Confidence Levels

| Level | Criteria | Evidence assessment |
|-------|----------|---------------------|
| **High_Precision*** | Alt depth >= 100 **and** Depth Score **> 0.00515** | Robust variant call with extensive k-mer support. |
| **High_Precision** | Alt depth >= 21 and < 100, **and** Depth Score **> 0.00515** | High-confidence call meeting both calibrated cutoffs. |
| **Low_Precision** | Depth Score in closed band `0.00469 <= Depth Score <= 0.00515`, at **any** alt depth; **or** alt depth <= 20 with Depth Score > 0.00515 | Marginal signal; requires orthogonal validation. |
| **Negative** | Depth Score < 0.00469, or undefined (active region depth 0 -> NaN) | Subthreshold noise; candidate is discarded. |

The High tiers require a strict inequality (`>`) against 0.00515. The mid-band demotion runs last and defines a closed interval at its upper limit, ensuring that a Depth Score of exactly 0.00515 receives `Low_Precision` (Issue #184).

### Decision Table

Depth Score is `alt depth / region depth`, making the metrics interdependent:

| Depth Score | Alt depth <= 20 | Alt depth 21 to 99 | Alt depth >= 100 |
|-------------|-----------------|--------------------|------------------|
| undefined (region depth 0) | Negative | Negative | Negative |
| < 0.00469 | Negative | Negative | Negative |
| **= 0.00469** exactly | Low_Precision | Low_Precision | Low_Precision |
| 0.00469 to 0.00515 (interior) | Low_Precision | Low_Precision | Low_Precision |
| **= 0.00515** exactly | Low_Precision | Low_Precision | Low_Precision |
| > 0.00515 | Low_Precision | High_Precision | High_Precision* |

The active-region depth threshold does not alter outcomes across this table.

!!! note "Alt depth interval between 20 and 21"
    Because `alt_depth_thresholds.low` is 20 and `mid_low` is 21, fractional values between them lack an explicit integer rule. Alternate depths come from Kestrel's `Sample` field as k-mer-path depths and are always whole numbers in practice, so this gap is unreachable in production. However, `confidence_assignment.py` performs no integer cast, so nothing in the module itself enforces it.

!!! warning "Active region depth threshold does not cap confidence labels"
    An active-region depth <= 200 demotes candidates to Low_Precision, but this demotion runs first and is superseded by High_Precision rules whenever applicable. A variant with an alt depth of 50 and a total active-region k-mer depth of 150 is reported as High_Precision today, not Low_Precision, exactly the same label it would receive at a total active-region k-mer depth of 5000.

    When alt depth is >= 21 and region depth is <= 200, the resulting Depth Score is at least 21/200 = 0.105 (twenty times the high threshold), causing rule 2 or rule 3 to fire. Where the demotion is not overwritten, the variant would be Low_Precision or Negative regardless.

    This precedence was confirmed on Issue #183: the sequential last-wins assignment is intentional, and the legacy v1.3 absolute cap (`region depth <= 200`) is deprecated.

### Assignment Logic

Confidence evaluation follows an ordered first-match table in `vntyper/scripts/confidence_rules.py`:

| Order | Condition | Confidence Label | Rationale |
|-------|-----------|------------------|-----------|
| 0 | `Depth_Score` is NaN or `< 0.00469` | `Negative` | Outer reporting floor: eliminates subthreshold candidates. |
| 1 | `0.00469 <= Depth_Score <= 0.00515` | `Low_Precision` | Mid-band demotion (Issue #184): covers closed calibration band across all alternate depths. |
| 2 | `Depth_Score > 0.00515` and alt depth `>= 100` | `High_Precision*` | High depth meeting high-confidence threshold. |
| 3 | `Depth_Score > 0.00515` and `21 <=` alt depth `< 100` | `High_Precision` | Standard depth meeting high-confidence threshold. |
| 4 | `Depth_Score > 0.00515` and alt depth `<= 20` | `Low_Precision` | Shallow alternate depth demotion. |
| 5 | `Depth_Score > 0.00515` and `20 <` alt depth `< 21` and region depth `<= 200` | `Low_Precision` | Fractional depth edge case with shallow active region. |
| 6 | Otherwise | `Negative` | Fallback for fractional depth gap with deep active region. |

Non-Negative variants set `depth_confidence_pass = True`, permitting downstream selection.

### Changes in Issue #184

Before Issue #184, the mid-band rule evaluated the open interval `0.00469 < Depth Score < 0.00515`. It now covers the **closed** interval `[0.00469, 0.00515]`. The only score affected is **exactly 0.00515**, which is demoted to `Low_Precision`:

| Alt depth at Depth Score = 0.00515 | Before #184 | After #184 |
|------------------------------------|-------------|------------|
| <= 20 | Low_Precision | Low_Precision (unchanged) |
| 21 to 99 | High_Precision | **Low_Precision** |
| >= 100 | High_Precision* | **Low_Precision** |

Under integer depths, `Depth Score == 0.00515` requires `alt = 103k` and `region = 20000k` (103/20000 is exact in IEEE 754 float). The minimum reachable alternate depth is 103, meaning the only tier shift on real data is `High_Precision*` to `Low_Precision`.

## Sample-Level Confidence Grade

VNtyper derives an overall sample-level confidence grade during screening summary generation. The grade synthesizes caller execution state, calling status, subthreshold signals, coverage metrics, and cross-caller concordance.

### Vocabulary and Pre-Emption Order

Evaluated against `confidence_grade_rules` in `report_config.json`, defaulting to `"not-established"`:

| Grade | Evaluated Conditions | Definition |
|-------|----------------------|------------|
| `not-established` | Execution failure (`kestrel_execution` or `advntr_execution == "failed"`), missing required caller, or invalid output token | Genotyping incomplete; outcome undetermined. |
| `no-finding-limited` | Negative call with subthreshold candidate, or negative call with coverage QC `FAIL` or `NOT_EVALUATED` | Negative call limited by shallow coverage or an unconfirmed subthreshold candidate. |
| `no-finding` | Negative call by both callers (or single caller when adVNTR off) with coverage QC `PASS` | Negative result with acceptable coverage depth across the array. |
| `finding-limited` | Variant call with coverage QC `FAIL` or `NOT_EVALUATED`, or a flagged Kestrel call | Variant detected, but confidence is reduced by coverage quality or caller flags. |
| `finding-corroborated` | Unflagged Kestrel call with coverage QC `PASS` and positive cross-caller concordance (`cross_match_is_positive == True`) | Variant detected with acceptable coverage and concordantly verified across callers. |
| `finding` | Unflagged call with coverage QC `PASS` lacking positive cross-match (or clean adVNTR finding) | Variant detected with passing coverage metrics. |

All outputs are for research use only.

## Reference

Saei H. et al., *iScience* 26, 107171 (2023).
