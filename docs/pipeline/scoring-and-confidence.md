# Scoring and Confidence Assignment

VNtyper uses two complementary scoring systems to evaluate variant calls: the **frame score** (biological validity) and the **depth score** (signal strength). Together, these determine the confidence label assigned to each variant.

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

| Level | Criteria | Clinical Interpretation |
|-------|----------|------------------------|
| **High_Precision*** | Alt depth >= 100 and Depth Score >= 0.00515 | Very high confidence call; strong supporting evidence |
| **High_Precision** | Alt depth >= 21 and < 100, Depth Score >= 0.00515, Region depth > 200 | High confidence call suitable for clinical consideration |
| **Low_Precision** | Alt depth <= 20, or Depth Score between 0.00469--0.00515, or Region depth <= 200 | Variant detected but with marginal evidence; requires independent validation |
| **Negative** | Depth Score < 0.00469 | Signal below noise threshold; variant is likely an artifact |

!!! warning "Empirically derived thresholds"
    These thresholds were calibrated on a cohort of known-positive and known-negative samples as described in Saei et al. (2023). They are specific to the MUC1 VNTR assay and should not be applied to other genomic regions or variant types without re-calibration.

### Assignment Logic

The confidence assignment follows a layered rule system where conditions are applied in sequence and later conditions can overwrite earlier assignments. All variants start as "Negative":

1. Variants with Depth Score **below** the low threshold (0.00469) remain Negative regardless of other metrics
2. Variants at exactly the low threshold or with region depth <= 200 are assigned Low_Precision
3. Variants with alt depth >= 100 **and** Depth Score >= 0.00515 are upgraded to High_Precision*
4. Variants with alt depth 21--100 **and** Depth Score >= 0.00515 are assigned High_Precision
5. Remaining variants with Depth Score between the low and high thresholds receive Low_Precision

A boolean column `depth_confidence_pass` is set to `True` for all non-Negative variants, enabling downstream filtering.

## Clinical Interpretation

The confidence level directly informs clinical decision-making:

- **High_Precision / High_Precision***: The variant call is supported by sufficient evidence for clinical reporting, though orthogonal validation (e.g., SNaPshot for dupC, long-read sequencing for complex variants) is recommended per clinical guidelines.
- **Low_Precision**: The variant signal is present but marginal. Independent validation is essential before clinical action.
- **Negative**: No evidence of a pathogenic VNTR variant above the noise floor.

## Reference

Saei H. et al., *iScience* 26, 107171 (2023).
