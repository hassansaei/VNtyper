"""
confidence_assignment.py

Module Purpose:
---------------
Contains logic for assigning confidence labels (e.g., Low_Precision, High_Precision)
to variants based on numeric scores such as Depth_Score. The thresholds and
confidence categories are derived from empirical cutoffs in Saei et al.,
iScience 26, 107171 (2023).

Typical Flow:
-------------
1. Coerce the two depth columns to numeric with ``pd.to_numeric(errors="coerce")``
   and fill the unparseable ones with 0. This is a float coercion, not an integer
   cast: a fractional depth is accepted and stays fractional. Production depths are
   whole numbers because they are read counts split out of Kestrel's ``Sample``
   field, but nothing in this module enforces that.
2. Compute Depth_Score = alt / active_region depth (float division; a zero region
   depth yields inf, which is then replaced with NaN).
3. Compare Depth_Score and alt-depth coverage to threshold intervals.
4. Assign textual confidence labels (e.g., 'Low_Precision', 'High_Precision').
5. Return updated DataFrame with 'Depth_Score' and 'Confidence' columns.

Rule ordering (see ``docs/pipeline/scoring-and-confidence.md`` for the full table):
the six conditions below are applied in sequence and later ones overwrite earlier
ones, so ``cond_midband`` -- the CLOSED interval [low, high], applied last -- has the
final word. A Depth_Score of exactly ``high`` is therefore Low_Precision even though
``cond2`` and ``cond5`` are written with ``>=``.

References:
-----------
- Saei et al., iScience 26, 107171 (2023)
- #6: Docs: Kestrel postprocessing heuristic (depth score usage)
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Confidence label for variants that fail depth-based filtering.
# Used as the default and for Depth_Score below the configured low_threshold.
NEGATIVE_LABEL = "Negative"


def calculate_depth_score_and_assign_confidence(df: pd.DataFrame, kestrel_config: dict) -> pd.DataFrame:
    """
    Calculates Depth_Score for each variant and assigns a confidence label
    based on coverage thresholds defined in kestrel_config['confidence_assignment'].

    (Refactored)
      - All rows remain in the DataFrame.
      - Variants with Depth_Score < low_threshold are assigned 'Negative' (filtered out).
      - A new boolean column 'depth_confidence_pass' is True if the row's
        final Confidence is not 'Negative'.
      - Depth_Score is computed as:
             Depth_Score = Estimated_Depth_AlternateVariant / Estimated_Depth_Variant_ActiveRegion
        Any infinite values (resulting from division by zero) are replaced with NaN.

    Args:
        df (pd.DataFrame):
            Must have numeric columns:
              - 'Estimated_Depth_AlternateVariant'
              - 'Estimated_Depth_Variant_ActiveRegion'
        kestrel_config (dict):
            Must contain a subdict 'confidence_assignment'. These calibration keys
            are required, not defaulted (shipped values in brackets):
              - 'depth_score_thresholds' -> 'low' [0.00469], 'high' [0.00515]
              - 'alt_depth_thresholds' -> 'low' [20], 'mid_low' [21],
                'mid_high' [100]
              - 'var_active_region_threshold' [200]
            'confidence_levels' is optional; its 'low_precision',
            'high_precision' and 'high_precision_star' labels fall back to their
            shipped strings because they are cosmetic, not calibration.

    Returns:
        pd.DataFrame:
            Same shape as input, with:
              - 'Depth_Score' (float)
              - 'Confidence' (str, e.g., 'Low_Precision', 'High_Precision', etc.)
              - 'depth_confidence_pass' (bool; True if Confidence != 'Negative')

    Raises:
        KeyError: If any calibration key is missing from a non-empty run. A
            partial config is not supported input: silently substituting a
            different threshold would change genotypes without a log line.
    """
    logger.debug("Entering calculate_depth_score_and_assign_confidence")
    logger.debug(f"Initial DataFrame shape: {df.shape}")

    if df.empty:
        logger.debug("DataFrame is empty. Exiting function without changes.")
        return df

    # Extract relevant config subdict.
    #
    # The calibration constants below are REQUIRED, deliberately. They used to be read
    # as `.get(key, <default>)`, and those defaults were not neutral placeholders: they
    # encoded a second, wrong calibration (0.4 for `high` where the shipped config
    # supplies 0.00515, a factor of 78). Dropping a single key was enough to genotype a
    # cohort against thresholds nobody chose, silently and with no log line -- an
    # ample-depth High_Precision* call came back as Low_Precision. `--config-path`
    # replaces the whole config rather than merging it, and a partial config already
    # fails with KeyError elsewhere in the pipeline, so it is not supported input and
    # the right behaviour here is to fail loudly rather than to recalibrate.
    conf_assign = kestrel_config["confidence_assignment"]
    thresholds = conf_assign["depth_score_thresholds"]
    alt_thresholds = conf_assign["alt_depth_thresholds"]
    var_region_threshold = conf_assign["var_active_region_threshold"]

    # Confidence LABELS keep their fallbacks: they are cosmetic strings, not
    # calibration, so a missing one cannot move a call.
    confidence_levels = conf_assign.get("confidence_levels", {})
    low_prec_label = confidence_levels.get("low_precision", "Low_Precision")
    high_prec_label = confidence_levels.get("high_precision", "High_Precision")
    high_prec_star_label = confidence_levels.get("high_precision_star", "High_Precision*")

    low_threshold = thresholds["low"]
    high_threshold = thresholds["high"]

    alt_low = alt_thresholds["low"]
    alt_mid_low = alt_thresholds["mid_low"]
    alt_mid_high = alt_thresholds["mid_high"]

    # Convert depth columns to numeric for arithmetic (using float)
    depth_cols = [
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
    ]
    for col in depth_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    # Step 2: Calculate Depth_Score (avoid division by zero -> np.nan)
    df["Depth_Score"] = df["Estimated_Depth_AlternateVariant"] / df["Estimated_Depth_Variant_ActiveRegion"]
    # Replace any infinite values (resulting from division by zero) with NaN
    df["Depth_Score"] = df["Depth_Score"].replace([np.inf, -np.inf], np.nan)

    # Step 3: Assign Confidence
    # Default all rows to NEGATIVE_LABEL, then upgrade based on conditions.
    # Variants with Depth_Score < low_threshold remain Negative (filtered out).
    df["Confidence"] = NEGATIVE_LABEL

    # Guard: only variants at or above the minimum depth threshold can receive
    # a non-Negative confidence label. This ensures all variants (GG and non-GG)
    # with insufficient depth are excluded from downstream analysis (fixes #147).
    above_min_threshold = df["Depth_Score"] >= low_threshold

    # Condition 1: Low Precision if Depth_Score == low_threshold OR region depth <= var_region_threshold
    cond1 = (df["Depth_Score"] == low_threshold) | (df["Estimated_Depth_Variant_ActiveRegion"] <= var_region_threshold)

    # Condition 2: High Precision STAR if alt depth >= alt_mid_high AND Depth_Score >= high_threshold
    cond2 = (df["Estimated_Depth_AlternateVariant"] >= alt_mid_high) & (df["Depth_Score"] >= high_threshold)

    # Condition 3: Low Precision if alt_depth is between mid_low (inclusive) and mid_high (exclusive),
    # and Depth_Score between low_threshold and high_threshold.
    # Retained as the expression that names the intent -- @hassansaei on #184: "do not
    # remove this intent". `cond_midband` below now subsumes it (same label, wider
    # alt-depth range), so this assignment no longer decides any row on its own.
    cond3 = df["Estimated_Depth_AlternateVariant"].between(
        alt_mid_low,
        alt_mid_high,
        inclusive="left",
    ) & df["Depth_Score"].between(low_threshold, high_threshold)

    # Condition 4: Low Precision if alt depth <= alt_low
    cond4 = df["Estimated_Depth_AlternateVariant"] <= alt_low

    # Condition 5: High Precision if alt_depth is between mid_low and mid_high and Depth_Score >= high_threshold
    cond5 = df["Estimated_Depth_AlternateVariant"].between(alt_mid_low, alt_mid_high, inclusive="left") & (
        df["Depth_Score"] >= high_threshold
    )

    # Mid-band demotion (#184): the whole CLOSED interval [low_threshold, high_threshold]
    # is Low_Precision, at every alternate depth. Applied last so it takes precedence
    # over the High tiers.
    #
    # @hassansaei on #184: "Any variant with Depth_Score between 0.00469 and 0.00515
    # (inclusive) must be Low_Precision, even when alternate depth is >= 21 (or higher).
    # That mid-range Depth_Score demotion from 1.3 is still the desired behaviour."
    #
    # This subsumes the former `cond6`, which covered the OPEN interval and was already
    # applied last -- so the interior of the band was correct and only Depth_Score ==
    # high_threshold exactly diverged, where cond2/cond5 promoted the row. Note the
    # alternative of changing cond2/cond5 to `> high_threshold` was rejected: on its own
    # it would leave a row at exactly high_threshold matching no condition at all, so it
    # would keep NEGATIVE_LABEL -- turning a reported call into a non-call.
    cond_midband = df["Depth_Score"].between(low_threshold, high_threshold, inclusive="both")

    # Apply conditions in order (later conditions can overwrite earlier ones).
    # The above_min_threshold guard ensures Depth_Score < low_threshold stays Negative.
    df.loc[cond1 & above_min_threshold, "Confidence"] = low_prec_label
    df.loc[cond2 & above_min_threshold, "Confidence"] = high_prec_star_label
    df.loc[cond3 & above_min_threshold, "Confidence"] = low_prec_label
    df.loc[cond4 & above_min_threshold, "Confidence"] = low_prec_label
    df.loc[cond5 & above_min_threshold, "Confidence"] = high_prec_label
    df.loc[cond_midband & above_min_threshold, "Confidence"] = low_prec_label

    # Step 4: Mark pass/fail: Passing means Confidence != NEGATIVE_LABEL
    df["depth_confidence_pass"] = df["Confidence"] != NEGATIVE_LABEL

    logger.debug("Exiting calculate_depth_score_and_assign_confidence")
    logger.debug(f"Final DataFrame shape: {df.shape}")
    return df
