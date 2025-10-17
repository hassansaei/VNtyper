#!/usr/bin/env python3
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
1. Convert depth columns to integers.
2. Compute Depth_Score = alt / active_region depth.
3. Compare Depth_Score and alt-depth coverage to threshold intervals.
4. Assign textual confidence labels (e.g., 'Low_Precision', 'High_Precision').
5. Return updated DataFrame with 'Depth_Score' and 'Confidence' columns.

References:
-----------
- Saei et al., iScience 26, 107171 (2023)
- #6: Docs: Kestrel postprocessing heuristic (depth score usage)
"""

import logging

import numpy as np  # Import NumPy directly
import pandas as pd


def calculate_depth_score_and_assign_confidence(df: pd.DataFrame, kestrel_config: dict) -> pd.DataFrame:
    """
    Calculates Depth_Score for each variant and assigns a confidence label
    based on coverage thresholds defined in kestrel_config['confidence_assignment'].

    (Refactored)
      - All rows remain in the DataFrame.
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
            Must contain a subdict 'confidence_assignment' with keys:
              - 'depth_score_thresholds' (e.g. {'low': 0.2, 'high': 0.4})
              - 'alt_depth_thresholds' (e.g. {'low': 5, 'mid_low': 10, 'mid_high': 20})
              - 'var_active_region_threshold' (e.g. 30)
              - 'confidence_levels' (dict with keys like
                'low_precision', 'high_precision', 'high_precision_star')

    Returns:
        pd.DataFrame:
            Same shape as input, with:
              - 'Depth_Score' (float)
              - 'Confidence' (str, e.g., 'Low_Precision', 'High_Precision', etc.)
              - 'depth_confidence_pass' (bool; True if Confidence != 'Negative')
    """
    logging.debug("Entering calculate_depth_score_and_assign_confidence")
    logging.debug(f"Initial DataFrame shape: {df.shape}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting function without changes.")
        return df

    # Extract relevant config subdict
    conf_assign = kestrel_config.get("confidence_assignment", {})
    thresholds = conf_assign.get("depth_score_thresholds", {})
    alt_thresholds = conf_assign.get("alt_depth_thresholds", {})
    var_region_threshold = conf_assign.get("var_active_region_threshold", 0)
    confidence_levels = conf_assign.get("confidence_levels", {})

    # Fallback for confidence levels
    low_prec_label = confidence_levels.get("low_precision", "Low_Precision")
    high_prec_label = confidence_levels.get("high_precision", "High_Precision")
    high_prec_star_label = confidence_levels.get("high_precision_star", "High_Precision*")

    # Default to 0.2 for 'low' and 0.4 for 'high' if not specified
    low_threshold = thresholds.get("low", 0.2)
    high_threshold = thresholds.get("high", 0.4)

    # Default alt-depth thresholds (in case they are missing)
    alt_low = alt_thresholds.get("low", 5)
    alt_mid_low = alt_thresholds.get("mid_low", 10)
    alt_mid_high = alt_thresholds.get("mid_high", 20)

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
    # Default all rows to 'Negative', then update based on conditions.
    df["Confidence"] = "Negative"

    # Condition 1: Low Precision if Depth_Score <= low_threshold OR region depth <= var_region_threshold
    cond1 = (df["Depth_Score"] <= low_threshold) | (df["Estimated_Depth_Variant_ActiveRegion"] <= var_region_threshold)

    # Condition 2: High Precision STAR if alt depth >= alt_mid_high AND Depth_Score >= high_threshold
    cond2 = (df["Estimated_Depth_AlternateVariant"] >= alt_mid_high) & (df["Depth_Score"] >= high_threshold)

    # Condition 3: Low Precision if alt_depth is between mid_low and mid_high,
    # and Depth_Score between low_threshold and high_threshold
    cond3 = df["Estimated_Depth_AlternateVariant"].between(alt_mid_low, alt_mid_high) & df["Depth_Score"].between(
        low_threshold, high_threshold
    )

    # Condition 4: Low Precision if alt depth <= alt_low
    cond4 = df["Estimated_Depth_AlternateVariant"] <= alt_low

    # Condition 5: High Precision if alt_depth is between mid_low and mid_high and Depth_Score >= high_threshold
    cond5 = df["Estimated_Depth_AlternateVariant"].between(alt_mid_low, alt_mid_high, inclusive="left") & (
        df["Depth_Score"] >= high_threshold
    )

    # Apply conditions in order (later conditions can overwrite earlier ones)
    df.loc[cond1, "Confidence"] = low_prec_label
    df.loc[cond2, "Confidence"] = high_prec_star_label
    df.loc[cond3, "Confidence"] = low_prec_label
    df.loc[cond4, "Confidence"] = low_prec_label
    df.loc[cond5, "Confidence"] = high_prec_label

    # Step 4: Mark pass/fail: Passing means Confidence != 'Negative'
    df["depth_confidence_pass"] = df["Confidence"] != "Negative"

    logging.debug("Exiting calculate_depth_score_and_assign_confidence")
    logging.debug(f"Final DataFrame shape: {df.shape}")
    return df
