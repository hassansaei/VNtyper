"""Module: confidence_assignment.py

Description: Calculates depth scores and assigns confidence levels to Kestrel variants
             using an ordered, first-match rule table (vntyper.scripts.confidence_rules).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from vntyper.scripts.confidence_rules import (
    NEGATIVE_LABEL,
    assign_confidence_labels,
)

logger = logging.getLogger(__name__)

HAS_REPORTING_FLOOR_SPLIT = True


def calculate_depth_score_and_assign_confidence(df: pd.DataFrame, kestrel_config: dict) -> pd.DataFrame:
    """Calculate Depth_Score and assign confidence levels to variants based on Kestrel config.

    Evaluates confidence using the ordered first-match rule table in confidence_rules.py.

    Args:
        df (pandas.DataFrame): DataFrame containing kmer results with columns:
                              'Estimated_Depth_AlternateVariant',
                              'Estimated_Depth_Variant_ActiveRegion'
        kestrel_config (dict): Kestrel component of the active decision profile,
                               including the 'confidence_assignment' section.

    Returns:
        pandas.DataFrame: DataFrame with added 'Depth_Score', 'Confidence', and
                          'depth_confidence_pass' columns.

    Raises:
        KeyError: If any required calibration threshold key is missing.
    """
    logger.debug("Entering calculate_depth_score_and_assign_confidence")
    logger.debug(f"Input DataFrame shape: {df.shape}")

    # Empty DataFrame check: short-circuit before accessing kestrel_config
    if df.empty:
        logger.warning("Empty DataFrame provided to calculate_depth_score_and_assign_confidence")
        return df

    # Enforce required calibration keys (fails loudly with KeyError if missing)
    conf_assign = kestrel_config["confidence_assignment"]
    _ = conf_assign["reporting_floor"]
    thresholds = conf_assign["depth_score_thresholds"]
    alt_thresholds = conf_assign["alt_depth_thresholds"]
    var_region_threshold = conf_assign["var_active_region_threshold"]

    _ = thresholds["low"]
    _ = thresholds["high"]
    _ = alt_thresholds["low"]
    _ = alt_thresholds["mid_low"]
    _ = alt_thresholds["mid_high"]
    _ = var_region_threshold

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

    # Step 3: Assign Confidence using the ordered rule table
    df["Confidence"] = assign_confidence_labels(df, conf_assign)

    # Step 4: Mark pass/fail: Passing means Confidence != NEGATIVE_LABEL
    df["depth_confidence_pass"] = df["Confidence"] != NEGATIVE_LABEL

    logger.debug("Exiting calculate_depth_score_and_assign_confidence")
    logger.debug(f"Final DataFrame shape: {df.shape}")
    return df
