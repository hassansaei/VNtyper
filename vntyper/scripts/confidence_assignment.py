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
import pandas as pd
import numpy as np  # Import NumPy directly


def calculate_depth_score_and_assign_confidence(df: pd.DataFrame, kestrel_config: dict) -> pd.DataFrame:
    """
    Calculates Depth_Score for each variant and assigns a confidence label
    based on coverage thresholds defined in kestrel_config['confidence_assignment'].

    Args:
        df (pd.DataFrame):
            DataFrame with columns:
              - 'Estimated_Depth_AlternateVariant' (alt depth)
              - 'Estimated_Depth_Variant_ActiveRegion' (region depth)
        kestrel_config (dict):
            Contains subdict 'confidence_assignment', which has:
              - 'depth_score_thresholds'
              - 'alt_depth_thresholds'
              - 'var_active_region_threshold'
              - 'confidence_levels'

    Returns:
        pd.DataFrame:
            Same DataFrame with two new columns:
              - 'Depth_Score' (float)
              - 'Confidence' (str: 'Low_Precision', 'High_Precision', etc.)
    """
    logging.debug("Entering calculate_depth_score_and_assign_confidence")
    logging.debug(f"Initial DataFrame shape: {df.shape}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting function.")
        return df

    # Step 1: Convert depth columns to integers
    depth_cols = ['Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']
    df[depth_cols] = df[depth_cols].astype(int)
    logging.debug("Converted depth columns to integers.")

    # Step 2: Calculate Depth_Score
    df['Depth_Score'] = df['Estimated_Depth_AlternateVariant'] / df['Estimated_Depth_Variant_ActiveRegion']
    logging.debug("Calculated 'Depth_Score'.")

    # Step 3: Assign confidence
    conf_assign = kestrel_config['confidence_assignment']
    thresholds = conf_assign['depth_score_thresholds']
    alt_thresholds = conf_assign['alt_depth_thresholds']
    var_region_threshold = conf_assign['var_active_region_threshold']
    confidence_levels = conf_assign['confidence_levels']

    conditions = [
        # Condition 1: Low Precision
        (df['Depth_Score'] <= thresholds['low']) | (df['Estimated_Depth_Variant_ActiveRegion'] <= var_region_threshold),

        # Condition 2: High Precision Star (Updated to include alt_depth >= mid_high)
        (df['Estimated_Depth_AlternateVariant'] >= alt_thresholds['mid_high']) &
        (df['Depth_Score'] >= thresholds['high']),

        # Condition 3: Low Precision
        (df['Estimated_Depth_AlternateVariant'].between(alt_thresholds['mid_low'], alt_thresholds['mid_high'])) &
        (df['Depth_Score'].between(thresholds['low'], thresholds['high'])),

        # Condition 4: Low Precision
        (df['Estimated_Depth_AlternateVariant'] <= alt_thresholds['low']),

        # Condition 5: High Precision
        (df['Estimated_Depth_AlternateVariant'].between(alt_thresholds['mid_low'], alt_thresholds['mid_high'], inclusive='left')) &
        (df['Depth_Score'] >= thresholds['high'])
    ]

    choices = [
        confidence_levels['low_precision'],
        confidence_levels['high_precision_star'],  # Updated to assign 'High_Precision*'
        confidence_levels['low_precision'],
        confidence_levels['low_precision'],
        confidence_levels['high_precision']
    ]

    # Ensure 'high_precision_star' exists in confidence_levels
    if 'high_precision_star' not in confidence_levels:
        logging.error("'high_precision_star' not found in confidence_levels.")
        raise KeyError("'high_precision_star' not found in confidence_levels.")

    df['Confidence'] = np.select(conditions, choices, default=confidence_levels.get('low_precision', 'Low_Precision'))
    logging.debug("Assigned 'Confidence' labels based on conditions.")

    logging.debug(f"Final DataFrame shape: {df.shape}")
    logging.debug("Exiting calculate_depth_score_and_assign_confidence")
    return df
