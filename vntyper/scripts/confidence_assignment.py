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

import pandas as pd


def calculate_depth_score_and_assign_confidence(df, kestrel_config):
    """
    Calculates Depth_Score for each variant and assigns a confidence label
    based on coverage thresholds defined in kestrel_config['confidence_assignment'].

    This step is crucial for distinguishing between likely real variants
    vs. potential false positives. Variants below a certain coverage ratio
    or having fewer than a certain alt-depth are labeled "Low_Precision".

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
    if df.empty:
        return df

    # Convert depth columns from string to integer
    df['Estimated_Depth_AlternateVariant'] = df['Estimated_Depth_AlternateVariant'].astype(int)
    df['Estimated_Depth_Variant_ActiveRegion'] = (
        df['Estimated_Depth_Variant_ActiveRegion'].astype(int)
    )

    # Depth_Score = alt_depth / total_region_depth
    df['Depth_Score'] = (
        df['Estimated_Depth_AlternateVariant'] / df['Estimated_Depth_Variant_ActiveRegion']
    )

    # Retrieve threshold settings from kestrel_config
    conf_assign = kestrel_config['confidence_assignment']
    depth_score_thresholds = conf_assign['depth_score_thresholds']
    alt_depth_thresholds = conf_assign['alt_depth_thresholds']
    var_active_region_threshold = conf_assign['var_active_region_threshold']
    confidence_levels = conf_assign['confidence_levels']

    def assign_confidence(row):
        depth_score = row['Depth_Score']
        alt_depth = row['Estimated_Depth_AlternateVariant']
        var_region = row['Estimated_Depth_Variant_ActiveRegion']

        # Compare depth_score, alt_depth, and var_region to thresholds
        if (
            depth_score <= depth_score_thresholds['low']
            or var_region <= var_active_region_threshold
        ):
            return confidence_levels['low_precision']
        elif (
            alt_depth_thresholds['mid_low']
            <= alt_depth
            <= alt_depth_thresholds['mid_high']
            and depth_score_thresholds['low']
            <= depth_score
            <= depth_score_thresholds['high']
        ):
            return confidence_levels['low_precision']
        elif alt_depth > alt_depth_thresholds['mid_high']:
            return confidence_levels['high_precision']
        elif alt_depth <= alt_depth_thresholds['low']:
            return confidence_levels['low_precision']
        elif (
            alt_depth_thresholds['mid_low']
            <= alt_depth
            < alt_depth_thresholds['mid_high']
            and depth_score >= depth_score_thresholds['high']
        ):
            return confidence_levels['high_precision']
        elif (
            alt_depth >= alt_depth_thresholds['mid_high']
            and depth_score >= depth_score_thresholds['high']
        ):
            return confidence_levels['high_precision_star']
        else:
            return confidence_levels['low_precision']

    df['Confidence'] = df.apply(assign_confidence, axis=1)
    return df
