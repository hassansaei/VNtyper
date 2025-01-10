#!/usr/bin/env python3
"""
scoring.py

Module Purpose:
---------------
Contains functions for frame-score calculations, depth splitting from the
"Sample" field, and other variant scoring or extraction logic.

Typical Flow (Refactored):
--------------------------
1. Split the 'Sample' column (colon-delimited) into relevant depths.
2. Calculate a numeric frame difference = (alt_len - ref_len)/3, 
   indicating potential frameshifts.
3. Mark each variant with booleans (e.g., is_frameshift, is_insertion_frameshift).
4. Preserve all rows, returning a DataFrame with extra annotation columns.
5. Perform final filtering outside these functions (e.g., in process_kmer_results).

References:
-----------
- Thresholds, heuristics from Saei et al., iScience 26, 107171 (2023).
"""

import logging
import pandas as pd

logger = logging.getLogger(__name__)


def split_depth_and_calculate_frame_score(df):
    """
    Splits the 'Sample' column to obtain alternate/depth coverage,
    then calculates a raw frame difference = (ALT length - REF length) / 3.
    Also flags potential frameshifts as a boolean (is_frameshift).

    Steps (refactored):
      1) df['Sample'].split(':') -> e.g., Del:AltDepth:ActiveDepth
      2) Map these parts to new columns:
         'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion'.
      3) Compute numeric frame difference: (alt_len - ref_len)/3,
         store as df['raw_frame_diff'].
      4) Create df['is_frameshift'] = True if raw_frame_diff is not a multiple of 3.
      5) Keep all rows; no final filtering here.

    Args:
        df (pd.DataFrame):
            Must contain columns: ['Sample', 'REF', 'ALT'].

    Returns:
        pd.DataFrame:
            Adds/retains columns:
              - Del
              - Estimated_Depth_AlternateVariant
              - Estimated_Depth_Variant_ActiveRegion
              - ref_len, alt_len
              - raw_frame_diff
              - is_frameshift
            Leaves filtering to later steps.
    """
    if df.empty:
        logger.debug("[split_depth_and_calculate_frame_score] DataFrame is empty. Returning as is.")
        return df

    logger.debug("[split_depth_and_calculate_frame_score] Starting with %d rows.", len(df))

    # 1) Split the Sample column
    split_cols = df['Sample'].str.split(':', expand=True)
    df['Del'] = split_cols[0]
    df['Estimated_Depth_AlternateVariant'] = split_cols[1]
    df['Estimated_Depth_Variant_ActiveRegion'] = split_cols[2]

    # 2) (Optionally) ensure we keep relevant columns but do NOT drop others yet.
    #    If you want to reorder or ensure columns exist, do so carefully:
    #    (No forced subset here—keep the entire df for debugging.)

    # 3) Compute lengths, then frame difference
    df['ref_len'] = df['REF'].str.len()
    df['alt_len'] = df['ALT'].str.len()
    df['raw_frame_diff'] = (df['alt_len'] - df['ref_len']) / 3.0

    # 4) Mark frameshift if raw_frame_diff is NOT a multiple of 3
    #    We can do a small rounding to avoid floating-precision issues.
    df['is_frameshift'] = (df['raw_frame_diff'].round(5) % 1 != 0)

    frameshift_count = df['is_frameshift'].sum()
    logger.debug(
        "[split_depth_and_calculate_frame_score] Found %d frameshift variants among %d rows.",
        frameshift_count, len(df)
    )

    return df


def split_frame_score(df):
    """
    Further annotates the frame difference by splitting it into "left" and "right"
    parts (e.g., converting numeric raw_frame_diff into a string like "-1.33"
    and then splitting).

    Steps:
      1) Convert df['raw_frame_diff'] to a string column 'Frame_Score_str' (rounded).
      2) df['Frame_Score_str'].split('.', 2) -> left, right
         (e.g., "-1.33" → ["-1", "33"]).
      3) If 'left' is '-0', replace with '-1' to standardize.
      4) Preserve all columns, adding 'left' and 'right' for later logic.

    Args:
        df (pd.DataFrame):
            Should contain 'raw_frame_diff' from previous step.

    Returns:
        pd.DataFrame:
            Now has columns 'Frame_Score_str', 'left', 'right'
            in addition to existing columns.
    """
    if df.empty:
        logger.debug("[split_frame_score] DataFrame is empty. Returning as is.")
        return df

    logger.debug("[split_frame_score] Starting with %d rows.", len(df))

    # 1) Make a string column for splitting
    #    Example: "-1.33"
    df['Frame_Score_str'] = df['raw_frame_diff'].apply(lambda x: f"{x:.2f}")

    # 2) Split into left/right
    df[['left', 'right']] = df['Frame_Score_str'].str.split('.', expand=True)

    # 3) Standardize '-0' to '-1'
    df['left'] = df['left'].replace('-0', '-1')

    logger.debug(
        "[split_frame_score] Completed splitting. Sample left/right values: %s / %s",
        df['left'].iloc[0] if len(df) > 0 else "N/A",
        df['right'].iloc[0] if len(df) > 0 else "N/A"
    )

    # Do NOT drop columns here; keep them for debugging.

    return df


def extract_frameshifts(df):
    """
    Annotates frameshift variants that follow known insertion/deletion patterns.

    We introduce:
      - df['is_insertion_frameshift']: True if (left >= 0) & right contains '33'
      - df['is_deletion_frameshift']:  True if (left < 0) & right contains '67'

    Steps:
      1) Evaluate conditions for insertion vs. deletion frameshifts.
      2) Store booleans, do NOT filter out rows.

    Args:
        df (pd.DataFrame):
            Must contain 'left' and 'right' columns (strings).

    Returns:
        pd.DataFrame:
            Adds columns:
              - is_insertion_frameshift
              - is_deletion_frameshift
            Leaves filtering to later steps.
    """
    if df.empty:
        logger.debug("[extract_frameshifts] DataFrame is empty. Returning as is.")
        return df

    logger.debug("[extract_frameshifts] Starting with %d rows.", len(df))

    # 1) Mark insertion frameshift
    #    (non-negative left + '33' in right)
    def is_insertion(l_str, r_str):
        # Example: left = "1", right = "33"
        try:
            return float(l_str) >= 0 and '33' in r_str
        except ValueError:
            # If l_str not convertible to float, fallback False
            return False

    # 2) Mark deletion frameshift
    #    (negative left + '67' in right)
    def is_deletion(l_str, r_str):
        try:
            return float(l_str) < 0 and '67' in r_str
        except ValueError:
            return False

    df['is_insertion_frameshift'] = df.apply(
        lambda row: is_insertion(row['left'], row['right']), axis=1
    )
    df['is_deletion_frameshift'] = df.apply(
        lambda row: is_deletion(row['left'], row['right']), axis=1
    )

    ins_count = df['is_insertion_frameshift'].sum()
    del_count = df['is_deletion_frameshift'].sum()
    logger.debug(
        "[extract_frameshifts] Found %d insertion frameshifts and %d deletion frameshifts among %d rows.",
        ins_count, del_count, len(df)
    )

    # 3) We do NOT filter. We keep all rows with flags for later filtering.
    return df
