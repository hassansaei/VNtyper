#!/usr/bin/env python3
"""
scoring.py

Module Purpose:
---------------
Contains functions for frame-score calculations, depth splitting from the
"Sample" field, and other variant scoring or extraction logic.

Typical Flow:
-------------
1. Split the 'Sample' column (colon-delimited) into relevant depths.
2. Calculate the frame score = (alt_len - ref_len) / 3, used to
   detect frameshift vs. non-frameshift.
3. (Refactored) Instead of filtering out non-frameshift variants,
   we now add a boolean column indicating whether it is frameshift.
4. Further split the frame score to identify insertion/deletion patterns
   (like 3n+1 vs 3n+2). Rather than dropping columns or filtering,
   we add boolean columns that mark relevant conditions.
5. Return the DataFrame with added columns. Final filtering
   is deferred to a later step.

References:
-----------
- Thresholds, heuristics from Saei et al., iScience 26, 107171 (2023).
"""

import logging
import pandas as pd
import numpy as np


def split_depth_and_calculate_frame_score(df: pd.DataFrame) -> pd.DataFrame:
    """
    Splits the 'Sample' column to obtain alternate/depth coverage,
    then calculates a 'Frame_Score' = (alt_len - ref_len) / 3.

    Steps:
      1) df['Sample'].split(':') -> e.g., Del:AltDepth:ActiveDepth
      2) Re-map to new columns (Estimated_Depth_AlternateVariant, etc.).
      3) Compute frame difference: (ALT length - REF length) / 3.
         Determine if it's a frameshift variant based on modulo operation.

    Args:
        df (pd.DataFrame):
            Must contain columns: ['Sample', 'REF', 'ALT'].

    Returns:
        pd.DataFrame:
            Adds 'Frame_Score' and a boolean column 'is_frameshift'.
            Does NOT filter out rows.
    """
    logging.debug("Entering split_depth_and_calculate_frame_score")
    logging.debug(f"Initial row count: {len(df)}, columns: {df.columns.tolist()}")

    if df.empty:
        logging.debug(
            "DataFrame is empty. Exiting split_depth_and_calculate_frame_score."
        )
        return df

    # Step 1) Split 'Sample' into 3 parts
    split_columns = df["Sample"].str.split(":", expand=True)
    split_columns.columns = [
        "Del",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
    ]
    df = pd.concat([df, split_columns], axis=1)

    # Step 2) Compute lengths and frame score
    df["ref_len"] = df["REF"].str.len()
    df["alt_len"] = df["ALT"].str.len()
    df["Frame_Score"] = (df["alt_len"] - df["ref_len"]) / 3

    # Step 3) Mark frameshift in a new boolean column
    df["is_frameshift"] = (df["alt_len"] - df["ref_len"]) % 3 != 0

    logging.debug("Exiting split_depth_and_calculate_frame_score")
    logging.debug(f"Final row count: {len(df)}, columns: {df.columns.tolist()}")
    return df


def split_frame_score(df: pd.DataFrame) -> pd.DataFrame:
    """
    Splits 'Frame_Score' into 'direction' and 'frameshift_amount'
    for further logic (3n+1 vs. 3n+2, etc.).

    Steps:
      1) Determine the direction based on the sign of (alt_len - ref_len).
      2) Calculate frameshift_amount as abs((alt_len - ref_len) % 3).
      3) (Refactored) We do not drop columns used for intermediate steps.

    Args:
        df (pd.DataFrame):
            Should contain 'Frame_Score', 'is_frameshift', 'ref_len', 'alt_len'.

    Returns:
        pd.DataFrame:
            Adds 'direction' and 'frameshift_amount'.
            Retains all rows and intermediate columns for debugging.
    """
    logging.debug("Entering split_frame_score")
    logging.debug(f"Initial row count: {len(df)}, columns: {df.columns.tolist()}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting split_frame_score.")
        return df

    # Step 1) Determine direction
    df["direction"] = np.sign(df["alt_len"] - df["ref_len"])

    # Step 2) Calculate frameshift_amount
    df["frameshift_amount"] = (df["alt_len"] - df["ref_len"]).abs() % 3

    logging.debug("Exiting split_frame_score")
    logging.debug(f"Final row count: {len(df)}, columns: {df.columns.tolist()}")
    return df


def extract_frameshifts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts frameshift variants that follow known insertion/deletion patterns.

    Steps:
      - For insertion frameshifts: direction > 0 & frameshift_amount == 1
      - For deletion frameshifts: direction < 0 & frameshift_amount == 2

    (Refactored) We do NOT remove non-matching rows here. Instead, we add a
    boolean column 'is_valid_frameshift' to indicate whether a row meets one
    of the frameshift patterns.

    Args:
        df (pd.DataFrame):
            Must contain 'direction' and 'frameshift_amount' columns.

    Returns:
        pd.DataFrame:
            Adds 'is_valid_frameshift' (bool).
            Keeps all rows.
    """
    logging.debug("Entering extract_frameshifts")
    logging.debug(f"Initial row count: {len(df)}, columns: {df.columns.tolist()}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting extract_frameshifts.")
        return df

    # Identify insertion vs deletion frameshifts
    condition_insertion = (df["direction"] > 0) & (df["frameshift_amount"] == 1)
    condition_deletion = (df["direction"] < 0) & (df["frameshift_amount"] == 2)

    df["is_valid_frameshift"] = condition_insertion | condition_deletion

    logging.debug("Exiting extract_frameshifts")
    logging.debug(f"Final row count: {len(df)}, columns: {df.columns.tolist()}")
    return df
