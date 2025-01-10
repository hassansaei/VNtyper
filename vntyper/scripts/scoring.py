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
3. Filter out non-frameshift variants.
4. Further split the frame score to identify insertion/deletion patterns
   (like 3n+1 vs 3n+2).
5. Return the filtered DataFrame that should then move on to confidence assignment.

References:
-----------
- Thresholds, heuristics from Saei et al., iScience 26, 107171 (2023).
"""

import logging
import pandas as pd


def split_depth_and_calculate_frame_score(df):
    """
    Splits the 'Sample' column to obtain alternate/depth coverage,
    then calculates a 'Frame_Score' = (alt_len - ref_len)/3.

    Steps:
      1) df['Sample'].split(':') -> e.g., Del:AltDepth:ActiveDepth
      2) Re-map to new columns (Estimated_Depth_AlternateVariant, etc.).
      3) Compute frame difference: (ALT length - REF length)/3.
         Round to two decimals, converting e.g. "3.0" -> "3C" to denote
         a multiple of 3.

    Args:
        df (pd.DataFrame):
            Must contain columns: ['Sample', 'REF', 'ALT'].

    Returns:
        pd.DataFrame:
            Adds 'Frame_Score' and keeps only frameshift
            (non-"C") variants. The others are filtered out.
    """
    logging.debug("Entering split_depth_and_calculate_frame_score")
    logging.debug(f"Initial row count: {len(df)}, columns: {df.columns.tolist()}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting split_depth_and_calculate_frame_score.")
        return df

    # Step 1) Split the Sample column into 3 parts
    pre_split_rows = len(df)
    pre_split_cols = df.columns.tolist()
    df[['Del', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']] = (
        df['Sample'].str.split(':', expand=True)
    )
    logging.debug("After splitting 'Sample':")
    logging.debug(f"Changed from {pre_split_rows} rows, {pre_split_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    # Step 2) Keep only necessary columns in a new DataFrame
    pre_select_rows = len(df)
    pre_select_cols = df.columns.tolist()
    df = df[
        [
            'Motifs',
            'Variant',
            'POS',
            'REF',
            'ALT',
            'Motif_sequence',
            'Estimated_Depth_AlternateVariant',
            'Estimated_Depth_Variant_ActiveRegion',
        ]
    ].copy()
    logging.debug("After selecting necessary columns:")
    logging.debug(f"Changed from {pre_select_rows} rows, {pre_select_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    # Step 3) Compute lengths, then frame score
    pre_frame_rows = len(df)
    pre_frame_cols = df.columns.tolist()
    df["ref_len"] = df["REF"].str.len()
    df["alt_len"] = df["ALT"].str.len()
    df["Frame_Score"] = (
        (df["alt_len"] - df["ref_len"]) / 3
    ).round(2).astype(str).apply(lambda x: x.replace('.0', 'C'))
    logging.debug("After computing 'Frame_Score':")
    logging.debug(f"Changed from {pre_frame_rows} rows, {pre_frame_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    # Step 4) Mark non-frameshift
    pre_truefalse_rows = len(df)
    pre_truefalse_cols = df.columns.tolist()
    df["TrueFalse"] = df['Frame_Score'].str.contains('C', regex=True)
    logging.debug("After marking non-frameshift variants:")
    logging.debug(f"Changed from {pre_truefalse_rows} rows, {pre_truefalse_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    # Keep only rows where Frame_Score does not contain 'C'
    pre_filter_rows = len(df)
    pre_filter_cols = df.columns.tolist()
    df = df[df["TrueFalse"] == False].copy()
    logging.debug("After filtering out non-frameshift variants:")
    logging.debug(f"Changed from {pre_filter_rows} rows, {pre_filter_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    logging.debug("Exiting split_depth_and_calculate_frame_score")
    logging.debug(f"Final row count: {len(df)}, columns: {df.columns.tolist()}")
    return df


def split_frame_score(df):
    """
    Splits 'Frame_Score' into 'left' and 'right' parts (e.g., "-1.33" → ["-1","33"])
    for further logic in identifying insertion vs. deletion frameshifts.

    Steps:
      1) df['Frame_Score'].split('.', 2)
      2) If 'left' is '-0', replace with '-1' to standardize
      3) Drop intermediate columns (TrueFalse, ref_len, alt_len)

    Args:
        df (pd.DataFrame):
            Should contain 'Frame_Score', 'TrueFalse', 'ref_len', 'alt_len'.

    Returns:
        pd.DataFrame: Now has columns 'left' and 'right' for analyzing frames.
    """
    logging.debug("Entering split_frame_score")
    logging.debug(f"Initial row count: {len(df)}, columns: {df.columns.tolist()}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting split_frame_score.")
        return df

    # Step 1) Split 'Frame_Score'
    pre_split_rows = len(df)
    pre_split_cols = df.columns.tolist()
    df[['left', 'right']] = df['Frame_Score'].str.split('.', expand=True)
    logging.debug("After splitting 'Frame_Score':")
    logging.debug(f"Changed from {pre_split_rows} rows, {pre_split_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    # Step 2) Replace '-0' with '-1'
    pre_replace_rows = len(df)
    pre_replace_cols = df.columns.tolist()
    df['left'] = df['left'].replace('-0', '-1')
    logging.debug("After replacing '-0' with '-1' in 'left':")
    logging.debug(f"Changed from {pre_replace_rows} rows, {pre_replace_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    # Step 3) Drop intermediate columns
    pre_drop_rows = len(df)
    pre_drop_cols = df.columns.tolist()
    df.drop(['TrueFalse', 'ref_len', 'alt_len'], axis=1, inplace=True)
    logging.debug("After dropping intermediate columns:")
    logging.debug(f"Changed from {pre_drop_rows} rows, {pre_drop_cols} columns")
    logging.debug(f"To {len(df)} rows, {df.columns.tolist()} columns")

    logging.debug("Exiting split_frame_score")
    logging.debug(f"Final row count: {len(df)}, columns: {df.columns.tolist()}")
    return df


def extract_frameshifts(df):
    """
    Extracts frameshift variants that follow known insertion/deletion patterns.

    Steps:
      - For insertion frameshifts: (left not negative) & right = "33"
        (3n+1 logic).
      - For deletion frameshifts: (left negative) & right = "67"
        (3n+2 logic).

    Args:
        df (pd.DataFrame):
            Must contain 'left' and 'right' columns.

    Returns:
        pd.DataFrame:
            Subset of frameshift variants meeting insertion or deletion patterns.
    """
    logging.debug("Entering extract_frameshifts")
    logging.debug(f"Initial row count: {len(df)}, columns: {df.columns.tolist()}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting extract_frameshifts.")
        return df

    # Identify insertion frameshifts
    pre_ins_rows = len(df)
    pre_ins_cols = df.columns.tolist()
    ins = df[
        df["left"].apply(lambda x: '-' not in x) &
        df["right"].apply(lambda y: '33' in y)
    ]
    logging.debug("After identifying insertion frameshifts:")
    logging.debug(f"Sliced from {pre_ins_rows} rows, {pre_ins_cols} columns to {len(ins)} rows.")

    # Identify deletion frameshifts
    pre_del_rows = len(df)
    pre_del_cols = df.columns.tolist()
    del_ = df[
        df["left"].apply(lambda x: '-' in x) &
        df["right"].apply(lambda y: '67' in y)
    ]
    logging.debug("After identifying deletion frameshifts:")
    logging.debug(f"Sliced from {pre_del_rows} rows, {pre_del_cols} columns to {len(del_)} rows.")

    # Combine
    combined = pd.concat([ins, del_], axis=0)
    logging.debug("After concatenating insertion and deletion frameshifts:")
    logging.debug(f"Resulting row count: {len(combined)}, columns: {combined.columns.tolist()}")

    logging.debug("Exiting extract_frameshifts")
    return combined
