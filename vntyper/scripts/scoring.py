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
    if df.empty:
        return df

    # 1) Split the Sample column into 3 parts
    df[['Del', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']] = (
        df['Sample'].str.split(':', expand=True)
    )

    # 2) Keep only necessary columns in a new DataFrame
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

    # 3) Compute lengths, then frame score
    df["ref_len"] = df["REF"].str.len()
    df["alt_len"] = df["ALT"].str.len()

    # (alt_len - ref_len) / 3
    df["Frame_Score"] = (
        (df["alt_len"] - df["ref_len"]) / 3
    ).round(2).astype(str).apply(lambda x: x.replace('.0', 'C'))

    # Mark non-frameshift
    df["TrueFalse"] = df['Frame_Score'].str.contains('C', regex=True)
    # Keep only rows where TrueFalse == False (meaning not a multiple of 3)
    df = df[df["TrueFalse"] == False].copy()

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
    if df.empty:
        return df

    df[['left', 'right']] = df['Frame_Score'].str.split('.', expand=True)
    df['left'] = df['left'].replace('-0', '-1')

    df.drop(['TrueFalse', 'ref_len', 'alt_len'], axis=1, inplace=True)

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
    if df.empty:
        return df

    ins = df[
        df["left"].apply(lambda x: '-' not in x)
        & df["right"].apply(lambda y: '33' in y)
    ]
    del_ = df[
        df["left"].apply(lambda x: '-' in x)
        & df["right"].apply(lambda y: '67' in y)
    ]

    return pd.concat([ins, del_], axis=0)