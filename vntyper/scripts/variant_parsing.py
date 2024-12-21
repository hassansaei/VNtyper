#!/usr/bin/env python3
"""
variant_parsing.py

Module Purpose:
---------------
Handles reading and parsing VCF files into DataFrames, along with final
ALT-based filtering rules (e.g., excluding certain ALT sequences or
checking coverage thresholds for 'GG').

Typical Flow:
-------------
1. read_vcf_without_comments():
   - Reads a VCF file (gzipped or not),
   - Ignores lines starting with "##" (meta),
   - Takes the "#CHROM" line as the header.

2. filter_by_alt_values_and_finalize():
   - If 'GG' is present in ALT, ensure Depth_Score >= threshold or remove it.
   - Exclude other undesired ALT sequences from kestrel_config['alt_filtering'].
   - Drop columns used only for intermediate steps like 'left', 'right'.

References:
-----------
- Saei et al., iScience 26, 107171 (2023)
"""

import logging
import pandas as pd


def read_vcf_without_comments(vcf_file):
    """
    Reads a VCF file (possibly gzipped) ignoring lines starting with '##'.
    The line starting with '#CHROM' is used to define DataFrame columns.

    Args:
        vcf_file (str):
            Path to the VCF file, e.g., "output_insertion.vcf" or "output.vcf.gz".

    Returns:
        pd.DataFrame:
            Contains the main variant records. May be empty if the file has no variants.
    """
    import gzip

    open_func = gzip.open if vcf_file.endswith(".gz") else open
    data = []
    header_line = None

    with open_func(vcf_file, 'rt') as f:
        for line in f:
            # The "##" lines are meta lines to be skipped
            if line.startswith("#CHROM"):
                header_line = line.strip().split('\t')
            elif not line.startswith("##"):
                data.append(line.strip().split('\t'))

    if data:
        return pd.DataFrame(data, columns=header_line)
    else:
        return pd.DataFrame()


def filter_by_alt_values_and_finalize(df, kestrel_config):
    """
    Applies final filtering rules based on ALT values, e.g., removing certain
    ALTs or requiring a minimal Depth_Score if ALT='GG'.

    Steps:
      1) If 'GG' is present in the ALT column, only keep those rows if
         Depth_Score >= 'gg_depth_score_threshold'.
      2) Exclude any ALTs in 'exclude_alts'.
      3) Drop 'left' and 'right' columns from earlier frame splitting steps
         to finalize the output.

    Args:
        df (pd.DataFrame):
            Must contain columns 'ALT', 'Depth_Score'.
        kestrel_config (dict):
            Must contain 'alt_filtering' subdict with keys:
              - 'gg_alt_value' (e.g., "GG")
              - 'gg_depth_score_threshold'
              - 'exclude_alts'

    Returns:
        pd.DataFrame: Filtered, finalized DataFrame ready for final steps.
    """
    if df.empty:
        return df

    alt_filter = kestrel_config['alt_filtering']
    gg_alt_value = alt_filter['gg_alt_value']
    gg_depth_score_threshold = alt_filter['gg_depth_score_threshold']
    exclude_alts = alt_filter['exclude_alts']

    # If 'GG' is in ALT, only retain it if Depth_Score >= threshold
    if df['ALT'].str.contains(r'\b' + gg_alt_value + r'\b').any():
        gg_condition = df['ALT'] == gg_alt_value
        df = pd.concat(
            [
                df[~gg_condition],
                df[gg_condition & (df['Depth_Score'] >= gg_depth_score_threshold)],
            ]
        )

    # Exclude specified ALTs
    df = df[~df['ALT'].isin(exclude_alts)]

    # Drop intermediate columns from frame splitting
    drop_cols = [col for col in ['left', 'right'] if col in df.columns]
    df.drop(drop_cols, axis=1, inplace=True)

    return df
