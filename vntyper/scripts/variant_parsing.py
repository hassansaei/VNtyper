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

import gzip
import logging
from typing import Optional

import pandas as pd


def read_vcf_without_comments(vcf_file: str) -> pd.DataFrame:
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
    open_func = gzip.open if vcf_file.endswith(".gz") else open
    data = []
    header: Optional[list] = None

    try:
        with open_func(vcf_file, "rt") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().split("\t")
                elif not line.startswith("##") and header:
                    data.append(line.strip().split("\t"))
    except FileNotFoundError:
        logging.error(f"VCF file not found: {vcf_file}")
        return pd.DataFrame()
    except Exception as e:
        logging.error(f"Error reading VCF file {vcf_file}: {e}")
        return pd.DataFrame()

    if data and header:
        logging.debug(f"VCF read successfully with {len(data)} records.")
        return pd.DataFrame(data, columns=header)
    logging.debug("No variant records found in VCF.")
    return pd.DataFrame()


def filter_by_alt_values_and_finalize(df: pd.DataFrame, kestrel_config: dict) -> pd.DataFrame:
    """
    Applies final filtering rules based on ALT values, e.g., removing certain
    ALTs or requiring a minimal Depth_Score if ALT='GG'.

    (Refactored)
      - Instead of removing rows, add a boolean column 'alt_filter_pass'
        that marks whether each row passes these checks.

    Args:
        df (pd.DataFrame):
            Must contain columns 'ALT', 'Depth_Score'.
        kestrel_config (dict):
            Must contain 'alt_filtering' subdict with keys:
              - 'gg_alt_value' (e.g., "GG")
              - 'gg_depth_score_threshold'
              - 'exclude_alts'

    Returns:
        pd.DataFrame:
            Same shape as input, with a new 'alt_filter_pass' column (bool).
            Intermediate columns like 'left', 'right' are no longer dropped;
            they remain for debugging.
    """
    logging.debug("Entering filter_by_alt_values_and_finalize")
    logging.debug(f"Initial DataFrame shape: {df.shape}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting filter_by_alt_values_and_finalize.")
        return df

    required_columns = {"ALT", "Depth_Score"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        logging.error(f"Missing required columns: {missing_columns}")
        raise KeyError(f"Missing required columns: {missing_columns}")

    alt_filter = kestrel_config.get("alt_filtering", {})
    gg_alt_value = alt_filter.get("gg_alt_value", "GG")
    gg_depth_threshold = alt_filter.get("gg_depth_score_threshold", 0.0)
    exclude_alts = alt_filter.get("exclude_alts", [])

    # Ensure Depth_Score is float
    df["Depth_Score"] = pd.to_numeric(df["Depth_Score"], errors="coerce")

    # Create the boolean mask
    is_gg = df["ALT"] == gg_alt_value
    meets_gg_threshold = df["Depth_Score"] >= gg_depth_threshold
    not_excluded_alt = ~df["ALT"].isin(exclude_alts)

    # Instead of filtering the DataFrame, we store a new boolean column
    df["alt_filter_pass"] = (~is_gg | meets_gg_threshold) & not_excluded_alt

    logging.debug("Exiting filter_by_alt_values_and_finalize")
    logging.debug(f"Final DataFrame shape: {df.shape}")
    return df
