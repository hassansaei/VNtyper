#!/usr/bin/env python3
# vntyper/scripts/motif_processing.py

"""
motif_processing.py

Module Purpose:
---------------
Handles all logic for working with MUC1 VNTR motifs:
  - Loading references from FASTA
  - Preprocessing insertion/deletion variants by merging them with motif info
  - Performing final motif correction & annotation based on the Kestrel
    config thresholds (position_threshold, exclude lists).

Typical Flow:
-------------
1. Load the main MUC1 reference to map 'Motif' IDs to sequences.
2. Preprocess insertion/deletion variants by renaming columns and merging
   them with the reference.
3. Possibly load additional motifs from 'MUC1_motifs_Rev_com.fa' for annotation.
4. Correct annotation by splitting left/right motifs, dropping duplicates,
   applying threshold-based logic from the config.

References:
-----------
- Saei et al., iScience 26, 107171 (2023).
"""

import logging

import pandas as pd
from Bio import SeqIO


def load_muc1_reference(reference_file):
    """
    Loads the MUC1 VNTR reference motifs from a FASTA file into a DataFrame.

    Args:
        reference_file (str):
            Path to the FASTA file containing MUC1 VNTR reference motifs.
            Each record ID is treated as 'Motifs', and the sequence as 'Motif_sequence'.

    Returns:
        pd.DataFrame:
            Columns: ['Motifs', 'Motif_sequence']
    """
    identifiers = []
    sequences = []

    with open(reference_file) as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            identifiers.append(seq_record.id)
            sequences.append(seq_record.seq)

    return pd.DataFrame({"Motifs": identifiers, "Motif_sequence": sequences})


def preprocessing_insertion(df, muc1_ref):
    """
    Preprocess insertion variants by merging them with the MUC1 reference motifs.

    Steps:
      - Rename '#CHROM' → 'Motifs'
      - Drop unused columns (ID, QUAL, FILTER, INFO, FORMAT)
      - Rename last column → 'Sample'
      - Merge with 'muc1_ref' to link motif IDs

    Args:
        df (pd.DataFrame): Insertion variants from a filtered VCF.
        muc1_ref (pd.DataFrame): MUC1 reference DataFrame with 'Motifs' & 'Motif_sequence'.

    Returns:
        pd.DataFrame: Updated with columns for 'Variant' = 'Insertion'.
    """
    df.rename(columns={"#CHROM": "Motifs"}, inplace=True)
    columns_to_drop = ["ID", "QUAL", "FILTER", "INFO", "FORMAT"]
    df.drop(columns=columns_to_drop, axis=1, inplace=True)
    last_column_name = df.columns[-1]
    df.rename(columns={last_column_name: "Sample"}, inplace=True)
    df = pd.merge(df, muc1_ref, on="Motifs", how="left")
    df["Variant"] = "Insertion"
    return df


def preprocessing_deletion(df, muc1_ref):
    """
    Preprocess deletion variants by merging them with the MUC1 reference motifs.

    Steps:
      - Rename '#CHROM' → 'Motifs'
      - Drop unused columns
      - Merge with MUC1 reference
      - Mark 'Variant' = 'Deletion'

    Args:
        df (pd.DataFrame): Deletion variants from a filtered VCF.
        muc1_ref (pd.DataFrame): MUC1 reference DataFrame with 'Motifs' & 'Motif_sequence'.

    Returns:
        pd.DataFrame: Updated with columns for 'Variant' = 'Deletion'.
    """
    df.rename(columns={"#CHROM": "Motifs"}, inplace=True)
    columns_to_drop = ["ID", "QUAL", "FILTER", "INFO", "FORMAT"]
    df.drop(columns=columns_to_drop, axis=1, inplace=True)
    last_column_name = df.columns[-1]
    df.rename(columns={last_column_name: "Sample"}, inplace=True)
    df = pd.merge(df, muc1_ref, on="Motifs", how="left")
    df["Variant"] = "Deletion"
    return df


def load_additional_motifs(config):
    """
    Load additional motifs from a FASTA file (e.g., 'MUC1_motifs_Rev_com.fa').

    Typically used to annotate the final results with reversed/complemented motifs
    or other expansions of the MUC1 motif set.

    Args:
        config (dict):
            The pipeline config containing
            ["reference_data"]["muc1_motifs_rev_com"].

    Returns:
        pd.DataFrame:
            Columns: ['Motif', 'Motif_sequence']
    """
    from Bio import SeqIO

    identifiers = []
    sequences = []

    muc1_motifs_rev_com_file = config["reference_data"]["muc1_motifs_rev_com"]

    with open(muc1_motifs_rev_com_file) as motif_file:
        for seq_record in SeqIO.parse(motif_file, "fasta"):
            identifiers.append(seq_record.id)
            sequences.append(str(seq_record.seq.upper()))

    return pd.DataFrame({"Motif": identifiers, "Motif_sequence": sequences})


def _apply_uniform_filtering_right_motif(motif_right, exclude_motifs_right, alt_for_motif_right_gg, motifs_for_alt_gg):
    """
    Apply uniform depth-score-based filtering for right motif variants.

    This function implements the fix for Issue #136 by using depth-score-first
    filtering instead of ALT-based filtering. This prevents silent deletion of
    non-GG frameshift variants (e.g., insG_pos54 mutations).

    Algorithm:
    1. Remove excluded motifs FIRST (conserved motifs where variants are likely artifacts)
    2. Sort by Depth_Score DESC, then POS DESC (highest coverage = highest confidence)
    3. Deduplicate by [POS, REF, ALT] - CRITICAL: position-specific to preserve different biological events
    4. Preserve motif-specific GG logic (keeps GG only in allowed motifs, preserves ALL non-GG)

    Args:
        motif_right (pd.DataFrame): Right motif variants (POS >= position_threshold)
        exclude_motifs_right (list): Motifs to exclude (e.g., conserved motifs)
        alt_for_motif_right_gg (str): ALT value for GG filtering (typically "GG")
        motifs_for_alt_gg (list): Allowed motifs for GG variants (typically ["X"])

    Returns:
        pd.DataFrame: Filtered variants with uniform depth-score-based filtering

    Biological Rationale:
        - insG_pos54: 176 GC variants + 85 GG variants → keep highest-depth GC (not delete all GC!)
        - dupC: G>GG, C>CG, T>TG at same position → keep highest-depth variant
        - Position-specific dedup: POS 60 vs POS 67 are DIFFERENT biological events

    References:
        - GitHub Issue #136
        - Hassan Saei's recommendation: "Sort by depth score, deduplicate, not ALT-based filtering"
    """
    if motif_right.empty:
        return motif_right

    # Step 1: Remove excluded motifs FIRST
    # Biological: Conserved motifs (Q, 8, 9, 7, etc.) are stable → variants likely artifacts
    motif_right = motif_right[~motif_right["Motif"].isin(exclude_motifs_right)]

    if motif_right.empty:
        return motif_right

    # Step 2: Sort by Depth_Score DESC, then POS DESC
    # Biological: Highest coverage first (signal > noise), then position for consistency
    motif_right = motif_right.sort_values(["Depth_Score", "POS"], ascending=[False, False])

    # Step 3: Deduplicate by [POS, REF, ALT] - CRITICAL FIX!
    # Why position-specific?
    #   - POS 60: C>GC (insG at pos 54) vs POS 67: G>GG (dupC) = DIFFERENT biological events
    #   - Old bug: subset=["REF", "ALT"] would lose one of them
    #   - New fix: subset=["POS", "REF", "ALT"] preserves both
    motif_right = motif_right.drop_duplicates(
        subset=["POS", "REF", "ALT"],
        keep="first",  # Keeps highest Depth_Score (sorted first)
    )

    # Step 4: Preserve motif-specific GG logic
    # Biological: GG variants only allowed in motif "X" (canonical repeat unit)
    #             All non-GG variants (GC, CG, CT, etc.) are PRESERVED
    if (motif_right["ALT"] == alt_for_motif_right_gg).any():
        # Filter GG: keep only if in allowed motifs
        gg_in_allowed = motif_right[
            (motif_right["ALT"] == alt_for_motif_right_gg) & (motif_right["Motif"].isin(motifs_for_alt_gg))
        ]
        # Preserve ALL non-GG variants (this is the FIX for Issue #136!)
        non_gg = motif_right[motif_right["ALT"] != alt_for_motif_right_gg]
        # Combine
        motif_right = pd.concat([gg_in_allowed, non_gg], ignore_index=True)

    return motif_right


def _prioritize_frameshift_and_dedupe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort by frameshift validity and depth score, then deduplicate on genomic locus.

    This helper function standardizes frameshift-aware deduplication across all
    motif filtering paths (left motifs, right motifs uniform, right motifs legacy).
    It ensures valid frameshift variants are prioritized over non-frameshift variants
    when multiple variants share the same genomic position.

    Algorithm:
        1. Fill missing `is_valid_frameshift` values with False (defensive programming)
        2. Sort by: is_valid_frameshift DESC, Depth_Score DESC, POS DESC
           - Frameshift-valid variants first (biologically significant)
           - Then highest depth score (signal > noise)
           - Then position for deterministic ordering
        3. Deduplicate on [POS, REF, ALT] to preserve distinct genomic events
           - Why not ALT-only? POS 60 C>GC vs POS 67 G>GG are DIFFERENT events

    Args:
        df: DataFrame with motif variants. Must have columns:
            - POS, REF, ALT (required for deduplication)
            - Depth_Score (required for sorting)
            - is_valid_frameshift (optional, filled with False if missing)

    Returns:
        DataFrame sorted and deduplicated with frameshift prioritization.

    Example:
        Two variants at same POS/REF/ALT:
        - Row 1: is_valid_frameshift=True, Depth_Score=0.010
        - Row 2: is_valid_frameshift=False, Depth_Score=0.015
        Result: Row 1 kept (frameshift=True prioritized despite lower depth)

    References:
        - GitHub Issue #136: Frameshift detection fixes
        - Sourcery AI review: DRY principle recommendation
    """
    if df.empty:
        return df

    # Create a copy to avoid modifying the original DataFrame
    result = df.copy()

    # Determine sort columns based on whether is_valid_frameshift exists
    if "is_valid_frameshift" in result.columns:
        # Fill missing values with False (variant is not a valid frameshift)
        result["is_valid_frameshift"] = result["is_valid_frameshift"].fillna(False)
        sort_cols = ["is_valid_frameshift", "Depth_Score", "POS"]
        sort_ascending = [False, False, False]  # True first, highest score, highest POS
    else:
        sort_cols = ["Depth_Score", "POS"]
        sort_ascending = [False, False]

    # Sort to prioritize valid frameshifts with highest depth scores
    result = result.sort_values(sort_cols, ascending=sort_ascending)

    # Deduplicate on genomic locus [POS, REF, ALT] - preserves distinct biological events
    # Why [POS, REF, ALT] instead of ALT-only?
    #   - POS 60: C>GC (insG at pos 54) vs POS 67: G>GG (dupC) = DIFFERENT events
    #   - ALT-only dedup would incorrectly merge these distinct variants
    result = result.drop_duplicates(subset=["POS", "REF", "ALT"], keep="first")

    return result


def motif_correction_and_annotation(df, merged_motifs, kestrel_config):
    """
    Final step of motif annotation: correct positions for left/right motifs,
    drop duplicates, handle special cases (e.g., 'GG' alt on the right motif).

    (Refactored):
      - Returns the same shape as `df`,
      - Adds 'motif_filter_pass' boolean,
      - Ensures 'Motif_fasta', 'POS_fasta', and 'Motif' columns
        exist in the final output, even for failing rows.

    Issue #136 Fix:
      - Supports new uniform filtering via 'use_uniform_filtering' flag
      - Default: false (legacy behavior for backward compatibility)
      - When true: applies depth-score-based filtering (fixes insG_pos54 detection)
    """
    logging.debug("Entering motif_correction_and_annotation")
    logging.debug(f"Initial row count: {len(df)}, columns: {df.columns.tolist()}")

    if df.empty:
        logging.debug("DataFrame is empty. Exiting motif_correction_and_annotation.")
        df["motif_filter_pass"] = False
        df["Motif_fasta"] = pd.NA
        df["POS_fasta"] = pd.NA
        df["Motif"] = pd.NA
        return df

    # Keep a copy to attach pass/fail
    original_df = df.copy(deep=True)
    original_df["original_index"] = original_df.index

    mf = kestrel_config["motif_filtering"]
    use_uniform_filtering = mf.get("use_uniform_filtering", False)  # Issue #136 fix flag
    position_threshold = mf.get("position_threshold", 60)
    exclude_motifs_right = mf.get("exclude_motifs_right", [])
    alt_for_motif_right_gg = mf.get("alt_for_motif_right_gg", "GG")
    motifs_for_alt_gg = mf.get("motifs_for_alt_gg", [])
    exclude_alts_combined = mf.get("exclude_alts_combined", [])
    exclude_motifs_combined = mf.get("exclude_motifs_combined", [])

    # =============== Original Logic ===============
    working_df = original_df.copy(deep=True)

    # Step 1) Ensure 'Motif_fasta'; check for dash
    if "Motifs" in working_df.columns:
        working_df["Motif_fasta"] = working_df["Motifs"]
    else:
        logging.error("Missing 'Motifs' column. Old code returns empty.")
        combined_df = pd.DataFrame(columns=working_df.columns)
        pass

    if "Motifs" not in working_df.columns or working_df["Motifs"].str.count("-").max() != 1:
        logging.error("Cannot split 'Motifs' into left-right. Old code returns empty.")
        combined_df = pd.DataFrame(columns=working_df.columns)
    else:
        working_df[["Motif_left", "Motif_right"]] = working_df["Motifs"].str.split("-", expand=True)
        working_df["POS"] = pd.to_numeric(working_df["POS"], errors="coerce").fillna(-1).astype(int)

        # Left vs. Right
        motif_left = working_df[working_df["POS"] < position_threshold].copy()
        motif_right = working_df[working_df["POS"] >= position_threshold].copy()

        # Merge + filter left
        if not motif_left.empty:
            motif_left.rename(columns={"Motif_right": "Motif"}, inplace=True)
            if "Motif_sequence" in motif_left.columns:
                motif_left.drop(columns=["Motif_sequence"], inplace=True, errors="ignore")
            motif_left = motif_left.merge(merged_motifs, on="Motif", how="left")

            keep_cols = [
                "Motif",
                "Motif_fasta",
                "Variant",
                "POS",
                "REF",
                "ALT",
                "Motif_sequence",
                "Estimated_Depth_AlternateVariant",
                "Estimated_Depth_Variant_ActiveRegion",
                "Depth_Score",
                "Confidence",
                "original_index",
            ]
            if "is_valid_frameshift" in motif_left.columns:
                keep_cols.append("is_valid_frameshift")
            motif_left = motif_left[keep_cols]

            # Apply frameshift-aware sorting and deduplication (DRY: uses shared helper)
            motif_left = _prioritize_frameshift_and_dedupe(motif_left)

        # Merge + filter right
        if not motif_right.empty:
            motif_right.rename(columns={"Motif_left": "Motif"}, inplace=True)
            if "Motif_sequence" in motif_right.columns:
                motif_right.drop(columns=["Motif_sequence"], inplace=True, errors="ignore")
            motif_right = motif_right.merge(merged_motifs, on="Motif", how="left")

            keep_cols = [
                "Motif",
                "Motif_fasta",
                "Variant",
                "POS",
                "REF",
                "ALT",
                "Motif_sequence",
                "Estimated_Depth_AlternateVariant",
                "Estimated_Depth_Variant_ActiveRegion",
                "Depth_Score",
                "Confidence",
                "original_index",
            ]
            if "is_valid_frameshift" in motif_right.columns:
                keep_cols.append("is_valid_frameshift")
            motif_right = motif_right[keep_cols]

            # Issue #136 Fix: Branch based on use_uniform_filtering flag
            if use_uniform_filtering:
                # NEW: Uniform depth-score-based filtering (fixes insG_pos54 detection)
                motif_right = _apply_uniform_filtering_right_motif(
                    motif_right, exclude_motifs_right, alt_for_motif_right_gg, motifs_for_alt_gg
                )
                # Re-apply frameshift prioritization after uniform filtering (DRY: uses shared helper)
                motif_right = _prioritize_frameshift_and_dedupe(motif_right)
            else:
                # IMPROVED LEGACY: Hassan's refactored GG logic (PR #140)
                # Better than old logic but still has limitations vs uniform filtering
                if motif_right["ALT"].str.contains(r"\b" + alt_for_motif_right_gg + r"\b").any():
                    motif_right = motif_right[~motif_right["Motif"].isin(exclude_motifs_right)]
                    # Apply frameshift-aware sorting and deduplication (DRY: uses shared helper)
                    motif_right = _prioritize_frameshift_and_dedupe(motif_right)
                    if motif_right["Motif"].isin(motifs_for_alt_gg).any():
                        motif_right = motif_right[motif_right["Motif"].isin(motifs_for_alt_gg)]

        # Combine
        combined_df = pd.concat([motif_right, motif_left], axis=0, ignore_index=True)
        combined_df = combined_df[~combined_df["ALT"].isin(exclude_alts_combined)]
        combined_df = combined_df[~combined_df["Motif"].isin(exclude_motifs_combined)]

        # Adjust POS => create POS_fasta
        combined_df["POS"] = pd.to_numeric(combined_df["POS"], errors="coerce").fillna(-1).astype(int)
        combined_df["POS_fasta"] = combined_df["POS"]
        combined_df.update(
            combined_df["POS"].mask(
                combined_df["POS"] >= position_threshold,
                lambda x: x - position_threshold,
            )
        )
    # =============== End Original Logic ===============

    # Mark pass/fail based on original_index
    pass_mask = original_df["original_index"].isin(combined_df.get("original_index", []))
    original_df["motif_filter_pass"] = pass_mask

    # Ensure final columns exist in the main DF even for failing rows
    # (They remain NaN if the row didn't pass.)
    if "Motif_fasta" not in original_df.columns:
        original_df["Motif_fasta"] = pd.NA
    if "POS_fasta" not in original_df.columns:
        original_df["POS_fasta"] = pd.NA
    if "Motif" not in original_df.columns:
        original_df["Motif"] = pd.NA

    # For rows that "passed", copy over the final POS_fasta, Motif_fasta, and Motif
    combined_df = combined_df.set_index("original_index", drop=False)
    for idx in combined_df.index:
        if idx in original_df.index:
            # Copy Motif_fasta / POS_fasta
            original_df.at[idx, "Motif_fasta"] = combined_df.at[idx, "Motif_fasta"]
            if "POS_fasta" in combined_df.columns:
                original_df.at[idx, "POS_fasta"] = combined_df.at[idx, "POS_fasta"]
            # Also copy 'Motif'
            if "Motif" in combined_df.columns:
                original_df.at[idx, "Motif"] = combined_df.at[idx, "Motif"]

    # Drop the temporary original_index column
    original_df.drop(columns=["original_index"], inplace=True, errors="ignore")

    logging.debug("Exiting motif_correction_and_annotation")
    logging.debug(f"Final row count: {len(original_df)}, columns: {original_df.columns.tolist()}")
    return original_df
