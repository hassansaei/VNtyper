#!/usr/bin/env python3
# tests/unit/test_scoring.py

"""
Unit tests for the scoring functionality in vntyper/scripts/scoring.py.
Validates frame-score calculations, depth splitting, and frameshift extraction.
"""

import pytest
import pandas as pd
import numpy as np

from vntyper.scripts.scoring import (
    split_depth_and_calculate_frame_score,
    split_frame_score,
    extract_frameshifts,
)


@pytest.mark.parametrize("df_input,expected_len", [
    (pd.DataFrame(), 0),
])
def test_split_depth_and_calculate_frame_score_empty_df(df_input, expected_len):
    """
    Verify that an empty input DataFrame remains empty.
    """
    out = split_depth_and_calculate_frame_score(df_input)
    assert len(out) == expected_len, (
        "Empty input should yield empty output after split_depth_and_calculate_frame_score."
    )


def test_split_depth_and_calculate_frame_score_no_frameshift():
    """
    If the difference (ALT length - REF length) is a multiple of 3,
    the variant should not be retained (non-frameshift).
    """
    df = pd.DataFrame({
        "Sample": ["Del:10:100"],   # Only the first 'Del' part is not used, but we keep format for test
        "REF": ["ATG"],            # length 3
        "ALT": ["ATGATG"],         # length 6  -> difference = 3 -> multiple of 3
        "Motifs": ["mock_motif"],
        "Variant": ["mock_variant"],
        "POS": [123],
        "Motif_sequence": ["mock_sequence"]
    })
    out = split_depth_and_calculate_frame_score(df)
    # Because it's a multiple of 3 difference, is_frameshift == False => filtered out
    assert out.empty, (
        "Variants with multiple-of-3 difference should be filtered out as non-frameshift."
    )


def test_split_depth_and_calculate_frame_score_frameshift():
    """
    If the difference (ALT length - REF length) is not a multiple of 3,
    the variant should be retained and a 'Frame_Score' should be added.
    """
    df = pd.DataFrame({
        "Sample": ["Del:50:500"],
        "REF": ["ATG"],            # length 3
        "ALT": ["ATGA"],           # length 4 -> difference = 1 -> frameshift
        "Motifs": ["mock_motif"],
        "Variant": ["mock_variant"],
        "POS": [456],
        "Motif_sequence": ["mock_sequence"]
    })
    out = split_depth_and_calculate_frame_score(df)
    assert not out.empty, "Expected to retain a frameshift variant (difference not multiple of 3)."
    assert "Frame_Score" in out.columns, "Output should have a 'Frame_Score' column."
    # Check that is_frameshift was True
    assert "is_frameshift" in out.columns, "Output should have 'is_frameshift' marking frameshift or not."
    assert all(out["is_frameshift"]), "All retained rows should be frameshift variants."


def test_split_frame_score_empty_df():
    """
    Verify that an empty input DataFrame remains empty when split_frame_score is called.
    """
    df = pd.DataFrame()
    out = split_frame_score(df)
    assert out.empty, "Empty input should yield empty output after split_frame_score."


def test_split_frame_score_basic():
    """
    Test basic splitting of frame score into 'direction' and 'frameshift_amount'.
    """
    df = pd.DataFrame({
        "Frame_Score": [1.0, -2.0],  # not directly used, but indicates frameshift
        "ref_len": [3, 6],
        "alt_len": [4, 4],          # alt_len - ref_len => [1, -2]
        "is_frameshift": [True, True]  # frameshift is assumed True from previous step
    })
    out = split_frame_score(df)

    # We drop 'is_frameshift', 'ref_len', 'alt_len'
    # We keep 'direction', 'frameshift_amount', 'Frame_Score', etc.
    expected_columns = {"Frame_Score", "direction", "frameshift_amount"}
    assert expected_columns.issubset(set(out.columns)), (
        f"Output must contain at least: {expected_columns}"
    )

    # direction = sign(alt_len - ref_len)
    # frameshift_amount = abs(alt_len - ref_len) % 3
    # For row0: alt_len - ref_len = 1 => direction=1, frameshift_amount=1
    # For row1: alt_len - ref_len = -2 => direction < 0 => -1, frameshift_amount=2
    assert out.loc[0, "direction"] == 1, "Expected direction=1 for alt_len-ref_len=1."
    assert out.loc[0, "frameshift_amount"] == 1, "Expected frameshift_amount=1 for difference=1."

    assert out.loc[1, "direction"] == -1, "Expected direction=-1 for alt_len-ref_len=-2."
    assert out.loc[1, "frameshift_amount"] == 2, "Expected frameshift_amount=2 for difference=-2."


def test_extract_frameshifts_empty_df():
    """
    Verify that an empty input DataFrame remains empty in extract_frameshifts.
    """
    df = pd.DataFrame()
    out = extract_frameshifts(df)
    assert out.empty, "Empty input should yield empty output after extract_frameshifts."


def test_extract_frameshifts_mixed():
    """
    Test that only frameshift rows meeting the 3n+1 insertion (direction>0) or
    3n+2 deletion (direction<0) are retained.
    """
    df = pd.DataFrame({
        "direction": [1, 1, -1, -1, 1],
        "frameshift_amount": [1, 2, 2, 1, 1],
        "Frame_Score": [0.33, 0.66, -0.66, -0.33, 0.33],
        "Variant": ["ins_ok", "ins_wrong", "del_ok", "del_wrong", "ins_ok2"]
    })
    # insertion frameshift => direction>0 and frameshift_amount=1
    # deletion frameshift => direction<0 and frameshift_amount=2
    out = extract_frameshifts(df)
    # Expect to retain rows: 0 (ins_ok), 2 (del_ok), 4 (ins_ok2)
    # Indices 1 (ins_wrong) and 3 (del_wrong) should be dropped
    assert len(out) == 3, "Expected to keep 3 rows that match frameshift patterns."
    assert sorted(out["Variant"].tolist()) == sorted(["ins_ok", "del_ok", "ins_ok2"]), (
        "Wrong set of frameshift variants retained."
    )
