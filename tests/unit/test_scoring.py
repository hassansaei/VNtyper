#!/usr/bin/env python3
# tests/unit/test_scoring.py

"""
Unit tests for the scoring functionality in vntyper/scripts/scoring.py.

Validates frame-score calculations, depth splitting, and frameshift extraction.
"""

import pytest
import pandas as pd

from vntyper.scripts.scoring import (
    split_depth_and_calculate_frame_score,
    split_frame_score,
    extract_frameshifts,
)


@pytest.mark.parametrize(
    "df_input,expected_len",
    [
        (pd.DataFrame(), 0),
    ],
)
def test_split_depth_and_calculate_frame_score_empty_df(df_input, expected_len):
    """Verify that an empty input DataFrame remains empty."""
    out = split_depth_and_calculate_frame_score(df_input)
    assert (
        len(out) == expected_len
    ), "Empty input should yield empty output after split_depth_and_calculate_frame_score."


def test_split_depth_and_calculate_frame_score_no_frameshift():
    """
    Test that non-frameshift variants (multiple of 3 difference) are correctly flagged.

    When the difference (ALT length - REF length) is a multiple of 3,
    the variant is marked with is_frameshift=False but retained in the DataFrame.
    This follows the flag-and-defer pattern where filtering is done downstream.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "Sample": [
                "Del:10:100"
            ],  # Only the first 'Del' part is not used, but we keep format for test
            "REF": ["ATG"],  # length 3
            "ALT": ["ATGATG"],  # length 6  -> difference = 3 -> multiple of 3
            "Motifs": ["mock_motif"],
            "Variant": ["mock_variant"],
            "POS": [123],
            "Motif_sequence": ["mock_sequence"],
        }
    )

    # Act
    out = split_depth_and_calculate_frame_score(df)

    # Assert - row is retained with is_frameshift flag
    assert len(out) == 1, "Row should be retained in the DataFrame"
    assert "is_frameshift" in out.columns, "Should have 'is_frameshift' column"
    assert not out.loc[0, "is_frameshift"], \
        "Variant with multiple-of-3 difference should be marked as is_frameshift=False"

    # Assert - Frame_Score calculation is correct
    assert "Frame_Score" in out.columns, "Should have 'Frame_Score' column"
    expected_frame_score = (6 - 3) / 3  # (alt_len - ref_len) / 3 = 1.0
    assert (
        out.loc[0, "Frame_Score"] == expected_frame_score
    ), f"Frame_Score should be {expected_frame_score}"


def test_split_depth_and_calculate_frame_score_frameshift():
    """
    Test frameshift variant detection and Frame_Score calculation.

    If the difference (ALT length - REF length) is not a multiple of 3,
    the variant should be retained and a 'Frame_Score' should be added.
    """
    df = pd.DataFrame(
        {
            "Sample": ["Del:50:500"],
            "REF": ["ATG"],  # length 3
            "ALT": ["ATGA"],  # length 4 -> difference = 1 -> frameshift
            "Motifs": ["mock_motif"],
            "Variant": ["mock_variant"],
            "POS": [456],
            "Motif_sequence": ["mock_sequence"],
        }
    )
    out = split_depth_and_calculate_frame_score(df)
    assert (
        not out.empty
    ), "Expected to retain a frameshift variant (difference not multiple of 3)."
    assert "Frame_Score" in out.columns, "Output should have a 'Frame_Score' column."
    # Check that is_frameshift was True
    assert (
        "is_frameshift" in out.columns
    ), "Output should have 'is_frameshift' marking frameshift or not."
    assert all(out["is_frameshift"]), "All retained rows should be frameshift variants."


def test_split_frame_score_empty_df():
    """Verify that an empty input DataFrame remains empty when split_frame_score is called."""
    df = pd.DataFrame()
    out = split_frame_score(df)
    assert out.empty, "Empty input should yield empty output after split_frame_score."


def test_split_frame_score_basic():
    """Test basic splitting of frame score into 'direction' and 'frameshift_amount'."""
    df = pd.DataFrame(
        {
            "Frame_Score": [1.0, -2.0],  # not directly used, but indicates frameshift
            "ref_len": [3, 6],
            "alt_len": [4, 4],  # alt_len - ref_len => [1, -2]
            "is_frameshift": [
                True,
                True,
            ],  # frameshift is assumed True from previous step
        }
    )
    out = split_frame_score(df)

    # We drop 'is_frameshift', 'ref_len', 'alt_len'
    # We keep 'direction', 'frameshift_amount', 'Frame_Score', etc.
    expected_columns = {"Frame_Score", "direction", "frameshift_amount"}
    assert expected_columns.issubset(
        set(out.columns)
    ), f"Output must contain at least: {expected_columns}"

    # direction = sign(alt_len - ref_len)
    # frameshift_amount = abs(alt_len - ref_len) % 3
    # For row0: alt_len - ref_len = 1 => direction=1, frameshift_amount=1
    # For row1: alt_len - ref_len = -2 => direction < 0 => -1, frameshift_amount=2
    assert out.loc[0, "direction"] == 1, "Expected direction=1 for alt_len-ref_len=1."
    assert (
        out.loc[0, "frameshift_amount"] == 1
    ), "Expected frameshift_amount=1 for difference=1."

    assert (
        out.loc[1, "direction"] == -1
    ), "Expected direction=-1 for alt_len-ref_len=-2."
    assert (
        out.loc[1, "frameshift_amount"] == 2
    ), "Expected frameshift_amount=2 for difference=-2."


def test_extract_frameshifts_empty_df():
    """Verify that an empty input DataFrame remains empty in extract_frameshifts."""
    df = pd.DataFrame()
    out = extract_frameshifts(df)
    assert out.empty, "Empty input should yield empty output after extract_frameshifts."


def test_extract_frameshifts_mixed():
    """
    Test that frameshift variants are correctly flagged based on insertion/deletion patterns.

    Valid frameshift patterns:
    - Insertion frameshift: direction > 0 AND frameshift_amount == 1 (3n+1)
    - Deletion frameshift: direction < 0 AND frameshift_amount == 2 (3n+2)

    All rows are retained but marked with is_valid_frameshift flag.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "direction": [1, 1, -1, -1, 1],
            "frameshift_amount": [1, 2, 2, 1, 1],
            "Frame_Score": [0.33, 0.66, -0.66, -0.33, 0.33],
            "Variant": ["ins_ok", "ins_wrong", "del_ok", "del_wrong", "ins_ok2"],
        }
    )

    # Act
    out = extract_frameshifts(df)

    # Assert - all rows are retained
    assert len(out) == 5, "All 5 rows should be retained in the DataFrame"

    # Assert - is_valid_frameshift column exists
    assert (
        "is_valid_frameshift" in out.columns
    ), "Should have 'is_valid_frameshift' column"

    # Assert - correct rows are marked as valid frameshifts
    # Row 0: direction=1, frameshift_amount=1 => valid insertion (3n+1)
    assert out.loc[0, "is_valid_frameshift"], \
        "Row 0 (ins_ok) should be valid: direction>0 and frameshift_amount=1"

    # Row 1: direction=1, frameshift_amount=2 => invalid (wrong frameshift amount for insertion)
    assert not out.loc[1, "is_valid_frameshift"], \
        "Row 1 (ins_wrong) should be invalid: direction>0 but frameshift_amount=2"

    # Row 2: direction=-1, frameshift_amount=2 => valid deletion (3n+2)
    assert out.loc[2, "is_valid_frameshift"], \
        "Row 2 (del_ok) should be valid: direction<0 and frameshift_amount=2"

    # Row 3: direction=-1, frameshift_amount=1 => invalid (wrong frameshift amount for deletion)
    assert not out.loc[3, "is_valid_frameshift"], \
        "Row 3 (del_wrong) should be invalid: direction<0 but frameshift_amount=1"

    # Row 4: direction=1, frameshift_amount=1 => valid insertion (3n+1)
    assert out.loc[4, "is_valid_frameshift"], \
        "Row 4 (ins_ok2) should be valid: direction>0 and frameshift_amount=1"

    # Assert - verify the count of valid frameshifts
    valid_count = out["is_valid_frameshift"].sum()
    assert valid_count == 3, "Should have exactly 3 valid frameshifts (rows 0, 2, 4)"
