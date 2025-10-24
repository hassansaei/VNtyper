#!/usr/bin/env python3
# tests/unit/test_variant_parsing.py

"""
Unit tests for variant_parsing.py, focusing on filter_by_alt_values_and_finalize().

Ensures that ALT-based filtering rules are correctly applied:
  - 'GG' alt requires a minimum Depth_Score.
  - exclude_alts removed from final DataFrame.
  - left/right columns dropped at the end.
"""

import pandas as pd
import pytest

from vntyper.scripts.variant_parsing import filter_by_alt_values_and_finalize

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


@pytest.fixture
def kestrel_config_mock():
    """
    Provide a mock kestrel_config dict specifically for ALT filtering tests.

    Real config may contain more fields, but we only need 'alt_filtering' here.
    """
    return {
        "alt_filtering": {
            "gg_alt_value": "GG",
            "gg_depth_score_threshold": 0.02,
            "exclude_alts": ["BAD_ALT", "ZZZ"],
        }
    }


def test_filter_by_alt_values_empty_df(kestrel_config_mock):
    """Test that an empty DataFrame simply returns empty without error."""
    df = pd.DataFrame()
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)
    assert out.empty, "Empty input should yield empty output."


def test_filter_by_alt_values_missing_columns(kestrel_config_mock):
    """Test that missing 'ALT' or 'Depth_Score' columns raises KeyError."""
    df = pd.DataFrame(
        {
            "ALT": ["GG", "ABC"],
            # 'Depth_Score' is missing here
        }
    )

    with pytest.raises(KeyError) as exc_info:
        filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    assert "Missing required columns" in str(
        exc_info.value
    ), "Expected KeyError due to missing required 'Depth_Score' column."


def test_filter_by_alt_values_gg_filter_below_threshold(kestrel_config_mock):
    """
    Test that GG variants below the depth score threshold are correctly flagged.

    When ALT='GG' and Depth_Score < threshold (0.02), the row is marked
    with alt_filter_pass=False but retained in the DataFrame.
    This follows the flag-and-defer pattern where filtering is done downstream.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "ALT": ["GG", "GG", "XYZ"],
            "Depth_Score": [
                0.019,  # Below threshold (0.02) for GG => should fail
                0.02,   # At threshold for GG => should pass
                0.5,    # Non-GG ALT => unaffected by GG filter
            ],
        }
    )

    # Act
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    # Assert - all rows are retained
    assert len(out) == 3, "All 3 rows should be retained in the DataFrame"

    # Assert - alt_filter_pass column exists
    assert "alt_filter_pass" in out.columns, "Should have 'alt_filter_pass' column"

    # Assert - row 0: GG with Depth_Score=0.019 should fail (below threshold)
    assert not out.loc[0, "alt_filter_pass"], \
        "Row 0 (GG with Depth_Score=0.019) should have alt_filter_pass=False (below threshold)"

    # Assert - row 1: GG with Depth_Score=0.02 should pass (at threshold)
    assert out.loc[1, "alt_filter_pass"], \
        "Row 1 (GG with Depth_Score=0.02) should have alt_filter_pass=True (meets threshold)"

    # Assert - row 2: XYZ should pass (not affected by GG-specific filter)
    assert out.loc[2, "alt_filter_pass"], \
        "Row 2 (XYZ) should have alt_filter_pass=True (non-GG ALT not affected)"

    # Assert - verify the count of passing rows
    pass_count = out["alt_filter_pass"].sum()
    assert pass_count == 2, "Should have exactly 2 rows with alt_filter_pass=True"


def test_filter_by_alt_values_exclude_alts(kestrel_config_mock):
    """
    Test that ALTs in 'exclude_alts' list are correctly flagged.

    ALTs in the exclude_alts list (e.g., 'BAD_ALT', 'ZZZ') are marked
    with alt_filter_pass=False but retained in the DataFrame.
    This follows the flag-and-defer pattern where filtering is done downstream.
    """
    # Arrange
    # Config has exclude_alts = ["BAD_ALT", "ZZZ"]
    df = pd.DataFrame(
        {
            "ALT": ["GG", "BAD_ALT", "OK_ALT", "ZZZ", "ANOTHER"],
            "Depth_Score": [0.5, 0.1, 0.3, 0.2, 0.6],  # just some placeholder scores
        }
    )

    # Act
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    # Assert - all rows are retained
    assert len(out) == 5, "All 5 rows should be retained in the DataFrame"

    # Assert - alt_filter_pass column exists
    assert "alt_filter_pass" in out.columns, "Should have 'alt_filter_pass' column"

    # Assert - check each row's filter status
    # Row 0: GG with good Depth_Score => should pass
    assert out.loc[0, "alt_filter_pass"], \
        "Row 0 (GG with Depth_Score=0.5) should have alt_filter_pass=True"

    # Row 1: BAD_ALT is in exclude_alts => should fail
    assert not out.loc[1, "alt_filter_pass"], \
        "Row 1 (BAD_ALT) should have alt_filter_pass=False (in exclude_alts)"

    # Row 2: OK_ALT is not in exclude_alts => should pass
    assert out.loc[2, "alt_filter_pass"], \
        "Row 2 (OK_ALT) should have alt_filter_pass=True (not in exclude_alts)"

    # Row 3: ZZZ is in exclude_alts => should fail
    assert not out.loc[3, "alt_filter_pass"], \
        "Row 3 (ZZZ) should have alt_filter_pass=False (in exclude_alts)"

    # Row 4: ANOTHER is not in exclude_alts => should pass
    assert out.loc[4, "alt_filter_pass"], \
        "Row 4 (ANOTHER) should have alt_filter_pass=True (not in exclude_alts)"

    # Assert - verify the count: 3 should pass (GG, OK_ALT, ANOTHER)
    pass_count = out["alt_filter_pass"].sum()
    assert (
        pass_count == 3
    ), "Should have exactly 3 rows with alt_filter_pass=True (GG, OK_ALT, ANOTHER)"

    # Assert - verify that excluded ALTs are still in the DataFrame
    all_alts = out["ALT"].tolist()
    assert "BAD_ALT" in all_alts, "'BAD_ALT' should still be present in DataFrame"
    assert "ZZZ" in all_alts, "'ZZZ' should still be present in DataFrame"


def test_filter_by_alt_values_drop_left_right(kestrel_config_mock):
    """
    Test that intermediate columns like 'left' and 'right' are retained for debugging.

    The refactored implementation keeps all columns including intermediate ones
    like 'left' and 'right' for debugging purposes. These columns may be dropped
    in a later finalization step downstream.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "ALT": ["GG", "ABC"],
            "Depth_Score": [0.05, 0.02],
            "left": ["some_left_data", "some_left_data"],
            "right": ["some_right_data", "some_right_data"],
        }
    )

    # Act
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    # Assert - both rows are retained
    assert len(out) == 2, "Expected 2 rows total."

    # Assert - alt_filter_pass column exists
    assert "alt_filter_pass" in out.columns, "Should have 'alt_filter_pass' column"

    # Assert - both rows pass the filter
    # Row 0: GG with Depth_Score=0.05 (above threshold) => should pass
    # Row 1: ABC (not in exclude_alts) => should pass
    assert out.loc[0, "alt_filter_pass"], "Row 0 (GG) should pass the filter"
    assert out.loc[1, "alt_filter_pass"], "Row 1 (ABC) should pass the filter"

    # Assert - 'left' and 'right' columns are RETAINED for debugging
    # This is a change from the old behavior where they were dropped
    assert (
        "left" in out.columns and "right" in out.columns
    ), "Intermediate columns 'left' and 'right' should be retained for debugging"

    # Assert - verify the data in left/right columns is unchanged
    assert all(
        out["left"] == "some_left_data"
    ), "'left' column data should be preserved"
    assert all(
        out["right"] == "some_right_data"
    ), "'right' column data should be preserved"
