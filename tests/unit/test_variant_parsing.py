#!/usr/bin/env python3
# tests/unit/test_variant_parsing.py

"""
Unit tests for variant_parsing.py, focusing on:
  filter_by_alt_values_and_finalize()

Ensures that ALT-based filtering rules are correctly applied:
  - 'GG' alt requires a minimum Depth_Score.
  - exclude_alts removed from final DataFrame.
  - left/right columns dropped at the end.
"""

import pytest
import pandas as pd
from vntyper.scripts.variant_parsing import filter_by_alt_values_and_finalize


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
            "exclude_alts": ["BAD_ALT", "ZZZ"]
        }
    }


def test_filter_by_alt_values_empty_df(kestrel_config_mock):
    """
    Test that an empty DataFrame simply returns empty without error.
    """
    df = pd.DataFrame()
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)
    assert out.empty, "Empty input should yield empty output."


def test_filter_by_alt_values_missing_columns(kestrel_config_mock):
    """
    Test that missing 'ALT' or 'Depth_Score' columns raises KeyError.
    """
    df = pd.DataFrame({
        "ALT": ["GG", "ABC"],
        # 'Depth_Score' is missing here
    })

    with pytest.raises(KeyError) as exc_info:
        filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    assert "Missing required columns" in str(exc_info.value), (
        "Expected KeyError due to missing required 'Depth_Score' column."
    )


def test_filter_by_alt_values_gg_filter_below_threshold(kestrel_config_mock):
    """
    If ALT='GG' but Depth_Score < threshold, that row should be removed.
    """
    df = pd.DataFrame({
        "ALT": ["GG", "GG", "XYZ"],
        "Depth_Score": [0.019, 0.02, 0.5]  # 0.019 < threshold => remove, 0.02 >= threshold => keep
    })

    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)
    # The first row has Depth_Score=0.019 => < 0.02 => removed
    # The second row has Depth_Score=0.02 => OK => keep
    # The third row has ALT=XYZ => unaffected by the GG filter => keep
    assert len(out) == 2, (
        "Expected only 2 rows to remain after removing GG with insufficient Depth_Score."
    )
    # Check that 'GG' row with Depth_Score=0.019 was removed
    assert (out["Depth_Score"] < 0.02).sum() == 0, "No row should have Depth_Score < 0.02 for 'GG' alt."


def test_filter_by_alt_values_exclude_alts(kestrel_config_mock):
    """
    Test that ALTs in 'exclude_alts' are removed from the DataFrame.
    """
    df = pd.DataFrame({
        "ALT": ["GG", "BAD_ALT", "OK_ALT", "ZZZ", "ANOTHER"],
        "Depth_Score": [0.5, 0.1, 0.3, 0.2, 0.6]  # just some placeholder scores
    })

    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)
    # Excluded: "BAD_ALT" and "ZZZ"
    # Keep: "GG", "OK_ALT", "ANOTHER"
    kept_alts = out["ALT"].tolist()
    assert len(out) == 3, "Expected 3 ALTs after excluding 'BAD_ALT' and 'ZZZ'."
    assert "BAD_ALT" not in kept_alts, "'BAD_ALT' should be removed."
    assert "ZZZ" not in kept_alts, "'ZZZ' should be removed."


def test_filter_by_alt_values_drop_left_right(kestrel_config_mock):
    """
    Test that 'left' and 'right' columns (if present) are dropped.
    """
    df = pd.DataFrame({
        "ALT": ["GG", "ABC"],
        "Depth_Score": [0.05, 0.02],
        "left": ["some_left_data", "some_left_data"],
        "right": ["some_right_data", "some_right_data"]
    })
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    # We keep both rows because 'GG' with Depth_Score=0.05 is fine,
    # and "ABC" is not in exclude_alts => all good.
    assert len(out) == 2, "Expected 2 rows total."
    assert "left" not in out.columns and "right" not in out.columns, (
        "Expected the 'left' and 'right' columns to be dropped."
    )
