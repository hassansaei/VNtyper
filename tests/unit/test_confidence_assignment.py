#!/usr/bin/env python3
# tests/unit/test_confidence_assignment.py

"""
Unit tests for confidence assignment functionality.
Validates confidence levels based on depth scores and thresholds
from vntyper/scripts/kestrel_config.json.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts.confidence_assignment import (
    NEGATIVE_LABEL,
    calculate_depth_score_and_assign_confidence,
)

pytestmark = pytest.mark.unit

# Fixed region depth used across tests to isolate Depth_Score effects.
# Must be above var_active_region_threshold (200) to avoid triggering cond1's region check.
_REGION_DEPTH = 10000


@pytest.fixture(scope="session")
def kestrel_config():
    """Load the Kestrel configuration from vntyper/scripts/kestrel_config.json."""
    this_file = Path(__file__).resolve()
    config_path = this_file.parents[2] / "vntyper" / "scripts" / "kestrel_config.json"
    if not config_path.exists():
        pytest.exit(f"kestrel_config.json not found at {config_path}", returncode=1)
    with config_path.open("r") as f:
        return json.load(f)


def _make_df(alt_depth: float, region_depth: float = _REGION_DEPTH) -> pd.DataFrame:
    """Helper to create a single-row DataFrame with the required depth columns."""
    return pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [alt_depth],
            "Estimated_Depth_Variant_ActiveRegion": [region_depth],
        }
    )


def test_calculate_depth_score_empty_df(kestrel_config):
    """Empty input should yield empty output."""
    df = pd.DataFrame()
    out = calculate_depth_score_and_assign_confidence(df, kestrel_config)
    assert out.empty


def test_depth_score_below_threshold_is_negative(kestrel_config):
    """Variants with Depth_Score < low_threshold get Negative (filtered out)."""
    low = kestrel_config["confidence_assignment"]["depth_score_thresholds"]["low"]
    # Derive alt_depth that produces Depth_Score just below low_threshold
    alt_depth = low * _REGION_DEPTH * 0.5  # half the threshold
    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == NEGATIVE_LABEL
    assert not out.loc[0, "depth_confidence_pass"]


def test_depth_score_at_threshold_boundary_not_negative(kestrel_config):
    """Variants with Depth_Score == low_threshold should NOT be Negative."""
    low = kestrel_config["confidence_assignment"]["depth_score_thresholds"]["low"]
    # Derive alt_depth that produces Depth_Score exactly at low_threshold
    alt_depth = low * _REGION_DEPTH
    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth), kestrel_config)
    assert out.loc[0, "Confidence"] != NEGATIVE_LABEL, (
        "Depth_Score at exactly low_threshold must not be Negative; "
        "the Negative filter applies only to strictly lower scores."
    )
    assert out.loc[0, "depth_confidence_pass"]


def test_depth_score_in_low_precision_band(kestrel_config):
    """Variants with Depth_Score in (low, high) band get Low_Precision."""
    conf = kestrel_config["confidence_assignment"]
    low = conf["depth_score_thresholds"]["low"]
    high = conf["depth_score_thresholds"]["high"]
    mid_low = conf["alt_depth_thresholds"]["mid_low"]
    mid_high = conf["alt_depth_thresholds"]["mid_high"]
    low_label = conf["confidence_levels"]["low_precision"]

    # Target Depth_Score in (low, high) and alt in [mid_low, mid_high] => cond3 Low_Precision
    depth_score = (low + high) / 2
    alt_depth = (mid_low + mid_high) // 2
    region_depth = int(alt_depth / depth_score) + 1

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == low_label


def test_depth_score_high_precision_star(kestrel_config):
    """High alt depth + high Depth_Score => High_Precision*."""
    conf = kestrel_config["confidence_assignment"]
    high = conf["depth_score_thresholds"]["high"]
    alt_mid_high = conf["alt_depth_thresholds"]["mid_high"]
    high_star_label = conf["confidence_levels"]["high_precision_star"]

    # alt >= mid_high AND Depth_Score >= high_threshold => cond2
    alt_depth = alt_mid_high
    region_depth = int(alt_depth / high)  # produces Depth_Score >= high

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == high_star_label


def test_low_region_depth_with_sufficient_depth_score(kestrel_config):
    """Low region depth (below var_active_region_threshold) should get Low_Precision, not Negative."""
    conf = kestrel_config["confidence_assignment"]
    low = conf["depth_score_thresholds"]["low"]
    var_region_threshold = conf["var_active_region_threshold"]
    low_label = conf["confidence_levels"]["low_precision"]

    # Region depth at threshold, but Depth_Score above low_threshold
    region_depth = var_region_threshold
    alt_depth = low * region_depth * 2  # Depth_Score = 2 * low_threshold (well above)

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == low_label


def test_low_region_depth_with_insufficient_depth_score(kestrel_config):
    """Low region depth AND Depth_Score below threshold should still be Negative."""
    conf = kestrel_config["confidence_assignment"]
    low = conf["depth_score_thresholds"]["low"]
    var_region_threshold = conf["var_active_region_threshold"]

    # Both below thresholds
    region_depth = var_region_threshold
    alt_depth = low * region_depth * 0.5  # Depth_Score = 0.5 * low_threshold

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == NEGATIVE_LABEL
    assert not out.loc[0, "depth_confidence_pass"]
