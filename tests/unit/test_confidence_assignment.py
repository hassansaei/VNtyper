#!/usr/bin/env python3
# tests/unit/test_confidence_assignment.py

"""
Unit tests for confidence assignment functionality.
Validates confidence levels based on depth scores and thresholds
from vntyper/scripts/kestrel_config.json.
"""

import json
import pytest
import pandas as pd
from pathlib import Path
from vntyper.scripts.confidence_assignment import (
    calculate_depth_score_and_assign_confidence
)


@pytest.fixture(scope="session")
def kestrel_config():
    """
    Loads the Kestrel configuration from vntyper/scripts/kestrel_config.json.
    Must contain a top-level "confidence_assignment" key.

    Example structure of kestrel_config.json:
    {
      "confidence_assignment": {
        "depth_score_thresholds": { ... },
        "alt_depth_thresholds": { ... },
        "var_active_region_threshold": 5000,
        "confidence_levels": { ... }
      },
      ...
    }
    """
    # Adjust path if needed; we go two levels up from tests/unit to
    # find vntyper/scripts/kestrel_config.json.
    # (Alternatively, use an absolute path or a more direct approach.)
    this_file = Path(__file__).resolve()
    scripts_dir = this_file.parents[2] / "vntyper" / "scripts"  # go up 2 dirs
    config_path = scripts_dir / "kestrel_config.json"

    if not config_path.exists():
        pytest.exit(
            f"kestrel_config.json not found at {config_path}",
            returncode=1
        )

    with config_path.open("r") as f:
        raw_config = json.load(f)

    # Return the entire dictionary. The confidence_assignment subdict
    # is used by calculate_depth_score_and_assign_confidence internally.
    # i.e., that function does kestrel_config["confidence_assignment"].
    return raw_config


def test_calculate_depth_score_empty_df(kestrel_config):
    """
    Test behavior when input DataFrame is empty.
    """
    df = pd.DataFrame()
    out = calculate_depth_score_and_assign_confidence(df, kestrel_config)
    assert out.empty, "Empty input should yield empty output."


def test_calculate_depth_score_low_precision(kestrel_config):
    """
    Test that low depth scores result in 'Low_Precision' confidence.
    """
    df = pd.DataFrame({
        "Estimated_Depth_AlternateVariant": [10],
        "Estimated_Depth_Variant_ActiveRegion": [10000]
    })
    out = calculate_depth_score_and_assign_confidence(df, kestrel_config)

    # If your code references:
    #   conf_assign = kestrel_config["confidence_assignment"]
    # then inside conf_assign["confidence_levels"]["low_precision"] ...
    low_label = kestrel_config["confidence_assignment"]["confidence_levels"]["low_precision"]

    assert out.loc[0, "Confidence"] == low_label, (
        "Expected Low_Precision for low depth score."
    )


def test_calculate_depth_score_high_precision(kestrel_config):
    df = pd.DataFrame({
        "Estimated_Depth_AlternateVariant": [100],
        "Estimated_Depth_Variant_ActiveRegion": [5000]
    })
    out = calculate_depth_score_and_assign_confidence(df, kestrel_config)

    # If your function returns "High_Precision*",
    # then let's just confirm the actual returned value is "High_Precision*".
    high_star_label = kestrel_config["confidence_assignment"]["confidence_levels"]["high_precision_star"]
    assert out.loc[0, "Confidence"] == high_star_label, (
        "Expected 'High_Precision*' for depth_score=0.02 with alt=100 >= mid_high=100"
    )
