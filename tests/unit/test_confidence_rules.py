"""Unit tests for the ordered confidence rule table (vntyper/scripts/confidence_rules.py).

Pins the pure decision logic:
- Row 0: Depth_Score is NaN or below the floor -> Negative
- Row 1: low <= Depth_Score <= high -> Low_Precision (closed interval)
- Row 2: Depth_Score > high and AD >= alt_mid_high -> High_Precision*
- Row 3: Depth_Score > high and alt_mid_low <= AD < alt_mid_high -> High_Precision
- Row 4: Depth_Score > high and AD <= alt_low -> Low_Precision
- Row 5: Depth_Score > high and alt_low < AD < alt_mid_low and RD <= var_region_threshold -> Low_Precision
- Row 6: otherwise -> Negative
- Outer precondition: sub-floor variants stay Negative
- Rule order: reversing adjacent overlapping rules alters classification
- Empty frame handling and threshold validation
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np
import pandas as pd
import pytest

from vntyper.scripts.confidence_rules import (
    CONFIDENCE_RULES,
    HIGH_PRECISION_LABEL,
    HIGH_PRECISION_STAR_LABEL,
    LOW_PRECISION_LABEL,
    NEGATIVE_LABEL,
    ConfidenceRule,
    assign_confidence_labels,
)

pytestmark = pytest.mark.unit

# Shipped threshold constants used for tests
_SHIPPED_THRESHOLDS: dict[str, Any] = {
    "reporting_floor": 0.00469,
    "depth_score_thresholds": {
        "low": 0.00469,
        "high": 0.00515,
    },
    "alt_depth_thresholds": {
        "low": 20,
        "mid_low": 21,
        "mid_high": 100,
    },
    "var_active_region_threshold": 200,
    "confidence_levels": {
        "low_precision": "Low_Precision",
        "high_precision": "High_Precision",
        "high_precision_star": "High_Precision*",
    },
}


def _make_frame(
    depth_score: float | None,
    alt_depth: float,
    region_depth: float,
) -> pd.DataFrame:
    """Construct a single-row DataFrame with the required depth and score columns."""
    ds_val = np.nan if depth_score is None else float(depth_score)
    return pd.DataFrame(
        {
            "Depth_Score": [ds_val],
            "Estimated_Depth_AlternateVariant": [float(alt_depth)],
            "Estimated_Depth_Variant_ActiveRegion": [float(region_depth)],
        }
    )


def test_confidence_rule_dataclass_structure() -> None:
    """ConfidenceRule must expose name, predicate, and label attributes and be frozen."""
    rule = ConfidenceRule(
        name="test_rule",
        predicate=lambda df, th: pd.Series(True, index=df.index),
        label="Test_Label",
    )
    assert rule.name == "test_rule"
    assert rule.label == "Test_Label"
    df = _make_frame(0.01, 50, 5000)
    assert bool(rule.predicate(df, _SHIPPED_THRESHOLDS).iloc[0]) is True
    with pytest.raises(dataclasses.FrozenInstanceError):
        rule.name = "other"  # type: ignore[misc]


def test_assign_confidence_labels_with_flat_thresholds() -> None:
    """Flat threshold configurations (without nested sections) must be supported."""
    flat_config = {
        "reporting_floor": 0.00469,
        "low": 0.00469,
        "high": 0.00515,
        "mid_low": 21,
        "mid_high": 100,
        "var_active_region_threshold": 200,
    }
    df = _make_frame(0.010, 50, 5000)
    assert assign_confidence_labels(df, flat_config).iloc[0] == HIGH_PRECISION_LABEL


def test_assign_confidence_labels_empty_frame() -> None:
    """An empty DataFrame must return an empty Series before accessing thresholds."""
    out = assign_confidence_labels(pd.DataFrame(), {})
    assert isinstance(out, pd.Series)
    assert out.empty
    assert out.name == "Confidence"


def test_row_0_subthreshold_or_nan_is_negative() -> None:
    """Row 0: Depth_Score NaN or below floor must be Negative."""
    df_nan = _make_frame(None, 50, 0)
    assert assign_confidence_labels(df_nan, _SHIPPED_THRESHOLDS).iloc[0] == NEGATIVE_LABEL

    df_sub = _make_frame(0.0040, 50, 12500)
    assert assign_confidence_labels(df_sub, _SHIPPED_THRESHOLDS).iloc[0] == NEGATIVE_LABEL


def test_row_1_midband_closed_interval_is_low_precision() -> None:
    """Row 1: Depth_Score in [low, high] inclusive must be Low_Precision across depths."""
    # Bottom endpoint: exactly low
    df_low = _make_frame(0.00469, 150, 31982)
    assert assign_confidence_labels(df_low, _SHIPPED_THRESHOLDS).iloc[0] == LOW_PRECISION_LABEL

    # Interior of band
    df_mid = _make_frame(0.0050, 150, 30000)
    assert assign_confidence_labels(df_mid, _SHIPPED_THRESHOLDS).iloc[0] == LOW_PRECISION_LABEL

    # Top endpoint: exactly high (#184 demotion)
    df_high = _make_frame(0.00515, 150, 29126)
    assert assign_confidence_labels(df_high, _SHIPPED_THRESHOLDS).iloc[0] == LOW_PRECISION_LABEL


def test_row_2_high_precision_star() -> None:
    """Row 2: Depth_Score > high and alt_depth >= alt_mid_high (100) -> High_Precision*."""
    df = _make_frame(0.010, 120, 12000)
    assert assign_confidence_labels(df, _SHIPPED_THRESHOLDS).iloc[0] == HIGH_PRECISION_STAR_LABEL


def test_row_3_high_precision() -> None:
    """Row 3: Depth_Score > high and alt_mid_low (21) <= alt_depth < alt_mid_high (100) -> High_Precision."""
    # Lower edge of alt depth (21)
    df_21 = _make_frame(0.010, 21, 2100)
    assert assign_confidence_labels(df_21, _SHIPPED_THRESHOLDS).iloc[0] == HIGH_PRECISION_LABEL

    # Middle of alt depth (50)
    df_50 = _make_frame(0.010, 50, 5000)
    assert assign_confidence_labels(df_50, _SHIPPED_THRESHOLDS).iloc[0] == HIGH_PRECISION_LABEL

    # Just below mid_high (99)
    df_99 = _make_frame(0.010, 99, 9900)
    assert assign_confidence_labels(df_99, _SHIPPED_THRESHOLDS).iloc[0] == HIGH_PRECISION_LABEL


def test_row_4_low_precision_alt_low() -> None:
    """Row 4: Depth_Score > high and alt_depth <= alt_low (20) -> Low_Precision."""
    df_10 = _make_frame(0.010, 10, 1000)
    assert assign_confidence_labels(df_10, _SHIPPED_THRESHOLDS).iloc[0] == LOW_PRECISION_LABEL

    df_20 = _make_frame(0.010, 20, 2000)
    assert assign_confidence_labels(df_20, _SHIPPED_THRESHOLDS).iloc[0] == LOW_PRECISION_LABEL


def test_row_5_low_precision_region_gap() -> None:
    """Row 5: Depth_Score > high, alt_depth in (20, 21) and region depth <= 200 -> Low_Precision."""
    df = _make_frame(0.1025, 20.5, 200)
    assert assign_confidence_labels(df, _SHIPPED_THRESHOLDS).iloc[0] == LOW_PRECISION_LABEL


def test_row_6_otherwise_negative() -> None:
    """Row 6: Depth_Score > high, alt_depth in (20, 21) and region depth > 200 -> Negative."""
    df = _make_frame(0.1020, 20.5, 201)
    assert assign_confidence_labels(df, _SHIPPED_THRESHOLDS).iloc[0] == NEGATIVE_LABEL


def test_outer_precondition_blocks_subfloor_promotion() -> None:
    """Variants with Depth_Score below low threshold must remain Negative even if alt <= 20 or region <= 200."""
    # alt_depth <= 20 but Depth_Score < 0.00469
    df_alt = _make_frame(0.0010, 10, 10000)
    assert assign_confidence_labels(df_alt, _SHIPPED_THRESHOLDS).iloc[0] == NEGATIVE_LABEL

    # region_depth <= 200 but Depth_Score < 0.00469
    df_reg = _make_frame(0.0010, 0.1, 100)
    assert assign_confidence_labels(df_reg, _SHIPPED_THRESHOLDS).iloc[0] == NEGATIVE_LABEL


def test_source_order_reversal_pins_high_precision_over_region_demotion() -> None:
    """Spec Section 5 & Issue #183: High tiers must be evaluated before region demotion.

    A variant with alt=50 on a 150-read region has Depth_Score = 50/150 = 0.333 > high.
    Even though region depth <= 200, it must be High_Precision, NOT Low_Precision.
    """
    df = _make_frame(50 / 150, 50, 150)
    assert assign_confidence_labels(df, _SHIPPED_THRESHOLDS).iloc[0] == HIGH_PRECISION_LABEL


def test_rule_table_length_and_names() -> None:
    """CONFIDENCE_RULES must contain exactly 7 rules in the verified order."""
    expected_names = (
        "subthreshold_or_nan",
        "midband_demotion",
        "high_precision_star",
        "high_precision",
        "low_precision_alt_low",
        "low_precision_region_gap",
        "fallback_negative",
    )
    assert len(CONFIDENCE_RULES) == 7
    assert tuple(r.name for r in CONFIDENCE_RULES) == expected_names


def test_reversing_adjacent_rules_changes_labels() -> None:
    """Reversing adjacent overlapping rule pairs alters classification outcomes.

    Proves that rule order is strictly load-bearing:
    - Pair (1, 2) [midband vs high*]: at Depth_Score == high (0.00515) and AD=120,
      Row 1 gives Low_Precision (#184 demotion). If Row 2 ran first, it would give High_Precision*.
    - Pair (5, 6) [gap region <=200 vs fallback]: for AD=20.5 and RD=200, Row 5 gives Low_Precision.
      If Row 6 (fallback) ran first, it would give Negative.
    """
    rules_list = list(CONFIDENCE_RULES)

    # Test swap (1, 2)
    swapped_1_2 = tuple(rules_list[0:1] + [rules_list[2], rules_list[1]] + rules_list[3:])
    df_boundary = _make_frame(0.00515, 120, 120 / 0.00515)
    # With normal rules: Low_Precision
    assert (
        assign_confidence_labels(df_boundary, _SHIPPED_THRESHOLDS, rules=CONFIDENCE_RULES).iloc[0]
        == LOW_PRECISION_LABEL
    )
    # With swapped rules (2 before 1): promoted to High_Precision*
    assert (
        assign_confidence_labels(df_boundary, _SHIPPED_THRESHOLDS, rules=swapped_1_2).iloc[0]
        == HIGH_PRECISION_STAR_LABEL
    )

    # Test swap (5, 6)
    swapped_5_6 = tuple(rules_list[0:5] + [rules_list[6], rules_list[5]])
    df_gap = _make_frame(0.1025, 20.5, 200)
    # With normal rules: Low_Precision
    assert assign_confidence_labels(df_gap, _SHIPPED_THRESHOLDS, rules=CONFIDENCE_RULES).iloc[0] == LOW_PRECISION_LABEL
    # With swapped rules (6 before 5): demoted to Negative
    assert assign_confidence_labels(df_gap, _SHIPPED_THRESHOLDS, rules=swapped_5_6).iloc[0] == NEGATIVE_LABEL


def test_missing_required_threshold_raises_keyerror() -> None:
    """Missing required calibration key must raise KeyError."""
    incomplete_config: dict[str, Any] = {
        "reporting_floor": 0.00469,
        "depth_score_thresholds": {"low": 0.00469},
        "alt_depth_thresholds": {"low": 20, "mid_low": 21, "mid_high": 100},
        "var_active_region_threshold": 200,
    }
    df = _make_frame(0.010, 50, 5000)
    with pytest.raises(KeyError, match="high"):
        assign_confidence_labels(df, incomplete_config)


def test_missing_reporting_floor_raises_keyerror() -> None:
    """Missing reporting_floor calibration key must raise KeyError."""
    incomplete_config = dict(_SHIPPED_THRESHOLDS)
    del incomplete_config["reporting_floor"]
    df = _make_frame(0.010, 50, 5000)
    with pytest.raises(KeyError, match="reporting_floor"):
        assign_confidence_labels(df, incomplete_config)


def test_cosmetic_confidence_levels_override() -> None:
    """Custom cosmetic labels in confidence_levels must be respected."""
    custom_config = dict(_SHIPPED_THRESHOLDS)
    custom_config["confidence_levels"] = {
        "low_precision": "Custom_Low",
        "high_precision": "Custom_High",
        "high_precision_star": "Custom_High*",
    }
    df_high = _make_frame(0.010, 50, 5000)
    assert assign_confidence_labels(df_high, custom_config).iloc[0] == "Custom_High"

    df_mid = _make_frame(0.0050, 50, 10000)
    assert assign_confidence_labels(df_mid, custom_config).iloc[0] == "Custom_Low"


def test_evaluation_breaks_early_when_all_rows_assigned() -> None:
    """Predicate evaluation must terminate as soon as all rows have been assigned."""

    class CountingTuple(tuple):  # type: ignore[type-arg]
        def __iter__(self):
            self.yielded = getattr(self, "yielded", 0)
            for item in super().__iter__():
                self.yielded += 1
                yield item

    counting_rules = CountingTuple(CONFIDENCE_RULES)
    counting_rules.yielded = 0
    df = _make_frame(None, 50, 0)
    assign_confidence_labels(df, _SHIPPED_THRESHOLDS, rules=counting_rules)
    assert counting_rules.yielded == 2


def test_non_mapping_confidence_levels_falls_back_to_empty_label_map() -> None:
    """Non-mapping confidence_levels in thresholds should use canonical rule labels."""
    cfg = dict(_SHIPPED_THRESHOLDS)
    cfg["confidence_levels"] = "invalid_not_a_mapping"
    df = _make_frame(0.010, 50, 5000)
    assert assign_confidence_labels(df, cfg).iloc[0] == HIGH_PRECISION_LABEL


def test_unassigned_rows_fallback_to_negative() -> None:
    """If custom rules leave unassigned rows, they must fall back to Negative."""
    empty_rules: tuple[ConfidenceRule, ...] = ()
    df = _make_frame(0.010, 50, 5000)
    assert assign_confidence_labels(df, _SHIPPED_THRESHOLDS, rules=empty_rules).iloc[0] == NEGATIVE_LABEL
