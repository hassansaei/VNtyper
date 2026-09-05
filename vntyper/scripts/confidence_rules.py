"""Ordered confidence rule table for Kestrel call confidence assignment.

This module provides the pure decision layer for assigning confidence levels to
Kestrel VNTR calls based on depth scores and active-region k-mer depths.

The confidence assignment evaluates an ordered, first-match rule table with the
reporting floor as an outer precondition:
- Row 0: Depth_Score is NaN or below the floor -> Negative
- Row 1: Depth_Score in [low, high] (closed interval) -> Low_Precision
- Row 2: Depth_Score > high and alt depth >= alt_mid_high -> High_Precision*
- Row 3: Depth_Score > high and alt depth in [alt_mid_low, alt_mid_high) -> High_Precision
- Row 4: Depth_Score > high and alt depth <= alt_low -> Low_Precision
- Row 5: Depth_Score > high and alt depth in (alt_low, alt_mid_low) and region depth <= var_region_threshold -> Low_Precision
- Row 6: otherwise -> Negative

Precedence note:
Row 1 covers the closed interval [low, high] at every alternate depth, subsuming
the former `cond3` (alt depth in [21, 100) and Depth_Score in [low, high]).
As documented by @hassansaei on #184: "do not remove this intent" -- the mid-range
Depth_Score demotion applies across all alternate depths, so `cond3`'s intent is
fully satisfied by Row 1 covering the entire closed interval.

A source-order first-match table is strictly forbidden: evaluating the active-region
depth check (cond1) before the High tiers would restore the region-depth cap rejected
on #183.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Final

import pandas as pd

logger = logging.getLogger(__name__)

#: Canonical confidence label names.
NEGATIVE_LABEL: Final[str] = "Negative"
LOW_PRECISION_LABEL: Final[str] = "Low_Precision"
HIGH_PRECISION_LABEL: Final[str] = "High_Precision"
HIGH_PRECISION_STAR_LABEL: Final[str] = "High_Precision*"


def _get_threshold(thresholds: Mapping[str, Any], key: str, section: str | None = None) -> Any:
    """Extract a threshold value supporting both nested config dicts and flat dictionaries.

    Args:
        thresholds: Dictionary containing threshold configurations.
        key: Specific threshold key to extract.
        section: Optional nested section name (e.g. 'depth_score_thresholds').

    Returns:
        The threshold value.

    Raises:
        KeyError: If the threshold cannot be resolved from the configuration.
    """
    if (
        section is not None
        and section in thresholds
        and isinstance(thresholds[section], Mapping)
        and key in thresholds[section]
    ):
        return thresholds[section][key]
    if key in thresholds:
        return thresholds[key]
    raise KeyError(f"Threshold '{key}' (section: '{section}') not found in configuration")


@dataclass(frozen=True)
class ConfidenceRule:
    """An immutable rule in the confidence assignment table.

    Attributes:
        name: A descriptive identifier for the rule.
        predicate: Callable taking (frame, thresholds) and returning a boolean Series.
        label: The confidence label assigned when this rule matches first.
    """

    name: str
    predicate: Callable[[pd.DataFrame, Mapping[str, Any]], pd.Series]
    label: str


def _predicate_subthreshold_or_nan(df: pd.DataFrame, thresholds: Mapping[str, Any]) -> pd.Series:
    floor = _get_threshold(thresholds, "reporting_floor", "confidence_assignment")
    ds = df["Depth_Score"]
    return ds.isna() | (ds < floor)


def _predicate_midband(df: pd.DataFrame, thresholds: Mapping[str, Any]) -> pd.Series:
    low = _get_threshold(thresholds, "low", "depth_score_thresholds")
    high = _get_threshold(thresholds, "high", "depth_score_thresholds")
    ds = df["Depth_Score"]
    return ds.between(low, high, inclusive="both")


def _predicate_high_precision_star(df: pd.DataFrame, thresholds: Mapping[str, Any]) -> pd.Series:
    high = _get_threshold(thresholds, "high", "depth_score_thresholds")
    alt_mid_high = _get_threshold(thresholds, "mid_high", "alt_depth_thresholds")
    ds = df["Depth_Score"]
    ad = df["Estimated_Depth_AlternateVariant"]
    return (ds >= high) & (ad >= alt_mid_high)


def _predicate_high_precision(df: pd.DataFrame, thresholds: Mapping[str, Any]) -> pd.Series:
    high = _get_threshold(thresholds, "high", "depth_score_thresholds")
    alt_mid_low = _get_threshold(thresholds, "mid_low", "alt_depth_thresholds")
    alt_mid_high = _get_threshold(thresholds, "mid_high", "alt_depth_thresholds")
    ds = df["Depth_Score"]
    ad = df["Estimated_Depth_AlternateVariant"]
    return (ds > high) & (ad >= alt_mid_low) & (ad < alt_mid_high)


def _predicate_low_precision_alt_low(df: pd.DataFrame, thresholds: Mapping[str, Any]) -> pd.Series:
    high = _get_threshold(thresholds, "high", "depth_score_thresholds")
    alt_low = _get_threshold(thresholds, "low", "alt_depth_thresholds")
    ds = df["Depth_Score"]
    ad = df["Estimated_Depth_AlternateVariant"]
    return (ds > high) & (ad <= alt_low)


def _predicate_low_precision_region_gap(df: pd.DataFrame, thresholds: Mapping[str, Any]) -> pd.Series:
    high = _get_threshold(thresholds, "high", "depth_score_thresholds")
    alt_low = _get_threshold(thresholds, "low", "alt_depth_thresholds")
    alt_mid_low = _get_threshold(thresholds, "mid_low", "alt_depth_thresholds")
    var_region = _get_threshold(thresholds, "var_active_region_threshold")
    ds = df["Depth_Score"]
    ad = df["Estimated_Depth_AlternateVariant"]
    rd = df["Estimated_Depth_Variant_ActiveRegion"]
    return (ds > high) & (ad > alt_low) & (ad < alt_mid_low) & (rd <= var_region)


def _predicate_otherwise(df: pd.DataFrame, _thresholds: Mapping[str, Any]) -> pd.Series:
    return pd.Series(True, index=df.index)


#: The canonical ordered confidence rule table (rows 0 to 6).
CONFIDENCE_RULES: Final[tuple[ConfidenceRule, ...]] = (
    ConfidenceRule("subthreshold_or_nan", _predicate_subthreshold_or_nan, NEGATIVE_LABEL),
    ConfidenceRule("midband_demotion", _predicate_midband, LOW_PRECISION_LABEL),
    ConfidenceRule("high_precision_star", _predicate_high_precision_star, HIGH_PRECISION_STAR_LABEL),
    ConfidenceRule("high_precision", _predicate_high_precision, HIGH_PRECISION_LABEL),
    ConfidenceRule("low_precision_alt_low", _predicate_low_precision_alt_low, LOW_PRECISION_LABEL),
    ConfidenceRule("low_precision_region_gap", _predicate_low_precision_region_gap, LOW_PRECISION_LABEL),
    ConfidenceRule("fallback_negative", _predicate_otherwise, NEGATIVE_LABEL),
)


def assign_confidence_labels(
    frame: pd.DataFrame,
    thresholds: Mapping[str, Any],
    rules: tuple[ConfidenceRule, ...] = CONFIDENCE_RULES,
) -> pd.Series:
    """Assign confidence labels to candidate rows using ordered, first-match rules.

    Args:
        frame: DataFrame containing Depth_Score, Estimated_Depth_AlternateVariant,
            and Estimated_Depth_Variant_ActiveRegion columns.
        thresholds: Dictionary or mapping of threshold configurations.
        rules: Ordered tuple of ConfidenceRule instances (defaults to CONFIDENCE_RULES).

    Returns:
        pd.Series of string confidence labels aligned with frame.index.
    """
    if frame.empty:
        return pd.Series(dtype=object, name="Confidence")

    # Cosmetic label override mapping if provided in thresholds
    confidence_levels = thresholds.get("confidence_levels", {})
    if isinstance(confidence_levels, Mapping):
        label_map = {
            NEGATIVE_LABEL: NEGATIVE_LABEL,
            LOW_PRECISION_LABEL: confidence_levels.get("low_precision", LOW_PRECISION_LABEL),
            HIGH_PRECISION_LABEL: confidence_levels.get("high_precision", HIGH_PRECISION_LABEL),
            HIGH_PRECISION_STAR_LABEL: confidence_levels.get("high_precision_star", HIGH_PRECISION_STAR_LABEL),
        }
    else:
        label_map = {}

    labels = pd.Series(index=frame.index, dtype=object, name="Confidence")
    unassigned = pd.Series(True, index=frame.index)

    for rule in rules:
        if not unassigned.any():
            break
        matched_mask = rule.predicate(frame, thresholds) & unassigned
        resolved_label = label_map.get(rule.label, rule.label)
        labels.loc[matched_mask] = resolved_label
        unassigned = unassigned & ~matched_mask

    if unassigned.any():
        labels.loc[unassigned] = NEGATIVE_LABEL

    return labels
