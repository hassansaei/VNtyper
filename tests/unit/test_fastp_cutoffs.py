"""Pure fastp cutoff configuration and presentation contracts."""

from __future__ import annotations

import importlib
from typing import Any

import pytest

from vntyper.scripts.report_formatting import threshold_icon

pytestmark = pytest.mark.unit


def _module() -> Any:
    """Import the extracted fastp cutoff module so its initial RED is collectable."""
    return importlib.import_module("vntyper.scripts.fastp_cutoffs")


def test_build_fastp_cutoffs_pairs_each_numeric_decision_with_its_label() -> None:
    """One validated representation owns each fastp decision cutoff and display label."""
    module = _module()

    cutoffs = module.build_fastp_cutoffs(
        {
            "duplication_rate": 0.1234,
            "q20_rate": 0.7555,
            "q30_rate": 0.6543,
            "passed_filter_reads_rate": 0.7765,
        }
    )

    assert cutoffs.duplication_rate.value == pytest.approx(0.1234)
    assert cutoffs.duplication_rate.label == "12.34%"
    assert cutoffs.q20_rate.value == pytest.approx(0.7555)
    assert cutoffs.q20_rate.label == "75.55%"
    assert cutoffs.q30_rate.value == pytest.approx(0.6543)
    assert cutoffs.q30_rate.label == "65.43%"
    assert cutoffs.passed_filter_rate.value == pytest.approx(0.7765)
    assert cutoffs.passed_filter_rate.label == "77.65%"


@pytest.mark.parametrize(("value", "label"), [(0.0, "0%"), (1.0, "100%")])
def test_build_fastp_cutoffs_accepts_inclusive_boundaries_with_exact_labels(value: float, label: str) -> None:
    """Catch either inclusive endpoint becoming invalid or rendering imprecisely."""
    module = _module()
    cutoffs = module.build_fastp_cutoffs(
        {
            "duplication_rate": value,
            "q20_rate": 0.8,
            "q30_rate": 0.7,
            "passed_filter_reads_rate": 0.8,
        }
    )

    assert cutoffs.duplication_rate.value == value
    assert cutoffs.duplication_rate.label == label


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("duplication_rate", None),
        ("duplication_rate", -0.01),
        ("q20_rate", "0.8"),
        ("q30_rate", True),
        ("passed_filter_reads_rate", float("nan")),
        ("passed_filter_reads_rate", 1.01),
    ],
)
def test_build_fastp_cutoffs_rejects_each_invalid_configured_fraction(key: str, value: object) -> None:
    """Missing, nonnumeric, nonfinite, and out-of-range cutoffs fail with their key."""
    module = _module()
    thresholds: dict[str, object] = {
        "duplication_rate": 0.1,
        "q20_rate": 0.8,
        "q30_rate": 0.7,
        "passed_filter_reads_rate": 0.8,
    }
    if value is None:
        del thresholds[key]
    else:
        thresholds[key] = value

    with pytest.raises(ValueError, match=key):
        module.build_fastp_cutoffs(thresholds)


@pytest.mark.parametrize(
    ("raw_rate", "expected"),
    [(None, None), (0.0, 0.0), (0.77654, 0.7765), (0.77655, 0.7766)],
)
def test_fastp_threshold_rate_matches_the_value_rendered_in_the_report(
    raw_rate: float | None, expected: float | None
) -> None:
    """Icon decisions use the same two-decimal percentage precision readers see."""
    module = _module()

    if expected is None:
        assert module.fastp_threshold_rate(raw_rate) is None
    else:
        assert module.fastp_threshold_rate(raw_rate) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("raw_rate", "cutoff", "higher_better", "expected_displayed_rate", "expected_colour"),
    [
        (0.10004, 0.1, False, 0.1, "green"),
        (0.10006, 0.1, False, 0.1001, "red"),
        (0.79996, 0.8, True, 0.8, "green"),
        (0.79994, 0.8, True, 0.7999, "red"),
        (0.00065, 0.0007, True, 0.0007, "green"),
    ],
)
def test_fastp_threshold_rate_is_the_icon_decision_value(
    raw_rate: float,
    cutoff: float,
    higher_better: bool,
    expected_displayed_rate: float,
    expected_colour: str,
) -> None:
    """Icon decisions compare the report's displayed precision to the paired cutoff."""
    module = _module()
    displayed_rate = module.fastp_threshold_rate(raw_rate)

    assert displayed_rate == pytest.approx(expected_displayed_rate)
    assert threshold_icon(displayed_rate, cutoff, higher_better=higher_better)[1] == expected_colour
