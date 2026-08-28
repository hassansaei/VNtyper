"""Pure fastp cutoff configuration and presentation contracts."""

from __future__ import annotations

import importlib
import json
from decimal import Decimal
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

    assert cutoffs.duplication_rate.value == Decimal("0.1234")
    assert cutoffs.duplication_rate.label == "12.34%"
    assert cutoffs.q20_rate.value == Decimal("0.7555")
    assert cutoffs.q20_rate.label == "75.55%"
    assert cutoffs.q30_rate.value == Decimal("0.6543")
    assert cutoffs.q30_rate.label == "65.43%"
    assert cutoffs.passed_filter_rate.value == Decimal("0.7765")
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

    assert cutoffs.duplication_rate.value == Decimal(str(value))
    assert cutoffs.duplication_rate.label == label


def test_build_fastp_cutoffs_derives_the_decision_and_label_at_display_precision() -> None:
    """A fractional cutoff shares the two-decimal percentage domain of its label."""
    module = _module()

    cutoffs = module.build_fastp_cutoffs(
        {
            "duplication_rate": 0.1,
            "q20_rate": 0.77654,
            "q30_rate": 0.7,
            "passed_filter_reads_rate": 0.8,
        }
    )

    assert cutoffs.q20_rate.value == Decimal("0.7765")
    assert cutoffs.q20_rate.label == "77.65%"


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


def test_build_fastp_cutoffs_normalizes_an_oversized_json_integer(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A valid-JSON integer too large for float still follows the config error contract."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")

    with pytest.raises(ValueError, match="duplication_rate"):
        module.build_fastp_cutoffs(
            {
                "duplication_rate": 10**400,
                "q20_rate": 0.8,
                "q30_rate": 0.7,
                "passed_filter_reads_rate": 0.8,
            }
        )

    assert "duplication_rate" in caplog.text


@pytest.mark.parametrize("metric_key", ("duplication_rate", "q20_rate", "q30_rate", "passed_filter_rate"))
@pytest.mark.parametrize(
    "raw_rate",
    (True, "0.8", float("nan"), float("inf"), float("-inf"), -0.0001, 1.0001, 10**400),
)
def test_build_fastp_measurement_rejects_every_malformed_nonmissing_rate(
    metric_key: str, raw_rate: object, caplog: pytest.LogCaptureFixture
) -> None:
    """Every displayed fastp metric fails closed with its output key in the log."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")

    with pytest.raises(ValueError, match=metric_key):
        module.build_fastp_measurement(raw_rate, metric_key)

    assert metric_key in caplog.text


@pytest.mark.parametrize("token", ("NaN", "Infinity", "-Infinity"))
def test_build_fastp_measurement_rejects_nonfinite_values_json_loads_accepts(token: str) -> None:
    """Bare nonfinite JSON tokens reach the report parser and must fail explicitly."""
    module = _module()
    raw_rate = json.loads(f'{{"rate": {token}}}')["rate"]

    with pytest.raises(ValueError, match="q20_rate"):
        module.build_fastp_measurement(raw_rate, "q20_rate")


@pytest.mark.parametrize(
    ("metric_key", "raw_rate", "cutoff", "higher_better", "expected_fraction", "expected_display", "expected_colour"),
    (
        ("duplication_rate", Decimal("0.05045"), Decimal("0.0505"), False, Decimal("0.0505"), "5.05%", "green"),
        ("q20_rate", Decimal("0.77645"), Decimal("0.7765"), True, Decimal("0.7765"), "77.65%", "green"),
        ("q30_rate", Decimal("0.70045"), Decimal("0.7005"), True, Decimal("0.7005"), "70.05%", "green"),
        (
            "passed_filter_rate",
            Decimal("0.80045"),
            Decimal("0.8005"),
            True,
            Decimal("0.8005"),
            "80.05%",
            "green",
        ),
    ),
)
def test_build_fastp_measurement_supplies_one_validated_display_and_icon_value(
    metric_key: str,
    raw_rate: Decimal,
    cutoff: Decimal,
    higher_better: bool,
    expected_fraction: Decimal,
    expected_display: str,
    expected_colour: str,
) -> None:
    """The same validated Decimal value feeds each metric's visible text and icon."""
    module = _module()

    measurement = module.build_fastp_measurement(raw_rate, metric_key)

    assert measurement.value == expected_fraction
    assert measurement.display == expected_display
    assert threshold_icon(measurement.value, cutoff, higher_better=higher_better)[1] == expected_colour


@pytest.mark.parametrize("metric_key", ("duplication_rate", "q20_rate", "q30_rate", "passed_filter_rate"))
def test_build_fastp_measurement_preserves_none_as_missing(metric_key: str) -> None:
    """A missing fastp measurement stays an empty display/icon representation."""
    module = _module()

    measurement = module.build_fastp_measurement(None, metric_key)

    assert measurement.value is None
    assert measurement.display is None


@pytest.mark.parametrize(
    ("raw_rate", "expected"),
    [
        (None, None),
        (0.0, Decimal(0)),
        (1.0, Decimal(1)),
        (0.6005, Decimal("0.6005")),
        (0.77654, Decimal("0.7765")),
        (0.77655, Decimal("0.7766")),
    ],
)
def test_fastp_threshold_rate_matches_the_value_rendered_in_the_report(
    raw_rate: float | None, expected: Decimal | None
) -> None:
    """Icon decisions use the same two-decimal percentage precision readers see."""
    module = _module()

    if expected is None:
        assert module.fastp_threshold_rate(raw_rate) is None
    else:
        assert module.fastp_threshold_rate(raw_rate) == expected


@pytest.mark.parametrize(
    ("raw_rate", "expected"),
    [
        (None, None),
        (0.1, "10.0%"),
        (0.05045, "5.05%"),
        (0.77645, "77.65%"),
        (0.80045, "80.05%"),
    ],
)
def test_fastp_display_rate_uses_the_icon_decision_rounding(raw_rate: float | None, expected: str | None) -> None:
    """Visible metric text is formatted from the same decimal decision rate."""
    module = _module()
    formatter = getattr(module, "fastp_display_rate", None)

    assert callable(formatter)
    assert formatter(raw_rate) == expected


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

    assert displayed_rate == Decimal(str(expected_displayed_rate))
    assert threshold_icon(displayed_rate, Decimal(str(cutoff)), higher_better=higher_better)[1] == expected_colour
