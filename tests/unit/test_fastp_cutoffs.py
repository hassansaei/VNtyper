"""Pure fastp cutoff configuration and presentation contracts."""

from __future__ import annotations

import importlib
import json
from decimal import Decimal
from typing import Any

import pytest

from vntyper.scripts.report_formatting import threshold_icon

pytestmark = pytest.mark.unit


class _UnstringifiableFloat(float):
    """A real number whose decimal conversion fails after type validation."""

    def __str__(self) -> str:
        raise ValueError("cannot stringify test value")


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


def test_build_fastp_cutoffs_logs_and_keys_a_decimal_conversion_failure(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The conversion handler, not only its preceding type guard, is executable."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")
    expected = "Config thresholds has invalid fastp cutoff 'duplication_rate': expected a finite fraction from 0 to 1."

    with pytest.raises(ValueError, match="duplication_rate") as exc_info:
        module.build_fastp_cutoffs(
            {
                "duplication_rate": _UnstringifiableFloat(0.1),
                "q20_rate": 0.8,
                "q30_rate": 0.7,
                "passed_filter_reads_rate": 0.8,
            }
        )

    assert str(exc_info.value) == expected
    assert isinstance(exc_info.value.__cause__, ValueError)
    assert str(exc_info.value.__cause__) == "cannot stringify test value"
    assert [record.getMessage() for record in caplog.records if record.name == "vntyper.scripts.fastp_cutoffs"] == [
        expected
    ]


@pytest.mark.parametrize(
    ("component", "value"),
    (
        ("passed_filter_reads", "not-a-count"),
        ("total_reads", "not-a-count"),
        ("passed_filter_reads", "80"),
        ("total_reads", "100"),
        ("passed_filter_reads", True),
        ("total_reads", True),
        ("passed_filter_reads", float("nan")),
        ("total_reads", float("nan")),
        ("passed_filter_reads", float("inf")),
        ("total_reads", float("inf")),
        ("passed_filter_reads", float("-inf")),
        ("total_reads", float("-inf")),
        ("passed_filter_reads", -1),
        ("total_reads", -1),
        ("passed_filter_reads", 80.5),
        ("total_reads", Decimal("100.5")),
    ),
)
def test_calculate_passed_filter_rate_rejects_each_invalid_source_count(
    component: str, value: object, caplog: pytest.LogCaptureFixture
) -> None:
    """Parsed source counts fail closed before arithmetic can lose their invalidity."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")
    counts: dict[str, object] = {"passed_filter_reads": 80, "total_reads": 100}
    counts[component] = value

    with pytest.raises(ValueError, match="passed_filter_rate"):
        module.calculate_passed_filter_rate(**counts)

    assert component in caplog.text


def test_calculate_passed_filter_rate_returns_the_valid_ratio() -> None:
    """A complete valid pair still produces its raw ratio before display rounding."""
    module = _module()

    rate = module.calculate_passed_filter_rate(passed_filter_reads=80, total_reads=100)

    assert rate == Decimal("0.8")
    assert type(rate) is Decimal


def test_calculate_passed_filter_rate_preserves_a_large_count_ratio_for_display_and_decision() -> None:
    """Large integer counts cannot cross a visible half-tie through binary float."""
    module = _module()

    rate = module.calculate_passed_filter_rate(
        passed_filter_reads=77644999999999999,
        total_reads=100000000000000000,
    )
    measurement = module.build_fastp_measurement(rate, "passed_filter_rate")

    assert rate == Decimal("0.77644999999999999")
    assert measurement.display == "77.64%"
    assert measurement.display != "77.65%"
    assert threshold_icon(measurement.value, Decimal("0.7765"), higher_better=True)[1] == "red"


@pytest.mark.parametrize(
    ("offset", "comparison", "displayed"),
    ((-1, "below", "60.04%"), (0, "equal", "60.05%"), (1, "above", "60.05%")),
)
def test_calculate_passed_filter_rate_does_not_preround_long_ratios_across_a_half_tie(
    offset: int, comparison: str, displayed: str
) -> None:
    """Count arithmetic must retain enough precision to decide the visible boundary."""
    module = _module()
    scale = 10**35

    rate = module.calculate_passed_filter_rate(
        passed_filter_reads=60045 * scale + offset,
        total_reads=100000 * scale,
    )
    measurement = module.build_fastp_measurement(rate, "passed_filter_rate")

    assert rate is not None
    if comparison == "below":
        assert rate < Decimal("0.60045")
    elif comparison == "equal":
        assert rate == Decimal("0.60045")
    else:
        assert rate > Decimal("0.60045")
    assert measurement.display == displayed


def test_calculate_passed_filter_rate_logs_and_keys_a_decimal_conversion_failure(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A source count that fails Decimal conversion follows the keyed error contract."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")
    expected = (
        "Fastp output has invalid passed_filter_rate source count 'passed_filter_reads': "
        "expected a finite non-negative integer count."
    )

    with pytest.raises(ValueError, match="passed_filter_reads") as exc_info:
        module.calculate_passed_filter_rate(
            passed_filter_reads=_UnstringifiableFloat(80),
            total_reads=100,
        )

    assert str(exc_info.value) == expected
    assert isinstance(exc_info.value.__cause__, ValueError)
    assert str(exc_info.value.__cause__) == "cannot stringify test value"
    assert [record.getMessage() for record in caplog.records if record.name == "vntyper.scripts.fastp_cutoffs"] == [
        expected
    ]


@pytest.mark.parametrize(
    ("missing_path", "before_filtering", "filtering_result"),
    (
        ("summary.before_filtering.total_reads", {}, {"passed_filter_reads": 80}),
        ("filtering_result.passed_filter_reads", {"total_reads": 100}, {}),
    ),
)
def test_calculate_passed_filter_rate_from_sources_rejects_each_missing_key(
    missing_path: str,
    before_filtering: dict[str, object],
    filtering_result: dict[str, object],
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Incomplete output cannot inherit plausible source-count defaults."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")
    expected = f"Fastp output is missing required passed_filter_rate source key {missing_path!r}."

    with pytest.raises(ValueError, match=missing_path.replace(".", r"\.")) as exc_info:
        module.calculate_passed_filter_rate_from_sources(before_filtering, filtering_result)

    assert str(exc_info.value) == expected
    assert isinstance(exc_info.value.__cause__, KeyError)
    assert [record.getMessage() for record in caplog.records if record.name == "vntyper.scripts.fastp_cutoffs"] == [
        expected
    ]


def test_calculate_passed_filter_rate_from_sources_extracts_both_required_keys() -> None:
    """The source helper reads the exact two paths used by shipped fastp output."""
    module = _module()

    rate = module.calculate_passed_filter_rate_from_sources(
        {"total_reads": 100},
        {"passed_filter_reads": 80},
    )

    assert rate == Decimal("0.8")
    assert type(rate) is Decimal


def test_calculate_passed_filter_rate_rejects_a_count_above_the_total(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A syntactically numeric count cannot imply an impossible rate above one."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")

    with pytest.raises(ValueError, match="passed_filter_rate"):
        module.calculate_passed_filter_rate(passed_filter_reads=101, total_reads=100)

    assert "passed_filter_reads" in caplog.text


def test_calculate_passed_filter_rate_preserves_zero_counts_as_missing() -> None:
    """An empty fastp measurement remains the supported missing-rate representation."""
    module = _module()

    assert module.calculate_passed_filter_rate(passed_filter_reads=0, total_reads=0) is None


@pytest.mark.parametrize("passed_filter_reads", (1, 80))
def test_calculate_passed_filter_rate_rejects_a_positive_count_with_zero_total(
    passed_filter_reads: int, caplog: pytest.LogCaptureFixture
) -> None:
    """A positive passed count cannot be missing when fastp says its total is zero."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")

    with pytest.raises(ValueError, match="passed_filter_rate"):
        module.calculate_passed_filter_rate(passed_filter_reads=passed_filter_reads, total_reads=0)

    assert "passed_filter_rate" in caplog.text
    assert "passed_filter_reads" in caplog.text


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


def test_build_fastp_measurement_logs_and_keys_a_decimal_conversion_failure(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A measured rate conversion error is wrapped and logged with its metric key."""
    module = _module()
    caplog.set_level("ERROR", logger="vntyper.scripts.fastp_cutoffs")
    expected = "Fastp output has invalid measured rate 'q20_rate': expected a finite numeric fraction from 0 to 1."

    with pytest.raises(ValueError, match="q20_rate") as exc_info:
        module.build_fastp_measurement(_UnstringifiableFloat(0.8), "q20_rate")

    assert str(exc_info.value) == expected
    assert isinstance(exc_info.value.__cause__, ValueError)
    assert str(exc_info.value.__cause__) == "cannot stringify test value"
    assert [record.getMessage() for record in caplog.records if record.name == "vntyper.scripts.fastp_cutoffs"] == [
        expected
    ]


def test_validated_fastp_fraction_preserves_an_exact_decimal_type_and_value() -> None:
    """The internal report summary cannot coerce an exact JSON rate back to float."""
    module = _module()

    fraction = module.validated_fastp_fraction(Decimal("0.60044999999999999"), "q20_rate")

    assert fraction == Decimal("0.60044999999999999")
    assert type(fraction) is Decimal


@pytest.mark.parametrize(
    ("rate", "displayed"),
    (
        ("0.600449999999999999999999999999999", "60.04%"),
        ("0.600450000000000000000000000000000", "60.05%"),
        ("0.600450000000000000000000000000001", "60.05%"),
    ),
)
def test_build_fastp_measurement_quantizes_long_decimals_without_ambient_prerounding(rate: str, displayed: str) -> None:
    """Digits beyond the default Decimal context still decide a half-up boundary."""
    module = _module()

    measurement = module.build_fastp_measurement(Decimal(rate), "q20_rate")

    assert measurement.display == displayed


def test_build_fastp_cutoffs_normalizes_signed_zero() -> None:
    """A valid negative-zero cutoff cannot leak its sign into the visible label."""
    module = _module()

    cutoff = module.build_fastp_cutoffs(
        {
            "duplication_rate": -0.0,
            "q20_rate": 0.8,
            "q30_rate": 0.7,
            "passed_filter_reads_rate": 0.8,
        }
    ).duplication_rate

    assert cutoff.value == Decimal(0)
    assert cutoff.value.as_tuple().sign == 0
    assert cutoff.label == "0%"


def test_build_fastp_measurement_normalizes_signed_zero() -> None:
    """A valid negative-zero measurement cannot render a negative percentage."""
    module = _module()

    measurement = module.build_fastp_measurement(-0.0, "duplication_rate")

    assert measurement.value == Decimal(0)
    assert measurement.value.as_tuple().sign == 0
    assert measurement.display == "0.0%"


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
def test_build_fastp_measurement_value_matches_the_value_rendered_in_the_report(
    raw_rate: float | None, expected: Decimal | None
) -> None:
    """Icon decisions use the same two-decimal percentage precision readers see."""
    module = _module()
    measurement = module.build_fastp_measurement(raw_rate, "q20_rate")

    assert measurement.value == expected


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
def test_build_fastp_measurement_display_uses_the_icon_decision_rounding(
    raw_rate: float | None, expected: str | None
) -> None:
    """Visible metric text is formatted from the same decimal decision rate."""
    module = _module()

    assert module.build_fastp_measurement(raw_rate, "q20_rate").display == expected


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
def test_build_fastp_measurement_is_the_icon_decision_value(
    raw_rate: float,
    cutoff: float,
    higher_better: bool,
    expected_displayed_rate: float,
    expected_colour: str,
) -> None:
    """Icon decisions compare the report's displayed precision to the paired cutoff."""
    module = _module()
    displayed_rate = module.build_fastp_measurement(raw_rate, "q20_rate").value

    assert displayed_rate == Decimal(str(expected_displayed_rate))
    assert threshold_icon(displayed_rate, Decimal(str(cutoff)), higher_better=higher_better)[1] == expected_colour


def test_fastp_cutoffs_does_not_expose_dead_split_value_or_display_wrappers() -> None:
    """Callers use the combined measurement so its visible value and decision cannot drift."""
    module = _module()

    assert not hasattr(module, "fastp_display_rate")
    assert not hasattr(module, "fastp_threshold_rate")
