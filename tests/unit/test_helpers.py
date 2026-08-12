"""Tests for the shared assertion helpers in ``tests/helpers.py``.

``tests/helpers.py`` is the assertion oracle for both slow tiers -- the pytest
integration tests and the Docker integration tests -- so a validator of its own
that cannot fail makes "no regression" unsubstantiable. Two of them are pinned
here for that reason:

* ``validate_coverage_output`` must read the column names production actually
  writes -- ``mean``, ``median``, ``percent_uncovered``. Reading any other name
  with a ``0`` default returns ``0.0`` for every metric, and the only caller
  asserts ``mean_cov >= 0``.
* ``assert_pattern_match`` must compare a trailing ``*`` literally, not as a
  glob. It is applied to ``Confidence``, and ``High_Precision*`` is a *literal*
  confidence label in this codebase, so glob semantics would let plain
  ``High_Precision`` satisfy an expectation of ``High_Precision*`` and erase the
  boundary the assertion exists to check.

Both are asserted from the failing side as well as the passing one, so a
validator that stops discriminating fails these tests.
"""

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

from tests.helpers import (  # noqa: E402
    COVERAGE_COLUMNS,
    assert_pattern_match,
    validate_coverage_output,
    validate_kestrel_output,
)

# The four literal confidence labels this codebase can emit. Three come from
# kestrel_config.json's `confidence_levels`; "Negative" is the placeholder row
# kestrel_genotyping.py writes when nothing was called. Not one is a pattern.
CONFIDENCE_LABELS = ("High_Precision*", "High_Precision", "Low_Precision", "Negative")


def _write_tsv(path: Path, header: list[str], row: list[str]) -> None:
    """Write a one-row TSV.

    Args:
        path: File to create; parents are created.
        header: Column names.
        row: The single data row.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")


def _production_coverage_header() -> list[str]:
    """Read the coverage header straight out of the production writer.

    Rendered by calling the production formatter rather than by regexing its
    source. The header used to be an inline literal inside
    ``fastq_bam_processing.calculate_vntr_coverage``, and this function used to
    scrape it out with a regex -- which broke the moment the pure half moved to
    ``coverage_stats``. Calling the formatter tracks it through any refactor, and
    ``calculate_vntr_coverage`` now writes exactly what it returns.

    Returns:
        list[str]: The column names ``calculate_vntr_coverage`` writes.
    """
    from vntyper.scripts.coverage_stats import COVERAGE_COLUMNS as PRODUCTION_COLUMNS
    from vntyper.scripts.coverage_stats import format_coverage_summary

    rendered = format_coverage_summary(dict.fromkeys(PRODUCTION_COLUMNS, 0))
    return rendered.splitlines()[0].split("\t")


# ---------------------------------------------------------------------------
# validate_coverage_output: the column names must be production's own
# ---------------------------------------------------------------------------


def test_the_helper_columns_are_the_ones_production_writes() -> None:
    """Derived from the production writer, so the two cannot drift apart.

    Nine since #172, which appended the ``coverage_qc`` verdict to the schema.
    """
    header = _production_coverage_header()
    assert len(header) == 16, f"the header scan produced {header}; it has drifted"
    assert list(COVERAGE_COLUMNS) == header


def test_coverage_output_returns_the_real_values(tmp_path) -> None:
    _write_tsv(
        tmp_path / "coverage" / "coverage_summary.tsv",
        list(COVERAGE_COLUMNS),
        # Eight measurements, then #222's seven build-comparable columns as NA - the
        # state a run with no `vntr_array_coords` produces - then the QC verdict.
        ["153.99", "150.00", "12.30", "7", "301", "5000", "12", "0.24"] + ["NA"] * 7 + ["PASS"],
    )

    metrics = validate_coverage_output(tmp_path)

    assert metrics["mean_cov"] == pytest.approx(153.99)
    assert metrics["median_cov"] == pytest.approx(150.00)
    assert metrics["uncovered_pct"] == pytest.approx(0.24)


def test_coverage_output_rejects_the_old_capitalised_column_names(tmp_path) -> None:
    """The regression this fix exists for: these names silently yielded 0.0."""
    _write_tsv(
        tmp_path / "coverage" / "coverage_summary.tsv",
        ["Mean", "Median", "StdDev", "Min", "Max", "RegionLength", "UncoveredBases", "Uncovered%"],
        ["153.99", "150.00", "12.30", "7", "301", "5000", "12", "0.24"],
    )

    with pytest.raises(AssertionError, match="missing"):
        validate_coverage_output(tmp_path)


@pytest.mark.parametrize("dropped", COVERAGE_COLUMNS)
def test_coverage_output_fails_loudly_on_any_missing_column(tmp_path, dropped: str) -> None:
    header = [column for column in COVERAGE_COLUMNS if column != dropped]
    _write_tsv(tmp_path / "coverage" / "coverage_summary.tsv", header, ["1"] * len(header))

    with pytest.raises(AssertionError, match=re.escape(dropped)):
        validate_coverage_output(tmp_path)


def test_coverage_output_fails_on_an_unparseable_value(tmp_path) -> None:
    row = ["not-a-number"] + ["1"] * (len(COVERAGE_COLUMNS) - 1)
    _write_tsv(tmp_path / "coverage" / "coverage_summary.tsv", list(COVERAGE_COLUMNS), row)

    with pytest.raises(AssertionError, match="mean"):
        validate_coverage_output(tmp_path)


def test_coverage_output_fails_when_there_are_no_data_rows(tmp_path) -> None:
    path = tmp_path / "coverage" / "coverage_summary.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(COVERAGE_COLUMNS) + "\n", encoding="utf-8")

    with pytest.raises(AssertionError, match="no data rows"):
        validate_coverage_output(tmp_path)


# ---------------------------------------------------------------------------
# assert_pattern_match: the star in High_Precision* is a character, not a glob
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", CONFIDENCE_LABELS)
def test_every_confidence_label_matches_itself(label: str) -> None:
    assert_pattern_match(label, label, "Confidence")


def test_high_precision_star_is_not_satisfied_by_high_precision() -> None:
    """The whole point: these are two different confidence levels.

    Under the old wildcard rule this passed, so no boundary between the two was
    ever asserted by either slow tier.
    """
    with pytest.raises(AssertionError):
        assert_pattern_match("High_Precision", "High_Precision*", "Confidence")


def test_high_precision_is_not_satisfied_by_high_precision_star() -> None:
    with pytest.raises(AssertionError):
        assert_pattern_match("High_Precision*", "High_Precision", "Confidence")


@pytest.mark.parametrize("actual", ["High_Precision_v2", "High_Precisionx", ""])
def test_a_trailing_star_is_not_a_wildcard_by_default(actual: str) -> None:
    with pytest.raises(AssertionError):
        assert_pattern_match(actual, "High_Precision*", "Confidence")


def test_a_bare_star_is_compared_literally() -> None:
    """`expected="*"` used to match everything; now it means the string "*"."""
    assert_pattern_match("*", "*", "Confidence")
    with pytest.raises(AssertionError):
        assert_pattern_match("anything at all", "*", "Confidence")


def test_a_bare_star_is_refused_even_when_wildcards_are_asked_for() -> None:
    """An opt-in wildcard that matches everything still asserts nothing."""
    with pytest.raises(AssertionError, match="asserts nothing"):
        assert_pattern_match("anything at all", "*", "Confidence", allow_wildcard=True)


def test_wildcard_matching_is_still_available_when_asked_for() -> None:
    assert_pattern_match("High_Precision_v2", "High_Precision*", "Confidence", allow_wildcard=True)
    with pytest.raises(AssertionError):
        assert_pattern_match("Low_Precision", "High_Precision*", "Confidence", allow_wildcard=True)


def test_the_failure_message_names_the_field_and_both_values() -> None:
    with pytest.raises(AssertionError) as excinfo:
        assert_pattern_match("High_Precision", "High_Precision*", "Confidence")
    message = str(excinfo.value)
    assert "Confidence" in message
    assert "High_Precision*" in message


# ---------------------------------------------------------------------------
# validate_kestrel_output
# ---------------------------------------------------------------------------


def _write_kestrel(tmp_path: Path, **columns: str) -> Path:
    """Write a minimal kestrel_result.tsv.

    Args:
        tmp_path: The output directory root.
        **columns: Column name to value.

    Returns:
        Path: ``tmp_path``, ready to hand to the validator.
    """
    _write_tsv(tmp_path / "kestrel" / "kestrel_result.tsv", list(columns), list(columns.values()))
    return tmp_path


def test_kestrel_output_rejects_an_empty_expectation_set(tmp_path) -> None:
    """An empty dict asserted nothing at all beyond the file being non-empty."""
    _write_kestrel(tmp_path, Confidence="High_Precision*")

    with pytest.raises(AssertionError, match="no expected values"):
        validate_kestrel_output(tmp_path, {})


def test_kestrel_output_compares_confidence_literally(tmp_path) -> None:
    """The star reaches the real assertion through the real call path."""
    _write_kestrel(tmp_path, Confidence="High_Precision")

    with pytest.raises(AssertionError):
        validate_kestrel_output(tmp_path, {"Confidence": "High_Precision*"})

    validate_kestrel_output(tmp_path, {"Confidence": "High_Precision"})


def test_kestrel_output_fails_on_a_missing_field(tmp_path) -> None:
    _write_kestrel(tmp_path, Confidence="High_Precision*")

    with pytest.raises(AssertionError, match="missing field: Depth_Score"):
        validate_kestrel_output(tmp_path, {"Depth_Score": 0.05})


def test_kestrel_output_applies_tolerance_to_the_depth_columns(tmp_path) -> None:
    _write_kestrel(
        tmp_path,
        Estimated_Depth_AlternateVariant="416",
        Estimated_Depth_Variant_ActiveRegion="7110",
        Depth_Score="0.0586",
    )

    validate_kestrel_output(
        tmp_path,
        {
            "Estimated_Depth_AlternateVariant": 416,
            "Estimated_Depth_Variant_ActiveRegion": {"value": 7000, "tolerance_percentage": 5},
            "Depth_Score": {"value": 0.0586, "tolerance_percentage": 5},
        },
    )

    with pytest.raises(AssertionError):
        validate_kestrel_output(tmp_path, {"Estimated_Depth_Variant_ActiveRegion": 7000})


def test_kestrel_output_fails_when_the_tsv_has_no_data_rows(tmp_path) -> None:
    path = tmp_path / "kestrel" / "kestrel_result.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# a comment\nConfidence\n", encoding="utf-8")

    with pytest.raises(AssertionError, match="no data rows"):
        validate_kestrel_output(tmp_path, {"Confidence": "High_Precision*"})
