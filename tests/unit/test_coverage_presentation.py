"""Presentation decisions for coverage quality-control verdicts."""

import pytest

from vntyper.scripts import screening_summary as ss
from vntyper.scripts.coverage_qc import (
    COVERAGE_QC_FAIL,
    COVERAGE_QC_NOT_EVALUATED,
    COVERAGE_QC_PASS,
    evaluate_coverage_qc,
)

pytestmark = pytest.mark.unit


@pytest.fixture(scope="module")
def report_config() -> dict:
    """The shipped report configuration."""
    return ss.load_report_config()


def test_the_note_renders_only_when_nothing_was_measured(report_config) -> None:
    """Catch deriving the note from ``passed``, which is true for NOT_EVALUATED."""
    unmeasured = evaluate_coverage_qc(None, None, 100.0, 50.0)
    measured_pass = evaluate_coverage_qc(250.0, 0.5, 100.0, 50.0)
    measured_fail = evaluate_coverage_qc(10.0, 80.0, 100.0, 50.0)

    note = ss.coverage_not_measured_note(report_config, unmeasured)

    assert note, "the shipped config must declare the note"
    assert ss.coverage_not_measured_note(report_config, measured_pass) == ""
    assert ss.coverage_not_measured_note(report_config, measured_fail) == ""


def test_a_config_written_before_the_note_renders_nothing(report_config) -> None:
    """Catch emitting new prose for an older deployment that did not opt in."""
    stripped = {key: value for key, value in report_config.items() if key != ss.COVERAGE_NOT_MEASURED_NOTE_KEY}
    unmeasured = evaluate_coverage_qc(None, None, 100.0, 50.0)

    assert ss.coverage_not_measured_note(stripped, unmeasured) == ""


@pytest.mark.parametrize(
    ("status", "expected"),
    [
        (COVERAGE_QC_PASS, "Pass"),
        (COVERAGE_QC_FAIL, "Fail"),
        (COVERAGE_QC_NOT_EVALUATED, "Not evaluated"),
        ("FUTURE_STATUS", "FUTURE_STATUS"),
    ],
)
def test_coverage_qc_word_translates_known_tokens_and_preserves_unknown_ones(status: str, expected: str) -> None:
    """Catch a future durable status being guessed at or discarded by the report."""
    assert ss.coverage_qc_word(status) == expected
