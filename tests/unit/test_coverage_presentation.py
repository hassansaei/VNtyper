"""Presentation decisions for coverage quality-control verdicts."""

import pytest

from vntyper.scripts import screening_summary as ss
from vntyper.scripts.coverage_presentation import coverage_notes, coverage_qc_tone, coverage_uncovered_exceeded
from vntyper.scripts.coverage_qc import (
    COVERAGE_QC_FAIL,
    COVERAGE_QC_NOT_EVALUATED,
    COVERAGE_QC_PASS,
    COVERAGE_QC_REDUCED,
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
        (COVERAGE_QC_PASS, "Adequate"),
        (COVERAGE_QC_REDUCED, "Reduced"),
        (COVERAGE_QC_FAIL, "Insufficient"),
        (COVERAGE_QC_NOT_EVALUATED, "Not evaluated"),
        ("FUTURE_STATUS", "FUTURE_STATUS"),
    ],
)
def test_coverage_qc_word_translates_known_tokens_and_preserves_unknown_ones(status: str, expected: str) -> None:
    """Catch a future durable status being guessed at or discarded by the report."""
    assert ss.coverage_qc_word(status) == expected


@pytest.mark.parametrize(
    ("status", "expected"),
    [
        (COVERAGE_QC_PASS, "ok"),
        (COVERAGE_QC_REDUCED, "none"),
        (COVERAGE_QC_FAIL, "caution"),
        (COVERAGE_QC_NOT_EVALUATED, "none"),
        ("FUTURE_STATUS", "none"),
    ],
)
def test_only_a_failing_depth_is_toned_as_a_caution(status: str, expected: str) -> None:
    """Catch a reduced or unmeasured depth being painted as a pass or as a caution."""
    assert coverage_qc_tone(status) == expected


def _notes(report_config, qc, *, is_positive=False, mean=75.0, percent=5.0, low=50.0) -> tuple[str, ...]:
    return coverage_notes(
        report_config,
        qc,
        is_positive=is_positive,
        mean_vntr_coverage=mean,
        percent_vntr_uncovered=percent,
        mean_threshold=100,
        low_mean_threshold=low,
        percent_threshold=50.0,
    )


def test_a_reduced_depth_adds_one_sentence_to_a_result_without_a_finding(report_config) -> None:
    """A real exome at 96x was told "quality metrics are below threshold" and to re-assess."""
    qc = evaluate_coverage_qc(96.24, 8.31, 100, 50.0, low_mean_threshold=50.0)

    notes = _notes(report_config, qc, mean=96.24, percent=8.31)

    assert notes == (
        "VNTR depth is reduced: mean 96.24x, below 100x. Detection is somewhat less sensitive at this depth.",
    )


def test_a_low_depth_names_the_threshold_that_was_applied_and_an_action(report_config) -> None:
    with_band = evaluate_coverage_qc(42.9, 22.79, 100, 50.0, low_mean_threshold=50.0)
    single_line = evaluate_coverage_qc(75.0, 0.0, 100, 50.0)

    banded = _notes(report_config, with_band, mean=42.9, percent=22.79)
    legacy = _notes(report_config, single_line, mean=75.0, percent=0.0, low=None)

    assert banded[0].startswith("VNTR depth is low: mean 42.9x, below 50x.")
    assert legacy[0].startswith("VNTR depth is low: mean 75x, below 100x."), "the single line is the threshold"
    assert banded[-1] == legacy[-1] and "if a variant is suspected" in banded[-1]
    assert len(banded) == len(legacy) == 2


@pytest.mark.parametrize("mean", [250.0, 75.0, 10.0])
def test_a_call_is_never_qualified_by_its_mean_depth(report_config, mean) -> None:
    """Four of fourteen real positives at 69-87x were graded limited for depth alone."""
    qc = evaluate_coverage_qc(mean, 5.0, 100, 50.0, low_mean_threshold=50.0)

    assert _notes(report_config, qc, is_positive=True, mean=mean) == ()


@pytest.mark.parametrize("is_positive", [True, False])
def test_a_mostly_uncovered_region_is_stated_for_every_result_and_never_as_low_depth(
    report_config, is_positive
) -> None:
    """FAIL also means a patchy region; a 250x mean must not be explained as low depth."""
    qc = evaluate_coverage_qc(250.0, 80.0, 100, 50.0, low_mean_threshold=50.0)

    notes = _notes(report_config, qc, is_positive=is_positive, mean=250.0, percent=80.0)

    assert coverage_uncovered_exceeded(qc) is True
    assert notes[0].startswith("80% of the VNTR region has no read coverage (limit 50%).")
    assert not any("VNTR depth is" in note for note in notes)
    assert len(notes) == (1 if is_positive else 2)


def test_both_failures_are_both_stated_in_declaration_order(report_config) -> None:
    qc = evaluate_coverage_qc(10.0, 90.0, 100, 50.0, low_mean_threshold=50.0)

    notes = _notes(report_config, qc, mean=10.0, percent=90.0)

    assert [note.split(" ")[0] for note in notes] == ["90%", "VNTR", "Re‐sequence"]


def test_a_passing_or_unmeasured_verdict_has_no_note(report_config) -> None:
    assert _notes(report_config, evaluate_coverage_qc(250.0, 5.0, 100, 50.0, low_mean_threshold=50.0)) == ()
    assert _notes(report_config, evaluate_coverage_qc(None, None, 100, 50.0), mean=None, percent=None) == ()
    assert coverage_uncovered_exceeded(evaluate_coverage_qc(10.0, 5.0, 100, 50.0)) is False


def test_a_config_without_the_notes_block_renders_none(report_config) -> None:
    stripped = {key: value for key, value in report_config.items() if key != "coverage_notes"}
    qc = evaluate_coverage_qc(10.0, 90.0, 100, 50.0, low_mean_threshold=50.0)

    assert _notes(stripped, qc, mean=10.0, percent=90.0) == ()


@pytest.mark.parametrize(
    "block",
    [
        [],
        {"reduced_no_finding": "a"},
        {"reduced_no_finding": "a", "low_no_finding": "b", "uncovered": "c", "action_no_finding": "d", "extra": "e"},
        {"reduced_no_finding": "a", "low_no_finding": "b", "uncovered": "c", "action_no_finding": " "},
        {"reduced_no_finding": "a", "low_no_finding": "b", "uncovered": 3, "action_no_finding": "d"},
    ],
)
def test_a_malformed_notes_block_is_refused_rather_than_half_rendered(report_config, block) -> None:
    qc = evaluate_coverage_qc(75.0, 5.0, 100, 50.0, low_mean_threshold=50.0)

    with pytest.raises(ValueError, match="coverage_notes"):
        _notes({**report_config, "coverage_notes": block}, qc)
