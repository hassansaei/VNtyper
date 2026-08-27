"""Pure presentation-state tests for the cross-match report section."""

from __future__ import annotations

import pytest

from vntyper.scripts.cross_match_presentation import build_cross_match_summary
from vntyper.scripts.screening_summary import load_report_config
from vntyper.scripts.summary_steps import STEP_CROSS_MATCH

pytestmark = pytest.mark.unit

NO_MATCH_MESSAGE = "No matches were found between Kestrel and adVNTR results."


def _step(parsed_result: object, **extra: object) -> dict[str, object]:
    """Build one recorded cross-match step."""
    return {"step": STEP_CROSS_MATCH, "parsed_result": parsed_result, **extra}


@pytest.mark.parametrize(
    ("rows", "expected_message", "expected_positive"),
    [
        ([{"Match": "Yes"}], "At least one match was found between Kestrel and adVNTR results.", True),
        ([{"Match": "No"}], "No matches were found between Kestrel and adVNTR results.", False),
        (
            [{"Match": "No"}, {"Match": "Yes"}],
            "At least one match was found between Kestrel and adVNTR results.",
            True,
        ),
        ([{"Match": "no"}], "No matches were found between Kestrel and adVNTR results.", False),
        ([{}], "No matches were found between Kestrel and adVNTR results.", False),
        ([], "No matches were found between Kestrel and adVNTR results.", False),
    ],
)
def test_read_cross_match_rows_keep_the_existing_wording_and_verdict(
    rows: list[dict[str, str]], expected_message: str, expected_positive: bool
) -> None:
    """Readable comparison wording and its structural state are independent contracts."""
    summary = {"steps": [_step({"data": rows})]}

    message, is_positive, is_assessable = build_cross_match_summary(summary, {})

    assert message == expected_message
    assert is_positive is expected_positive
    assert is_assessable is True


def test_absent_cross_match_has_no_section_or_verdict() -> None:
    """A stage that did not run produces no section and no assessability claim."""
    assert build_cross_match_summary({"steps": []}, {}) == ("", False, False)


def test_unreadable_cross_match_uses_the_configured_not_assessable_message() -> None:
    """Unreadable data establishes neither agreement nor disagreement."""
    summary = {"steps": [_step(None)]}
    config = {"cross_match": {"not_assessable_message": "Comparison result unavailable."}}

    assert build_cross_match_summary(summary, config) == ("Comparison result unavailable.", False, False)


def test_shipped_config_declares_the_unreadable_cross_match_message() -> None:
    """The shipped report opts into explicit unreadable-comparison wording."""
    config = load_report_config()
    summary = {"steps": [_step(None)]}

    message, is_positive, is_assessable = build_cross_match_summary(summary, config)

    assert message == config["cross_match"]["not_assessable_message"]
    assert "not assessed" in message
    assert (is_positive, is_assessable) == (False, False)


@pytest.mark.parametrize("cross_match_config", [None, {}, "invalid", {"not_assessable_message": ""}])
def test_legacy_or_malformed_config_keeps_the_old_sentence_but_never_claims_assessability(
    cross_match_config: object,
) -> None:
    """Compatibility may preserve text, but it cannot restore the false No-match chip."""
    summary = {"steps": [_step({"error": "unreadable"})]}
    config = {} if cross_match_config is None else {"cross_match": cross_match_config}

    assert build_cross_match_summary(summary, config) == (NO_MATCH_MESSAGE, False, False)
