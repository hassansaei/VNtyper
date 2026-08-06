"""Rolling row-level verdicts up into the three categories the donut charts count.

`cohort_rules.compute_algorithm_result` answers a question about one *row*. A cohort
report answers a question about one *sample*, and a sample can have several rows. This
module is the reduction in between: map each row's verdict onto Positive /
Positive_Flagged / Negative, then take the highest category any of the sample's rows
reached.

It is the step where a diagnosis is made, so it is also the step where a refactor can
change one without changing a line of report markup. Everything here is
**characterisation** - it records the reduction as it stands - with the single exception
of `test_the_reduction_leaves_the_caller_s_frame_untouched`, which is a specification:
it states the contract that closed the internal-column leak, and it is named as such in
its own docstring. No other test in this file has been ratified.
"""

from __future__ import annotations

import pandas as pd
import pytest

from vntyper.scripts.cohort_categories import (
    aggregate_sample_category,
    category_counts,
    sample_categories,
    unify_advntr_result,
    unify_kestrel_result,
)

pytestmark = pytest.mark.unit

#: The shipped Kestrel rules' four possible results, plus the default.
_KESTREL_LOGIC = {
    "rules": [
        {"conditions": {"Confidence": {"operator": "in", "value": ["High_Precision"]}}, "result": "High_Precision"},
        {"conditions": {"Confidence": {"operator": "in", "value": ["Low_Precision"]}}, "result": "Low_Precision"},
        {"conditions": {"Confidence": {"operator": "in", "value": ["Flagged"]}}, "result": "Low_Precision_flagged"},
    ],
    "default": "negative",
}


# ---------------------------------------------------------------------------
# Row-level verdict -> category
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "row_result,category",
    [
        ("High_Precision", "Positive"),
        ("Low_Precision", "Positive"),
        ("High_Precision_flagged", "Positive_Flagged"),
        ("Low_Precision_flagged", "Positive_Flagged"),
        ("negative", "Negative"),
        ("", "Negative"),
        (None, "Negative"),
    ],
)
def test_each_kestrel_verdict_maps_to_its_category(row_result: str | None, category: str) -> None:
    assert unify_kestrel_result(row_result) == category


def test_a_starred_confidence_label_is_not_a_kestrel_verdict() -> None:
    """`High_Precision*` is a *confidence* label, and this function takes *verdicts*.

    The shipped rules already fold `High_Precision*` into the `High_Precision` verdict,
    so nothing reaches here with a star. Pinned because passing one straight in - the
    obvious mistake when wiring a new caller - is silently counted as Negative rather
    than rejected.
    """
    assert unify_kestrel_result("High_Precision*") == "Negative"


@pytest.mark.parametrize(
    "row_result,category",
    [
        ("positive", "Positive"),
        ("positive flagged", "Positive_Flagged"),
        ("negative", "Negative"),
        ("Positive", "Negative"),
        (None, "Negative"),
    ],
)
def test_each_advntr_verdict_maps_to_its_category(row_result: str | None, category: str) -> None:
    """The adVNTR verdicts are lower-case and space-separated where Kestrel's are
    capitalised and underscored, so the two mappings cannot be shared."""
    assert unify_advntr_result(row_result) == category


# ---------------------------------------------------------------------------
# Row categories -> one sample category
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "results,category",
    [
        ([], "Negative"),
        (["Negative"], "Negative"),
        (["Negative", "Negative"], "Negative"),
        (["Positive_Flagged", "Negative"], "Positive_Flagged"),
        (["Positive", "Negative"], "Positive"),
        (["Negative", "Positive_Flagged", "Positive"], "Positive"),
        (["Positive_Flagged", "Positive_Flagged"], "Positive_Flagged"),
    ],
)
def test_a_sample_takes_the_highest_category_any_of_its_rows_reached(results: list[str], category: str) -> None:
    assert aggregate_sample_category(results) == category


def test_a_sample_with_no_rows_at_all_is_negative() -> None:
    """A sample whose every row was filtered out counts as Negative, not as absent."""
    assert aggregate_sample_category([]) == "Negative"


# ---------------------------------------------------------------------------
# The frame-level reduction
# ---------------------------------------------------------------------------


def test_an_empty_frame_produces_an_empty_series() -> None:
    result = sample_categories(pd.DataFrame(), _KESTREL_LOGIC, unify_kestrel_result)

    assert result.empty


def test_a_frame_with_no_sample_column_produces_an_empty_series() -> None:
    """Without `Sample` there is nothing to group by, so the whole cohort drops out
    rather than being counted as one anonymous sample."""
    frame = pd.DataFrame([{"Confidence": "High_Precision"}])

    result = sample_categories(frame, _KESTREL_LOGIC, unify_kestrel_result)

    assert result.empty


def test_rows_are_grouped_by_sample_and_reduced() -> None:
    frame = pd.DataFrame(
        [
            {"Sample": "s1", "Confidence": "Flagged"},
            {"Sample": "s1", "Confidence": "High_Precision"},
            {"Sample": "s2", "Confidence": "Flagged"},
            {"Sample": "s3", "Confidence": "Nothing"},
        ]
    )

    result = sample_categories(frame, _KESTREL_LOGIC, unify_kestrel_result)

    assert result.to_dict() == {"s1": "Positive", "s2": "Positive_Flagged", "s3": "Negative"}


def test_the_result_is_indexed_by_sample_name() -> None:
    frame = pd.DataFrame([{"Sample": "s1", "Confidence": "High_Precision"}])

    result = sample_categories(frame, _KESTREL_LOGIC, unify_kestrel_result)

    assert list(result.index) == ["s1"]


def test_the_reduction_leaves_the_caller_s_frame_untouched() -> None:
    """**Specification**, and the contract this module's whole job statement implies.

    A function documented as reducing a frame to categories has no business writing to
    it. It used to: the two working columns `__row_result` and `__unified` went onto the
    frame the caller passed rather than onto a copy, and because `aggregate_cohort`
    renders the report *before* it exports CSV/TSV/JSON from the same frame, both
    reached every machine-readable output. Fixed by annotating a copy, so the caller's
    frame comes back exactly as it went in - columns, order and values. See
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-internal-columns-leak-into-exports.md`.

    The frame is asserted whole rather than by column name, because "no new columns" is
    the visible half of the contract and "no changed cells" is the other half.
    """
    frame = pd.DataFrame([{"Sample": "s1", "Confidence": "High_Precision"}])
    before = frame.copy(deep=True)

    result = sample_categories(frame, _KESTREL_LOGIC, unify_kestrel_result)

    assert list(frame.columns) == ["Sample", "Confidence"]
    pd.testing.assert_frame_equal(frame, before)
    assert result.to_dict() == {"s1": "Positive"}


# ---------------------------------------------------------------------------
# Categories -> the three donut segments
# ---------------------------------------------------------------------------


def test_an_empty_series_counts_as_three_zeroes_and_a_zero_total() -> None:
    assert category_counts(pd.Series(dtype=str)) == (0, 0, 0, 0)


def test_the_counts_come_back_in_segment_order_with_their_total() -> None:
    """Positive, Positive_Flagged, Negative, total - the order the donut expects."""
    series = pd.Series(["Positive", "Negative", "Positive_Flagged", "Positive", "Negative"])

    assert category_counts(series) == (2, 1, 2, 5)


def test_a_category_the_reduction_never_produces_is_not_counted() -> None:
    """The total is the sum of the three known categories, not the length of the
    series, so an unexpected label is dropped from the donut rather than mis-attributed
    to one of the three segments."""
    series = pd.Series(["Positive", "Unexpected"])

    assert category_counts(series) == (1, 0, 0, 1)
