"""Unit tests for the golden-cohort gate's declared-delta waiver policy.

The policy is small and load-bearing: it decides which deltas fail the gate, and a
mistake here does not fail loudly - it turns a gate off. These tests pin the two halves
that must not drift apart. :func:`golden_cohort.waiver.fatal_deltas` decides fatality and
:func:`golden_cohort.waiver.waived_cases` decides what the report says was waived; if the
second stopped agreeing with the first, a run could reach ``IDENTICAL`` by waiving a delta
and never name the case whose delta it waived.

The fatality rules themselves are exercised against ``compare`` in
``test_golden_cohort_compare.py``, through the re-export this module's first test pins.
"""

import sys
from pathlib import Path
from typing import Any

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from golden_cohort import compare, waiver  # noqa: E402


def _result(cases: dict[str, tuple[list[str], list[str]]]) -> dict[str, Any]:
    """Build a comparison document carrying only what the waiver query reads.

    Args:
        cases: Case id -> (all deltas, fatal deltas).

    Returns:
        dict[str, Any]: A document in ``compare_sides`` shape, cases only.
    """
    return {"cases": {case_id: {"deltas": deltas, "fatal_deltas": fatal} for case_id, (deltas, fatal) in cases.items()}}


def test_the_waiver_policy_is_defined_once_and_re_exported_by_compare() -> None:
    """``compare`` must re-export the policy, not carry a second copy of it.

    The names are part of ``compare``'s surface and callers use them from there. Binding
    the same objects is what keeps the extraction a move rather than a fork: two
    definitions of "which delta may be waived" would be free to drift.
    """
    assert compare.fatal_deltas is waiver.fatal_deltas
    assert compare.DECLARABLE_DELTA is waiver.DECLARABLE_DELTA
    assert compare.DECLARED_DELTA_HEADING is waiver.DECLARED_DELTA_HEADING


def test_waived_cases_names_a_case_whose_command_delta_was_waived() -> None:
    """A waived delta is still a delta, and the reader has to be told which case had one."""
    result = _result({"a5c1_hg19_subset": ([waiver.DECLARABLE_DELTA], [])})

    assert waiver.waived_cases(result) == ["a5c1_hg19_subset"]


def test_waived_cases_ignores_a_command_delta_that_stayed_fatal() -> None:
    """Without the declaration nothing is waived, so the section must not appear."""
    result = _result({"a5c1_hg19_subset": ([waiver.DECLARABLE_DELTA], [waiver.DECLARABLE_DELTA])})

    assert waiver.waived_cases(result) == []


def test_waived_cases_ignores_a_case_whose_only_delta_is_someone_else() -> None:
    """A fatal genotype delta is not a waiver, and must never be reported as one."""
    result = _result({"a5c1_hg19_subset": (["kestrel_result"], ["kestrel_result"])})

    assert waiver.waived_cases(result) == []


def test_waived_cases_are_sorted_so_two_runs_read_the_same() -> None:
    """The rendered report is diffed between runs; case order cannot be dict order."""
    result = _result(
        {
            "b178_hg19_subset": ([waiver.DECLARABLE_DELTA], []),
            "a5c1_hg19_subset": ([waiver.DECLARABLE_DELTA], []),
        }
    )

    assert waiver.waived_cases(result) == ["a5c1_hg19_subset", "b178_hg19_subset"]
