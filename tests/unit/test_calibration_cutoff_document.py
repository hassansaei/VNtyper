"""Report-document assembly rules that need no optimize run: axis documents and curve availability."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import pytest

from tests.unit.test_calibration_cutoff_curves import scenario
from tests.unit.test_calibration_cutoff_optimize import _advntr_baseline_policy
from vntyper.scripts.calibration_cutoff_axes import ADVNTR_CUTOFF, ADVNTR_MIN_SUPPORT, declared_axis
from vntyper.scripts.calibration_cutoff_document import _CURVE_UNAVAILABLE, _axis_documents, _curve

pytestmark = pytest.mark.unit

#: One sample (p2) is a no-call at the strictest breakpoint and called at the others, so
#: the no-call set changes across the axis.
CHANGING_NO_CALLS: dict[int, dict[str, bool | None]] = {
    0: {"p1": True, "p2": True, "n1": False, "n2": False},
    1: {"p1": True, "p2": True, "n1": False, "n2": False},
    2: {"p1": True, "p2": None, "n1": False, "n2": False},
}


def _inputs(**fields: Any) -> Any:
    return SimpleNamespace(**fields)


def test_an_advntr_cutoff_axis_without_its_search_record_is_refused() -> None:
    """``unrejectable_samples`` comes from the probe search; a cutoff axis without one is a wiring defect."""
    axis = declared_axis(ADVNTR_CUTOFF, [0.001, 0.004], baseline=_advntr_baseline_policy())

    with pytest.raises(ValueError, match="adVNTR axis advntr_cutoff has no adVNTR search record"):
        _axis_documents(_inputs(derived=[(axis, ())], advntr_search=None))


def test_only_the_advntr_cutoff_axis_carries_unrejectable_samples() -> None:
    """Spec section 8: a p-value of 0 defeats every cutoff; read support has no such floor."""
    baseline = _advntr_baseline_policy()
    cutoff = declared_axis(ADVNTR_CUTOFF, [0.001, 0.004], baseline=baseline)
    support = declared_axis(ADVNTR_MIN_SUPPORT, [3, 4], baseline=baseline)
    search = SimpleNamespace(unrejectable={ADVNTR_CUTOFF: 2, ADVNTR_MIN_SUPPORT: 0})

    documents = _axis_documents(_inputs(derived=[(cutoff, ()), (support, ())], advntr_search=search))

    assert documents[0]["unrejectable_samples"] == 2
    assert "unrejectable_samples" not in documents[1]


@pytest.mark.parametrize("caller", ["kestrel", "advntr"])
def test_a_single_caller_curve_with_a_changing_no_call_set_is_still_refused(caller: str) -> None:
    """Only the either-caller union may change its no-call set; for one caller it is a replay defect."""
    axis, candidates, arms = scenario(CHANGING_NO_CALLS)

    with pytest.raises(ValueError, match="cannot hide changes in assessability"):
        _curve(axis, candidates, arms, ">=", caller)


def test_a_union_curve_with_a_changing_no_call_set_is_published_unavailable() -> None:
    axis, candidates, arms = scenario(CHANGING_NO_CALLS)

    assert _curve(axis, candidates, arms, ">=", "both") == {
        "axis": axis.axis,
        "status": "unavailable",
        "reason": _CURVE_UNAVAILABLE,
    }


@pytest.mark.parametrize("caller", ["kestrel", "both"])
def test_a_curve_with_a_fixed_no_call_set_is_built_and_marked_available(caller: str) -> None:
    fixed = {index: {**calls, "p2": True} for index, calls in CHANGING_NO_CALLS.items()}
    axis, candidates, arms = scenario(fixed)

    document = _curve(axis, candidates, arms, ">=", caller)

    assert document["status"] == "available"
    assert document["axis"] == axis.axis
    assert len(document["points"]) == 3
