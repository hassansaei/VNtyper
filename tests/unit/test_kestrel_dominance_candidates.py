"""Complete generated-dominance candidate frame projections."""

import pandas as pd
import pytest

from vntyper.scripts.identity_candidate_persistence import IDENTITY_SELECTION_COLUMNS
from vntyper.scripts.kestrel_dominance_candidates import (
    merge_candidate_annotations,
    passing_candidate_frame,
    selected_candidate_frame,
)

pytestmark = pytest.mark.unit


def _frame() -> pd.DataFrame:
    rows = []
    for ordinal, passing in enumerate((True, True, False)):
        row: dict[str, object] = {
            "__Identity_Observation_Ordinal": str(ordinal),
            "is_frameshift": str(passing),
        }
        row.update(dict.fromkeys(IDENTITY_SELECTION_COLUMNS, "owned" if passing else ""))
        rows.append(row)
    return pd.DataFrame(rows)


def test_passing_candidate_frame_keeps_all_positive_hypotheses_before_selection() -> None:
    result = passing_candidate_frame(_frame(), ("is_frameshift",))

    assert result["__Identity_Observation_Ordinal"].tolist() == ["0", "1"]


@pytest.mark.parametrize("mutation", ["missing", "malformed"])
def test_passing_candidate_frame_fails_closed_on_invalid_prerequisites(mutation: str) -> None:
    frame = _frame()
    if mutation == "missing":
        frame = frame.drop(columns=[IDENTITY_SELECTION_COLUMNS[0]])
    else:
        frame.loc[0, "is_frameshift"] = "yes"

    with pytest.raises(ValueError, match="required columns|Boolean"):
        passing_candidate_frame(frame, ("is_frameshift",))


def test_merge_candidate_annotations_retains_failed_evidence_rows() -> None:
    frame = _frame()
    annotated = frame.iloc[:2].copy()
    annotated["Nomenclature"] = ["A", "B"]

    merged = merge_candidate_annotations(frame, annotated)

    assert len(merged) == 3
    assert merged["Nomenclature"].tolist() == ["A", "B", ""]


def test_merge_candidate_annotations_rejects_foreign_rows() -> None:
    annotated = _frame().iloc[:1].copy()
    annotated.index = [99]

    with pytest.raises(ValueError, match="belong"):
        merge_candidate_annotations(_frame(), annotated)


def test_selected_candidate_frame_uses_the_legacy_ordinal_without_reranking() -> None:
    result = selected_candidate_frame(_frame().iloc[:2], 1)

    assert result["__Identity_Observation_Ordinal"].tolist() == ["1"]


@pytest.mark.parametrize("ordinal", [-1, 8])
def test_selected_candidate_frame_rejects_invalid_or_absent_ordinal(ordinal: int) -> None:
    with pytest.raises(ValueError, match="ordinal"):
        selected_candidate_frame(_frame(), ordinal)
