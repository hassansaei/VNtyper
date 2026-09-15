"""The extracted Kestrel evaluator preserves complete decision evidence."""

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.scripts.kestrel_decision_config import KestrelSelection
from vntyper.scripts.kestrel_genotyping import _resolve_selection, add_haplo_count, select_single_best_variant
from vntyper.scripts.kestrel_postprocessing import (
    evaluate_kestrel_candidates,
    filter_and_select_kestrel_candidates,
)

pytestmark = pytest.mark.unit


def _motifs() -> pd.DataFrame:
    return pd.DataFrame({"Motif": ["X", "5"], "Motif_sequence": ["SEQ1", "SEQ2"]})


def _select(frame: pd.DataFrame, selection: KestrelSelection) -> pd.DataFrame:
    return select_single_best_variant(frame, selection=selection)


def test_evaluator_retains_an_artifact_row_that_the_final_result_suppresses() -> None:
    frame = kestrel_stage_frame("raw", ref="C", alt="CGGCA")
    selection = _resolve_selection(kestrel_config())

    result = evaluate_kestrel_candidates(
        frame,
        _motifs(),
        kestrel_config(),
        selection=selection,
        add_haplo_count_fn=add_haplo_count,
        select_single_best_variant_fn=_select,
    )

    assert result.reached_final_filter is True
    assert len(result.prefilter) == 1
    assert bool(result.prefilter.iloc[0]["flag_filter_pass"]) is False
    assert result.prefilter.iloc[0]["Flag"] == "False_Positive_4bp_Insertion"
    assert result.selected.empty


def test_evaluator_preserves_ordered_duplicate_source_ordinals() -> None:
    frame = kestrel_stage_frame("raw", rows=2)
    frame.loc[1, "POS"] = frame.loc[0, "POS"]
    frame.loc[1, "Sample"] = "Del:121:12000"
    selection = _resolve_selection(kestrel_config())

    result = evaluate_kestrel_candidates(
        frame,
        _motifs(),
        kestrel_config(),
        selection=selection,
        add_haplo_count_fn=add_haplo_count,
        select_single_best_variant_fn=_select,
    )

    assert result.prefilter["Sample"].tolist() == ["Del:120:12000", "Del:121:12000"]
    assert result.prefilter["haplo_count"].tolist() == [2, 2]
    assert len(result.selected) == 1


def test_pure_final_filter_requires_every_gate_before_selection() -> None:
    frame = kestrel_stage_frame("flagged").drop(columns=["motif_filter_pass"])
    selection = _resolve_selection(kestrel_config())

    with pytest.raises(ValueError, match="motif_filter_pass"):
        filter_and_select_kestrel_candidates(
            frame,
            selection=selection,
            select_single_best_variant_fn=_select,
        )


def test_empty_evaluator_never_claims_to_have_reached_the_final_gate() -> None:
    empty = kestrel_stage_frame("raw").iloc[0:0]
    selection = _resolve_selection(kestrel_config())

    result = evaluate_kestrel_candidates(
        empty,
        _motifs(),
        kestrel_config(),
        selection=selection,
        add_haplo_count_fn=add_haplo_count,
        select_single_best_variant_fn=_select,
    )

    assert result.reached_final_filter is False
    assert result.prefilter.empty
    assert result.selected.empty
