"""The extracted Kestrel evaluator preserves complete decision evidence."""

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.scripts.identity_candidates import translation_component_from_config
from vntyper.scripts.kestrel_decision_config import KestrelSelection
from vntyper.scripts.kestrel_genotyping import (
    _resolve_selection,
    add_haplo_count,
    process_kmer_results,
    select_single_best_variant,
)
from vntyper.scripts.kestrel_postprocessing import (
    evaluate_kestrel_candidates,
    filter_and_select_kestrel_candidates,
)
from vntyper.scripts.nomenclature import nomenclature_config

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


def test_live_postprocessing_observer_receives_complete_defensive_raw_inputs_before_decisions(tmp_path) -> None:
    first = kestrel_stage_frame("raw")
    suppressed = kestrel_stage_frame("raw", ref="C", alt="CGGCA")
    suppressed.loc[0, "POS"] = 68
    raw = pd.concat([first, suppressed], ignore_index=True)
    identity_motifs = nomenclature_config["motifs"]
    assert isinstance(identity_motifs, dict)
    raw["Motifs"] = "S-C"
    raw["Motif_sequence"] = identity_motifs["C"] + identity_motifs["S"]
    motifs = pd.DataFrame({"Motif": ["S"], "Motif_sequence": [identity_motifs["S"]]})
    config = kestrel_config()
    identity = translation_component_from_config(nomenclature_config)
    observed: dict[str, object] = {}

    def observe(raw_frame, motif_frame, frozen_config, selection, identity_component) -> None:
        observed["raw"] = raw_frame.copy(deep=True)
        observed["motifs"] = motif_frame.copy(deep=True)
        observed["selection"] = selection
        observed["identity"] = identity_component
        raw_frame.loc[0, "ALT"] = "MUTATED"
        motif_frame.loc[0, "Motif"] = "MUTATED"
        frozen_config["artifact_flags"] = []

    result = process_kmer_results(
        raw,
        motifs,
        str(tmp_path),
        config,
        identity_component=identity,
        raw_capture_observer=observe,
    )

    observed_raw = observed["raw"]
    observed_motifs = observed["motifs"]
    assert isinstance(observed_raw, pd.DataFrame) and isinstance(observed_motifs, pd.DataFrame)
    pd.testing.assert_frame_equal(observed_raw, raw)
    pd.testing.assert_frame_equal(observed_motifs, motifs)
    assert observed["identity"] is identity
    assert observed["selection"] == _resolve_selection(config)
    assert result.iloc[0]["ALT"] != "MUTATED"
    prefilter = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
    assert len(prefilter) == 2
    assert bool(prefilter.iloc[1]["flag_filter_pass"]) is False


def test_live_raw_capture_observer_requires_explicit_frozen_identity(tmp_path) -> None:
    with pytest.raises(ValueError, match="identity"):
        process_kmer_results(
            kestrel_stage_frame("raw"),
            _motifs(),
            str(tmp_path),
            kestrel_config(),
            raw_capture_observer=lambda *_args: None,
        )
