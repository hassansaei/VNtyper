"""Unit tests for Kestrel selection configuration projection and strategy dispatch."""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from vntyper.scripts.kestrel_decision_config import (
    KestrelSortField,
    project_kestrel_selection,
)
from vntyper.scripts.kestrel_genotyping import (
    filter_final_dataframe,
    select_single_best_variant,
)

pytestmark = pytest.mark.unit


def _valid_selection_mapping() -> dict[str, Any]:
    return {
        "confidence_priority": {
            "High_Precision*": 3,
            "High_Precision": 2,
            "Low_Precision": 1,
            "Negative": 0,
        },
        "final_filter_columns": [
            "is_valid_frameshift",
            "depth_confidence_pass",
            "alt_filter_pass",
            "motif_filter_pass",
            "flag_filter_pass",
        ],
        "frameshift": {
            "modulus": 3,
            "insertion_remainder": 2,
            "deletion_remainder": 1,
        },
        "sort_order": [
            {"column": "confidence_priority", "ascending": False},
            {"column": "is_unflagged", "ascending": False},
            {"column": "Depth_Score", "ascending": False},
            {"column": "haplo_count", "ascending": False},
            {"column": "POS", "ascending": True},
        ],
        "unflagged_value": "Not flagged",
    }


def test_project_kestrel_selection_default_is_identity_dominance() -> None:
    raw = _valid_selection_mapping()
    projected = project_kestrel_selection(raw)
    assert projected.strategy == "identity_dominance"
    assert projected.modulus == 3
    assert projected.insertion_remainder == 2
    assert projected.deletion_remainder == 1
    assert projected.unflagged_value == "Not flagged"
    assert len(projected.sort_order) == 5
    assert projected.sort_order[0] == KestrelSortField("confidence_priority", False)


def test_project_kestrel_selection_explicit_legacy() -> None:
    raw = _valid_selection_mapping()
    raw["strategy"] = "legacy"
    projected = project_kestrel_selection(raw)
    assert projected.strategy == "legacy"


def test_project_kestrel_selection_identity_dominance() -> None:
    raw = _valid_selection_mapping()
    raw["strategy"] = "identity_dominance"
    projected = project_kestrel_selection(raw)
    assert projected.strategy == "identity_dominance"


def test_project_kestrel_selection_invalid_strategy() -> None:
    raw = _valid_selection_mapping()
    raw["strategy"] = "unsupported_xyz"
    with pytest.raises(ValueError, match="unsupported Kestrel selection strategy: 'unsupported_xyz'"):
        project_kestrel_selection(raw)


@pytest.mark.parametrize(
    "missing_key",
    ["confidence_priority", "final_filter_columns", "frameshift", "sort_order", "unflagged_value"],
)
def test_project_kestrel_selection_missing_keys(missing_key: str) -> None:
    raw = _valid_selection_mapping()
    del raw[missing_key]
    with pytest.raises(ValueError, match=f"Kestrel selection is missing '{missing_key}'"):
        project_kestrel_selection(raw)


def test_project_kestrel_selection_malformed_priority() -> None:
    raw = _valid_selection_mapping()
    raw["confidence_priority"] = "not-a-mapping"
    with pytest.raises(ValueError, match="priority and frameshift values must be mappings"):
        project_kestrel_selection(raw)


def test_project_kestrel_selection_malformed_columns() -> None:
    raw = _valid_selection_mapping()
    raw["final_filter_columns"] = "string_is_not_sequence"
    with pytest.raises(ValueError, match="final_filter_columns must be an ordered sequence"):
        project_kestrel_selection(raw)


def test_project_kestrel_selection_malformed_sort_order() -> None:
    raw = _valid_selection_mapping()
    raw["sort_order"] = "not_a_sequence"
    with pytest.raises(ValueError, match="sort_order must be an ordered sequence"):
        project_kestrel_selection(raw)

    raw["sort_order"] = ["not_a_mapping"]
    with pytest.raises(ValueError, match="sort_order entries must be mappings"):
        project_kestrel_selection(raw)

    raw["sort_order"] = [{"column": "col"}]  # missing ascending
    with pytest.raises(ValueError, match="require string column and boolean ascending"):
        project_kestrel_selection(raw)


def test_project_kestrel_selection_malformed_unflagged() -> None:
    raw = _valid_selection_mapping()
    raw["unflagged_value"] = 123
    with pytest.raises(ValueError, match="unflagged_value must be a string"):
        project_kestrel_selection(raw)


def test_select_single_best_variant_dispatches_between_legacy_and_identity_dominance() -> None:
    # Row 0: high raw Depth_Score, but single isolated assembly with low AltDepth
    # Row 1 & 2: lower Depth_Score each, but form a 2-assembly consensus group with higher peak AltDepth
    frame = pd.DataFrame(
        {
            "Confidence": ["High_Precision", "High_Precision", "High_Precision"],
            "Flag": ["Not flagged", "Not flagged", "Not flagged"],
            "Depth_Score": [0.090, 0.040, 0.038],
            "haplo_count": [10, 5, 5],
            "POS": [60, 53, 53],
            "REF": ["C", "C", "C"],
            "ALT": ["CG", "CG", "CG"],
            "Motifs": ["A-B", "C-D", "E-F"],
            "Estimated_Depth_AlternateVariant": [5.0, 25.0, 22.0],
            "Molecular_Identity": [
                "repeat:1:58:ins:G",
                "repeat:1:53:ins:G",
                "repeat:1:53:ins:G",
            ],
        }
    )

    # Legacy strategy selects Row 0 based on raw Depth_Score
    raw_legacy = _valid_selection_mapping()
    raw_legacy["strategy"] = "legacy"
    legacy_sel = project_kestrel_selection(raw_legacy)
    assert legacy_sel.strategy == "legacy"
    legacy_best = select_single_best_variant(frame, selection=legacy_sel)
    assert legacy_best.iloc[0]["Molecular_Identity"] == "repeat:1:58:ins:G"
    assert legacy_best.iloc[0]["Depth_Score"] == 0.090

    # Identity dominance strategy (now default) groups Rows 1 & 2 into repeat:1:53:ins:G (2 assemblies, peak AltDepth 25.0)
    # dominating Row 0 (1 assembly, peak AltDepth 5.0)
    dom_sel = project_kestrel_selection(_valid_selection_mapping())
    assert dom_sel.strategy == "identity_dominance"
    dom_best = select_single_best_variant(frame, selection=dom_sel)
    assert dom_best.iloc[0]["Molecular_Identity"] == "repeat:1:53:ins:G"
    assert dom_best.iloc[0]["Depth_Score"] == 0.040

    # Test explicit strategy override parameter
    override_best = select_single_best_variant(frame, strategy="legacy")
    assert override_best.iloc[0]["Molecular_Identity"] == "repeat:1:58:ins:G"


def test_filter_final_dataframe_with_identity_dominance(tmp_path: Any) -> None:
    frame = pd.DataFrame(
        {
            "is_frameshift": [True, True],
            "is_valid_frameshift": [True, True],
            "depth_confidence_pass": [True, True],
            "alt_filter_pass": [True, True],
            "motif_filter_pass": [True, True],
            "flag_filter_pass": [True, True],
            "Confidence": ["High_Precision", "High_Precision"],
            "Flag": ["Not flagged", "Not flagged"],
            "Depth_Score": [0.080, 0.030],
            "haplo_count": [10, 5],
            "POS": [60, 53],
            "REF": ["C", "C"],
            "ALT": ["CG", "CG"],
            "Motifs": ["A-B", "C-D"],
            "Estimated_Depth_AlternateVariant": [3.0, 30.0],
            "Molecular_Identity": [
                "repeat:1:58:ins:G",
                "repeat:1:53:ins:G",
            ],
        }
    )

    # Dominance strategy via kwarg
    filtered = filter_final_dataframe(frame, str(tmp_path), strategy="identity_dominance")
    assert len(filtered) == 1
    assert filtered.iloc[0]["Molecular_Identity"] == "repeat:1:53:ins:G"


def test_resolve_selection_handles_nested_and_strategy_replacement() -> None:
    from vntyper.scripts.kestrel_genotyping import _resolve_selection

    raw = {"selection": _valid_selection_mapping()}
    resolved = _resolve_selection(raw, strategy="identity_dominance")
    assert resolved.strategy == "identity_dominance"

    # Strategy change on existing KestrelSelection
    replaced = _resolve_selection(resolved, strategy="legacy")
    assert replaced.strategy == "legacy"


def test_identity_translation_component_validates_permit_boundary_insertions() -> None:
    from vntyper.scripts.identity_candidates import IdentityTranslationComponent

    with pytest.raises(ValueError, match="permit_boundary_insertions must be a boolean"):
        IdentityTranslationComponent(
            kestrel_motifs={"X": "ACGT"},
            advntr_repeat_unit_motifs={"2": "ACGT"},
            advntr_rotation_offset=1,
            permit_boundary_insertions="not-a-bool",  # type: ignore[arg-type]
        )
