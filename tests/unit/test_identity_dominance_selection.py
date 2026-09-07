"""Unit tests for pure candidate grouping, metric aggregation, and identity dominance selection."""

from __future__ import annotations

import pandas as pd
import pytest

from vntyper.scripts.identity_dominance_selection import (
    IdentityCandidateGroup,
    group_candidates_by_identity,
    rank_candidate_groups,
    select_group_representative_row,
    select_identity_dominance_variant,
)

pytestmark = pytest.mark.unit


def test_empty_dataframe_returns_empty_groups_and_empty_selection() -> None:
    """An empty candidate pool yields no groups and returns an empty DataFrame."""
    empty_df = pd.DataFrame()
    assert group_candidates_by_identity(empty_df) == []
    assert select_identity_dominance_variant(empty_df).empty


def test_single_row_dataframe_returns_immediately() -> None:
    """A candidate pool with a single row short-circuits and returns unchanged."""
    single = pd.DataFrame([{"POS": 60, "REF": "C", "ALT": "CT", "Confidence": "High_Precision*"}])
    result = select_identity_dominance_variant(single)
    assert len(result) == 1
    assert result.iloc[0]["POS"] == 60


def test_promiscuous_artifact_cannot_defeat_true_variant_via_summing() -> None:
    """Peak alternate depth (max) ensures a high-depth true mutation defeats promiscuous low-depth stutter."""
    # True variant: present in only 2 motif pairs, but each has high depth (120)
    true_rows = [
        {
            "Molecular_Identity": "MUC1-X-60-coding-v1|55|54|-|G",
            "Molecular_Identity_Context_Diverges": False,
            "Confidence": "High_Precision*",
            "Estimated_Depth_AlternateVariant": 120,
            "Depth_Score": 0.05,
            "haplo_count": 2,
            "POS": 66,
            "REF": "G",
            "ALT": "GC",
            "Flag": "Not flagged",
        },
        {
            "Molecular_Identity": "MUC1-X-60-coding-v1|55|54|-|G",
            "Molecular_Identity_Context_Diverges": False,
            "Confidence": "High_Precision*",
            "Estimated_Depth_AlternateVariant": 115,
            "Depth_Score": 0.048,
            "haplo_count": 2,
            "POS": 66,
            "REF": "G",
            "ALT": "GC",
            "Flag": "Not flagged",
        },
    ]

    # Stutter artifact: present across 10 motif pairs, each with low depth (20)
    artifact_rows = [
        {
            "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|C",
            "Molecular_Identity_Context_Diverges": False,
            "Confidence": "High_Precision*",
            "Estimated_Depth_AlternateVariant": 20,
            "Depth_Score": 0.01,
            "haplo_count": 10,
            "POS": 65,
            "REF": "G",
            "ALT": "GG",
            "Flag": "Not flagged",
        }
        for _ in range(10)
    ]

    df = pd.DataFrame(true_rows + artifact_rows)
    groups = group_candidates_by_identity(df)
    assert len(groups) == 2

    # Group for true variant should have peak_alternate_depth = 120, not 235
    true_group = next(g for g in groups if "55|54" in g.identity_key)
    artifact_group = next(g for g in groups if "60|59" in g.identity_key)

    assert true_group.peak_alternate_depth == 120.0
    assert artifact_group.peak_alternate_depth == 20.0

    # If peak depth is prioritized or when assembly counts are evaluated,
    # peak alternate depth prevents the 10x artifact sum (200) from beating true variant (120)
    ranked = rank_candidate_groups([true_group, artifact_group])
    assert ranked[0].identity_key == true_group.identity_key
    assert true_group.peak_alternate_depth > artifact_group.peak_alternate_depth


def test_context_fidelity_prefers_canonical_x_over_divergent_motif() -> None:
    """When confidence is tied, canonical X repeat context beats divergent motif context."""
    divergent_row = {
        "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|C",
        "Molecular_Identity_Context_Diverges": True,
        "Confidence": "High_Precision*",
        "Estimated_Depth_AlternateVariant": 276,
        "Depth_Score": 0.0258,
        "haplo_count": 26,
        "POS": 67,
        "REF": "G",
        "ALT": "GG",
        "Flag": "Not flagged",
    }
    canonical_row = {
        "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|A",
        "Molecular_Identity_Context_Diverges": False,
        "Confidence": "High_Precision*",
        "Estimated_Depth_AlternateVariant": 276,
        "Depth_Score": 0.0251,
        "haplo_count": 132,
        "POS": 60,
        "REF": "C",
        "ALT": "CT",
        "Flag": "Not flagged",
    }

    df = pd.DataFrame([divergent_row, canonical_row])
    groups = group_candidates_by_identity(df)
    assert len(groups) == 2

    div_group = next(g for g in groups if g.identity_key.endswith("|-|C"))
    can_group = next(g for g in groups if g.identity_key.endswith("|-|A"))

    assert div_group.context_fidelity == 0
    assert can_group.context_fidelity == 1

    ranked = rank_candidate_groups(groups)
    assert ranked[0].identity_key == can_group.identity_key

    best = select_identity_dominance_variant(df)
    assert len(best) == 1
    assert best.iloc[0]["Molecular_Identity"] == "MUC1-X-60-coding-v1|60|59|-|A"


def test_unresolved_rows_group_distinctly_by_coordinates() -> None:
    """Rows with missing molecular identities group safely by distinct coordinates."""
    row_1 = {"POS": 60, "REF": "C", "ALT": "CT", "Confidence": "Low_Precision"}
    row_2 = {"POS": 60, "REF": "C", "ALT": "CT", "Confidence": "Low_Precision"}
    row_3 = {"POS": 54, "REF": "C", "ALT": "CA", "Confidence": "Low_Precision"}

    df = pd.DataFrame([row_1, row_2, row_3])
    groups = group_candidates_by_identity(df)
    assert len(groups) == 2
    keys = {g.identity_key for g in groups}
    assert "unresolved:60:C:CT" in keys
    assert "unresolved:54:C:CA" in keys

    grp_60 = next(g for g in groups if g.identity_key == "unresolved:60:C:CT")
    assert grp_60.assembly_cardinality == 2


def test_representative_row_selection_prefers_unflagged_and_high_haplo_count() -> None:
    """Within the winning group, representative selection applies the deterministic hierarchy."""
    rows = [
        {
            "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|A",
            "Confidence": "High_Precision*",
            "Flag": "Low_Depth_Conserved_Motifs",
            "haplo_count": 100,
            "Depth_Score": 0.03,
            "POS": 60,
        },
        {
            "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|A",
            "Confidence": "High_Precision*",
            "Flag": "Not flagged",
            "haplo_count": 50,
            "Depth_Score": 0.02,
            "POS": 60,
        },
        {
            "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|A",
            "Confidence": "High_Precision*",
            "Flag": "Not flagged",
            "haplo_count": 80,
            "Depth_Score": 0.02,
            "POS": 60,
        },
    ]
    df = pd.DataFrame(rows)
    groups = group_candidates_by_identity(df)
    assert len(groups) == 1

    selected = select_group_representative_row(groups[0])
    assert len(selected) == 1
    # Unflagged preferred over flagged; among unflagged, haplo_count 80 preferred over 50
    assert selected.iloc[0]["Flag"] == "Not flagged"
    assert selected.iloc[0]["haplo_count"] == 80


def test_candidate_group_invariants_raise_on_invalid_construction() -> None:
    """Invalid parameters passed to IdentityCandidateGroup violate safety invariants."""
    with pytest.raises(ValueError, match="requires a non-empty identity_key"):
        IdentityCandidateGroup(
            identity_key="",
            identity=None,
            max_confidence_priority=1,
            context_fidelity=1,
            assembly_cardinality=1,
            peak_alternate_depth=10.0,
            peak_depth_score=0.01,
            candidates=pd.DataFrame([{"A": 1}]),
        )

    with pytest.raises(ValueError, match="requires at least one candidate row"):
        IdentityCandidateGroup(
            identity_key="valid-key",
            identity=None,
            max_confidence_priority=1,
            context_fidelity=1,
            assembly_cardinality=0,
            peak_alternate_depth=10.0,
            peak_depth_score=0.01,
            candidates=pd.DataFrame(),
        )

    with pytest.raises(ValueError, match="context_fidelity must be binary"):
        IdentityCandidateGroup(
            identity_key="valid-key",
            identity=None,
            max_confidence_priority=1,
            context_fidelity=2,
            assembly_cardinality=1,
            peak_alternate_depth=10.0,
            peak_depth_score=0.01,
            candidates=pd.DataFrame([{"A": 1}]),
        )
