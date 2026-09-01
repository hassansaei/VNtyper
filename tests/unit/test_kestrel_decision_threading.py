"""Resolved Kestrel selection fields reach their exact decision points."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

from vntyper.scripts.kestrel_genotyping import select_single_best_variant
from vntyper.scripts.run_configuration import resolve_run_configuration
from vntyper.scripts.scoring import extract_frameshifts

pytestmark = pytest.mark.unit


def test_explicit_confidence_priority_changes_the_selected_row() -> None:
    raw_selection = resolve_run_configuration().kestrel["selection"]
    assert isinstance(raw_selection, Mapping)
    selection = dict(raw_selection)
    selection["confidence_priority"] = {
        "High_Precision*": 0,
        "High_Precision": 1,
        "Low_Precision": 9,
        "Negative": 0,
    }
    frame = pd.DataFrame(
        {
            "Confidence": ["High_Precision*", "Low_Precision"],
            "Flag": ["Not flagged", "Not flagged"],
            "Depth_Score": [0.9, 0.1],
            "haplo_count": [20, 1],
            "POS": [1, 2],
        }
    )

    result = select_single_best_variant(frame, selection=selection)

    assert result.iloc[0]["Confidence"] == "Low_Precision"


def test_explicit_frameshift_projection_controls_both_arms() -> None:
    frame = pd.DataFrame({"direction": [1, -1], "frameshift_amount": [2, 1]})
    frameshift = {"modulus": 3, "insertion_remainder": 2, "deletion_remainder": 1}

    result = extract_frameshifts(frame, frameshift=frameshift)

    assert result["is_valid_frameshift"].tolist() == [True, True]
