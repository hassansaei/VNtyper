"""The emitted motif sequence must be the 60 bp half named by ``Motif``.

Pair records are ``seq(R) + seq(L)``: for a record named ``<L>-<R>``,
positions 1-60 hold the R-named motif and 61-120 the L-named motif.
"""

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from vntyper.scripts.motif_processing import motif_correction_and_annotation
from vntyper.scripts.report_formatting import COLUMN_HELP

pytestmark = pytest.mark.unit

CONFIG_PATH = Path("vntyper/scripts/kestrel_config.json")

SIXTY = {
    "X": "A" * 60,
    "3": "G" * 60,
    "1": "T" * 60,
    "2": "C" * 60,
    "8": "AC" * 30,
    "7": "AG" * 30,
}


def _shipped_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def _merged_motifs() -> pd.DataFrame:
    return pd.DataFrame({"Motif": list(SIXTY), "Motif_sequence": list(SIXTY.values())})


def _frame(motifs_id: str, pos: int, ref: str, alt: str) -> pd.DataFrame:
    left_name, right_name = motifs_id.split("-")
    pair = SIXTY[right_name] + SIXTY[left_name]
    return pd.DataFrame(
        {
            "Motifs": [motifs_id],
            "Variant": ["Insertion"],
            "POS": [pos],
            "REF": [ref],
            "ALT": [alt],
            "Motif_sequence": [pair],
            "Estimated_Depth_AlternateVariant": [50],
            "Estimated_Depth_Variant_ActiveRegion": [5000],
            "Depth_Score": [0.010],
            "Confidence": ["High_Precision"],
        }
    )


def test_a_right_half_row_emits_the_60bp_half_named_before_the_dash() -> None:
    out = motif_correction_and_annotation(_frame("X-3", 67, "G", "GC"), _merged_motifs(), _shipped_config())

    assert out["motif_filter_pass"].tolist() == [True]
    assert out["Motif"].tolist() == ["X"]
    assert out["Motif_sequence"].tolist() == [SIXTY["X"]]
    assert len(out["Motif_sequence"].iloc[0]) == 60


def test_a_left_half_row_emits_the_60bp_half_named_after_the_dash() -> None:
    out = motif_correction_and_annotation(_frame("1-2", 54, "C", "CC"), _merged_motifs(), _shipped_config())

    assert out["motif_filter_pass"].tolist() == [True]
    assert out["Motif"].tolist() == ["2"]
    assert out["Motif_sequence"].tolist() == [SIXTY["2"]]
    assert len(out["Motif_sequence"].iloc[0]) == 60


def test_a_passing_row_without_an_input_pair_column_still_emits_the_selected_half() -> None:
    frame = _frame("X-3", 67, "G", "GC").drop(columns=["Motif_sequence"])

    out = motif_correction_and_annotation(frame, _merged_motifs(), _shipped_config())

    assert out["Motif_sequence"].tolist() == [SIXTY["X"]]


def test_a_failing_row_keeps_its_input_pair_sequence() -> None:
    """The new sequence copy-back reaches final-passing rows only."""
    frame = _frame("8-7", 70, "C", "CGGCA")
    pair = frame["Motif_sequence"].iloc[0]

    out = motif_correction_and_annotation(frame, _merged_motifs(), _shipped_config())

    assert out["motif_filter_pass"].tolist() == [False]
    assert out["Motif_sequence"].tolist() == [pair]


def test_an_annotation_surviving_nonframeshift_row_keeps_legacy_annotations_but_not_the_new_half() -> None:
    """The final gate applies to the new sequence copy-back, not legacy annotations."""
    frame = _frame("X-3", 67, "G", "GC")
    frame["is_valid_frameshift"] = False
    pair = frame["Motif_sequence"].iloc[0]

    out = motif_correction_and_annotation(frame, _merged_motifs(), _shipped_config())

    assert out["motif_filter_pass"].tolist() == [False]
    assert out["Motif_sequence"].tolist() == [pair]
    assert out["Motif"].tolist() == ["X"]
    assert out["Motif_fasta"].tolist() == ["X-3"]
    assert out["POS_fasta"].tolist() == [67]


def test_the_report_explains_that_the_emitted_value_is_the_named_60bp_half() -> None:
    help_text = COLUMN_HELP["Motif Sequence"]

    assert "60 bp" in help_text
    assert "Motif column" in help_text
    assert "121 bp" not in help_text
