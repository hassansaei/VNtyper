# tests/unit/test_kestrel_subthreshold_note.py

"""The subthreshold note reaches ``kestrel_result.tsv`` as a comment, and only there (#266).

The note is a ``##`` banner line, never a row and never a column: ``summary.parse_tsv``
routes ``#`` lines into ``comments`` and ``data`` never sees them, so no consumer that reads
the table can mistake a suppressed candidate for a call.

It is written from exactly one of ``process_kestrel_output``'s three empty-result branches
-- the one with a scored frame behind it. The other two run before anything is scored, so a
pre-result read there could only pick up a stale file from an earlier run into the same
output directory.
"""

import csv
from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts import kestrel_genotyping as kg
from vntyper.scripts import subthreshold as st
from vntyper.scripts.summary import parse_tsv

pytestmark = pytest.mark.unit

#: The 10 columns the negative placeholder has always carried. #266 must not touch them.
PLACEHOLDER_COLUMNS = [
    "Motif",
    "Variant",
    "POS",
    "REF",
    "ALT",
    "Motif_sequence",
    "Estimated_Depth_AlternateVariant",
    "Estimated_Depth_Variant_ActiveRegion",
    "Depth_Score",
    "Confidence",
]

NOTE = f"{st.NOTE_MARKER} 1 candidate variant in the pathogenic frame, best Depth_Score 0.0031."


def read_back(path: Path) -> tuple[list[str], list[dict]]:
    """Split a written result into its comment lines and its data rows."""
    comments: list[str] = []
    body: list[str] = []
    for line in path.read_text().splitlines():
        (comments if line.startswith("#") else body).append(line)
    rows = list(csv.DictReader(body, delimiter="\t")) if body else []
    return comments, rows


def pre_result(path: Path, rows: list[dict]) -> None:
    """Write a ``kestrel_pre_result.tsv`` whose rows default to the eligible state."""
    built = []
    for row in rows:
        base: dict[str, object] = dict.fromkeys((gate for gate in kg.FILTER_COLUMNS if gate != st.DEPTH_GATE), True)
        base[st.DEPTH_GATE] = False
        base |= {"POS": 60, "REF": "C", "ALT": "CC", "Depth_Score": 0.0031, "Motifs": "5-A"}
        base |= row
        built.append(base)
    pd.DataFrame(built).to_csv(path, sep="\t", index=False)


class TestOutputEmptyResult:
    def test_without_a_note_the_file_is_what_it_always_was(self, tmp_path: Path):
        kg.output_empty_result(str(tmp_path), ["## VNtyper Kestrel result"])

        comments, rows = read_back(tmp_path / "kestrel_result.tsv")

        assert comments == ["## VNtyper Kestrel result"]
        assert len(rows) == 1
        assert rows[0]["Confidence"] == "Negative"
        assert list(rows[0]) == PLACEHOLDER_COLUMNS

    def test_the_note_is_appended_as_a_banner_line(self, tmp_path: Path):
        kg.output_empty_result(str(tmp_path), ["## VNtyper Kestrel result"], note=NOTE)

        comments, rows = read_back(tmp_path / "kestrel_result.tsv")

        assert comments[-1] == f"## {NOTE}"
        assert len(rows) == 1, "the note must not add a row"
        assert list(rows[0]) == PLACEHOLDER_COLUMNS, "the note must not add a column"
        assert rows[0]["Confidence"] == "Negative", "the sample is still a negative"

    def test_the_note_follows_the_header_it_was_given(self, tmp_path: Path):
        header = ["## VNtyper Kestrel result", "## VNtyper Version: 9.9.9"]

        kg.output_empty_result(str(tmp_path), header, note=NOTE)

        comments, _ = read_back(tmp_path / "kestrel_result.tsv")

        assert comments == [*header, f"## {NOTE}"]

    @pytest.mark.parametrize("note", [None, ""])
    def test_an_absent_note_adds_no_line(self, tmp_path: Path, note):
        kg.output_empty_result(str(tmp_path), ["## VNtyper Kestrel result"], note=note)

        comments, _ = read_back(tmp_path / "kestrel_result.tsv")

        assert comments == ["## VNtyper Kestrel result"]

    def test_parse_tsv_routes_the_note_to_comments_and_not_to_data(self, tmp_path: Path):
        kg.output_empty_result(str(tmp_path), ["## VNtyper Kestrel result"], note=NOTE)

        parsed = parse_tsv(str(tmp_path / "kestrel_result.tsv"))

        assert st.find_note(parsed["comments"]) == NOTE
        assert len(parsed["data"]) == 1
        assert all(st.NOTE_MARKER not in str(value) for value in parsed["data"][0].values())


class TestTheNoteHelper:
    """``_subthreshold_note`` is the seam between the evidence file and the banner."""

    def config(self, **overrides) -> dict:
        settings = {"enabled": True, "template": "{marker} {events} {noun} at {best_depth_score}/{floor}"}
        settings.update(overrides)
        return {
            "confidence_assignment": {"depth_score_thresholds": {"low": 0.00469}},
            "subthreshold_note": settings,
        }

    def test_it_describes_an_eligible_row(self, tmp_path: Path):
        pre_result(tmp_path / "kestrel_pre_result.tsv", [{}])

        note = kg._subthreshold_note(str(tmp_path), self.config())

        assert note is not None
        assert note.startswith(st.NOTE_MARKER)
        assert "0.0031" in note

    def test_it_is_silent_when_nothing_is_eligible(self, tmp_path: Path):
        pre_result(tmp_path / "kestrel_pre_result.tsv", [{st.DEPTH_GATE: True}])

        assert kg._subthreshold_note(str(tmp_path), self.config()) is None

    def test_it_is_silent_for_an_artifact_row_however_deep(self, tmp_path: Path):
        """A row killed by the artifact gate is strong signal deliberately discarded (#174),
        not weak signal, and must never be described as subthreshold."""
        pre_result(tmp_path / "kestrel_pre_result.tsv", [{"flag_filter_pass": False, "Depth_Score": 0.9}])

        assert kg._subthreshold_note(str(tmp_path), self.config()) is None

    def test_disabling_it_restores_the_previous_output_exactly(self, tmp_path: Path):
        pre_result(tmp_path / "kestrel_pre_result.tsv", [{}])

        assert kg._subthreshold_note(str(tmp_path), self.config(enabled=False)) is None

    def test_a_missing_template_yields_no_note(self, tmp_path: Path):
        pre_result(tmp_path / "kestrel_pre_result.tsv", [{}])
        config = self.config()
        del config["subthreshold_note"]["template"]

        assert kg._subthreshold_note(str(tmp_path), config) is None

    def test_a_config_without_the_floor_loses_the_note_not_the_run(self, tmp_path: Path):
        """A partial config raises in `calculate_depth_score_and_assign_confidence`, which is
        right for a gate. This is an annotation, so it degrades instead."""
        pre_result(tmp_path / "kestrel_pre_result.tsv", [{}])
        config = self.config()
        config["confidence_assignment"] = {}

        assert kg._subthreshold_note(str(tmp_path), config) is None

    def test_a_missing_pre_result_yields_no_note(self, tmp_path: Path):
        assert kg._subthreshold_note(str(tmp_path), self.config()) is None

    def test_the_shipped_config_renders_a_note(self):
        """The template that actually ships must be renderable; a `{typo}` in it would
        silently disable the feature."""
        signal = st.SubthresholdSignal(events=1, rows=1, best_depth_score=0.0031, floor=0.00469)
        template = kg.kestrel_config["subthreshold_note"]["template"]

        note = st.format_note(signal, template)

        assert note is not None
        assert note.startswith(st.NOTE_MARKER)
        assert "0.0031" in note and "0.00469" in note
        assert "NOT a call" in note

    def test_the_shipped_config_enables_the_feature(self):
        assert kg.kestrel_config["subthreshold_note"]["enabled"] is True


class TestFilterColumns:
    def test_the_constant_is_the_six_gates_in_order(self):
        assert kg.FILTER_COLUMNS == (
            "is_frameshift",
            "is_valid_frameshift",
            "depth_confidence_pass",
            "alt_filter_pass",
            "motif_filter_pass",
            "flag_filter_pass",
        )

    def test_the_depth_gate_the_note_keys_on_is_one_of_them(self):
        assert st.DEPTH_GATE in kg.FILTER_COLUMNS

    def test_filter_final_dataframe_still_requires_every_one_of_them(self, tmp_path: Path):
        """A gate dropped from a non-empty frame is an error, not a permit (#185). The
        extraction to a constant must not have loosened that."""
        frame = pd.DataFrame([dict.fromkeys(kg.FILTER_COLUMNS, True) | {"POS": 60, "Depth_Score": 0.9}]).drop(
            columns=["flag_filter_pass"]
        )

        with pytest.raises(ValueError, match="flag_filter_pass"):
            kg.filter_final_dataframe(frame, str(tmp_path))
