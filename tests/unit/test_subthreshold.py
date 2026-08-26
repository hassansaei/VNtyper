# tests/unit/test_subthreshold.py

"""The below-reporting-floor signal (#266).

A subthreshold candidate is a row that fails the depth gate and *nothing else*. That
strictness is the whole design: a row killed by ``flag_filter_pass`` carries a declared
artifact flag (#174) and its ``Depth_Score`` may be excellent, so calling it "subthreshold"
would say *weak signal* where the truth is *strong signal, deliberately discarded*.

The second half of the design is that a gate verdict must be **explicit**. The frame is
read back from a TSV, where one missing value turns a whole boolean column into strings,
and ``bool("False")`` is ``True`` -- so a cast would turn a failing gate into a passing one.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from vntyper.scripts import subthreshold as st

pytestmark = pytest.mark.unit

#: The six gates ``filter_final_dataframe`` ANDs, in source order.
GATES = (
    "is_frameshift",
    "is_valid_frameshift",
    "depth_confidence_pass",
    "alt_filter_pass",
    "motif_filter_pass",
    "flag_filter_pass",
)

STRUCTURAL = tuple(gate for gate in GATES if gate != st.DEPTH_GATE)

FLOOR = 0.00469

#: Values a gate cell can carry that assert nothing. ``pd.NA`` and ``NaN`` arrive from a
#: TSV round trip; the strings from a hand-edited or foreign file.
UNREADABLE = [None, float("nan"), "", "maybe", 2, "NA", pd.NA]


def frame(rows: list[dict]) -> pd.DataFrame:
    """Build a pre-result-shaped frame, defaulting every gate to the eligible state."""
    built = []
    for row in rows:
        base: dict[str, object] = dict.fromkeys(STRUCTURAL, True)
        base[st.DEPTH_GATE] = False
        base |= {"POS": 60, "REF": "C", "ALT": "CC", "Depth_Score": 0.001, "Motifs": "5-A"}
        base |= row
        built.append(base)
    return pd.DataFrame(built)


class TestEligibility:
    def test_a_row_failing_only_the_depth_gate_is_subthreshold(self):
        signal = st.detect(frame([{}]), GATES, FLOOR)

        assert signal is not None
        assert signal.events == 1
        assert signal.rows == 1
        assert signal.best_depth_score == pytest.approx(0.001)
        assert signal.floor == pytest.approx(FLOOR)

    def test_a_row_that_passes_the_depth_gate_is_not_subthreshold(self):
        assert st.detect(frame([{st.DEPTH_GATE: True}]), GATES, FLOOR) is None

    @pytest.mark.parametrize("gate", STRUCTURAL)
    def test_a_row_failing_any_structural_gate_is_not_subthreshold(self, gate):
        """The decisive case is ``flag_filter_pass``: a declared artifact is not weak
        signal, whatever its depth score says."""
        assert st.detect(frame([{gate: False, "Depth_Score": 0.9}]), GATES, FLOOR) is None

    def test_an_eligible_row_beside_an_artifact_row_reports_only_the_eligible_one(self):
        signal = st.detect(
            frame(
                [
                    {"Depth_Score": 0.002},
                    {"flag_filter_pass": False, "Depth_Score": 0.9, "POS": 61},
                ]
            ),
            GATES,
            FLOOR,
        )

        assert signal is not None
        assert signal.rows == 1
        assert signal.best_depth_score == pytest.approx(0.002)


class TestCounting:
    def test_events_counts_distinct_pos_ref_alt_not_rows(self):
        """The same event against several motif contexts is one event."""
        signal = st.detect(
            frame(
                [
                    {"Motifs": "5-A", "Depth_Score": 0.001},
                    {"Motifs": "6-B", "Depth_Score": 0.002},
                    {"POS": 61, "Motifs": "5-A", "Depth_Score": 0.003},
                ]
            ),
            GATES,
            FLOOR,
        )

        assert signal is not None
        assert signal.rows == 3
        assert signal.events == 2
        assert signal.best_depth_score == pytest.approx(0.003)

    def test_booleans_that_survived_a_tsv_round_trip_as_strings_are_honoured(self):
        raw = frame([{}])
        for gate in GATES:
            raw[gate] = raw[gate].astype(str)

        signal = st.detect(raw, GATES, FLOOR)

        assert signal is not None
        assert signal.rows == 1

    def test_a_string_false_does_not_pass_a_structural_gate(self):
        raw = frame([{}])
        raw["flag_filter_pass"] = "False"

        assert st.detect(raw, GATES, FLOOR) is None

    def test_a_non_numeric_depth_score_is_excluded(self):
        assert st.detect(frame([{"Depth_Score": "n/a"}]), GATES, FLOOR) is None

    def test_a_non_numeric_row_does_not_hide_a_numeric_one(self):
        signal = st.detect(frame([{"Depth_Score": "n/a"}, {"POS": 61, "Depth_Score": 0.004}]), GATES, FLOOR)

        assert signal is not None
        assert signal.rows == 1
        assert signal.best_depth_score == pytest.approx(0.004)


class TestVerdictsAreThreeValued:
    """ "Not True" and "False" are different claims, and eligibility needs both."""

    @pytest.mark.parametrize("value", UNREADABLE)
    def test_an_unrecognisable_structural_verdict_disqualifies_the_row(self, value):
        raw = frame([{}])
        raw["motif_filter_pass"] = pd.Series([value], dtype=object)

        assert st.detect(raw, GATES, FLOOR) is None

    @pytest.mark.parametrize("value", UNREADABLE)
    def test_an_unrecognisable_depth_verdict_disqualifies_the_row(self, value):
        """The dangerous direction: reading "not True" as "failed the depth gate" would
        make a row with no recorded verdict look like a suppressed candidate."""
        raw = frame([{}])
        raw[st.DEPTH_GATE] = pd.Series([value], dtype=object)

        assert st.detect(raw, GATES, FLOOR) is None

    def test_numpy_booleans_are_read_as_booleans(self):
        """``np.bool_`` does not subclass ``bool``; an isinstance test alone misses it."""
        raw = frame([{}])
        raw["motif_filter_pass"] = pd.Series([np.True_], dtype=object)
        raw[st.DEPTH_GATE] = pd.Series([np.False_], dtype=object)

        assert st.detect(raw, GATES, FLOOR) is not None

    def test_zero_and_one_are_read_as_the_verdicts_they_encode(self):
        raw = frame([{}])
        for gate in STRUCTURAL:
            raw[gate] = 1
        raw[st.DEPTH_GATE] = 0

        assert st.detect(raw, GATES, FLOOR) is not None

    def test_one_unreadable_row_does_not_hide_a_readable_one(self):
        raw = frame([{"Depth_Score": 0.002}, {"POS": 61, "Depth_Score": 0.004}])
        raw["motif_filter_pass"] = pd.Series([True, "?"], dtype=object)

        signal = st.detect(raw, GATES, FLOOR)

        assert signal is not None
        assert signal.rows == 1
        assert signal.best_depth_score == pytest.approx(0.002)

    def test_verdict_reads_every_shipped_spelling(self):
        """Pinned directly, because every mask above is built from this one function."""
        assert st._verdict(True) is True
        assert st._verdict(False) is False
        assert st._verdict("True") is True
        assert st._verdict("false") is False
        assert st._verdict(" TRUE ") is True
        assert st._verdict(np.True_) is True
        assert st._verdict(np.False_) is False
        assert st._verdict(1) is True
        assert st._verdict(0) is False
        for value in UNREADABLE:
            assert st._verdict(value) is None, value


class TestDegenerateInput:
    def test_an_empty_frame_yields_no_signal(self):
        assert st.detect(pd.DataFrame(), GATES, FLOOR) is None

    def test_a_none_frame_yields_no_signal(self):
        assert st.detect(None, GATES, FLOOR) is None

    def test_a_missing_gate_column_yields_no_signal_and_does_not_raise(self):
        raw = frame([{}]).drop(columns=["flag_filter_pass"])

        assert st.detect(raw, GATES, FLOOR) is None

    def test_a_missing_depth_score_column_yields_no_signal(self):
        raw = frame([{}]).drop(columns=["Depth_Score"])

        assert st.detect(raw, GATES, FLOOR) is None

    def test_a_gate_list_without_the_depth_gate_yields_no_signal(self):
        """Guard the guard: with the depth gate absent from the list there is no
        magnitude judgement left, so nothing may be called subthreshold."""
        assert st.detect(frame([{}]), STRUCTURAL, FLOOR) is None


class TestFileEntry:
    def test_it_reads_a_written_pre_result(self, tmp_path: Path):
        path = tmp_path / "kestrel_pre_result.tsv"
        frame([{}]).to_csv(path, sep="\t", index=False)

        signal = st.detect_from_file(path, GATES, FLOOR)

        assert signal is not None
        assert signal.events == 1

    def test_it_reads_a_pre_result_whose_gate_column_carries_a_missing_value(self, tmp_path: Path):
        """The round trip that turns booleans into strings: one NA in any gate column."""
        raw = frame([{"Depth_Score": 0.002}, {"POS": 61, "Depth_Score": 0.004}])
        raw["flag_filter_pass"] = pd.Series([True, None], dtype=object)
        path = tmp_path / "kestrel_pre_result.tsv"
        raw.to_csv(path, sep="\t", index=False)

        signal = st.detect_from_file(path, GATES, FLOOR)

        assert signal is not None
        assert signal.rows == 1, "the row whose artifact gate says nothing is not eligible"
        assert signal.best_depth_score == pytest.approx(0.002)

    def test_a_missing_file_yields_no_signal(self, tmp_path: Path):
        assert st.detect_from_file(tmp_path / "absent.tsv", GATES, FLOOR) is None

    def test_a_directory_in_place_of_the_file_yields_no_signal(self, tmp_path: Path):
        assert st.detect_from_file(tmp_path, GATES, FLOOR) is None

    @pytest.mark.parametrize(
        ("name", "payload"),
        [
            # A ragged file: pandas.errors.ParserError, which subclasses ValueError.
            ("ragged", b"a\tb\tc\n1\t2\n3\t4\t5\t6\t7\n"),
            # Not UTF-8: UnicodeDecodeError, which also subclasses ValueError.
            ("undecodable", b"a\tb\n\xff\xfe\t2\n"),
        ],
    )
    def test_an_unparseable_file_yields_no_signal(self, tmp_path: Path, name, payload):
        """These two really do raise. A file of NUL bytes does **not** -- pandas reads it
        as a one-column empty frame -- so testing with that would exercise the empty-frame
        path and leave the read handler uncovered while looking like it did not."""
        path = tmp_path / f"{name}.tsv"
        path.write_bytes(payload)

        assert st.detect_from_file(path, GATES, FLOOR) is None

    def test_a_file_of_nul_bytes_is_read_as_empty_rather_than_raising(self, tmp_path: Path):
        """Pinned so the parametrisation above is not quietly reverted to this input."""
        path = tmp_path / "nulls.tsv"
        path.write_bytes(b"\x00\x01\x02")

        assert st.detect_from_file(path, GATES, FLOOR) is None

    def test_an_empty_file_yields_no_signal(self, tmp_path: Path):
        path = tmp_path / "kestrel_pre_result.tsv"
        path.write_text("")

        assert st.detect_from_file(path, GATES, FLOOR) is None


class TestRendering:
    TEMPLATE = (
        "{marker} {events} {noun} in the pathogenic frame identified below the "
        "reporting floor (best Depth_Score {best_depth_score}, floor {floor}); "
        "filtered out, not a call."
    )

    def signal(self, **overrides) -> st.SubthresholdSignal:
        values = {"events": 1, "rows": 1, "best_depth_score": 0.00312, "floor": 0.00469}
        values.update(overrides)
        return st.SubthresholdSignal(**values)

    def test_the_note_begins_with_the_marker(self):
        assert st.format_note(self.signal(), self.TEMPLATE).startswith(st.NOTE_MARKER)

    def test_the_note_carries_the_score_the_floor_and_the_count(self):
        note = st.format_note(self.signal(rows=4), self.TEMPLATE)

        assert "0.00312" in note
        assert "0.00469" in note
        assert " 1 candidate variant " in note

    def test_the_noun_agrees_in_number(self):
        one = st.format_note(self.signal(), self.TEMPLATE)
        many = st.format_note(self.signal(events=3, rows=9), self.TEMPLATE)

        assert " 1 candidate variant " in one
        assert " 3 candidate variants " in many

    def test_the_note_never_carries_a_newline_or_a_tab(self):
        """It is written as one ``##`` line of a TSV; either would split it."""
        note = st.format_note(self.signal(), "{marker}\n{events}\t{noun}")

        assert "\n" not in note
        assert "\t" not in note

    def test_a_template_naming_an_unknown_field_yields_no_note(self):
        assert st.format_note(self.signal(), "{marker} {not_a_field}") is None

    def test_a_malformed_template_yields_no_note(self):
        assert st.format_note(self.signal(), "{marker} {") is None

    def test_the_rendered_note_is_findable(self):
        """The writer's output and the reader's matcher agree, which is the contract
        between `kestrel_genotyping` and `generate_report`."""
        note = st.format_note(self.signal(), self.TEMPLATE)

        assert st.find_note([f"## {note}"]) == note


class TestFindNote:
    def test_it_finds_the_marked_line_among_the_banner(self):
        comments = [
            "VNtyper Kestrel result",
            "VNtyper Version: 2.0.22",
            f"{st.NOTE_MARKER} 1 candidate variant ...",
        ]

        assert st.find_note(comments) == comments[2]

    def test_it_returns_none_when_no_line_is_marked(self):
        assert st.find_note(["VNtyper Kestrel result", "Analysis date: 2026-08-26"]) is None

    def test_it_ignores_leading_hashes_and_whitespace(self):
        assert st.find_note([f"##  {st.NOTE_MARKER} x"]) == f"{st.NOTE_MARKER} x"

    def test_it_returns_none_for_no_comments(self):
        assert st.find_note([]) is None

    def test_it_returns_none_for_none(self):
        assert st.find_note(None) is None

    def test_it_returns_the_first_marked_line(self):
        first = f"{st.NOTE_MARKER} first"
        assert st.find_note([first, f"{st.NOTE_MARKER} second"]) == first
