"""Tests for fail-closed idxstats parsing and scan selection."""

import pytest

from vntyper.scripts.idxstats_parsing import SCAN_INDEXED, SCAN_STREAM, choose_scan, parse_idxstats

pytestmark = pytest.mark.unit

GOOD = "chr1\t20000\t600\t50\n*\t0\t0\t80\n"
CLEAN = "chr1\t20000\t600\t0\n*\t0\t0\t80\n"


class TestParse:
    def test_it_returns_placed_and_unplaced_counts(self):
        assert parse_idxstats(GOOD) == (50, 80)

    def test_a_clean_file_has_no_placed_unmapped_reads(self):
        assert parse_idxstats(CLEAN) == (0, 80)

    @pytest.mark.parametrize(
        "bad",
        [
            "",
            "   \n",
            "chr1\t20000\t600\n",
            "chr1\t20000\t600\t50\textra\n",
            "chr1\t20000\tsix\t50\n*\t0\t0\t80\n",
            "chr1\t20000\t600\t-1\n*\t0\t0\t80\n",
            "chr1\t20000\t600\t50\n",
            "*\t0\t0\t1\n*\t0\t0\t2\n",
        ],
    )
    def test_anything_malformed_is_rejected_rather_than_guessed(self, bad):
        assert parse_idxstats(bad) is None

    @pytest.mark.parametrize(
        ("column", "spelling"),
        [
            (1, " 20000"),
            (2, " 600"),
            (3, " 50"),
            (1, "20000 "),
            (2, "600 "),
            (3, "50 "),
            (1, "+20000"),
            (2, "+600"),
            (3, "+50"),
            (1, "-0"),
            (2, "-0"),
            (3, "-0"),
        ],
    )
    def test_numeric_fields_must_be_ascii_decimal_tokens(self, column, spelling):
        fields = ["chr1", "20000", "600", "50"]
        fields[column] = spelling
        text = "\t".join(fields) + "\n*\t0\t0\t80\n"
        assert parse_idxstats(text) is None

    def test_an_oversized_decimal_count_is_rejected_without_raising(self):
        oversized_count = "9" * 5000
        text = f"chr1\t{oversized_count}\t600\t50\n*\t0\t0\t80\n"
        assert parse_idxstats(text) is None


class TestChooseScan:
    def test_zero_placed_counts_do_not_authorize_an_incomplete_literal_star_fetch(self):
        """The CRAM index can expose fewer literal-star records than idxstats reports."""
        scan, reason = choose_scan(
            "auto",
            CLEAN,
            exit_ok=True,
            indexed_count_text="20\n",
            indexed_count_exit_ok=True,
        )

        assert scan == SCAN_STREAM
        assert reason == "indexed '*' count 20 differs from idxstats unplaced count 80; using lossless stream scan"

    def test_forced_indexed_rejects_an_incomplete_literal_star_fetch(self):
        """An explicit indexed policy cannot override evidence that its consumer is lossy."""
        with pytest.raises(ValueError, match=r"indexed '\*' count 20 differs from idxstats unplaced count 80"):
            choose_scan(
                "indexed",
                CLEAN,
                exit_ok=True,
                indexed_count_text="20\n",
                indexed_count_exit_ok=True,
            )

    def test_matching_literal_star_count_authorizes_indexed_recovery(self):
        """Indexed recovery is safe only when the exact consumer count matches idxstats."""
        assert choose_scan(
            "auto",
            CLEAN,
            exit_ok=True,
            indexed_count_text="80\n",
            indexed_count_exit_ok=True,
        ) == (SCAN_INDEXED, "idxstats and indexed '*' count agree on 80 unplaced reads")

    @pytest.mark.parametrize(("count_text", "count_exit_ok"), [("not-a-count\n", True), (None, False)])
    def test_missing_or_malformed_literal_star_evidence_fails_closed(self, count_text, count_exit_ok):
        """A count command failure cannot authorize the consumer it was meant to prove."""
        assert (
            choose_scan(
                "auto",
                CLEAN,
                exit_ok=True,
                indexed_count_text=count_text,
                indexed_count_exit_ok=count_exit_ok,
            )[0]
            == SCAN_STREAM
        )

    @pytest.mark.parametrize(("count_text", "count_exit_ok"), [("not-a-count\n", True), (None, False)])
    def test_forced_indexed_rejects_missing_or_malformed_literal_star_evidence(self, count_text, count_exit_ok):
        """A forced mode must raise when its exact-consumer evidence is unusable."""
        with pytest.raises(ValueError, match=r"indexed '\*' count"):
            choose_scan(
                "indexed",
                CLEAN,
                exit_ok=True,
                indexed_count_text=count_text,
                indexed_count_exit_ok=count_exit_ok,
            )

    def test_oversized_literal_star_count_fails_closed_without_an_unbounded_diagnostic(self):
        """Hostile captured output cannot escape count parsing or the message bound."""
        scan, reason = choose_scan(
            "auto",
            CLEAN,
            exit_ok=True,
            indexed_count_text="9" * 5000,
            indexed_count_exit_ok=True,
        )

        assert scan == SCAN_STREAM
        assert len(reason) <= 256

    def test_auto_picks_indexed_only_when_nothing_would_be_lost(self):
        assert (
            choose_scan("auto", CLEAN, exit_ok=True, indexed_count_text="80\n", indexed_count_exit_ok=True)[0]
            == SCAN_INDEXED
        )

    def test_auto_falls_back_to_stream_when_placed_unmapped_reads_exist(self):
        scan, reason = choose_scan("auto", GOOD, exit_ok=True)
        assert scan == SCAN_STREAM and "50" in reason

    def test_auto_falls_back_to_stream_on_a_malformed_table(self):
        assert choose_scan("auto", "garbage", exit_ok=True)[0] == SCAN_STREAM

    def test_malformed_reason_names_the_exact_offending_line(self):
        malformed = "chr1\t20000\tsix\t50\n*\t0\t0\t80\n"

        scan, reason = choose_scan("auto", malformed, exit_ok=True)

        assert scan == SCAN_STREAM
        assert reason == (
            "idxstats output is malformed at line 1: 'chr1\\t20000\\tsix\\t50'; using lossless stream scan"
        )

    def test_oversized_malformed_line_is_escaped_and_bounded(self):
        oversized_count = "9" * 5000
        malformed = f"chr1\t{oversized_count}\t600\t50\n*\t0\t0\t80\n"

        scan, reason = choose_scan("auto", malformed, exit_ok=True)

        assert scan == SCAN_STREAM
        assert reason.startswith("idxstats output is malformed at line 1: 'chr1\\t999")
        assert "...<truncated>" in reason
        assert "\\t600\\t50" not in reason
        assert len(reason) <= 256

    def test_empty_output_reason_names_the_missing_evidence(self):
        assert choose_scan("auto", "", exit_ok=True) == (
            SCAN_STREAM,
            "idxstats output is empty; using lossless stream scan",
        )

    def test_missing_star_reason_names_the_required_terminal_row(self):
        without_star = "chr1\t20000\t600\t0\n"

        assert choose_scan("auto", without_star, exit_ok=True) == (
            SCAN_STREAM,
            "idxstats output is missing its terminal '*' row; using lossless stream scan",
        )

    def test_auto_falls_back_to_stream_when_idxstats_itself_failed(self):
        assert choose_scan("auto", CLEAN, exit_ok=False)[0] == SCAN_STREAM

    def test_forcing_indexed_where_reads_would_be_lost_raises_rather_than_dropping(self):
        with pytest.raises(ValueError, match="50"):
            choose_scan("indexed", GOOD, exit_ok=True)

    def test_stream_is_always_allowed_because_it_is_never_lossy(self):
        assert choose_scan("stream", GOOD, exit_ok=True)[0] == SCAN_STREAM

    def test_an_unknown_scan_mode_is_rejected(self):
        with pytest.raises(ValueError, match="bogus"):
            choose_scan("bogus", CLEAN, exit_ok=True)

    def test_auto_falls_back_to_stream_when_idxstats_output_is_missing(self):
        assert choose_scan("auto", None, exit_ok=True)[0] == SCAN_STREAM

    @pytest.mark.parametrize("idxstats_text", [None, "garbage"])
    def test_forced_indexed_falls_back_to_stream_without_valid_evidence(self, idxstats_text):
        assert choose_scan("indexed", idxstats_text, exit_ok=True)[0] == SCAN_STREAM
