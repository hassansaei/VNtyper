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


class TestChooseScan:
    def test_auto_picks_indexed_only_when_nothing_would_be_lost(self):
        assert choose_scan("auto", CLEAN, exit_ok=True)[0] == SCAN_INDEXED

    def test_auto_falls_back_to_stream_when_placed_unmapped_reads_exist(self):
        scan, reason = choose_scan("auto", GOOD, exit_ok=True)
        assert scan == SCAN_STREAM and "50" in reason

    def test_auto_falls_back_to_stream_on_a_malformed_table(self):
        assert choose_scan("auto", "garbage", exit_ok=True)[0] == SCAN_STREAM

    def test_auto_falls_back_to_stream_when_idxstats_itself_failed(self):
        assert choose_scan("auto", CLEAN, exit_ok=False)[0] == SCAN_STREAM

    def test_forcing_indexed_where_reads_would_be_lost_raises_rather_than_dropping(self):
        with pytest.raises(ValueError, match="50"):
            choose_scan("indexed", GOOD, exit_ok=True)

    def test_stream_is_always_allowed_because_it_is_never_lossy(self):
        assert choose_scan("stream", GOOD, exit_ok=True)[0] == SCAN_STREAM
