#!/usr/bin/env python3
"""
tests/unit/test_summary_parsers.py

Unit tests for the TSV/CSV parsers in vntyper.scripts.summary.

These parsers feed pipeline_summary.json, which generate_report.py and
cohort_summary.py then read, so a silently mis-parsed row propagates into the report a
clinician sees. The ragged-row behaviour is therefore pinned explicitly: a row whose
field count disagrees with the header is skipped and logged, never truncated into a
plausible-looking partial record, and never allowed to discard the rest of the file.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.summary import md5sum, parse_csv, parse_json_file, parse_tsv

pytestmark = pytest.mark.unit
SUMMARY_LOGGER = "vntyper.scripts.summary"


def _write(path: Path, text: str) -> str:
    """Write text to a file and return its path as a string.

    Args:
        path: Destination path.
        text: File contents.

    Returns:
        str: The path, as the parsers expect.
    """
    path.write_text(text, encoding="utf-8")
    return str(path)


# --- well-formed input -------------------------------------------------------


def test_parse_tsv_reads_comments_and_rows(tmp_path: Path) -> None:
    """Comment lines become comments; the first non-comment line is the header."""
    path = _write(
        tmp_path / "ok.tsv",
        "# VNtyper Kestrel result\n# Version: 2.0.4\nMotifs\tPOS\tREF\nD-C\t67\tG\nD-A\t12\tT\n",
    )
    result = parse_tsv(path)

    assert result["comments"] == ["VNtyper Kestrel result", "Version: 2.0.4"]
    assert result["data"] == [
        {"Motifs": "D-C", "POS": "67", "REF": "G"},
        {"Motifs": "D-A", "POS": "12", "REF": "T"},
    ]


def test_parse_csv_reads_comments_and_rows(tmp_path: Path) -> None:
    """The CSV parser treats a row whose first cell starts with '#' as a comment."""
    path = _write(tmp_path / "ok.csv", "# header comment\nMotifs,POS\nD-C,67\n")
    result = parse_csv(path)

    assert result["comments"] == ["header comment"]
    assert result["data"] == [{"Motifs": "D-C", "POS": "67"}]


def test_parse_tsv_skips_blank_lines(tmp_path: Path) -> None:
    """Blank lines are ignored rather than becoming empty records."""
    path = _write(tmp_path / "blank.tsv", "A\tB\n\n1\t2\n\n")
    assert parse_tsv(path)["data"] == [{"A": "1", "B": "2"}]


def test_parse_tsv_keeps_a_row_whose_last_column_is_empty(tmp_path: Path) -> None:
    """A trailing empty field is a value, not absent whitespace.

    Stripping the line removes the final tab along with it, so the row arrives one
    field short of the header and is discarded as ragged -- silently, because the row
    itself is perfectly well formed. Every nullable column is written as the empty
    string by contract, so whenever one of them is last this drops real data: it took
    out the whole Kestrel result whenever adVNTR had not run, leaving the summary with
    no Kestrel records and the cross-match step raising on a run that had succeeded.
    """
    path = _write(tmp_path / "trailing.tsv", "A\tB\tC\n1\t2\t\n")

    assert parse_tsv(path)["data"] == [{"A": "1", "B": "2", "C": ""}]


def test_parse_tsv_keeps_a_row_that_is_empty_but_for_its_first_column(tmp_path: Path) -> None:
    """The degenerate case of the same bug: every field after the first is empty."""
    path = _write(tmp_path / "mostly-empty.tsv", "A\tB\tC\n1\t\t\n")

    assert parse_tsv(path)["data"] == [{"A": "1", "B": "", "C": ""}]


def test_parse_tsv_decodes_writer_quoted_json_cells(tmp_path: Path) -> None:
    """TSV quoting added by pandas must not become part of a summary cell."""
    path = _write(
        tmp_path / "quoted.tsv",
        'Name\tMetadata\ncall\t"{""source"":""kestrel"",""values"":[67,""G""]}"\n',
    )

    assert parse_tsv(path)["data"] == [{"Name": "call", "Metadata": '{"source":"kestrel","values":[67,"G"]}'}]


# --- ragged rows: the behaviour this module pins -----------------------------


@pytest.mark.parametrize(
    ("row", "found"),
    [("1\t2", 2), ("1\t2\t3\t4", 4)],
    ids=["too-few-fields", "too-many-fields"],
)
def test_parse_tsv_skips_ragged_row_and_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, row: str, found: int
) -> None:
    """A ragged row is dropped with a warning naming the file and line.

    Truncating it instead would produce a record that looks valid but is missing
    columns - corruption that reaches the report unnoticed.

    Args:
        tmp_path: Pytest temporary directory.
        caplog: Log capture fixture.
        row: The malformed row.
        found: Number of fields that row actually has.
    """
    path = _write(tmp_path / "ragged.tsv", f"A\tB\tC\n{row}\n7\t8\t9\n")

    with caplog.at_level(logging.WARNING):
        result = parse_tsv(path)

    # The malformed row is gone; the good row after it survives.
    assert result["data"] == [{"A": "7", "B": "8", "C": "9"}]
    assert f"expected 3 fields, found {found}" in caplog.text
    assert "ragged.tsv:2" in caplog.text


def test_parse_csv_skips_ragged_row_and_warns(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """The CSV parser applies the same rule as the TSV parser."""
    path = _write(tmp_path / "ragged.csv", "A,B,C\n1,2\n7,8,9\n")

    with caplog.at_level(logging.WARNING):
        result = parse_csv(path)

    assert result["data"] == [{"A": "7", "B": "8", "C": "9"}]
    assert "expected 3 fields, found 2" in caplog.text


def test_one_bad_row_does_not_discard_the_file(tmp_path: Path) -> None:
    """A single malformed line must not cost every other row in the file.

    This is why the field count is checked per row instead of relying on
    zip(strict=True), whose ValueError would escape to the module's broad handler.
    """
    lines = "\n".join(f"{i}\t{i}\t{i}" for i in range(20))
    path = _write(tmp_path / "mostly-good.tsv", f"A\tB\tC\n{lines}\nBROKEN\n")

    result = parse_tsv(path)

    assert len(result["data"]) == 20
    assert not any("Error parsing" in comment for comment in result["comments"])


# --- error handling ----------------------------------------------------------


def test_parse_tsv_reports_unreadable_file(tmp_path: Path) -> None:
    """An unreadable file is reported as a comment, not raised."""
    result = parse_tsv(str(tmp_path / "does-not-exist.tsv"))

    assert result["data"] == []
    assert any("Error parsing TSV file" in comment for comment in result["comments"])


def test_parse_csv_reports_unreadable_file(tmp_path: Path) -> None:
    """The CSV parser degrades the same way."""
    result = parse_csv(str(tmp_path / "does-not-exist.csv"))

    assert result["data"] == []
    assert any("Error parsing CSV file" in comment for comment in result["comments"])


def test_md5_failure_returns_none(caplog: pytest.LogCaptureFixture) -> None:
    """An unreadable artifact has no fabricated checksum or ERROR log."""
    caplog.set_level(logging.ERROR, logger=SUMMARY_LOGGER)

    with patch("vntyper.scripts.summary.open", side_effect=OSError("hash open failed")):
        assert md5sum("artifact.tsv") is None

    errors = [record for record in caplog.records if record.name == SUMMARY_LOGGER and record.levelno == logging.ERROR]
    assert errors == []


def test_csv_failure_returns_error_comment() -> None:
    """An unreadable CSV remains visible through its exact structured comment."""
    with patch("vntyper.scripts.summary.open", side_effect=OSError("csv open failed")):
        result = parse_csv("broken.csv")

    assert result == {"comments": ["Error parsing CSV file: csv open failed"], "data": []}


def test_json_failure_returns_error_mapping() -> None:
    """An invalid JSON artifact remains visible through its exact error mapping."""
    error = json.JSONDecodeError("json decode failed", "not-json", 0)

    with patch("vntyper.scripts.summary.open", side_effect=error):
        result = parse_json_file("broken.json")

    assert result == {"error": "Error reading JSON file: json decode failed: line 1 column 1 (char 0)"}


def test_tsv_failure_returns_error_comment() -> None:
    """An unreadable TSV remains visible through its exact structured comment."""
    with patch("vntyper.scripts.summary.open", side_effect=OSError("tsv open failed")):
        result = parse_tsv("broken.tsv")

    assert result == {"comments": ["Error parsing TSV file: tsv open failed"], "data": []}


def test_header_only_file_yields_no_rows(tmp_path: Path) -> None:
    """A file with only a header is valid and contributes no data."""
    assert parse_tsv(_write(tmp_path / "hdr.tsv", "A\tB\n"))["data"] == []
