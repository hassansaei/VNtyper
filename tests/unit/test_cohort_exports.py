"""The machine-readable cohort outputs, and the table that de-pseudonymises them.

`vntyper cohort --output-format csv,tsv,json` writes the aggregated Kestrel and adVNTR
frames beside the HTML report. These are the files a downstream analysis actually reads,
so what is in them and what they are called is a contract; none of it had a test.

Characterisation throughout, including one `..._today` recording a defect the split was
not authorised to fix.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts.cohort_exports import (
    parse_output_formats,
    write_cohort_frame,
    write_pseudonymization_table,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Parsing the --output-format argument
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "argument,formats",
    [
        ("", []),
        ("csv", ["csv"]),
        ("csv,tsv,json", ["csv", "tsv", "json"]),
        ("CSV, TSV", ["csv", "tsv"]),
        ("  csv  ", ["csv"]),
        ("csv,,tsv", ["csv", "tsv"]),
        (",", []),
        ("csv,csv", ["csv", "csv"]),
    ],
)
def test_the_format_list_is_split_lower_cased_and_stripped(argument: str, formats: list[str]) -> None:
    assert parse_output_formats(argument) == formats


def test_an_unrecognised_format_is_carried_through_and_simply_never_matches() -> None:
    """Nothing validates the list; an unsupported name is silently inert rather than
    an error, so `--output-format xlsx` produces the HTML report and nothing else."""
    assert parse_output_formats("xlsx") == ["xlsx"]


# ---------------------------------------------------------------------------
# Writing one frame
# ---------------------------------------------------------------------------


def _frame() -> pd.DataFrame:
    return pd.DataFrame([{"Sample": "s1", "Motif": "5"}, {"Sample": "s2", "Motif": "8"}])


def test_no_requested_format_writes_no_file(tmp_path) -> None:
    write_cohort_frame(_frame(), tmp_path, "cohort_kestrel", "Kestrel", [])

    assert list(tmp_path.iterdir()) == []


def test_an_empty_frame_writes_no_file(tmp_path) -> None:
    """A cohort in which no sample produced adVNTR results gets no `cohort_advntr.csv`
    rather than a header-only one."""
    write_cohort_frame(pd.DataFrame(), tmp_path, "cohort_advntr", "adVNTR", ["csv", "tsv", "json"])

    assert list(tmp_path.iterdir()) == []


def test_the_csv_is_written_without_the_index(tmp_path) -> None:
    write_cohort_frame(_frame(), tmp_path, "cohort_kestrel", "Kestrel", ["csv"])

    assert (tmp_path / "cohort_kestrel.csv").read_text(encoding="utf-8").splitlines() == [
        "Sample,Motif",
        "s1,5",
        "s2,8",
    ]


def test_the_tsv_is_tab_separated(tmp_path) -> None:
    write_cohort_frame(_frame(), tmp_path, "cohort_kestrel", "Kestrel", ["tsv"])

    assert (tmp_path / "cohort_kestrel.tsv").read_text(encoding="utf-8").splitlines()[0] == "Sample\tMotif"


def test_the_json_is_a_list_of_row_objects(tmp_path) -> None:
    write_cohort_frame(_frame(), tmp_path, "cohort_kestrel", "Kestrel", ["json"])

    assert json.loads((tmp_path / "cohort_kestrel.json").read_text(encoding="utf-8")) == [
        {"Sample": "s1", "Motif": "5"},
        {"Sample": "s2", "Motif": "8"},
    ]


def test_every_requested_format_is_written(tmp_path) -> None:
    write_cohort_frame(_frame(), tmp_path, "cohort_advntr", "adVNTR", ["csv", "tsv", "json"])

    assert {p.name for p in tmp_path.iterdir()} == {
        "cohort_advntr.csv",
        "cohort_advntr.tsv",
        "cohort_advntr.json",
    }


def test_each_written_file_is_announced_under_its_algorithm_name(tmp_path, caplog) -> None:
    """The log line is how a user finds the file; the label distinguishes the two
    frames, which otherwise write identical-looking messages."""
    caplog.set_level(logging.INFO, logger="vntyper.scripts.cohort_exports")

    write_cohort_frame(_frame(), tmp_path, "cohort_advntr", "adVNTR", ["csv"])

    assert f"Cohort adVNTR CSV written to: {tmp_path / 'cohort_advntr.csv'}" in caplog.text


def test_a_frame_carrying_the_reductions_working_columns_exports_them_today(tmp_path) -> None:
    """Characterisation of a live defect, at the boundary where it is visible.

    Nothing here strips columns: whatever the frame carries is exported. `aggregate_cohort`
    renders the report first, and rendering annotates the frame with `__row_result` and
    `__unified`, so those two internal columns are in every export. See
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-internal-columns-leak-into-exports.md`.
    """
    frame = pd.DataFrame([{"Sample": "s1", "__row_result": "High_Precision", "__unified": "Positive"}])

    write_cohort_frame(frame, tmp_path, "cohort_kestrel", "Kestrel", ["csv"])

    assert (tmp_path / "cohort_kestrel.csv").read_text(encoding="utf-8").splitlines()[0] == (
        "Sample,__row_result,__unified"
    )


def test_the_output_directory_may_be_given_as_a_string(tmp_path) -> None:
    write_cohort_frame(_frame(), str(tmp_path), "cohort_kestrel", "Kestrel", ["csv"])

    assert (tmp_path / "cohort_kestrel.csv").is_file()


# ---------------------------------------------------------------------------
# The pseudonymization table
# ---------------------------------------------------------------------------


def test_the_pseudonymization_table_is_a_two_column_tsv(tmp_path) -> None:
    write_pseudonymization_table(tmp_path, {"anon_65622": "sample_one", "anon_a1b2c": "sample_two"})

    assert (tmp_path / "pseudonymization_table.tsv").read_text(encoding="utf-8").splitlines() == [
        "Pseudonym\tOriginal",
        "anon_65622\tsample_one",
        "anon_a1b2c\tsample_two",
    ]


def test_an_empty_mapping_still_writes_the_header(tmp_path) -> None:
    """`aggregate_cohort` only calls this when the mapping is non-empty, so the
    header-only file is what a direct caller gets, not what a cohort run produces."""
    write_pseudonymization_table(tmp_path, {})

    assert (tmp_path / "pseudonymization_table.tsv").read_text(encoding="utf-8") == "Pseudonym\tOriginal\n"


def test_a_table_that_cannot_be_written_is_reported_rather_than_raised(tmp_path, caplog) -> None:
    """This runs after the report is already on disk, so failing here must not turn a
    completed cohort into a traceback."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_exports")

    write_pseudonymization_table(Path(tmp_path) / "no-such-directory", {"a": "b"})

    assert "Failed to write pseudonymization table" in caplog.text
