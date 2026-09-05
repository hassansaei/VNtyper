"""The operator-facing summary tables ``summary.py`` writes (#119).

``pipeline_summary.<csv|tsv>`` is one row per step with the run's provenance first;
``pipeline_summary_rows.<csv|tsv>`` is the long form of every result row. Both are
computed by ``summary_flattening`` and written here; these tests pin what lands on disk.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pytest

from vntyper.scripts import summary
from vntyper.scripts.summary import (
    convert_summary_rows_to_delimited,
    convert_summary_to_csv,
    convert_summary_to_tsv,
)
from vntyper.scripts.summary_flattening import LONG_COLUMNS, STEP_COLUMNS, long_rows, run_columns

pytestmark = pytest.mark.unit

PROFILE_SHA256 = "5d41402abc4b2a76b9719d911017c592" * 2

KESTREL_ROW = {
    "Motifs": "X-X",
    "POS": "67",
    "REF": "G",
    "ALT": "GG",
    "Depth_Score": "0.01",
    "Confidence": "High_Precision",
    "Nomenclature": "59dupC",
}

ADVNTR_ROWS = [
    {"VID": "25561", "Variant": "I22_2_G_LEN1", "NumberOfSupportingReads": "7", "Pvalue": "0.001"},
    {"VID": "25561", "Variant": "D58_2&D59_2", "NumberOfSupportingReads": "3", "Pvalue": "0.04"},
]


def _summary(steps: list[dict[str, Any]] | None = None, kestrel_row: dict[str, str] | None = None) -> dict[str, Any]:
    kestrel = {
        "step": "Kestrel Genotyping",
        "start": "2026-09-05T10:15:01.000000",
        "end": "2026-09-05T10:15:02.000000",
        "command": "run_kestrel(...)",
        "result_file": "/work/results/kestrel/kestrel_result.tsv",
        "file_type": "tsv",
        "md5sum": "d41d8cd98f00b204e9800998ecf8427e",
        "parsed_result": {
            "comments": ["VNtyper Kestrel result"],
            "data": [KESTREL_ROW if kestrel_row is None else kestrel_row],
        },
    }
    advntr = {
        "step": "adVNTR Genotyping",
        "start": "2026-09-05T10:16:01.000000",
        "end": "2026-09-05T10:17:02.000000",
        "command": "run_advntr(...)",
        "result_file": "/work/results/advntr/output_adVNTR_result.tsv",
        "file_type": "tsv",
        "md5sum": "9e107d9d372bb6826bd81d3542a419d6",
        "parsed_result": {"comments": [], "data": ADVNTR_ROWS},
    }
    return {
        "schema_version": 3,
        "decision_policy": "legacy-selection-v1",
        "advntr_evidence_digest": None,
        "decision_profile_id": "vntyper-default",
        "decision_profile_revision": "1",
        "decision_profile_kind": "packaged",
        "decision_profile_source": "package",
        "decision_profile_sha256": PROFILE_SHA256,
        "decision_profile_snapshot": "provenance/decision_profile.json",
        "pipeline_start": "2026-09-05T10:15:00.000001",
        "version": "2.0.27",
        "input_files": {"bam": "/data/sample.bam"},
        "sample_name": "sample",
        "sample_name_is_explicit": True,
        "reference_assembly_requested": "hg19",
        "reference_key_used": None,
        "reference_path": None,
        "reference_source_effective": "not-required",
        "region_resolved": "chr1:155158000-155163000",
        "steps": [kestrel, advntr] if steps is None else steps,
        "kestrel_counting_mode": "split",
        "pipeline_end": "2026-09-05T10:18:30.000002",
    }


def _read(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = [{str(field): str(value) for field, value in row.items()} for row in reader]
        return [str(name) for name in reader.fieldnames or []], rows


EXPECTED_PARSED_COLUMNS = sorted(
    ["parsed_result_comments", "parsed_result_n_rows", *[f"parsed_result_data_{field}" for field in KESTREL_ROW]]
)


# ---------------------------------------------------------------------------
# The one-row-per-step table
# ---------------------------------------------------------------------------


def test_the_tsv_puts_run_columns_first_then_step_columns_then_sorted_result_columns(tmp_path: Path) -> None:
    table = _summary()
    convert_summary_to_tsv(table, tmp_path / "pipeline_summary.tsv")

    header, _rows = _read(tmp_path / "pipeline_summary.tsv", "\t")

    assert header == [*run_columns(table), *STEP_COLUMNS, *EXPECTED_PARSED_COLUMNS]
    assert header[0] == "run_schema_version"


def test_every_row_carries_the_run_provenance(tmp_path: Path) -> None:
    convert_summary_to_csv(_summary(), tmp_path / "pipeline_summary.csv")

    _header, rows = _read(tmp_path / "pipeline_summary.csv", ",")

    assert len(rows) == 2
    assert [row["run_decision_profile_sha256"] for row in rows] == [PROFILE_SHA256, PROFILE_SHA256]
    assert [row["run_version"] for row in rows] == ["2.0.27", "2.0.27"]
    assert [row["run_sample_name_is_explicit"] for row in rows] == ["True", "True"]
    assert [row["step"] for row in rows] == ["Kestrel Genotyping", "adVNTR Genotyping"]


def test_cells_that_do_not_apply_to_a_step_are_blank(tmp_path: Path) -> None:
    convert_summary_to_csv(_summary(), tmp_path / "pipeline_summary.csv")

    _header, (kestrel, advntr) = _read(tmp_path / "pipeline_summary.csv", ",")

    assert kestrel["parsed_result_data_Motifs"] == "X-X"
    assert kestrel["parsed_result_n_rows"] == ""
    assert advntr["parsed_result_n_rows"] == "2"
    assert advntr["parsed_result_data_Motifs"] == ""
    assert advntr["run_reference_key_used"] == ""


def test_no_cell_contains_embedded_json(tmp_path: Path) -> None:
    convert_summary_to_tsv(_summary(), tmp_path / "pipeline_summary.tsv")
    convert_summary_to_csv(_summary(), tmp_path / "pipeline_summary.csv")

    for name in ("pipeline_summary.tsv", "pipeline_summary.csv"):
        text = (tmp_path / name).read_text(encoding="utf-8")
        assert '{"' not in text, name
        assert "; {" not in text, name


def test_csv_and_tsv_carry_the_same_table(tmp_path: Path) -> None:
    convert_summary_to_csv(_summary(), tmp_path / "pipeline_summary.csv")
    convert_summary_to_tsv(_summary(), tmp_path / "pipeline_summary.tsv")

    assert _read(tmp_path / "pipeline_summary.csv", ",") == _read(tmp_path / "pipeline_summary.tsv", "\t")


def test_a_summary_with_no_steps_writes_the_header_only(tmp_path: Path) -> None:
    empty = _summary(steps=[])
    convert_summary_to_tsv(empty, tmp_path / "pipeline_summary.tsv")

    lines = (tmp_path / "pipeline_summary.tsv").read_text(encoding="utf-8").splitlines()

    assert len(lines) == 1
    assert lines[0].split("\t") == [*run_columns(empty), *STEP_COLUMNS]


def test_a_cell_containing_the_delimiter_round_trips(tmp_path: Path) -> None:
    row = {**KESTREL_ROW, "Nomenclature_Note": "comma, and\ttab"}
    convert_summary_to_csv(_summary(kestrel_row=row), tmp_path / "pipeline_summary.csv")
    convert_summary_to_tsv(_summary(kestrel_row=row), tmp_path / "pipeline_summary.tsv")

    _h, (csv_kestrel, _a) = _read(tmp_path / "pipeline_summary.csv", ",")
    _h, (tsv_kestrel, _a) = _read(tmp_path / "pipeline_summary.tsv", "\t")

    assert csv_kestrel["parsed_result_data_Nomenclature_Note"] == "comma, and\ttab"
    assert tsv_kestrel["parsed_result_data_Nomenclature_Note"] == "comma, and\ttab"


# ---------------------------------------------------------------------------
# The long rows file
# ---------------------------------------------------------------------------


def test_the_rows_file_lists_every_result_row(tmp_path: Path) -> None:
    table = _summary()
    convert_summary_rows_to_delimited(table, tmp_path / "pipeline_summary_rows.tsv", delimiter="\t")

    header, rows = _read(tmp_path / "pipeline_summary_rows.tsv", "\t")

    assert header == list(LONG_COLUMNS)
    assert rows == long_rows(table)
    assert len(rows) == len(KESTREL_ROW) + len(ADVNTR_ROWS) * len(ADVNTR_ROWS[0])
    assert {"step": "adVNTR Genotyping", "row_index": "1", "field": "Variant", "value": "D58_2&D59_2"} in rows


def test_the_rows_file_defaults_to_a_comma_delimiter(tmp_path: Path) -> None:
    table = _summary()
    convert_summary_rows_to_delimited(table, tmp_path / "pipeline_summary_rows.csv")

    header, rows = _read(tmp_path / "pipeline_summary_rows.csv", ",")

    assert header == ["step", "row_index", "field", "value"]
    assert rows == long_rows(table)


def test_the_rows_file_of_a_summary_with_no_rows_is_the_header_only(tmp_path: Path) -> None:
    convert_summary_rows_to_delimited(_summary(steps=[]), tmp_path / "pipeline_summary_rows.csv")

    assert (tmp_path / "pipeline_summary_rows.csv").read_text(encoding="utf-8").splitlines() == [
        "step,row_index,field,value"
    ]


# ---------------------------------------------------------------------------
# What is gone
# ---------------------------------------------------------------------------


def test_the_json_embedding_flattener_and_the_demo_are_gone() -> None:
    source = Path(summary.__file__).read_text(encoding="utf-8")

    assert not hasattr(summary, "flatten_dict")
    assert '__name__ == "__main__"' not in source
    assert "json.dumps" not in source
