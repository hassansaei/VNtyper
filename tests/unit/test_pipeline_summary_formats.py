"""``--summary-formats`` writes a per-step table and a rows file per format (#119).

The flattening itself is pinned in ``test_summary_flattening.py`` and the writers in
``test_summary_writers.py``; these tests prove ``run_pipeline`` calls them with the
summary it wrote, under the names the docs promise, and only for the formats requested.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts import summary_steps
from vntyper.scripts.summary_flattening import LONG_COLUMNS, long_rows, run_columns

pytestmark = pytest.mark.unit

KESTREL_RESULT = (
    "## VNtyper Kestrel result\n"
    "## VNtyper Version: 2.0.27\n"
    "Motifs\tPOS\tREF\tALT\tConfidence\n"
    "X-X\t67\tG\tGG\tHigh_Precision\n"
)


def _write_kestrel_result(*args: Any, **kwargs: Any) -> None:
    """Stand in for ``run_kestrel``: leave a real result where ``record_step`` reads it."""
    del args
    result = Path(kwargs["output_dir"]) / "kestrel_result.tsv"
    result.parent.mkdir(parents=True, exist_ok=True)
    result.write_text(KESTREL_RESULT, encoding="utf-8")


def _read(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = [{str(field): str(value) for field, value in row.items()} for row in reader]
        return [str(name) for name in reader.fieldnames or []], rows


def _written_summary(output_dir: Path) -> dict[str, Any]:
    return json.loads((output_dir / "pipeline_summary.json").read_text(encoding="utf-8"))


def test_a_requested_format_writes_its_table_and_rows_file_and_no_other(tmp_path: Path) -> None:
    output_dir = tmp_path / "results"

    harness = run_pipeline_under_harness(output_dir, summary_formats=["csv"])

    assert harness.error is None
    assert (output_dir / "pipeline_summary.csv").exists()
    assert (output_dir / "pipeline_summary_rows.csv").exists()
    assert not (output_dir / "pipeline_summary.tsv").exists()
    assert not (output_dir / "pipeline_summary_rows.tsv").exists()


def test_no_requested_format_writes_no_table(tmp_path: Path) -> None:
    output_dir = tmp_path / "results"

    harness = run_pipeline_under_harness(output_dir)

    assert harness.error is None
    assert (output_dir / "pipeline_summary.json").exists()
    assert sorted(path.name for path in output_dir.glob("pipeline_summary*")) == ["pipeline_summary.json"]


def test_the_table_starts_with_the_run_provenance_of_the_written_json(tmp_path: Path) -> None:
    output_dir = tmp_path / "results"

    harness = run_pipeline_under_harness(output_dir, summary_formats=["csv"])

    assert harness.error is None
    written = _written_summary(output_dir)
    run = run_columns(written)
    header, rows = _read(output_dir / "pipeline_summary.csv", ",")
    assert header[: len(run)] == list(run)
    assert len(rows) == len(written["steps"])
    assert {row["run_decision_profile_sha256"] for row in rows} == {written["decision_profile_sha256"]}
    assert [row["step"] for row in rows] == [step["step"] for step in written["steps"]]


def test_the_rows_file_carries_the_kestrel_result_row(tmp_path: Path) -> None:
    output_dir = tmp_path / "results"

    harness = run_pipeline_under_harness(
        output_dir,
        summary_formats=["tsv"],
        stage_side_effects={"run_kestrel": _write_kestrel_result},
    )

    assert harness.error is None
    written = _written_summary(output_dir)
    header, rows = _read(output_dir / "pipeline_summary_rows.tsv", "\t")
    assert header == list(LONG_COLUMNS)
    assert rows == long_rows(written)
    assert {
        "step": summary_steps.STEP_KESTREL,
        "row_index": "0",
        "field": "Confidence",
        "value": "High_Precision",
    } in rows
    table_header, table_rows = _read(output_dir / "pipeline_summary.tsv", "\t")
    kestrel = next(row for row in table_rows if row["step"] == summary_steps.STEP_KESTREL)
    assert kestrel["parsed_result_data_Confidence"] == "High_Precision"
    assert kestrel["parsed_result_comments"] == "VNtyper Kestrel result | VNtyper Version: 2.0.27"
    assert '{"' not in (output_dir / "pipeline_summary.tsv").read_text(encoding="utf-8")
