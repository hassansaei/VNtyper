"""Unit tests for ``record_step``'s handling of an absent result file (#212).

``record_step`` is the only place that sees the path a stage was supposed to have
written. When that file is missing, everything downstream degrades quietly:
``md5sum`` swallows the ``FileNotFoundError`` and returns ``None``, and ``parse_tsv``
turns it into a comment and returns ``{"comments": [...], "data": []}``. An empty
``data`` list is exactly what a stage that legitimately found no variant produces, so
"the step produced nothing" and "the step found nothing" become indistinguishable --
and both the HTML report and cohort mode render the latter as a negative.

The flag is set **only** when the file is absent. That is deliberate and is pinned
below: the golden-cohort harness compares whole step records, so an unconditionally
present key would diff every step of every run for no behavioural reason.
"""

import json
import logging
from datetime import datetime, timezone
from unittest.mock import patch

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts.summary import record_step, refresh_step

pytestmark = pytest.mark.unit

#: Fixed timestamps. `record_step` only calls `.isoformat()` on them, so the values are
#: irrelevant beyond being stable -- but the naive shape is not: `pipeline.py` passes
#: `datetime.now(timezone.utc).replace(tzinfo=None)`, and these are built the same way so
#: the recorded `start`/`end` strings match what a real run writes.
_START = datetime(2026, 1, 1, 12, 0, 0, tzinfo=timezone.utc).replace(tzinfo=None)
_END = datetime(2026, 1, 1, 12, 0, 30, tzinfo=timezone.utc).replace(tzinfo=None)


def test_record_step_flags_a_missing_result_file(tmp_path, caplog):
    """#212's second half: a step that produced nothing must not look like a negative."""
    summary = {"steps": []}
    missing = tmp_path / "absent.tsv"

    with caplog.at_level(logging.ERROR):
        record_step(summary, "Kestrel Genotyping", str(missing), "tsv", "cmd", _START, _END)

    assert summary["steps"][0]["result_file_missing"] is True
    errors = [record.getMessage() for record in caplog.records if record.levelno == logging.ERROR]
    assert any("absent.tsv" in message for message in errors), errors


def test_record_step_adds_no_key_when_the_file_exists(tmp_path):
    """Pins gate-neutrality: the golden harness compares whole step records, so an
    unconditional key would diff every case of every run."""
    path = tmp_path / "present.tsv"
    path.write_text("a\tb\n1\t2\n", encoding="utf-8")
    summary = {"steps": []}

    record_step(summary, "Kestrel Genotyping", str(path), "tsv", "cmd", _START, _END)

    assert "result_file_missing" not in summary["steps"][0]


def test_record_step_logs_nothing_at_error_when_the_file_exists(tmp_path, caplog):
    """The ERROR is the operator-facing half of the flag; a clean run must stay quiet."""
    path = tmp_path / "present.tsv"
    path.write_text("a\tb\n1\t2\n", encoding="utf-8")
    summary = {"steps": []}

    with caplog.at_level(logging.ERROR):
        record_step(summary, "Kestrel Genotyping", str(path), "tsv", "cmd", _START, _END)

    errors = [record.getMessage() for record in caplog.records if record.levelno == logging.ERROR]
    assert errors == []


def test_record_step_still_records_the_rest_of_the_step_when_the_file_is_missing(tmp_path):
    """The flag is additive: it must not replace or suppress the existing record.

    The missing file is still worth recording -- the path is the operator's only clue to
    which stage failed -- so `md5sum` staying None and `parsed_result` staying an empty
    parse is the documented shape, not a second bug to fix here.
    """
    summary = {"steps": []}
    missing = tmp_path / "absent.tsv"

    record_step(summary, "Kestrel Genotyping", str(missing), "tsv", "run_kestrel(...)", _START, _END)

    record = summary["steps"][0]
    assert record["step"] == "Kestrel Genotyping"
    assert record["result_file"] == str(missing)
    assert record["md5sum"] is None
    assert record["parsed_result"]["data"] == []


def test_parser_failure_is_recorded_and_step_is_appended(tmp_path):
    """A parser failure remains in its step and does not block later summary records."""
    path = tmp_path / "result.tsv"
    path.write_text("A\n", encoding="utf-8")
    summary = {"steps": []}

    with patch("vntyper.scripts.summary.parse_tsv", side_effect=ValueError("parser failed")):
        record_step(summary, "Broken parser", str(path), "tsv", "first command", _START, _END)

    assert len(summary["steps"]) == 1
    failed_step = summary["steps"][0]
    assert failed_step["step"] == "Broken parser"
    assert failed_step["result_file"] == str(path)
    assert failed_step["parsed_result"] == {"error": "Error parsing file: parser failed"}

    record_step(summary, "Later step", str(path), "tsv", "later command", _START, _END)

    assert [step["step"] for step in summary["steps"]] == ["Broken parser", "Later step"]


def test_conversion_summary_uses_the_first_routed_fastq_from_the_nonempty_tuple(tmp_path):
    """Tuple migration keeps the conversion artifact stable without a removed ``fastq1`` local."""
    routed = ("/routed/R1.fastq.gz", "/routed/R2.fastq.gz", "/routed/single.fastq.gz")
    output_dir = tmp_path / "out"
    harness = run_pipeline_under_harness(
        output_dir,
        stage_side_effects={"route_converted_fastqs": lambda *args, **kwargs: routed},
    )

    assert harness.error is None
    summary = json.loads((output_dir / "pipeline_summary.json").read_text(encoding="utf-8"))
    conversion = next(step for step in summary["steps"] if step["step"] == "BAM to FASTQ Conversion")
    assert conversion["result_file"] == routed[0]


def test_refresh_step_picks_up_a_result_file_rewritten_after_recording(tmp_path) -> None:
    """A step recorded before its file is rewritten must not keep the stale copy.

    The cross-caller nomenclature stage rewrites ``kestrel_result.tsv`` after the
    Kestrel step has been recorded. Without refreshing, the summary -- and the HTML
    report and cohort tables built from it, which read the summary rather than the TSV
    -- keep the pre-reconciliation row: the tier before promotion and an empty adVNTR
    column. The stored checksum no longer matched the file on disk either.
    """
    result = tmp_path / "kestrel_result.tsv"
    result.write_text("A\tB\n1\told\n")
    summary: dict = {"steps": []}
    record_step(summary, "Kestrel Genotyping", str(result), "tsv", "cmd", _START, _END)
    stale_md5 = summary["steps"][0]["md5sum"]
    assert summary["steps"][0]["parsed_result"]["data"] == [{"A": "1", "B": "old"}]

    result.write_text("A\tB\n1\tnew\n")
    assert refresh_step(summary, "Kestrel Genotyping") is True

    assert summary["steps"][0]["parsed_result"]["data"] == [{"A": "1", "B": "new"}]
    assert summary["steps"][0]["md5sum"] != stale_md5


def test_refresh_step_reports_an_unknown_step_rather_than_raising() -> None:
    """A summary without the step is a no-op, not an exception."""
    assert refresh_step({"steps": []}, "Kestrel Genotyping") is False
