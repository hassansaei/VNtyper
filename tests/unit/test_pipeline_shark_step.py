"""Pipeline wiring for the readable SHARK step summary."""

from __future__ import annotations

import json
from datetime import datetime as real_datetime
from pathlib import Path
from unittest import mock

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.modules.shark import shark_filtering
from vntyper.scripts import pipeline as pipeline_module
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def fake_shark_filter(kept_r1: int, kept_r2: int):
    """Return a filter stub that writes the same paths as real SHARK."""

    def _filter(*, output_dir, sample_name, **_kwargs):
        out1 = Path(output_dir) / f"{sample_name}_shark_R1.fastq"
        out2 = Path(output_dir) / f"{sample_name}_shark_R2.fastq"
        out1.write_text("".join(f"@r{i}\nACGT\n+\nFFFF\n" for i in range(kept_r1)), encoding="utf-8")
        out2.write_text("".join(f"@r{i}\nACGT\n+\nFFFF\n" for i in range(kept_r2)), encoding="utf-8")
        return str(out1), str(out2)

    return _filter


def shark_step(output_dir: Path) -> dict:
    """Return the sole recorded SHARK step."""
    summary = json.loads((output_dir / "pipeline_summary.json").read_text(encoding="utf-8"))
    steps = [step for step in summary["steps"] if step["step"] == "SHARK Filtering"]
    assert len(steps) == 1
    return steps[0]


def run_with_shark(tmp_path: Path, kept_r1: int = 2, kept_r2: int = 2) -> tuple[Path, dict]:
    """Run the real orchestration with only the SHARK binary call replaced."""
    out = tmp_path / "out"
    with mock.patch(
        "vntyper.modules.shark.shark_filtering.run_shark_filter",
        side_effect=fake_shark_filter(kept_r1, kept_r2),
    ):
        run_pipeline_under_harness(
            out,
            bam=None,
            fastq1="/in/r1.fastq.gz",
            fastq2="/in/r2.fastq.gz",
            extra_modules=["shark"],
        )
    return out, shark_step(out)


def test_the_step_records_a_readable_file_the_run_wrote(tmp_path) -> None:
    """The old path and unsupported ``fastq`` type made every successful step unreadable."""
    _out, step = run_with_shark(tmp_path)

    assert "result_file_missing" not in step
    assert Path(step["result_file"]).exists()
    assert Path(step["result_file"]).name == "sample_shark_step.json"
    assert step["file_type"] == "json"
    assert step["md5sum"] is not None


def test_the_step_records_exact_fastq_paths_and_string_counts(tmp_path) -> None:
    """The metadata points at the actual uncompressed SHARK outputs."""
    out, step = run_with_shark(tmp_path, kept_r1=3, kept_r2=3)

    fastq_dir = out / "fastq_bam_processing"
    assert step["parsed_result"] == {
        "filtered_fastq_1": str(fastq_dir / "sample_shark_R1.fastq"),
        "filtered_fastq_2": str(fastq_dir / "sample_shark_R2.fastq"),
        "kept_reads_r1": "3",
        "kept_reads_r2": "3",
    }


def test_a_filter_that_kept_nothing_records_readable_string_zero(tmp_path) -> None:
    """Zero kept reads remains distinct from an unreadable result."""
    _out, step = run_with_shark(tmp_path, kept_r1=0, kept_r2=0)

    assert "result_file_missing" not in step
    assert step["parsed_result"]["kept_reads_r1"] == "0"
    assert step["parsed_result"]["kept_reads_r2"] == "0"


def test_a_count_failure_stops_before_the_step_is_recorded(tmp_path) -> None:
    """Counting and sidecar creation must finish before the step end is captured and recorded."""
    out = tmp_path / "out"

    def truncated_filter(*, output_dir, sample_name, **_kwargs):
        r1 = Path(output_dir) / f"{sample_name}_shark_R1.fastq"
        r2 = Path(output_dir) / f"{sample_name}_shark_R2.fastq"
        r1.write_text("@r0\nACGT\n+\n", encoding="utf-8")
        r2.write_text("@r0\nACGT\n+\nFFFF\n", encoding="utf-8")
        return str(r1), str(r2)

    with mock.patch("vntyper.modules.shark.shark_filtering.run_shark_filter", side_effect=truncated_filter):
        harness = run_pipeline_under_harness(
            out,
            bam=None,
            fastq1="/in/r1.fastq.gz",
            fastq2="/in/r2.fastq.gz",
            extra_modules=["shark"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not (out / "fastq_bam_processing" / "sample_shark_step.json").exists()
    assert not (out / "pipeline_summary.json").exists()


def test_the_step_end_is_captured_after_counting_and_writing(tmp_path, monkeypatch) -> None:
    """The recorded duration includes the work required to make the step readable."""
    events: list[str] = []
    real_write_summary = shark_filtering.write_shark_step_summary
    write_filter_outputs = fake_shark_filter(1, 1)

    class RecordingDateTime:
        @classmethod
        def now(cls, tz=None):
            events.append("now")
            return real_datetime.now(tz)

    def recording_filter(**kwargs):
        events.append("filter")
        return write_filter_outputs(**kwargs)

    def recording_summary(*args, **kwargs):
        events.append("summary_start")
        payload = real_write_summary(*args, **kwargs)
        events.append("summary_done")
        return payload

    monkeypatch.setattr(pipeline_module, "datetime", RecordingDateTime)
    monkeypatch.setattr(shark_filtering, "run_shark_filter", recording_filter)
    monkeypatch.setattr(shark_filtering, "write_shark_step_summary", recording_summary)

    run_pipeline_under_harness(
        tmp_path / "out",
        bam=None,
        fastq1="/in/r1.fastq.gz",
        fastq2="/in/r2.fastq.gz",
        extra_modules=["shark"],
    )

    filter_index = events.index("filter")
    assert events[filter_index - 1 : filter_index + 4] == [
        "now",
        "filter",
        "summary_start",
        "summary_done",
        "now",
    ]


def test_pipeline_forwards_the_resolved_empty_policy_and_runtime_references(tmp_path) -> None:
    seen: dict[str, object] = {}
    write_filter_outputs = fake_shark_filter(1, 1)

    def recording_filter(**kwargs):
        seen.update(kwargs)
        return write_filter_outputs(**kwargs)

    with mock.patch("vntyper.modules.shark.shark_filtering.run_shark_filter", side_effect=recording_filter):
        run_pipeline_under_harness(
            tmp_path / "out",
            bam=None,
            fastq1="/in/r1.fastq.gz",
            fastq2="/in/r2.fastq.gz",
            extra_modules=["shark"],
        )

    run = resolve_run_configuration()
    assert seen["resolved_component"] == run.shark == {}
    assert seen["config"] == run.shark_runtime
    assert seen["custom_context_active"] is False
