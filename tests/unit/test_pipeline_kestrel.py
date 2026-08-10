"""The extracted pipeline-to-Kestrel stage preserves its atomic boundary."""

import logging
from pathlib import Path
from typing import Any, get_type_hints
from unittest import mock

import pytest

from vntyper.scripts import pipeline_kestrel
from vntyper.scripts.pipeline_kestrel import run_kestrel_stage

pytestmark = pytest.mark.unit


def _kwargs(tmp_path: Path) -> dict[str, Any]:
    return {
        "fastq_files": ("R1.fastq.gz", "R2.fastq.gz", "single.fastq.gz"),
        "dirs": {"kestrel": str(tmp_path / "kestrel")},
        "config": {
            "tools": {"kestrel": "kestrel.jar"},
            "reference_data": {"muc1_reference_vntr": "/refs/muc1.fa"},
        },
        "sample_name": "patient-7",
        "log_level": logging.DEBUG,
        "cwd": "/project/root",
        "summary": {"steps": []},
        "summary_file_path": str(tmp_path / "pipeline_summary.json"),
        "runner": mock.Mock(),
    }


def test_stage_forwards_the_complete_tuple_and_records_the_exact_summary(tmp_path: Path, monkeypatch) -> None:
    kwargs = _kwargs(tmp_path)
    record = mock.Mock()
    monkeypatch.setattr(pipeline_kestrel, "record_step", record)

    run_kestrel_stage(**kwargs)

    runner = kwargs["runner"]
    assert isinstance(runner, mock.Mock)
    runner.assert_called_once_with(
        vcf_path=tmp_path / "kestrel" / "output.vcf",
        output_dir=tmp_path / "kestrel",
        fastq_files=("R1.fastq.gz", "R2.fastq.gz", "single.fastq.gz"),
        reference_vntr="/refs/muc1.fa",
        kestrel_path="kestrel.jar",
        config=kwargs["config"],
        sample_name="patient-7",
        log_level=logging.DEBUG,
        cwd="/project/root",
    )
    positional = record.call_args.args
    assert positional[1:5] == (
        "Kestrel Genotyping",
        str(tmp_path / "kestrel" / "kestrel_result.tsv"),
        "tsv",
        "run_kestrel(...)",
    )
    assert record.call_args.kwargs == {"write_summary_path": str(tmp_path / "pipeline_summary.json")}


def test_stage_summary_annotation_matches_the_real_pipeline_payload() -> None:
    assert get_type_hints(run_kestrel_stage)["summary"] == dict[str, Any]


@pytest.mark.parametrize("fastq_files", [(), ("same.fastq.gz", "same.fastq.gz")])
def test_stage_rejects_empty_or_duplicate_inputs_before_runner(
    fastq_files: tuple[str, ...], tmp_path: Path, monkeypatch
) -> None:
    kwargs = _kwargs(tmp_path)
    kwargs["fastq_files"] = fastq_files
    record = mock.Mock()
    monkeypatch.setattr(pipeline_kestrel, "record_step", record)

    with pytest.raises(ValueError, match="FASTQ"):
        run_kestrel_stage(**kwargs)

    runner = kwargs["runner"]
    assert isinstance(runner, mock.Mock)
    runner.assert_not_called()
    record.assert_not_called()


def test_runner_failure_cannot_create_a_successful_summary_step(tmp_path: Path, monkeypatch) -> None:
    kwargs = _kwargs(tmp_path)
    kwargs["runner"] = mock.Mock(side_effect=RuntimeError("Kestrel failed"))
    record = mock.Mock()
    monkeypatch.setattr(pipeline_kestrel, "record_step", record)

    with pytest.raises(RuntimeError, match="Kestrel failed"):
        run_kestrel_stage(**kwargs)

    record.assert_not_called()


def test_summary_failure_rolls_back_a_partially_appended_success_step(tmp_path: Path, monkeypatch) -> None:
    kwargs = _kwargs(tmp_path)
    summary = kwargs["summary"]
    assert isinstance(summary, dict)

    def append_then_fail(*args: object, **kwargs: object) -> None:
        del args, kwargs
        summary["steps"].append({"step": "Kestrel Genotyping"})
        raise OSError("summary write failed")

    monkeypatch.setattr(pipeline_kestrel, "record_step", append_then_fail)

    with pytest.raises(OSError, match="summary write failed"):
        run_kestrel_stage(**kwargs)

    assert summary["steps"] == []
