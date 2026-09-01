"""The extracted pipeline-to-Kestrel stage preserves its atomic boundary."""

import logging
from pathlib import Path
from typing import Any, get_type_hints
from unittest import mock

import pytest

from vntyper.scripts import pipeline_kestrel
from vntyper.scripts.pipeline_kestrel import run_kestrel_stage
from vntyper.scripts.run_configuration import resolve_run_configuration

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
        "threads": 4,
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
        threads=4,
    )
    positional = record.call_args.args
    assert positional[1:5] == (
        "Kestrel Genotyping",
        str(tmp_path / "kestrel" / "kestrel_result.tsv"),
        "tsv",
        "run_kestrel(...)",
    )
    assert record.call_args.kwargs == {"write_summary_path": str(tmp_path / "pipeline_summary.json")}


def test_stage_forwards_one_explicit_decision_and_runtime_context(tmp_path: Path, monkeypatch) -> None:
    kwargs = _kwargs(tmp_path)
    run = resolve_run_configuration()
    kwargs.update(
        resolved_component=run.kestrel,
        runtime_component=run.kestrel_runtime,
        custom_context_active=True,
    )
    monkeypatch.setattr(pipeline_kestrel, "record_step", mock.Mock())

    run_kestrel_stage(**kwargs)

    runner = kwargs["runner"]
    assert isinstance(runner, mock.Mock)
    assert runner.call_args.kwargs["resolved_component"] is run.kestrel
    assert runner.call_args.kwargs["runtime_component"] is run.kestrel_runtime
    assert runner.call_args.kwargs["custom_context_active"] is True
    assert kwargs["summary"]["kestrel_counting_mode"] == "split"


def test_custom_stage_cannot_fall_back_to_packaged_decisions(tmp_path: Path, monkeypatch) -> None:
    kwargs = _kwargs(tmp_path)
    kwargs["custom_context_active"] = True
    monkeypatch.setattr(pipeline_kestrel, "record_step", mock.Mock())

    with pytest.raises(ValueError, match="custom Kestrel run context requires an explicit resolved component"):
        run_kestrel_stage(**kwargs)

    runner = kwargs["runner"]
    assert isinstance(runner, mock.Mock)
    runner.assert_not_called()


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


# ---------------------------------------------------------------------------
# The thread budget reaches Kestrel (#262)
# ---------------------------------------------------------------------------


def test_the_stage_forwards_the_runs_thread_count(tmp_path: Path, monkeypatch) -> None:
    """The KAnalyze counting step is what needs it; nothing else in the stage does."""
    kwargs = _kwargs(tmp_path)
    kwargs["threads"] = 11
    monkeypatch.setattr(pipeline_kestrel, "record_step", mock.Mock())

    run_kestrel_stage(**kwargs)

    runner = kwargs["runner"]
    assert isinstance(runner, mock.Mock)
    assert runner.call_args.kwargs["threads"] == 11


def test_the_stage_has_a_thread_default_so_existing_callers_still_bind(tmp_path: Path, monkeypatch) -> None:
    """A required parameter would break every caller that has not been updated yet.

    The default matches ``config.json``'s ``default_values.threads``, so a caller that
    omits it gets the same allocation the CLI would have chosen.
    """
    kwargs = _kwargs(tmp_path)
    del kwargs["threads"]
    monkeypatch.setattr(pipeline_kestrel, "record_step", mock.Mock())

    run_kestrel_stage(**kwargs)

    runner = kwargs["runner"]
    assert isinstance(runner, mock.Mock)
    assert runner.call_args.kwargs["threads"] == 4


def test_threads_is_appended_after_every_existing_run_kestrel_parameter() -> None:
    """Inserting it earlier would silently rebind arguments for a positional caller.

    ``run_kestrel``'s parameters are not keyword-only. Every caller in this repository
    uses keywords, but an external one passing ``log_level`` and ``cwd`` positionally
    would start passing ``cwd`` as ``threads`` if the new parameter went in front of
    them -- which type-checks, runs, and produces a wrong thread allocation.
    """
    import inspect

    from vntyper.scripts import kestrel_genotyping

    names = list(inspect.signature(kestrel_genotyping.run_kestrel).parameters)

    assert names[-6:-3] == ["log_level", "cwd", "threads"]
    assert names[-3:] == ["resolved_component", "runtime_component", "custom_context_active"]
