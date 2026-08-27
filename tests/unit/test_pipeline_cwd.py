"""The ``cwd=`` plumbing every external tool depends on (AGENTS.md trap 7).

Every tool and reference path in ``config.json`` is relative to the process working
directory. ``run_pipeline`` therefore pins ``project_root = os.getcwd()`` on entry and
threads it as ``cwd=`` into the Java and samtools calls, so that a stage which changes
directory - or a Celery worker whose CWD is not the repo - still resolves
``vntyper/dependencies/kestrel/kestrel.jar`` and the reference FASTAs.

Nothing guarded that. Removing any one ``cwd=`` argument left the whole unit tier
green, and the failure it causes in production is a ``FileNotFoundError`` on a
reference path - or, worse, a *different* file that happens to exist relative to the
new directory.

These tests run the real ``run_pipeline`` from a directory that is not the repo root
and assert each tool call carries that directory. They also pin the two properties
that make the plumbing worth having: the value is captured **once, on entry**, so a
stage that chdirs mid-run cannot poison later stages; and an unreadable working
directory falls back to the package root rather than raising.
"""

import copy
import json
import logging
import os
from pathlib import Path
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import pipeline as pipeline_module

pytestmark = pytest.mark.unit


#: Stage name -> how ``run_pipeline`` is expected to pass the working directory.
#: Every entry is a tool that reads a path out of ``config.json``.
CWD_CARRYING_STAGES = ("validate_bam_file", "run_kestrel", "run_advntr")


def _cwd_of(harness, stage: str) -> str:
    """Return the ``cwd`` keyword ``stage`` was called with.

    Args:
        harness: The harness returned by ``run_pipeline_under_harness``.
        stage: The stage name.

    Returns:
        str: The recorded ``cwd`` argument.

    Raises:
        AssertionError: If the stage was never called or carried no ``cwd``.
    """
    kwargs = harness.kwargs(stage)
    assert "cwd" in kwargs, f"{stage}() was called without cwd=; AGENTS.md trap 7"
    return kwargs["cwd"]


@pytest.mark.parametrize("stage", CWD_CARRYING_STAGES)
def test_every_tool_call_carries_the_working_directory(stage: str, tmp_path: Path, monkeypatch) -> None:
    """Run from a directory that is not the repo root and follow it into each tool.

    Args:
        stage: The stage under test.
        tmp_path: Pytest temporary directory.
        monkeypatch: Pytest monkeypatch fixture.
    """
    run_dir = tmp_path / "elsewhere"
    run_dir.mkdir()
    monkeypatch.chdir(run_dir)
    expected = os.getcwd()
    assert expected != str(Path(__file__).resolve().parents[2]), "the test did not actually leave the repo root"

    harness = run_pipeline_under_harness(tmp_path / "out", extra_modules=["advntr"])

    assert _cwd_of(harness, stage) == expected


def test_the_working_directory_is_captured_once_at_entry(tmp_path: Path, monkeypatch) -> None:
    """A stage that changes directory must not change what later stages receive.

    ``process_bam_to_fastq`` runs before Kestrel and before adVNTR. If the pipeline
    re-read ``os.getcwd()`` at each call site instead of pinning it on entry, a stage
    that chdirs would silently repoint every later tool at the new directory.

    Args:
        tmp_path: Pytest temporary directory.
        monkeypatch: Pytest monkeypatch fixture.
    """
    run_dir = tmp_path / "start"
    run_dir.mkdir()
    wandered = tmp_path / "wandered"
    wandered.mkdir()
    monkeypatch.chdir(run_dir)
    expected = os.getcwd()

    def _chdir_mid_run(*args, **kwargs):
        os.chdir(wandered)
        return mock.DEFAULT

    harness = run_pipeline_under_harness(
        tmp_path / "out",
        extra_modules=["advntr"],
        stage_side_effects={"process_bam_to_fastq": _chdir_mid_run},
    )

    assert os.getcwd() != expected, "the side effect did not actually change directory"
    for stage in ("run_kestrel", "run_advntr"):
        assert _cwd_of(harness, stage) == expected


def test_an_unreadable_working_directory_falls_back_to_the_package_root(tmp_path: Path, monkeypatch) -> None:
    """``os.getcwd()`` raises when the CWD has been deleted; the run must continue.

    Args:
        tmp_path: Pytest temporary directory.
        monkeypatch: Pytest monkeypatch fixture.
    """
    monkeypatch.setattr(pipeline_module.os, "getcwd", mock.Mock(side_effect=FileNotFoundError("cwd is gone")))

    harness = run_pipeline_under_harness(tmp_path / "out")

    fallback = _cwd_of(harness, "run_kestrel")
    assert fallback == str(Path(pipeline_module.__file__).parent.parent.parent)


def test_the_fastq_path_also_carries_the_working_directory(tmp_path: Path, monkeypatch) -> None:
    """FASTQ input reaches Kestrel through a different branch, with the same rule.

    Args:
        tmp_path: Pytest temporary directory.
        monkeypatch: Pytest monkeypatch fixture.
    """
    run_dir = tmp_path / "elsewhere"
    run_dir.mkdir()
    monkeypatch.chdir(run_dir)
    expected = os.getcwd()

    harness = run_pipeline_under_harness(
        tmp_path / "out",
        bam=None,
        fastq1="/in/r1.fastq.gz",
        fastq2="/in/r2.fastq.gz",
    )

    assert _cwd_of(harness, "run_kestrel") == expected


def test_kestrel_stage_boundary_preserves_inputs_configuration_and_summary(tmp_path: Path, monkeypatch) -> None:
    """Extracting the Kestrel stage must preserve every caller-visible argument and artifact.

    The expected dict is **exhaustive on purpose**: it fails if a parameter is added as
    well as if one is dropped, so growing ``run_kestrel``'s call contract has to be a
    decision rather than a side effect. ``threads`` was added that way in #262, for the
    KAnalyze counting step.
    """
    run_dir = tmp_path / "elsewhere"
    run_dir.mkdir()
    monkeypatch.chdir(run_dir)
    out = tmp_path / "out"

    harness = run_pipeline_under_harness(
        out,
        bam=None,
        fastq1="/input/R1.fastq.gz",
        fastq2="/input/R2.fastq.gz",
        sample_name="patient-7",
        log_level=logging.DEBUG,
    )

    assert harness.kwargs("run_kestrel") == {
        "vcf_path": out / "kestrel" / "output.vcf",
        "output_dir": out / "kestrel",
        "fastq_files": (
            str(out / "fastq_bam_processing" / "output_R1.fastq.gz"),
            str(out / "fastq_bam_processing" / "output_R2.fastq.gz"),
            str(out / "fastq_bam_processing" / "output_single.fastq.gz"),
        ),
        "reference_vntr": "/refs/muc1.fa",
        "kestrel_path": "kestrel.jar",
        "config": MINIMAL_CONFIG,
        "sample_name": "patient-7",
        "log_level": logging.DEBUG,
        "cwd": str(run_dir),
        "threads": 4,
    }
    summary = json.loads((out / "pipeline_summary.json").read_text(encoding="utf-8"))
    kestrel_steps = [step for step in summary["steps"] if step["step"] == "Kestrel Genotyping"]
    assert len(kestrel_steps) == 1
    assert kestrel_steps[0]["result_file"] == str(out / "kestrel" / "kestrel_result.tsv")


def test_a_cram_input_validates_with_the_working_directory(tmp_path: Path, monkeypatch) -> None:
    """CRAM goes through the same ``validate_bam_file``, and needs the same cwd.

    Args:
        tmp_path: Pytest temporary directory.
        monkeypatch: Pytest monkeypatch fixture.
    """
    run_dir = tmp_path / "elsewhere"
    run_dir.mkdir()
    monkeypatch.chdir(run_dir)
    expected = os.getcwd()
    patient_dir = tmp_path / "patient"
    patient_dir.mkdir()
    cram = patient_dir / "in.cram"
    cram.touch()

    run_root = tmp_path / "run"
    harness = run_pipeline_under_harness(run_root / "out", bam=None, cram=str(cram))

    assert _cwd_of(harness, "validate_bam_file") == expected


def test_validate_bam_file_receives_the_configured_samtools(tmp_path: Path) -> None:
    """F1: the pipeline binds config's samtools path into the validator callback."""
    config = copy.deepcopy(MINIMAL_CONFIG)
    config["tools"]["samtools"] = "/opt/conda/envs/vntyper/bin/samtools"

    harness = run_pipeline_under_harness(tmp_path / "out", config=config)

    assert harness.kwargs("validate_bam_file")["samtools_path"] == "/opt/conda/envs/vntyper/bin/samtools"


# ---------------------------------------------------------------------------
# The quickcheck log destination (#201, #162)
#
# ``cwd`` is not the only path this call site has to get right. The quickcheck
# log used to be derived from the *input* alignment, and ``run_command`` opens
# its log before it spawns anything, so a read-only input mount failed every BAM
# and CRAM run before quickcheck even executed. The pipeline now names its own
# output directory, and both branches have to -- the BAM one being fixed and the
# CRAM one left alone would still fail every CRAM run.
# ---------------------------------------------------------------------------


def _log_dir_of(harness) -> Path:
    """Return the ``log_dir`` ``validate_bam_file`` was called with.

    Args:
        harness: The harness returned by ``run_pipeline_under_harness``.

    Returns:
        Path: The recorded ``log_dir`` argument.

    Raises:
        AssertionError: If the call carried no ``log_dir``.
    """
    kwargs = harness.kwargs("validate_bam_file")
    assert "log_dir" in kwargs, "validate_bam_file() was called without log_dir=; the log lands beside the input (#201)"
    return Path(kwargs["log_dir"])


def test_a_bam_input_validates_with_the_log_dir_set_to_the_output_directory(tmp_path: Path) -> None:
    """The BAM branch names the run's output directory as the log destination.

    Args:
        tmp_path: Pytest temporary directory.
    """
    out = tmp_path / "out"
    bam = tmp_path / "in" / "sample.bam"
    bam.parent.mkdir()
    bam.touch()

    harness = run_pipeline_under_harness(out, bam=str(bam))

    assert _log_dir_of(harness) == out
    assert _log_dir_of(harness) != bam.parent


def test_a_cram_input_validates_with_the_log_dir_set_to_the_output_directory(tmp_path: Path) -> None:
    """The CRAM branch is a second call site and gets the same argument.

    Args:
        tmp_path: Pytest temporary directory.
    """
    out = tmp_path / "out"
    cram = tmp_path / "in" / "sample.cram"
    cram.parent.mkdir()
    cram.touch()

    harness = run_pipeline_under_harness(out, bam=None, cram=str(cram))

    assert _log_dir_of(harness) == out
    assert _log_dir_of(harness) != cram.parent
