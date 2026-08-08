"""Orchestration tests for the alignment preflight boundary in ``run_pipeline``."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import pipeline
from vntyper.scripts.alignment_contract import AlignmentPlan

pytestmark = pytest.mark.unit


def _plan(output_dir: Path, file_format: str, *, reference_path: str | None = None) -> AlignmentPlan:
    """Build the exact object a mocked successful preflight returns.

    Args:
        output_dir: Pipeline output directory.
        file_format: Alignment format proven by preflight.
        reference_path: Proven CRAM reference, when required.

    Returns:
        A frozen alignment plan rooted in the run directory.
    """
    stage_dir = output_dir / "fastq_bam_processing"
    return AlignmentPlan(
        input_path=str(output_dir / f"in.{file_format}"),
        view_path=str(stage_dir / f"input.{file_format}"),
        file_format=file_format,
        index_path=str(stage_dir / f"input.{file_format}.{'bai' if file_format == 'bam' else 'crai'}"),
        reference_path=reference_path,
        reference_source="test",
        uncovered_contigs=(),
        unmapped_scan="indexed",
    )


def test_alignment_validation_and_target_resolution_precede_preflight_and_all_stages(tmp_path: Path) -> None:
    out = tmp_path / "out"
    plan = _plan(out, "bam")
    events: list[str] = []

    def preflight(*args, **kwargs):
        events.append("preflight")
        return plan

    def record(event: str) -> object:
        events.append(event)
        return mock.DEFAULT

    with (
        mock.patch.object(
            pipeline, "pin_reference_resolution", side_effect=lambda config: events.append("pin"), create=True
        ),
        mock.patch.object(pipeline, "restore_reference_resolution", create=True),
    ):
        harness = run_pipeline_under_harness(
            out,
            stage_side_effects={
                "validate_bam_file": lambda *args, **kwargs: events.append("validate"),
                "read_alignment_header": lambda *args, **kwargs: record("header"),
                "get_region_string_with_fallback": lambda *args, **kwargs: record("region"),
                "run_preflight": preflight,
                "get_tool_versions": lambda *args, **kwargs: record("tools"),
                "process_bam_to_fastq": lambda *args, **kwargs: record("process"),
            },
        )

    assert harness.error is None
    assert events.index("validate") < events.index("header") < events.index("region")
    assert events.index("region") < events.index("pin") < events.index("preflight")
    assert events.index("preflight") < events.index("tools") < events.index("process")
    preflight_kwargs = harness.kwargs("run_preflight")
    process_kwargs = harness.kwargs("process_bam_to_fastq")
    assert preflight_kwargs["bed_file"] == process_kwargs["bed_file"]
    assert preflight_kwargs["region"] is None
    assert process_kwargs["plan"] is plan


def test_fresh_output_root_exists_for_validation_but_stage_directories_wait_for_preflight(tmp_path: Path) -> None:
    out = tmp_path / "fresh-output"
    bam = tmp_path / "patient.bam"
    bam.touch()

    def validate(*args: Any, **kwargs: Any) -> None:
        assert out.is_dir()
        assert not (out / "kestrel").exists()
        assert not (out / "coverage").exists()

    harness = run_pipeline_under_harness(
        out,
        create_output_dir=False,
        bam=str(bam),
        stage_side_effects={"validate_bam_file": validate},
    )

    assert harness.error is None
    assert harness.stages["run_preflight"].called


def test_patient_input_tree_cannot_be_used_as_the_output_root(tmp_path: Path) -> None:
    input_root = tmp_path / "patient-input"
    input_root.mkdir()
    bam = input_root / "patient.bam"
    bam.write_bytes(b"patient-bytes")

    harness = run_pipeline_under_harness(
        input_root,
        create_output_dir=False,
        bam=str(bam),
        expect_failure=True,
    )

    assert isinstance(harness.error, SystemExit)
    assert set(input_root.iterdir()) == {bam}
    assert bam.read_bytes() == b"patient-bytes"
    assert not harness.stages["validate_bam_file"].called
    assert not harness.stages["run_preflight"].called


def test_cram_plan_drives_conversion_and_reference_aware_coverage(tmp_path: Path) -> None:
    out = tmp_path / "out"
    cram = tmp_path / "patient.cram"
    cram.touch()
    reference = str(tmp_path / "reference genome.fa")
    plan = _plan(out, "cram", reference_path=reference)
    plan = AlignmentPlan(**{**plan.__dict__, "input_path": str(cram)})

    with (
        mock.patch.object(pipeline, "pin_reference_resolution", return_value=None, create=True),
        mock.patch.object(pipeline, "restore_reference_resolution", create=True),
    ):
        harness = run_pipeline_under_harness(
            out,
            bam=None,
            cram=str(cram),
            stage_side_effects={"run_preflight": lambda *args, **kwargs: plan},
        )

    assert harness.kwargs("process_bam_to_fastq")["plan"] is plan
    coverage = harness.kwargs("calculate_vntr_coverage")
    assert coverage["bam_file"] == plan.view_path
    assert coverage["reference_path"] == reference


def test_explicit_cram_reference_reaches_preflight_as_the_cli_candidate(tmp_path: Path) -> None:
    """Dropping this value would let config or htslib outrank the CLI flag.

    Args:
        tmp_path: Pytest temporary directory.
    """
    out = tmp_path / "out"
    cram = tmp_path / "patient.cram"
    cram.touch()
    reference = tmp_path / "full reference.fa"

    with (
        mock.patch.object(pipeline, "pin_reference_resolution", return_value=None, create=True),
        mock.patch.object(pipeline, "restore_reference_resolution", create=True),
    ):
        harness = run_pipeline_under_harness(
            out,
            bam=None,
            cram=str(cram),
            reference_fasta=reference,
        )

    assert harness.kwargs("run_preflight")["reference_fasta"] == str(reference)
    assert harness.kwargs("run_preflight")["error_output_dir"] == str(out)


def test_post_alignment_bam_is_preflighted_and_its_returned_plan_is_consumed(tmp_path: Path) -> None:
    out = tmp_path / "out"
    plan = _plan(out, "bam")
    events: list[str] = []

    def preflight(*args, **kwargs):
        events.append("preflight")
        return plan

    def record(event: str) -> object:
        events.append(event)
        return mock.DEFAULT

    with (
        mock.patch.object(pipeline, "pin_reference_resolution", return_value=None, create=True),
        mock.patch.object(pipeline, "restore_reference_resolution", create=True),
    ):
        harness = run_pipeline_under_harness(
            out,
            bam=None,
            fastq1="/in/r1.fastq.gz",
            fastq2="/in/r2.fastq.gz",
            stage_side_effects={
                "align_and_sort_fastq": lambda *args, **kwargs: record("align"),
                "run_preflight": preflight,
                "process_bam_to_fastq": lambda *args, **kwargs: record("process"),
            },
        )

    assert events == ["align", "preflight", "process"]
    preflight_kwargs = harness.kwargs("run_preflight")
    assert preflight_kwargs["in_path"] == harness.stages["align_and_sort_fastq"].return_value
    assert preflight_kwargs["file_format"] == "bam"
    assert harness.kwargs("process_bam_to_fastq")["plan"] is plan
    assert harness.kwargs("calculate_vntr_coverage")["bam_file"] == plan.view_path


@pytest.mark.parametrize("initial_ref_path", [None, "http://operator.example/%s"])
@pytest.mark.parametrize("outcome", ["return", "system_exit", "base_exception"])
def test_ref_path_is_pinned_through_the_run_and_restored_exactly(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    initial_ref_path: str | None,
    outcome: str,
) -> None:
    out = tmp_path / "out"
    pinned = "/local/cache/%2s/%2s/%s"
    if initial_ref_path is None:
        monkeypatch.delenv("REF_PATH", raising=False)
    else:
        monkeypatch.setenv("REF_PATH", initial_ref_path)
    config = {**MINIMAL_CONFIG, "cram": {"local_ref_path": pinned}}
    plan = _plan(out, "bam")

    def observe_pin(*args, **kwargs):
        assert os.environ["REF_PATH"] == pinned
        return mock.DEFAULT

    stage_side_effects: dict[str, Any] = {"process_bam_to_fastq": observe_pin}
    if outcome == "system_exit":
        stage_side_effects["run_kestrel"] = SystemExit(7)
    elif outcome == "base_exception":
        stage_side_effects["run_kestrel"] = KeyboardInterrupt("operator interrupted")

    stage_side_effects["run_preflight"] = lambda *args, **kwargs: plan
    harness = run_pipeline_under_harness(
        out,
        config=config,
        expect_failure=outcome != "return",
        stage_side_effects=stage_side_effects,
    )

    if outcome == "system_exit":
        assert isinstance(harness.error, SystemExit)
        assert harness.error.code == 7
    elif outcome == "base_exception":
        assert isinstance(harness.error, KeyboardInterrupt)
    else:
        assert harness.error is None
    if initial_ref_path is None:
        assert "REF_PATH" not in os.environ
    else:
        assert os.environ["REF_PATH"] == initial_ref_path


def test_preflight_failure_stops_before_any_processing_stage(tmp_path: Path) -> None:
    out = tmp_path / "out"

    harness = run_pipeline_under_harness(
        out,
        expect_failure=True,
        stage_side_effects={"run_preflight": RuntimeError("probe failed")},
    )

    assert isinstance(harness.error, SystemExit)
    assert not harness.stages["get_tool_versions"].called
    assert not harness.stages["process_bam_to_fastq"].called
    assert not harness.stages["run_kestrel"].called
