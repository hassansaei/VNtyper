"""Orchestration tests for the alignment preflight boundary in ``run_pipeline``."""

from __future__ import annotations

import json
import os
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import pipeline, pipeline_alignment
from vntyper.scripts.alignment_binding import AlignmentBinding
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


def test_owned_alignment_preflight_boundary_immediately_follows_validation_and_precedes_all_stages(
    tmp_path: Path,
) -> None:
    """Header and target preparation belong inside the post-validation seam."""
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
        mock.patch.object(
            pipeline_alignment,
            "pin_reference_resolution",
            side_effect=lambda config: events.append("pin"),
            create=True,
        ),
        mock.patch.object(pipeline, "restore_reference_resolution", create=True),
        mock.patch.object(pipeline_alignment, "restore_reference_resolution", create=True),
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
    assert events[:7] == ["validate", "header", "region", "region", "preflight", "tools", "process"]
    preflight_kwargs = harness.kwargs("run_preflight")
    process_kwargs = harness.kwargs("process_bam_to_fastq")
    assert preflight_kwargs["bed_file"] == process_kwargs["bed_file"]
    assert preflight_kwargs["region"] is None
    assert process_kwargs["plan"] is plan


def test_coverage_region_is_resolved_once_before_preflight_and_reused_by_the_consumer(tmp_path: Path) -> None:
    """The later depth stage must consume the exact region whose reference was proven."""
    out = tmp_path / "out"
    events: list[str] = []

    def resolve_region(*args: object, **kwargs: object) -> str:
        del args
        region_type = str(kwargs["region_type"])
        events.append(region_type)
        return "chr1:10-20" if region_type == "bam_region" else "chr2:30-40"

    def record(event: str) -> object:
        events.append(event)
        return mock.DEFAULT

    harness = run_pipeline_under_harness(
        out,
        stage_side_effects={
            "get_region_string_with_fallback": resolve_region,
            "run_preflight": lambda *args, **kwargs: record("preflight"),
            "get_tool_versions": lambda *args, **kwargs: record("tools"),
        },
    )

    assert harness.error is None
    assert events == ["bam_region", "vntr_region", "preflight", "tools"]
    assert harness.kwargs("run_preflight")["coverage_region"] == "chr2:30-40"
    assert harness.kwargs("calculate_vntr_coverage")["region"] == "chr2:30-40"


def test_alignment_binding_is_released_after_coverage_and_before_archiving(tmp_path: Path) -> None:
    """The descriptor and exact owned view are gone before archiving can reuse the FD."""
    out = tmp_path / "out"
    input_dir = tmp_path / "patient-input"
    input_dir.mkdir()
    alignment = input_dir / "patient.bam"
    alignment.write_bytes(b"patient alignment")
    bindings: list[AlignmentBinding] = []
    open_state_at_archive: list[bool] = []
    view_state_at_archive: list[bool] = []

    def return_bound_plan(*args: object, **kwargs: object) -> AlignmentPlan:
        del args
        binding = kwargs["binding"]
        assert isinstance(binding, AlignmentBinding)
        bindings.append(binding)
        return replace(
            _plan(out, "bam"),
            input_path=str(alignment),
            view_path=str(kwargs["bound_view_path"]),
            binding=binding,
        )

    def observe_archive(*args: object, **kwargs: object) -> str:
        del args, kwargs
        binding = bindings[0]
        open_state_at_archive.append(binding.is_open)
        view_state_at_archive.append(os.path.lexists(binding.view_path or ""))
        return str(tmp_path / "out.zip")

    with mock.patch.object(pipeline, "create_safe_archive", side_effect=observe_archive):
        harness = run_pipeline_under_harness(
            out,
            archive_results=True,
            stage_side_effects={"run_preflight": return_bound_plan},
            bam=str(alignment),
        )

    assert harness.error is None
    assert open_state_at_archive == [False]
    assert view_state_at_archive == [False]


@pytest.mark.parametrize("input_type", ["BAM", "FASTQ"])
def test_non_cram_runs_neither_read_nor_mutate_cram_ref_path_policy(
    input_type: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """CRAM-only environment policy cannot abort unrelated input modes."""
    out = tmp_path / "output"
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    selected: dict[str, Any]
    if input_type == "BAM":
        alignment = inputs / "input.bam"
        alignment.write_bytes(b"BAM")
        selected = {"bam": str(alignment)}
    else:
        fastq = inputs / "input.fastq.gz"
        fastq.write_bytes(b"FASTQ")
        selected = {"bam": None, "fastq1": str(fastq)}
    config = deepcopy(MINIMAL_CONFIG)
    config["cram"] = {
        "allow_ambient_reference_resolution": "false",
        "local_ref_path": "https://refget.invalid/%s",
    }
    original = "https://operator.invalid/%s"
    monkeypatch.setenv("REF_PATH", original)
    real_pin = pipeline.pin_reference_resolution
    real_restore = pipeline.restore_reference_resolution

    with (
        mock.patch.object(pipeline, "pin_reference_resolution", wraps=real_pin) as pin,
        mock.patch.object(pipeline, "restore_reference_resolution", wraps=real_restore) as restore,
    ):
        harness = run_pipeline_under_harness(out, config=config, expect_failure=True, **selected)

    assert harness.error is None
    pin.assert_not_called()
    restore.assert_not_called()
    assert os.environ["REF_PATH"] == original


@pytest.mark.parametrize("failure", [ValueError("target resolution failed"), KeyboardInterrupt("interrupted")])
def test_alignment_boundary_restores_ref_path_when_target_preparation_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: BaseException,
) -> None:
    """A failure after pinning cannot leak run-local REF_PATH into the process."""
    out = tmp_path / "out"
    original = "http://operator.example/%s"
    pinned = "/local/cache/%2s/%2s/%s"
    observed_ref_paths: list[str | None] = []
    monkeypatch.setenv("REF_PATH", original)
    config = {**MINIMAL_CONFIG, "cram": {"local_ref_path": pinned}}
    cram = tmp_path / "patient-input" / "input.cram"
    cram.parent.mkdir()
    cram.touch()
    out.mkdir(parents=True)
    artifact = out / "preflight_error.json"
    artifact.write_text('{"code":"stale"}\n', encoding="utf-8")

    def fail_target_resolution(*args: object, **kwargs: object) -> str:
        del args, kwargs
        observed_ref_paths.append(os.environ.get("REF_PATH"))
        raise failure

    harness = run_pipeline_under_harness(
        out,
        config=config,
        bam=None,
        cram=str(cram),
        expect_failure=True,
        stage_side_effects={"get_region_string_with_fallback": fail_target_resolution},
    )

    assert isinstance(harness.error, SystemExit if isinstance(failure, Exception) else KeyboardInterrupt)
    assert observed_ref_paths == [pinned]
    assert os.environ["REF_PATH"] == original
    assert not harness.stages["run_preflight"].called
    assert not harness.stages["get_tool_versions"].called
    if isinstance(failure, Exception):
        assert json.loads(artifact.read_text(encoding="utf-8")) == {
            "code": "alignment_target_invalid",
            "message": "Alignment preflight could not prepare the requested target; verify the BED file or "
            "configured regions.",
            "candidates": [],
        }
    else:
        assert not artifact.exists()


def test_fresh_output_root_exists_for_validation_but_stage_directories_wait_for_preflight(tmp_path: Path) -> None:
    out = tmp_path / "run-output" / "fresh-output"
    bam = tmp_path / "patient-input" / "patient.bam"
    bam.parent.mkdir()
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


def test_pipeline_forwards_a_cram_header_to_remote_uri_policy_before_region_or_preflight(tmp_path: Path) -> None:
    """The top-level orchestration cannot bypass the owned post-header policy check."""
    out = tmp_path / "run-output"
    cram = tmp_path / "patient-input" / "patient.cram"
    cram.parent.mkdir()
    cram.touch()
    header = "@SQ\tSN:chr1\tLN:100\tUR:https://refget.example/private/reference.fa\n"

    harness = run_pipeline_under_harness(
        out,
        bam=None,
        cram=str(cram),
        config=MINIMAL_CONFIG,
        expect_failure=True,
        stage_side_effects={"read_alignment_header": lambda *args, **kwargs: header},
    )

    assert isinstance(harness.error, SystemExit)
    assert json.loads((out / "preflight_error.json").read_text(encoding="utf-8"))["code"] == (
        "reference_policy_invalid"
    )
    assert not harness.stages["get_region_string_with_fallback"].called
    assert not harness.stages["run_preflight"].called
    assert not harness.stages["get_tool_versions"].called


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
    out = tmp_path / "run-output"
    cram = tmp_path / "patient-input" / "patient.cram"
    cram.parent.mkdir()
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
    assert harness.kwargs("run_kestrel")["fastq_files"] == harness.stages["route_converted_fastqs"].return_value


def test_pipeline_threads_the_explicit_cram_reference_to_preflight(tmp_path: Path) -> None:
    """The parser value must survive pipeline orchestration without being dropped.

    Args:
        tmp_path: Pytest temporary directory.
    """
    out = tmp_path / "run-output"
    cram = tmp_path / "patient-input" / "patient.cram"
    cram.parent.mkdir()
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
    cram = tmp_path / "patient-input" / "input.cram"
    cram.parent.mkdir()
    cram.touch()
    plan = _plan(out, "cram", reference_path="/refs/hg19.fa")

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
        bam=None,
        cram=str(cram),
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


@pytest.mark.parametrize("initial_ref_path", [None, "http://operator.example/%s"])
@pytest.mark.parametrize("outcome", ["return", "system_exit", "exception"])
def test_cleanup_failure_never_masks_a_primary_outcome_or_skips_exact_ref_path_restoration(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    initial_ref_path: str | None,
    outcome: str,
) -> None:
    """CRAM cleanup failure preserves a primary outcome and always restores REF_PATH."""
    out = tmp_path / "out"
    pinned = "/local/cache/%2s/%2s/%s"
    if initial_ref_path is None:
        monkeypatch.delenv("REF_PATH", raising=False)
    else:
        monkeypatch.setenv("REF_PATH", initial_ref_path)
    config = {**MINIMAL_CONFIG, "cram": {"local_ref_path": pinned}}
    cram = tmp_path / "patient-input" / "input.cram"
    cram.parent.mkdir()
    cram.touch()
    binding = AlignmentBinding(str(cram))
    plan = replace(_plan(out, "cram", reference_path="/refs/hg19.fa"), input_path=str(cram), binding=binding)
    stage_side_effects: dict[str, Any] = {"run_preflight": lambda *args, **kwargs: plan}
    if outcome == "system_exit":
        stage_side_effects["process_bam_to_fastq"] = SystemExit(7)
    elif outcome == "exception":
        stage_side_effects["process_bam_to_fastq"] = RuntimeError("stage failed")
    cleanup_effect: object = (
        [None, RuntimeError("alignment cleanup failed")]
        if outcome == "return"
        else RuntimeError("alignment cleanup failed")
    )

    try:
        with mock.patch.object(binding, "close", side_effect=cleanup_effect):
            harness = run_pipeline_under_harness(
                out,
                config=config,
                bam=None,
                cram=str(cram),
                expect_failure=True,
                stage_side_effects=stage_side_effects,
            )

        if outcome == "return":
            assert isinstance(harness.error, RuntimeError)
            assert str(harness.error) == "alignment cleanup failed"
        elif outcome == "system_exit":
            assert isinstance(harness.error, SystemExit)
            assert harness.error.code == 7
        else:
            assert isinstance(harness.error, SystemExit)
            assert harness.error.code == 1
            assert isinstance(harness.error.__context__, RuntimeError)
            assert str(harness.error.__context__) == "stage failed"
        if initial_ref_path is None:
            assert "REF_PATH" not in os.environ
        else:
            assert os.environ["REF_PATH"] == initial_ref_path
    finally:
        binding.close()
        if initial_ref_path is None:
            monkeypatch.delenv("REF_PATH", raising=False)
        else:
            monkeypatch.setenv("REF_PATH", initial_ref_path)


def test_caller_exception_state_does_not_hide_successful_run_cleanup_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Only this run's outcome can decide whether its cleanup failure is suppressed."""
    out = tmp_path / "out"
    original_ref_path = "http://operator.example/%s"
    pinned_ref_path = "/local/cache/%2s/%2s/%s"
    monkeypatch.setenv("REF_PATH", original_ref_path)
    config = {**MINIMAL_CONFIG, "cram": {"local_ref_path": pinned_ref_path}}
    cram = tmp_path / "patient-input" / "input.cram"
    cram.parent.mkdir()
    cram.touch()
    binding = AlignmentBinding(str(cram))
    plan = replace(_plan(out, "cram", reference_path="/refs/hg19.fa"), input_path=str(cram), binding=binding)

    try:
        try:
            raise LookupError("unrelated caller failure")
        except LookupError:
            with mock.patch.object(
                binding,
                "close",
                side_effect=[None, RuntimeError("alignment cleanup failed")],
            ):
                harness = run_pipeline_under_harness(
                    out,
                    config=config,
                    bam=None,
                    cram=str(cram),
                    expect_failure=True,
                    stage_side_effects={"run_preflight": lambda *args, **kwargs: plan},
                )

        assert isinstance(harness.error, RuntimeError)
        assert str(harness.error) == "alignment cleanup failed"
        assert os.environ["REF_PATH"] == original_ref_path
    finally:
        binding.close()


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
