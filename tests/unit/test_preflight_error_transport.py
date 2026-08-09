"""Public failure artifacts emitted by the complete preflight boundary."""

from __future__ import annotations

import json
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_preflight import run_preflight
from vntyper.scripts.preflight_error_io import PreflightErrorContext, persist_preflight_failure

pytestmark = pytest.mark.unit


def _alignment(tmp_path: Path, file_format: str, *, indexed: bool = True) -> Path:
    """Create a minimal alignment and, by default, its conventional index."""
    input_dir = tmp_path / "private-worker-input"
    input_dir.mkdir()
    alignment = input_dir / f"patient.{file_format}"
    alignment.write_bytes(b"alignment")
    if indexed:
        suffix = "bai" if file_format == "bam" else "crai"
        (input_dir / f"patient.{file_format}.{suffix}").write_bytes(b"index")
    return alignment


def _artifact(run_root: Path) -> dict:
    """Load the exact artifact installed before the preflight exception escapes."""
    return json.loads((run_root / "preflight_error.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("failure", "expected_type", "original_fragment", "code", "public_message", "candidate_count"),
    [
        (
            "unsafe_output",
            ValueError,
            "single basename",
            "preflight_output_unsafe",
            "Alignment preflight rejected an unsafe output or log destination; use a dedicated output directory "
            "and remove conflicting entries.",
            0,
        ),
        (
            "missing_index",
            RuntimeError,
            "samtools index",
            "alignment_index_unusable",
            "Alignment preflight could not prepare a safe local view and fresh index; remove conflicting "
            "run-output entries and verify samtools can index the alignment.",
            0,
        ),
        (
            "candidate_policy",
            ValueError,
            "exactly one terminal",
            "reference_policy_invalid",
            "CRAM reference candidate policy is invalid; configure a list ending in exactly one terminal "
            "htslib-resolved candidate.",
            0,
        ),
        (
            "scan_policy",
            ValueError,
            "invalid unmapped scan mode",
            "unmapped_scan_invalid",
            "Alignment unmapped-read scan selection failed; use auto or stream mode to avoid losing "
            "placed-unmapped reads.",
            0,
        ),
        (
            "forced_indexed",
            ValueError,
            "3 placed-unmapped reads",
            "unmapped_scan_invalid",
            "Alignment unmapped-read scan selection failed; use auto or stream mode to avoid losing "
            "placed-unmapped reads.",
            0,
        ),
        (
            "bam_probe",
            RuntimeError,
            "private target",
            "alignment_probe_failed",
            "BAM preflight could not retrieve the requested target; verify the index and target coordinates.",
            0,
        ),
        (
            "reference",
            ValueError,
            "contig=chr1",
            "reference_unresolved",
            "Unable to resolve reference for CRAM input: contig=chr1, M5=digest. Candidates: source=cli, "
            "path=no path, reason=not supplied; source=config_cram_reference, path=no path, reason=not supplied; "
            "source=config_bwa_reference, path=no path, reason=not supplied; source=htslib-resolved "
            "(header UR: or REF_PATH), path=no path, reason=probe exited non-zero",
            4,
        ),
    ],
)
def test_every_actionable_preflight_failure_writes_one_curated_artifact_before_reraising(
    tmp_path: Path,
    failure: str,
    expected_type: type[Exception],
    original_fragment: str,
    code: str,
    public_message: str,
    candidate_count: int,
) -> None:
    """Each preflight stage crosses the same path-free three-field boundary."""
    run_root = tmp_path / "run-output"
    run_root.mkdir()
    work_dir = run_root / "alignment"
    file_format = "bam" if failure in {"unsafe_output", "missing_index", "bam_probe"} else "cram"
    alignment = _alignment(tmp_path, file_format, indexed=failure != "missing_index")
    output_name = "../private-view" if failure == "unsafe_output" else "sample"
    config: dict = {}
    command_results: list[tuple[bool, str]] = []
    binding: AlignmentBinding | None = None
    if failure == "candidate_policy":
        config = {"cram": {"reference_candidate_order": ["cli"]}}
    elif failure == "scan_policy":
        config = {"cram": {"unmapped_scan": "/private/invalid-mode"}}
        command_results = [(True, "chr1\t4\t1\t0\n*\t0\t0\t2\n")]
    elif failure == "forced_indexed":
        config = {"cram": {"unmapped_scan": "indexed"}}
        command_results = [(True, "chr1\t4\t1\t3\n*\t0\t0\t2\n")]
    elif failure == "reference":
        command_results = [
            (True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"),
            (False, "failed to decode /private/worker/patient.cram"),
        ]
    elif failure == "missing_index":
        command_results = [(False, "cannot build /private/worker/patient.bam")]
    elif failure == "bam_probe":
        command_results = [
            (True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"),
            (False, "private target /private/worker/patient.bam is stale"),
        ]

    with ExitStack() as stack:
        if failure in {"scan_policy", "forced_indexed", "bam_probe", "reference"}:
            view = work_dir / f"sample.{file_format}"
            index = Path(f"{view}.{'bai' if file_format == 'bam' else 'crai'}")
            binding = AlignmentBinding(str(alignment))
            stack.enter_context(
                patch(
                    "vntyper.scripts.alignment_preflight.build_alignment_view",
                    return_value=(str(view), str(index), binding),
                )
            )
        if command_results:
            stack.enter_context(
                patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=command_results)
            )
        raised = stack.enter_context(pytest.raises(expected_type))
        run_preflight(
            str(alignment),
            str(work_dir),
            output_name,
            file_format,
            config,
            1,
            region="chr1:1-2",
            header_contigs=("chr1",),
            m5="digest",
            error_output_dir=run_root,
        )

    assert original_fragment in str(raised.value)
    payload = _artifact(run_root)
    assert set(payload) == {"code", "message", "candidates"}
    assert payload["code"] == code
    assert payload["message"] == public_message
    assert len(payload["candidates"]) == candidate_count
    serialized = json.dumps(payload)
    assert "/private/" not in serialized
    assert "\\private\\" not in serialized
    if binding is not None:
        assert binding.is_open is False


@pytest.mark.parametrize(
    "candidate_order",
    ["cli", ["cli"]],
)
def test_reference_policy_type_and_order_share_one_stable_public_failure(
    tmp_path: Path,
    candidate_order: str | list[str],
) -> None:
    """Policy transport is stable across validation prose and never copies the invalid value."""
    run_root = tmp_path / "run-output"
    run_root.mkdir()
    alignment = _alignment(tmp_path, "cram")

    with pytest.raises(ValueError):
        run_preflight(
            str(alignment),
            str(run_root / "alignment"),
            "sample",
            "cram",
            {"cram": {"reference_candidate_order": candidate_order}},
            1,
            region="chr1:1-2",
            error_output_dir=run_root,
        )

    assert _artifact(run_root) == {
        "code": "reference_policy_invalid",
        "message": "CRAM reference candidate policy is invalid; configure a list ending in exactly one terminal "
        "htslib-resolved candidate.",
        "candidates": [],
    }


@pytest.mark.parametrize(
    "original",
    [
        ValueError("Regular reserved index has no generated-index provenance: /private/view.bai"),
        ValueError("Unable to inspect generated-index provenance /private/view.bai: private failure"),
        RuntimeError("Unable to allocate a temporary symlink beside /private/view.bai"),
    ],
)
def test_every_view_index_provenance_exception_uses_the_same_stable_phase_payload(
    tmp_path: Path,
    original: Exception,
) -> None:
    """View/index internals cannot change the public code or leak their destination."""
    run_root = tmp_path / "run-output"
    run_root.mkdir()
    alignment = _alignment(tmp_path, "bam")

    with (
        patch("vntyper.scripts.alignment_preflight.build_alignment_view", side_effect=original),
        pytest.raises(type(original)) as raised,
    ):
        run_preflight(
            str(alignment),
            str(run_root / "alignment"),
            "sample",
            "bam",
            {},
            1,
            region="chr1:1-2",
            error_output_dir=run_root,
        )

    assert raised.value is original
    assert _artifact(run_root) == {
        "code": "alignment_index_unusable",
        "message": "Alignment preflight could not prepare a safe local view and fresh index; remove conflicting "
        "run-output entries and verify samtools can index the alignment.",
        "candidates": [],
    }


def test_an_unexpected_reference_probe_exception_uses_the_reference_phase_payload(tmp_path: Path) -> None:
    """Reference probe internals cannot fall back to scan semantics or expose their prose."""
    run_root = tmp_path / "run-output"
    run_root.mkdir()
    alignment = _alignment(tmp_path, "cram")
    original = RuntimeError("unexpected reference probe /private/worker failure")
    binding = AlignmentBinding(str(alignment))

    with (
        patch(
            "vntyper.scripts.alignment_preflight.build_alignment_view",
            return_value=(
                str(run_root / "alignment" / "sample.cram"),
                str(run_root / "alignment" / "sample.cram.crai"),
                binding,
            ),
        ),
        patch(
            "vntyper.scripts.alignment_preflight.capture_command",
            return_value=(True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"),
        ),
        patch("vntyper.scripts.alignment_preflight.resolve_reference", side_effect=original),
        pytest.raises(RuntimeError) as raised,
    ):
        run_preflight(
            str(alignment),
            str(run_root / "alignment"),
            "sample",
            "cram",
            {},
            1,
            region="chr1:1-2",
            error_output_dir=run_root,
        )

    assert raised.value is original
    assert binding.is_open is False
    assert _artifact(run_root) == {
        "code": "reference_unresolved",
        "message": "CRAM reference preflight could not prove the requested target; verify the reference FASTA and "
        "target contigs.",
        "candidates": [],
    }


def test_an_unclassified_boundary_exception_uses_only_the_generic_public_failure(tmp_path: Path) -> None:
    """Unexpected internals retain their identity while exposing no exception prose."""
    original = RuntimeError("unexpected /private/worker detail")
    context = PreflightErrorContext(tmp_path)

    with pytest.raises(RuntimeError) as raised, persist_preflight_failure(context):
        raise original

    assert raised.value is original
    assert _artifact(tmp_path) == {
        "code": "alignment_preflight_failed",
        "message": "Alignment preflight failed before processing; inspect the server logs for the job.",
        "candidates": [],
    }


def test_artifact_write_failure_does_not_mask_the_original_preflight_exception(tmp_path: Path) -> None:
    """The primary diagnostic survives every ordinary persistence failure."""
    run_root = tmp_path / "run-output"
    run_root.mkdir()
    alignment = _alignment(tmp_path, "bam")
    original = RuntimeError("BAM preflight probe failed: stale index")

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=original),
        patch(
            "vntyper.scripts.preflight_error_io.write_preflight_error",
            side_effect=RuntimeError("artifact storage unavailable"),
        ) as writer,
        pytest.raises(RuntimeError) as raised,
    ):
        run_preflight(
            str(alignment),
            str(run_root / "alignment"),
            "sample",
            "bam",
            {},
            1,
            region="chr1:1-2",
            error_output_dir=run_root,
        )

    assert raised.value is original
    writer.assert_called_once()


def test_preflight_boundary_does_not_catch_base_exceptions(tmp_path: Path) -> None:
    """Process-control exceptions bypass both serialization and persistence."""
    run_root = tmp_path / "run-output"
    run_root.mkdir()
    alignment = _alignment(tmp_path, "bam")

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=KeyboardInterrupt),
        patch("vntyper.scripts.preflight_error_io.write_preflight_error") as writer,
        pytest.raises(KeyboardInterrupt),
    ):
        run_preflight(
            str(alignment),
            str(run_root / "alignment"),
            "sample",
            "bam",
            {},
            1,
            region="chr1:1-2",
            error_output_dir=run_root,
        )

    writer.assert_not_called()


def test_shared_failure_context_is_persisted_only_by_its_outer_owner(tmp_path: Path) -> None:
    """Passing a context suppresses run_preflight's otherwise self-owned persistence."""
    run_root = tmp_path / "run-output"
    run_root.mkdir()
    alignment = _alignment(tmp_path, "bam")
    context = PreflightErrorContext(run_root)
    original = RuntimeError("index preparation failed")

    with (
        patch("vntyper.scripts.alignment_preflight.build_alignment_view", side_effect=original),
        patch("vntyper.scripts.preflight_error_io.write_preflight_error") as writer,
        pytest.raises(RuntimeError) as raised,
        persist_preflight_failure(context),
    ):
        run_preflight(
            str(alignment),
            str(run_root / "alignment"),
            "sample",
            "bam",
            {},
            1,
            region="chr1:1-2",
            error_output_dir=run_root,
            failure_context=context,
        )

    assert raised.value is original
    writer.assert_called_once()
