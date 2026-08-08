"""Public failure artifacts emitted by the complete preflight boundary."""

from __future__ import annotations

import json
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_preflight import run_preflight

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
    ("failure", "expected_type", "original_fragment", "code", "public_fragment", "candidate_count"),
    [
        ("unsafe_output", ValueError, "single basename", "preflight_output_unsafe", "output name", 0),
        ("missing_index", RuntimeError, "samtools index", "alignment_index_unusable", "BAM index", 0),
        (
            "candidate_policy",
            ValueError,
            "exactly one terminal",
            "reference_policy_invalid",
            "candidate policy",
            0,
        ),
        ("scan_policy", ValueError, "invalid unmapped scan mode", "unmapped_scan_invalid", "scan policy", 0),
        (
            "forced_indexed",
            ValueError,
            "3 placed-unmapped reads",
            "unmapped_scan_lossy",
            "3 placed-unmapped reads",
            0,
        ),
        ("bam_probe", RuntimeError, "private target", "alignment_probe_failed", "requested target", 0),
        ("reference", ValueError, "contig=chr1", "reference_unresolved", "contig=chr1", 4),
    ],
)
def test_every_actionable_preflight_failure_writes_one_curated_artifact_before_reraising(
    tmp_path: Path,
    failure: str,
    expected_type: type[Exception],
    original_fragment: str,
    code: str,
    public_fragment: str,
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
        command_results = [(False, "private target /private/worker/patient.bam is stale")]

    with ExitStack() as stack:
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
    assert public_fragment in payload["message"]
    assert len(payload["candidates"]) == candidate_count
    serialized = json.dumps(payload)
    assert "/private/" not in serialized
    assert "\\private\\" not in serialized


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
