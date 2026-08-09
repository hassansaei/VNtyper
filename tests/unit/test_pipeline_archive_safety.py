"""Archive safety at the completed CLI pipeline boundary."""

import zipfile
from pathlib import Path

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness

pytestmark = pytest.mark.unit

PATIENT_BYTES = b"patient-alignment-bytes-that-must-never-enter-an-archive"


def test_cli_archive_refuses_a_run_local_alignment_symlink_before_install(tmp_path: Path, caplog) -> None:
    """A live preflight view must fail the completed run rather than leak its target.

    Args:
        tmp_path: Scratch directory standing in for the input and output trees.
        caplog: Pytest log capture.
    """
    patient_alignment = tmp_path / "patient" / "patient.bam"
    patient_alignment.parent.mkdir()
    patient_alignment.write_bytes(PATIENT_BYTES)
    output_dir = tmp_path / "results"

    def leave_alignment_view(*args, **kwargs) -> None:
        del args, kwargs
        (output_dir / "alignment_view.bam").symlink_to(patient_alignment)

    harness = run_pipeline_under_harness(
        output_dir,
        archive_results=True,
        expect_failure=True,
        stage_side_effects={"generate_summary_report": leave_alignment_view},
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not Path(f"{output_dir}.zip").exists(), "an unsafe archive must never be installed"
    assert (output_dir / "alignment_view.bam").is_symlink(), "failure must retain the diagnosable result tree"
    assert patient_alignment.read_bytes() == PATIENT_BYTES
    assert "Pipeline finished successfully." not in caplog.messages


def test_cli_archive_still_packages_regular_results(tmp_path: Path) -> None:
    """Rejecting aliases must not disable the ordinary archive success path.

    Args:
        tmp_path: Scratch directory standing in for the output tree.
    """
    output_dir = tmp_path / "results"
    result_bytes = b"ordinary-result"

    def write_result(*args, **kwargs) -> None:
        del args, kwargs
        (output_dir / "result.txt").write_bytes(result_bytes)

    harness = run_pipeline_under_harness(
        output_dir,
        archive_results=True,
        stage_side_effects={"generate_summary_report": write_result},
    )

    assert harness.error is None
    archive_path = Path(f"{output_dir}.zip")
    with zipfile.ZipFile(archive_path) as archive:
        assert archive.read("result.txt") == result_bytes
