"""Archive safety at the completed CLI pipeline boundary."""

import os
import tarfile
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts import pipeline as pipeline_module

pytestmark = pytest.mark.unit

PATIENT_BYTES = b"patient-alignment-bytes-that-must-never-enter-an-archive"


@pytest.mark.parametrize(
    ("archive_format", "suffix"),
    [pytest.param("zip", ".zip", id="zip"), pytest.param("tar.gz", ".tar.gz", id="tar")],
)
def test_cli_releases_the_h1_owned_view_before_archiving(tmp_path: Path, archive_format: str, suffix: str) -> None:
    """The combined H1/H3 boundary closes the exact view before ZIP or tar."""
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    output_dir = tmp_path / "results"
    plans: list[SimpleNamespace] = []

    def owned_plan(*args, **kwargs):
        del args
        view = output_dir / f"input.{kwargs['file_format']}"
        view.symlink_to(patient)
        plan = SimpleNamespace(
            input_path=str(patient),
            view_path=str(view),
            file_format=kwargs["file_format"],
            index_path=f"{view}.bai",
            reference_path=None,
            reference_source="test",
            uncovered_contigs=(),
            unmapped_scan="indexed",
            close_calls=0,
        )

        def close() -> None:
            plan.close_calls += 1
            if not os.path.lexists(view):
                return
            assert os.readlink(view) == str(patient)
            view.unlink()

        plan.close = close
        plans.append(plan)
        return plan

    harness = run_pipeline_under_harness(
        output_dir,
        archive_results=True,
        archive_format=archive_format,
        stage_side_effects={"run_preflight": owned_plan},
    )

    assert harness.error is None
    assert plans and not Path(plans[0].view_path).exists()
    assert plans[0].close_calls == 2  # final coverage plus the outer idempotent fallback
    archive = Path(f"{output_dir}{suffix}")
    assert PATIENT_BYTES not in archive.read_bytes()
    if archive_format == "zip":
        with zipfile.ZipFile(archive) as zip_package:
            assert "input.bam" not in zip_package.namelist()
    else:
        with tarfile.open(archive, "r:gz") as tar_package:
            assert "./input.bam" not in tar_package.getnames()


def test_cli_archive_runs_after_all_final_summary_formats_are_written(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The direct CLI archive owns the complete, final result tree.

    Args:
        tmp_path: Scratch directory standing in for the input and output trees.
    """
    output_dir = tmp_path / "results"

    def assert_final_outputs(base_name, archive_format, root_dir, **kwargs) -> str:
        del base_name, archive_format, root_dir, kwargs
        assert (output_dir / "pipeline_summary.json").exists()
        assert (output_dir / "pipeline_summary.csv").exists()
        assert (output_dir / "pipeline_summary.tsv").exists()
        return str(output_dir) + ".zip"

    monkeypatch.setattr(pipeline_module, "create_safe_archive", assert_final_outputs)
    harness = run_pipeline_under_harness(
        output_dir,
        archive_results=True,
        summary_formats=["csv", "tsv"],
    )

    assert harness.error is None


def test_cli_malicious_result_symlink_fails_with_completed_outcome_one(tmp_path: Path) -> None:
    """An alias not owned by H1 still fails the CLI without exposing patient bytes."""
    patient = tmp_path / "external-patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    output_dir = tmp_path / "results"

    def add_malicious_alias(*args, **kwargs) -> None:
        del args, kwargs
        (output_dir / "malicious.bam").symlink_to(patient)

    harness = run_pipeline_under_harness(
        output_dir,
        archive_results=True,
        expect_failure=True,
        stage_side_effects={"generate_summary_report": add_malicious_alias},
    )

    assert isinstance(harness.error, SystemExit) and harness.error.code == 1
    assert not Path(f"{output_dir}.zip").exists()
    assert patient.read_bytes() == PATIENT_BYTES


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
