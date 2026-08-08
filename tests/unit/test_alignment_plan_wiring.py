"""Wiring tests for consumers of the proven alignment preflight plan."""

from __future__ import annotations

import hashlib
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vntyper.scripts import fastq_bam_processing
from vntyper.scripts.alignment_contract import AlignmentPlan

pytestmark = pytest.mark.unit


CONFIG = {"tools": {"samtools": "samtools"}}
REGION = "chr1:1-2"


def _plan(tmp_path: Path, file_format: str, *, unmapped_scan: str = "indexed") -> AlignmentPlan:
    """Return a plan whose writable view is separate from its patient input.

    Args:
        tmp_path: Test directory holding the input and view trees.
        file_format: Alignment format for the plan.
        unmapped_scan: Proven CRAM scan selected by preflight.

    Returns:
        A complete alignment plan for the stage under test.
    """
    input_path = tmp_path / "input" / f"sample.{file_format}"
    view_path = tmp_path / "run" / f"input.{file_format}"
    input_path.parent.mkdir(exist_ok=True)
    view_path.parent.mkdir(exist_ok=True)
    input_path.write_bytes(f"patient-{file_format}".encode())
    return AlignmentPlan(
        input_path=str(input_path),
        view_path=str(view_path),
        file_format=file_format,
        index_path=f"{view_path}.{'bai' if file_format == 'bam' else 'crai'}",
        reference_path=str(tmp_path / "reference genome.fa") if file_format == "cram" else None,
        reference_source="test",
        uncovered_contigs=(),
        unmapped_scan=unmapped_scan,
    )


def _run_conversion(tmp_path: Path, plan: AlignmentPlan, *, fast_mode: bool) -> tuple[list[str], MagicMock]:
    """Run conversion with shell and filesystem effects recorded.

    Args:
        tmp_path: Test directory receiving generated output names.
        plan: Proven alignment plan consumed by the stage.
        fast_mode: Whether the slice is final and therefore needs an index.

    Returns:
        Emitted commands and the mocked offset extractor.
    """
    commands: list[str] = []

    def record(command, log_file, critical=False, cwd=None):
        commands.append(command)
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        Path(log_file).write_text("")
        return True

    with (
        patch.object(fastq_bam_processing, "run_command", record),
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value=REGION),
        patch.object(fastq_bam_processing, "extract_unmapped_reads_from_offset") as extractor,
        patch.object(fastq_bam_processing.os, "replace"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            in_bam=plan.input_path,
            output=str(tmp_path / "run"),
            output_name="output",
            threads=4,
            config=CONFIG,
            fast_mode=fast_mode,
            plan=plan,
        )
    assert isinstance(extractor, MagicMock)
    return commands, extractor


@pytest.mark.parametrize(("fast_mode", "expects_index"), [(False, False), (True, True)])
def test_slice_uses_the_view_reference_threads_and_mode_specific_indexing(
    tmp_path: Path, fast_mode: bool, expects_index: bool
) -> None:
    plan = _plan(tmp_path, "cram")

    commands, _ = _run_conversion(tmp_path, plan, fast_mode=fast_mode)

    slice_command = commands[0]
    assert plan.view_path in slice_command
    assert plan.input_path not in slice_command
    assert f"-T '{plan.reference_path}'" in slice_command
    assert "-@ 4" in slice_command
    assert ("&& samtools index" in slice_command) is expects_index


def test_non_fast_bam_uses_the_plan_view_and_exact_plan_index_for_offsets(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")

    commands, extractor = _run_conversion(tmp_path, plan, fast_mode=False)

    assert plan.view_path in commands[0]
    assert all("input.bam.bai" not in command or "output_sliced.bam" in command for command in commands)
    extractor.assert_called_once_with(
        bam_file=plan.view_path,
        bai_file=plan.index_path,
        output_bam=str(tmp_path / "run" / "output_unmapped.bam"),
    )


@pytest.mark.parametrize(
    ("scan", "expected_fragment", "unexpected_fragment"),
    [
        ("indexed", "'*' -o", "set -o pipefail"),
        ("stream", "set -o pipefail", " '*' -o"),
    ],
)
def test_cram_unmapped_command_is_selected_by_the_proven_plan(
    tmp_path: Path, scan: str, expected_fragment: str, unexpected_fragment: str
) -> None:
    plan = _plan(tmp_path, "cram", unmapped_scan=scan)

    commands, _ = _run_conversion(tmp_path, plan, fast_mode=False)

    (unmapped_command,) = [command for command in commands if "-f 12" in command]
    assert expected_fragment in unmapped_command
    assert unexpected_fragment not in unmapped_command
    assert plan.view_path in unmapped_command
    assert f"-T '{plan.reference_path}'" in unmapped_command


def test_coverage_passes_the_proven_reference_to_samtools_depth(tmp_path: Path) -> None:
    depth_file = tmp_path / "cov_vntr_coverage.txt"
    seen: list[str] = []

    def run_depth(command, log_file, critical=False, cwd=None):
        seen.append(command)
        depth_file.write_text("chr1\t1\t10\nchr1\t2\t10\n")
        return True

    reference_path = tmp_path / "reference genome.fa"
    with patch.object(fastq_bam_processing, "run_command", run_depth):
        fastq_bam_processing.calculate_vntr_coverage(
            bam_file=str(tmp_path / "run" / "input.cram"),
            region=REGION,
            threads=4,
            config=CONFIG,
            output_dir=str(tmp_path),
            output_name="cov",
            reference_path=str(reference_path),
        )

    assert seen == [
        f"samtools depth -a -@ 4 --reference '{reference_path}' -r {REGION} {tmp_path}/run/input.cram > {depth_file}"
    ]


def _tree_digest(root: Path) -> dict[str, str]:
    """Hash every regular file below a test input directory.

    Args:
        root: Input tree to snapshot.

    Returns:
        Relative file names mapped to SHA-256 digests.
    """
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


@pytest.mark.parametrize("file_format", ["bam", "cram"])
def test_conversion_leaves_the_input_tree_byte_identical(tmp_path: Path, file_format: str) -> None:
    plan = _plan(tmp_path, file_format)
    input_root = Path(plan.input_path).parent
    Path(plan.index_path.replace(str(tmp_path / "run"), str(input_root))).write_bytes(b"patient-index")
    before = _tree_digest(input_root)

    _run_conversion(tmp_path, plan, fast_mode=False)

    assert _tree_digest(input_root) == before
