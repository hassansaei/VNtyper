"""Wiring tests for consumers of the proven alignment preflight plan."""

from __future__ import annotations

import hashlib
import inspect
import os
from dataclasses import replace
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vntyper.scripts import fastq_bam_processing
from vntyper.scripts.alignment_contract import AlignmentPlan, index_candidate_names

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


def _run_conversion(
    tmp_path: Path,
    plan: AlignmentPlan,
    *,
    fast_mode: bool,
    keep_intermediates: bool = False,
) -> tuple[list[str], MagicMock]:
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
            output=str(tmp_path / "run"),
            output_name="output",
            threads=4,
            config=CONFIG,
            fast_mode=fast_mode,
            keep_intermediates=keep_intermediates,
            plan=plan,
        )
    assert isinstance(extractor, MagicMock)
    return commands, extractor


def test_conversion_api_has_one_authoritative_alignment_input() -> None:
    assert "in_bam" not in inspect.signature(fastq_bam_processing.process_bam_to_fastq).parameters


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


def test_keep_intermediates_regenerates_stale_artifacts_for_the_new_plan(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")
    run_dir = tmp_path / "run"
    stale_paths = [
        run_dir / "output_sliced.bam",
        run_dir / "output_R1.fastq.gz",
        run_dir / "output_R2.fastq.gz",
        run_dir / "output_other.fastq.gz",
        run_dir / "output_single.fastq.gz",
    ]
    for stale_path in stale_paths:
        stale_path.write_bytes(b"prior-patient-artifact")

    commands, _ = _run_conversion(tmp_path, plan, fast_mode=True, keep_intermediates=True)

    assert len(commands) == 2
    assert plan.view_path in commands[0]
    assert "samtools view" in commands[0]
    assert "samtools fastq" in commands[1]
    assert str(run_dir / "output_sliced.bam") in commands[1]


def _assert_rejected_before_conversion_work(
    tmp_path: Path,
    plan: AlignmentPlan,
    *,
    fast_mode: bool,
    error_match: str = "derived alignment-conversion destination",
) -> None:
    """Assert an unsafe derived entry stops both region and command execution."""
    region_resolver = MagicMock(return_value=REGION)
    command_runner = MagicMock(return_value=True)
    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", region_resolver),
        patch.object(fastq_bam_processing, "run_command", command_runner),
        patch.object(fastq_bam_processing, "extract_unmapped_reads_from_offset"),
        patch.object(fastq_bam_processing.os, "replace"),
        pytest.raises(ValueError, match=error_match),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=4,
            config=CONFIG,
            fast_mode=fast_mode,
            plan=plan,
        )

    region_resolver.assert_not_called()
    command_runner.assert_not_called()


def test_stale_slice_symlink_to_patient_alignment_is_rejected_before_work(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")
    input_root = Path(plan.input_path).parent
    stale_slice = tmp_path / "run" / "output_sliced.bam"
    stale_slice.symlink_to(plan.input_path)
    before = _tree_digest(input_root)

    _assert_rejected_before_conversion_work(tmp_path, plan, fast_mode=True)

    assert _tree_digest(input_root) == before
    assert stale_slice.is_symlink()


def test_stale_r1_hardlink_to_patient_index_is_rejected_before_work(tmp_path: Path) -> None:
    initial_plan = _plan(tmp_path, "bam")
    source_index = Path(initial_plan.input_path).with_suffix(".bam.bai")
    source_index.write_bytes(b"patient-index")
    plan = replace(initial_plan, index_path=str(source_index))
    stale_r1 = tmp_path / "run" / "output_R1.fastq.gz"
    stale_r1.hardlink_to(source_index)
    before = _tree_digest(source_index.parent)

    _assert_rejected_before_conversion_work(tmp_path, plan, fast_mode=True)

    assert _tree_digest(source_index.parent) == before
    assert stale_r1.stat().st_ino == source_index.stat().st_ino


@pytest.mark.parametrize(
    ("stale_name", "file_format", "fast_mode"),
    [
        ("output_sliced.bam.bai", "bam", True),
        ("output_sliced.bam.csi", "bam", True),
        ("output_sliced.bai", "bam", True),
        ("output_sliced.csi", "bam", True),
        ("output_unmapped.bam", "bam", False),
        ("output_sliced_unmapped.bam", "bam", False),
        ("output_R2.fastq.gz", "bam", True),
        ("output_other.fastq.gz", "bam", True),
        ("output_single.fastq.gz", "bam", True),
        ("output_slice.log", "bam", True),
        ("output_filter.log", "cram", False),
        ("output_merge.log", "bam", False),
        ("output_index.log", "bam", False),
        ("output_sort_fastq.log", "bam", True),
    ],
)
def test_stale_index_and_log_symlinks_are_rejected_before_work(
    tmp_path: Path,
    stale_name: str,
    file_format: str,
    fast_mode: bool,
) -> None:
    plan = _plan(tmp_path, file_format)
    protected = Path(plan.input_path)
    stale_path = tmp_path / "run" / stale_name
    stale_path.symlink_to(protected)
    before = _tree_digest(protected.parent)

    _assert_rejected_before_conversion_work(tmp_path, plan, fast_mode=fast_mode)

    assert _tree_digest(protected.parent) == before
    assert stale_path.is_symlink()


def test_stale_non_regular_derived_destination_is_rejected_before_work(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")
    stale_unmapped = tmp_path / "run" / "output_unmapped.bam"
    stale_unmapped.mkdir()

    _assert_rejected_before_conversion_work(tmp_path, plan, fast_mode=False)

    assert stale_unmapped.is_dir()


def test_derived_hardlink_to_the_plan_view_is_reported_as_a_protected_alias(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")
    view = Path(plan.view_path)
    view.write_bytes(b"proven-view")
    stale_r2 = tmp_path / "run" / "output_R2.fastq.gz"
    stale_r2.hardlink_to(view)

    _assert_rejected_before_conversion_work(
        tmp_path,
        plan,
        fast_mode=True,
        error_match="aliases protected input",
    )


def test_multiply_linked_derived_artifact_is_rejected_even_when_unrelated_to_the_plan(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")
    unrelated = tmp_path / "unrelated-artifact"
    unrelated.write_bytes(b"unrelated")
    stale_other = tmp_path / "run" / "output_other.fastq.gz"
    stale_other.hardlink_to(unrelated)

    _assert_rejected_before_conversion_work(
        tmp_path,
        plan,
        fast_mode=True,
        error_match="multiple hard links",
    )

    assert unrelated.read_bytes() == b"unrelated"


def test_single_link_regular_derived_artifacts_remain_replaceable(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")
    run_dir = tmp_path / "run"
    (run_dir / "output_sliced.bam").write_bytes(b"stale-slice")
    (run_dir / "output_R1.fastq.gz").write_bytes(b"stale-fastq")

    commands, _ = _run_conversion(tmp_path, plan, fast_mode=True)

    assert len(commands) == 2


@pytest.mark.parametrize(
    "stale_index_name",
    [
        "output_sliced.bam.bai",
        "output_sliced.bai",
        "output_sliced.bam.csi",
        "output_sliced.csi",
    ],
)
def test_successful_rerun_removes_every_stale_slice_index_before_regenerating_the_default_bai(
    tmp_path: Path,
    stale_index_name: str,
) -> None:
    plan = _plan(tmp_path, "bam")
    run_dir = tmp_path / "run"
    sliced_bam = run_dir / "output_sliced.bam"
    index_candidates = tuple(Path(path) for path in index_candidate_names(str(sliced_bam), "bam"))
    stale_index = run_dir / stale_index_name
    stale_index.write_bytes(b"prior-patient-index")
    stale_inode = stale_index.stat().st_ino
    commands: list[str] = []

    def regenerate_default_index(command, log_file, critical=False, cwd=None):
        commands.append(command)
        if len(commands) == 1:
            assert all(not os.path.lexists(candidate) for candidate in index_candidates)
            sliced_bam.write_bytes(b"fresh-slice")
            index_candidates[0].write_bytes(b"fresh-default-bai")
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value=REGION),
        patch.object(fastq_bam_processing, "run_command", regenerate_default_index),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=run_dir,
            output_name="output",
            threads=4,
            config=CONFIG,
            fast_mode=True,
            plan=plan,
        )

    assert len(commands) == 2
    assert index_candidates[0].read_bytes() == b"fresh-default-bai"
    if stale_index == index_candidates[0]:
        assert index_candidates[0].stat().st_ino != stale_inode
    assert all(not candidate.exists() for candidate in index_candidates[1:])


def test_unsafe_fastq_destination_leaves_validated_stale_csi_untouched(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "bam")
    run_dir = tmp_path / "run"
    stale_csi = run_dir / "output_sliced.bam.csi"
    stale_csi.write_bytes(b"prior-patient-csi")
    stale_inode = stale_csi.stat().st_ino
    unsafe_fastq = run_dir / "output_R1.fastq.gz"
    unsafe_fastq.symlink_to(plan.input_path)
    before = hashlib.sha256(stale_csi.read_bytes()).hexdigest()

    _assert_rejected_before_conversion_work(tmp_path, plan, fast_mode=True)

    assert stale_csi.stat().st_ino == stale_inode
    assert hashlib.sha256(stale_csi.read_bytes()).hexdigest() == before
    assert unsafe_fastq.is_symlink()
