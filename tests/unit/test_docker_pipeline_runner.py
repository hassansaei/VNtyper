"""Unit contracts for Docker's path-only pipeline transport."""

from pathlib import Path
from typing import Any

import pytest

from tests.docker import conftest as docker_support
from tests.support.orchestration import PipelineRequest, PipelineRunResult

pytestmark = pytest.mark.unit


def _request(data_root: Path, output_dir: Path, **changes: Any) -> PipelineRequest:
    values: dict[str, Any] = {
        "input_kind": "bam",
        "input_paths": (data_root / "nested" / "sample.bam",),
        "reference_assembly": "hg19",
        "output_dir": output_dir,
        "threads": 7,
        "log_level": "WARNING",
        "cli_options": ("--fast-mode", "--archive-results"),
        "reference_fasta": None,
    }
    values.update(changes)
    return PipelineRequest(**values)


def test_docker_mapper_preserves_nested_input_and_case_output_paths(tmp_path: Path) -> None:
    """Basename mapping or a global output path would collide distinct registered resources/cases."""
    data_root = tmp_path / "data"
    input_path = data_root / "nested" / "sample.bam"
    input_path.parent.mkdir(parents=True)
    input_path.write_bytes(b"bam")
    output_root = tmp_path / "docker_output0"
    output_dir = output_root / "suite" / "case"
    output_dir.mkdir(parents=True)
    request = _request(data_root, output_dir)

    assert (
        docker_support.map_docker_request_path(
            input_path,
            request=request,
            test_data_root=data_root,
            output_mount_root=output_root,
        )
        == "/opt/vntyper/input/nested/sample.bam"
    )
    assert (
        docker_support.map_docker_request_path(
            output_dir,
            request=request,
            test_data_root=data_root,
            output_mount_root=output_root,
        )
        == "/opt/vntyper/output/suite/case"
    )


def test_docker_mapper_rejects_unregistered_escape_symlink_and_nonempty_output(tmp_path: Path) -> None:
    """Only registered read-only inputs and a fresh isolated case directory may cross the mount boundary."""
    data_root = tmp_path / "data"
    data_root.mkdir()
    output_root = tmp_path / "docker_output0"
    output_root.mkdir()
    output_dir = output_root / "case"
    output_dir.mkdir()
    request = _request(data_root, output_dir)
    outside = tmp_path / "outside.bam"
    outside.write_bytes(b"bam")
    linked = data_root / "linked.bam"
    linked.symlink_to(outside)

    for path in (outside, data_root / "missing.bam", linked):
        with pytest.raises(ValueError, match="registered test data"):
            docker_support.map_docker_request_path(
                path,
                request=request,
                test_data_root=data_root,
                output_mount_root=output_root,
            )

    with pytest.raises(ValueError, match="mount root"):
        docker_support.map_docker_request_path(
            output_root,
            request=_request(data_root, output_root),
            test_data_root=data_root,
            output_mount_root=output_root,
        )
    (output_dir / "stale.txt").write_text("stale", encoding="utf-8")
    with pytest.raises(ValueError, match="initially empty"):
        docker_support.map_docker_request_path(
            output_dir,
            request=request,
            test_data_root=data_root,
            output_mount_root=output_root,
        )


@pytest.mark.parametrize(
    ("reference", "expected"),
    [
        ("mounted", "/opt/vntyper/input/references/chr1.fa"),
        ("image", "/opt/vntyper/reference/alignment/chr1.hg19.fa"),
    ],
)
def test_docker_mapper_supports_mounted_or_image_cram_reference(tmp_path: Path, reference: str, expected: str) -> None:
    """An explicit CRAM reference remains read-only whether fixture-mounted or shipped in the image."""
    data_root = tmp_path / "data"
    cram = data_root / "cram" / "sample.cram"
    mounted_reference = data_root / "references" / "chr1.fa"
    cram.parent.mkdir(parents=True)
    mounted_reference.parent.mkdir()
    cram.write_bytes(b"cram")
    mounted_reference.write_text(">chr1\nA\n", encoding="utf-8")
    output_root = tmp_path / "docker_output0"
    output_dir = output_root / "case"
    output_dir.mkdir(parents=True)
    reference_path = mounted_reference if reference == "mounted" else Path("reference/alignment/chr1.hg19.fa")
    request = _request(
        data_root,
        output_dir,
        input_kind="cram",
        input_paths=(cram,),
        reference_fasta=reference_path,
    )

    assert (
        docker_support.map_docker_request_path(
            reference_path,
            request=request,
            test_data_root=data_root,
            output_mount_root=output_root,
        )
        == expected
    )


class _ExecResult:
    exit_code = 1
    output = b"captured combined docker output"


class _Container:
    def __init__(self) -> None:
        self.commands: list[list[str]] = []

    def exec(self, command: list[str]) -> _ExecResult:
        self.commands.append(command)
        return _ExecResult()


def test_generic_docker_runner_builds_exact_request_argv_and_captures_output(tmp_path: Path) -> None:
    """Docker may map paths, but it may not omit options, append defaults, or discard exec output."""
    data_root = tmp_path / "data"
    input_path = data_root / "nested" / "sample.bam"
    input_path.parent.mkdir(parents=True)
    input_path.write_bytes(b"bam")
    output_root = tmp_path / "docker_output0"
    output_dir = output_root / "suite" / "case"
    output_dir.mkdir(parents=True)
    request = _request(data_root, output_dir)
    container = _Container()

    result = docker_support.run_vntyper_pipeline(
        container,
        request,
        test_data_root=data_root,
        output_mount_root=output_root,
    )

    assert result == PipelineRunResult(1, "captured combined docker output", "")
    assert container.commands == [
        [
            "/bin/bash",
            "-c",
            (
                "source /opt/conda/etc/profile.d/conda.sh && conda run --no-capture-output -n vntyper "
                "vntyper -l WARNING pipeline --bam /opt/vntyper/input/nested/sample.bam --threads 7 "
                "--reference-assembly hg19 --output-dir /opt/vntyper/output/suite/case --fast-mode "
                "--archive-results"
            ),
        ]
    ]


def test_generic_docker_runner_maps_fastq_output_to_its_case_suffix(tmp_path: Path) -> None:
    """The old FASTQ runner's basename/global-output mapping cannot recur."""
    data_root = tmp_path / "data"
    fastq1 = data_root / "nested" / "R1.fastq.gz"
    fastq2 = data_root / "nested" / "R2.fastq.gz"
    fastq1.parent.mkdir(parents=True)
    fastq1.write_bytes(b"r1")
    fastq2.write_bytes(b"r2")
    output_root = tmp_path / "docker_output0"
    output_dir = output_root / "fastq-case"
    output_dir.mkdir(parents=True)
    request = _request(
        data_root,
        output_dir,
        input_kind="fastq",
        input_paths=(fastq1, fastq2),
        cli_options=(),
    )
    container = _Container()

    docker_support.run_vntyper_pipeline(
        container,
        request,
        test_data_root=data_root,
        output_mount_root=output_root,
    )

    command = container.commands[0][2]
    assert "--fastq1 /opt/vntyper/input/nested/R1.fastq.gz" in command
    assert "--fastq2 /opt/vntyper/input/nested/R2.fastq.gz" in command
    assert "--output-dir /opt/vntyper/output/fastq-case" in command
    assert "--output-dir /opt/vntyper/output " not in command
