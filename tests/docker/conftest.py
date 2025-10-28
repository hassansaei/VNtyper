"""
Docker test fixtures using testcontainers.

This module provides pytest fixtures for Docker-based integration testing.
Uses testcontainers for automatic container lifecycle management.

Fixtures:
- vntyper_image: Session-scoped Docker image build
- vntyper_container: Module-scoped container with volume mounts
"""

from collections.abc import Generator
from pathlib import Path

import pytest
from testcontainers.core.container import DockerContainer


@pytest.fixture(scope="session")
def vntyper_image() -> Generator[str, None, None]:
    """
    Build VNtyper Docker image once per test session.

    Yields:
        str: Image tag

    Notes:
        - Built from project root Dockerfile
        - Cached for entire test session
        - Automatically cleaned up after session
    """
    import subprocess

    # Build image
    image_tag = "vntyper:test"
    project_root = Path(__file__).parent.parent.parent

    result = subprocess.run(
        ["docker", "build", "-t", image_tag, "."],
        cwd=str(project_root),
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(f"Docker build failed: {result.stderr}")

    yield image_tag

    # Cleanup (optional - Docker will handle)
    # subprocess.run(["docker", "rmi", image_tag], check=False)


@pytest.fixture(scope="module")
def vntyper_container(
    vntyper_image: str, tmp_path_factory: pytest.TempPathFactory
) -> Generator[DockerContainer, None, None]:
    """
    Create VNtyper container with volume mounts.

    Args:
        vntyper_image: Docker image tag from vntyper_image fixture
        tmp_path_factory: Pytest temp directory factory

    Yields:
        DockerContainer: Running container with mounts configured

    Notes:
        - Module-scoped: One container per test module
        - Mounts test data and output directories
        - Automatically cleaned up after module
    """
    # Get project paths
    project_root = Path(__file__).parent.parent.parent
    test_data_dir = project_root / "tests" / "test_data"

    # Create temp output directory for this module
    output_dir = tmp_path_factory.mktemp("docker_output")

    # Create container with volume mounts
    container = DockerContainer(vntyper_image)

    # Mount test data (read-only)
    container.with_volume_mapping(
        str(test_data_dir.absolute()),
        "/data",
        mode="ro",
    )

    # Mount output directory (read-write)
    container.with_volume_mapping(
        str(output_dir.absolute()),
        "/output",
        mode="rw",
    )

    # Keep container running
    container.with_command("tail -f /dev/null")

    # Start container
    container.start()

    # Verify container is running
    if not container.get_wrapped_container():
        raise RuntimeError("Container failed to start")

    yield container

    # Cleanup
    container.stop()


def run_vntyper_pipeline(
    container: DockerContainer,
    bam_file: Path,
    reference: str,
    output_dir: Path,
    extra_modules: list[str] | None = None,
) -> int:
    """
    Execute VNtyper pipeline inside Docker container.

    Args:
        container: Running Docker container
        bam_file: BAM file path (will be mapped to /data)
        reference: Reference assembly (hg19, hg38)
        output_dir: Output directory path (will be mapped to /output)
        extra_modules: Optional list of extra modules (e.g., ["advntr"])

    Returns:
        int: Exit code from pipeline execution

    Example:
        >>> exit_code = run_vntyper_pipeline(
        ...     container,
        ...     Path("tests/test_data/test.bam"),
        ...     "hg19",
        ...     Path("/tmp/output"),
        ...     extra_modules=["advntr"],
        ... )
        >>> assert exit_code == 0
    """
    # Build command
    bam_name = bam_file.name
    output_name = output_dir.name

    cmd = [
        "vntyper",
        "pipeline",
        "--bam",
        f"/data/{bam_name}",
        "--reference",
        reference,
        "--output",
        f"/output/{output_name}",
        "--threads",
        "2",
    ]

    if extra_modules:
        cmd.extend(["--extra-modules", ",".join(extra_modules)])

    # Execute command
    result = container.exec(cmd)

    return result.exit_code


def run_vntyper_fastq_pipeline(
    container: DockerContainer,
    fastq1: Path,
    fastq2: Path | None,
    reference: str,
    output_dir: Path,
    extra_modules: list[str] | None = None,
) -> int:
    """
    Execute VNtyper FASTQ pipeline inside Docker container.

    Args:
        container: Running Docker container
        fastq1: FASTQ file 1 path (will be mapped to /data)
        fastq2: Optional FASTQ file 2 path
        reference: Reference assembly
        output_dir: Output directory path
        extra_modules: Optional list of extra modules

    Returns:
        int: Exit code from pipeline execution
    """
    # Build command
    fastq1_name = fastq1.name
    output_name = output_dir.name

    cmd = [
        "vntyper",
        "pipeline",
        "--fastq",
        f"/data/{fastq1_name}",
        "--reference",
        reference,
        "--output",
        f"/output/{output_name}",
        "--threads",
        "2",
    ]

    if fastq2:
        cmd.extend(["--fastq2", f"/data/{fastq2.name}"])

    if extra_modules:
        cmd.extend(["--extra-modules", ",".join(extra_modules)])

    # Execute command
    result = container.exec(cmd)

    return result.exit_code
