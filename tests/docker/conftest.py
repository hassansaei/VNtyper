"""
Docker test fixtures using testcontainers.

This module provides pytest fixtures for Docker-based integration testing.
Uses testcontainers for automatic container lifecycle management.

Fixtures:
- vntyper_image: Session-scoped Docker image build
- vntyper_container: Module-scoped container with volume mounts
"""

import subprocess
from collections.abc import Generator
from pathlib import Path
from typing import Optional

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
        ["docker", "build", "-f", "docker/Dockerfile", "-t", image_tag, "."],
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

    # Set permissions on output directory for container user (appuser, UID 1001)
    # This is necessary because the container runs as non-root user
    subprocess.run(["chmod", "777", str(output_dir)], check=True)

    # Create container with volume mounts
    container = DockerContainer(vntyper_image)

    # Mount test data (read-only) - using correct path from Docker README
    container.with_volume_mapping(
        str(test_data_dir.absolute()),
        "/opt/vntyper/input",
        mode="ro",
    )

    # Mount output directory (read-write) - using correct path from Docker README
    container.with_volume_mapping(
        str(output_dir.absolute()),
        "/opt/vntyper/output",
        mode="rw",
    )

    # Override entrypoint and keep container running for testing
    # The production entrypoint validates commands, but for testing we need
    # direct shell access to run arbitrary commands via exec()
    # We bypass the entrypoint but maintain the correct user (appuser, UID 1001)
    # to preserve file permissions on mounted volumes
    # Use bash (not sh) for `source` command support
    container.with_kwargs(
        entrypoint=["/bin/bash", "-c", "tail -f /dev/null"],
        user="1001:1001",  # Match the appuser UID:GID from Dockerfile
        working_dir="/opt/vntyper",  # Set working directory as per Docker README
    )

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
    extra_modules: Optional[list[str]] = None,
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
    # Since we bypassed the entrypoint for testing, we need to manually
    # activate the conda environment and run vntyper
    bam_name = bam_file.name
    output_name = output_dir.name

    # Build vntyper command arguments
    # Use correct paths as per Docker README documentation
    vntyper_args = [
        "pipeline",
        "--bam",
        f"/opt/vntyper/input/{bam_name}",
        "--reference",
        reference,
        "--output",
        f"/opt/vntyper/output/{output_name}",
        "--threads",
        "2",
    ]

    if extra_modules:
        vntyper_args.extend(["--extra-modules", ",".join(extra_modules)])

    # Execute via conda run since we bypassed the entrypoint
    # Use --no-capture-output to stream stdout/stderr properly
    cmd = [
        "/bin/bash",
        "-c",
        f"source /opt/conda/etc/profile.d/conda.sh && conda run --no-capture-output -n vntyper vntyper {' '.join(vntyper_args)}",
    ]

    # Execute command
    result = container.exec(cmd)

    return result.exit_code


def run_vntyper_fastq_pipeline(
    container: DockerContainer,
    fastq1: Path,
    fastq2: Optional[Path],
    reference: str,
    output_dir: Path,
    extra_modules: Optional[list[str]] = None,
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
    # Since we bypassed the entrypoint for testing, we need to manually
    # activate the conda environment and run vntyper
    fastq1_name = fastq1.name
    output_name = output_dir.name

    # Build vntyper command arguments
    # Use correct paths as per Docker README documentation
    vntyper_args = [
        "pipeline",
        "--fastq",
        f"/opt/vntyper/input/{fastq1_name}",
        "--reference",
        reference,
        "--output",
        f"/opt/vntyper/output/{output_name}",
        "--threads",
        "2",
    ]

    if fastq2:
        vntyper_args.extend(["--fastq2", f"/opt/vntyper/input/{fastq2.name}"])

    if extra_modules:
        vntyper_args.extend(["--extra-modules", ",".join(extra_modules)])

    # Execute via conda run since we bypassed the entrypoint
    # Use --no-capture-output to stream stdout/stderr properly
    cmd = [
        "/bin/bash",
        "-c",
        f"source /opt/conda/etc/profile.d/conda.sh && conda run --no-capture-output -n vntyper vntyper {' '.join(vntyper_args)}",
    ]

    # Execute command
    result = container.exec(cmd)

    return result.exit_code
