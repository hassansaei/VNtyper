"""
Docker test fixtures using testcontainers.

This module provides pytest fixtures for Docker-based integration testing.
Uses testcontainers for automatic container lifecycle management.

Fixtures:
- vntyper_image: Session-scoped Docker image build
- vntyper_container: Module-scoped container with volume mounts
"""

import os
import shlex
import subprocess
from collections.abc import Generator
from pathlib import Path
from typing import Any, Protocol

import pytest
from testcontainers.core.container import DockerContainer

from tests.support.orchestration import PipelineRequest, PipelineRunResult, build_pipeline_argv


def _exit_code(result: Any) -> int:
    """Return an exec result's exit code, treating an unknown code as failure.

    testcontainers types `ExecResult.exit_code` as `int | None`; None means the daemon
    did not report one (e.g. the exec was interrupted). Callers compare against 0, so
    mapping None to a non-zero code keeps "unknown" from reading as success.

    Args:
        result: An exec result carrying an `exit_code` attribute.

    Returns:
        int: The exit code, or 1 when it is unknown.
    """
    code = getattr(result, "exit_code", None)
    return 1 if code is None else int(code)


@pytest.fixture(scope="session")
def vntyper_image() -> Generator[str, None, None]:
    """
    Provide a VNtyper Docker image for the session, building one only if needed.

    Yields:
        str: Image tag

    Raises:
        RuntimeError: If no image was supplied and the local build fails.

    Notes:
        - Set VNTYPER_TEST_IMAGE to an existing tag to skip building entirely. CI does
          this with the image the build job already produced and pushed; without it the
          fixture rebuilt the whole image from scratch (measured: 1042s of fixture setup
          to run 17s of assertions), which is why the Docker test jobs took ~20 minutes.
        - The local build path is retained so `make test-docker` still works on a
          developer machine with no image present.
    """
    preexisting = os.environ.get("VNTYPER_TEST_IMAGE")
    if preexisting:
        yield preexisting
        return

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
) -> Generator[tuple[DockerContainer, Path], None, None]:
    """
    Create VNtyper container with volume mounts.

    Args:
        vntyper_image: Docker image tag from vntyper_image fixture
        tmp_path_factory: Pytest temp directory factory

    Yields:
        tuple: (DockerContainer, Path) - Running container and mounted output directory

    Notes:
        - Module-scoped: One container per test module
        - Mounts test data and output directories
        - Automatically cleaned up after module
    """
    # Get project paths
    project_root = Path(__file__).parent.parent.parent
    test_data_dir = project_root / "tests" / "data"  # BAM files are in tests/data/, not tests/test_data/

    # Create temp output directory for this module
    output_dir = tmp_path_factory.mktemp("docker_output")

    # Set permissions on output directory for container user (appuser, UID 1001)
    # This is necessary because the container runs as non-root user
    # Using chmod 777 is acceptable here because:
    # 1. This is an isolated pytest tmpdir that gets cleaned up after tests
    # 2. The container user (UID 1001) needs write access
    # 3. This approach is simpler than chown/chgrp which requires sudo
    # Note: For production environments, use chmod 775 with group ownership instead
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

    yield container, output_dir

    # Cleanup
    container.stop()


class ContainerExecutor(Protocol):
    """Structural subset of ``DockerContainer`` used by the pipeline runner."""

    def exec(self, command: list[str]) -> Any:
        """Execute one command in the running container."""


def _map_registered_input(path: Path, test_data_root: Path) -> str:
    lexical_root = test_data_root.absolute()
    lexical_path = path.absolute()
    try:
        relative = lexical_path.relative_to(lexical_root)
    except ValueError as exc:
        raise ValueError(f"Path is not registered test data: {path}") from exc
    current = lexical_root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"Path is not registered test data because it contains a symlink: {path}")

    root = test_data_root.resolve(strict=True)
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ValueError(f"Path is not registered test data: {path}") from exc
    if not resolved.is_file() or not resolved.is_relative_to(root):
        raise ValueError(f"Path is not registered test data: {path}")
    return f"/opt/vntyper/input/{resolved.relative_to(root).as_posix()}"


def _map_case_output(path: Path, output_mount_root: Path) -> str:
    root = output_mount_root.resolve(strict=True)
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ValueError(f"Docker case output directory does not exist: {path}") from exc
    if resolved == root:
        raise ValueError("Parameterized Docker cases may not use the output mount root")
    if not resolved.is_dir() or not resolved.is_relative_to(root):
        raise ValueError(f"Docker case output directory escapes its mount root: {path}")
    if any(resolved.iterdir()):
        raise ValueError(f"Docker case output directory must be initially empty: {path}")
    return f"/opt/vntyper/output/{resolved.relative_to(root).as_posix()}"


def map_docker_request_path(
    path: Path,
    *,
    request: PipelineRequest,
    test_data_root: Path,
    output_mount_root: Path,
) -> str:
    """Map one canonical request path into a registered read-only or isolated mount.

    Args:
        path: Host or image-relative path from ``request``.
        request: Canonical request owning the path roles.
        test_data_root: Host directory mounted read-only at ``/opt/vntyper/input``.
        output_mount_root: Host directory mounted read-write at ``/opt/vntyper/output``.

    Returns:
        Container path preserving the resource's relative identity.

    Raises:
        ValueError: If a path is unregistered, escapes a mount, or targets shared output.
    """
    if path == request.output_dir:
        return _map_case_output(path, output_mount_root)
    if path == request.reference_fasta and not path.is_absolute():
        if not path.parts or path.parts[0] != "reference" or ".." in path.parts:
            raise ValueError(f"Image reference path is not registered: {path}")
        return f"/opt/vntyper/{path.as_posix()}"
    return _map_registered_input(path, test_data_root)


def _exec_output(result: Any) -> str:
    output = getattr(result, "output", b"")
    if isinstance(output, bytes):
        return output.decode(errors="replace")
    return str(output or "")


def run_vntyper_pipeline(
    container: ContainerExecutor,
    request: PipelineRequest,
    *,
    test_data_root: Path,
    output_mount_root: Path,
) -> PipelineRunResult:
    """Execute one canonical pipeline request through Docker's path transport.

    Args:
        container: Running container executor.
        request: Transport-independent pipeline request.
        test_data_root: Host input mount root.
        output_mount_root: Host output mount root.

    Returns:
        Captured Docker exec exit and combined output.
    """

    def mapper(path: Path) -> str:
        return map_docker_request_path(
            path,
            request=request,
            test_data_root=test_data_root,
            output_mount_root=output_mount_root,
        )

    argv = build_pipeline_argv(request, mapper)
    command = [
        "/bin/bash",
        "-c",
        "source /opt/conda/etc/profile.d/conda.sh && conda run --no-capture-output -n vntyper " + shlex.join(argv),
    ]
    result = container.exec(command)
    return PipelineRunResult(exit_code=_exit_code(result), stdout=_exec_output(result), stderr="")
