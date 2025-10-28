"""
Docker integration tests for VNtyper pipeline.

These tests use the SAME test orchestration and validation logic as local tests,
guaranteeing 100% test identity. The ONLY difference is the runner function.

Architecture:
- Shared test cases: tests.parametrization (get_bam_test_cases, etc.)
- Shared test logic: tests.test_orchestration (run_bam_test_case, etc.)
- Shared validation: tests.helpers (validate_kestrel_output, etc.)
- Docker-specific: Runner functions (this file) + fixtures (conftest.py)

This ensures local and Docker tests are identical except for execution environment.
"""

from pathlib import Path

import pytest

from tests.parametrization import get_advntr_test_cases, get_advntr_test_ids, get_bam_test_cases, get_bam_test_ids
from tests.test_orchestration import run_advntr_test_case, run_bam_test_case

from .conftest import run_vntyper_pipeline

# ============================================================================
# BAM Pipeline Tests
# ============================================================================


@pytest.mark.docker
@pytest.mark.parametrize("test_case", get_bam_test_cases(), ids=get_bam_test_ids())
def test_docker_bam_pipeline(test_case: dict, vntyper_container, tmp_path) -> None:
    """
    Test BAM pipeline in Docker container.

    This test uses the SAME orchestration logic as local tests.
    100% test identity guaranteed by shared modules.

    Args:
        test_case: Test configuration from test_data_config.json
        vntyper_container: Docker container fixture
        tmp_path: Pytest temp directory
    """

    # Define Docker-specific runner
    def docker_runner(bam_file: Path, reference: str, output_dir: Path) -> int:
        """Execute pipeline in Docker container."""
        return run_vntyper_pipeline(
            vntyper_container,
            bam_file,
            reference,
            output_dir,
            extra_modules=None,
        )

    # Create output directory
    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run test using SHARED orchestration logic
    # This is the SAME function used by local tests!
    run_bam_test_case(test_case, docker_runner, output_dir)


# ============================================================================
# adVNTR Module Tests
# ============================================================================


@pytest.mark.docker
@pytest.mark.slow
@pytest.mark.parametrize("test_case", get_advntr_test_cases(), ids=get_advntr_test_ids())
def test_docker_advntr_pipeline(test_case: dict, vntyper_container, tmp_path) -> None:
    """
    Test adVNTR module pipeline in Docker container.

    This test uses the SAME orchestration logic as local tests.
    100% test identity guaranteed by shared modules.

    Args:
        test_case: Test configuration from test_data_config.json
        vntyper_container: Docker container fixture
        tmp_path: Pytest temp directory
    """

    # Define Docker-specific runner
    def docker_runner(bam_file: Path, reference: str, output_dir: Path, extra_modules: list[str]) -> int:
        """Execute pipeline with adVNTR in Docker container."""
        return run_vntyper_pipeline(
            vntyper_container,
            bam_file,
            reference,
            output_dir,
            extra_modules=extra_modules,
        )

    # Create output directory
    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run test using SHARED orchestration logic
    run_advntr_test_case(test_case, docker_runner, output_dir)


# ============================================================================
# Container Health Tests
# ============================================================================


@pytest.mark.docker
def test_docker_container_health(vntyper_container) -> None:
    """
    Test Docker container is healthy and vntyper is accessible.

    Args:
        vntyper_container: Docker container fixture
    """
    # Test vntyper command exists
    result = vntyper_container.exec(["which", "vntyper"])
    assert result.exit_code == 0, "vntyper command not found in container"

    # Test vntyper version
    result = vntyper_container.exec(["vntyper", "--version"])
    assert result.exit_code == 0, "vntyper --version failed"

    # Test vntyper help
    result = vntyper_container.exec(["vntyper", "--help"])
    assert result.exit_code == 0, "vntyper --help failed"


@pytest.mark.docker
def test_docker_volume_mounts(vntyper_container) -> None:
    """
    Test Docker container has correct volume mounts.

    Args:
        vntyper_container: Docker container fixture
    """
    # Test /data mount (test data)
    result = vntyper_container.exec(["ls", "/data"])
    assert result.exit_code == 0, "/data mount not accessible"

    # Test /output mount (output directory)
    result = vntyper_container.exec(["ls", "/output"])
    assert result.exit_code == 0, "/output mount not accessible"

    # Test write permission in /output
    result = vntyper_container.exec(["touch", "/output/test_write.txt"])
    assert result.exit_code == 0, "/output mount not writable"


@pytest.mark.docker
def test_docker_dependencies(vntyper_container) -> None:
    """
    Test Docker container has required dependencies.

    Args:
        vntyper_container: Docker container fixture
    """
    # Test Java (required for Kestrel)
    result = vntyper_container.exec(["java", "-version"])
    assert result.exit_code == 0, "Java not found in container"

    # Test Python
    result = vntyper_container.exec(["python", "--version"])
    assert result.exit_code == 0, "Python not found in container"

    # Test samtools
    result = vntyper_container.exec(["samtools", "--version"])
    assert result.exit_code == 0, "samtools not found in container"

    # Test BWA
    result = vntyper_container.exec(["bwa"])
    # BWA returns non-zero when no args, but should be found
    # Just check it doesn't return 127 (command not found)
    assert result.exit_code != 127, "BWA not found in container"
