"""
End-to-end integration tests for the VNtyper pipeline.

Every scenario is read from ``tests/test_data_config.json`` and validated through the
shared orchestration in ``tests/support/orchestration.py``, so this tier and the Docker
tier assert exactly the same things about exactly the same samples. The only difference
between the two is the runner: here the CLI is a subprocess on this machine, there it is
``container.exec``.

What this module deliberately does **not** define:

* ``test_config`` and ``ensure_test_data`` - both live in ``tests/conftest.py``. The
  copies that used to sit here ignored ``VNTYPER_TEST_DATA_SKIP_DOWNLOAD``, so the CI
  fast-fail path was silently bypassed for this tier.
* ``compute_md5`` and ``download_file`` - both live in ``tests/support/data_utils.py``.
* the test-case lists - ``tests/parametrization.py`` is the single source, shared with
  ``tests/docker/test_docker_pipeline.py``.

These tests perform real I/O (subprocesses, ~1.2 GB of sample data, reference genomes)
and are therefore **not** part of the unit tier. They carry ``integration`` only.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.parametrization import (
    get_advntr_test_cases,
    get_advntr_test_ids,
    get_bam_test_cases,
    get_bam_test_ids,
    get_fastq_test_cases,
    get_fastq_test_ids,
)
from tests.support.orchestration import (
    ADVNTR_TIMEOUT_SECONDS,
    PipelineRequest,
    PipelineRunResult,
    build_pipeline_argv,
    run_advntr_test_case,
    run_bam_test_case,
    run_fastq_test_case,
)

logger = logging.getLogger(__name__)


def _fresh_output_dir(tmp_path: Path, test_name: str) -> Path:
    """Return an empty, per-case output directory under ``tmp_path``.

    Args:
        tmp_path: The test's temporary directory.
        test_name: The case's ``test_name``, used as the directory name.

    Returns:
        Path: The created directory.
    """
    output_dir = tmp_path / test_name
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _run_cli(command: list[str]) -> subprocess.CompletedProcess[str]:
    """Run the ``vntyper`` CLI and return its complete captured result.

    stdout and stderr are logged rather than swallowed: the pipeline writes its whole
    stage log to stderr, and it is the only diagnosis available when a case fails.

    Args:
        command: The full argv to execute.

    Returns:
        Captured process result, including the diagnostic stderr contract.
    """
    logger.info("Command to execute: %s", " ".join(command))
    result = subprocess.run(command, capture_output=True, text=True, check=False)

    logger.info("Return code: %d", result.returncode)
    logger.info("STDOUT:\n%s", result.stdout)
    if result.returncode != 0:
        logger.error("STDERR:\n%s", result.stderr)
    else:
        logger.debug("STDERR:\n%s", result.stderr)

    return result


def _run_local_pipeline(request: PipelineRequest) -> PipelineRunResult:
    """Execute a canonical request locally, mapping paths without changing identity.

    Args:
        request: Transport-independent pipeline request.

    Returns:
        Complete captured local process result.
    """
    pipeline_argv = build_pipeline_argv(request, str)
    completed = _run_cli([sys.executable, "-m", "vntyper.cli", *pipeline_argv[1:]])
    return PipelineRunResult(completed.returncode, completed.stdout, completed.stderr)


#
# 1) FASTQ Tests
#
@pytest.mark.integration
@pytest.mark.parametrize("fastq_case", get_fastq_test_cases(), ids=get_fastq_test_ids())
def test_fastq_input(tmp_path: Path, ensure_test_data: None, fastq_case: dict) -> None:
    """Run the pipeline from paired FASTQ input and validate the declared outputs.

    Args:
        tmp_path: Pytest temporary directory.
        ensure_test_data: Session fixture that verifies (and, outside CI, fetches) the data.
        fastq_case: One entry of ``integration_tests.fastq_tests``.
    """
    output_dir = _fresh_output_dir(tmp_path, fastq_case["test_name"])

    run_fastq_test_case(fastq_case, _run_local_pipeline, output_dir)


#
# 2) BAM Tests
#
@pytest.mark.integration
@pytest.mark.parametrize("bam_case", get_bam_test_cases(), ids=get_bam_test_ids())
def test_bam_input_with_kestrel_checks(tmp_path: Path, ensure_test_data: None, bam_case: dict) -> None:
    """Run the pipeline from BAM input and validate Kestrel, coverage and artefacts.

    Validation is delegated wholesale to ``run_bam_test_case``. That matters for the
    expected-negative cases: this test used to skip the confidence assertion whenever the
    expectation was a negative call, which is precisely the class most likely to move
    under a filtering change. The shared validator asserts every expectation, positive or
    negative, with no exemption.

    Args:
        tmp_path: Pytest temporary directory.
        ensure_test_data: Session fixture that verifies (and, outside CI, fetches) the data.
        bam_case: One entry of ``integration_tests.bam_tests``.
    """
    output_dir = _fresh_output_dir(tmp_path, bam_case["test_name"])
    run_bam_test_case(bam_case, _run_local_pipeline, output_dir)


#
# 3) adVNTR Tests
#
@pytest.mark.integration
@pytest.mark.timeout(ADVNTR_TIMEOUT_SECONDS)
@pytest.mark.parametrize("advntr_case", get_advntr_test_cases(), ids=get_advntr_test_ids())
def test_advntr_input(tmp_path: Path, ensure_test_data: None, advntr_case: dict) -> None:
    """Run the pipeline with ``--extra-modules advntr`` and validate the adVNTR call.

    Args:
        tmp_path: Pytest temporary directory.
        ensure_test_data: Session fixture that verifies (and, outside CI, fetches) the data.
        advntr_case: One entry of ``integration_tests.advntr_tests``.
    """
    output_dir = _fresh_output_dir(tmp_path, advntr_case["test_name"])

    run_advntr_test_case(advntr_case, _run_local_pipeline, output_dir)
