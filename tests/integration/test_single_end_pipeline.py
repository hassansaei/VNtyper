"""A derived single-end BAM is routed to Kestrel without losing its reads."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pytest

from scripts.single_end_fixture import parse_single_end_fixture
from tests.parametrization import load_test_config
from tests.support.orchestration import PipelineRequest, PipelineRunResult, build_pipeline_argv, run_bam_test_case

logger = logging.getLogger(__name__)

pytestmark = pytest.mark.integration

_CONFIG = load_test_config()
_CASES = _CONFIG.get("integration_tests", {}).get("single_end_bam_tests", [])


@pytest.mark.parametrize("single_end_case", _CASES, ids=lambda case: case["test_name"])
def test_single_end_bam_produces_a_genotype(
    tmp_path: Path,
    ensure_test_data: None,
    single_end_case: dict,
) -> None:
    """Exercise A-161-1 with the single-end fixture declared in test data config."""
    declarations = {
        entry["name"]: parse_single_end_fixture(entry)
        for entry in _CONFIG.get("derived_fixtures", [])
        if entry.get("kind") == "single_end_bam"
    }
    spec = declarations[single_end_case["fixture_name"]]
    assert spec.output_bam.is_file(), f"Missing {spec.output_bam}; run `make cram-fixtures` to derive it."

    output_dir = tmp_path / single_end_case["test_name"]
    output_dir.mkdir()
    case = {**single_end_case, "bam": str(spec.output_bam)}

    def local_runner(request: PipelineRequest) -> PipelineRunResult:
        command = build_pipeline_argv(request, str)
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        logger.info("single-end pipeline stdout:\n%s", result.stdout)
        logger.info("single-end pipeline stderr:\n%s", result.stderr)
        return PipelineRunResult(result.returncode, result.stdout, result.stderr)

    run_bam_test_case(case, local_runner, output_dir)
