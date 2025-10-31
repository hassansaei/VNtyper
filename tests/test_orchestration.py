"""
Shared test orchestration logic for VNtyper.

This module contains the CORE test logic that must be identical
for both local and Docker tests.

The ONLY difference between local and Docker tests is the
"runner" function that executes the pipeline.

This architecture guarantees 100% test identity.
"""

from pathlib import Path
from typing import Callable, Optional

from tests.helpers import (
    assert_required_files,
    validate_advntr_output,
    validate_coverage_output,
    validate_kestrel_output,
)


def run_bam_test_case(
    test_case: dict,
    runner: Callable[[Path, str, Path], int],
    output_dir: Path,
) -> None:
    """
    Universal BAM test case orchestration.

    This function contains the COMPLETE test logic for BAM tests.
    It's used by BOTH local and Docker tests to guarantee 100% identity.

    Args:
        test_case: Test configuration dict from test_data_config.json
        runner: Function that executes the pipeline
                Signature: runner(bam_file, reference, output_dir) -> exit_code
        output_dir: Output directory path

    Raises:
        AssertionError: If any validation fails

    Example (Local):
        def local_runner(bam, ref, out):
            result = subprocess.run(["vntyper", "pipeline", "--bam", str(bam), ...])
            return result.returncode

        run_bam_test_case(test_case, local_runner, output_dir)

    Example (Docker):
        def docker_runner(bam, ref, out):
            return container.exec(["vntyper", "pipeline", "--bam", f"/input/{bam.name}", ...]).exit_code

        run_bam_test_case(test_case, docker_runner, output_dir)
    """
    # 1. Extract test configuration
    bam_file = Path(test_case["bam"])
    reference = test_case["reference_assembly"]

    # 2. Run pipeline via runner (ONLY difference between local and Docker)
    exit_code = runner(bam_file, reference, output_dir)

    # 3. Assert success
    assert exit_code == 0, f"Pipeline failed with exit code {exit_code}"

    # 4. Validate required files exist
    required_files = [
        "summary_report.html",
        "kestrel/kestrel_result.tsv",
        "coverage/coverage_summary.tsv",
        "pipeline_summary.json",
    ]

    if test_case.get("check_igv_report"):
        required_files.append("igv_report.html")

    assert_required_files(output_dir, required_files)

    # 5. Validate Kestrel output
    validate_kestrel_output(output_dir, test_case["kestrel_assertions"])

    # 6. Validate coverage output
    coverage_metrics = validate_coverage_output(output_dir)
    assert coverage_metrics["mean_cov"] >= 0, "Coverage mean should be non-negative"


def run_advntr_test_case(
    test_case: dict,
    runner: Callable[[Path, str, Path, list[str]], int],
    output_dir: Path,
) -> None:
    """
    Universal adVNTR test case orchestration.

    This function contains the COMPLETE test logic for adVNTR tests.
    Identical for local and Docker tests.

    Args:
        test_case: Test configuration dict
        runner: Function that executes the pipeline with extra modules
                Signature: runner(bam_file, reference, output_dir, extra_modules) -> exit_code
        output_dir: Output directory path

    Raises:
        AssertionError: If any validation fails
    """
    # 1. Extract test configuration
    bam_file = Path(test_case["bam"])
    reference = test_case["reference_assembly"]
    extra_modules = ["advntr"]

    # 2. Run pipeline via runner
    exit_code = runner(bam_file, reference, output_dir, extra_modules)

    # 3. Assert success
    assert exit_code == 0, f"adVNTR pipeline failed with exit code {exit_code}"

    # 4. Validate required files
    required_files = [
        "summary_report.html",
        "kestrel/kestrel_result.tsv",
        "advntr/output_adVNTR_result.tsv",  # adVNTR outputs go in advntr/ subdirectory
    ]
    assert_required_files(output_dir, required_files)

    # 5. Validate adVNTR output
    validate_advntr_output(output_dir, test_case["advntr_assertions"])


def run_fastq_test_case(
    test_case: dict,
    runner: Callable[[Path, Optional[Path], str, Path, list[str]], int],
    output_dir: Path,
) -> None:
    """
    Universal FASTQ test case orchestration.

    Args:
        test_case: Test configuration dict
        runner: Function that executes the pipeline
                Signature: runner(fastq1, fastq2, reference, output_dir, extra_modules) -> exit_code
        output_dir: Output directory path

    Raises:
        AssertionError: If any validation fails
    """
    # 1. Extract test configuration
    fastq1 = Path(test_case["fastq1"])
    fastq2 = Path(test_case.get("fastq2", "")) if test_case.get("fastq2") else None
    reference = test_case.get("reference_assembly", "hg19")
    extra_modules: list[str] = []

    # Check for extra modules in cli_options
    cli_options = test_case.get("cli_options", [])
    if "--extra-modules" in cli_options:
        idx = cli_options.index("--extra-modules")
        if idx + 1 < len(cli_options):
            extra_modules = cli_options[idx + 1].split(",")

    # 2. Run pipeline via runner
    exit_code = runner(fastq1, fastq2, reference, output_dir, extra_modules)

    # 3. Assert success
    assert exit_code == 0, f"FASTQ pipeline failed with exit code {exit_code}"

    # 4. Validate expected files
    expected_files = test_case.get("expected_files", [])
    if expected_files:
        assert_required_files(output_dir, expected_files)
