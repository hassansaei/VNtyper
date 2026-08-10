"""
Shared test orchestration logic for VNtyper.

This module contains the CORE test logic that must be identical
for both local and Docker tests.

The ONLY difference between local and Docker tests is the
"runner" function that executes the pipeline.

This architecture guarantees 100% test identity.
"""

import os
from collections.abc import Callable
from pathlib import Path

from tests.helpers import (
    assert_required_files,
    validate_advntr_output,
    validate_coverage_output,
    validate_kestrel_output,
)

# Seconds allowed for a single adVNTR test, overriding pytest.ini's global 600 s.
#
# adVNTR's HMM genotyping is the slowest thing in the suite: ~5-6 min on a 32-core
# workstation and ~9 min locally, but GitHub's runners give 2 cores, so it needs several
# times that. 600 s passed locally and timed out on CI. The global timeout stays tight
# for everything else - this is the one genuine outlier - and 45 min still bounds a hang
# far below the 6 h job limit.
#
# Both adVNTR tests (docker and integration) import this rather than repeating the
# literal, so recalibrating for new CI hardware is a one-line change here.
ADVNTR_TIMEOUT_SECONDS = 2700

_FASTQ_OUTPUTS = (
    ("r1", "output_R1.fastq.gz"),
    ("r2", "output_R2.fastq.gz"),
    ("other", "output_other.fastq.gz"),
    ("single", "output_single.fastq.gz"),
)


def _declared_artifact_paths(test_case: dict, field: str, output_dir: Path) -> list[Path]:
    """Validate and resolve one declared artifact list beneath an output directory.

    Args:
        test_case: Integration case carrying the artifact declaration.
        field: Declaration field to validate.
        output_dir: Per-case pipeline output directory.

    Returns:
        Validated artifact paths rooted beneath ``output_dir``.

    Raises:
        ValueError: If the declaration is malformed, duplicated, or escapes the output directory.
    """
    test_name = str(test_case.get("test_name", "<unnamed>"))
    raw_paths = test_case.get(field, [])
    if not isinstance(raw_paths, list):
        raise ValueError(f"case={test_name} field={field} must be a list of relative artifact paths")
    root = output_dir.resolve()
    result: list[Path] = []
    seen: set[str] = set()
    for raw in raw_paths:
        if not isinstance(raw, str) or not raw:
            raise ValueError(f"case={test_name} field={field} contains an empty or non-string artifact path")
        relative = Path(raw)
        if relative.is_absolute() or ".." in relative.parts or raw in seen:
            raise ValueError(f"case={test_name} field={field} invalid artifact: {raw}")
        seen.add(raw)
        candidate = root / relative
        resolved_parent = candidate.parent.resolve(strict=False)
        if not resolved_parent.is_relative_to(root):
            raise ValueError(f"case={test_name} field={field} artifact escapes output_dir: {raw}")
        if field == "expected_present" and not candidate.resolve(strict=False).is_relative_to(root):
            raise ValueError(f"case={test_name} field={field} artifact escapes output_dir: {raw}")
        result.append(candidate)
    return result


def assert_declared_artifacts(test_case: dict, output_dir: Path) -> None:
    """Assert every case-declared output-relative artifact state.

    Args:
        test_case: Integration case carrying optional presence and absence declarations.
        output_dir: Per-case pipeline output directory.

    Raises:
        AssertionError: If a declared artifact has the wrong filesystem state.
        ValueError: If a declaration is malformed or escapes the output directory.
    """
    test_name = str(test_case.get("test_name", "<unnamed>"))
    present = _declared_artifact_paths(test_case, "expected_present", output_dir)
    absent = _declared_artifact_paths(test_case, "expected_absent", output_dir)
    overlap = set(present).intersection(absent)
    if overlap:
        raise ValueError(
            f"case={test_name} fields=expected_present,expected_absent overlap: {sorted(map(str, overlap))}"
        )
    failures = [f"case={test_name} field=expected_present missing: {path}" for path in present if not path.exists()]
    failures.extend(
        f"case={test_name} field=expected_absent unexpectedly present: {path}"
        for path in absent
        if os.path.lexists(path)
    )
    assert not failures, "Declared artifact mismatch:\n" + "\n".join(failures)


def mixed_layout_diagnostic(test_case: dict, output_dir: Path) -> str:
    """Render the exact fail-closed diagnostic declared for one mixed fixture.

    Args:
        test_case: Integration case carrying ``expected_mixed_fastq_records``.
        output_dir: Per-case pipeline output directory.

    Returns:
        The complete diagnostic emitted by ``route_converted_fastqs``.

    Raises:
        KeyError: If the case omits a required FASTQ record count.
    """
    counts = test_case["expected_mixed_fastq_records"]
    fastq_dir = output_dir / "fastq_bam_processing"
    details = ", ".join(f"{fastq_dir / filename}: {counts[key]} records" for key, filename in _FASTQ_OUTPUTS)
    return f"FASTQ layout 'mixed' cannot be consumed without dropping reads. Produced FASTQs: {details}"


def _assert_expected_exit(test_case: dict, exit_code: int, *, label: str) -> bool:
    """Assert the configured application exit and report whether success artefacts apply."""
    expected_exit = test_case.get("expected_exit_code", 0)
    assert exit_code == expected_exit, f"{label} pipeline exit: expected {expected_exit}, got {exit_code}"
    return expected_exit == 0


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

    # 3. Assert the declared outcome. Fail-closed cases stop before success artefacts.
    if not _assert_expected_exit(test_case, exit_code, label="BAM"):
        return

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
    runner: Callable[[Path, str, Path, list[str], list[str]], int],
    output_dir: Path,
) -> None:
    """
    Universal adVNTR test case orchestration.

    This function contains the COMPLETE test logic for adVNTR tests.
    Identical for local and Docker tests.

    Args:
        test_case: Test configuration dict
        runner: Function that executes the pipeline with extra modules
                Signature: runner(bam_file, reference, output_dir, extra_modules, extra_cli_options) -> exit_code
        output_dir: Output directory path

    Raises:
        AssertionError: If any validation fails
    """
    # 1. Extract test configuration
    bam_file = Path(test_case["bam"])
    reference = test_case["reference_assembly"]

    # Parse cli_options: separate --extra-modules from other CLI options
    cli_options = list(test_case.get("cli_options", []))
    extra_modules = ["advntr"]
    extra_cli_options: list[str] = []

    i = 0
    while i < len(cli_options):
        if cli_options[i] == "--extra-modules":
            # Extract module names, merging with default "advntr"
            if i + 1 < len(cli_options):
                for mod in cli_options[i + 1].split(","):
                    if mod not in extra_modules:
                        extra_modules.append(mod)
                i += 2
            else:
                i += 1
        else:
            extra_cli_options.append(cli_options[i])
            i += 1

    # 2. Run pipeline via runner
    exit_code = runner(bam_file, reference, output_dir, extra_modules, extra_cli_options)

    # 3. Assert the declared outcome. Fail-closed cases stop before success artefacts.
    if not _assert_expected_exit(test_case, exit_code, label="adVNTR"):
        return

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
    runner: Callable[[Path, Path | None, str, Path, list[str]], int],
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
