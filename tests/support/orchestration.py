"""Transport-independent VNtyper integration-test orchestration."""

import csv
import inspect
import json
import os
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from tests.helpers import (
    COVERAGE_COLUMNS,
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


@dataclass(frozen=True)
class PipelineRequest:
    """Canonical semantic request shared by local and Docker pipeline runners."""

    input_kind: Literal["bam", "cram", "fastq"]
    input_paths: tuple[Path, ...]
    reference_assembly: str
    output_dir: Path
    threads: int
    log_level: str
    cli_options: tuple[str, ...]
    reference_fasta: Path | None


@dataclass(frozen=True)
class PipelineRunResult:
    """Captured application outcome returned by every pipeline runner."""

    exit_code: int
    stdout: str
    stderr: str


_FASTQ_OUTPUTS = (
    ("r1", "output_R1.fastq.gz"),
    ("r2", "output_R2.fastq.gz"),
    ("other", "output_other.fastq.gz"),
    ("single", "output_single.fastq.gz"),
)

_LOG_LEVELS = frozenset({"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"})
_REQUEST_OWNED_OPTIONS = frozenset(
    {
        "-l",
        "--log-level",
        "--bam",
        "--cram",
        "--fastq1",
        "--fastq2",
        "--threads",
        "--reference-assembly",
        "-o",
        "--output",
        "--output-dir",
        "--reference-fasta",
    }
)


def build_pipeline_argv(request: PipelineRequest, path_mapper: Callable[[Path], str]) -> list[str]:
    """Build the only supported pipeline argv from a canonical request.

    Args:
        request: Transport-independent pipeline request.
        path_mapper: Function mapping each request path for one transport.

    Returns:
        Complete argv beginning with the ``vntyper`` executable.

    Raises:
        ValueError: If the request is ambiguous or invalid.
    """
    if isinstance(request.threads, bool) or request.threads <= 0:
        raise ValueError("PipelineRequest threads must be a positive integer")
    if request.log_level not in _LOG_LEVELS:
        raise ValueError(f"Unsupported PipelineRequest log level: {request.log_level}")
    expected_arity = {"bam": (1, 1), "cram": (1, 1), "fastq": (1, 2)}[request.input_kind]
    if not expected_arity[0] <= len(request.input_paths) <= expected_arity[1]:
        if request.input_kind == "fastq":
            raise ValueError("FASTQ requests require one or two input paths")
        raise ValueError(f"{request.input_kind.upper()} requests require exactly one input path")

    flags = [option for option in request.cli_options if option.startswith("-")]
    owned = sorted(_REQUEST_OWNED_OPTIONS.intersection(flags))
    if owned:
        raise ValueError(f"Options {owned} are owned by PipelineRequest fields")
    duplicate_flags = sorted(flag for flag in set(flags) if flags.count(flag) > 1)
    if duplicate_flags:
        raise ValueError(f"Duplicate CLI options are ambiguous: {duplicate_flags}")

    argv = ["vntyper", "-l", request.log_level, "pipeline"]
    if request.input_kind == "fastq":
        argv.extend(["--fastq1", path_mapper(request.input_paths[0])])
        if len(request.input_paths) == 2:
            argv.extend(["--fastq2", path_mapper(request.input_paths[1])])
    else:
        argv.extend([f"--{request.input_kind}", path_mapper(request.input_paths[0])])
    argv.extend(
        [
            "--threads",
            str(request.threads),
            "--reference-assembly",
            request.reference_assembly,
            "--output-dir",
            path_mapper(request.output_dir),
        ]
    )
    if request.reference_fasta is not None:
        argv.extend(["--reference-fasta", path_mapper(request.reference_fasta)])
    argv.extend(request.cli_options)
    return argv


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
    seen: set[Path] = set()
    for raw in raw_paths:
        if not isinstance(raw, str) or not raw:
            raise ValueError(f"case={test_name} field={field} contains an empty or non-string artifact path")
        relative = Path(raw)
        if relative.is_absolute() or ".." in relative.parts or relative in seen:
            raise ValueError(f"case={test_name} field={field} invalid artifact: {raw}")
        seen.add(relative)
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


def assert_declared_archive(test_case: dict, output_dir: Path) -> None:
    """Assert the optional default sibling ZIP state for a successful case.

    Args:
        test_case: Integration case carrying an optional archive declaration.
        output_dir: Per-case pipeline output directory.

    Raises:
        AssertionError: If the sibling archive has the wrong filesystem state.
        ValueError: If ``expected_archive`` is not a Boolean.
    """
    if test_case.get("expected_exit_code", 0) != 0 or "expected_archive" not in test_case:
        return
    case_id = str(test_case.get("test_name", "<unnamed>"))
    expected = test_case["expected_archive"]
    if not isinstance(expected, bool):
        raise ValueError(f"case={case_id} field=expected_archive must be a Boolean")
    archive = Path(f"{output_dir}.zip")
    if expected:
        assert archive.is_file(), f"case={case_id} field=expected_archive expected file is missing: {archive}"
    else:
        assert not os.path.lexists(archive), (
            f"case={case_id} field=expected_archive unexpected entry is present: {archive}"
        )


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


def current_declared_failure_diagnostic_adapter(test_case: dict[str, Any]) -> tuple[str, ...]:
    """Adapt only the current valid-mixed failure declarations to stable fragments.

    Task 3 cannot change the live manifest while the production routing fix is developed
    independently. Task 5 removes this seam when those declarations become successes.

    Args:
        test_case: Current live integration declaration.

    Returns:
        Stable causal and per-output count fragments expected in captured process output.

    Raises:
        ValueError: If the declaration is not one of the temporary nonzero mixed cases.
    """
    if test_case.get("expected_exit_code", 0) == 0 or "expected_mixed_fastq_records" not in test_case:
        raise ValueError("temporary current-failure adapter accepts only nonzero mixed-read declarations")
    counts = test_case["expected_mixed_fastq_records"]
    fragments = ["FASTQ layout 'mixed' cannot be consumed without dropping reads."]
    fragments.extend(f"{filename}: {counts[key]} records" for key, filename in _FASTQ_OUTPUTS)
    return tuple(fragments)


def _failure_diagnostic_fragments(test_case: dict[str, Any]) -> tuple[str, ...]:
    diagnostic = test_case.get("expected_diagnostic")
    if diagnostic is not None:
        if not isinstance(diagnostic, str) or not diagnostic:
            raise ValueError("expected_diagnostic must be a non-empty string")
        return (diagnostic,)
    return current_declared_failure_diagnostic_adapter(test_case)


def _coerce_run_result(result: PipelineRunResult | int) -> tuple[PipelineRunResult, bool]:
    """Keep old pure harness tests readable while real runners migrate atomically."""
    if isinstance(result, PipelineRunResult):
        return result, False
    return PipelineRunResult(exit_code=result, stdout="", stderr=""), True


def _assert_expected_exit(
    test_case: dict[str, Any], result: PipelineRunResult | int, *, label: str
) -> tuple[bool, PipelineRunResult]:
    """Assert the configured application exit and report whether success artefacts apply."""
    captured, legacy_int = _coerce_run_result(result)
    expected_exit = test_case.get("expected_exit_code", 0)
    assert captured.exit_code == expected_exit, (
        f"{label} pipeline exit: expected {expected_exit}, got {captured.exit_code}"
    )
    if expected_exit != 0 and not legacy_int:
        combined = f"{captured.stdout}\n{captured.stderr}"
        for fragment in _failure_diagnostic_fragments(test_case):
            assert fragment in combined, f"{label} pipeline missing declared causal diagnostic: {fragment}"
    return expected_exit == 0, captured


def _request_from_case(
    test_case: dict[str, Any],
    input_kind: Literal["bam", "cram", "fastq"],
    input_paths: tuple[Path, ...],
    output_dir: Path,
) -> PipelineRequest:
    """Create the canonical request while live declarations await Task 5 fields."""
    return PipelineRequest(
        input_kind=input_kind,
        input_paths=input_paths,
        reference_assembly=str(test_case.get("reference_assembly", "hg19")),
        output_dir=output_dir,
        threads=int(test_case.get("threads", 2)),
        log_level=str(test_case.get("log_level", "DEBUG")),
        cli_options=tuple(test_case.get("cli_options", ())),
        reference_fasta=Path(test_case["reference_fasta"]) if test_case.get("reference_fasta") else None,
    )


def _invoke_runner(
    runner: Callable[..., PipelineRunResult | int],
    request: PipelineRequest,
    legacy_args: tuple[Any, ...],
) -> PipelineRunResult | int:
    """Invoke typed runners while retaining compatibility with unowned legacy callers."""
    if isinstance(getattr(runner, "return_value", None), int):
        return runner(*legacy_args)
    try:
        inspect.signature(runner).bind(request)
    except TypeError:
        return runner(*legacy_args)
    except ValueError:
        pass
    return runner(request)


def assert_read_set_routing(
    result: PipelineRunResult,
    *,
    expected_counts: dict[str, int],
    expected_selected: tuple[str, ...],
) -> dict[str, Any]:
    """Assert exactly one canonical routing record and its complete declared values.

    Args:
        result: Captured pipeline output.
        expected_counts: Declared r1/r2/other/single record counts.
        expected_selected: Declared selected FASTQ basenames in execution order.

    Returns:
        Parsed canonical routing record.

    Raises:
        AssertionError: If evidence is missing, duplicated, malformed, noncanonical, or different.
    """
    marker = "READ_SET_ROUTING "
    payloads = [
        line.partition(marker)[2] for line in f"{result.stdout}\n{result.stderr}".splitlines() if marker in line
    ]
    assert len(payloads) == 1, f"Expected exactly one READ_SET_ROUTING record, found {len(payloads)}"
    payload = payloads[0]
    try:
        record = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise AssertionError(f"Malformed READ_SET_ROUTING JSON: {payload}") from exc
    assert json.dumps(record, sort_keys=True, separators=(",", ":")) == payload, (
        "READ_SET_ROUTING record is not canonical JSON"
    )
    assert set(record) == {"counts", "layout", "selected"}, "READ_SET_ROUTING fields differ from the contract"
    assert record["counts"] == expected_counts, "READ_SET_ROUTING counts differ from the declaration"
    assert record["selected"] == list(expected_selected), "READ_SET_ROUTING selected paths differ from the declaration"
    assert all(Path(path).name == path for path in record["selected"]), (
        "READ_SET_ROUTING selected paths must be stable basenames"
    )
    return dict(record)


_STRICT_FASTQ_FIELDS = frozenset(
    {
        "expected_files",
        "kestrel_assertions",
        "coverage_assertions",
        "pipeline_summary_assertions",
        "report_assertions",
        "expected_archive",
    }
)
_STRICT_ADVNTR_FIELDS = _STRICT_FASTQ_FIELDS | {"advntr_assertions", "cross_match_assertions"}


def _require_oracle_fields(test_case: dict[str, Any], required: frozenset[str]) -> None:
    missing = sorted(required.difference(test_case))
    if missing:
        raise ValueError(f"case={test_case.get('test_name', '<unnamed>')} missing strict oracle field(s): {missing}")


def _read_tsv_row(path: Path) -> dict[str, str]:
    assert path.is_file() and path.stat().st_size > 0, f"Result TSV is missing or empty: {path}"
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    assert len(rows) == 1, f"Expected exactly one result row in {path}, found {len(rows)}"
    return {key: value for key, value in rows[0].items() if key is not None and value is not None}


def _assert_coverage_values(test_case: dict[str, Any], output_dir: Path) -> None:
    expected = test_case["coverage_assertions"]
    assert set(expected) == set(COVERAGE_COLUMNS), (
        f"coverage_assertions must declare every coverage column: {list(COVERAGE_COLUMNS)}"
    )
    actual = _read_tsv_row(output_dir / "coverage" / "coverage_summary.tsv")
    assert actual == expected, f"Coverage values differ: expected {expected}, got {actual}"


def _summary_steps(output_dir: Path) -> dict[str, Any]:
    summary_path = output_dir / "pipeline_summary.json"
    assert summary_path.is_file() and summary_path.stat().st_size > 0, f"Pipeline summary is missing: {summary_path}"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    steps: dict[str, Any] = {}
    for entry in summary.get("steps", []):
        name = entry.get("step")
        assert isinstance(name, str), "Pipeline summary step is missing its name"
        assert name not in steps, f"Pipeline summary contains duplicate step: {name}"
        steps[name] = entry.get("parsed_result")
    return steps


def _assert_summary_values(test_case: dict[str, Any], output_dir: Path, *, advntr: bool) -> None:
    expected = test_case["pipeline_summary_assertions"]
    required = {"Coverage Calculation", "Kestrel Genotyping"}
    if advntr:
        required.add("adVNTR Genotyping")
    assert set(expected) == required, f"pipeline_summary_assertions must declare exactly {sorted(required)}"
    steps = _summary_steps(output_dir)
    for name, expected_result in expected.items():
        assert name in steps, f"Pipeline summary is missing step: {name}"
        assert steps[name] == expected_result, f"Pipeline summary step {name} differs"
    if advntr:
        cross_match = "Cross-Match Variant Comparison"
        assert cross_match in steps, f"Pipeline summary is missing step: {cross_match}"
        assert steps[cross_match] == test_case["cross_match_assertions"], f"Pipeline summary step {cross_match} differs"


def _assert_report_values(test_case: dict[str, Any], output_dir: Path) -> None:
    expected = test_case["report_assertions"]
    if not isinstance(expected, list) or not expected or not all(isinstance(item, str) and item for item in expected):
        raise ValueError("report_assertions must be a non-empty list of text fragments")
    report_path = output_dir / "summary_report.html"
    assert report_path.is_file() and report_path.stat().st_size > 0, f"Summary report is missing: {report_path}"
    report = report_path.read_text(encoding="utf-8")
    for fragment in expected:
        assert fragment in report, f"Summary report is missing declared text: {fragment}"


def validate_strict_fastq_success(test_case: dict[str, Any], output_dir: Path) -> None:
    """Validate the complete value-bearing FASTQ success schema.

    Args:
        test_case: Complete strict FASTQ outcome declaration.
        output_dir: Per-case pipeline output directory.

    Raises:
        ValueError: If the declaration omits a required oracle.
        AssertionError: If any declared artifact or value differs.
    """
    _require_oracle_fields(test_case, _STRICT_FASTQ_FIELDS)
    assert_required_files(output_dir, test_case["expected_files"])
    validate_kestrel_output(output_dir, test_case["kestrel_assertions"])
    _assert_coverage_values(test_case, output_dir)
    _assert_summary_values(test_case, output_dir, advntr=False)
    _assert_report_values(test_case, output_dir)
    assert_declared_archive(test_case, output_dir)


def validate_strict_advntr_success(test_case: dict[str, Any], output_dir: Path) -> None:
    """Validate Kestrel, adVNTR, cross-match, report, and archive as one outcome.

    Args:
        test_case: Complete strict adVNTR outcome declaration.
        output_dir: Per-case pipeline output directory.

    Raises:
        ValueError: If the declaration omits a required oracle.
        AssertionError: If any declared artifact or value differs.
    """
    _require_oracle_fields(test_case, _STRICT_ADVNTR_FIELDS)
    assert set(test_case["advntr_assertions"]) == {
        "VID",
        "State",
        "NumberOfSupportingReads",
        "MeanCoverage",
        "Pvalue",
    }, "advntr_assertions must declare every supported adVNTR field"
    assert_required_files(output_dir, test_case["expected_files"])
    validate_kestrel_output(output_dir, test_case["kestrel_assertions"])
    _assert_coverage_values(test_case, output_dir)
    validate_advntr_output(output_dir, test_case["advntr_assertions"])
    _assert_summary_values(test_case, output_dir, advntr=True)
    _assert_report_values(test_case, output_dir)
    assert_declared_archive(test_case, output_dir)


def run_bam_test_case(
    test_case: dict,
    runner: Callable[[PipelineRequest], PipelineRunResult | int] | Callable[[Path, str, Path], int],
    output_dir: Path,
) -> None:
    """Run and validate one BAM case identically across transports.

    Args:
        test_case: Test configuration from ``test_data_config.json``.
        runner: Typed local or Docker runner.
        output_dir: Per-case output directory.

    Raises:
        AssertionError: If any declared outcome differs.
    """
    request = _request_from_case(test_case, "bam", (Path(test_case["bam"]),), output_dir)
    result = _invoke_runner(runner, request, (request.input_paths[0], request.reference_assembly, output_dir))
    success, captured = _assert_expected_exit(test_case, result, label="BAM")
    if not success:
        return

    if "expected_fastq_records" in test_case:
        assert_read_set_routing(
            captured,
            expected_counts=test_case["expected_fastq_records"],
            expected_selected=tuple(test_case["expected_selected_fastqs"]),
        )

    required_files = [
        "summary_report.html",
        "kestrel/kestrel_result.tsv",
        "coverage/coverage_summary.tsv",
        "pipeline_summary.json",
    ]

    if test_case.get("check_igv_report"):
        required_files.append("igv_report.html")

    assert_required_files(output_dir, required_files)

    validate_kestrel_output(output_dir, test_case["kestrel_assertions"])
    coverage_metrics = validate_coverage_output(output_dir)
    assert coverage_metrics["mean_cov"] >= 0, "Coverage mean should be non-negative"
    assert_declared_artifacts(test_case, output_dir)
    assert_declared_archive(test_case, output_dir)


def run_advntr_test_case(
    test_case: dict,
    runner: Callable[[PipelineRequest], PipelineRunResult | int],
    output_dir: Path,
) -> None:
    """Run and validate one adVNTR case identically across transports.

    Args:
        test_case: Test configuration from ``test_data_config.json``.
        runner: Typed local or Docker runner.
        output_dir: Per-case output directory.

    Raises:
        AssertionError: If any declared outcome differs.
    """
    bam_file = Path(test_case["bam"])
    reference = test_case["reference_assembly"]

    cli_options = list(test_case.get("cli_options", []))
    extra_modules = ["advntr"]
    extra_cli_options: list[str] = []

    i = 0
    while i < len(cli_options):
        if cli_options[i] == "--extra-modules":
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

    request = _request_from_case(test_case, "bam", (bam_file,), output_dir)
    if "--extra-modules" not in request.cli_options:
        request = PipelineRequest(
            input_kind=request.input_kind,
            input_paths=request.input_paths,
            reference_assembly=request.reference_assembly,
            output_dir=request.output_dir,
            threads=request.threads,
            log_level=request.log_level,
            cli_options=(*request.cli_options, "--extra-modules", ",".join(extra_modules)),
            reference_fasta=request.reference_fasta,
        )

    result = _invoke_runner(runner, request, (bam_file, reference, output_dir, extra_modules, extra_cli_options))
    success, captured = _assert_expected_exit(test_case, result, label="adVNTR")
    if not success:
        return

    if "expected_fastq_records" in test_case:
        assert_read_set_routing(
            captured,
            expected_counts=test_case["expected_fastq_records"],
            expected_selected=tuple(test_case["expected_selected_fastqs"]),
        )

    required_files = [
        "summary_report.html",
        "kestrel/kestrel_result.tsv",
        "advntr/output_adVNTR_result.tsv",  # adVNTR outputs go in advntr/ subdirectory
    ]
    assert_required_files(output_dir, required_files)

    validate_advntr_output(output_dir, test_case["advntr_assertions"])
    assert_declared_archive(test_case, output_dir)


def run_fastq_test_case(
    test_case: dict,
    runner: Callable[[PipelineRequest], PipelineRunResult | int],
    output_dir: Path,
) -> None:
    """Run and validate one FASTQ case identically across transports.

    Args:
        test_case: Test configuration from ``test_data_config.json``.
        runner: Typed local or Docker runner.
        output_dir: Per-case output directory.

    Raises:
        AssertionError: If any declared outcome differs.
    """
    fastq1 = Path(test_case["fastq1"])
    fastq2 = Path(test_case.get("fastq2", "")) if test_case.get("fastq2") else None
    reference = test_case.get("reference_assembly", "hg19")
    extra_modules: list[str] = []

    cli_options = test_case.get("cli_options", [])
    if "--extra-modules" in cli_options:
        idx = cli_options.index("--extra-modules")
        if idx + 1 < len(cli_options):
            extra_modules = cli_options[idx + 1].split(",")

    request = _request_from_case(
        test_case,
        "fastq",
        (fastq1,) if fastq2 is None else (fastq1, fastq2),
        output_dir,
    )

    result = _invoke_runner(runner, request, (fastq1, fastq2, reference, output_dir, extra_modules))
    success, _captured = _assert_expected_exit(test_case, result, label="FASTQ")
    assert success, "FASTQ declarations currently require a successful outcome"

    expected_files = test_case.get("expected_files", [])
    if expected_files:
        assert_required_files(output_dir, expected_files)
