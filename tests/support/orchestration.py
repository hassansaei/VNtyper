"""Transport-independent VNtyper integration-test orchestration."""

import csv
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
    validate_kestrel_output,
)
from vntyper.scripts.coverage_stats import _BUILD_COMPARABLE_COLUMNS

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

    flags = [option.partition("=")[0] for option in request.cli_options if option.startswith("-")]
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


def _failure_diagnostic_fragments(test_case: dict[str, Any]) -> tuple[str, ...]:
    diagnostic = test_case.get("expected_diagnostic")
    if not isinstance(diagnostic, str) or not diagnostic:
        raise ValueError("expected_diagnostic must be a non-empty string")
    return (diagnostic,)


def _assert_expected_exit(
    test_case: dict[str, Any], result: PipelineRunResult, *, label: str
) -> tuple[bool, PipelineRunResult]:
    """Assert the configured application exit and report whether success artefacts apply."""
    if not isinstance(result, PipelineRunResult):
        raise TypeError("Pipeline runners must return PipelineRunResult with captured diagnostics")
    expected_exit = test_case.get("expected_exit_code", 0)
    assert result.exit_code == expected_exit, f"{label} pipeline exit: expected {expected_exit}, got {result.exit_code}"
    if expected_exit != 0:
        combined = f"{result.stdout}\n{result.stderr}"
        for fragment in _failure_diagnostic_fragments(test_case):
            assert fragment in combined, f"{label} pipeline missing declared causal diagnostic: {fragment}"
    return expected_exit == 0, result


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


def _merge_required_advntr(cli_options: tuple[str, ...]) -> tuple[str, ...]:
    """Return one canonical extra-modules option with adVNTR first."""
    retained: list[str] = []
    declared_modules: str | None = None
    index = 0
    while index < len(cli_options):
        option = cli_options[index]
        if option == "--extra-modules":
            if declared_modules is not None:
                raise ValueError("--extra-modules may be declared only once")
            if index + 1 >= len(cli_options) or cli_options[index + 1].startswith("-"):
                raise ValueError("--extra-modules requires a non-empty value")
            declared_modules = cli_options[index + 1]
            index += 2
            continue
        if option.startswith("--extra-modules="):
            if declared_modules is not None:
                raise ValueError("--extra-modules may be declared only once")
            declared_modules = option.partition("=")[2]
            index += 1
            continue
        retained.append(option)
        index += 1

    extras = [] if declared_modules is None else [module.strip() for module in declared_modules.split(",")]
    if any(not module for module in extras):
        raise ValueError("--extra-modules requires non-empty comma-separated module names")
    modules = ["advntr"]
    for module in extras:
        if module not in modules:
            modules.append(module)
    return (*retained, "--extra-modules", ",".join(modules))


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
_STRICT_CASE_FIELDS = _STRICT_FASTQ_FIELDS | {
    "expected_absent",
    "expected_fastq_records",
    "expected_selected_fastqs",
}
_STRICT_ADVNTR_CASE_FIELDS = _STRICT_ADVNTR_FIELDS | {
    "expected_absent",
    "expected_fastq_records",
    "expected_selected_fastqs",
}
_STRICT_KESTREL_FIELDS = frozenset(
    {
        "Confidence",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "Depth_Score",
    }
)


def _require_oracle_fields(test_case: dict[str, Any], required: frozenset[str]) -> None:
    missing = sorted(required.difference(test_case))
    if missing:
        raise ValueError(f"case={test_case.get('test_name', '<unnamed>')} missing strict oracle field(s): {missing}")


def _require_kestrel_fields(test_case: dict[str, Any]) -> None:
    missing = sorted(_STRICT_KESTREL_FIELDS.difference(test_case["kestrel_assertions"]))
    if missing:
        raise ValueError(f"kestrel_assertions missing characterized field(s): {missing}")


def _read_tsv_row(path: Path) -> dict[str, str]:
    assert path.is_file() and path.stat().st_size > 0, f"Result TSV is missing or empty: {path}"
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    assert len(rows) == 1, f"Expected exactly one result row in {path}, found {len(rows)}"
    return {key: value for key, value in rows[0].items() if key is not None and value is not None}


#: The coverage columns an integration case declares expected values for.
#:
#: Deliberately *not* every column in the TSV. #222 appended seven build-comparable columns,
#: and `scripts/integration_compatibility.py` treats each case's declarations as an
#: append-only contract that may not be mutated - so requiring the new columns to be declared
#: would invalidate every historical case rather than extend it. The new columns are checked
#: for presence here and asserted on properly in `tests/unit/test_coverage_stats.py`.
_DECLARED_COVERAGE_COLUMNS = tuple(column for column in COVERAGE_COLUMNS if column not in _BUILD_COMPARABLE_COLUMNS)


def _assert_coverage_values(test_case: dict[str, Any], output_dir: Path) -> None:
    expected = test_case["coverage_assertions"]
    assert set(expected) == set(_DECLARED_COVERAGE_COLUMNS), (
        f"coverage_assertions must declare every asserted coverage column: {list(_DECLARED_COVERAGE_COLUMNS)}"
    )
    actual = _read_tsv_row(output_dir / "coverage" / "coverage_summary.tsv")
    declared = {key: value for key, value in actual.items() if key in expected}
    assert declared == expected, f"Coverage values differ: expected {expected}, got {declared}"
    missing = [column for column in _BUILD_COMPARABLE_COLUMNS if column not in actual]
    assert not missing, f"coverage_summary.tsv is missing #222's columns: {missing}"


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


def _summary_expected_result(test_case: dict[str, Any], output_dir: Path, step: str) -> dict[str, str]:
    """Read one result-bearing summary oracle from its independently checked TSV."""
    paths = {
        "Coverage Calculation": output_dir / "coverage" / "coverage_summary.tsv",
        "Kestrel Genotyping": output_dir / "kestrel" / "kestrel_result.tsv",
        "adVNTR Genotyping": output_dir / "advntr" / "output_adVNTR_result.tsv",
        "Cross-Match Variant Comparison": output_dir / "advntr" / "cross_match_results.tsv",
    }
    expected = _read_tsv_row(paths[step])
    if "cross_match_assertions" in test_case:
        declared = test_case["cross_match_assertions"]
        assert set(declared) == {"comments", "data"}, "cross_match_assertions fields differ"
        assert isinstance(declared["comments"], list), "cross_match_assertions comments must be a list"
        assert isinstance(declared["data"], list) and len(declared["data"]) == 1
        cross_match = declared["data"][0]
        assert isinstance(cross_match, dict), "cross_match_assertions data must contain one result object"
        enrichment_fields = {
            "Kestrel Genotyping": {
                "Allele_Change": "Kestrel_Allele_Change",
                "Variant_Type": "Kestrel_Variant_Type",
            },
            "adVNTR Genotyping": {
                "Allele_Change": "Advntr_Allele_Change",
                "Variant_Type": "Advntr_Variant_Type",
            },
        }
        if step in enrichment_fields:
            for summary_field, cross_match_field in enrichment_fields[step].items():
                assert cross_match_field in cross_match, (
                    f"cross_match_assertions missing summary enrichment field: {cross_match_field}"
                )
                expected[summary_field] = cross_match[cross_match_field]
        elif step == "Cross-Match Variant Comparison":
            assert declared["data"] == [expected], "Cross-Match Variant Comparison TSV differs from its declared oracle"
    return expected


def _assert_summary_values(test_case: dict[str, Any], output_dir: Path, *, advntr: bool) -> None:
    expected = test_case["pipeline_summary_assertions"]
    required = {"Coverage Calculation", "Kestrel Genotyping"}
    if advntr:
        required.update({"adVNTR Genotyping", "Cross-Match Variant Comparison"})
    steps = _summary_steps(output_dir)
    assert isinstance(expected, dict), "pipeline_summary_assertions must be an explicit result contract"
    assert set(expected) == {"steps", "parsed_results"}, (
        "pipeline_summary_assertions must declare exactly steps and parsed_results"
    )
    expected_steps = expected["steps"]
    parsed_results = expected["parsed_results"]
    if (
        not isinstance(expected_steps, list)
        or not expected_steps
        or not all(isinstance(name, str) and name for name in expected_steps)
    ):
        raise ValueError("pipeline_summary_assertions.steps must be a non-empty list of step names")
    if not isinstance(parsed_results, list) or not all(isinstance(name, str) and name for name in parsed_results):
        raise ValueError("pipeline_summary_assertions.parsed_results must be a list of step names")
    assert list(steps) == expected_steps, (
        f"Pipeline summary step sequence differs: expected {expected_steps}, got {list(steps)}"
    )
    expected_result_order = [name for name in expected_steps if name in required]
    assert parsed_results == expected_result_order and set(parsed_results) == required, (
        f"pipeline_summary_assertions.parsed_results must declare exactly {sorted(required)}"
    )
    for name in parsed_results:
        parsed = steps[name]
        assert isinstance(parsed, dict) and set(parsed) == {"comments", "data"}, (
            f"Pipeline summary step {name} parsed_result must contain only comments and data"
        )
        assert isinstance(parsed["comments"], list), f"Pipeline summary step {name} comments must be a list"
        expected_result = _summary_expected_result(test_case, output_dir, name)
        data = parsed["data"]
        assert isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict), (
            f"Pipeline summary step {name} data must contain exactly one result object"
        )
        actual_result = data[0]
        assert actual_result == expected_result, f"Pipeline summary step {name} parsed_result differs"


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
    _require_kestrel_fields(test_case)
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
    _require_kestrel_fields(test_case)
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
    runner: Callable[[PipelineRequest], PipelineRunResult],
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
    result = runner(request)
    success, captured = _assert_expected_exit(test_case, result, label="BAM")
    if not success:
        return

    _require_oracle_fields(test_case, _STRICT_CASE_FIELDS)
    assert_read_set_routing(
        captured,
        expected_counts=test_case["expected_fastq_records"],
        expected_selected=tuple(test_case["expected_selected_fastqs"]),
    )
    validate_strict_fastq_success(test_case, output_dir)
    assert_declared_artifacts(test_case, output_dir)


def run_cram_test_case(
    test_case: dict,
    runner: Callable[[PipelineRequest], PipelineRunResult],
    output_dir: Path,
) -> None:
    """Run and validate one CRAM case identically across transports.

    Args:
        test_case: Test configuration from ``test_data_config.json``.
        runner: Typed local or Docker runner.
        output_dir: Per-case output directory.

    Raises:
        AssertionError: If any declared outcome differs.
        ValueError: If the strict success declaration is incomplete.
    """
    request = _request_from_case(test_case, "cram", (Path(test_case["cram"]),), output_dir)
    result = runner(request)
    success, captured = _assert_expected_exit(test_case, result, label="CRAM")
    if not success:
        return
    _require_oracle_fields(test_case, _STRICT_CASE_FIELDS)
    assert_read_set_routing(
        captured,
        expected_counts=test_case["expected_fastq_records"],
        expected_selected=tuple(test_case["expected_selected_fastqs"]),
    )
    validate_strict_fastq_success(test_case, output_dir)
    assert_declared_artifacts(test_case, output_dir)


def run_advntr_test_case(
    test_case: dict,
    runner: Callable[[PipelineRequest], PipelineRunResult],
    output_dir: Path,
) -> None:
    """Run and validate one adVNTR case identically across transports.

    Args:
        test_case: Test configuration from ``test_data_config.json``.
        runner: Typed local or Docker runner.
        output_dir: Per-case output directory.

    Raises:
        AssertionError: If any declared outcome differs.
        ValueError: If the extra-module declaration is ambiguous or malformed.
    """
    bam_file = Path(test_case["bam"])

    request = _request_from_case(test_case, "bam", (bam_file,), output_dir)
    request = PipelineRequest(
        input_kind=request.input_kind,
        input_paths=request.input_paths,
        reference_assembly=request.reference_assembly,
        output_dir=request.output_dir,
        threads=request.threads,
        log_level=request.log_level,
        cli_options=_merge_required_advntr(request.cli_options),
        reference_fasta=request.reference_fasta,
    )

    result = runner(request)
    success, captured = _assert_expected_exit(test_case, result, label="adVNTR")
    if not success:
        return

    _require_oracle_fields(test_case, _STRICT_ADVNTR_CASE_FIELDS)
    assert_read_set_routing(
        captured,
        expected_counts=test_case["expected_fastq_records"],
        expected_selected=tuple(test_case["expected_selected_fastqs"]),
    )
    validate_strict_advntr_success(test_case, output_dir)
    assert_declared_artifacts(test_case, output_dir)


def run_fastq_test_case(
    test_case: dict,
    runner: Callable[[PipelineRequest], PipelineRunResult],
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

    request = _request_from_case(
        test_case,
        "fastq",
        (fastq1,) if fastq2 is None else (fastq1, fastq2),
        output_dir,
    )

    result = runner(request)
    success, captured = _assert_expected_exit(test_case, result, label="FASTQ")
    assert success, "FASTQ declarations currently require a successful outcome"
    _require_oracle_fields(test_case, _STRICT_CASE_FIELDS)
    assert_read_set_routing(
        captured,
        expected_counts=test_case["expected_fastq_records"],
        expected_selected=tuple(test_case["expected_selected_fastqs"]),
    )
    validate_strict_fastq_success(test_case, output_dir)
    assert_declared_artifacts(test_case, output_dir)
