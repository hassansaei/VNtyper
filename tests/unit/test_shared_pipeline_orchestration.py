"""Unit contracts for transport-independent pipeline orchestration."""

import json
from dataclasses import FrozenInstanceError
from pathlib import Path
from typing import Any

import pytest

from tests.support import orchestration

pytestmark = pytest.mark.unit


def _request(tmp_path: Path, **changes: object) -> orchestration.PipelineRequest:
    values: dict[str, Any] = {
        "input_kind": "bam",
        "input_paths": (Path("tests/data/nested/sample.bam"),),
        "reference_assembly": "hg19",
        "output_dir": tmp_path / "case",
        "threads": 2,
        "log_level": "DEBUG",
        "cli_options": ("--fast-mode", "--archive-results"),
        "reference_fasta": None,
    }
    values.update(changes)
    return orchestration.PipelineRequest(**values)


def test_pipeline_request_and_result_are_frozen_transport_values(tmp_path: Path) -> None:
    """A runner must not be able to mutate its canonical request or captured result."""
    assert hasattr(orchestration, "PipelineRequest")
    assert hasattr(orchestration, "PipelineRunResult")

    request = orchestration.PipelineRequest(
        input_kind="bam",
        input_paths=(Path("tests/data/nested/sample.bam"),),
        reference_assembly="hg19",
        output_dir=tmp_path / "case",
        threads=2,
        log_level="DEBUG",
        cli_options=("--fast-mode", "--archive-results"),
        reference_fasta=None,
    )
    result = orchestration.PipelineRunResult(exit_code=0, stdout="out", stderr="err")

    with pytest.raises(FrozenInstanceError):
        request.threads = 4  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        result.stderr = "changed"  # type: ignore[misc]


def test_build_pipeline_argv_owns_every_non_path_argument(tmp_path: Path) -> None:
    """Dropping fast/archive options or restoring a hidden thread default must change this literal argv."""
    assert hasattr(orchestration, "build_pipeline_argv")
    request = _request(tmp_path)

    argv = orchestration.build_pipeline_argv(request, lambda path: f"mapped:{path}")

    assert argv == [
        "vntyper",
        "-l",
        "DEBUG",
        "pipeline",
        "--bam",
        "mapped:tests/data/nested/sample.bam",
        "--threads",
        "2",
        "--reference-assembly",
        "hg19",
        "--output-dir",
        f"mapped:{tmp_path / 'case'}",
        "--fast-mode",
        "--archive-results",
    ]


def test_build_pipeline_argv_maps_fastq_and_explicit_cram_reference(tmp_path: Path) -> None:
    """The sole builder must cover both FASTQ arities and CRAM's optional reference path."""
    fastq = _request(
        tmp_path,
        input_kind="fastq",
        input_paths=(Path("reads/R1.fastq.gz"), Path("reads/R2.fastq.gz")),
        cli_options=(),
    )
    cram = _request(
        tmp_path,
        input_kind="cram",
        input_paths=(Path("reads/sample.cram"),),
        reference_fasta=Path("reference/alignment/chr1.hg19.fa"),
        cli_options=(),
    )

    assert orchestration.build_pipeline_argv(fastq, str)[4:8] == [
        "--fastq1",
        "reads/R1.fastq.gz",
        "--fastq2",
        "reads/R2.fastq.gz",
    ]
    assert orchestration.build_pipeline_argv(cram, str)[4:] == [
        "--cram",
        "reads/sample.cram",
        "--threads",
        "2",
        "--reference-assembly",
        "hg19",
        "--output-dir",
        str(tmp_path / "case"),
        "--reference-fasta",
        "reference/alignment/chr1.hg19.fa",
    ]


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"threads": 0}, "positive"),
        ({"threads": True}, "positive"),
        ({"log_level": "TRACE"}, "log level"),
        ({"input_paths": ()}, "exactly one"),
        ({"input_kind": "fastq", "input_paths": (Path("a"), Path("b"), Path("c"))}, "one or two"),
        ({"cli_options": ("--threads", "9")}, "owned by PipelineRequest"),
        ({"cli_options": ("--threads=9",)}, "owned by PipelineRequest"),
        ({"cli_options": ("--output-dir=/tmp/shared",)}, "owned by PipelineRequest"),
        ({"cli_options": ("--log-level=INFO",)}, "owned by PipelineRequest"),
        ({"cli_options": ("--fast-mode", "--fast-mode")}, "Duplicate"),
    ],
)
def test_build_pipeline_argv_rejects_ambiguous_requests(
    tmp_path: Path,
    changes: dict[str, object],
    message: str,
) -> None:
    """Invalid arity and transport-owned CLI flags must fail before either transport executes."""
    with pytest.raises(ValueError, match=message):
        orchestration.build_pipeline_argv(_request(tmp_path, **changes), str)


def test_bam_orchestration_sends_one_canonical_request_and_checks_failure_text(tmp_path: Path) -> None:
    """Exit one alone must not pass, and the runner must receive all semantic options in one value."""
    case = {
        "test_name": "refusal",
        "bam": "tests/data/nested/sample.bam",
        "reference_assembly": "hg19",
        "threads": 7,
        "log_level": "WARNING",
        "cli_options": ["--fast-mode", "--archive-results"],
        "expected_exit_code": 1,
        "expected_diagnostic": "mate outputs are inconsistent",
    }
    received: list[orchestration.PipelineRequest] = []

    def runner(request: orchestration.PipelineRequest) -> orchestration.PipelineRunResult:
        received.append(request)
        return orchestration.PipelineRunResult(1, "", "ERROR mate outputs are inconsistent")

    orchestration.run_bam_test_case(case, runner, tmp_path / "refusal")

    assert received == [
        orchestration.PipelineRequest(
            input_kind="bam",
            input_paths=(Path("tests/data/nested/sample.bam"),),
            reference_assembly="hg19",
            output_dir=tmp_path / "refusal",
            threads=7,
            log_level="WARNING",
            cli_options=("--fast-mode", "--archive-results"),
            reference_fasta=None,
        )
    ]


def test_bam_orchestration_rejects_a_bare_integer_result(tmp_path: Path) -> None:
    """An exit code without captured diagnostics can never satisfy a declared failure."""
    case = {
        "test_name": "legacy-single-end",
        "bam": "tests/data/single.bam",
        "reference_assembly": "hg19",
        "expected_exit_code": 1,
        "expected_diagnostic": "known refusal",
    }

    def runner(_request: orchestration.PipelineRequest) -> int:
        return 1

    with pytest.raises(TypeError, match="PipelineRunResult"):
        orchestration.run_bam_test_case(case, runner, tmp_path)  # type: ignore[arg-type]


def test_expected_failure_diagnostic_is_found_across_captured_streams(tmp_path: Path) -> None:
    """Local stderr and Docker's combined stdout must be equivalent diagnostic sources."""
    case = {
        "test_name": "refusal",
        "bam": "tests/data/sample.bam",
        "reference_assembly": "hg19",
        "expected_exit_code": 1,
        "expected_diagnostic": "causal phrase",
    }
    for result in (
        orchestration.PipelineRunResult(1, "causal phrase", ""),
        orchestration.PipelineRunResult(1, "", "causal phrase"),
    ):

        def runner(
            _request: orchestration.PipelineRequest,
            result: orchestration.PipelineRunResult = result,
        ) -> orchestration.PipelineRunResult:
            return result

        orchestration.run_bam_test_case(case, runner, tmp_path)

    with pytest.raises(AssertionError, match="causal phrase"):
        orchestration.run_bam_test_case(
            case,
            lambda _request: orchestration.PipelineRunResult(1, "", "unrelated"),
            tmp_path,
        )


@pytest.mark.parametrize(
    "module_options",
    [["--extra-modules", "shark"], ["--extra-modules=shark,advntr,shark"]],
)
def test_advntr_orchestration_merges_an_additional_module_into_the_typed_request(
    tmp_path: Path, module_options: list[str]
) -> None:
    """An additional module must not displace the adVNTR stage owned by this orchestration path."""
    case = {
        "test_name": "advntr-with-shark",
        "bam": "tests/data/sample.bam",
        "reference_assembly": "hg19",
        "cli_options": ["--fast-mode", *module_options, "--archive-results"],
        "expected_exit_code": 1,
        "expected_diagnostic": "known refusal",
    }
    received: list[orchestration.PipelineRequest] = []

    def runner(request: orchestration.PipelineRequest) -> orchestration.PipelineRunResult:
        received.append(request)
        return orchestration.PipelineRunResult(1, "", "known refusal")

    orchestration.run_advntr_test_case(case, runner, tmp_path)

    assert received[0].cli_options == (
        "--fast-mode",
        "--archive-results",
        "--extra-modules",
        "advntr,shark",
    )


@pytest.mark.parametrize(
    "cli_options",
    [
        ["--extra-modules"],
        ["--extra-modules", "--fast-mode"],
        ["--extra-modules", "shark", "--extra-modules=advntr"],
        ["--extra-modules="],
    ],
)
def test_advntr_orchestration_rejects_ambiguous_module_declarations(tmp_path: Path, cli_options: list[str]) -> None:
    """A missing, empty, or duplicate module value must fail before transport execution."""
    case = {
        "bam": "tests/data/sample.bam",
        "reference_assembly": "hg19",
        "cli_options": cli_options,
    }

    with pytest.raises(ValueError, match="extra-modules"):
        orchestration.run_advntr_test_case(
            case,
            lambda _request: pytest.fail("invalid request reached the transport"),
            tmp_path,
        )


def _routing_line(record: dict[str, Any]) -> str:
    import json

    return "READ_SET_ROUTING " + json.dumps(record, sort_keys=True, separators=(",", ":"))


def test_read_set_routing_requires_one_canonical_exact_record() -> None:
    """Every non-empty converted FASTQ must be represented once and in Kestrel execution order."""
    record = {
        "counts": {"other": 0, "r1": 7, "r2": 7, "single": 1},
        "layout": "mixed",
        "selected": ["output_R1.fastq.gz", "output_R2.fastq.gz", "output_single.fastq.gz"],
    }
    result = orchestration.PipelineRunResult(0, f"prefix INFO {_routing_line(record)}\n", "")

    parsed = orchestration.assert_read_set_routing(
        result,
        expected_counts={"r1": 7, "r2": 7, "other": 0, "single": 1},
        expected_selected=("output_R1.fastq.gz", "output_R2.fastq.gz", "output_single.fastq.gz"),
    )

    assert parsed == record


@pytest.mark.parametrize(
    "output",
    [
        "no routing record",
        'READ_SET_ROUTING {"counts":',
        (
            'READ_SET_ROUTING {"counts":{"other":0,"r1":7,"r2":7,"single":1},"layout":"mixed",'
            '"selected":["output_R1.fastq.gz","output_R2.fastq.gz","output_single.fastq.gz"]}\n'
        )
        * 2,
        (
            'READ_SET_ROUTING {"layout":"mixed","counts":{"other":0,"r1":7,"r2":7,"single":1},'
            '"selected":["output_R1.fastq.gz","output_R2.fastq.gz","output_single.fastq.gz"]}'
        ),
    ],
)
def test_read_set_routing_rejects_missing_duplicate_malformed_or_noncanonical_records(output: str) -> None:
    """The evidence parser must fail closed instead of accepting a plausible log fragment."""
    with pytest.raises(AssertionError, match="READ_SET_ROUTING"):
        orchestration.assert_read_set_routing(
            orchestration.PipelineRunResult(0, output, ""),
            expected_counts={"r1": 7, "r2": 7, "other": 0, "single": 1},
            expected_selected=("output_R1.fastq.gz", "output_R2.fastq.gz", "output_single.fastq.gz"),
        )


@pytest.mark.parametrize(
    ("counts", "selected"),
    [
        ({"r1": 7, "r2": 7, "other": 0, "single": 0}, ("output_R1.fastq.gz", "output_R2.fastq.gz")),
        (
            {"r1": 7, "r2": 7, "other": 0, "single": 1},
            ("output_R2.fastq.gz", "output_R1.fastq.gz", "output_single.fastq.gz"),
        ),
        (
            {"r1": 7, "r2": 7, "other": 0, "single": 1},
            ("output_R1.fastq.gz", "output_R2.fastq.gz"),
        ),
        (
            {"r1": 7, "r2": 7, "other": 0, "single": 1},
            ("output_R1.fastq.gz", "output_R2.fastq.gz", "output_single.fastq.gz", "extra.fastq.gz"),
        ),
    ],
)
def test_read_set_routing_rejects_wrong_counts_dropped_reordered_or_extra_selection(
    counts: dict[str, int], selected: tuple[str, ...]
) -> None:
    """A genotype cannot compensate for lost, reordered, or invented sequence inputs."""
    actual = {
        "counts": {"other": 0, "r1": 7, "r2": 7, "single": 1},
        "layout": "mixed",
        "selected": ["output_R1.fastq.gz", "output_R2.fastq.gz", "output_single.fastq.gz"],
    }
    with pytest.raises(AssertionError):
        orchestration.assert_read_set_routing(
            orchestration.PipelineRunResult(0, _routing_line(actual), ""),
            expected_counts=counts,
            expected_selected=selected,
        )


def _write_strict_success_tree(output_dir: Path, *, advntr: bool = False) -> None:
    (output_dir / "kestrel").mkdir(parents=True)
    (output_dir / "coverage").mkdir()
    (output_dir / "kestrel" / "kestrel_result.tsv").write_text(
        (
            "Confidence\tEstimated_Depth_AlternateVariant\t"
            "Estimated_Depth_Variant_ActiveRegion\tDepth_Score\n"
            "High_Precision*\t125\t1024\t0.1220703125\n"
        ),
        encoding="utf-8",
    )
    (output_dir / "coverage" / "coverage_summary.tsv").write_text(
        "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\t"
        # #222's build-comparable columns, as a run with no `vntr_array_coords` records them.
        "vntr_array_length\tvntr_array_depth_sum\tvntr_array_depth_sum_per_unit_length\t"
        "depth_sum_reference_length\tvntr_flank_bases\tvntr_flank_mean_depth\tdepth_counting_policy\tcoverage_qc\n"
        "566.92\t593.00\t300.71\t23\t1062\t1501\t0\t0.00\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tPASS\n",
        encoding="utf-8",
    )
    coverage_result = {
        "mean": "566.92",
        "median": "593.00",
        "stdev": "300.71",
        "min": "23",
        "max": "1062",
        "region_length": "1501",
        "uncovered_bases": "0",
        "percent_uncovered": "0.00",
        # #222's columns as a run with no `vntr_array_coords` records them: the literal
        # not-measured token, never 0 - zero would read as "no coverage was seen".
        "vntr_array_length": "NA",
        "vntr_array_depth_sum": "NA",
        "vntr_array_depth_sum_per_unit_length": "NA",
        "depth_sum_reference_length": "NA",
        "vntr_flank_bases": "NA",
        "vntr_flank_mean_depth": "NA",
        "depth_counting_policy": "NA",
        "coverage_qc": "PASS",
    }
    kestrel_result = {
        "Confidence": "High_Precision*",
        "Estimated_Depth_AlternateVariant": "125",
        "Estimated_Depth_Variant_ActiveRegion": "1024",
        "Depth_Score": "0.1220703125",
    }
    steps: list[dict[str, Any]] = [
        {"step": "Coverage Calculation", "parsed_result": {"comments": [], "data": [coverage_result]}},
        {"step": "Kestrel Genotyping", "parsed_result": {"comments": [], "data": [kestrel_result]}},
    ]
    if advntr:
        steps[1]["parsed_result"]["data"][0].update({"Variant_Type": "Insertion", "Allele_Change": "G"})
        (output_dir / "advntr").mkdir()
        (output_dir / "advntr" / "output_adVNTR_result.tsv").write_text(
            "VID\tVariant\tNumberOfSupportingReads\tMeanCoverage\tPvalue\n"
            "25561\tI22_4_G_LEN1\t39\t70.3333333333\t5.774455097259999e-59\n",
            encoding="utf-8",
        )
        (output_dir / "advntr" / "cross_match_results.tsv").write_text(
            "Kestrel_Allele_Change\tKestrel_Variant_Type\tAdvntr_Allele_Change\tAdvntr_Variant_Type\tMatch\n"
            "G\tInsertion\tG\tInsertion\tYes\n",
            encoding="utf-8",
        )
        steps.extend(
            [
                {
                    "step": "adVNTR Genotyping",
                    "parsed_result": {
                        "comments": [],
                        "data": [
                            {
                                "VID": "25561",
                                "Variant": "I22_4_G_LEN1",
                                "NumberOfSupportingReads": "39",
                                "MeanCoverage": "70.3333333333",
                                "Pvalue": "5.774455097259999e-59",
                                "Variant_Type": "Insertion",
                                "Allele_Change": "G",
                            }
                        ],
                    },
                },
                {
                    "step": "Cross-Match Variant Comparison",
                    "parsed_result": {
                        "comments": [],
                        "data": [
                            {
                                "Kestrel_Allele_Change": "G",
                                "Kestrel_Variant_Type": "Insertion",
                                "Advntr_Allele_Change": "G",
                                "Advntr_Variant_Type": "Insertion",
                                "Match": "Yes",
                            }
                        ],
                    },
                },
            ]
        )
    (output_dir / "pipeline_summary.json").write_text(json.dumps({"steps": steps}), encoding="utf-8")
    (output_dir / "summary_report.html").write_text(
        "<html>Kestrel detected a high-precision pathogenic variant.</html>", encoding="utf-8"
    )


def _strict_fastq_case() -> dict[str, Any]:
    return {
        "test_name": "strict-fastq",
        "expected_files": [
            "summary_report.html",
            "kestrel/kestrel_result.tsv",
            "coverage/coverage_summary.tsv",
            "pipeline_summary.json",
        ],
        "kestrel_assertions": {
            "Confidence": "High_Precision*",
            "Estimated_Depth_AlternateVariant": 125,
            "Estimated_Depth_Variant_ActiveRegion": 1024,
            "Depth_Score": 0.1220703125,
        },
        "coverage_assertions": {
            "mean": "566.92",
            "median": "593.00",
            "stdev": "300.71",
            "min": "23",
            "max": "1062",
            "region_length": "1501",
            "uncovered_bases": "0",
            "percent_uncovered": "0.00",
            "coverage_qc": "PASS",
        },
        "pipeline_summary_assertions": {
            "steps": ["Coverage Calculation", "Kestrel Genotyping"],
            "parsed_results": ["Coverage Calculation", "Kestrel Genotyping"],
        },
        "report_assertions": ["Kestrel detected a high-precision pathogenic variant."],
        "expected_archive": False,
    }


def test_strict_fastq_success_validates_values_summary_report_and_archive(tmp_path: Path) -> None:
    """A FASTQ success needs value-bearing artifacts, not merely an exit code and two files."""
    output_dir = tmp_path / "fastq"
    _write_strict_success_tree(output_dir)

    orchestration.validate_strict_fastq_success(_strict_fastq_case(), output_dir)

    changed = _strict_fastq_case()
    changed["kestrel_assertions"] = {
        "Confidence": "Negative",
        "Estimated_Depth_AlternateVariant": 125,
        "Estimated_Depth_Variant_ActiveRegion": 1024,
        "Depth_Score": 0.1220703125,
    }
    with pytest.raises(AssertionError, match="Confidence"):
        orchestration.validate_strict_fastq_success(changed, output_dir)


@pytest.mark.parametrize("missing", ["Estimated_Depth_Variant_ActiveRegion", "Depth_Score"])
def test_strict_fastq_success_requires_every_characterized_kestrel_field(tmp_path: Path, missing: str) -> None:
    """The strict oracle must not silently weaken either characterized depth value."""
    output_dir = tmp_path / "fastq"
    _write_strict_success_tree(output_dir)
    case = _strict_fastq_case()
    del case["kestrel_assertions"][missing]

    with pytest.raises(ValueError, match=missing):
        orchestration.validate_strict_fastq_success(case, output_dir)


@pytest.mark.parametrize(
    "missing",
    [
        "expected_files",
        "kestrel_assertions",
        "coverage_assertions",
        "pipeline_summary_assertions",
        "report_assertions",
        "expected_archive",
    ],
)
def test_strict_fastq_success_rejects_incomplete_oracles(tmp_path: Path, missing: str) -> None:
    """No missing value or artifact field may degrade the strict path into a smoke check."""
    output_dir = tmp_path / "fastq"
    _write_strict_success_tree(output_dir)
    case = _strict_fastq_case()
    del case[missing]

    with pytest.raises(ValueError, match=missing):
        orchestration.validate_strict_fastq_success(case, output_dir)


def test_strict_advntr_success_adds_exact_module_and_cross_match_oracles(tmp_path: Path) -> None:
    """adVNTR success must bind its TSV and cross-match result to the same pipeline summary."""
    output_dir = tmp_path / "advntr"
    _write_strict_success_tree(output_dir, advntr=True)
    case = _strict_fastq_case()
    case.update(
        {
            "test_name": "strict-advntr",
            "expected_files": [
                *case["expected_files"],
                "advntr/output_adVNTR_result.tsv",
            ],
            "pipeline_summary_assertions": {
                "steps": [
                    *case["pipeline_summary_assertions"]["steps"],
                    "adVNTR Genotyping",
                    "Cross-Match Variant Comparison",
                ],
                "parsed_results": [
                    *case["pipeline_summary_assertions"]["parsed_results"],
                    "adVNTR Genotyping",
                    "Cross-Match Variant Comparison",
                ],
            },
            "advntr_assertions": {
                "VID": "25561",
                "State": "I22_4_G_LEN1",
                "NumberOfSupportingReads": 39,
                "MeanCoverage": 70.3333333333,
                "Pvalue": 5.774455097259999e-59,
            },
            "cross_match_assertions": {
                "comments": [],
                "data": [
                    {
                        "Kestrel_Allele_Change": "G",
                        "Kestrel_Variant_Type": "Insertion",
                        "Advntr_Allele_Change": "G",
                        "Advntr_Variant_Type": "Insertion",
                        "Match": "Yes",
                    }
                ],
            },
        }
    )

    orchestration.validate_strict_advntr_success(case, output_dir)

    case["cross_match_assertions"]["data"][0]["Match"] = "No"
    with pytest.raises(AssertionError, match="Cross-Match Variant Comparison"):
        orchestration.validate_strict_advntr_success(case, output_dir)


@pytest.mark.parametrize(
    ("step", "field"),
    [
        ("Kestrel Genotyping", "Variant_Type"),
        ("Kestrel Genotyping", "Allele_Change"),
        ("Kestrel Genotyping", "Unexpected_Stable_Field"),
        ("adVNTR Genotyping", "Variant_Type"),
        ("adVNTR Genotyping", "Allele_Change"),
        ("adVNTR Genotyping", "Unexpected_Stable_Field"),
    ],
)
def test_strict_advntr_summary_rejects_mutated_enriched_cross_match_values(
    tmp_path: Path, step: str, field: str
) -> None:
    """Cross-match enrichment must remain bound to its independent declared result."""
    output_dir = tmp_path / "advntr"
    _write_strict_success_tree(output_dir, advntr=True)
    case = _strict_fastq_case()
    case.update(
        {
            "expected_files": [*case["expected_files"], "advntr/output_adVNTR_result.tsv"],
            "pipeline_summary_assertions": {
                "steps": [
                    *case["pipeline_summary_assertions"]["steps"],
                    "adVNTR Genotyping",
                    "Cross-Match Variant Comparison",
                ],
                "parsed_results": [
                    *case["pipeline_summary_assertions"]["parsed_results"],
                    "adVNTR Genotyping",
                    "Cross-Match Variant Comparison",
                ],
            },
            "advntr_assertions": {
                "VID": "25561",
                "State": "I22_4_G_LEN1",
                "NumberOfSupportingReads": 39,
                "MeanCoverage": 70.3333333333,
                "Pvalue": 5.774455097259999e-59,
            },
            "cross_match_assertions": {
                "comments": [],
                "data": [
                    {
                        "Kestrel_Allele_Change": "G",
                        "Kestrel_Variant_Type": "Insertion",
                        "Advntr_Allele_Change": "G",
                        "Advntr_Variant_Type": "Insertion",
                        "Match": "Yes",
                    }
                ],
            },
        }
    )
    summary_path = output_dir / "pipeline_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    result = next(entry["parsed_result"]["data"][0] for entry in summary["steps"] if entry["step"] == step)
    result[field] = "mutated"
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(AssertionError, match=step):
        orchestration.validate_strict_advntr_success(case, output_dir)
