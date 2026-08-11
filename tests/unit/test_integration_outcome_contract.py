"""Pure checks for strict, value-bearing real-integration outcome declarations."""

from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests import parametrization
from tests.integration import test_pipeline_integration as local_integration
from tests.parametrization import get_advntr_test_cases, get_bam_test_cases, get_fastq_test_cases, load_test_config
from tests.support import orchestration

pytestmark = pytest.mark.unit


EXPECTED_BAM_ROUTING = {
    "example_b178_hg19_subset_fast": (14_690, 14_690, 0, 1),
    "example_a5c1_hg19_subset_fast": (20_888, 20_888, 0, 3),
    "example_66bf_hg19_subset_fast": (19_841, 19_841, 0, 3),
    "example_7a61_hg19_subset_fast": (2_359, 2_359, 0, 8),
    "example_dfc3_hg19_subset_fast": (31_596, 31_596, 0, 53),
    "example_dfc3_GRCh37_fast": (31_593, 31_593, 0, 2),
    "example_dfc3_GRCh38_fast": (31_603, 31_603, 0, 2),
    "example_dfc3_hg38_ensembl_fast": (31_603, 31_603, 0, 2),
    "example_40cf_hg38_subset_fast_gdp_guard": (3_474, 3_474, 0, 93),
    "example_b178_hg19_subset_default": (16_929, 16_929, 0, 1),
    "example_40cf_hg38_subset_default": (19_492, 19_492, 0, 93),
}


def _negative_case() -> dict[str, Any]:
    return {
        "bam": "fixture.bam",
        "reference_assembly": "hg19",
        "expected_exit_code": 1,
        "expected_diagnostic": "mate outputs are inconsistent",
        "expected_mixed_fastq_records": {"r1": 7, "r2": 7, "other": 0, "single": 1},
    }


def _negative_result() -> orchestration.PipelineRunResult:
    return orchestration.PipelineRunResult(1, "", "mate outputs are inconsistent")


def test_every_bam_fixture_is_a_real_success_with_exact_measured_routing() -> None:
    """Catch restoring the temporary mixed-layout refusal or losing default-mode controls."""
    configured = {}
    for case in get_bam_test_cases():
        assert case["expected_exit_code"] == 0
        assert case["threads"] == 2
        assert case["log_level"] == "DEBUG"
        counts = case["expected_fastq_records"]
        configured[case["test_name"]] = (counts["r1"], counts["r2"], counts["other"], counts["single"])
        assert case["expected_selected_fastqs"] == [
            "output_R1.fastq.gz",
            "output_R2.fastq.gz",
            "output_single.fastq.gz",
        ]

    assert configured == EXPECTED_BAM_ROUTING


def test_declared_failures_do_not_retain_unreachable_success_expectations() -> None:
    """A case that stops at routing must not promise downstream reports or genotype fields."""
    success_only = {
        "expected_archive",
        "expected_absent",
        "expected_files",
        "expected_present",
        "kestrel_assertions",
        "check_igv_report",
        "expected_vcf",
        "advntr_assertions",
    }
    stale = {}
    single_end_cases = load_test_config()["integration_tests"]["single_end_bam_tests"]
    for case in [*get_bam_test_cases(), *get_advntr_test_cases(), *get_fastq_test_cases(), *single_end_cases]:
        if case.get("expected_exit_code", 0) != 0:
            retained = sorted(success_only.intersection(case))
            if retained:
                stale[case["test_name"]] = retained

    assert stale == {}


def test_clean_remapped_paired_bam_retains_a_real_advntr_success_contract() -> None:
    """MED11: early-return negative tests cannot replace downstream paired/adVNTR proof."""
    successful = [case for case in get_advntr_test_cases() if case.get("expected_exit_code", 0) == 0]

    assert [case["test_name"] for case in successful] == [
        "example_a5c1_hg19_subset_advntr",
        "example_b178_hg19_bwa_advntr",
    ]
    case = successful[1]
    assert case["bam"] == "tests/data/remapped/bwa/hg19/example_b178_hg19_bwa.bam"
    assert case["expected_fastq_records"] == {"r1": 14689, "r2": 14689, "other": 0, "single": 0}
    assert case["expected_selected_fastqs"] == ["output_R1.fastq.gz", "output_R2.fastq.gz"]
    assert case["advntr_assertions"] == {
        "VID": "25561",
        "State": "I22_4_G_LEN1",
        "NumberOfSupportingReads": 39,
        "MeanCoverage": {"value": 70.3333333333, "tolerance_percentage": 10},
        "Pvalue": {"value": 5.774455097259999e-59, "log10_tolerance": 2},
    }


def test_direct_single_fastq_has_an_end_to_end_integration_contract() -> None:
    """A-161-3 must exercise the single-input fastp/BWA/Kestrel path, not only argument parsing."""
    single_fastq_cases = [case for case in get_fastq_test_cases() if not case.get("fastq2")]

    assert len(single_fastq_cases) == 1
    case = single_fastq_cases[0]
    assert case["test_name"] == "example_6449_hg19_subset_single_fastq"
    assert case["fastq1"] == "tests/data/example_6449_hg19_subset_R1.fastq.gz"
    assert case["expected_fastq_records"] == {"r1": 0, "r2": 0, "other": 40203, "single": 0}
    assert case["expected_selected_fastqs"] == ["output_other.fastq.gz"]
    assert case["kestrel_assertions"]["Depth_Score"] == {"value": 0.15457227138643068, "tolerance_percentage": 0}
    assert case["coverage_assertions"] == {
        "mean": "2255.56",
        "median": "1613.00",
        "stdev": "1575.22",
        "min": "216",
        "max": "6110",
        "region_length": "1501",
        "uncovered_bases": "0",
        "percent_uncovered": "0.00",
        "coverage_qc": "PASS",
    }


def test_paired_shark_fastq_has_the_deterministic_post_dedup_contract() -> None:
    """Pin the exact SHARK output after deduplicating fastp was made deterministic."""
    case = {case["test_name"]: case for case in get_fastq_test_cases()}["example_6449_hg19_subset_fastq_shark"]

    assert case["expected_fastq_records"] == {"r1": 40954, "r2": 40954, "other": 0, "single": 0}
    assert case["coverage_assertions"] == {
        "mean": "3709.29",
        "median": "4138.00",
        "stdev": "1704.10",
        "min": "222",
        "max": "6803",
        "region_length": "1501",
        "uncovered_bases": "0",
        "percent_uncovered": "0.00",
        "coverage_qc": "PASS",
    }


def test_fastq_cases_name_the_exact_quality_control_summary_step() -> None:
    """Catch the display-name shorthand that never appears in pipeline_summary.json."""
    for case in get_fastq_test_cases():
        steps = case["pipeline_summary_assertions"]["steps"]
        assert "FASTQ Quality Control" in steps
        assert "FASTQ QC" not in steps


def test_alternate_paired_fastq_uses_b178_and_omits_shark() -> None:
    case = {case["test_name"]: case for case in get_fastq_test_cases()}[
        "example_b178_hg19_subset_paired_fastq_no_shark"
    ]
    assert case["fastq1"] == "tests/data/example_b178_hg19_subset_R1.fastq.gz"
    assert case["fastq2"] == "tests/data/example_b178_hg19_subset_R2.fastq.gz"
    assert case["reference_assembly"] == "hg19"
    assert case["expected_files"] == ["summary_report.html", "kestrel/kestrel_result.tsv"]
    assert case["coverage_assertions"] == {
        "mean": "566.90",
        "median": "573.00",
        "stdev": "299.18",
        "min": "21",
        "max": "1062",
        "region_length": "1501",
        "uncovered_bases": "0",
        "percent_uncovered": "0.00",
        "coverage_qc": "PASS",
    }
    routing = 'READ_SET_ROUTING {"counts":{"other":0,"r1":11452,"r2":11452,"single":0},"layout":"paired","selected":["output_R1.fastq.gz","output_R2.fastq.gz"]}'
    runner = mock.Mock(return_value=orchestration.PipelineRunResult(0, "", routing))
    with (
        mock.patch.object(orchestration, "validate_strict_fastq_success"),
        mock.patch.object(orchestration, "assert_declared_artifacts"),
    ):
        orchestration.run_fastq_test_case(case, runner, Path("output"))
    request = runner.call_args.args[0]
    assert isinstance(request, orchestration.PipelineRequest)
    assert request.cli_options == ()


def test_local_runner_uses_the_current_worktree_module_not_an_ambient_console_script(tmp_path: Path) -> None:
    """Catch an editable install from another worktree silently executing stale pipeline code."""
    request = orchestration.PipelineRequest(
        input_kind="bam",
        input_paths=(Path("tests/data/input.bam"),),
        reference_assembly="hg19",
        output_dir=tmp_path,
        threads=2,
        log_level="DEBUG",
        cli_options=(),
        reference_fasta=None,
    )
    completed = __import__("subprocess").CompletedProcess([], 0, "stdout", "stderr")

    with mock.patch.object(local_integration, "_run_cli", return_value=completed) as run_cli:
        result = local_integration._run_local_pipeline(request)

    argv = run_cli.call_args.args[0]
    assert argv[:3] == [__import__("sys").executable, "-m", "vntyper.cli"]
    assert argv[3:] == orchestration.build_pipeline_argv(request, str)[1:]
    assert result == orchestration.PipelineRunResult(0, "stdout", "stderr")


def test_single_end_keep_case_names_the_real_unmapped_artifact() -> None:
    cases = {case["test_name"]: case for case in load_test_config()["integration_tests"]["single_end_bam_tests"]}
    keep = cases["example_b178_hg19_single_end_keep"]
    assert keep["cli_options"] == ["--keep-intermediates"]
    assert "fastq_bam_processing/output_unmapped.bam" in keep["expected_files"]
    assert keep["expected_archive"] is False


def test_single_end_cases_declare_the_negative_matrix_after_keep() -> None:
    cases = {case["test_name"]: case for case in load_test_config()["integration_tests"]["single_end_bam_tests"]}
    artifact = "fastq_bam_processing/output_unmapped.bam"
    assert cases["example_b178_hg19_single_end_delete"]["cli_options"] == ["--delete-intermediates"]
    assert cases["example_b178_hg19_single_end_delete"]["expected_absent"] == [artifact]
    assert cases["example_b178_hg19_single_end_delete_overrides_keep"]["cli_options"] == [
        "--keep-intermediates",
        "--delete-intermediates",
    ]
    assert cases["example_b178_hg19_single_end_delete_overrides_keep"]["expected_absent"] == [artifact]
    assert cases["example_b178_hg19_single_end_archive"]["cli_options"] == ["--archive-results"]
    assert cases["example_b178_hg19_single_end_archive"]["expected_archive"] is True


def test_single_end_fast_and_default_cases_pin_their_distinct_measured_read_sets() -> None:
    """Catch copying the default-mode oracle onto the smaller fast-mode region slice."""
    cases = {case["test_name"]: case for case in load_test_config()["integration_tests"]["single_end_bam_tests"]}
    fast = cases["example_b178_hg19_single_end"]
    default = cases["example_b178_hg19_single_end_keep"]

    assert fast["expected_fastq_records"] == {"r1": 0, "r2": 0, "other": 14683, "single": 0}
    assert fast["kestrel_assertions"] == {
        "Estimated_Depth_AlternateVariant": {"value": 414, "tolerance_percentage": 5},
        "Estimated_Depth_Variant_ActiveRegion": {"value": 6111, "tolerance_percentage": 5},
        "Depth_Score": {"value": 0.0677466863033873, "tolerance_percentage": 5},
        "Confidence": "High_Precision*",
    }
    assert default["expected_fastq_records"] == {"r1": 0, "r2": 0, "other": 16922, "single": 0}
    assert default["kestrel_assertions"] == {
        "Estimated_Depth_AlternateVariant": {"value": 416, "tolerance_percentage": 5},
        "Estimated_Depth_Variant_ActiveRegion": {"value": 6168, "tolerance_percentage": 5},
        "Depth_Score": {"value": 0.06744487678339818, "tolerance_percentage": 5},
        "Confidence": "High_Precision*",
    }


def test_mixed_layout_diagnostic_renders_the_exact_dynamic_paths_and_counts(tmp_path: Path) -> None:
    """The contract must match the pipeline's full no-discard diagnostic, not a vague prefix."""
    expected = orchestration.mixed_layout_diagnostic(_negative_case(), tmp_path / "case-output")

    fastq_dir = tmp_path / "case-output" / "fastq_bam_processing"
    assert expected == (
        "FASTQ layout 'mixed' cannot be consumed without dropping reads. Produced FASTQs: "
        f"{fastq_dir / 'output_R1.fastq.gz'}: 7 records, "
        f"{fastq_dir / 'output_R2.fastq.gz'}: 7 records, "
        f"{fastq_dir / 'output_other.fastq.gz'}: 0 records, "
        f"{fastq_dir / 'output_single.fastq.gz'}: 1 records"
    )


def test_mixed_layout_diagnostic_rejects_an_incomplete_count_contract(tmp_path: Path) -> None:
    """No produced FASTQ may disappear from a malformed negative-case declaration."""
    case = _negative_case()
    del case["expected_mixed_fastq_records"]["single"]

    with pytest.raises(KeyError, match="single"):
        orchestration.mixed_layout_diagnostic(case, tmp_path)


def test_declared_artifacts_assert_both_presence_and_absence(tmp_path: Path) -> None:
    kept = tmp_path / "kept.txt"
    kept.write_text("kept\n", encoding="utf-8")
    case = {
        "test_name": "artifact-case",
        "expected_present": ["kept.txt"],
        "expected_absent": ["removed.txt"],
    }
    orchestration.assert_declared_artifacts(case, tmp_path)

    kept.unlink()
    (tmp_path / "removed.txt").write_text("unexpected\n", encoding="utf-8")
    with pytest.raises(AssertionError) as exc_info:
        orchestration.assert_declared_artifacts(case, tmp_path)
    assert f"case=artifact-case field=expected_present missing: {kept.resolve()}" in str(exc_info.value)
    assert (
        f"case=artifact-case field=expected_absent unexpectedly present: {(tmp_path / 'removed.txt').resolve()}"
        in str(exc_info.value)
    )


@pytest.mark.parametrize(
    "case",
    [
        {"expected_present": "same.txt"},
        {"expected_present": [None]},
        {"expected_present": [""]},
        {"expected_present": ["/absolute.txt"]},
        {"expected_present": ["../escape.txt"]},
        {"expected_present": ["nested/../../escape.txt"]},
        {"expected_present": ["same.txt", "same.txt"]},
        {"expected_present": ["same.txt"], "expected_absent": ["same.txt"]},
    ],
)
def test_declared_artifacts_reject_invalid_paths(tmp_path: Path, case: dict) -> None:
    with pytest.raises(ValueError, match="artifact|expected_"):
        orchestration.assert_declared_artifacts(case, tmp_path)


def test_declared_artifacts_reject_normalized_duplicate_paths(tmp_path: Path) -> None:
    case = {"expected_present": ["same.txt", "./same.txt"]}
    with pytest.raises(ValueError, match="invalid artifact"):
        orchestration.assert_declared_artifacts(case, tmp_path)


def test_declared_artifacts_reject_leaf_and_intermediate_symlink_escapes(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "leaf.txt").write_text("outside\n", encoding="utf-8")
    (output / "leaf.txt").symlink_to(outside / "leaf.txt")
    (output / "linked-dir").symlink_to(outside, target_is_directory=True)

    for declared in ("leaf.txt", "linked-dir/leaf.txt"):
        with pytest.raises(ValueError, match="escapes output_dir"):
            orchestration.assert_declared_artifacts({"expected_present": [declared]}, output)


def test_declared_absence_rejects_a_dangling_symlink_entry(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    (output / "broken.txt").symlink_to(tmp_path / "missing.txt")
    case = {"test_name": "broken-link", "expected_absent": ["broken.txt"]}
    with pytest.raises(AssertionError, match="case=broken-link field=expected_absent unexpectedly present"):
        orchestration.assert_declared_artifacts(case, output)


def test_successful_bam_enforces_declared_artifacts(tmp_path: Path) -> None:
    case = dict(get_bam_test_cases()[0])
    case.update(test_name="bam-artifact", expected_present=["must-exist.txt"])
    counts = case["expected_fastq_records"]
    selected = case["expected_selected_fastqs"]
    routing = "READ_SET_ROUTING " + __import__("json").dumps(
        {"counts": counts, "layout": "mixed", "selected": selected},
        sort_keys=True,
        separators=(",", ":"),
    )
    with (
        mock.patch.object(orchestration, "validate_strict_fastq_success"),
        pytest.raises(AssertionError, match="case=bam-artifact field=expected_present"),
    ):
        orchestration.run_bam_test_case(
            case,
            mock.Mock(return_value=orchestration.PipelineRunResult(0, "", routing)),
            tmp_path,
        )


def test_declared_archive_distinguishes_omitted_false_and_true(tmp_path: Path) -> None:
    output = tmp_path / "result"
    output.mkdir()
    archive = Path(f"{output}.zip")
    orchestration.assert_declared_archive({"test_name": "omitted"}, output)
    orchestration.assert_declared_archive({"test_name": "absent", "expected_archive": False}, output)
    with pytest.raises(AssertionError, match="case=present field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "present", "expected_archive": True}, output)
    archive.write_bytes(b"zip")
    orchestration.assert_declared_archive({"test_name": "present", "expected_archive": True}, output)
    with pytest.raises(AssertionError, match="case=absent field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "absent", "expected_archive": False}, output)


def test_declared_archive_rejects_invalid_boolean_and_broken_symlink(tmp_path: Path) -> None:
    output = tmp_path / "result"
    output.mkdir()
    with pytest.raises(ValueError, match="case=invalid field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "invalid", "expected_archive": "false"}, output)
    Path(f"{output}.zip").symlink_to(tmp_path / "missing.zip")
    with pytest.raises(AssertionError, match="case=absent field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "absent", "expected_archive": False}, output)


def test_declared_archive_skips_every_nonzero_outcome(tmp_path: Path) -> None:
    output = tmp_path / "result"
    output.mkdir()
    orchestration.assert_declared_archive(
        {"test_name": "usage", "expected_exit_code": 2, "expected_archive": True},
        output,
    )
    Path(f"{output}.zip").write_bytes(b"zip")
    orchestration.assert_declared_archive(
        {"test_name": "failure", "expected_exit_code": 1, "expected_archive": False},
        output,
    )


def test_bam_orchestration_accepts_the_declared_exit_one_without_success_artifact_checks(tmp_path: Path) -> None:
    """An expected fail-closed BAM case must stop before Kestrel/coverage success assertions."""
    case = _negative_case()
    runner = mock.Mock(return_value=_negative_result())

    with (
        mock.patch.object(orchestration, "validate_strict_fastq_success") as strict_success,
        mock.patch.object(orchestration, "assert_declared_artifacts") as declared_artifacts,
    ):
        orchestration.run_bam_test_case(case, runner, tmp_path)

    request = runner.call_args.args[0]
    assert isinstance(request, orchestration.PipelineRequest)
    assert request.input_paths == (Path("fixture.bam"),)
    assert request.output_dir == tmp_path
    strict_success.assert_not_called()
    declared_artifacts.assert_not_called()


def test_advntr_orchestration_accepts_the_declared_exit_one_without_success_artifact_checks(tmp_path: Path) -> None:
    """The mixed adVNTR fixture has the same early fail-closed contract as ordinary BAM cases."""
    case = _negative_case()
    runner = mock.Mock(return_value=_negative_result())

    with (
        mock.patch.object(orchestration, "assert_required_files") as required_files,
        mock.patch.object(orchestration, "validate_advntr_output") as advntr,
    ):
        orchestration.run_advntr_test_case(case, runner, tmp_path)

    request = runner.call_args.args[0]
    assert isinstance(request, orchestration.PipelineRequest)
    assert request.input_paths == (Path("fixture.bam"),)
    assert request.cli_options == ("--extra-modules", "advntr")
    required_files.assert_not_called()
    advntr.assert_not_called()


def test_clean_bam_cases_still_require_exit_zero(tmp_path: Path) -> None:
    """Omitting the negative contract must retain the existing success requirement."""
    case = {"bam": "clean.bam", "reference_assembly": "hg19"}

    with pytest.raises(AssertionError, match="expected 0, got 1"):
        orchestration.run_bam_test_case(
            case,
            mock.Mock(return_value=orchestration.PipelineRunResult(1, "", "")),
            tmp_path,
        )


@pytest.mark.parametrize(
    ("runner_name", "case"),
    [
        (
            "run_bam_test_case",
            {"test_name": "bam", "bam": "sample.bam", "reference_assembly": "hg19"},
        ),
        (
            "run_fastq_test_case",
            {"test_name": "fastq", "fastq1": "reads.fastq.gz", "reference_assembly": "hg19"},
        ),
        (
            "run_advntr_test_case",
            {"test_name": "advntr", "bam": "sample.bam", "reference_assembly": "hg19"},
        ),
    ],
)
def test_every_success_runner_rejects_a_missing_strict_oracle(
    tmp_path: Path, runner_name: str, case: dict[str, Any]
) -> None:
    """Catch a transport-specific success path bypassing the shared strict schema."""
    runner = getattr(orchestration, runner_name)

    with pytest.raises(ValueError, match="missing strict oracle field"):
        runner(case, lambda _request: orchestration.PipelineRunResult(0, "", ""), tmp_path)


def test_cram_parametrization_and_runner_preserve_the_exact_request(tmp_path: Path) -> None:
    """Catch CRAM declarations being absent or translated outside PipelineRequest."""
    cases = parametrization.get_cram_test_cases()
    assert parametrization.get_cram_test_ids() == [case["test_name"] for case in cases]
    case = {
        "test_name": "cram-case",
        "cram": "tests/data/cram/sample.cram",
        "reference_assembly": "hg19",
        "reference_fasta": "reference/alignment/chr1.hg19.fa",
        "threads": 2,
        "log_level": "DEBUG",
        "cli_options": ["--fast-mode"],
    }
    received: list[orchestration.PipelineRequest] = []

    def runner(request: orchestration.PipelineRequest) -> orchestration.PipelineRunResult:
        received.append(request)
        return orchestration.PipelineRunResult(0, "", "")

    with pytest.raises(ValueError, match="missing strict oracle field"):
        orchestration.run_cram_test_case(case, runner, tmp_path)

    assert received == [
        orchestration.PipelineRequest(
            input_kind="cram",
            input_paths=(Path("tests/data/cram/sample.cram"),),
            reference_assembly="hg19",
            output_dir=tmp_path,
            threads=2,
            log_level="DEBUG",
            cli_options=("--fast-mode",),
            reference_fasta=Path("reference/alignment/chr1.hg19.fa"),
        )
    ]


def test_nonzero_case_requires_an_explicit_causal_diagnostic(tmp_path: Path) -> None:
    """Catch retaining the temporary adapter after valid mixed cases become successes."""
    case = _negative_case()
    case.pop("expected_diagnostic")
    diagnostic = "\n".join(
        [
            "FASTQ layout 'mixed' cannot be consumed without dropping reads.",
            "output_R1.fastq.gz: 7 records",
            "output_R2.fastq.gz: 7 records",
            "output_other.fastq.gz: 0 records",
            "output_single.fastq.gz: 1 records",
        ]
    )

    with pytest.raises(ValueError, match="expected_diagnostic"):
        orchestration.run_bam_test_case(
            case,
            lambda _request: orchestration.PipelineRunResult(1, "", diagnostic),
            tmp_path,
        )


def test_summary_step_list_is_an_exact_stable_oracle(tmp_path: Path) -> None:
    """Catch dropping, reordering, or requiring volatile timestamps in summary validation."""
    summary = {
        "steps": [
            {"step": "BAM Header Parsing", "parsed_result": {"date": "volatile"}},
            {
                "step": "Coverage Calculation",
                "parsed_result": {"comments": [], "data": [{"mean": "878.21"}]},
            },
            {
                "step": "Kestrel Genotyping",
                "parsed_result": {"comments": ["volatile"], "data": [{"mean": "878.21"}]},
            },
        ]
    }
    (tmp_path / "pipeline_summary.json").write_text(__import__("json").dumps(summary), encoding="utf-8")
    case: dict[str, Any] = {
        "pipeline_summary_assertions": {
            "steps": ["BAM Header Parsing", "Coverage Calculation", "Kestrel Genotyping"],
            "parsed_results": ["Coverage Calculation", "Kestrel Genotyping"],
        },
        "coverage_assertions": {"mean": "878.21"},
    }

    with mock.patch.object(
        orchestration,
        "_summary_expected_result",
        return_value={"mean": "878.21"},
        create=True,
    ):
        orchestration._assert_summary_values(case, tmp_path, advntr=False)
    case["pipeline_summary_assertions"]["steps"] = ["Coverage Calculation", "Kestrel Genotyping"]
    with pytest.raises(AssertionError, match="step sequence"):
        orchestration._assert_summary_values(case, tmp_path, advntr=False)


@pytest.mark.parametrize(
    ("step", "expected_error"),
    [
        ("Coverage Calculation", "Coverage Calculation parsed_result differs"),
        ("Kestrel Genotyping", "Kestrel Genotyping parsed_result differs"),
        ("adVNTR Genotyping", "adVNTR Genotyping parsed_result differs"),
        ("Cross-Match Variant Comparison", "Cross-Match Variant Comparison parsed_result differs"),
    ],
)
def test_summary_required_parsed_result_mutations_fail(
    tmp_path: Path,
    step: str,
    expected_error: str,
) -> None:
    """Every result-bearing summary step must agree with its independently checked artifact."""
    required = ["Coverage Calculation", "Kestrel Genotyping"]
    if step in {"adVNTR Genotyping", "Cross-Match Variant Comparison"}:
        required.extend(["adVNTR Genotyping", "Cross-Match Variant Comparison"])
    summary_steps = [
        {
            "step": name,
            "parsed_result": {"comments": [], "data": [{"value": "mutated" if name == step else name}]},
        }
        for name in required
    ]
    (tmp_path / "pipeline_summary.json").write_text(
        __import__("json").dumps({"steps": summary_steps}),
        encoding="utf-8",
    )
    case = {
        "pipeline_summary_assertions": {"steps": required, "parsed_results": required},
        "cross_match_assertions": {
            "comments": [],
            "data": [{"value": "Cross-Match Variant Comparison"}],
        },
    }

    def expected_result(_case: dict, _output: Path, name: str) -> dict[str, str]:
        return {"value": name}

    with (
        mock.patch.object(orchestration, "_summary_expected_result", side_effect=expected_result, create=True),
        pytest.raises(AssertionError, match=expected_error),
    ):
        orchestration._assert_summary_values(case, tmp_path, advntr=len(required) == 4)
