"""Pure checks for integration cases that intentionally fail closed on mixed reads."""

from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests.parametrization import get_advntr_test_cases, get_bam_test_cases, get_fastq_test_cases, load_test_config
from tests.support import orchestration

pytestmark = pytest.mark.unit


EXPECTED_MIXED_COUNTS = {
    "example_b178_hg19_subset_fast": (14_690, 14_690, 0, 1),
    "example_a5c1_hg19_subset_fast": (20_888, 20_888, 0, 3),
    "example_66bf_hg19_subset_fast": (19_841, 19_841, 0, 3),
    "example_7a61_hg19_subset_fast": (2_359, 2_359, 0, 8),
    "example_dfc3_hg19_subset_fast": (31_596, 31_596, 0, 53),
    "example_dfc3_GRCh37_fast": (31_593, 31_593, 0, 2),
    "example_dfc3_GRCh38_fast": (31_603, 31_603, 0, 2),
    "example_dfc3_hg38_ensembl_fast": (31_603, 31_603, 0, 2),
    "example_40cf_hg38_subset_fast_gdp_guard": (3_474, 3_474, 0, 93),
    "example_a5c1_hg19_subset_advntr": (20_888, 20_888, 0, 3),
}


def _negative_case() -> dict[str, Any]:
    return {
        "bam": "fixture.bam",
        "reference_assembly": "hg19",
        "expected_exit_code": 1,
        "expected_mixed_fastq_records": {"r1": 7, "r2": 7, "other": 0, "single": 1},
    }


def test_every_known_mixed_fixture_declares_exit_one_and_its_exact_record_counts() -> None:
    """The ten measured mixed fixtures must be explicit rather than silently treated as successes."""
    configured = {}
    for case in [*get_bam_test_cases(), *get_advntr_test_cases()]:
        if case.get("expected_exit_code") == 1:
            counts = case["expected_mixed_fastq_records"]
            configured[case["test_name"]] = (counts["r1"], counts["r2"], counts["other"], counts["single"])

    assert configured == EXPECTED_MIXED_COUNTS


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

    assert successful == [
        {
            "test_name": "example_b178_hg19_bwa_advntr",
            "bam": "tests/data/remapped/bwa/hg19/example_b178_hg19_bwa.bam",
            "reference_assembly": "hg19",
            "cli_options": [
                "--fast-mode",
                "--keep-intermediates",
                "--extra-modules",
                "advntr",
                "--advntr-max-coverage",
                "300",
            ],
            "advntr_assertions": {
                "VID": "25561",
                "State": "I22_4_G_LEN1",
                "NumberOfSupportingReads": 39,
                "MeanCoverage": {"value": 70.3333333333, "tolerance_percentage": 10},
                "Pvalue": {"value": 5.774455097259999e-59, "log10_tolerance": 2},
            },
        }
    ]


def test_direct_single_fastq_has_an_end_to_end_integration_contract() -> None:
    """A-161-3 must exercise the single-input fastp/BWA/Kestrel path, not only argument parsing."""
    single_fastq_cases = [case for case in get_fastq_test_cases() if not case.get("fastq2")]

    assert single_fastq_cases == [
        {
            "test_name": "example_6449_hg19_subset_single_fastq",
            "fastq1": "tests/data/example_6449_hg19_subset_R1.fastq.gz",
            "reference_assembly": "hg19",
            "expected_files": ["summary_report.html", "kestrel/kestrel_result.tsv"],
        }
    ]


def test_single_end_keep_case_names_the_real_unmapped_artifact() -> None:
    cases = {case["test_name"]: case for case in load_test_config()["integration_tests"]["single_end_bam_tests"]}
    keep = cases["example_b178_hg19_single_end_keep"]
    assert keep["cli_options"] == ["--keep-intermediates"]
    assert keep["expected_present"] == ["fastq_bam_processing/output_unmapped.bam"]
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
    case = {
        "test_name": "bam-artifact",
        "bam": "clean.bam",
        "reference_assembly": "hg19",
        "kestrel_assertions": {},
        "expected_present": ["must-exist.txt"],
    }
    with (
        mock.patch.object(orchestration, "assert_required_files"),
        mock.patch.object(orchestration, "validate_kestrel_output"),
        mock.patch.object(orchestration, "validate_coverage_output", return_value={"mean_cov": 0}),
        pytest.raises(AssertionError, match="case=bam-artifact field=expected_present"),
    ):
        orchestration.run_bam_test_case(case, mock.Mock(return_value=0), tmp_path)


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
    runner = mock.Mock(return_value=1)

    with (
        mock.patch.object(orchestration, "assert_required_files") as required_files,
        mock.patch.object(orchestration, "validate_kestrel_output") as kestrel,
        mock.patch.object(orchestration, "validate_coverage_output") as coverage,
    ):
        orchestration.run_bam_test_case(case, runner, tmp_path)

    runner.assert_called_once_with(Path("fixture.bam"), "hg19", tmp_path)
    required_files.assert_not_called()
    kestrel.assert_not_called()
    coverage.assert_not_called()


def test_advntr_orchestration_accepts_the_declared_exit_one_without_success_artifact_checks(tmp_path: Path) -> None:
    """The mixed adVNTR fixture has the same early fail-closed contract as ordinary BAM cases."""
    case = _negative_case()
    runner = mock.Mock(return_value=1)

    with (
        mock.patch.object(orchestration, "assert_required_files") as required_files,
        mock.patch.object(orchestration, "validate_advntr_output") as advntr,
    ):
        orchestration.run_advntr_test_case(case, runner, tmp_path)

    runner.assert_called_once_with(Path("fixture.bam"), "hg19", tmp_path, ["advntr"], [])
    required_files.assert_not_called()
    advntr.assert_not_called()


def test_clean_bam_cases_still_require_exit_zero(tmp_path: Path) -> None:
    """Omitting the negative contract must retain the existing success requirement."""
    case = {"bam": "clean.bam", "reference_assembly": "hg19"}

    with pytest.raises(AssertionError, match="expected 0, got 1"):
        orchestration.run_bam_test_case(case, mock.Mock(return_value=1), tmp_path)
