"""Focused unit tests for the golden-cohort CRAM matrix policy.

CRAM declarations are policy rather than filesystem-derived base cases. A missing
derived fixture must therefore be reported as ordinary matrix drift instead of silently
shrinking the attested cohort.
"""

import logging
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from golden_cohort import admissibility, cram_cases, cram_evidence, matrix  # noqa: E402
from golden_cohort.admissibility import PIPELINE_REQUIRED_ARTIFACTS  # noqa: E402

from tests.unit import test_golden_cohort_matrix as matrix_fixtures  # noqa: E402

_data_dir = matrix_fixtures._data_dir
_derive_cram_fixtures = matrix_fixtures._derive_cram_fixtures
_documented_data_dir = matrix_fixtures._documented_data_dir
_build = matrix_fixtures._build


# --------------------------------------------------------------------------------------
# the CRAM group (#188)
# --------------------------------------------------------------------------------------


def _cram_cases(built) -> dict:
    """Index a built matrix's CRAM cases by id.

    Args:
        built: The matrix.

    Returns:
        dict: The CRAM cases, keyed on case id.
    """
    return {case["case_id"]: case for case in built["cases"] if case["group"] == "cram"}


def test_a_cram_case_is_built_for_each_declared_id_when_the_fixture_exists(tmp_path: Path) -> None:
    """The whole point of #188: the gate finally feeds VNtyper a CRAM.

    Args:
        tmp_path: pytest's per-test directory.
    """
    root = _documented_data_dir(tmp_path)
    built = _build(root)
    cram = _cram_cases(built)
    assert sorted(cram) == [
        "7a61_hg38_ensembl_indexed_cram",
        "7a61_hg38_ensembl_stream_cram",
        "b178_hg19_indexed_cram",
        "b178_hg19_stream_cram",
        "indexed_safe_indexed_cram",
        "indexed_safe_stream_cram",
    ]
    case = cram["b178_hg19_indexed_cram"]
    assert case["kind"] == "pipeline"
    assert case["alignment_kind"] == "cram"
    assert case["expect_exit"] == "nonzero"
    assert case["required_artifacts"] == []
    assert case["expected_stderr_contains"] == "idxstats reports 329 placed-unmapped reads; using stream scan"
    assert case["repeat_of"] == "b178_hg19_subset"


def test_cram_cases_pin_forced_index_refusal_and_stream_success_evidence(tmp_path: Path) -> None:
    """Unsafe forced-index cases fail causally while their lossless stream twins succeed."""
    cram = _cram_cases(_build(_documented_data_dir(tmp_path)))
    indexed = cram["7a61_hg38_ensembl_indexed_cram"]
    stream = cram["7a61_hg38_ensembl_stream_cram"]
    expected_raw = {
        "count": 2_690,
        "sorted_read_name_sha256": "c64be739cd6344b8b62004fc9ea568779b3cc06ff1d472ac0e5d97c343130d7d",
    }
    expected_stream = {
        "count": 634_261,
        "sorted_read_name_sha256": "b7f75d19497698f12d6dbbc38afc12702b2d262670a4c893b39f95967ebf7b7b",
    }

    assert indexed.get("side_expectations") is None
    assert indexed["expect_exit"] == "nonzero"
    assert indexed["required_artifacts"] == []
    assert indexed["expected_stderr_contains"] == "idxstats reports 11571 placed-unmapped reads; using stream scan"
    assert stream.get("side_expectations") is None
    assert stream["expect_exit"] == "zero"
    assert stream["required_artifacts"] == list(PIPELINE_REQUIRED_ARTIFACTS)
    assert "expected_stderr_contains" not in stream
    for case in (indexed, stream):
        assert case["cram_evidence_expectation"] == {
            "placed_unmapped_guard_count": 11_571,
            "raw_indexed_read_set": expected_raw,
            "stream_read_set": expected_stream,
        }

    b178_expected_raw = {
        "count": 4_478,
        "sorted_read_name_sha256": "dad9a81a4e8cf30d1d938717459614f7d8ac6decb84978a5bc23c090b4d90a8b",
    }
    b178_expected_stream = {
        "count": 4_807,
        "sorted_read_name_sha256": "d3aa88fe91c8964b2f9a1b053a672f2bc3d1896b71de986f5cde02999d552591",
    }
    b178_indexed = cram["b178_hg19_indexed_cram"]
    b178_stream = cram["b178_hg19_stream_cram"]
    assert b178_indexed["expect_exit"] == "nonzero"
    assert b178_indexed["required_artifacts"] == []
    assert b178_indexed["expected_stderr_contains"] == ("idxstats reports 329 placed-unmapped reads; using stream scan")
    assert b178_stream["expect_exit"] == "zero"
    assert b178_stream["required_artifacts"] == list(PIPELINE_REQUIRED_ARTIFACTS)
    assert "expected_stderr_contains" not in b178_stream
    for case in (b178_indexed, b178_stream):
        assert case.get("side_expectations") is None
        assert case["cram_evidence_expectation"] == {
            "placed_unmapped_guard_count": 329,
            "raw_indexed_read_set": b178_expected_raw,
            "stream_read_set": b178_expected_stream,
        }


@pytest.mark.parametrize(
    "case_id",
    ["b178_hg19_indexed_cram", "7a61_hg38_ensembl_indexed_cram"],
)
def test_forced_indexed_outcome_agrees_with_its_read_set_evidence(tmp_path: Path, case_id: str) -> None:
    """The ordinary outcome oracle and A-178-2 evidence must accept the same guard refusal."""
    case = _cram_cases(_build(_documented_data_dir(tmp_path)))[case_id]
    expectation = case["cram_evidence_expectation"]
    guard_count = expectation["placed_unmapped_guard_count"]
    record = {
        "exit_code": 1,
        "aborted": False,
        "timed_out": False,
        "raw_indexed_read_set": expectation["raw_indexed_read_set"],
        "unmapped_read_set": None,
        "raw_indexed_loss": None,
        "placed_unmapped_guard_count": guard_count,
        "observed_unmapped_scan": None,
        "observed_unmapped_command": None,
        "scan_observation_problems": ["A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"],
    }

    ordinary = admissibility.check_case(
        case,
        record,
        tmp_path,
        stderr=f"idxstats reports {guard_count} placed-unmapped reads; using stream scan",
    )
    evidence_case = {**case, "effective_unmapped_scan": "indexed"}

    assert ordinary["expectation_met"] is True
    assert ordinary["expectation_problems"] == []
    assert cram_evidence.validate_cram_evidence(evidence_case, record) == []


def test_indexed_safe_cram_cases_require_nonempty_authorized_scan_equivalence(tmp_path: Path) -> None:
    """HIGH2: A-178-2 needs one genuine indexed execution, not only forced rejections."""
    cram = _cram_cases(_build(_documented_data_dir(tmp_path)))
    expected = {
        "count": 20,
        "sorted_read_name_sha256": "16a0efa7785630c3d80716d9a386ddaa24f4933b5671f4ecd221b42a8dffe740",
    }

    for scan in ("indexed", "stream"):
        case = cram[f"indexed_safe_{scan}_cram"]
        assert case["expect_exit"] == "zero"
        assert case["required_artifacts"] == [
            artifact for artifact in PIPELINE_REQUIRED_ARTIFACTS if artifact != "kestrel/kestrel_pre_result.tsv"
        ]
        assert case["cram_evidence_expectation"] == {
            "indexed_authorized": True,
            "raw_indexed_read_set": expected,
            "stream_read_set": expected,
        }


def test_a_cram_case_never_runs_in_fast_mode(tmp_path: Path) -> None:
    """``--fast-mode`` skips the unmapped-read extraction, which is the code under test.

    ``fastq_bam_processing.process_bam_to_fastq`` reaches
    ``build_cram_unmapped_filter_command`` only inside ``if not fast_mode:``. A fast-mode
    CRAM case would exercise the slice and the FASTQ conversion and none of the branch the
    fixtures exist for, so it would be theatre rather than coverage.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(_documented_data_dir(tmp_path))
    cram = _cram_cases(built)
    assert cram, "the fixture must build the CRAM cases"
    for case_id, case in cram.items():
        assert case["fast_mode"] is False, case_id


def test_cram_cases_exercise_both_lossless_scan_strategies(tmp_path: Path) -> None:
    """Every declared CRAM fixture is run once through indexed and stream extraction."""
    cram = _cram_cases(_build(_documented_data_dir(tmp_path)))

    assert len(cram) == 6
    assert {case["unmapped_scan"] for case in cram.values()} == {"indexed", "stream"}
    assert {case["repeat_of"] for case in cram.values()} == {
        "b178_hg19_subset",
        "7a61_hg38_ensembl_bwa",
        "indexed_safe",
    }


def test_a_cram_case_points_at_the_derived_fixture_and_not_at_the_bam(tmp_path: Path) -> None:
    """The fixture path mirrors the source layout, on both shapes the cohort uses.

    ``b178_hg19_subset`` is a top-level subset BAM and ``7a61_hg38_ensembl_bwa`` sits under
    ``remapped/<aligner>/<assembly>/``; the mirrored path has to be right for both. The
    ``bam`` key is dropped so that nothing downstream can read it and quietly launch the BAM.

    Args:
        tmp_path: pytest's per-test directory.
    """
    root = _documented_data_dir(tmp_path)
    cram = _cram_cases(_build(root))

    subset = cram["b178_hg19_indexed_cram"]
    assert Path(subset["cram"]) == root / "cram" / "example_b178_hg19_subset.cram"
    assert Path(subset["source_bam"]) == root / "example_b178_hg19_subset.bam"

    remapped = cram["7a61_hg38_ensembl_indexed_cram"]
    expected = root / "cram" / "remapped" / "bwa" / "hg38_ensembl" / "example_7a61_hg38_ensembl.cram"
    assert Path(remapped["cram"]) == expected

    for case_id, case in cram.items():
        assert "bam" not in case, case_id


def test_a_missing_cram_fixture_is_logged_and_leaves_the_group_short(tmp_path: Path) -> None:
    """A cohort on which ``make cram-fixtures`` never ran must say so, not run 58 cases.

    Args:
        tmp_path: pytest's per-test directory.
    """
    root = _documented_data_dir(tmp_path, with_cram=False)
    built = _build(root, strict=False)
    assert _cram_cases(built) == {}
    assert built["check"]["counts"]["cram"] == 0
    assert "cram: derived 0, page records 6" in built["check"]["mismatches"]
    assert built["check"]["attestation_grade"] is False
    assert any("CRAM fixtures missing for" in line for line in built["derivation_log"])


def test_a_missing_cram_fixture_is_reported_at_error_level(tmp_path: Path, caplog) -> None:
    """Being loud means an error a reader sees, not only a line in the derivation log.

    Args:
        tmp_path: pytest's per-test directory.
        caplog: pytest's log capture.
    """
    root = _documented_data_dir(tmp_path, with_cram=False)
    with caplog.at_level(logging.ERROR, logger="golden_cohort.matrix"):
        _build(root, strict=False)
    messages = [record.getMessage() for record in caplog.records if record.levelno >= logging.ERROR]
    assert any("no CRAM fixture for b178_hg19_subset" in message for message in messages), messages
    assert any("make cram-fixtures" in message for message in messages), messages


def test_a_missing_cram_fixture_refuses_a_strict_build(tmp_path: Path) -> None:
    """A short CRAM group is ordinary drift, so it is refused like any other drift.

    There is deliberately no "0 or 2 is fine" rule. A run without the CRAM cases covers
    strictly less than the contract records and must not earn an attestation-grade verdict;
    ``--allow-matrix-drift`` remains the knowing way to run it anyway.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If a matrix missing its CRAM cases builds in strict mode.
    """
    root = _documented_data_dir(tmp_path, with_cram=False)
    with pytest.raises(ValueError) as excinfo:
        _build(root)
    message = str(excinfo.value)
    assert "cram: derived 0, page records 6" in message
    assert "--allow-matrix-drift" in message


def test_check_matrix_counts_the_cram_group(tmp_path: Path) -> None:
    """``check_matrix`` has to count the group, or the drift check cannot see it shrink.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(_documented_data_dir(tmp_path))
    recount = matrix.check_matrix(built["cases"], built["probes"])
    assert recount["counts"]["cram"] == 6
    assert recount["counts"]["total"] == 64
    assert recount["documented"]["cram"] == matrix.DOCUMENTED_TOTALS["cram"]
    assert recount["mismatches"] == []

    without_cram = [case for case in built["cases"] if case["group"] != "cram"]
    reduced = matrix.check_matrix(without_cram, built["probes"])
    assert reduced["counts"]["cram"] == 0
    assert "cram: derived 0, page records 6" in reduced["mismatches"]
    assert reduced["attestation_grade"] is False


def test_the_cram_cases_move_the_documented_per_assembly_counts(tmp_path: Path) -> None:
    """Adding a case adds it to its assembly's count too, so the two must move together.

    ``b178_hg19_subset`` plus ``indexed_safe`` take hg19 from the 20 runs 1-5 measured to 24, and
    ``7a61_hg38_ensembl_bwa`` takes hg38_ensembl from 7 to 9. Leaving
    ``DOCUMENTED_ASSEMBLY_COUNTS`` alone would make every full run drift on two assemblies.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(_documented_data_dir(tmp_path))
    assert built["check"]["counts"]["by_assembly"]["hg19"] == 24
    assert built["check"]["counts"]["by_assembly"]["hg38_ensembl"] == 9
    assert matrix.DOCUMENTED_ASSEMBLY_COUNTS["hg19"] == 24
    assert matrix.DOCUMENTED_ASSEMBLY_COUNTS["hg38_ensembl"] == 9


def test_a_cram_policy_naming_a_case_the_data_does_not_provide_is_refused(tmp_path: Path) -> None:
    """A stale CRAM id is the same failure as a stale non-fast id, and fails the same way.

    This is distinct from a missing *fixture*: the base case itself is gone, so there is
    nothing to derive a fixture path from and no honest smaller matrix to run.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If a stale CRAM id is dropped instead of refused.
    """
    root = _data_dir(tmp_path, ("a5c1",), ("hg19",))
    _derive_cram_fixtures(root)
    with pytest.raises(ValueError, match="CRAM policy names 1 case"):
        _build(root, non_fast_ids=(), advntr_ids=(), cram_ids=("no_such_case",), include_probes=False, strict=False)


def test_a_fixture_path_cannot_be_derived_for_a_bam_outside_the_data_directory(tmp_path: Path) -> None:
    """The fixture tree mirrors the data tree, so a BAM outside it has no mirrored position.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If an unmirrorable BAM silently produces a path.
    """
    case = {"case_id": "stray", "bam": str(tmp_path / "elsewhere" / "example_a5c1_hg19_subset.bam")}
    with pytest.raises(ValueError, match="is not under the data directory"):
        cram_cases.cram_fixture_for(case, tmp_path / "data", tmp_path / "data" / "cram")
