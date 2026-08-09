"""Unit tests for the golden-cohort gate's matrix derivation and its drift check.

The check against the documented 50/5/3/3 contract existed before these tests and was
advisory in every direction: ``build_matrix`` logged deviations as warnings, ``cmd_matrix``
returned 0 regardless, and the comparison's verdict ignored them. A silently reduced run
therefore earned the same ``IDENTICAL`` as a full one, and a ``--case`` filter that matched
nothing produced a zero-case matrix that every ``all()`` in the harness then agreed was
verified - because ``all()`` over an empty mapping is True.

These tests also pin what "derived" means. Only the base cases are derived from
``tests/data``; the non-fast ids, the adVNTR ids and the probes are declared policy.
"""

import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from golden_cohort import matrix  # noqa: E402
from golden_cohort.admissibility import COHORT_REQUIRED_ARTIFACTS, PIPELINE_REQUIRED_ARTIFACTS  # noqa: E402


def _data_dir(tmp_path: Path, samples, assemblies) -> Path:
    """Build a fake ``tests/data`` holding one subset BAM per sample per assembly.

    Args:
        tmp_path: pytest's per-test directory.
        samples: Short sample ids, e.g. ``("a5c1",)``.
        assemblies: Assembly names to create a BAM for.

    Returns:
        Path: The data directory.
    """
    root = tmp_path / "data"
    root.mkdir()
    for sample in samples:
        for assembly in assemblies:
            (root / f"example_{sample}_{assembly}_subset.bam").write_bytes(b"")
    return root


def _derive_cram_fixtures(root: Path) -> Path:
    """Mirror what ``make cram-fixtures`` writes: one CRAM beside every BAM.

    ``scripts/make_cram_fixtures.py`` writes ``(fixture_root / relative).with_suffix(
    ".cram")`` for every BAM under the data root, so this reproduces the *layout* it
    produces. The content is irrelevant here - no test in this file opens one - but the
    paths are exactly what ``matrix.cram_fixture_for`` has to compute.

    Args:
        root: The data directory.

    Returns:
        Path: The fixture root.
    """
    cram_root = root / matrix.CRAM_FIXTURE_DIRNAME
    for bam in sorted(root.rglob("*.bam")):
        if cram_root in bam.parents:
            continue
        fixture = (cram_root / bam.relative_to(root)).with_suffix(".cram")
        fixture.parent.mkdir(parents=True, exist_ok=True)
        fixture.write_bytes(b"")
    indexed_safe = cram_root / "indexed-safe" / "indexed-safe.cram"
    indexed_safe.parent.mkdir(parents=True, exist_ok=True)
    indexed_safe.write_bytes(b"")
    return cram_root


def _documented_data_dir(tmp_path: Path, *, with_cram: bool = True) -> Path:
    """Build a fake data directory that reproduces the documented contract.

    Seven multi-reference samples at six assemblies is 42; plus their hg19 subsets already
    counted, plus the hg38-only ``example_40cf``. The shape below reaches the documented
    per-assembly counts (24 hg19, 9 hg38, 8 GRCh38, 9 hg38_ensembl, 7 each of the rest)
    using the two filename shapes ``derive_base_cases`` recognises, plus the six CRAM
    repeats.

    Args:
        tmp_path: pytest's per-test directory.
        with_cram: Whether to derive the CRAM fixtures. ``False`` models a cohort on which
            ``make cram-fixtures`` has never been run.

    Returns:
        Path: The data directory.
    """
    samples = ("6449", "66bf", "6c28", "7a61", "a5c1", "b178", "dfc3")
    root = _data_dir(tmp_path, samples, ("hg19",))
    (root / "example_40cf_hg38_subset.bam").write_bytes(b"")
    remapped = root / "remapped"
    for assembly in ("hg19", "hg38", "GRCh37", "GRCh38", "hg19_ensembl", "hg38_ensembl"):
        target = remapped / "bwa" / assembly
        target.mkdir(parents=True)
        for sample in samples:
            (target / f"example_{sample}_{assembly}.bam").write_bytes(b"")
    # 8 subset + 42 remapped = 50 base cases, and the per-assembly totals land on the
    # documented ones once the 5 non-fast, 3 adVNTR and 6 CRAM repeats are added.
    if with_cram:
        _derive_cram_fixtures(root)
    return root


def _probe_capable_data_dir(tmp_path: Path) -> Path:
    """Build the smallest data directory the three declared probes resolve against.

    ``PROBE_SPECS`` names ``dfc3_hg19_subset``, ``40cf_hg38_subset`` and
    ``dfc3_GRCh38_bwa``, and ``_resolve`` refuses a policy naming a case the data does not
    provide - so a probe fixture has to supply all three, including the remapped one.

    Args:
        tmp_path: pytest's per-test directory.

    Returns:
        Path: The data directory.
    """
    root = _data_dir(tmp_path, ("dfc3",), ("hg19",))
    (root / "example_40cf_hg38_subset.bam").write_bytes(b"")
    remapped = root / "remapped" / "bwa" / "GRCh38"
    remapped.mkdir(parents=True)
    (remapped / "example_dfc3_GRCh38.bam").write_bytes(b"")
    return root


def _build(root: Path, **kwargs):
    """Build a matrix over a fake data directory with the real policies.

    Args:
        root: The data directory.
        **kwargs: Passed through to ``build_matrix``.

    Returns:
        dict: The matrix.
    """
    return matrix.build_matrix(root, **kwargs)


# --------------------------------------------------------------------------------------
# required_artifacts
# --------------------------------------------------------------------------------------


def test_every_derived_base_case_declares_the_pipeline_artefact_requirement(tmp_path: Path) -> None:
    """A case expected to exit zero must say what "exited zero and produced output" means.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(
        _data_dir(tmp_path, ("a5c1", "b178", "dfc3", "40cf"), ("hg19",)),
        non_fast_ids=("a5c1_hg19_subset",),
        advntr_ids=("b178_hg19_subset",),
        cram_ids=(),
        include_probes=False,
        strict=False,
    )
    assert built["cases"], "the fixture must derive at least one case"
    for case in built["cases"]:
        assert case["expect_exit"] == "zero"
        assert case["required_artifacts"] == list(PIPELINE_REQUIRED_ARTIFACTS), case["case_id"]


def test_mixed_base_cases_declare_side_specific_fail_closed_outcomes(tmp_path: Path) -> None:
    """The baseline may succeed only because it drops reads; the candidate must reject them."""
    built = _build(_documented_data_dir(tmp_path))
    by_id = {case["case_id"]: case for case in built["cases"]}
    mixed_base_ids = {
        "40cf_hg38_subset",
        "6449_hg19_subset",
        "66bf_GRCh37_bwa",
        "66bf_GRCh38_bwa",
        "66bf_hg19_bwa",
        "66bf_hg19_ensembl_bwa",
        "66bf_hg19_subset",
        "66bf_hg38_bwa",
        "66bf_hg38_ensembl_bwa",
        "6c28_hg19_subset",
        "7a61_GRCh37_bwa",
        "7a61_GRCh38_bwa",
        "7a61_hg19_bwa",
        "7a61_hg19_ensembl_bwa",
        "7a61_hg19_subset",
        "7a61_hg38_bwa",
        "7a61_hg38_ensembl_bwa",
        "a5c1_GRCh37_bwa",
        "a5c1_GRCh38_bwa",
        "a5c1_hg19_bwa",
        "a5c1_hg19_ensembl_bwa",
        "a5c1_hg19_subset",
        "a5c1_hg38_bwa",
        "a5c1_hg38_ensembl_bwa",
        "b178_hg19_subset",
        "dfc3_GRCh37_bwa",
        "dfc3_GRCh38_bwa",
        "dfc3_hg19_bwa",
        "dfc3_hg19_ensembl_bwa",
        "dfc3_hg19_subset",
        "dfc3_hg38_bwa",
        "dfc3_hg38_ensembl_bwa",
    }

    declared = {case_id for case_id, case in by_id.items() if case.get("side_expectations")}
    base_ids = {case["case_id"] for case in built["cases"] if case["group"] == "base"}
    assert declared.intersection(base_ids) == mixed_base_ids
    for case_id in mixed_base_ids:
        assert by_id[case_id]["side_expectations"]["before"] == {
            "expect_exit": "zero",
            "required_artifacts": list(PIPELINE_REQUIRED_ARTIFACTS),
        }
        assert by_id[case_id]["side_expectations"]["after"] == {
            "expect_exit": "nonzero",
            "required_artifacts": [],
        }


def test_repeats_do_not_blindly_inherit_the_base_layout_outcome(tmp_path: Path) -> None:
    """Non-fast extraction repairs layout, while fast adVNTR repeats retain the refusal."""
    built = _build(_documented_data_dir(tmp_path))
    by_id = {case["case_id"]: case for case in built["cases"]}

    assert by_id["b178_hg19_nonfast"].get("side_expectations") is None
    assert by_id["b178_hg19_nonfast"]["expect_exit"] == "zero"
    assert by_id["a5c1_hg19_advntr"]["side_expectations"]["after"] == {
        "expect_exit": "nonzero",
        "required_artifacts": [],
    }


def test_a_probe_expected_to_fail_requires_no_artefacts(tmp_path: Path) -> None:
    """The two mismatch probes refuse before writing anything; requiring output inverts them.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(_probe_capable_data_dir(tmp_path), non_fast_ids=(), advntr_ids=(), cram_ids=(), strict=False)
    by_id = {probe["case_id"]: probe for probe in built["probes"]}
    assert by_id["probe_mismatch_hg19_as_hg38"]["required_artifacts"] == []
    assert by_id["probe_mismatch_hg38_as_hg19"]["required_artifacts"] == []


def test_the_naming_probe_is_expected_to_succeed_and_carries_the_full_requirement(tmp_path: Path) -> None:
    """It exits 0 on both sides in every recorded run, so it must produce output.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(_probe_capable_data_dir(tmp_path), non_fast_ids=(), advntr_ids=(), cram_ids=(), strict=False)
    naming = next(probe for probe in built["probes"] if probe["case_id"] == "probe_naming_ncbi_as_hg38")
    assert naming["expect_exit"] == "zero"
    assert naming["required_artifacts"] == list(PIPELINE_REQUIRED_ARTIFACTS)
    assert naming.get("side_expectations") is None


def test_cohorts_only_consume_cases_expected_to_write_candidate_summaries(tmp_path: Path) -> None:
    """Intentional candidate refusals must not block every downstream cohort case."""
    built = _build(_documented_data_dir(tmp_path))
    candidate_successes = {
        case["case_id"]
        for case in built["cases"]
        if case.get("side_expectations", {}).get("after", {}).get("expect_exit", case["expect_exit"]) == "zero"
    }

    assert candidate_successes
    assert len(candidate_successes) == 26
    for cohort in built["cohort_cases"]:
        if cohort.get("empty_input_dir"):
            continue
        assert set(cohort["inputs"]).issubset(candidate_successes)
    assert set(built["cohort_cases"][0]["inputs"]) == candidate_successes


def test_exporting_cohorts_receive_a_candidate_successful_advntr_case(tmp_path: Path) -> None:
    """Cohort adVNTR exports are impossible when every upstream adVNTR case refuses."""
    built = _build(_documented_data_dir(tmp_path))
    successful_advntr = [
        case
        for case in built["cases"]
        if case["group"] == "advntr" and matrix.materialize_side_expectation(case, "after")["expect_exit"] == "zero"
    ]

    assert [(case["case_id"], case["repeat_of"]) for case in successful_advntr] == [
        ("b178_hg19_advntr", "b178_hg19_bwa")
    ]
    for cohort in built["cohort_cases"]:
        if "cohort_advntr.tsv" in cohort["required_artifacts"]:
            assert "b178_hg19_advntr" in cohort["inputs"]


def test_the_cohort_cases_declare_their_exports_and_the_empty_one_declares_none() -> None:
    """``cohort_empty`` writes only its log and exits 0 by design, and says so explicitly.

    It is the one legitimate "exited zero, produced no substantive output" case in the
    matrix. Every other cohort case must have written its seven exports.
    """
    cases = {case["case_id"]: case for case in matrix.build_cohort_cases(["a", "b"])}
    assert cases["cohort_multi"]["required_artifacts"] == list(COHORT_REQUIRED_ARTIFACTS)
    assert cases["cohort_single"]["required_artifacts"] == list(COHORT_REQUIRED_ARTIFACTS)
    assert cases["cohort_empty"]["required_artifacts"] == []
    assert "pseudonymization_table.tsv" in cases["cohort_multi_pseudonymized"]["required_artifacts"]


def test_only_the_pseudonymized_cohort_case_requires_the_pseudonymization_table() -> None:
    """It is the only case that writes it, so it is the only case that can require it."""
    cases = {case["case_id"]: case for case in matrix.build_cohort_cases(["a"])}
    requiring = [
        case_id for case_id, case in cases.items() if "pseudonymization_table.tsv" in case["required_artifacts"]
    ]
    assert requiring == ["cohort_multi_pseudonymized"]


# --------------------------------------------------------------------------------------
# drift and emptiness
# --------------------------------------------------------------------------------------


def test_a_matrix_matching_the_documented_contract_builds_and_is_attestation_grade(tmp_path: Path) -> None:
    """The fixture reproduces 50/5/3/6 and the per-assembly counts, so strict mode passes.

    The total is 64 rather than the 58 runs 1-5 measured: the CRAM group is six cases wider.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(_documented_data_dir(tmp_path))
    assert built["check"]["mismatches"] == []
    assert built["check"]["attestation_grade"] is True
    assert built["check"]["counts"]["total"] == 64
    assert built["check"]["counts"]["cram"] == 6
    assert built["check"]["counts"]["probes"] == 3


def test_a_drifted_matrix_is_refused_before_anything_launches(tmp_path: Path) -> None:
    """A run over a matrix that is not the documented one must not reach the pipeline.

    This is the whole of C4. Deriving 4 cases where the page records 58 used to log six
    warnings, return the matrix anyway, and then earn the same ``IDENTICAL`` verdict.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If a drifted matrix is built.
    """
    root = _data_dir(tmp_path, ("a5c1", "b178", "dfc3", "40cf"), ("hg19",))
    with pytest.raises(ValueError, match="deviates from the contract the gate page records"):
        _build(
            root,
            non_fast_ids=("a5c1_hg19_subset",),
            advntr_ids=("b178_hg19_subset",),
            cram_ids=(),
            include_probes=False,
        )


def test_the_refusal_names_every_deviation_so_the_fix_is_obvious(tmp_path: Path) -> None:
    """A refusal that does not say what drifted is a refusal nobody can act on.

    Args:
        tmp_path: pytest's per-test directory.
    """
    root = _data_dir(tmp_path, ("a5c1", "b178", "dfc3", "40cf"), ("hg19",))
    with pytest.raises(ValueError) as excinfo:
        _build(
            root,
            non_fast_ids=("a5c1_hg19_subset",),
            advntr_ids=("b178_hg19_subset",),
            cram_ids=(),
            include_probes=False,
        )
    message = str(excinfo.value)
    assert "base: derived 4, page records 50" in message
    assert "--allow-matrix-drift" in message


def test_a_drifted_matrix_builds_when_drift_is_allowed_but_is_not_attestation_grade(tmp_path: Path) -> None:
    """A deliberately reduced smoke run stays possible, and is marked as what it is.

    Args:
        tmp_path: pytest's per-test directory.
    """
    root = _data_dir(tmp_path, ("a5c1", "b178", "dfc3", "40cf"), ("hg19",))
    built = _build(
        root,
        non_fast_ids=("a5c1_hg19_subset",),
        advntr_ids=("b178_hg19_subset",),
        cram_ids=(),
        include_probes=False,
        strict=False,
    )
    assert built["check"]["mismatches"]
    assert built["check"]["attestation_grade"] is False


def test_a_case_filter_makes_the_matrix_not_attestation_grade_without_refusing_it(tmp_path: Path) -> None:
    """``--case`` is for smoke runs, so it narrows rather than fails - but it is recorded.

    Args:
        tmp_path: pytest's per-test directory.
    """
    built = _build(_documented_data_dir(tmp_path), case_filter=["a5c1"])
    assert built["check"]["skipped"] is True
    assert built["check"]["attestation_grade"] is False
    assert all("a5c1" in case["case_id"] for case in built["cases"])


def test_a_filter_matching_nothing_is_a_hard_failure_even_though_a_filter_is_advisory(tmp_path: Path) -> None:
    """A zero-case matrix is the worst one, and being filtered does not excuse it.

    ``all()`` over an empty mapping is True, so a side that launched nothing reported
    ``launch_verified: true`` and compared clean against another side that also launched
    nothing.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If an empty matrix is returned.
    """
    with pytest.raises(ValueError, match="no pipeline cases"):
        _build(_documented_data_dir(tmp_path), case_filter=["no_such_sample"])


def test_a_zero_case_matrix_is_refused_even_with_drift_allowed(tmp_path: Path) -> None:
    """``--allow-matrix-drift`` permits a smaller matrix, never an absent one.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If ``strict=False`` lets an empty matrix through.
    """
    with pytest.raises(ValueError, match="no pipeline cases"):
        _build(_documented_data_dir(tmp_path), case_filter=["no_such_sample"], strict=False)


def test_an_empty_data_directory_is_refused(tmp_path: Path) -> None:
    """No BAMs at all derives no cases, which is the same failure by another route.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If an empty data directory yields a usable matrix.
    """
    empty = tmp_path / "data"
    empty.mkdir()
    with pytest.raises(ValueError, match="no pipeline cases"):
        _build(empty, non_fast_ids=(), advntr_ids=(), cram_ids=(), include_probes=False, strict=False)


# --------------------------------------------------------------------------------------
# what "derived" means
# --------------------------------------------------------------------------------------


def test_a_policy_naming_a_case_the_data_does_not_provide_is_refused(tmp_path: Path) -> None:
    """Silently dropping it would shrink the matrix without changing a line of the page.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If a stale policy id is dropped instead of refused.
    """
    root = _data_dir(tmp_path, ("a5c1",), ("hg19",))
    with pytest.raises(ValueError, match="non-fast policy names 1 case"):
        _build(root, non_fast_ids=("no_such_case",), advntr_ids=(), cram_ids=(), include_probes=False, strict=False)


def test_the_non_fast_and_advntr_selections_are_policy_and_not_derived(tmp_path: Path) -> None:
    """The matrix records the policy it applied, which is what makes it overridable.

    The gate page describes the whole matrix as "derived from tests/data". Only the base
    cases are; this records the four selections that are not.

    Args:
        tmp_path: pytest's per-test directory.
    """
    root = _probe_capable_data_dir(tmp_path)
    _derive_cram_fixtures(root)
    built = _build(
        root,
        non_fast_ids=("dfc3_hg19_subset",),
        advntr_ids=("40cf_hg38_subset",),
        cram_ids=("dfc3_GRCh38_bwa",),
        include_probes=True,
        strict=False,
    )
    policies = built["policies"]
    assert policies["non_fast_case_ids"] == ["dfc3_hg19_subset"]
    assert policies["advntr_case_ids"] == ["40cf_hg38_subset"]
    assert policies["cram_case_ids"] == ["dfc3_GRCh38_bwa"]
    assert policies["probe_specs"] == [list(spec) for spec in matrix.PROBE_SPECS]
    # The base cases come from the filesystem; the four selections above do not.
    assert {case["case_id"] for case in built["cases"] if case["group"] == "base"} == {
        "dfc3_hg19_subset",
        "40cf_hg38_subset",
        "dfc3_GRCh38_bwa",
    }
