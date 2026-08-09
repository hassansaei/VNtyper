"""Unit tests for the golden-cohort gate's admissibility checks.

The gate harness had no unit test of any kind, and ``scripts/`` is not in
``[tool.coverage.run] source``, so none of the ~2600 lines backing every genotype claim on
this project was measured by anything. These tests cover the three checks that decide
whether a gate run is evidence at all:

* a case did what its matrix entry declared (exit code and required artefacts);
* the two sides are genuinely opposed, rather than one tree compared against itself;
* the revision each side ran is recorded rather than asserted.

The first is the one that matters most. Before it existed, two runs that both died at
exit 1 before writing a single genotype artefact compared *equal* and the gate exited 0 -
see ``test_two_sides_that_both_produced_nothing_are_not_identical`` in
``test_golden_cohort_compare.py`` for that failure end to end.
"""

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from golden_cohort import admissibility  # noqa: E402


def _record(exit_code, *, aborted=False, timed_out=False):
    """Build a minimal runner record.

    Args:
        exit_code: The recorded exit code, or None.
        aborted: Whether the launcher refused to start.
        timed_out: Whether the runner killed the case.

    Returns:
        dict: The record shape ``check_case`` reads.
    """
    return {"exit_code": exit_code, "aborted": aborted, "timed_out": timed_out}


# --------------------------------------------------------------------------------------
# check_exit
# --------------------------------------------------------------------------------------


def test_a_zero_expected_case_that_exits_zero_is_admissible() -> None:
    """The ordinary success path reports no problem."""
    assert admissibility.check_exit("zero", 0, aborted=False, timed_out=False) == ""


def test_a_zero_expected_case_that_exits_one_is_reported() -> None:
    """A case declared to succeed and did not is named, with both codes."""
    problem = admissibility.check_exit("zero", 1, aborted=False, timed_out=False)
    assert problem == "expected exit 0, got 1"


def test_a_nonzero_expected_probe_that_exits_one_is_admissible() -> None:
    """The two deliberate-mismatch probes are supposed to fail, and exit 1 when they do."""
    assert admissibility.check_exit("nonzero", 1, aborted=False, timed_out=False) == ""


def test_a_nonzero_expected_probe_that_succeeds_is_reported() -> None:
    """A probe that stops refusing is a finding, not a pass."""
    problem = admissibility.check_exit("nonzero", 0, aborted=False, timed_out=False)
    assert problem == "expected a nonzero exit, got 0"


@pytest.mark.parametrize("code", [97, 98, 99, 2])
def test_an_exit_code_vntyper_never_produces_never_satisfies_an_expectation(code: int) -> None:
    """97/98/99 are the harness's own codes and 2 is nobody's; none is a pipeline result.

    A probe that is supposed to fail must fail the way the product fails. Accepting the
    launcher's own refusal (97) as "nonzero, as expected" would let an abort masquerade as
    a successful probe.

    Args:
        code: An exit code outside ``vntyper``'s own {0, 1}.
    """
    problem = admissibility.check_exit("nonzero", code, aborted=False, timed_out=False)
    assert "not one of vntyper's own exit codes" in problem


def test_an_aborted_run_fails_its_expectation_whichever_way_it_points() -> None:
    """The launcher refusing to start is not a pipeline outcome."""
    problem = admissibility.check_exit("nonzero", 97, aborted=True, timed_out=False)
    assert "the launcher refused to start" in problem


def test_a_timed_out_run_fails_its_expectation() -> None:
    """A killed run's exit code says nothing about what the code decided."""
    problem = admissibility.check_exit("zero", 99, aborted=False, timed_out=True)
    assert "killed on the harness timeout" in problem


def test_a_case_with_no_recorded_exit_code_fails_its_expectation() -> None:
    """None means the case never ran, which no expectation can be met by."""
    problem = admissibility.check_exit("zero", None, aborted=False, timed_out=False)
    assert problem == "no exit code was recorded, so this case never ran"


def test_an_unknown_expect_exit_raises_rather_than_exempting_the_case() -> None:
    """A typo in ``expect_exit`` must not silently skip the whole check.

    Raises:
        AssertionError: If a mistyped expectation is accepted.
    """
    with pytest.raises(ValueError, match="Unknown expect_exit"):
        admissibility.check_exit("nonzeo", 1, aborted=False, timed_out=False)


# --------------------------------------------------------------------------------------
# required artefacts
# --------------------------------------------------------------------------------------


def test_a_case_that_wrote_every_required_artefact_is_admissible(tmp_path: Path) -> None:
    """All five declared files present means nothing is reported.

    Args:
        tmp_path: pytest's per-test directory.
    """
    for name in admissibility.PIPELINE_REQUIRED_ARTIFACTS:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("x", encoding="utf-8")
    assert admissibility.missing_artifacts(tmp_path, admissibility.PIPELINE_REQUIRED_ARTIFACTS) == []


def test_a_case_that_exited_zero_and_wrote_nothing_is_reported(tmp_path: Path) -> None:
    """Exit 0 with an empty output directory is the silent-nothing failure.

    Args:
        tmp_path: pytest's per-test directory.
    """
    case = {
        "case_id": "x",
        "expect_exit": "zero",
        "required_artifacts": list(admissibility.PIPELINE_REQUIRED_ARTIFACTS),
    }
    detail = admissibility.check_case(case, _record(0), tmp_path)
    assert detail["expectation_met"] is False
    assert detail["missing_artifacts"] == list(admissibility.PIPELINE_REQUIRED_ARTIFACTS)
    assert "did not write 5 of 5 required artefact(s)" in detail["expectation_problems"][0]


def test_a_directory_named_like_a_required_artefact_does_not_count(tmp_path: Path) -> None:
    """The check is ``is_file()``, so an empty directory of the right name is still missing.

    Args:
        tmp_path: pytest's per-test directory.
    """
    (tmp_path / "pipeline_summary.json").mkdir()
    assert "pipeline_summary.json" in admissibility.missing_artifacts(tmp_path, ("pipeline_summary.json",))


def test_a_nonzero_probe_requires_no_artefacts(tmp_path: Path) -> None:
    """The mismatch probes refuse before writing anything, which is the point of them.

    Args:
        tmp_path: pytest's per-test directory.
    """
    case = {"case_id": "probe_mismatch_hg19_as_hg38", "expect_exit": "nonzero", "required_artifacts": []}
    assert admissibility.check_case(case, _record(1), tmp_path)["expectation_met"] is True


@pytest.mark.parametrize("stderr", ["", "reference could not be decoded", "FASTQ layout 'single' was rejected"])
def test_a_declared_mixed_layout_refusal_rejects_an_unrelated_exit_one(tmp_path: Path, stderr: str) -> None:
    """Missing tools, references, or another layout cannot satisfy the measured refusal."""
    expected = "FASTQ layout 'mixed' cannot be consumed without dropping reads."
    case = {
        "case_id": "measured-mixed",
        "expect_exit": "nonzero",
        "required_artifacts": [],
        "expected_stderr_contains": expected,
    }

    detail = admissibility.check_case(case, _record(1), tmp_path, stderr=stderr)

    assert detail["expectation_met"] is False
    assert detail["stderr_problem"] == f"expected stderr to contain {expected!r}"
    assert detail["expectation_problems"] == [detail["stderr_problem"]]


def test_a_declared_mixed_layout_refusal_accepts_the_causal_diagnostic(tmp_path: Path) -> None:
    """The stable prefix proves the exit came from fail-closed read routing."""
    expected = "FASTQ layout 'mixed' cannot be consumed without dropping reads."
    case = {
        "case_id": "measured-mixed",
        "expect_exit": "nonzero",
        "required_artifacts": [],
        "expected_stderr_contains": expected,
    }

    detail = admissibility.check_case(
        case,
        _record(1),
        tmp_path,
        stderr=f"ERROR: {expected} Produced FASTQs: r1=1, r2=1, single=1",
    )

    assert detail["expectation_met"] is True
    assert detail["stderr_problem"] == ""


@pytest.mark.parametrize("invalid", ["", "   ", [], ["mixed"], False, True])
def test_a_legacy_case_cannot_bypass_causal_diagnostic_validation(tmp_path: Path, invalid: object) -> None:
    """Top-level expectations must be as strict as materialized side expectations."""
    case = {
        "case_id": "legacy-mixed",
        "expect_exit": "nonzero",
        "required_artifacts": [],
        "expected_stderr_contains": invalid,
    }

    with pytest.raises(ValueError) as raised:
        admissibility.check_case(case, _record(1), tmp_path, stderr="an unrelated exit")

    assert str(raised.value) == "Case legacy-mixed has invalid expected_stderr_contains"


def test_artefacts_are_not_checked_when_the_exit_code_already_failed(tmp_path: Path) -> None:
    """One failure, one message. A crashed run's missing files are a consequence, not a finding.

    Args:
        tmp_path: pytest's per-test directory.
    """
    case = {"case_id": "x", "expect_exit": "zero", "required_artifacts": ["pipeline_summary.json"]}
    detail = admissibility.check_case(case, _record(1), tmp_path)
    assert detail["missing_artifacts"] == []
    assert detail["expectation_problems"] == ["expected exit 0, got 1"]


def test_the_pipeline_requirement_does_not_include_advntr() -> None:
    """adVNTR is optional, so requiring its output would fail 56 legitimate cases.

    The consequence is named in ``admissibility.PIPELINE_REQUIRED_ARTIFACTS``'s own note:
    an adVNTR case that silently writes no adVNTR output is caught by the row-set
    comparison, not by this rule.
    """
    assert not any("advntr" in name.lower() for name in admissibility.PIPELINE_REQUIRED_ARTIFACTS)


# --------------------------------------------------------------------------------------
# side opposition
# --------------------------------------------------------------------------------------


def _side(label, tree, expect_marker, *, marker="vntyper.scripts.cohort_rules", head=None):
    """Build a minimal side record.

    Args:
        label: ``before`` or ``after``.
        tree: The recorded source tree.
        expect_marker: ``present`` or ``absent``.
        marker: The marker module name.
        head: A recorded commit, or None for a side with no revision.

    Returns:
        dict: The side record shape ``check_sides_opposed`` reads.
    """
    record = {
        "side": label,
        "tree": tree,
        "marker": marker,
        "expect_marker": expect_marker,
        # One result, because a side with none is refused on that ground alone and would
        # make every other assertion in this section pass for the wrong reason.
        "pipeline_results": {"a5c1_hg19_subset": {"exit_code": 0}},
        "cohort_results": {},
    }
    if head is not None:
        record["revision"] = {"head": head}
    return record


def test_two_properly_opposed_sides_raise_no_problem() -> None:
    """The configuration every recorded run used is accepted."""
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/base", "absent", head="a" * 40),
        _side("after", "/trees/cand", "present", head="b" * 40),
    )
    assert problems == []


def test_the_same_run_root_for_both_sides_is_refused() -> None:
    """Comparing a directory with itself yields IDENTICAL and means nothing."""
    problems = admissibility.check_sides_opposed(
        Path("/gate/one"),
        Path("/gate/one"),
        _side("before", "/trees/base", "absent"),
        _side("after", "/trees/cand", "present"),
    )
    assert any("same directory" in problem for problem in problems)


def test_the_same_source_tree_on_both_sides_is_refused() -> None:
    """This is the failure the whole resolution wrapper exists to prevent.

    The launcher proves each process ran the tree its caller named. It cannot notice that
    the caller named the same tree twice.
    """
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/same", "absent"),
        _side("after", "/trees/same", "present"),
    )
    assert any("same source tree" in problem for problem in problems)


def test_the_same_commit_on_both_sides_is_refused() -> None:
    """Two different paths can be checkouts of one commit; the revision catches that."""
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/base", "absent", head="c" * 40),
        _side("after", "/trees/cand", "present", head="c" * 40),
    )
    assert any("same commit" in problem for problem in problems)


def test_a_side_labelled_before_in_the_after_slot_is_refused() -> None:
    """``--after-root`` pointing at a baseline run inverts every attribution."""
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/base", "absent"),
        _side("before", "/trees/cand", "present"),
    )
    assert any("labelled 'before', not 'after'" in problem for problem in problems)


def test_two_sides_expecting_the_same_marker_state_are_refused() -> None:
    """The marker distinguishes the trees only when the two expectations are opposed."""
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/base", "present"),
        _side("after", "/trees/cand", "present"),
    )
    assert any("both sides expected the marker" in problem for problem in problems)


def test_two_sides_using_different_marker_modules_are_refused() -> None:
    """Neither side's resolution proof says anything about the other's tree."""
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/base", "absent", marker="vntyper.scripts.pipeline_guards"),
        _side("after", "/trees/cand", "present", marker="vntyper.scripts.cohort_rules"),
    )
    assert any("different marker modules" in problem for problem in problems)


def test_a_side_that_recorded_no_case_results_is_refused() -> None:
    """Two sides that measured nothing compare clean, for lack of anything to differ.

    ``all()`` over an empty mapping is True, so a zero-case side used to report itself
    launch-verified as well.
    """
    empty = _side("after", "/trees/cand", "present")
    empty["pipeline_results"] = {}
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/base", "absent"),
        empty,
    )
    assert problems == ["after: this side has no case results at all, so it measured nothing"]


def test_two_sides_with_no_recorded_revision_are_not_refused_for_that_alone() -> None:
    """A run root written by the 1.0.0 harness has no revision, and must still compare.

    Refusing it would make every existing run root unreadable. The absence is warned about
    at the call site instead, and ``--expect-*-sha`` turns it into an error when the caller
    cares.
    """
    problems = admissibility.check_sides_opposed(
        Path("/gate/before"),
        Path("/gate/after"),
        _side("before", "/trees/base", "absent"),
        _side("after", "/trees/cand", "present"),
    )
    assert problems == []


# --------------------------------------------------------------------------------------
# revision recording
# --------------------------------------------------------------------------------------


def _fake_git(responses):
    """Build a subprocess runner that replays canned git output.

    Args:
        responses: Mapping of the git subcommand tuple to ``(returncode, stdout, stderr)``.

    Returns:
        Callable: A stand-in for ``subprocess.run``.
    """

    def run(command, **_kwargs):
        key = tuple(command[3:])
        code, out, err = responses.get(key, (1, "", f"unexpected git call: {key}"))
        return subprocess.CompletedProcess(command, code, out, err)

    return run


def test_a_clean_tree_records_its_head_branch_and_cleanliness() -> None:
    """The happy path records everything an attestation needs to name a commit."""
    run = _fake_git(
        {
            ("rev-parse", "HEAD"): (0, "9b5486e" + "0" * 33 + "\n", ""),
            ("rev-parse", "--abbrev-ref", "HEAD"): (0, "fix/issue-181-197-followups\n", ""),
            ("status", "--porcelain"): (0, "", ""),
            ("status", "--porcelain", "--", "vntyper", "docker", "scripts"): (0, "", ""),
        }
    )
    info = admissibility.describe_tree(Path("/trees/cand"), run=run)
    assert info["head"].startswith("9b5486e")
    assert info["branch"] == "fix/issue-181-197-followups"
    assert info["dirty"] is False
    assert info["dirty_relevant"] is False
    assert info["error"] is None


def test_a_tree_dirty_only_outside_the_genotype_paths_is_not_relevant_dirty() -> None:
    """A modified README does not change a genotype, and must not fail a run.

    The scoped ``git status -- vntyper docker scripts`` is what ``--require-clean`` reads;
    the whole-tree state is recorded beside it so the distinction stays visible.
    """
    run = _fake_git(
        {
            ("rev-parse", "HEAD"): (0, "a" * 40 + "\n", ""),
            ("rev-parse", "--abbrev-ref", "HEAD"): (0, "main\n", ""),
            ("status", "--porcelain"): (0, " M README.md\n", ""),
            ("status", "--porcelain", "--", "vntyper", "docker", "scripts"): (0, "", ""),
        }
    )
    info = admissibility.describe_tree(Path("/trees/cand"), run=run)
    assert info["dirty"] is True
    assert info["dirty_relevant"] is False
    assert info["dirty_paths"] == []


def test_a_tree_dirty_inside_the_genotype_paths_records_the_lines() -> None:
    """An edited stage module is exactly what must not be attested as a commit."""
    run = _fake_git(
        {
            ("rev-parse", "HEAD"): (0, "a" * 40 + "\n", ""),
            ("rev-parse", "--abbrev-ref", "HEAD"): (0, "main\n", ""),
            ("status", "--porcelain"): (0, " M vntyper/scripts/scoring.py\n", ""),
            ("status", "--porcelain", "--", "vntyper", "docker", "scripts"): (
                0,
                " M vntyper/scripts/scoring.py\n",
                "",
            ),
        }
    )
    info = admissibility.describe_tree(Path("/trees/cand"), run=run)
    assert info["dirty_relevant"] is True
    assert info["dirty_paths"] == [" M vntyper/scripts/scoring.py"]


def test_a_tree_that_is_not_a_git_checkout_records_the_error_and_does_not_raise() -> None:
    """Losing a completed 90-minute side to a missing .git would be the worse failure."""
    run = _fake_git({})
    info = admissibility.describe_tree(Path("/trees/nope"), run=run)
    assert info["head"] is None
    assert info["error"]


def test_a_missing_git_binary_records_the_error_and_does_not_raise() -> None:
    """``subprocess.run`` raising OSError is caught, not propagated."""

    def run(_command, **_kwargs):
        raise FileNotFoundError("no git here")

    info = admissibility.describe_tree(Path("/trees/cand"), run=run)
    assert info["head"] is None
    assert "FileNotFoundError" in info["error"]


def test_an_expected_sha_matching_an_abbreviated_head_is_accepted() -> None:
    """Callers name commits the way git log prints them, seven characters at a time."""
    assert admissibility.verify_revision({"head": "ec67fff" + "0" * 33}, "ec67fff", side="after") == []


def test_an_expected_sha_that_disagrees_with_the_recorded_head_is_reported() -> None:
    """This is the check that makes "attests commit X" falsifiable."""
    problems = admissibility.verify_revision({"head": "b" * 40}, "ec67fff", side="after")
    assert problems == [f"after: expected ec67fff, but the side ran at {'b' * 40}"]


def test_expecting_a_sha_from_a_side_that_recorded_none_is_reported() -> None:
    """A run root from the 1.0.0 harness cannot satisfy an ``--expect-after-sha``."""
    problems = admissibility.verify_revision(None, "ec67fff", side="after")
    assert len(problems) == 1
    assert "recorded no revision" in problems[0]


def test_no_expected_sha_means_no_check() -> None:
    """Recording is unconditional; checking is opt-in."""
    assert admissibility.verify_revision(None, None, side="before") == []
