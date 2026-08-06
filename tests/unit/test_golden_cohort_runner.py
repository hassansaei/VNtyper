"""Unit tests for the golden-cohort gate's side runner and its command-line entry point.

Two things are pinned here that no other test reaches:

* what a side records about itself - its revision, how many cases it launched, and whether
  they did what the matrix declared. ``side.json`` used to record a *path* and a
  ``launch_verified`` computed as ``all()`` over the launched cases, which is True for a
  side that launched none;
* what ``compare`` refuses. The launcher proves each process ran the tree its caller named;
  nothing proved the caller named two different, oppositely-labelled trees, so "compared
  against itself" was reachable by misconfiguration and returned exit 0.

Nothing here launches a pipeline: ``_run_one`` and the git call are both replaced.
"""

import importlib.util
import json
import sys
from pathlib import Path
from unittest import mock

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from golden_cohort import runner  # noqa: E402


def _load_gate_entry_point():
    """Import ``scripts/golden_cohort_gate.py`` as a module without running it.

    It is a script rather than a package member, so it is loaded by path. The module-level
    ``sys.path`` insert it performs is the same one this file already did.

    Returns:
        ModuleType: The loaded entry point.
    """
    spec = importlib.util.spec_from_file_location("golden_cohort_gate", REPO_ROOT / "scripts" / "golden_cohort_gate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = _load_gate_entry_point()

CLEAN_REVISION = {
    "head": "e" * 40,
    "branch": "main",
    "dirty": False,
    "dirty_paths": [],
    "dirty_relevant": False,
    "revision_paths": ["vntyper", "docker", "scripts"],
    "error": None,
}


def _case(case_id, **extra):
    """Build one pipeline case in the shape ``pipeline_argv`` reads.

    Args:
        case_id: The case id.
        **extra: Further keys to set.

    Returns:
        dict: The case.
    """
    return {
        "case_id": case_id,
        "kind": "pipeline",
        "group": "base",
        "bam": f"/data/example_{case_id}.bam",
        "assembly": "hg19",
        "fast_mode": True,
        "advntr": False,
        "expect_exit": "zero",
        "required_artifacts": [],
        **extra,
    }


def _matrix(cases):
    """Build the minimum matrix ``run_side`` reads.

    Args:
        cases: The pipeline cases.

    Returns:
        dict: The matrix.
    """
    return {"cases": list(cases), "probes": [], "cohort_cases": [], "check": {}}


def _run_side(tmp_path: Path, cases, records, revision=CLEAN_REVISION):
    """Run one side with the launch replaced by canned records.

    Args:
        tmp_path: pytest's per-test directory.
        cases: The matrix's pipeline cases.
        records: Mapping of case id to the record ``_run_one`` would have returned.
        revision: What ``describe_tree`` reports.

    Returns:
        dict: The side record.
    """

    def fake_run_one(*, case, **_kwargs):
        return records[case["case_id"]]

    with (
        mock.patch.object(runner, "_run_one", side_effect=fake_run_one),
        mock.patch.object(runner.admissibility, "describe_tree", return_value=revision),
    ):
        return runner.run_side(
            matrix=_matrix(cases),
            tree=tmp_path / "tree",
            run_root=tmp_path / "run",
            side="after",
            marker="vntyper.scripts.cohort_rules",
            expect_marker=True,
            threads=4,
            advntr_threads=8,
            jobs=1,
            timeout=60,
            skip_cohort=True,
        )


def _ok_record(case_id):
    """A record for a case that ran and met its expectation.

    Args:
        case_id: The case id.

    Returns:
        dict: The record.
    """
    return {
        "case_id": case_id,
        "exit_code": 0,
        "launch_verified": True,
        "expectation_met": True,
        "expectation_problems": [],
    }


# --------------------------------------------------------------------------------------
# what a side records
# --------------------------------------------------------------------------------------


def test_a_side_records_the_revision_it_ran(tmp_path: Path) -> None:
    """Record the revision, because "this run attests commit X" must not be an assertion.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases = [_case("a")]
    record = _run_side(tmp_path, cases, {"a": _ok_record("a")})
    assert record["revision"]["head"] == "e" * 40
    on_disk = json.loads((tmp_path / "run" / "side.json").read_text(encoding="utf-8"))
    assert on_disk["revision"]["head"] == "e" * 40
    assert on_disk["revision"]["dirty_relevant"] is False


def test_a_side_whose_cases_all_met_their_expectations_says_so(tmp_path: Path) -> None:
    """The healthy path, which run 4 took on all 65 cases on both sides.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases = [_case("a")]
    record = _run_side(tmp_path, cases, {"a": _ok_record("a")})
    assert record["expectations_met"] is True
    assert record["expectations_unmet"] == []
    assert record["cases_launched"] == 1


def test_a_side_with_an_unmet_expectation_reports_it(tmp_path: Path) -> None:
    """And ``cmd_run`` turns that into exit 1, so a broken side never reaches ``compare``.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases = [_case("a")]
    failed = {**_ok_record("a"), "expectation_met": False, "expectation_problems": ["expected exit 0, got 1"]}
    record = _run_side(tmp_path, cases, {"a": failed})
    assert record["expectations_met"] is False
    assert record["expectations_unmet"] == ["a"]


def test_a_side_that_launched_nothing_is_not_verified(tmp_path: Path) -> None:
    """``all()`` over an empty mapping is True, so an empty side used to report itself verified.

    ``build_matrix`` now refuses a zero-case matrix, which is the first lock on this door.
    This is the second: even reached some other way, a side that ran nothing reports
    ``launch_verified: false``, so the comparison's verdict is ``UNVERIFIED``.

    Args:
        tmp_path: pytest's per-test directory.
    """
    record = _run_side(tmp_path, [], {})
    assert record["cases_launched"] == 0
    assert record["launch_verified"] is False


def test_a_side_records_a_failed_revision_lookup_rather_than_losing_the_run(tmp_path: Path) -> None:
    """A tree that is not a git checkout is a fact about the run, not a reason to lose it.

    Args:
        tmp_path: pytest's per-test directory.
    """
    broken = {**CLEAN_REVISION, "head": None, "branch": None, "error": "not a git repository"}
    cases = [_case("a")]
    record = _run_side(tmp_path, cases, {"a": _ok_record("a")}, revision=broken)
    assert record["revision"]["head"] is None
    assert record["revision"]["error"] == "not a git repository"
    assert record["launch_verified"] is True


# --------------------------------------------------------------------------------------
# what compare refuses
# --------------------------------------------------------------------------------------


def _write_minimal_side(root: Path, label: str, tree: str, expect_marker: str, head=None) -> None:
    """Write just enough of a run root for ``cmd_compare`` to load and admit it.

    One case result is included deliberately: a side with none is refused for that on its
    own, which would make every refusal test below pass for the wrong reason.

    Args:
        root: The run root.
        label: ``before`` or ``after``.
        tree: The recorded source tree.
        expect_marker: ``present`` or ``absent``.
        head: A recorded commit, or None.
    """
    root.mkdir(parents=True, exist_ok=True)
    (root / "cases" / "a").mkdir(parents=True, exist_ok=True)
    (root / "logs" / "a").mkdir(parents=True, exist_ok=True)
    result = {"case_id": "a", "exit_code": 0, "launch_verified": True, "expectation_problems": []}
    (root / "logs" / "a" / "result.json").write_text(json.dumps(result), encoding="utf-8")
    (root / "matrix.json").write_text(json.dumps({"cases": [], "probes": [], "check": {}}), encoding="utf-8")
    side = {
        "side": label,
        "tree": tree,
        "marker": "vntyper.scripts.cohort_rules",
        "expect_marker": expect_marker,
        "launch_verified": True,
        "pipeline_results": {"a": result},
        "cohort_results": {},
    }
    if head:
        side["revision"] = {"head": head, "branch": "b", "dirty_relevant": False, "dirty_paths": []}
    (root / "side.json").write_text(json.dumps(side), encoding="utf-8")


def _compare_args(before: Path, after: Path, **overrides):
    """Build the argparse namespace ``cmd_compare`` reads.

    Args:
        before: The baseline run root.
        after: The candidate run root.
        **overrides: Fields to replace.

    Returns:
        argparse.Namespace: The arguments.
    """
    import argparse

    defaults = {
        "before_root": before,
        "after_root": after,
        "json_out": None,
        "text_out": None,
        "expect_before_sha": None,
        "expect_after_sha": None,
        "require_clean": False,
    }
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def test_comparing_a_run_root_with_itself_is_refused(tmp_path: Path) -> None:
    """It used to load the same side twice and report a confident IDENTICAL.

    Args:
        tmp_path: pytest's per-test directory.
    """
    root = tmp_path / "only"
    _write_minimal_side(root, "before", "/trees/base", "absent", head="a" * 40)
    assert gate.cmd_compare(_compare_args(root, root)) == 1


def test_comparing_two_sides_that_ran_the_same_tree_is_refused(tmp_path: Path) -> None:
    """The failure the whole resolution wrapper exists to prevent, reached by misconfiguration.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before, after = tmp_path / "before", tmp_path / "after"
    _write_minimal_side(before, "before", "/trees/same", "absent", head="a" * 40)
    _write_minimal_side(after, "after", "/trees/same", "present", head="b" * 40)
    assert gate.cmd_compare(_compare_args(before, after)) == 1


def test_a_recorded_revision_that_disagrees_with_the_expected_sha_is_refused(tmp_path: Path) -> None:
    """This is what makes "attests ``ec67fff`` and nothing after it" falsifiable.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before, after = tmp_path / "before", tmp_path / "after"
    _write_minimal_side(before, "before", "/trees/base", "absent", head="a" * 40)
    _write_minimal_side(after, "after", "/trees/cand", "present", head="b" * 40)
    assert gate.cmd_compare(_compare_args(before, after, expect_after_sha="ec67fff")) == 1


def test_a_dirty_side_is_refused_when_the_caller_requires_a_clean_tree(tmp_path: Path) -> None:
    """A run over uncommitted edits attests that commit plus edits, which is not a commit.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before, after = tmp_path / "before", tmp_path / "after"
    _write_minimal_side(before, "before", "/trees/base", "absent", head="a" * 40)
    _write_minimal_side(after, "after", "/trees/cand", "present", head="b" * 40)
    side = json.loads((after / "side.json").read_text(encoding="utf-8"))
    side["revision"]["dirty_relevant"] = True
    side["revision"]["dirty_paths"] = [" M vntyper/scripts/scoring.py"]
    side["revision"]["revision_paths"] = ["vntyper", "docker", "scripts"]
    (after / "side.json").write_text(json.dumps(side), encoding="utf-8")
    assert gate.cmd_compare(_compare_args(before, after, require_clean=True)) == 1


def test_a_dirty_side_still_compares_when_the_caller_does_not_require_a_clean_tree(tmp_path: Path) -> None:
    """Dirtiness is recorded and warned about unconditionally; refusing it is opt-in.

    Making it fatal by default would break every run from a worktree with a scratch file
    in it, which is most of them.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before, after = tmp_path / "before", tmp_path / "after"
    _write_minimal_side(before, "before", "/trees/base", "absent", head="a" * 40)
    _write_minimal_side(after, "after", "/trees/cand", "present", head="b" * 40)
    side = json.loads((after / "side.json").read_text(encoding="utf-8"))
    side["revision"]["dirty_relevant"] = True
    side["revision"]["dirty_paths"] = [" M vntyper/scripts/scoring.py"]
    (after / "side.json").write_text(json.dumps(side), encoding="utf-8")
    assert gate.cmd_compare(_compare_args(before, after, require_clean=False)) == 0


def test_a_side_that_ran_no_cases_at_all_is_refused(tmp_path: Path) -> None:
    """Two sides that measured nothing compare clean, for lack of anything to differ.

    ``build_matrix`` refuses a zero-case matrix and ``run_side`` reports a zero-case side as
    unverified, but neither applies to a run root written by something other than this
    harness. This is the check that does.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before, after = tmp_path / "before", tmp_path / "after"
    _write_minimal_side(before, "before", "/trees/base", "absent", head="a" * 40)
    _write_minimal_side(after, "after", "/trees/cand", "present", head="b" * 40)
    side = json.loads((after / "side.json").read_text(encoding="utf-8"))
    side["pipeline_results"] = {}
    (after / "side.json").write_text(json.dumps(side), encoding="utf-8")
    assert gate.cmd_compare(_compare_args(before, after)) == 1


def test_two_properly_opposed_sides_get_past_the_refusal(tmp_path: Path) -> None:
    """The control: an admissible pair reaches the comparison and produces a verdict.

    Without this the refusal tests above would all pass against a ``cmd_compare`` that
    refused everything.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before, after = tmp_path / "before", tmp_path / "after"
    _write_minimal_side(before, "before", "/trees/base", "absent", head="a" * 40)
    _write_minimal_side(after, "after", "/trees/cand", "present", head="b" * 40)
    json_out = tmp_path / "result.json"
    assert gate.cmd_compare(_compare_args(before, after, json_out=json_out)) == 0
    assert json.loads(json_out.read_text(encoding="utf-8"))["verdict"] == "IDENTICAL"
