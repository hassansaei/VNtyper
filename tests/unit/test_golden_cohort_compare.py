"""Unit tests for the golden-cohort gate's comparison and its verdict.

The headline is ``test_two_sides_that_both_produced_nothing_are_not_identical``. Before the
expectation check existed, a case where both sides died at exit 1 before writing a single
genotype artefact produced nothing but ``absent_both`` statuses; ``diff_case`` excludes
those from the delta list, ``_summarise`` excludes them from the compared count, and
``_verdict`` therefore returned ``IDENTICAL`` over two runs that measured nothing. The gate
exited 0 and the page recorded a PASS.

The rest pin the verdict ordering, the ``REDUCED`` word that keeps a smoke run from earning
an attestation verdict, and the two comparison narrowings this lane repaired: a changed
provenance banner that was computed and then discarded, and a cohort-order justification
that cited a test which is not in the repository.
"""

import json
import sys
from pathlib import Path
from typing import Any

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from golden_cohort import HARNESS_VERSION, artifacts, compare, normalise  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]

#: One side's cases in ``_write_side`` shape: case id -> (exit code, expectation
#: problems, relative path -> file contents). Named because the empty-list and
#: empty-dict literals below are otherwise ambiguous to mypy at every call site.
CaseSpec = dict[str, tuple[int, list[str], dict[str, str]]]


def _table(rows, columns=("Motifs", "POS"), provenance=("## VNtyper Version: <VERSION>",)):
    """Build a parsed TSV as ``read_tsv`` returns one.

    Args:
        rows: The data rows.
        columns: The header.
        provenance: The ``##`` banner lines, already normalised.

    Returns:
        dict: The parsed table.
    """
    return {"columns": list(columns), "rows": list(rows), "provenance": list(provenance)}


# --------------------------------------------------------------------------------------
# absent_both, and why it needs the expectation check above it
# --------------------------------------------------------------------------------------


def test_an_artefact_absent_on_both_sides_is_not_a_delta() -> None:
    """Pinned deliberately: 56 of 59 cases legitimately have no adVNTR output.

    Turning that into a delta would drown the signal, so this behaviour is correct and
    stays. What it cannot do is tell a legitimate absence from a catastrophic one - that
    judgement is made per case from the matrix, not here.
    """
    assert compare.diff_table(None, None, ())["status"] == "absent_both"
    assert compare.diff_sequence(None, None)["status"] == "absent_both"
    assert compare.diff_opaque(None, None)["status"] == "absent_both"


def test_a_case_where_every_artefact_is_absent_on_both_sides_has_no_deltas() -> None:
    """The exact shape of the failure: two empty runs produce an empty delta list."""
    empty = dict.fromkeys(compare.PIPELINE_ARTIFACTS)
    result = compare.diff_case(empty, empty, compare.PIPELINE_ARTIFACTS)
    assert result["deltas"] == []


def test_an_unmet_expectation_becomes_a_delta_in_its_own_right() -> None:
    """This is what makes the case above visible."""
    case = {"deltas": [], "expectation_problems": {"before": ["expected exit 0, got 1"]}}
    compare._fold_expectation_into_deltas(case)
    assert case["deltas"] == ["EXPECTATION"]


def test_a_met_expectation_adds_nothing() -> None:
    """The 65 healthy cases of run 4 must not each grow a spurious delta."""
    case: dict[str, Any] = {"deltas": [], "expectation_problems": {}}
    compare._fold_expectation_into_deltas(case)
    assert case["deltas"] == []


def test_the_expectation_delta_is_not_added_twice() -> None:
    """Folding is idempotent, so a re-render cannot inflate the delta list."""
    case = {"deltas": ["EXPECTATION"], "expectation_problems": {"after": ["boom"]}}
    compare._fold_expectation_into_deltas(case)
    assert case["deltas"] == ["EXPECTATION"]


# --------------------------------------------------------------------------------------
# the verdict
# --------------------------------------------------------------------------------------


def _summary(cases_with_delta=()):
    """Build the rolled-up summary shape ``_verdict`` reads.

    Args:
        cases_with_delta: Case ids carrying a delta.

    Returns:
        dict: The summary.
    """
    return {"per_artifact": {}, "cases_total": 1, "cases_with_any_delta": list(cases_with_delta)}


def test_a_clean_attestation_run_is_identical() -> None:
    """The verdict every recorded PASS rests on."""
    assert compare._verdict(_summary(), [], True, [], True) == "IDENTICAL"


def test_a_clean_run_over_a_reduced_matrix_is_not_identical() -> None:
    """A smoke run must not earn the attestation word. This is C4's second half.

    ``REDUCED`` is exactly as free of deltas as ``IDENTICAL`` and attests strictly less,
    and the two used to be indistinguishable in both the verdict and the exit code.
    """
    assert compare._verdict(_summary(), [], True, [], False) == "REDUCED"


def test_an_unmet_expectation_outranks_a_clean_comparison() -> None:
    """Every "same" beneath a case that produced nothing is comparing two absences."""
    assert compare._verdict(_summary(), [], True, ["probe_x"], True) == "EXPECTATIONS_UNMET"


def test_an_unmet_expectation_outranks_a_delta() -> None:
    """A run whose cases did not do what they declared is not a run with interesting deltas."""
    assert compare._verdict(_summary(["a"]), [], True, ["b"], True) == "EXPECTATIONS_UNMET"


def test_an_unverified_launch_outranks_an_unmet_expectation() -> None:
    """If the harness cannot say which package ran, nothing below it means anything."""
    assert compare._verdict(_summary(), [], False, ["b"], True) == "UNVERIFIED"


def test_a_blocked_case_outranks_everything() -> None:
    """A case that never ran is the worst finding, and stays the reported one."""
    assert compare._verdict(_summary(), ["cohort_multi"], False, ["b"], True) == "BLOCKED"


def test_a_delta_over_a_full_matrix_is_reported_as_deltas() -> None:
    """Run 4's own verdict, which the page correctly records as mechanical."""
    assert compare._verdict(_summary(["a"]), [], True, [], True) == "DELTAS"


# --------------------------------------------------------------------------------------
# provenance
# --------------------------------------------------------------------------------------


def test_a_changed_provenance_banner_makes_the_table_differ() -> None:
    """It was computed and then discarded: ``status`` stayed ``same`` and never reached the verdict.

    The banner is compared *after* normalisation, which has already replaced the analysis
    timestamp and the VNtyper version, so anything still differing in it is real.
    """
    before = _table([{"Motifs": "4-5", "POS": "1"}], provenance=["## VNtyper Kestrel result"])
    after = _table([{"Motifs": "4-5", "POS": "1"}], provenance=["## Something else entirely"])
    detail = compare.diff_table(before, after, ("Motifs", "POS"))
    assert detail["status"] == "differ"
    assert detail["provenance_changed"] is True
    assert detail["provenance_before"] == ["## VNtyper Kestrel result"]


def test_an_unchanged_provenance_banner_leaves_the_table_the_same() -> None:
    """Run 4 has zero provenance changes, so folding this in must cost it nothing."""
    rows = [{"Motifs": "4-5", "POS": "1"}]
    detail = compare.diff_table(_table(rows), _table(rows), ("Motifs", "POS"))
    assert detail["status"] == "same"
    assert detail["provenance_changed"] is False


def test_a_row_change_still_differs_independently_of_the_banner() -> None:
    """The pre-existing comparison is unaffected by the new term."""
    before = _table([{"Motifs": "4-5", "POS": "1"}])
    after = _table([{"Motifs": "5C-Q", "POS": "1"}])
    detail = compare.diff_table(before, after, ("Motifs", "POS"))
    assert detail["status"] == "differ"
    assert detail["provenance_changed"] is False


# --------------------------------------------------------------------------------------
# the cohort-order justification
# --------------------------------------------------------------------------------------


def test_the_cohort_order_justification_cites_only_tests_that_exist() -> None:
    """It cited ``test_discovery_returns_an_unordered_set_today``, which is not in the repo.

    That test was replaced when the determinism fix landed. A justification citing a test
    that does not exist reads as verified and is not; this reads the real test files and
    checks every name the note mentions is defined in one of them.

    Two files rather than one since #205: the directory-input evidence is in
    ``test_cohort_inputs.py`` and the ZIP evidence - which is what reason 2 of the note
    now turns on - is in ``test_cohort_zip_identity.py``. The whole cohort-identity suite
    is read rather than the two files that happen to be cited today, so splitting a module
    again moves a citation without breaking this guard; ``test_cohort_identity.py`` was
    1,210 lines and became four modules, and the citation it carried moved with the test.
    """
    sources = "\n".join(
        path.read_text(encoding="utf-8") for path in sorted((REPO_ROOT / "tests" / "unit").glob("test_cohort_*.py"))
    )
    why = compare.COHORT_ORDER_WHY
    cited = [token.strip(".,;()") for token in why.split() if token.startswith("::") or "test_" in token]
    names = {token.split("::")[-1] for token in cited if "test_" in token}
    names = {name for name in names if name.startswith("test_")}
    assert names, "the justification must cite at least one test"
    for name in names:
        assert f"def {name}(" in sources, (
            f"{name} is cited by the normalisation note but is in neither cohort test file"
        )


def test_the_cohort_order_justification_quotes_the_sort_key_discovery_actually_uses() -> None:
    """The note has to name the key discovery sorts on, and #205 changed that key.

    It used to grep for ``return sorted(processed_dirs), temp_dirs``, which is gone: the
    de-duplicating ``set`` became a mapping keyed on the sample directory and the sort
    moved onto ``DiscoveredSample.sort_key``. The claim being defended is unchanged -
    discovery is ordered, so the note must not say it is unordered - and what changed is
    which line of source establishes it.

    It moved a second time when identity qualification landed: the sorted list is now bound
    and handed to ``qualify_colliding_identities`` on the way out, so the ``sorted(...)``
    call and the ``return`` are two statements rather than one. The sort **expression** is
    what carries the claim, so that is what is grepped for, and the return is asserted
    separately - qualification rewrites identities and never touches ``sort_key``, so the
    order the note describes survives it. Matching the whole old line would have failed
    here for a change that cannot affect ordering at all.
    """
    source = (REPO_ROOT / "vntyper" / "scripts" / "cohort_inputs.py").read_text(encoding="utf-8")
    assert "sorted(processed_dirs.values(), key=lambda sample: sample.sort_key)" in source
    assert "return qualify_colliding_identities(ordered), temp_dirs" in source
    assert (
        "returns sorted() on an effective-path key of the parts of the input path followed by the path "
        "relative to that input's root today" in compare.COHORT_ORDER_WHY
    )


def test_the_cohort_order_justification_bounds_the_zip_reason_to_baselines_before_the_fix() -> None:
    """Reason 2 was "ZIP inputs, on any version", and #205 made that false.

    Extraction still goes to ``tempfile.mkdtemp(prefix="cohort_zip_")`` - the temporary
    directory did not go away, only its place in the sort key did - so the reason survives
    for a baseline predating the fix and must be stated with that boundary. A stale
    justification here is not cosmetic: it is the recorded reason this gate does not
    compare cohort order at all, so leaving it would keep normalising away a difference
    that no longer exists.
    """
    source = (REPO_ROOT / "vntyper" / "scripts" / "cohort_inputs.py").read_text(encoding="utf-8")
    why = compare.COHORT_ORDER_WHY
    assert 'tempfile.mkdtemp(prefix="cohort_zip_")' in source
    assert "a baseline predating #205 sorts ZIP samples on their extracted path" in why
    assert "version-bounded" in why
    assert "Neither reason describes the current tree" in why
    assert "on any version" not in why


def test_the_gate_states_that_it_does_not_attest_the_ordering_fix() -> None:
    """Normalising a difference away also makes the fix to that difference invisible."""
    assert "NOT by this gate" in compare.COHORT_ORDER_WHY


# --------------------------------------------------------------------------------------
# compare_sides, end to end on disk
# --------------------------------------------------------------------------------------


def _write_side(root: Path, label: str, tree: str, expect_marker: str, cases, *, head=None, attestation=True):
    """Write a complete, readable run root for one side.

    Args:
        root: The side's run root.
        label: ``before`` or ``after``.
        tree: The recorded source tree.
        expect_marker: ``present`` or ``absent``.
        cases: Mapping of case id to ``(exit_code, expectation_problems, artefact_files)``,
            where ``artefact_files`` maps a relative path to its content.
        head: A recorded commit.
        attestation: Whether the matrix records itself as attestation-grade.

    Returns:
        dict: The side record, as ``load_side`` would return it.
    """
    (root / "cases").mkdir(parents=True, exist_ok=True)
    (root / "logs").mkdir(parents=True, exist_ok=True)
    pipeline_results = {}
    matrix_cases = []
    for case_id, (exit_code, problems, files) in cases.items():
        case_dir = root / "cases" / case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        for name, content in files.items():
            target = case_dir / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content, encoding="utf-8")
        log_dir = root / "logs" / case_id
        log_dir.mkdir(parents=True, exist_ok=True)
        record = {
            "case_id": case_id,
            "exit_code": exit_code,
            "launch_line": f"GATE-LAUNCH side={label}",
            "launch_verified": True,
            "aborted": False,
            "timed_out": False,
            "expectation_met": not problems,
            "expectation_problems": list(problems),
            "missing_artifacts": [],
        }
        (log_dir / "result.json").write_text(json.dumps(record), encoding="utf-8")
        pipeline_results[case_id] = record
        matrix_cases.append({"case_id": case_id, "group": "base", "expect_exit": "zero"})

    (root / "matrix.json").write_text(
        json.dumps(
            {
                "data_dir": "/shared/tests/data",
                "cases": matrix_cases,
                "probes": [],
                "cohort_cases": [],
                "check": {"counts": {"total": len(matrix_cases)}, "mismatches": [], "attestation_grade": attestation},
            }
        ),
        encoding="utf-8",
    )
    side = {
        "harness_version": "1.1.0",
        "side": label,
        "tree": tree,
        "run_root": str(root),
        "marker": "vntyper.scripts.cohort_rules",
        "expect_marker": expect_marker,
        "launch_verified": True,
        "cases_launched": len(cases),
        "expectations_unmet": sorted(cid for cid, (_e, problems, _f) in cases.items() if problems),
        "expectations_met": not any(problems for _e, problems, _f in cases.values()),
        "pipeline_results": pipeline_results,
        "cohort_results": {},
    }
    if head:
        side["revision"] = {"head": head, "branch": "b", "dirty_relevant": False, "dirty_paths": []}
    (root / "side.json").write_text(json.dumps(side), encoding="utf-8")
    return side


def _write_commands(root: Path, case_id: str, commands) -> None:
    """Write one side's recorded command stream for a case.

    Args:
        root: The side's run root.
        case_id: The case whose commands these are.
        commands: The command strings, in execution order.
    """
    log_dir = root / "logs" / case_id
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "commands.jsonl").write_text(
        "".join(json.dumps({"command": command, "shell": True}) + "\n" for command in commands),
        encoding="utf-8",
    )


def _compare(tmp_path: Path, before_cases, after_cases, **kwargs):
    """Run ``compare_sides`` over two written run roots.

    Args:
        tmp_path: pytest's per-test directory.
        before_cases: The baseline side's cases, in ``_write_side`` shape.
        after_cases: The candidate side's cases.
        **kwargs: Passed to both ``_write_side`` calls.

    Returns:
        dict: The comparison document.
    """
    before_root, after_root = tmp_path / "before", tmp_path / "after"
    before = _write_side(before_root, "before", "/trees/base", "absent", before_cases, head="a" * 40, **kwargs)
    after = _write_side(after_root, "after", "/trees/cand", "present", after_cases, head="b" * 40, **kwargs)
    rules = normalise.build_rules(source_root=Path("/trees/x"), run_root=before_root)
    return compare.compare_sides(before_root, after_root, before, after, normalise.manifest(rules), rules, rules)


_GOOD_FILES = {
    "pipeline_summary.json": json.dumps({"steps": []}),
    "kestrel/kestrel_result.tsv": "## VNtyper\nMotifs\tPOS\n4-5\t1\n",
}


def test_two_healthy_sides_compare_identical(tmp_path: Path) -> None:
    """The control: the fixture must be able to produce a PASS.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    result = _compare(tmp_path, cases, cases)
    assert result["verdict"] == "IDENTICAL"


def test_cram_read_set_and_raw_loss_evidence_are_compared(tmp_path: Path) -> None:
    """H7: A-178-2 evidence must not remain an ignored result.json field."""
    before_root, after_root = tmp_path / "before", tmp_path / "after"
    cases: CaseSpec = {"cram": (1, [], {})}
    before = _write_side(before_root, "before", "/trees/base", "absent", cases, head="a" * 40)
    after = _write_side(after_root, "after", "/trees/cand", "present", cases, head="b" * 40)
    before_result = before_root / "logs" / "cram" / "result.json"
    after_result = after_root / "logs" / "cram" / "result.json"
    before_payload = json.loads(before_result.read_text())
    after_payload = json.loads(after_result.read_text())
    before_payload.update(
        {
            "unmapped_read_set": None,
            "raw_indexed_read_set": {"count": 1},
            "raw_indexed_loss": None,
            "placed_unmapped_guard_count": 11_571,
        }
    )
    after_payload.update(
        {
            "unmapped_read_set": {"count": 4},
            "raw_indexed_read_set": {"count": 1},
            "raw_indexed_loss": 3,
            "placed_unmapped_guard_count": 329,
        }
    )
    before_result.write_text(json.dumps(before_payload))
    after_result.write_text(json.dumps(after_payload))
    rules = normalise.build_rules(source_root=Path("/trees/x"), run_root=before_root)

    result = compare.compare_sides(before_root, after_root, before, after, normalise.manifest(rules), rules, rules)

    assert result["cases"]["cram"]["artefacts"]["unmapped_read_set"]["status"] == "added_after"
    assert result["cases"]["cram"]["artefacts"]["raw_indexed_read_set"]["status"] == "same"
    assert result["cases"]["cram"]["artefacts"]["raw_indexed_loss"]["status"] == "added_after"
    guard_diff = result["cases"]["cram"]["artefacts"]["placed_unmapped_guard_count"]
    assert guard_diff == {"status": "differ", "before": 11_571, "after": 329}


def test_pipeline_artifact_reader_retains_the_causal_guard_count(tmp_path: Path) -> None:
    """The parsed guard cannot disappear between result.json and comparison."""
    output_dir = tmp_path / "output"
    log_dir = tmp_path / "logs"
    output_dir.mkdir()
    log_dir.mkdir()
    (log_dir / "result.json").write_text(json.dumps({"placed_unmapped_guard_count": 11_571}))

    parsed = artifacts.read_pipeline_case(output_dir, log_dir, [])

    assert parsed["placed_unmapped_guard_count"] == 11_571


def test_evidence_validation_changes_bump_the_harness_minor_version() -> None:
    """Tightening which CRAM runs can earn attestation changes the instrument.

    1.5.0 widens the compared set with ``output_bed`` and adds the declared
    command-delta mode (#262). Both change what a recorded result means, so a result
    document must be able to say which instrument produced it.
    """
    assert HARNESS_VERSION == "1.5.0"


def test_two_sides_that_both_produced_nothing_are_not_identical(tmp_path: Path) -> None:
    """The gate's worst failure, end to end.

    Both sides exit 1 and write no artefact at all. Every compared artefact is
    ``absent_both``, which ``diff_case`` correctly excludes from the delta list - so before
    the expectation check the verdict was ``IDENTICAL``, ``cmd_compare`` returned 0, and a
    run that measured nothing was recorded as a PASS.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases: CaseSpec = {"a5c1_hg19_subset": (1, ["expected exit 0, got 1"], {})}
    result = _compare(tmp_path, cases, cases)
    assert result["verdict"] == "EXPECTATIONS_UNMET"
    assert result["expectations_unmet"] == ["a5c1_hg19_subset"]
    assert result["cases"]["a5c1_hg19_subset"]["deltas"] == ["EXPECTATION"]


def test_a_case_that_exited_zero_and_wrote_nothing_on_one_side_is_caught(tmp_path: Path) -> None:
    """The asymmetric version, which the artefact comparison would also have caught.

    Args:
        tmp_path: pytest's per-test directory.
    """
    result = _compare(
        tmp_path,
        {"a5c1_hg19_subset": (0, [], _GOOD_FILES)},
        {"a5c1_hg19_subset": (0, ["exited as expected but did not write 5 of 5 required artefact(s)"], {})},
    )
    assert result["verdict"] == "EXPECTATIONS_UNMET"
    assert "after" in result["cases"]["a5c1_hg19_subset"]["expectation_problems"]


def test_the_result_records_each_side_s_revision(tmp_path: Path) -> None:
    """``side.json`` used to record a path, which is not a revision.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    result = _compare(tmp_path, cases, cases)
    assert result["before"]["revision"]["head"] == "a" * 40
    assert result["after"]["revision"]["head"] == "b" * 40


def test_a_side_written_by_the_older_harness_reports_no_revision_rather_than_a_blank(tmp_path: Path) -> None:
    """None is honest; an empty string would render as a commit nobody can check.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before_root, after_root = tmp_path / "before", tmp_path / "after"
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    before = _write_side(before_root, "before", "/trees/base", "absent", cases)
    after = _write_side(after_root, "after", "/trees/cand", "present", cases)
    rules = normalise.build_rules(source_root=Path("/trees/x"), run_root=before_root)
    result = compare.compare_sides(before_root, after_root, before, after, normalise.manifest(rules), rules, rules)
    assert result["before"]["revision"] is None
    assert "revision not recorded" in compare.render_text(result)


def test_a_clean_run_over_a_reduced_matrix_reports_reduced_end_to_end(tmp_path: Path) -> None:
    """The matrix's own ``attestation_grade`` flag reaches the verdict.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    result = _compare(tmp_path, cases, cases, attestation=False)
    assert result["verdict"] == "REDUCED"
    assert result["attestation_grade"] is False
    assert "Not attestation-grade" in compare.render_text(result)


def test_a_run_root_from_before_attestation_grade_existed_is_not_downgraded(tmp_path: Path) -> None:
    """A 1.0.0 matrix has no such key, and must not be retroactively marked reduced.

    Args:
        tmp_path: pytest's per-test directory.
    """
    before_root, after_root = tmp_path / "before", tmp_path / "after"
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    before = _write_side(before_root, "before", "/trees/base", "absent", cases)
    after = _write_side(after_root, "after", "/trees/cand", "present", cases)
    for root in (before_root, after_root):
        data = json.loads((root / "matrix.json").read_text())
        data["check"].pop("attestation_grade")
        (root / "matrix.json").write_text(json.dumps(data), encoding="utf-8")
    rules = normalise.build_rules(source_root=Path("/trees/x"), run_root=before_root)
    result = compare.compare_sides(before_root, after_root, before, after, normalise.manifest(rules), rules, rules)
    assert result["verdict"] == "IDENTICAL"


def test_the_rendered_report_names_the_unmet_expectations(tmp_path: Path) -> None:
    """A verdict word nobody can act on is half a finding.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases: CaseSpec = {"a5c1_hg19_subset": (1, ["expected exit 0, got 1"], {})}
    text = compare.render_text(_compare(tmp_path, cases, cases))
    assert "Cases that did not do what the matrix declared" in text
    assert "`a5c1_hg19_subset` (before): expected exit 0, got 1" in text


def test_the_rendered_report_does_not_call_the_whole_matrix_derived(tmp_path: Path) -> None:
    """Only the base cases are derived; the rest is declared policy.

    Args:
        tmp_path: pytest's per-test directory.
    """
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    text = compare.render_text(_compare(tmp_path, cases, cases))
    assert "selected by declared policy" in text


# ---------------------------------------------------------------------------
# The cohort statistics exports (#172), and why widening the gate needs its own test
# ---------------------------------------------------------------------------


def test_the_cohort_statistics_exports_are_compared() -> None:
    """``cohort_stats_{csv,tsv,json}`` must be in the compared set, not merely written.

    #172 adds a third cohort export carrying every ``cov_*`` column, ``coverage_qc``
    among them. A new output the harness writes but never reads is a gate that narrowed
    relative to the product: the three files would show up in ``cohort_output_files`` and
    their *contents* would never be diffed. That is the same failure this module's
    docstring describes for an artefact that stops being written.
    """
    for name in ("cohort_stats_csv", "cohort_stats_tsv", "cohort_stats_json"):
        assert name in compare.COHORT_ARTIFACTS, f"{name} is written but not compared"

    assert compare.COHORT_ARTIFACTS["cohort_stats_csv"] == ("table", ("Sample",))
    assert compare.COHORT_ARTIFACTS["cohort_stats_json"] == ("opaque", ())


def test_the_cohort_call_frequency_exports_are_compared() -> None:
    """``cohort_call_frequency_{csv,tsv,json}`` must be in the compared set, not merely written.

    Issue #33 adds the fourth cohort export carrying call frequencies across the roster.
    A new output the harness writes but never reads is a gate that narrowed relative to the
    product.
    """
    for name in ("cohort_call_frequency_csv", "cohort_call_frequency_tsv", "cohort_call_frequency_json"):
        assert name in compare.COHORT_ARTIFACTS, f"{name} is written but not compared"

    assert compare.COHORT_ARTIFACTS["cohort_call_frequency_csv"] == ("table", ("Grouping_Key",))
    assert compare.COHORT_ARTIFACTS["cohort_call_frequency_tsv"] == ("table", ("Grouping_Key",))
    assert compare.COHORT_ARTIFACTS["cohort_call_frequency_json"] == ("opaque", ())


def test_a_statistics_export_missing_from_one_side_is_a_delta_not_a_skip() -> None:
    """The candidate failing to write the export must fail the gate, not pass it quietly.

    ``absent_both`` is legitimately excluded from the delta list, so a widening that is
    never exercised can report a pass it did not measure. One-sided absence is the case
    that has to bite.
    """
    before = [{"Sample": "s1", "cov_coverage_qc": "PASS"}]

    assert compare.diff_table(before, None, ("Sample",))["status"] != "absent_both"
    assert compare.diff_table(before, None, ("Sample",))["status"] != "identical"
    assert compare.diff_table(None, before, ("Sample",))["status"] != "identical"


# ---------------------------------------------------------------------------
# The declared command delta (#262)
# ---------------------------------------------------------------------------


def _command_case(before_commands, after_commands, **artefacts):
    """Diff one case built from two command streams plus optional artefacts.

    Args:
        before_commands: The baseline side's ``executed_commands``.
        after_commands: The candidate side's ``executed_commands``.
        **artefacts: ``name=(before_value, after_value)`` for any other artefact.

    Returns:
        dict: The :func:`compare.diff_case` result.
    """
    before: dict[str, Any] = {"executed_commands": before_commands}
    after: dict[str, Any] = {"executed_commands": after_commands}
    for name, (left, right) in artefacts.items():
        before[name] = left
        after[name] = right
    return compare.diff_case(before, after, compare.PIPELINE_ARTIFACTS)


def test_a_declared_command_delta_is_reported_but_not_fatal() -> None:
    """Performance work changes the command stream on purpose.

    The delta stays in ``deltas`` so the report always names it; only its fatality
    changes. A mode that hid the delta would be worse than no mode at all, because the
    whole point is that a human reads what changed.
    """
    case = _command_case(["a", "b"], ["a", "x", "b"])

    assert "executed_commands" in case["deltas"]
    assert compare.fatal_deltas(case, expect_command_delta=True) == []


def test_an_undeclared_command_delta_is_still_fatal() -> None:
    """Without the declaration the gate keeps its existing meaning exactly."""
    case = _command_case(["a"], ["a", "b"])

    assert compare.fatal_deltas(case, expect_command_delta=False) == ["executed_commands"]


def test_a_result_delta_is_always_fatal_even_when_a_command_delta_is_declared() -> None:
    """The exemption is one artefact wide; a genotype change is never excusable."""
    tables = (_table([{"Motifs": "4-5", "POS": "1"}]), _table([{"Motifs": "4-5", "POS": "2"}]))
    case = _command_case(["a"], ["a", "b"], kestrel_result=tables)

    assert compare.fatal_deltas(case, expect_command_delta=True) == ["kestrel_result"]


def test_a_non_genotype_delta_is_still_fatal_under_the_declaration() -> None:
    """``report_tables`` and ``pipeline_steps`` are deliberate gate failures today.

    An exemption phrased as "fatal only if genotype-bearing" would silently weaken
    every one of them, which is why the rule is the inverse: everything stays fatal
    except ``executed_commands``.
    """
    case = _command_case(["a"], ["a"], report_tables=({"x": 1}, {"x": 2}))

    assert compare.fatal_deltas(case, expect_command_delta=True) == ["report_tables"]


def test_the_expectation_delta_is_fatal_under_the_declaration() -> None:
    """A case that did not do what it declared is not excused by a command delta."""
    case = _command_case(["a", "b"], ["a", "x", "b"])
    case["expectation_problems"] = {"before": ["expected exit 0, got 1"]}
    compare._fold_expectation_into_deltas(case)

    assert compare.fatal_deltas(case, expect_command_delta=True) == ["EXPECTATION"]


def test_output_bed_is_compared_as_an_ordered_line_sequence() -> None:
    """BED is headerless, so ``read_tsv`` would eat the first variant row as a header."""
    assert compare.PIPELINE_ARTIFACTS["output_bed"] == ("sequence", ())


def test_a_changed_bed_interval_is_a_delta() -> None:
    """#203 moved every BED interval by one base and no gate would have seen it."""
    case = _command_case(
        ["a"],
        ["a"],
        output_bed=(["motif-1\t66\t67"], ["motif-1\t67\t68"]),
    )

    assert compare.fatal_deltas(case, expect_command_delta=True) == ["output_bed"]


def test_the_bed_collector_reads_the_lines_it_finds(tmp_path: Path) -> None:
    """Collected as lines, in file order, with no header interpretation."""
    case_dir = tmp_path / "case"
    (case_dir / "kestrel").mkdir(parents=True)
    (case_dir / "kestrel" / "output.bed").write_text("m-1\t66\t67\nm-2\t10\t11\n", encoding="utf-8")
    rules = normalise.build_rules(source_root=tmp_path, run_root=tmp_path)

    collected = artifacts.read_pipeline_case(case_dir, tmp_path / "logs", rules)

    assert collected["output_bed"] == ["m-1\t66\t67", "m-2\t10\t11"]


def test_an_absent_bed_collects_as_none_rather_than_an_empty_list(tmp_path: Path) -> None:
    """A negative case writes no BED at all, and ``None`` is what makes that absent_both.

    An empty list would compare equal to a *present* but empty file, which is a
    different run.
    """
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    rules = normalise.build_rules(source_root=tmp_path, run_root=tmp_path)

    assert artifacts.read_pipeline_case(case_dir, tmp_path / "logs", rules)["output_bed"] is None


def test_an_empty_bed_collects_as_an_empty_list_not_none(tmp_path: Path) -> None:
    """Present-and-empty must be distinguishable from never-written."""
    case_dir = tmp_path / "case"
    (case_dir / "kestrel").mkdir(parents=True)
    (case_dir / "kestrel" / "output.bed").write_text("", encoding="utf-8")
    rules = normalise.build_rules(source_root=tmp_path, run_root=tmp_path)

    assert artifacts.read_pipeline_case(case_dir, tmp_path / "logs", rules)["output_bed"] == []


def test_a_declared_command_delta_makes_the_whole_run_identical(tmp_path: Path) -> None:
    """End to end: the flag is what lets a performance change be validated at all."""
    before: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    after: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    before_root, after_root = tmp_path / "before", tmp_path / "after"
    before_side = _write_side(before_root, "before", "/trees/base", "absent", before, head="a" * 40)
    after_side = _write_side(after_root, "after", "/trees/cand", "present", after, head="b" * 40)
    _write_commands(before_root, "a5c1_hg19_subset", ["samtools view -h in.bam | samtools view -b -f 4 -"])
    _write_commands(after_root, "a5c1_hg19_subset", ["samtools view -b -f 4 -u in.bam -o u.bam"])
    rules = normalise.build_rules(source_root=Path("/trees/x"), run_root=before_root)

    undeclared = compare.compare_sides(
        before_root, after_root, before_side, after_side, normalise.manifest(rules), rules, rules
    )
    declared = compare.compare_sides(
        before_root,
        after_root,
        before_side,
        after_side,
        normalise.manifest(rules),
        rules,
        rules,
        expect_command_delta=True,
    )

    assert undeclared["verdict"] == "DELTAS"
    assert declared["verdict"] == "IDENTICAL"
    assert declared["expect_command_delta"] is True
    assert declared["cases"]["a5c1_hg19_subset"]["deltas"] == ["executed_commands"]


def test_a_result_delta_still_fails_the_run_when_the_command_delta_is_declared(tmp_path: Path) -> None:
    """The declaration must not become a blanket pass."""
    good: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    changed: CaseSpec = {
        "a5c1_hg19_subset": (
            0,
            [],
            {**_GOOD_FILES, "kestrel/kestrel_result.tsv": "## VNtyper\nMotifs\tPOS\n4-5\t2\n"},
        )
    }
    before_root, after_root = tmp_path / "before", tmp_path / "after"
    before_side = _write_side(before_root, "before", "/trees/base", "absent", good, head="a" * 40)
    after_side = _write_side(after_root, "after", "/trees/cand", "present", changed, head="b" * 40)
    _write_commands(before_root, "a5c1_hg19_subset", ["samtools view -h in.bam"])
    _write_commands(after_root, "a5c1_hg19_subset", ["samtools view -b in.bam"])
    rules = normalise.build_rules(source_root=Path("/trees/x"), run_root=before_root)

    result = compare.compare_sides(
        before_root,
        after_root,
        before_side,
        after_side,
        normalise.manifest(rules),
        rules,
        rules,
        expect_command_delta=True,
    )

    assert result["verdict"] == "DELTAS"
    assert result["cases"]["a5c1_hg19_subset"]["fatal_deltas"] == ["kestrel_result"]


def test_the_report_says_which_cases_had_a_waived_command_delta(tmp_path: Path) -> None:
    """``IDENTICAL`` must never be printed without saying what was waived to earn it.

    This is the same failure ``REDUCED`` exists to prevent: a clean word over a run that
    attests strictly less than the word implies. The declaration is only worth anything
    if the reader is told which cases used it and can go and read the diff.
    """
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    before_root, after_root = tmp_path / "before", tmp_path / "after"
    before_side = _write_side(before_root, "before", "/trees/base", "absent", cases, head="a" * 40)
    after_side = _write_side(after_root, "after", "/trees/cand", "present", cases, head="b" * 40)
    _write_commands(before_root, "a5c1_hg19_subset", ["samtools view -h in.bam"])
    _write_commands(after_root, "a5c1_hg19_subset", ["samtools view -b in.bam"])
    rules = normalise.build_rules(source_root=Path("/trees/x"), run_root=before_root)

    text = compare.render_text(
        compare.compare_sides(
            before_root,
            after_root,
            before_side,
            after_side,
            normalise.manifest(rules),
            rules,
            rules,
            expect_command_delta=True,
        )
    )

    assert compare.DECLARED_DELTA_HEADING in text
    assert "a5c1_hg19_subset" in text


def test_the_report_does_not_mention_a_waiver_when_none_was_declared(tmp_path: Path) -> None:
    """The default run's report must read exactly as it did before this mode existed."""
    cases: CaseSpec = {"a5c1_hg19_subset": (0, [], _GOOD_FILES)}
    result = _compare(tmp_path, cases, cases)

    assert compare.DECLARED_DELTA_HEADING not in compare.render_text(result)
