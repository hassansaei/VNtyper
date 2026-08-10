"""Unit tests for the golden-cohort gate's side runner and its command-line entry point.

Three things are pinned here that no other test reaches:

* what a side records about itself - its revision, how many cases it launched, and whether
  they did what the matrix declared. ``side.json`` used to record a *path* and a
  ``launch_verified`` computed as ``all()`` over the launched cases, which is True for a
  side that launched none;
* what ``compare`` refuses. The launcher proves each process ran the tree its caller named;
  nothing proved the caller named two different, oppositely-labelled trees, so "compared
  against itself" was reachable by misconfiguration and returned exit 0;
* which input flag each case is launched with. ``pipeline_argv`` hardcoded ``--bam``, so a
  CRAM case added to the matrix would have been launched against a path that is not a BAM,
  or - worse, once the case also carried a BAM - launched as a BAM and compared clean while
  attesting nothing about the CRAM path.

Nothing here launches a pipeline: ``_run_one`` and the git call are both replaced.
"""

import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path
from unittest import mock

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from golden_cohort import runner  # noqa: E402
from golden_cohort.cram_cases import CRAM_READ_SET_EXPECTATIONS  # noqa: E402

from vntyper.scripts.cli_parser import build_parser  # noqa: E402


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
# which input flag a case is launched with
# --------------------------------------------------------------------------------------


def _argv(case, tmp_path: Path):
    """Build one case's ``vntyper pipeline`` argument list.

    Args:
        case: The case from the matrix.
        tmp_path: pytest's per-test directory, used as the output directory.

    Returns:
        list[str]: The argument list.
    """
    return runner.pipeline_argv(case, tmp_path / "out", threads=4, advntr_threads=8)


@pytest.mark.parametrize("group", ["base", "nonfast", "advntr"])
def test_a_non_cram_case_is_still_launched_with_bam(group: str, tmp_path: Path) -> None:
    """Every case predating the CRAM group must be launched exactly as before.

    ``alignment_kind`` is absent from all of them, so the default is what runs, and the
    default has to be ``--bam`` or the whole existing matrix changes meaning.

    Args:
        group: The case group.
        tmp_path: pytest's per-test directory.
    """
    argv = _argv(_case("a", group=group), tmp_path)
    assert "--bam" in argv
    assert "--cram" not in argv
    assert argv[argv.index("--bam") + 1] == "/data/example_a.bam"


def test_a_cram_case_is_launched_with_cram_and_the_fixture_path(tmp_path: Path) -> None:
    """``pipeline_argv`` hardcoded ``--bam``, so a CRAM case could not have been run at all.

    Args:
        tmp_path: pytest's per-test directory.
    """
    case = _case(
        "b178_hg19_cram",
        group="cram",
        alignment_kind="cram",
        cram="/data/cram/example_b178_hg19_subset.cram",
        source_bam="/data/example_b178_hg19_subset.bam",
        fast_mode=False,
    )
    del case["bam"]
    argv = _argv(case, tmp_path)
    assert "--cram" in argv
    assert "--bam" not in argv
    assert argv[argv.index("--cram") + 1] == "/data/cram/example_b178_hg19_subset.cram"
    assert "/data/example_b178_hg19_subset.bam" not in argv


def test_a_cram_case_is_launched_without_fast_mode(tmp_path: Path) -> None:
    """The CRAM unmapped-read extraction only runs when ``--fast-mode`` is off.

    Args:
        tmp_path: pytest's per-test directory.
    """
    case = _case("b178_hg19_cram", group="cram", alignment_kind="cram", cram="/data/x.cram", fast_mode=False)
    del case["bam"]
    assert "--fast-mode" not in _argv(case, tmp_path)


def test_a_cram_case_launches_with_a_complete_per_case_scan_config(tmp_path: Path) -> None:
    """The case's scan policy must reach VNtyper through replacement-safe config input."""
    tree = tmp_path / "tree"
    config = tree / "vntyper" / "config.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({"cram": {"unmapped_scan": "auto"}, "tools": {"samtools": "samtools"}}))
    case = _case("b178_indexed_cram", alignment_kind="cram", cram="/data/x.cram", unmapped_scan="indexed")
    del case["bam"]

    config_path, effective_mode = runner.materialize_case_config(tree, case, tmp_path / "logs")
    argv = runner.pipeline_argv(case, tmp_path / "out", threads=4, advntr_threads=8, config_path=config_path)

    assert effective_mode == "indexed"
    assert json.loads(config_path.read_text())["cram"]["unmapped_scan"] == "indexed"
    assert argv[argv.index("--config-path") + 1] == str(config_path)


def test_a_cram_case_places_global_config_before_pipeline_and_current_cli_accepts_it(tmp_path: Path) -> None:
    """Generated CRAM argv must obey the CLI's top-level-only config option contract."""
    config_path = tmp_path / "case-config.json"
    case = _case("b178_indexed_cram", alignment_kind="cram", cram="/data/x.cram", unmapped_scan="indexed")
    del case["bam"]

    argv = runner.pipeline_argv(case, tmp_path / "out", threads=4, advntr_threads=8, config_path=config_path)
    parsed = build_parser().parse_args(argv)

    assert argv[:3] == ["--config-path", str(config_path), "pipeline"]
    assert parsed.command == "pipeline"
    assert parsed.config_path == config_path


def test_side_expectation_is_materialized_before_admissibility() -> None:
    """A candidate refusal and a baseline success are both intentional outcomes."""
    case = {
        "case_id": "mixed",
        "expect_exit": "zero",
        "required_artifacts": ["pipeline_summary.json"],
        "side_expectations": {
            "before": {"expect_exit": "zero", "required_artifacts": ["pipeline_summary.json"]},
            "after": {"expect_exit": "nonzero", "required_artifacts": []},
        },
    }

    assert runner.materialize_side_expectation(case, "before") == {
        **case,
        "expect_exit": "zero",
        "required_artifacts": ["pipeline_summary.json"],
    }
    assert runner.materialize_side_expectation(case, "after") == {
        **case,
        "expect_exit": "nonzero",
        "required_artifacts": [],
    }


def test_side_expectation_rejects_a_missing_side() -> None:
    """A malformed differential declaration must not fall back to the legacy default."""
    case = {
        "case_id": "mixed",
        "side_expectations": {"before": {"expect_exit": "zero", "required_artifacts": []}},
    }

    with pytest.raises(ValueError, match="has no 'after' expectation"):
        runner.materialize_side_expectation(case, "after")


def test_a_cram_case_adds_a_missing_cram_section_to_its_complete_config_copy(tmp_path: Path) -> None:
    """Older target trees must still receive a replacement-safe complete config."""
    tree = tmp_path / "tree"
    config = tree / "vntyper" / "config.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({"tools": {"samtools": "samtools"}}))
    case = _case("b178_indexed_cram", alignment_kind="cram", cram="/data/x.cram", unmapped_scan="indexed")
    del case["bam"]

    config_path, effective_mode = runner.materialize_case_config(tree, case, tmp_path / "logs")

    assert effective_mode == "indexed"
    assert json.loads(config_path.read_text())["cram"] == {"unmapped_scan": "indexed"}


@pytest.mark.parametrize("override", ["auto", "indexed", "stream"])
def test_the_harness_environment_override_beats_a_cram_cases_mode(tmp_path: Path, monkeypatch, override: str) -> None:
    """Wave 3 can compare one override mode across all CRAM cases without touching product config."""
    tree = tmp_path / "tree"
    config = tree / "vntyper" / "config.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({"cram": {"unmapped_scan": "auto"}}))
    case = _case("b178_stream_cram", alignment_kind="cram", cram="/data/x.cram", unmapped_scan="stream")
    del case["bam"]
    monkeypatch.setenv("VNTYPER_CRAM_UNMAPPED_SCAN", override)

    config_path, effective_mode = runner.materialize_case_config(tree, case, tmp_path / "logs")

    assert effective_mode == override
    assert json.loads(config_path.read_text())["cram"]["unmapped_scan"] == override


@pytest.mark.parametrize("declared_mode", ["indexed", "stream"])
def test_the_harness_override_cannot_collapse_declared_a178_scan_cases(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    declared_mode: str,
) -> None:
    """The indexed and stream evidence cases must remain distinct measurements."""
    tree = tmp_path / "tree"
    config = tree / "vntyper" / "config.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({"cram": {"unmapped_scan": "auto"}}))
    case = _case(
        f"indexed_safe_{declared_mode}_cram",
        alignment_kind="cram",
        cram="/data/x.cram",
        unmapped_scan=declared_mode,
        cram_evidence_expectation={"indexed_authorized": True},
    )
    del case["bam"]
    monkeypatch.setenv("VNTYPER_CRAM_UNMAPPED_SCAN", "stream" if declared_mode == "indexed" else "indexed")

    config_path, effective_mode = runner.materialize_case_config(tree, case, tmp_path / "logs")

    assert effective_mode == declared_mode
    assert json.loads(config_path.read_text())["cram"]["unmapped_scan"] == declared_mode


@pytest.mark.parametrize(
    ("source_id", "case_prefix"),
    [("7a61_hg38_ensembl_bwa", "7a61_hg38_ensembl"), ("b178_hg19_subset", "b178_hg19")],
)
@pytest.mark.parametrize("declared_mode", ["indexed", "stream"])
def test_runner_preserves_side_wrapped_a178_scan_declarations_before_applying_an_override(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    source_id: str,
    case_prefix: str,
    declared_mode: str,
) -> None:
    """The real lossy CRAM cases cannot be collapsed before their after-side policy is selected."""
    tree = tmp_path / "tree"
    config = tree / "vntyper" / "config.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({"cram": {"unmapped_scan": "auto"}}))
    case_id = f"{case_prefix}_{declared_mode}_cram"
    case = _case(
        case_id,
        alignment_kind="cram",
        cram=f"/data/{source_id}.cram",
        unmapped_scan=declared_mode,
        side_expectations={
            "before": {"expect_exit": "zero", "required_artifacts": []},
            "after": {
                "expect_exit": "nonzero",
                "required_artifacts": [],
                "cram_evidence_expectation": CRAM_READ_SET_EXPECTATIONS[source_id],
            },
        },
    )
    del case["bam"]
    monkeypatch.setenv("VNTYPER_CRAM_UNMAPPED_SCAN", "stream" if declared_mode == "indexed" else "indexed")
    captured: dict[str, object] = {}

    def fake_run_one(*, case: dict[str, object], **_kwargs: object) -> dict[str, object]:
        captured.update(case)
        return _ok_record(str(case["case_id"]))

    with (
        mock.patch.object(runner, "_run_one", side_effect=fake_run_one),
        mock.patch.object(runner.admissibility, "describe_tree", return_value=CLEAN_REVISION),
    ):
        runner.run_side(
            matrix=_matrix([case]),
            tree=tree,
            run_root=tmp_path / "run",
            side="after",
            marker="vntyper.scripts.alignment_contract",
            expect_marker=True,
            threads=4,
            advntr_threads=8,
            jobs=1,
            timeout=60,
            skip_cohort=True,
        )

    assert captured["cram_evidence_expectation"] == CRAM_READ_SET_EXPECTATIONS[source_id]
    assert captured["effective_unmapped_scan"] == declared_mode
    case_config = json.loads(Path(str(captured["case_config_path"])).read_text())
    assert case_config["cram"]["unmapped_scan"] == declared_mode


def test_an_invalid_harness_scan_override_is_refused_before_launch(tmp_path: Path, monkeypatch) -> None:
    """A typo must not become a silently auto-selected scan."""
    tree = tmp_path / "tree"
    config = tree / "vntyper" / "config.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({"cram": {"unmapped_scan": "auto"}}))
    case = _case("b178_stream_cram", alignment_kind="cram", cram="/data/x.cram", unmapped_scan="stream")
    del case["bam"]
    monkeypatch.setenv("VNTYPER_CRAM_UNMAPPED_SCAN", "lossy")

    with pytest.raises(ValueError, match="VNTYPER_CRAM_UNMAPPED_SCAN"):
        runner.materialize_case_config(tree, case, tmp_path / "logs")


def test_an_unrecognised_alignment_kind_is_refused_rather_than_defaulted(tmp_path: Path) -> None:
    """Falling back to ``--bam`` would launch a different input from the declared one.

    Args:
        tmp_path: pytest's per-test directory.

    Raises:
        AssertionError: If an unknown ``alignment_kind`` silently launches something.
    """
    case = _case("odd", alignment_kind="sam")
    with pytest.raises(ValueError, match="alignment_kind='sam'"):
        _argv(case, tmp_path)


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


def test_load_side_rejects_a_non_object_record(tmp_path: Path) -> None:
    """A syntactically valid JSON array cannot masquerade as a side attestation."""
    (tmp_path / "side.json").write_text("[]", encoding="utf-8")

    with pytest.raises(ValueError, match="must contain a JSON object"):
        runner.load_side(tmp_path)


def test_run_one_timeout_records_a_failed_expectation(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    case = _case("timeout")
    timeout = subprocess.TimeoutExpired(
        ["golden-cohort-gate"],
        60,
        output=f"{runner.launcher.LAUNCH_PREFIX} side=after\n".encode(),
        stderr=b"partial timeout diagnostic",
    )

    def raise_timeout(*args: object, **kwargs: object) -> None:
        del args, kwargs
        raise timeout

    monkeypatch.setattr(runner.subprocess, "run", raise_timeout)
    record = runner._run_one(
        case=case,
        argv=["pipeline", "--bam", case["bam"]],
        tree=tmp_path / "tree",
        side="after",
        marker="vntyper.scripts.marker",
        expect_marker=True,
        output_dir=tmp_path / "output",
        log_dir=tmp_path / "logs",
        timeout=60,
    )

    assert record["exit_code"] == 99
    assert record["timed_out"] is True
    assert record["expectation_met"] is False
    assert record["expectation_problems"] == [
        "the run was killed on the harness timeout, so its exit code is not a pipeline result"
    ]
    assert (tmp_path / "logs/stdout.log").read_text(encoding="utf-8") == "GATE-LAUNCH side=after\n"
    assert (tmp_path / "logs/stderr.log").read_text(encoding="utf-8") == "partial timeout diagnostic"


def test_run_one_names_every_missing_required_artifact(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    case = _case("missing", required_artifacts=["pipeline_summary.json", "kestrel/kestrel_result.tsv"])
    completed = subprocess.CompletedProcess(
        ["golden-cohort-gate"],
        0,
        stdout=f"{runner.launcher.LAUNCH_PREFIX} side=after\n",
        stderr="",
    )
    monkeypatch.setattr(runner.subprocess, "run", lambda *args, **kwargs: completed)

    record = runner._run_one(
        case=case,
        argv=["pipeline", "--bam", case["bam"]],
        tree=tmp_path / "tree",
        side="after",
        marker="vntyper.scripts.marker",
        expect_marker=True,
        output_dir=tmp_path / "output",
        log_dir=tmp_path / "logs",
        timeout=60,
    )

    assert record["missing_artifacts"] == ["pipeline_summary.json", "kestrel/kestrel_result.tsv"]
    assert record["expectation_problems"] == [
        "exited as expected but did not write 2 of 2 required artefact(s): "
        "pipeline_summary.json, kestrel/kestrel_result.tsv"
    ]
    assert record["expectation_met"] is False


def test_missing_launch_marker_makes_the_side_unverified_even_when_exit_is_expected(tmp_path: Path) -> None:
    case = _case("expected-refusal", expect_exit="nonzero")
    record = {**_ok_record(case["case_id"]), "exit_code": 1, "launch_verified": False}

    side = _run_side(tmp_path, [case], {case["case_id"]: record})

    assert side["expectations_met"] is True
    assert side["launch_verified"] is False


def test_parallel_results_retain_matrix_case_order(tmp_path: Path) -> None:
    cases = [_case("slow"), _case("fast")]

    def fake_run_one(*, case: dict, **_kwargs: object) -> dict:
        if case["case_id"] == "slow":
            time.sleep(0.02)
        return _ok_record(case["case_id"])

    with (
        mock.patch.object(runner, "_run_one", side_effect=fake_run_one),
        mock.patch.object(runner.admissibility, "describe_tree", return_value=CLEAN_REVISION),
    ):
        record = runner.run_side(
            matrix=_matrix(cases),
            tree=tmp_path,
            run_root=tmp_path / "run",
            side="after",
            marker="vntyper.scripts.cohort_rules",
            expect_marker=True,
            threads=4,
            advntr_threads=8,
            jobs=2,
            timeout=60,
            skip_cohort=True,
        )
    assert list(record["pipeline_results"]) == ["slow", "fast"]


def test_cohort_execution_is_skipped_only_when_requested(tmp_path: Path) -> None:
    case = _case("pipeline")
    cohort_record = {**_ok_record("cohort"), "blocked": False}
    cohort_calls: list[bool] = []

    def fake_cohorts(**_kwargs: object) -> dict[str, dict]:
        cohort_calls.append(True)
        return {"cohort": cohort_record}

    records: list[dict] = []
    for skip_cohort in (False, True):
        with (
            mock.patch.object(runner, "_run_one", return_value=_ok_record(case["case_id"])),
            mock.patch.object(runner, "_run_cohorts", side_effect=fake_cohorts),
            mock.patch.object(runner.admissibility, "describe_tree", return_value=CLEAN_REVISION),
        ):
            records.append(
                runner.run_side(
                    matrix=_matrix([case]),
                    tree=tmp_path,
                    run_root=tmp_path / f"run-{skip_cohort}",
                    side="after",
                    marker="vntyper.scripts.cohort_rules",
                    expect_marker=True,
                    threads=4,
                    advntr_threads=8,
                    jobs=1,
                    timeout=60,
                    skip_cohort=skip_cohort,
                )
            )

    assert cohort_calls == [True]
    assert records[0]["cohort_results"] == {"cohort": cohort_record}
    assert records[1]["cohort_results"] == {}


@pytest.mark.parametrize("malformed", [False, True])
def test_load_side_rejects_missing_and_malformed_records(tmp_path: Path, malformed: bool) -> None:
    if malformed:
        (tmp_path / "side.json").write_text("{broken", encoding="utf-8")

    with pytest.raises(
        ValueError,
        match=r"No side\.json under .* That side has not been run, so there is nothing to compare\.",
    ):
        runner.load_side(tmp_path)


def test_run_cohorts_blocks_each_missing_pipeline_summary(tmp_path: Path) -> None:
    cases_root = tmp_path / "cases"
    matrix = {
        "cohort_cases": [
            {
                "case_id": "cohort",
                "inputs": ["sample-a", "sample-b"],
                "required_artifacts": [],
                "expect_exit": "zero",
            }
        ]
    }

    result = runner._run_cohorts(
        matrix=matrix,
        cases_root=cases_root,
        run_root=tmp_path / "run",
        logs_root=tmp_path / "logs",
        tree=tmp_path / "tree",
        side="after",
        marker="vntyper.scripts.marker",
        expect_marker=True,
        timeout=60,
    )["cohort"]

    assert result["blocked"] is True
    assert result["missing_inputs"] == [
        str(cases_root / "sample-a"),
        str(cases_root / "sample-b"),
    ]
    assert result["blocked_reason"].startswith("2 of 2 input directories have no pipeline_summary.json")


def test_run_cohorts_can_measure_an_explicitly_allowed_partial_cohort(tmp_path: Path) -> None:
    cases_root = tmp_path / "cases"
    matrix = {
        "cohort_cases": [
            {
                "case_id": "cohort",
                "inputs": ["sample-a", "sample-b"],
                "allow_missing_inputs": True,
                "required_artifacts": [],
                "expect_exit": "zero",
            }
        ]
    }

    with mock.patch.object(runner, "_run_one", return_value=_ok_record("cohort")):
        result = runner._run_cohorts(
            matrix=matrix,
            cases_root=cases_root,
            run_root=tmp_path / "run",
            logs_root=tmp_path / "logs",
            tree=tmp_path / "tree",
            side="after",
            marker="vntyper.scripts.marker",
            expect_marker=True,
            timeout=60,
        )["cohort"]

    assert result["blocked"] is False
    assert result["missing_inputs"] == [
        str(cases_root / "sample-a"),
        str(cases_root / "sample-b"),
    ]
    assert result["input_count"] == 2


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
