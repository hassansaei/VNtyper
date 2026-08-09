"""Unit tests for A-178-2 read-set evidence collection and runner wiring."""

import io
import json
import sys
from pathlib import Path
from typing import cast
from unittest import mock

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from golden_cohort import cram_evidence, read_set_commands, runner  # noqa: E402
from golden_cohort.read_sets import summarize_sorted_read_names, summarize_unmapped_read_set  # noqa: E402


def _case(case_id: str, config_path: Path) -> dict[str, object]:
    """Build one successful non-fast CRAM case for the runner seam."""
    return {
        "case_id": case_id,
        "kind": "pipeline",
        "group": "cram",
        "alignment_kind": "cram",
        "cram": "/data/x.cram",
        "assembly": "hg19",
        "fast_mode": False,
        "advntr": False,
        "expect_exit": "zero",
        "required_artifacts": [],
        "case_config_path": str(config_path),
    }


def test_unmapped_read_set_evidence_counts_records_and_hashes_sorted_names() -> None:
    """The evidence must be order-independent while retaining duplicate QNAMEs."""
    records = "read-b\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t!!\nread-a\t4\t*\t0\t0\t*\t*\t0\t0\tTG\t!!\nread-b\t4\t*\t0\t0\t*\t*\t0\t0\tGT\t!!\n"

    evidence = summarize_unmapped_read_set("3\n", records)

    assert evidence == {
        "count": 3,
        "sorted_read_name_sha256": "87d9644834b735d999f234c38234d23d78df482bc818a1d29e7a1a1d73696e1a",
    }


def test_unmapped_read_set_evidence_is_identical_for_a_different_record_order() -> None:
    """Indexed and stream output ordering must not change the comparison digest."""
    first = "read-b\t4\nread-a\t4\n"
    second = "read-a\t4\nread-b\t4\n"

    assert summarize_unmapped_read_set("2", first) == summarize_unmapped_read_set("2", second)


@pytest.mark.parametrize("count_output", ["", "not-a-count", "-1", "1 2"])
def test_unmapped_read_set_evidence_rejects_a_malformed_count(count_output: str) -> None:
    """A broken ``samtools view -c`` result must not become plausible evidence."""
    with pytest.raises(ValueError, match="samtools view -c"):
        summarize_unmapped_read_set(count_output, "")


def test_unmapped_read_set_evidence_rejects_count_and_record_disagreement() -> None:
    """The separately counted and digested record streams must describe the same reads."""
    with pytest.raises(ValueError, match="reported 2 records but emitted 1"):
        summarize_unmapped_read_set("2", "read-a\t4\n")


def test_unmapped_read_set_evidence_rejects_a_record_without_a_qname_field() -> None:
    """Malformed view output must fail instead of hashing an unverifiable string."""
    with pytest.raises(ValueError, match="malformed SAM record"):
        summarize_unmapped_read_set("1", "read-a\n")


def test_sorted_read_name_evidence_accepts_a_zero_record_stream() -> None:
    """An empty unmapped BAM has a stable count and the SHA-256 of empty input."""
    assert summarize_sorted_read_names("0\n", iter(())) == {
        "count": 0,
        "sorted_read_name_sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    }


def test_sorted_read_name_evidence_rejects_an_empty_qname_field() -> None:
    """A tab-only line cannot be accepted as a read name by the streaming seam."""
    with pytest.raises(ValueError, match="invalid QNAME"):
        summarize_sorted_read_names("1", iter(["\t\n"]))


def test_a_successful_cram_case_records_unmapped_read_set_evidence(tmp_path: Path) -> None:
    """The runner must persist A-178-2 evidence beside each successful CRAM result."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "/opt/samtools"}}), encoding="utf-8")
    unmapped_bam = output_dir / "fastq_bam_processing" / "output_unmapped.bam"
    unmapped_bam.parent.mkdir(parents=True)
    unmapped_bam.touch()
    case = _case("b178_indexed_cram", config_path)
    pipeline_result = mock.Mock(stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n", stderr="", returncode=0)
    count_result = mock.Mock(stdout="1\n", stderr="", returncode=0)
    view_process = mock.Mock(stdout=io.StringIO("read-a\t4\n"))
    view_process.wait.return_value = 0

    def run_evidence_command(command: list[str], **_kwargs: object) -> mock.Mock:
        if command[0] == sys.executable:
            return pipeline_result
        if command[:3] == ["/opt/samtools", "view", "-c"]:
            return count_result
        if command[0] == "sort":
            sorted_path = Path(command[command.index("-o") + 1])
            names_path = Path(command[-1])
            sorted_path.write_text("".join(sorted(names_path.read_text().splitlines(keepends=True))))
            return mock.Mock(stdout="", stderr="", returncode=0)
        raise AssertionError(f"unexpected command: {command}")

    with (
        mock.patch.object(read_set_commands.subprocess, "run", side_effect=run_evidence_command) as evidence_run,
        mock.patch.object(read_set_commands.subprocess, "Popen", return_value=view_process) as popen,
    ):
        record = runner._run_one(
            case=case,
            argv=["--config-path", str(config_path), "pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.cram_preflight",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    assert record["unmapped_read_set"] == {
        "count": 1,
        "sorted_read_name_sha256": "dc9ecb288105192549da23a9bd863fa8bbcf34bf47ddc3e7728031cd14afb044",
    }
    assert evidence_run.call_args_list[1].args[0] == ["/opt/samtools", "view", "-c", str(unmapped_bam)]
    assert popen.call_args.args[0] == ["/opt/samtools", "view", str(unmapped_bam)]
    assert evidence_run.call_args_list[2].args[0][0:2] == ["sort", "-o"]
    assert evidence_run.call_args_list[2].kwargs["env"]["LC_ALL"] == "C"
    assert json.loads((log_dir / "result.json").read_text())["unmapped_read_set"] == record["unmapped_read_set"]


def test_a_nonzero_cram_case_with_an_unmapped_bam_still_records_read_set_evidence(tmp_path: Path) -> None:
    """A later fail-closed routing exit must not erase the extraction evidence already written."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "/case/samtools"}}), encoding="utf-8")
    unmapped_bam = output_dir / "fastq_bam_processing" / "output_unmapped.bam"
    unmapped_bam.parent.mkdir(parents=True)
    unmapped_bam.touch()
    case = _case("7a61_stream_cram", config_path)
    case["expect_exit"] = "nonzero"
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n", stderr="mixed layout", returncode=1
    )
    evidence = {
        "count": 622_690,
        "sorted_read_name_sha256": "a" * 64,
    }

    with (
        mock.patch.object(read_set_commands, "collect_read_set_evidence", return_value=evidence) as collect,
        mock.patch.object(runner.subprocess, "run", return_value=pipeline_result),
        mock.patch.object(runner.time, "monotonic", side_effect=[100.0, 102.345]),
    ):
        record = runner._run_one(
            case=case,
            argv=["--config-path", str(config_path), "pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.cram_preflight",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    collect.assert_called_once_with(
        unmapped_bam,
        "/case/samtools",
        cwd=tree,
        temporary_parent=log_dir,
        timeout=60,
    )
    assert record["unmapped_read_set"] == evidence
    assert record["exit_code"] == 1
    assert record["seconds"] == 2.34
    assert record["expectation_met"] is True
    assert record["expectation_problems"] == []


def test_rejected_indexed_cram_records_raw_star_evidence_without_bypassing_the_guard(tmp_path: Path) -> None:
    """A-178-2 requires the rejected strategy's would-be read set and no produced BAM."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "/case/samtools"}}), encoding="utf-8")
    case = _case("7a61_indexed_cram", config_path)
    case.update(
        {
            "expect_exit": "nonzero",
            "effective_unmapped_scan": "indexed",
            "cram_evidence_expectation": {
                "placed_unmapped_guard_count": 11_571,
                "raw_indexed_read_set": {
                    "count": 2_690,
                    "sorted_read_name_sha256": "c64be739cd6344b8b62004fc9ea568779b3cc06ff1d472ac0e5d97c343130d7d",
                },
                "stream_read_set": {
                    "count": 622_690,
                    "sorted_read_name_sha256": "6677ba2931466a1519bf0f8b9527d783e2a73049462fe9a72cf9326c83cb8b91",
                },
            },
        }
    )
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n",
        stderr="idxstats reports 11571 placed-unmapped reads; using stream scan",
        returncode=1,
    )
    expectation = cast(dict[str, object], case["cram_evidence_expectation"])
    raw_evidence = cast(dict[str, object], expectation["raw_indexed_read_set"])

    with (
        mock.patch.object(read_set_commands, "collect_read_set_evidence", return_value=raw_evidence) as collect,
        mock.patch.object(runner.subprocess, "run", return_value=pipeline_result),
    ):
        record = runner._run_one(
            case=case,
            argv=["--config-path", str(config_path), "pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.cram_preflight",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    collect.assert_called_once_with(
        Path(str(case["cram"])),
        "/case/samtools",
        cwd=tree,
        temporary_parent=log_dir,
        timeout=60,
        regions=("*",),
    )
    assert record["unmapped_read_set"] is None
    assert record["raw_indexed_read_set"] == raw_evidence
    assert record["raw_indexed_loss"] is None
    assert record["placed_unmapped_guard_count"] == 11_571
    assert record["expectation_met"] is True
    assert record["expectation_problems"] == []


@pytest.mark.parametrize(
    ("stderr", "expected"),
    [
        ("idxstats reports 329 placed-unmapped reads; using stream scan", 329),
        (
            "2026-08-08 23:46:51,964 - vntyper.scripts.idxstats_parsing - ERROR - "
            "idxstats reports 11571 placed-unmapped reads; using stream scan",
            11_571,
        ),
        (
            "ValueError: idxstats reports 329 placed-unmapped reads; using stream scan",
            329,
        ),
    ],
)
def test_placed_unmapped_guard_parser_accepts_only_the_stable_diagnostic(stderr: str, expected: int) -> None:
    """The causal indexed rejection is recovered from the pipeline's stable diagnostic."""
    assert cram_evidence.parse_placed_unmapped_guard_count(stderr) == expected


@pytest.mark.parametrize(
    "stderr",
    [
        "",
        "missing CRAM reference",
        "idxstats reports many placed-unmapped reads; using stream scan",
        "idxstats reports -1 placed-unmapped reads; using stream scan",
        "idxstats reports 329 placed unmapped reads; using stream scan",
        "user message embeds idxstats reports 329 placed-unmapped reads; using stream scan here",
        "/tmp/idxstats reports 329 placed-unmapped reads; using stream scan.log",
        "WARNING - idxstats reports 329 placed-unmapped reads; using stream scan",
        (
            "idxstats reports 329 placed-unmapped reads; using stream scan\n"
            "idxstats reports 631571 placed-unmapped reads; using stream scan"
        ),
    ],
)
def test_placed_unmapped_guard_parser_rejects_missing_malformed_unrelated_or_conflicting_diagnostics(
    stderr: str,
) -> None:
    """Near matches and conflicting counts cannot authorize the golden rejection."""
    assert cram_evidence.parse_placed_unmapped_guard_count(stderr) is None


def test_placed_unmapped_guard_parser_rejects_a_hostile_huge_count_without_raising() -> None:
    """Digit length is bounded before conversion so hostile stderr fails closed cheaply."""
    stderr = f"ValueError: idxstats reports {'9' * 100_000} placed-unmapped reads; using stream scan"

    assert cram_evidence.parse_placed_unmapped_guard_count(stderr) is None


@pytest.mark.parametrize(("actual", "expected_problem_value"), [(None, "None"), (328, "328")])
def test_forced_indexed_decision_rejects_missing_or_wrong_guard_count(
    actual: int | None,
    expected_problem_value: str,
) -> None:
    """An unrelated exit 1 or the wrong guard count cannot satisfy A-178-2."""
    raw = {"count": 4_478, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": "b178_indexed_cram",
        "effective_unmapped_scan": "indexed",
        "cram_evidence_expectation": {
            "placed_unmapped_guard_count": 329,
            "raw_indexed_read_set": raw,
            "stream_read_set": {"count": 4_807, "sorted_read_name_sha256": "b" * 64},
        },
    }
    record = {
        "exit_code": 1,
        "raw_indexed_read_set": raw,
        "unmapped_read_set": None,
        "raw_indexed_loss": None,
        "placed_unmapped_guard_count": actual,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        f"A-178-2 indexed guard count differs: expected 329, got {expected_problem_value}"
    ]


def test_forced_indexed_decision_rejects_a_boolean_guard_declaration() -> None:
    """Boolean is an ``int`` subclass but cannot be accepted as a measured read count."""
    raw = {"count": 1, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": "malformed_indexed_cram",
        "effective_unmapped_scan": "indexed",
        "cram_evidence_expectation": {
            "placed_unmapped_guard_count": True,
            "raw_indexed_read_set": raw,
            "stream_read_set": raw,
        },
    }
    record = {
        "exit_code": 1,
        "raw_indexed_read_set": raw,
        "unmapped_read_set": None,
        "raw_indexed_loss": None,
        "placed_unmapped_guard_count": 1,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        "A-178-2 indexed guard expectation is missing or malformed: True"
    ]


@pytest.mark.parametrize("exit_code", [0, 2])
def test_forced_indexed_decision_requires_the_guard_to_cause_exit_one(exit_code: int) -> None:
    """The stable guard is evidence only when it produces the declared refusal exit."""
    raw = {"count": 4_478, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": "b178_indexed_cram",
        "effective_unmapped_scan": "indexed",
        "cram_evidence_expectation": {
            "placed_unmapped_guard_count": 329,
            "raw_indexed_read_set": raw,
            "stream_read_set": {"count": 4_807, "sorted_read_name_sha256": "b" * 64},
        },
    }
    record = {
        "exit_code": exit_code,
        "raw_indexed_read_set": raw,
        "unmapped_read_set": None,
        "raw_indexed_loss": None,
        "placed_unmapped_guard_count": 329,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        f"A-178-2 indexed rejection exit differs: expected 1, got {exit_code}"
    ]


def test_forced_indexed_runner_rejects_an_unrelated_preflight_exit(tmp_path: Path) -> None:
    """Exit 1 with no output is insufficient without the measured guard diagnostic."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "/case/samtools"}}), encoding="utf-8")
    case = _case("b178_indexed_cram", config_path)
    raw = {"count": 4_478, "sorted_read_name_sha256": "a" * 64}
    case.update(
        {
            "expect_exit": "nonzero",
            "effective_unmapped_scan": "indexed",
            "cram_evidence_expectation": {
                "placed_unmapped_guard_count": 329,
                "raw_indexed_read_set": raw,
                "stream_read_set": {"count": 4_807, "sorted_read_name_sha256": "b" * 64},
            },
        }
    )
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n",
        stderr="CRAM reference probe failed before idxstats decision",
        returncode=1,
    )

    with (
        mock.patch.object(read_set_commands, "collect_read_set_evidence", return_value=raw),
        mock.patch.object(runner.subprocess, "run", return_value=pipeline_result),
    ):
        record = runner._run_one(
            case=case,
            argv=["--config-path", str(config_path), "pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.cram_preflight",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    assert record["exit_code"] == 1
    assert record["unmapped_read_set"] is None
    assert record["placed_unmapped_guard_count"] is None
    assert record["expectation_met"] is False
    assert record["expectation_problems"] == ["A-178-2 indexed guard count differs: expected 329, got None"]


def test_stream_cram_decision_records_exact_read_set_and_raw_loss() -> None:
    """The stream decision fails closed unless both measurements match the declared evidence."""
    raw = {
        "count": 2_690,
        "sorted_read_name_sha256": "c64be739cd6344b8b62004fc9ea568779b3cc06ff1d472ac0e5d97c343130d7d",
    }
    stream = {
        "count": 634_261,
        "sorted_read_name_sha256": "b7f75d19497698f12d6dbbc38afc12702b2d262670a4c893b39f95967ebf7b7b",
    }
    case = {
        "case_id": "7a61_stream_cram",
        "effective_unmapped_scan": "stream",
        "cram_evidence_expectation": {"raw_indexed_read_set": raw, "stream_read_set": stream},
    }
    record = {
        "observed_unmapped_scan": "stream",
        "observed_unmapped_command": "set -o pipefail; samtools view ... | samtools view -f 4 ...",
        "raw_indexed_read_set": raw,
        "unmapped_read_set": stream,
        "raw_indexed_loss": 631_571,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == []


def test_stream_cram_decision_fails_when_raw_loss_evidence_is_missing() -> None:
    """A successful exit expectation cannot hide an unmeasured A-178-2 decision."""
    stream = {"count": 4_478, "sorted_read_name_sha256": "d" * 64}
    case = {
        "case_id": "b178_stream_cram",
        "effective_unmapped_scan": "stream",
        "cram_evidence_expectation": {"raw_indexed_read_set": stream, "stream_read_set": stream},
    }
    record = {
        "observed_unmapped_scan": "stream",
        "observed_unmapped_command": "set -o pipefail; samtools view ... | samtools view -f 4 ...",
        "raw_indexed_read_set": stream,
        "unmapped_read_set": stream,
        "raw_indexed_loss": None,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        "A-178-2 stream evidence did not record the raw indexed loss"
    ]


def test_stream_cram_decision_rejects_a_raw_loss_inconsistent_with_the_read_sets() -> None:
    """The recorded loss is a decision input, not an unchecked decorative number."""
    raw = {"count": 10, "sorted_read_name_sha256": "a" * 64}
    stream = {"count": 14, "sorted_read_name_sha256": "b" * 64}
    case = {
        "case_id": "x_stream_cram",
        "effective_unmapped_scan": "stream",
        "cram_evidence_expectation": {"raw_indexed_read_set": raw, "stream_read_set": stream},
    }
    record = {
        "observed_unmapped_scan": "stream",
        "observed_unmapped_command": "set -o pipefail; samtools view ... | samtools view -f 4 ...",
        "raw_indexed_read_set": raw,
        "unmapped_read_set": stream,
        "raw_indexed_loss": 3,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == ["A-178-2 raw indexed loss differs: expected 4, got 3"]


def test_raw_star_evidence_uses_the_same_bounded_collector_with_a_region(tmp_path: Path) -> None:
    """The raw diagnostic must select ``*`` without restoring full-SAM buffering."""
    alignment = tmp_path / "input.cram"
    alignment.touch()
    temporary_parent = tmp_path / "logs"
    temporary_parent.mkdir()
    count_result = mock.Mock(stdout="1\n", stderr="", returncode=0)
    view_process = mock.Mock(stdout=io.StringIO("read-a\t4\n"))
    view_process.wait.return_value = 0

    def run_command(command: list[str], **_kwargs: object) -> mock.Mock:
        if command[:3] == ["samtools", "view", "-c"]:
            return count_result
        sorted_path = Path(command[command.index("-o") + 1])
        names_path = Path(command[-1])
        sorted_path.write_text(names_path.read_text())
        return mock.Mock(stdout="", stderr="", returncode=0)

    with (
        mock.patch.object(read_set_commands.subprocess, "run", side_effect=run_command) as run,
        mock.patch.object(read_set_commands.subprocess, "Popen", return_value=view_process) as popen,
    ):
        evidence = read_set_commands.collect_read_set_evidence(
            alignment,
            "samtools",
            cwd=tmp_path,
            temporary_parent=temporary_parent,
            timeout=30,
            regions=("*",),
        )

    assert evidence["count"] == 1
    assert run.call_args_list[0].args[0] == ["samtools", "view", "-c", str(alignment), "*"]
    assert popen.call_args.args[0] == ["samtools", "view", str(alignment), "*"]


def test_a_cram_preflight_failure_without_an_unmapped_bam_preserves_the_original_expectation(tmp_path: Path) -> None:
    """No extraction artefact means no evidence diagnostic may mask the pipeline's exit failure."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "/case/samtools"}}), encoding="utf-8")
    case = _case("cram_preflight_failure", config_path)
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n", stderr="missing reference", returncode=1
    )

    with (
        mock.patch.object(read_set_commands, "collect_read_set_evidence") as collect,
        mock.patch.object(runner.subprocess, "run", return_value=pipeline_result),
        mock.patch.object(runner.time, "monotonic", side_effect=[200.0, 203.456]),
    ):
        record = runner._run_one(
            case=case,
            argv=["--config-path", str(config_path), "pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.cram_preflight",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    collect.assert_not_called()
    assert record["unmapped_read_set"] is None
    assert record["seconds"] == 3.46
    assert record["expectation_met"] is False
    assert record["expectation_problems"] == ["expected exit 0, got 1"]


def test_an_expected_nonzero_cram_case_fails_closed_when_its_existing_bam_cannot_be_evidenced(tmp_path: Path) -> None:
    """Evidence failure must make an otherwise expected routing exit inadmissible."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "samtools"}}), encoding="utf-8")
    unmapped_bam = output_dir / "fastq_bam_processing" / "output_unmapped.bam"
    unmapped_bam.parent.mkdir(parents=True)
    unmapped_bam.touch()
    case = _case("b178_stream_cram", config_path)
    case["expect_exit"] = "nonzero"
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n", stderr="mixed layout", returncode=1
    )
    failed_count = mock.Mock(stdout="", stderr="cannot decode", returncode=7)

    with mock.patch.object(read_set_commands.subprocess, "run", side_effect=[pipeline_result, failed_count]):
        record = runner._run_one(
            case=case,
            argv=["--config-path", str(config_path), "pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.cram_preflight",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    assert record["unmapped_read_set"] is None
    assert record["exit_code"] == 1
    assert record["exit_problem"] == ""
    assert record["expectation_met"] is False
    assert record["expectation_problems"] == [
        "could not collect CRAM read-set evidence: samtools read-set command exited 7: cannot decode"
    ]


def test_a_stalled_read_set_view_is_killed_reaped_and_its_temporary_files_are_removed(tmp_path: Path) -> None:
    """The evidence watchdog must bound streaming and leave neither a child nor QNAME files."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "samtools"}}), encoding="utf-8")
    unmapped_bam = output_dir / "fastq_bam_processing" / "output_unmapped.bam"
    unmapped_bam.parent.mkdir(parents=True)
    unmapped_bam.touch()
    case = _case("b178_stream_cram", config_path)
    pipeline_result = mock.Mock(stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n", stderr="", returncode=0)
    count_result = mock.Mock(stdout="1\n", stderr="", returncode=0)
    view_process = mock.Mock(stdout=io.StringIO("read-a\t4\n"))
    view_process.wait.return_value = 0
    view_future = mock.Mock()
    view_future.result.side_effect = read_set_commands.FutureTimeoutError
    executor = mock.Mock()
    executor.submit.return_value = view_future

    with (
        mock.patch.object(read_set_commands.subprocess, "run", side_effect=[pipeline_result, count_result]),
        mock.patch.object(read_set_commands.subprocess, "Popen", return_value=view_process),
        mock.patch.object(read_set_commands, "ThreadPoolExecutor", return_value=executor),
    ):
        record = runner._run_one(
            case=case,
            argv=["--config-path", str(config_path), "pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.cram_preflight",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=1,
        )

    view_process.kill.assert_called_once_with()
    view_process.wait.assert_called_once_with()
    view_future.result.assert_called_once_with(timeout=1)
    executor.shutdown.assert_called_once_with(wait=True, cancel_futures=True)
    assert record["expectation_met"] is False
    assert "timed out after 1 seconds" in record["expectation_problems"][0]
    assert list(log_dir.glob("vntyper-read-set-*")) == []
