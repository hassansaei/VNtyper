"""Unit tests for A-178-2 read-set evidence collection and runner wiring."""

import io
import json
import sys
from pathlib import Path
from unittest import mock

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from golden_cohort import read_set_commands, runner  # noqa: E402
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


def test_a_successful_cram_case_fails_closed_when_read_set_evidence_cannot_be_collected(tmp_path: Path) -> None:
    """A zero-exit CRAM run cannot be admissible when its samtools evidence command fails."""
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
