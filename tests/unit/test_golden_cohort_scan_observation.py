"""Focused evidence tests for the golden cohort's executed CRAM scan."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest import mock

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from golden_cohort import cram_evidence, read_set_commands, runner  # noqa: E402

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("scan", ["indexed", "stream"])
def test_runner_records_the_observed_scan_that_decides_a_successful_case(tmp_path: Path, scan: str) -> None:
    """Deleting the runner's command-observer wiring must fail this contract."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "/case/samtools"}}), encoding="utf-8")
    unmapped_bam = output_dir / "fastq_bam_processing" / "output_unmapped.bam"
    unmapped_bam.parent.mkdir(parents=True)
    unmapped_bam.touch()
    expected = {"count": 20, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": f"indexed_safe_{scan}_cram",
        "alignment_kind": "cram",
        "cram": str(tmp_path / "input.cram"),
        "expect_exit": "zero",
        "required_artifacts": [],
        "case_config_path": str(config_path),
        "effective_unmapped_scan": scan,
        "cram_evidence_expectation": {
            "indexed_authorized": True,
            "raw_indexed_read_set": expected,
            "stream_read_set": expected,
        },
    }
    command = f"observed {scan} extraction command"
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n",
        stderr="",
        returncode=0,
    )

    with (
        mock.patch.object(runner.subprocess, "run", return_value=pipeline_result),
        mock.patch.object(
            cram_evidence,
            "observe_unmapped_scan",
            return_value=(scan, command, []),
        ) as observe,
        mock.patch.object(read_set_commands, "collect_read_set_evidence", return_value=expected),
    ):
        record = runner._run_one(
            case=case,
            argv=["pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.alignment_contract",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    observe.assert_called_once_with(log_dir / "commands.jsonl", unmapped_bam)
    assert record["observed_unmapped_scan"] == scan
    assert record["observed_unmapped_command"] == command
    assert record["scan_observation_problems"] == []
    assert record["expectation_met"] is True
    assert record["expectation_problems"] == []


@pytest.mark.parametrize("scan", ["indexed", "stream"])
def test_authorized_scan_decisions_require_the_same_nonempty_read_set(scan: str) -> None:
    """Both genuine extraction paths must satisfy the existing A-178-2 oracle."""
    expected = {"count": 20, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": f"indexed_safe_{scan}_cram",
        "effective_unmapped_scan": scan,
        "cram_evidence_expectation": {
            "indexed_authorized": True,
            "raw_indexed_read_set": expected,
            "stream_read_set": expected,
        },
    }
    record = {
        "observed_unmapped_scan": scan,
        "observed_unmapped_command": f"observed {scan} extraction command",
        "raw_indexed_read_set": expected,
        "unmapped_read_set": expected,
        "raw_indexed_loss": 0,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == []


def test_forced_indexed_decision_rejects_an_observed_stream_command() -> None:
    """A stream extraction cannot satisfy the required forced-indexed rejection."""
    raw = {"count": 10, "sorted_read_name_sha256": "a" * 64}
    stream = {"count": 14, "sorted_read_name_sha256": "b" * 64}
    case = {
        "case_id": "b178_indexed_cram",
        "effective_unmapped_scan": "indexed",
        "cram_evidence_expectation": {
            "placed_unmapped_guard_count": 329,
            "raw_indexed_read_set": raw,
            "stream_read_set": stream,
        },
    }
    record = {
        "exit_code": 1,
        "observed_unmapped_scan": "stream",
        "observed_unmapped_command": "set -o pipefail; samtools view ... | samtools view -f 4 ...",
        "raw_indexed_read_set": raw,
        "unmapped_read_set": stream,
        "raw_indexed_loss": 4,
        "placed_unmapped_guard_count": 329,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        "A-178-2 indexed rejection observed a stream unmapped-extraction command instead of failing before work",
        "A-178-2 indexed rejection executed an unmapped-extraction command instead of failing before work",
        "A-178-2 indexed rejection produced an unmapped BAM instead of failing before work",
    ]


@pytest.mark.parametrize(("declared", "observed"), [("indexed", "stream"), ("stream", "indexed")])
def test_authorized_scan_decision_rejects_an_observed_mode_different_from_its_declaration(
    declared: str,
    observed: str,
) -> None:
    """The two authorized arms are distinct measurements, not interchangeable labels."""
    expected = {"count": 20, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": f"indexed_safe_{declared}_cram",
        "effective_unmapped_scan": declared,
        "cram_evidence_expectation": {
            "indexed_authorized": True,
            "raw_indexed_read_set": expected,
            "stream_read_set": expected,
        },
    }
    record = {
        "observed_unmapped_scan": observed,
        "observed_unmapped_command": f"observed {observed} extraction command",
        "raw_indexed_read_set": expected,
        "unmapped_read_set": expected,
        "raw_indexed_loss": 0,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        f"A-178-2 observed {observed} extraction, but declared {declared}"
    ]


def test_runner_rejects_an_authorized_case_when_the_observed_mode_differs_from_its_declaration(tmp_path: Path) -> None:
    """The observer's mismatch must reach the runner's admissibility result."""
    tree = tmp_path / "tree"
    output_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    config_path = log_dir / "config.json"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(json.dumps({"tools": {"samtools": "/case/samtools"}}), encoding="utf-8")
    unmapped_bam = output_dir / "fastq_bam_processing" / "output_unmapped.bam"
    unmapped_bam.parent.mkdir(parents=True)
    unmapped_bam.touch()
    expected = {"count": 20, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": "indexed_safe_indexed_cram",
        "alignment_kind": "cram",
        "cram": str(tmp_path / "input.cram"),
        "expect_exit": "zero",
        "required_artifacts": [],
        "case_config_path": str(config_path),
        "effective_unmapped_scan": "indexed",
        "cram_evidence_expectation": {
            "indexed_authorized": True,
            "raw_indexed_read_set": expected,
            "stream_read_set": expected,
        },
    }
    pipeline_result = mock.Mock(
        stdout=f"{runner.launcher.LAUNCH_PREFIX} verified\n",
        stderr="",
        returncode=0,
    )

    with (
        mock.patch.object(runner.subprocess, "run", return_value=pipeline_result),
        mock.patch.object(
            cram_evidence,
            "observe_unmapped_scan",
            return_value=("stream", "observed stream extraction command", []),
        ),
        mock.patch.object(read_set_commands, "collect_read_set_evidence", return_value=expected),
    ):
        record = runner._run_one(
            case=case,
            argv=["pipeline", "--cram", str(case["cram"])],
            tree=tree,
            side="after",
            marker="vntyper.scripts.alignment_contract",
            expect_marker=True,
            output_dir=output_dir,
            log_dir=log_dir,
            timeout=60,
        )

    assert record["expectation_met"] is False
    assert record["expectation_problems"] == ["A-178-2 observed stream extraction, but declared indexed"]


def test_decision_fails_closed_without_an_observed_extraction_command() -> None:
    """A declared stream mode alone is not evidence that its command executed."""
    expected = {"count": 20, "sorted_read_name_sha256": "a" * 64}
    case = {
        "case_id": "indexed_safe_stream_cram",
        "effective_unmapped_scan": "stream",
        "cram_evidence_expectation": {
            "indexed_authorized": True,
            "raw_indexed_read_set": expected,
            "stream_read_set": expected,
        },
    }
    record = {
        "observed_unmapped_scan": None,
        "observed_unmapped_command": None,
        "raw_indexed_read_set": expected,
        "unmapped_read_set": expected,
        "raw_indexed_loss": 0,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        "A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"
    ]


@pytest.mark.parametrize(
    ("command", "expected_mode"),
    [
        (
            "samtools view -@ 4 -b -f 4 -T reference.fa input.cram '*' "
            "-o /run/fastq_bam_processing/output_unmapped.bam",
            "indexed",
        ),
        (
            "set -o pipefail; samtools view -@ 4 -h -T reference.fa input.cram | "
            "samtools view -b -f 4 -@ 4 - -o /run/fastq_bam_processing/output_unmapped.bam",
            "stream",
        ),
        # The single-process stream shape (#262): no shell pipe, no `set`, no '*'.
        # Keyed on the absence of the index query, which is the property that actually
        # distinguishes a whole-file decode from an indexed fetch.
        (
            "samtools view -b -f 4 -u -T reference.fa -@ 4 input.cram -o /run/fastq_bam_processing/output_unmapped.bam",
            "stream",
        ),
        # Without a reference, and without a thread flag at one thread.
        (
            "samtools view -b -f 4 -u input.bam -o /run/fastq_bam_processing/output_unmapped.bam",
            "stream",
        ),
        # The -X custom-index form of the indexed fetch still classifies as indexed.
        (
            "samtools view -b -f 4 -u -@ 4 -X input.cram /proc/9/fd/3 '*' "
            "-o /run/fastq_bam_processing/output_unmapped.bam",
            "indexed",
        ),
    ],
)
def test_command_log_records_the_executed_scan(tmp_path: Path, command: str, expected_mode: str) -> None:
    """Only the subprocess command log can attest the extraction strategy."""
    commands_log = tmp_path / "commands.jsonl"
    commands_log.write_text(json.dumps({"command": command, "shell": True}) + "\n", encoding="utf-8")

    mode, observed_command, problems = cram_evidence.observe_unmapped_scan(
        commands_log,
        Path("/run/fastq_bam_processing/output_unmapped.bam"),
    )

    assert mode == expected_mode
    assert observed_command == command
    assert problems == []


def test_command_log_rejects_contradictory_extraction_modes(tmp_path: Path) -> None:
    """Both strategies in one case are ambiguous evidence, never a passing mode."""
    target = "/run/fastq_bam_processing/output_unmapped.bam"
    commands = [
        f"samtools view -b -f 4 input.cram '*' -o {target}",
        f"set -o pipefail; samtools view -h input.cram | samtools view -b -f 4 - -o {target}",
    ]
    commands_log = tmp_path / "commands.jsonl"
    commands_log.write_text(
        "".join(json.dumps({"command": command, "shell": True}) + "\n" for command in commands),
        encoding="utf-8",
    )

    mode, observed_command, problems = cram_evidence.observe_unmapped_scan(commands_log, Path(target))

    assert mode is None
    assert observed_command is None
    assert problems == ["A-178-2 observed contradictory CRAM unmapped-extraction modes: indexed, stream"]


def test_command_log_rejects_an_unknown_command_that_references_the_unmapped_bam(tmp_path: Path) -> None:
    """An unrecognized extraction mutation is not equivalent to proof that no extraction ran."""
    target = "/run/fastq_bam_processing/output_unmapped.bam"
    commands_log = tmp_path / "commands.jsonl"
    commands_log.write_text(
        json.dumps({"command": f"samtools view --write-index input.cram >{target}", "shell": True}) + "\n",
        encoding="utf-8",
    )

    mode, observed_command, problems = cram_evidence.observe_unmapped_scan(commands_log, Path(target))

    assert mode is None
    assert observed_command is None
    assert problems == [
        "A-178-2 command log references the unmapped BAM but no recognized extraction mode was observed"
    ]


def test_a_command_that_writes_the_target_without_filtering_is_still_unrecognised(tmp_path: Path) -> None:
    """Widening the stream branch must not turn every write into a recognised scan.

    Keying `stream` on "not indexed" is only safe because the `-f 4` and `-o <target>`
    guards sit above it. Dropping either would make an arbitrary command that touches
    the unmapped BAM look like a proven whole-file extraction.
    """
    target = "/run/fastq_bam_processing/output_unmapped.bam"
    commands_log = tmp_path / "commands.jsonl"
    commands_log.write_text(
        json.dumps({"command": f"samtools view -b -u input.cram -o {target}", "shell": True}) + "\n",
        encoding="utf-8",
    )

    mode, observed_command, problems = cram_evidence.observe_unmapped_scan(commands_log, Path(target))

    assert mode is None
    assert observed_command is None
    assert problems == [
        "A-178-2 command log references the unmapped BAM but no recognized extraction mode was observed"
    ]


def test_a_flag_four_command_writing_somewhere_else_is_not_observed(tmp_path: Path) -> None:
    """The target guard is what keeps an unrelated flag-4 command out of the evidence."""
    target = "/run/fastq_bam_processing/output_unmapped.bam"
    commands_log = tmp_path / "commands.jsonl"
    commands_log.write_text(
        json.dumps({"command": "samtools view -b -f 4 -u input.cram -o /somewhere/else.bam", "shell": True}) + "\n",
        encoding="utf-8",
    )

    mode, observed_command, problems = cram_evidence.observe_unmapped_scan(commands_log, Path(target))

    assert mode is None
    assert problems == ["A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode"]


@pytest.mark.parametrize(
    ("content", "expected_problem"),
    [
        (
            json.dumps({"command": "samtools idxstats input.cram", "shell": False}) + "\n",
            "A-178-2 did not observe exactly one executed CRAM unmapped-extraction mode",
        ),
        ("not-json\n", "A-178-2 command log is malformed at line 1"),
    ],
)
def test_command_log_fails_closed_on_missing_or_malformed_evidence(
    tmp_path: Path,
    content: str,
    expected_problem: str,
) -> None:
    """A missing or corrupt recorder entry cannot attest a scan mode."""
    commands_log = tmp_path / "commands.jsonl"
    commands_log.write_text(content, encoding="utf-8")

    mode, observed_command, problems = cram_evidence.observe_unmapped_scan(
        commands_log,
        Path("/run/fastq_bam_processing/output_unmapped.bam"),
    )

    assert mode is None
    assert observed_command is None
    assert problems == [expected_problem]


def test_authorized_indexed_decision_fails_closed_on_a_different_produced_read_set() -> None:
    """Authorization never permits a lossy indexed result to pass the gate."""
    expected = {"count": 20, "sorted_read_name_sha256": "a" * 64}
    different = {"count": 19, "sorted_read_name_sha256": "b" * 64}
    case = {
        "case_id": "indexed_safe_indexed_cram",
        "effective_unmapped_scan": "indexed",
        "cram_evidence_expectation": {
            "indexed_authorized": True,
            "raw_indexed_read_set": expected,
            "stream_read_set": expected,
        },
    }
    record = {
        "observed_unmapped_scan": "indexed",
        "observed_unmapped_command": "samtools view -f 4 input.cram '*' -o output.bam",
        "raw_indexed_read_set": expected,
        "unmapped_read_set": different,
        "raw_indexed_loss": 1,
    }

    assert cram_evidence.validate_cram_evidence(case, record) == [
        f"A-178-2 stream evidence differs: expected {expected}, got {different}",
        "A-178-2 raw indexed loss differs: expected 0, got 1",
    ]
