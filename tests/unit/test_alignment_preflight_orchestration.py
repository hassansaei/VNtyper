"""Orchestration and shipped-policy tests for alignment preflight."""

from __future__ import annotations

import json
import shlex
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_preflight import choose_unmapped_scan, run_preflight
from vntyper.scripts.alignment_preflight_logs import preflight_log_paths

pytestmark = pytest.mark.unit


def test_fast_cram_preflight_logs_include_only_the_probes_that_can_run(tmp_path: Path) -> None:
    """Destination validation must mirror fast mode's target and coverage probes exactly."""
    paths = preflight_log_paths(
        tmp_path,
        "sample",
        "cram",
        candidate_count=2,
        fast_mode=True,
        coverage_region="chr2:30-40",
    )

    assert paths == (
        tmp_path / "sample_reference_probe_1.log",
        tmp_path / "sample_reference_probe_2.log",
        tmp_path / "sample_reference_probe_3.log",
        tmp_path / "sample_reference_coverage_probe_1.log",
        tmp_path / "sample_reference_coverage_probe_2.log",
        tmp_path / "sample_reference_coverage_probe_3.log",
    )


def test_clean_idxstats_selects_the_indexed_scan(tmp_path: Path) -> None:
    """Complete evidence with no placed-unmapped reads permits the indexed scan."""
    clean = "chr1\t20000\t600\t0\n*\t0\t0\t80\n"
    config = {"tools": {"samtools": "samtools"}, "cram": {"unmapped_scan": "auto"}}

    with patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, clean)):
        scan = choose_unmapped_scan("/run/view.cram", config, 4, str(tmp_path), "sample")

    assert scan == "indexed"


def test_failed_idxstats_selects_the_lossless_stream_scan(tmp_path: Path) -> None:
    """A command failure fails closed to the lossless scan."""
    with patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "bad index")):
        scan = choose_unmapped_scan("/run/view.cram", {}, 2, str(tmp_path), "sample")

    assert scan == "stream"


def test_malformed_idxstats_logs_the_offending_line(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Fail-closed selection still tells the operator which evidence line was invalid."""
    malformed = "chr1\t20000\tsix\t0\n*\t0\t0\t80\n"

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, malformed)),
        caplog.at_level("INFO"),
    ):
        scan = choose_unmapped_scan("/run/view.cram", {}, 2, str(tmp_path), "sample")

    assert scan == "stream"
    assert "idxstats output is malformed at line 1: 'chr1\\t20000\\tsix\\t0'" in caplog.text


def test_oversized_malformed_idxstats_log_is_bounded(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Untrusted idxstats output cannot expand a single operator diagnostic without bound."""
    malformed = f"chr1\t{'9' * 5000}\t600\t0\n*\t0\t0\t80\n"

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, malformed)),
        caplog.at_level("INFO"),
    ):
        scan = choose_unmapped_scan("/run/view.cram", {}, 2, str(tmp_path), "sample")

    messages = [
        record.getMessage() for record in caplog.records if "idxstats output is malformed" in record.getMessage()
    ]
    assert scan == "stream"
    assert len(messages) == 1
    assert "...<truncated>" in messages[0]
    assert len(messages[0]) <= 300


def test_bam_idxstats_with_placed_unmapped_reads_selects_the_complete_stream(tmp_path: Path) -> None:
    """The BAI tail shortcut cannot recover placed unmapped records."""
    placed = "chr1\t20000\t600\t329\n*\t0\t0\t4478\n"
    config = {"tools": {"samtools": "samtools"}, "bam": {"unmapped_scan": "auto"}}

    with patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, placed)):
        scan = choose_unmapped_scan("/run/view.bam", config, 4, str(tmp_path), "sample", file_format="bam")

    assert scan == "stream"


def test_an_invalid_configured_scan_is_rejected(tmp_path: Path) -> None:
    """Only the shipped automatic and explicit safe scan modes are accepted."""
    config = {"cram": {"unmapped_scan": "quick"}}

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "")),
        pytest.raises(ValueError, match="invalid unmapped scan mode"),
    ):
        choose_unmapped_scan("/run/view.cram", config, 2, str(tmp_path), "sample")


def test_shipped_config_exposes_the_actual_safe_preflight_policy() -> None:
    """Replacement configs can copy every new key from the shipped full config."""
    config_path = Path(__file__).parents[2] / "vntyper" / "config.json"
    with config_path.open() as config_handle:
        config = json.load(config_handle)

    assert config["cram"] == {
        "allow_ambient_reference_resolution": False,
        "local_ref_path": "%2s/%2s/%s",
        "reference_probe_timeout_seconds": 120,
        "unmapped_scan": "auto",
        "reference_candidate_order": [
            "cli",
            "config_cram_reference",
            "config_bwa_reference",
            "htslib_resolved",
        ],
    }
    assert config["bam"] == {"unmapped_scan": "auto"}
    assert config["reference_data"]["cram_reference_hg19"] is None
    assert config["reference_data"]["cram_reference_hg38"] is None
    assert config["assembly_detection"]["naming_convention_threshold"] == 0.5
    assert config["assembly_detection"]["primary_contig_patterns"] == [
        "^chr[0-9XYM]+$",
        r"^NC_\d{6}\.\d+$",
        "^([0-9]+|X|Y|MT?)$",
    ]
    assert "read_layout" not in config


def test_a_cram_plan_uses_configured_candidate_order_and_idxstats_evidence(tmp_path: Path) -> None:
    """The returned immutable plan records the first proven configured reference."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.cram"
    alignment.write_bytes(b"CRAM")
    (input_dir / "sample.cram.crai").write_bytes(b"CRAI")
    cli_reference = tmp_path / "cli.fa"
    bwa_reference = tmp_path / "bwa.fa"
    cram_reference = tmp_path / "cram.fa"
    for reference in (cli_reference, bwa_reference, cram_reference):
        reference.write_text(">chr1\nACGT\n")
        Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n")
    config = {
        "tools": {"samtools": "samtools"},
        "cram": {
            "unmapped_scan": "auto",
            "reference_candidate_order": [
                "config_bwa_reference",
                "cli",
                "config_cram_reference",
                "htslib_resolved",
            ],
        },
        "reference_data": {
            "cram_reference_hg19": str(cram_reference),
            "bwa_reference_hg19": str(bwa_reference),
        },
    }

    def successful_commands(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        if arguments[1] == "index":
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"CRAI")
            return True, ""
        if " idxstats " in f" {command} ":
            return True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"
        return True, "decoded"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=successful_commands):
        plan = run_preflight(
            str(alignment),
            str(output),
            "sample",
            "cram",
            config,
            4,
            region="chr1:1-2",
            reference_assembly="hg19",
            reference_fasta=str(cli_reference),
            header_contigs=("chr1",),
            m5="abc",
        )

    assert Path(plan.view_path).samefile(alignment)
    assert plan.reference_path == str(bwa_reference)
    assert plan.reference_source == "config_bwa_reference"
    assert plan.uncovered_contigs == ()
    assert plan.unmapped_scan == "indexed"


def test_default_policy_probes_the_cli_reference_before_configured_candidates(tmp_path: Path) -> None:
    """The real orchestrator must issue the CLI ``-T`` probe first by default."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.cram"
    alignment.write_bytes(b"CRAM")
    (input_dir / "sample.cram.crai").write_bytes(b"CRAI")
    cli_reference = tmp_path / "cli reference.fa"
    cram_reference = tmp_path / "configured-cram.fa"
    bwa_reference = tmp_path / "configured-bwa.fa"
    for reference in (cli_reference, cram_reference, bwa_reference):
        reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    config = {
        "reference_data": {
            "cram_reference_hg19": str(cram_reference),
            "bwa_reference_hg19": str(bwa_reference),
        }
    }
    reference_probes: list[str] = []

    def commands(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        if arguments[1] == "index":
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"CRAI")
            return True, ""
        if " idxstats " in f" {command} ":
            return True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"
        reference_probes.append(command)
        return str(cram_reference) in shlex.split(command), "decode result"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=commands):
        plan = run_preflight(
            str(alignment),
            str(output),
            "sample",
            "cram",
            config,
            1,
            region="chr1:1-2",
            reference_fasta=str(cli_reference),
            header_contigs=("chr1",),
        )

    assert [shlex.split(command)[shlex.split(command).index("-T") + 1] for command in reference_probes] == [
        str(cli_reference),
        str(cram_reference),
    ]
    assert plan.reference_source == "config_cram_reference"


def test_a_bam_plan_builds_an_index_selects_a_scan_and_probes_the_target(tmp_path: Path) -> None:
    """BAM preflight proves its own index, unmapped scan, and target retrieval."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"BAM")
    (input_dir / "sample.bam.bai").write_bytes(b"BAI")
    commands: list[str] = []

    def successful_probe(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        commands.append(command)
        arguments = shlex.split(command)
        if arguments[1] == "index":
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"BAI")
        if arguments[1] == "idxstats":
            return True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=successful_probe):
        plan = run_preflight(
            str(alignment),
            str(output),
            "sample",
            "bam",
            {},
            2,
            region="chr1:1-2",
            fast_mode=False,
        )

    assert len(commands) == 3
    assert any(" index " in f" {command} " for command in commands)
    assert any(" idxstats " in f" {command} " for command in commands)
    assert " -T " not in f" {commands[-1]} "
    assert plan.reference_path is None
    assert plan.reference_source == "not-required"
    assert plan.unmapped_scan == "indexed"
    assert plan.index_path.endswith(".bai")


def test_a_bam_plan_records_stream_when_idxstats_reports_placed_unmapped(tmp_path: Path) -> None:
    """BAM preflight must carry the lossless scan decision into conversion."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"BAM")

    def commands(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        if arguments[1] == "index":
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"BAI")
            return True, ""
        if arguments[1] == "idxstats":
            return True, "chr1\t100\t2\t3\n*\t0\t0\t4\n"
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=commands):
        plan = run_preflight(
            str(alignment),
            str(output),
            "sample",
            "bam",
            {"bam": {"unmapped_scan": "auto"}},
            2,
            region="chr1:1-2",
        )

    assert plan.unmapped_scan == "stream"


@pytest.mark.parametrize("file_format", ["bam", "cram"])
def test_fast_mode_ignores_an_unused_invalid_unmapped_scan(file_format: str, tmp_path: Path) -> None:
    """Historical target-only mode cannot be rejected by a scan it never consumes."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / f"sample.{file_format}"
    alignment.write_bytes(file_format.upper().encode())
    reference = tmp_path / "reference.fa" if file_format == "cram" else None
    if reference is not None:
        reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    commands: list[str] = []

    def successful_target(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        commands.append(command)
        arguments = shlex.split(command)
        if arguments[1] == "index":
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"INDEX")
        return True, "decoded"

    config = {file_format: {"unmapped_scan": "not-a-real-scan"}}
    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=successful_target):
        plan = run_preflight(
            str(alignment),
            str(output),
            "sample",
            file_format,
            config,
            2,
            region="chr1:1-2",
            reference_fasta=str(reference) if reference is not None else None,
            fast_mode=True,
        )

    assert plan.unmapped_scan == "not-required"
    assert not any(" idxstats " in f" {command} " for command in commands)


def test_fast_cram_mode_does_not_whole_stream_probe_an_unused_stream_scan(tmp_path: Path) -> None:
    """Fast CRAM proves its target but never decodes the whole file for skipped recovery."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.cram"
    alignment.write_bytes(b"CRAM")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    commands: list[str] = []

    def successful_target(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        commands.append(command)
        arguments = shlex.split(command)
        if arguments[1] == "index":
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"INDEX")
        return True, "decoded"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=successful_target):
        plan = run_preflight(
            str(alignment),
            str(output),
            "sample",
            "cram",
            {"cram": {"unmapped_scan": "stream"}},
            2,
            region="chr1:1-2",
            reference_fasta=str(reference),
            fast_mode=True,
        )

    assert plan.unmapped_scan == "not-required"
    reference_commands = [command for command in commands if " view " in f" {command} "]
    assert len(reference_commands) == 1
    assert " -P " in f" {reference_commands[0]} "
    assert " -h " not in f" {reference_commands[0]} "


def test_a_bam_whose_index_cannot_retrieve_the_target_fails_preflight(tmp_path: Path) -> None:
    """A stale or incompatible BAM index cannot reach genotyping stages."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"BAM")
    (input_dir / "sample.bam.bai").write_bytes(b"BAI")

    def fail_probe_only(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        if arguments[1] == "index":
            Path(arguments[arguments.index("-o") + 1]).write_bytes(b"BAI")
            return True, ""
        if arguments[1] == "idxstats":
            return True, "chr1\t4\t1\t0\n*\t0\t0\t2\n"
        return False, "stale index"

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=fail_probe_only),
        pytest.raises(RuntimeError, match="stale index"),
    ):
        run_preflight(str(alignment), str(output), "sample", "bam", {}, 2, region="chr1:1-2")
