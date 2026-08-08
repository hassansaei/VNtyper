"""Orchestration and shipped-policy tests for alignment preflight."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_preflight import choose_unmapped_scan, run_preflight

pytestmark = pytest.mark.unit


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
        "unmapped_scan": "auto",
        "reference_candidate_order": [
            "cli",
            "config_cram_reference",
            "config_bwa_reference",
            "htslib_resolved",
        ],
    }
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

    def successful_commands(command: str, log_file: str, cwd: str | None = None) -> tuple[bool, str]:
        del log_file, cwd
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


def test_a_bam_plan_probes_the_index_without_reference_or_idxstats_resolution(tmp_path: Path) -> None:
    """BAM needs one retrieval proof but no CRAM reference decision."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"BAM")
    (input_dir / "sample.bam.bai").write_bytes(b"BAI")
    commands: list[str] = []

    def successful_probe(command: str, log_file: str, cwd: str | None = None) -> tuple[bool, str]:
        del log_file, cwd
        commands.append(command)
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

    assert len(commands) == 1
    assert " -T " not in f" {commands[0]} "
    assert plan.reference_path is None
    assert plan.reference_source == "not-required"
    assert plan.unmapped_scan == "indexed"
    assert plan.index_path.endswith(".bai")


def test_a_bam_whose_index_cannot_retrieve_the_target_fails_preflight(tmp_path: Path) -> None:
    """A stale or incompatible BAM index cannot reach genotyping stages."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"BAM")
    (input_dir / "sample.bam.bai").write_bytes(b"BAI")

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "stale index")),
        pytest.raises(RuntimeError, match="stale index"),
    ):
        run_preflight(str(alignment), str(output), "sample", "bam", {}, 2, region="chr1:1-2")
