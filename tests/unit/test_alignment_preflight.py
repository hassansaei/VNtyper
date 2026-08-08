"""Unit tests for the alignment preflight I/O seam."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_preflight import (
    HTSLIB_REFERENCE_SOURCE,
    build_alignment_view,
    capture_command,
    choose_unmapped_scan,
    resolve_reference,
    run_preflight,
)

pytestmark = pytest.mark.unit


def _config() -> dict:
    return {"tools": {"samtools": "samtools"}}


class TestAlignmentView:
    """The run-local view keeps every write inside the output tree."""

    def test_an_existing_index_is_symlinked_not_rebuilt(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.cram"
        index = input_dir / "sample.cram.crai"
        alignment.touch()
        index.touch()

        with patch("vntyper.scripts.alignment_preflight.run_command") as run_command:
            view_path, index_path = build_alignment_view(
                str(alignment), str(output_dir), "sample", "cram", _config(), threads=4
            )

        assert Path(view_path).samefile(alignment)
        assert Path(index_path).is_symlink()
        assert Path(index_path).samefile(index)
        run_command.assert_not_called()

    def test_a_missing_index_is_built_beside_the_view_never_beside_the_input(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.bam"
        alignment.touch()
        expected_index = output_dir / "sample.bam.bai"

        def build_index(command: str, log_file: str, critical: bool = False) -> bool:
            assert str(output_dir / "sample.bam") in command
            assert str(alignment) not in command
            assert Path(log_file).parent == output_dir
            assert critical is False
            expected_index.touch()
            return True

        with patch("vntyper.scripts.alignment_preflight.run_command", side_effect=build_index):
            view_path, index_path = build_alignment_view(
                str(alignment), str(output_dir), "sample", "bam", _config(), threads=3
            )

        assert Path(view_path).samefile(alignment)
        assert Path(index_path) == expected_index
        assert not (input_dir / "sample.bam.bai").exists()
        assert not (input_dir / "sample.bai").exists()

    def test_the_input_directory_is_untouched_even_when_it_is_read_only(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.cram"
        alignment.write_bytes(b"CRAM")
        before = {path.name: path.read_bytes() for path in input_dir.iterdir()}
        expected_index = output_dir / "readonly.cram.crai"
        input_dir.chmod(0o555)

        def build_index(command: str, log_file: str, critical: bool = False) -> bool:
            del command, log_file, critical
            expected_index.touch()
            return True

        try:
            with patch("vntyper.scripts.alignment_preflight.run_command", side_effect=build_index):
                build_alignment_view(str(alignment), str(output_dir), "readonly", "cram", _config(), threads=2)
        finally:
            input_dir.chmod(0o755)

        after = {path.name: path.read_bytes() for path in input_dir.iterdir()}
        assert after == before

    def test_a_stale_view_pointing_at_a_different_alignment_is_replaced(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        current = input_dir / "current.cram"
        stale = input_dir / "stale.cram"
        current.touch()
        stale.touch()
        (input_dir / "current.cram.crai").touch()
        view = output_dir / "sample.cram"
        view.symlink_to(stale)

        with patch("vntyper.scripts.alignment_preflight.run_command"):
            view_path, _ = build_alignment_view(str(current), str(output_dir), "sample", "cram", _config(), threads=2)

        assert Path(view_path).samefile(current)
        assert os.readlink(view_path) == str(current.resolve())

    def test_a_view_already_pointing_at_this_same_input_is_reused(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.cram"
        alignment.touch()
        (input_dir / "sample.cram.crai").touch()
        view = output_dir / "sample.cram"
        view.symlink_to(os.path.relpath(alignment, output_dir))

        with (
            patch("vntyper.scripts.alignment_preflight.os.unlink", side_effect=AssertionError("view was replaced")),
            patch("vntyper.scripts.alignment_preflight.run_command"),
        ):
            view_path, _ = build_alignment_view(str(alignment), str(output_dir), "sample", "cram", _config(), threads=2)

        assert Path(view_path).samefile(alignment)

    def test_a_regular_file_at_the_view_path_is_replaced_without_touching_the_input(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.cram"
        alignment.write_bytes(b"input")
        (input_dir / "sample.cram.crai").touch()
        view = output_dir / "sample.cram"
        view.write_bytes(b"stale output")

        view_path, _ = build_alignment_view(str(alignment), str(output_dir), "sample", "cram", _config(), threads=2)

        assert Path(view_path).is_symlink()
        assert Path(view_path).samefile(alignment)
        assert alignment.read_bytes() == b"input"

    def test_a_non_fast_bam_builds_the_required_bai_when_only_a_csi_exists(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.bam"
        alignment.touch()
        (input_dir / "sample.bam.csi").touch()
        expected_index = output_dir / "sample.bam.bai"

        def build_index(command: str, log_file: str, critical: bool = False) -> bool:
            del command, log_file, critical
            expected_index.touch()
            return True

        with patch("vntyper.scripts.alignment_preflight.run_command", side_effect=build_index):
            _, index_path = build_alignment_view(
                str(alignment),
                str(output_dir),
                "sample",
                "bam",
                _config(),
                threads=2,
                bai_only=True,
            )

        assert Path(index_path) == expected_index
        assert not Path(index_path).is_symlink()

    def test_a_stale_index_symlink_is_removed_before_building_for_a_new_input(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "current.cram"
        alignment.touch()
        old_index = input_dir / "old.cram.crai"
        old_index.write_bytes(b"old index")
        stale_view_index = output_dir / "sample.cram.crai"
        stale_view_index.symlink_to(old_index)

        def build_index(command: str, log_file: str, critical: bool = False) -> bool:
            del command, log_file, critical
            assert not stale_view_index.exists()
            assert not stale_view_index.is_symlink()
            stale_view_index.write_bytes(b"current index")
            return True

        with patch("vntyper.scripts.alignment_preflight.run_command", side_effect=build_index):
            _, index_path = build_alignment_view(
                str(alignment), str(output_dir), "sample", "cram", _config(), threads=2
            )

        assert Path(index_path).read_bytes() == b"current index"
        assert old_index.read_bytes() == b"old index"

    def test_a_stale_higher_priority_csi_cannot_override_the_current_bai(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "current.bam"
        alignment.touch()
        current_bai = input_dir / "current.bam.bai"
        current_bai.touch()
        old_csi = input_dir / "old.bam.csi"
        old_csi.touch()
        stale_view_csi = output_dir / "sample.bam.csi"
        stale_view_csi.symlink_to(old_csi)

        with patch("vntyper.scripts.alignment_preflight.run_command"):
            _, index_path = build_alignment_view(str(alignment), str(output_dir), "sample", "bam", _config(), threads=2)

        assert Path(index_path).samefile(current_bai)
        assert not stale_view_csi.exists()
        assert not stale_view_csi.is_symlink()

    @pytest.mark.parametrize("command_succeeds", [False, True])
    def test_a_failed_or_missing_index_build_raises_an_actionable_error(
        self, tmp_path: Path, command_succeeds: bool
    ) -> None:
        alignment = tmp_path / "sample.cram"
        alignment.touch()
        output_dir = tmp_path / "output"

        with (
            patch("vntyper.scripts.alignment_preflight.run_command", return_value=command_succeeds),
            pytest.raises(RuntimeError, match="samtools index"),
        ):
            build_alignment_view(str(alignment), str(output_dir), "sample", "cram", _config(), threads=2)

    def test_an_unknown_alignment_format_is_rejected_before_creating_a_view(self, tmp_path: Path) -> None:
        alignment = tmp_path / "sample.sam"
        alignment.touch()
        output_dir = tmp_path / "output"

        with pytest.raises(ValueError, match="unknown alignment format"):
            build_alignment_view(str(alignment), str(output_dir), "sample", "sam", _config(), threads=2)

        assert not output_dir.exists()

    def test_an_input_at_the_exact_view_path_is_rejected_without_deleting_it(self, tmp_path: Path) -> None:
        alignment = tmp_path / "sample.cram"
        alignment.write_bytes(b"patient alignment")
        index = tmp_path / "sample.cram.crai"
        index.write_bytes(b"patient index")

        with pytest.raises(ValueError, match="same path"):
            build_alignment_view(str(alignment), str(tmp_path), "sample", "cram", _config(), threads=2)

        assert alignment.is_file()
        assert not alignment.is_symlink()
        assert alignment.read_bytes() == b"patient alignment"
        assert index.is_file()
        assert not index.is_symlink()
        assert index.read_bytes() == b"patient index"

    def test_the_alignment_view_cannot_escape_the_output_directory(self, tmp_path: Path) -> None:
        alignment = tmp_path / "input.cram"
        alignment.touch()
        (tmp_path / "input.cram.crai").touch()
        output_dir = tmp_path / "output"
        escaped_view = tmp_path / "escaped.cram"

        with pytest.raises(ValueError, match="inside output directory"):
            build_alignment_view(str(alignment), str(output_dir), "../escaped", "cram", _config(), threads=2)

        assert not output_dir.exists()
        assert not escaped_view.exists()


class TestCaptureCommand:
    """Captured stdout is both returned to parsers and retained in a run log."""

    def test_a_success_returns_stdout_and_writes_the_log(self, tmp_path: Path) -> None:
        completed = subprocess.CompletedProcess(args="command", returncode=0, stdout="chr1\t1\t2\t0\n")
        log_file = tmp_path / "capture.log"

        with patch("vntyper.scripts.alignment_preflight.subprocess.run", return_value=completed):
            result = capture_command("command", str(log_file), cwd=str(tmp_path))

        assert result == (True, "chr1\t1\t2\t0\n")
        assert log_file.read_text() == "chr1\t1\t2\t0\n"

    def test_a_failure_returns_false_with_the_diagnostic_instead_of_raising(self, tmp_path: Path) -> None:
        completed = subprocess.CompletedProcess(args="command", returncode=1, stdout="bad index\n")
        log_file = tmp_path / "capture.log"

        with patch("vntyper.scripts.alignment_preflight.subprocess.run", return_value=completed):
            result = capture_command("command", str(log_file))

        assert result == (False, "bad index\n")
        assert log_file.read_text() == "bad index\n"


class TestChooseUnmappedScan:
    """The I/O seam feeds complete idxstats evidence to the fail-closed decision."""

    def test_clean_idxstats_selects_the_indexed_scan(self, tmp_path: Path) -> None:
        clean = "chr1\t20000\t600\t0\n*\t0\t0\t80\n"
        config = {"tools": {"samtools": "samtools"}, "cram": {"unmapped_scan": "auto"}}

        with patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, clean)):
            scan = choose_unmapped_scan("/run/view.cram", config, 4, str(tmp_path), "sample")

        assert scan == "indexed"

    def test_failed_idxstats_selects_the_lossless_stream_scan(self, tmp_path: Path) -> None:
        with patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "bad index")):
            scan = choose_unmapped_scan("/run/view.cram", {}, 2, str(tmp_path), "sample")

        assert scan == "stream"

    def test_an_invalid_configured_scan_is_rejected(self, tmp_path: Path) -> None:
        config = {"cram": {"unmapped_scan": "quick"}}

        with (
            patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "")),
            pytest.raises(ValueError, match="invalid unmapped scan mode"),
        ):
            choose_unmapped_scan("/run/view.cram", config, 2, str(tmp_path), "sample")


class TestResolveReference:
    """CRAM references are authorised only by the real-shape decode probe."""

    def test_explicit_candidates_are_probed_before_the_no_reference_one(self, tmp_path: Path) -> None:
        commands: list[str] = []

        def fail_until_htslib(command: str, log_file: str, cwd: str | None = None) -> tuple[bool, str]:
            del log_file, cwd
            commands.append(command)
            return (" -T " not in f" {command} ", "did not decode")

        candidates = (
            ("cli", "/refs/cli.fa"),
            ("config_cram_reference", "/refs/whole.fa"),
            ("config_bwa_reference", "/refs/chr1.fa"),
        )
        with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=fail_until_htslib):
            resolve_reference(
                "/run/view.cram",
                candidates,
                "chr1:1-2",
                None,
                {},
                4,
                str(tmp_path),
                "sample",
                ("chr1",),
                "abc",
            )

        assert len(commands) == 4
        assert all(" -T " in f" {command} " for command in commands[:3])
        assert " -T " not in f" {commands[-1]} "

    def test_the_no_reference_candidate_is_last_and_is_recorded_as_htslib_resolved(self, tmp_path: Path) -> None:
        outcomes = [(False, "explicit failed"), (True, "")]

        with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=outcomes):
            reference, source, uncovered = resolve_reference(
                "/run/view.cram",
                (("cli", "/refs/cli.fa"),),
                "chr1:1-2",
                None,
                {},
                1,
                str(tmp_path),
                "sample",
                ("chr1",),
                "abc",
            )

        assert reference is None
        assert source == HTSLIB_REFERENCE_SOURCE == "htslib-resolved (header UR: or REF_PATH)"
        assert uncovered == ()

    def test_the_first_candidate_that_decodes_wins(self, tmp_path: Path) -> None:
        outcomes = [(False, "wrong digest"), (True, "decoded")]
        candidates = (("cli", "/refs/first.fa"), ("config_cram_reference", "/refs/second.fa"))

        with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=outcomes) as capture:
            reference, source, _ = resolve_reference(
                "/run/view.cram",
                candidates,
                "chr1:1-2",
                None,
                {},
                1,
                str(tmp_path),
                "sample",
                (),
                None,
            )

        assert reference == "/refs/second.fa"
        assert source == "config_cram_reference"
        assert capture.call_count == 2

    def test_no_candidate_decoding_raises_naming_every_one_with_its_reason(self, tmp_path: Path) -> None:
        candidates = (
            ("cli", "/refs/cli.fa"),
            ("config_cram_reference", None),
            ("config_bwa_reference", "/refs/chr1.fa"),
        )
        outcomes = [(False, "checksum mismatch"), (False, "missing chr2"), (False, "UR unavailable")]

        with (
            patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=outcomes),
            pytest.raises(ValueError) as error,
        ):
            resolve_reference(
                "/run/view.cram",
                candidates,
                "chr1:1-2",
                None,
                {},
                1,
                str(tmp_path),
                "sample",
                ("chr1", "chr2"),
                "012345",
            )

        message = str(error.value)
        for expected in (
            "contig=chr1",
            "M5=012345",
            "source=cli",
            "path=/refs/cli.fa",
            "checksum mismatch",
            "source=config_cram_reference",
            "not supplied",
            "source=config_bwa_reference",
            "missing chr2",
            HTSLIB_REFERENCE_SOURCE,
            "UR unavailable",
        ):
            assert expected in message

    def test_a_reference_error_names_the_region_contig_not_the_first_header_contig(self, tmp_path: Path) -> None:
        with (
            patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "not found")),
            pytest.raises(ValueError) as error,
        ):
            resolve_reference(
                "/run/view.cram",
                (),
                "chr1:100-200",
                None,
                {},
                1,
                str(tmp_path),
                "sample",
                ("chrM", "chr1"),
                "chr1-digest",
            )

        assert "contig=chr1" in str(error.value)
        assert "contig=chrM" not in str(error.value)

    def test_a_reference_error_names_the_bed_contig_not_the_first_header_contig(self, tmp_path: Path) -> None:
        bed_file = tmp_path / "target.bed"
        bed_file.write_text("# target\n\nchr2\t100\t200\n")

        with (
            patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "not found")),
            pytest.raises(ValueError) as error,
        ):
            resolve_reference(
                "/run/view.cram",
                (),
                "chr1:100-200",
                bed_file,
                {},
                1,
                str(tmp_path),
                "sample",
                ("chrM", "chr2"),
                "chr2-digest",
            )

        assert "contig=chr2" in str(error.value)
        assert "contig=chrM" not in str(error.value)

    def test_the_probe_uses_the_runs_own_region_and_bed_not_a_hardcoded_one(self, tmp_path: Path) -> None:
        commands: list[str] = []

        def succeed(command: str, log_file: str, cwd: str | None = None) -> tuple[bool, str]:
            del log_file, cwd
            commands.append(command)
            return True, ""

        with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=succeed):
            resolve_reference(
                "/run/view.cram",
                (("cli", "/refs/full.fa"),),
                "chr9:90-99",
                None,
                {},
                1,
                str(tmp_path),
                "region",
                (),
                None,
            )
            resolve_reference(
                "/run/view.cram",
                (("cli", "/refs/full.fa"),),
                "chr9:90-99",
                "/run/actual regions.bed",
                {},
                1,
                str(tmp_path),
                "bed",
                (),
                None,
            )

        assert "chr9:90-99" in commands[0]
        assert "-L '/run/actual regions.bed'" in commands[1]
        assert "chr9:90-99" not in commands[1]

    def test_a_reference_not_covering_every_header_contig_logs_a_warning(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        reference = tmp_path / "chr1.fa"
        reference.write_text(">chr1\nACGT\n")
        Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n")

        with (
            patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "")),
            caplog.at_level("WARNING"),
        ):
            _, _, uncovered = resolve_reference(
                "/run/view.cram",
                (("config_bwa_reference", str(reference)),),
                "chr1:1-2",
                None,
                {},
                1,
                str(tmp_path),
                "sample",
                ("chr1", "chr2"),
                "abc",
            )

        assert uncovered == ("chr2",)
        assert "chr2" in caplog.text


def test_shipped_config_exposes_the_actual_safe_preflight_policy() -> None:
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


class TestRunPreflight:
    """The orchestration returns one frozen plan backed by proven prerequisites."""

    def test_a_cram_plan_uses_configured_candidate_order_and_idxstats_evidence(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.cram"
        alignment.touch()
        (input_dir / "sample.cram.crai").touch()
        cli_reference = tmp_path / "cli.fa"
        bwa_reference = tmp_path / "bwa.fa"
        cram_reference = tmp_path / "cram.fa"
        for reference in (cli_reference, bwa_reference, cram_reference):
            reference.write_text(">chr1\nACGT\n")
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
                str(output_dir),
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

        assert plan.input_path == str(alignment)
        assert Path(plan.view_path).samefile(alignment)
        assert plan.file_format == "cram"
        assert plan.reference_path == str(bwa_reference)
        assert plan.reference_source == "config_bwa_reference"
        assert plan.uncovered_contigs == ()
        assert plan.unmapped_scan == "indexed"

    def test_a_bam_plan_probes_the_index_but_skips_reference_and_idxstats_resolution(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.bam"
        alignment.touch()
        (input_dir / "sample.bam.bai").touch()
        commands: list[str] = []

        def successful_probe(command: str, log_file: str, cwd: str | None = None) -> tuple[bool, str]:
            del log_file, cwd
            commands.append(command)
            return True, ""

        with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=successful_probe):
            plan = run_preflight(
                str(alignment),
                str(output_dir),
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

    def test_a_bam_whose_index_cannot_retrieve_the_region_fails_preflight(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        alignment = input_dir / "sample.bam"
        alignment.touch()
        (input_dir / "sample.bam.bai").touch()

        with (
            patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "stale index")),
            pytest.raises(RuntimeError, match="stale index"),
        ):
            run_preflight(
                str(alignment),
                str(output_dir),
                "sample",
                "bam",
                {},
                2,
                region="chr1:1-2",
            )

    def test_an_unknown_reference_candidate_source_is_rejected(self, tmp_path: Path) -> None:
        alignment = tmp_path / "sample.cram"
        alignment.touch()
        (tmp_path / "sample.cram.crai").touch()
        output_dir = tmp_path / "output"
        config = {"cram": {"reference_candidate_order": ["mystery", "htslib_resolved"]}}

        clean = "chr1\t4\t1\t0\n*\t0\t0\t2\n"
        with (
            patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, clean)),
            pytest.raises(ValueError, match="unknown reference candidate source: mystery"),
        ):
            run_preflight(
                str(alignment),
                str(output_dir),
                "sample",
                "cram",
                config,
                1,
                region="chr1:1-2",
            )
