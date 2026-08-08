"""Filesystem-safety regressions for the run-local alignment preflight view."""

from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_preflight import build_alignment_view, capture_command, run_preflight

pytestmark = pytest.mark.unit


def _config() -> dict:
    """Return the minimal shipped-default-compatible tool configuration."""
    return {"tools": {"samtools": "samtools"}}


def _indexed_cram(directory: Path, name: str = "input") -> tuple[Path, Path]:
    """Create a literal CRAM and CRAI fixture in ``directory``."""
    directory.mkdir(parents=True, exist_ok=True)
    alignment = directory / f"{name}.cram"
    index = directory / f"{name}.cram.crai"
    alignment.write_bytes(b"patient alignment")
    index.write_bytes(b"patient index")
    return alignment, index


def _write_built_index(command: str) -> None:
    """Materialize the explicit temporary samtools-index destination."""
    arguments = shlex.split(command)
    Path(arguments[arguments.index("-o") + 1]).write_bytes(b"current index")


def test_an_output_name_cannot_traverse_an_in_output_symlink_component(tmp_path: Path) -> None:
    """A nested output name must not replace a victim outside the run directory."""
    alignment, _ = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    outside = tmp_path / "outside"
    output.mkdir()
    outside.mkdir()
    (output / "escape").symlink_to(outside, target_is_directory=True)
    victim = outside / "sample.cram"
    victim.write_bytes(b"outside victim")

    with pytest.raises(ValueError, match="single basename"):
        build_alignment_view(str(alignment), str(output), "escape/sample", "cram", _config(), threads=1)

    assert victim.is_file()
    assert not victim.is_symlink()
    assert victim.read_bytes() == b"outside victim"


def test_an_input_symlink_resolving_to_the_view_cannot_delete_that_target(tmp_path: Path) -> None:
    """Replacing the view must not unlink the file that an input symlink resolves to."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    output.mkdir()
    view = output / "sample.cram"
    view.write_bytes(b"patient target")
    input_link = input_dir / "input.cram"
    input_link.symlink_to(view)
    (input_dir / "input.cram.crai").write_bytes(b"patient index")

    with pytest.raises(ValueError, match="same file"):
        build_alignment_view(str(input_link), str(output), "sample", "cram", _config(), threads=1)

    assert view.is_file()
    assert not view.is_symlink()
    assert view.read_bytes() == b"patient target"
    assert input_link.samefile(view)


def test_a_hardlinked_view_collision_is_rejected_without_mutation(tmp_path: Path) -> None:
    """A regular view entry sharing the input inode is patient data, not stale output."""
    alignment, _ = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    view = output / "sample.cram"
    os.link(alignment, view)

    with pytest.raises(ValueError, match="same file"):
        build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    assert alignment.read_bytes() == b"patient alignment"
    assert view.read_bytes() == b"patient alignment"
    assert not view.is_symlink()


@pytest.mark.parametrize("protected_source", ["alignment", "index"])
def test_a_derived_index_collision_with_a_source_inode_is_rejected_before_creating_the_view(
    tmp_path: Path, protected_source: str
) -> None:
    """A derived index name must not replace the alignment or its selected source index."""
    alignment, source_index = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    view_index = output / "sample.cram.crai"
    source = alignment if protected_source == "alignment" else source_index
    os.link(source, view_index)

    with pytest.raises(ValueError, match="same file"):
        build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    assert not (output / "sample.cram").exists()
    assert source.read_bytes() == view_index.read_bytes()
    assert not view_index.is_symlink()


def test_a_dangling_view_link_is_atomically_replaced_without_unlinking_the_final_entry(tmp_path: Path) -> None:
    """A dangling stale symlink remains safely replaceable without an unlink window."""
    alignment, _ = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    view = output / "sample.cram"
    view.symlink_to(tmp_path / "missing.cram")

    with patch("vntyper.scripts.alignment_preflight.os.unlink", side_effect=AssertionError("blind unlink")):
        view_path, _ = build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    assert Path(view_path).is_symlink()
    assert Path(view_path).samefile(alignment)


@pytest.mark.parametrize("entry_kind", ["file", "directory"])
def test_a_wrong_type_view_collision_is_rejected_without_replacement(tmp_path: Path, entry_kind: str) -> None:
    """Only a symlink is an owned reusable view entry; files and directories are rejected."""
    alignment, _ = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    view = output / "sample.cram"
    if entry_kind == "file":
        view.write_bytes(b"unowned output")
    else:
        view.mkdir()

    with pytest.raises(ValueError, match="wrong type"):
        build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    if entry_kind == "file":
        assert view.read_bytes() == b"unowned output"
    else:
        assert view.is_dir()


@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_capture_command_rejects_a_log_alias_of_a_protected_input(tmp_path: Path, alias_kind: str) -> None:
    """A final log entry must never truncate a protected patient-data inode."""
    patient = tmp_path / "patient.cram"
    patient.write_bytes(b"patient bytes")
    log_file = patient if alias_kind == "direct" else tmp_path / "command.log"
    if alias_kind == "symlink":
        log_file.symlink_to(patient)
    elif alias_kind == "hardlink":
        os.link(patient, log_file)

    with patch("vntyper.scripts.alignment_preflight.subprocess.run") as run:
        success, diagnostic = capture_command("command", str(log_file), protected_paths=(patient,))

    assert success is False
    assert "protected" in diagnostic
    assert patient.read_bytes() == b"patient bytes"
    run.assert_not_called()


def test_capture_command_atomically_replaces_a_stale_log_symlink_without_following_it(tmp_path: Path) -> None:
    """Rerun logging replaces the directory entry and preserves its former target."""
    victim = tmp_path / "outside.log"
    victim.write_bytes(b"outside victim")
    log_file = tmp_path / "command.log"
    log_file.symlink_to(victim)
    completed = subprocess.CompletedProcess(args="command", returncode=0, stdout="new log\n")

    with patch("vntyper.scripts.alignment_preflight.subprocess.run", return_value=completed):
        success, output = capture_command("command", str(log_file))

    assert (success, output) == (True, "new log\n")
    assert victim.read_bytes() == b"outside victim"
    assert log_file.is_file()
    assert not log_file.is_symlink()
    assert log_file.read_text() == "new log\n"


def test_capture_command_turns_an_oserror_into_a_safe_logged_diagnostic(tmp_path: Path) -> None:
    """An OS-level spawn failure is a parseable failure, not an escaped exception."""
    log_file = tmp_path / "command.log"

    with patch("vntyper.scripts.alignment_preflight.subprocess.run", side_effect=OSError("cannot spawn")):
        success, diagnostic = capture_command("command", str(log_file))

    assert success is False
    assert "cannot spawn" in diagnostic
    assert "cannot spawn" in log_file.read_text()


@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_an_index_log_alias_of_the_input_is_rejected_before_the_view_is_created(
    tmp_path: Path, alias_kind: str
) -> None:
    """All derived index-log paths are checked before any run-local entry is mutated."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    output.mkdir()
    alignment = input_dir / "sample.cram"
    alignment.write_bytes(b"patient alignment")
    log_file = output / "sample_index.log"
    if alias_kind == "symlink":
        log_file.symlink_to(alignment)
    else:
        os.link(alignment, log_file)

    with (
        patch("vntyper.scripts.alignment_preflight.run_command", return_value=False, create=True),
        pytest.raises(ValueError, match="protected"),
    ):
        build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    assert alignment.read_bytes() == b"patient alignment"
    assert not (output / "sample.cram").exists()


@pytest.mark.parametrize(
    ("file_format", "log_name"),
    [
        ("cram", "sample_idxstats.log"),
        ("cram", "sample_reference_probe_4.log"),
        ("bam", "sample_alignment_probe.log"),
    ],
)
def test_every_later_preflight_log_is_validated_before_the_alignment_view_is_created(
    tmp_path: Path, file_format: str, log_name: str
) -> None:
    """A later-stage log collision must fail before build_alignment_view mutates output."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    output.mkdir()
    alignment = input_dir / f"sample.{file_format}"
    alignment.write_bytes(b"patient alignment")
    suffix = "crai" if file_format == "cram" else "bai"
    (input_dir / f"sample.{file_format}.{suffix}").write_bytes(b"patient index")
    (output / log_name).symlink_to(alignment)
    clean_idxstats = "chr1\t4\t1\t0\n*\t0\t0\t2\n"

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, clean_idxstats)),
        pytest.raises(ValueError, match="protected"),
    ):
        run_preflight(
            str(alignment),
            str(output),
            "sample",
            file_format,
            {},
            1,
            region="chr1:1-2",
        )

    assert alignment.read_bytes() == b"patient alignment"
    assert not (output / f"sample.{file_format}").exists()


def test_an_existing_cram_csi_is_reused_beside_the_view_without_reindexing(tmp_path: Path) -> None:
    """A CRAM CSI proven usable by htslib is a reusable general preflight index."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    output.mkdir()
    alignment = input_dir / "sample.cram"
    source_csi = input_dir / "sample.cram.csi"
    alignment.write_bytes(b"CRAM\x02")
    source_csi.write_bytes(b"CSI\x01")

    with patch("vntyper.scripts.alignment_preflight.capture_command") as capture:
        _, index_path = build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    assert Path(index_path).is_symlink()
    assert Path(index_path).samefile(source_csi)
    capture.assert_not_called()


def test_a_stale_cram_csi_cannot_override_the_selected_crai(tmp_path: Path) -> None:
    """Selecting CRAI removes a stale co-located CSI before downstream htslib calls."""
    alignment, source_crai = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    stale_source = tmp_path / "stale.cram.csi"
    stale_source.write_bytes(b"stale CSI")
    stale_view_csi = output / "sample.cram.csi"
    stale_view_csi.symlink_to(stale_source)

    _, index_path = build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    assert Path(index_path).samefile(source_crai)
    assert not os.path.lexists(stale_view_csi)


def test_an_existing_index_is_symlinked_without_invoking_the_index_builder(tmp_path: Path) -> None:
    """A selected source index is represented by an owned run-local symlink."""
    alignment, source_index = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"

    with patch("vntyper.scripts.alignment_preflight.capture_command") as capture:
        view_path, index_path = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=4
        )

    assert Path(view_path).samefile(alignment)
    assert Path(index_path).is_symlink()
    assert Path(index_path).samefile(source_index)
    capture.assert_not_called()


def test_a_missing_index_is_built_beside_the_view_and_never_beside_the_input(tmp_path: Path) -> None:
    """Index construction targets an atomic temporary path in the output directory."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"BAM")

    def build_index(
        command: str,
        log_file: str,
        cwd: str | None = None,
        *,
        protected_paths: tuple[str | Path, ...] = (),
    ) -> tuple[bool, str]:
        assert str(output / "sample.bam") in command
        assert str(alignment) not in command
        assert Path(log_file).parent == output
        assert cwd is None
        assert protected_paths
        _write_built_index(command)
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_index):
        _, index_path = build_alignment_view(str(alignment), str(output), "sample", "bam", _config(), threads=3)

    assert Path(index_path) == output / "sample.bam.bai"
    assert Path(index_path).read_bytes() == b"current index"
    assert not (input_dir / "sample.bam.bai").exists()
    assert not (input_dir / "sample.bai").exists()


@pytest.mark.parametrize("command_succeeds", [False, True])
def test_a_failed_or_empty_index_build_raises_an_actionable_error(tmp_path: Path, command_succeeds: bool) -> None:
    """Both command failure and a missing promised artifact stop preflight."""
    alignment = tmp_path / "sample.cram"
    alignment.write_bytes(b"CRAM")

    with (
        patch(
            "vntyper.scripts.alignment_preflight.capture_command",
            return_value=(command_succeeds, ""),
        ),
        pytest.raises(RuntimeError, match="samtools index"),
    ):
        build_alignment_view(str(alignment), str(tmp_path / "output"), "sample", "cram", _config(), threads=2)


def test_an_unknown_alignment_format_is_rejected_before_output_creation(tmp_path: Path) -> None:
    """Unsupported formats cannot create even a partial run-local view."""
    alignment = tmp_path / "sample.sam"
    alignment.touch()
    output = tmp_path / "output"

    with pytest.raises(ValueError, match="unknown alignment format"):
        build_alignment_view(str(alignment), str(output), "sample", "sam", _config(), threads=2)

    assert not output.exists()


def test_a_view_already_pointing_at_the_same_input_is_reused(tmp_path: Path) -> None:
    """A correct owned view survives a rerun without replacement."""
    alignment, source_index = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    view = output / "sample.cram"
    view.symlink_to(os.path.relpath(alignment, output))
    (output / "sample.cram.crai").symlink_to(os.path.relpath(source_index, output))

    with patch("vntyper.scripts.alignment_preflight.os.replace", side_effect=AssertionError("view replaced")):
        view_path, _ = build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=2)

    assert Path(view_path).samefile(alignment)


@pytest.mark.parametrize(
    ("return_code", "expected_success", "output"),
    [(0, True, "chr1\t1\t2\t0\n"), (1, False, "bad index\n")],
)
def test_capture_command_returns_and_persists_complete_output(
    tmp_path: Path, return_code: int, expected_success: bool, output: str
) -> None:
    """The parser result and durable diagnostic contain the same complete text."""
    completed = subprocess.CompletedProcess(args="command", returncode=return_code, stdout=output)
    log_file = tmp_path / "capture.log"

    with patch("vntyper.scripts.alignment_preflight.subprocess.run", return_value=completed):
        result = capture_command("command", str(log_file), cwd=str(tmp_path))

    assert result == (expected_success, output)
    assert log_file.read_text() == output
