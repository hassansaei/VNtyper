"""Filesystem-safety regressions for the run-local alignment preflight view."""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
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


def _write_successful_index(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
    """Materialize a valid index and return the capture-command success contract."""
    _write_built_index(command)
    return True, ""


def test_an_existing_source_index_is_rebuilt_from_the_alignment_before_use(tmp_path: Path) -> None:
    """A filename-compatible but wrong-sample index must never be trusted."""
    alignment, source_index = _indexed_cram(tmp_path / "input")
    source_index.write_bytes(b"valid index for a different sample")
    output = tmp_path / "output"

    def build_current_index(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        _write_built_index(command)
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_current_index) as capture:
        view_path, view_index, binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=3
        )

    assert Path(view_path).is_symlink()
    assert Path(view_path).samefile(alignment)
    assert Path(view_index).read_bytes() == b"current index"
    assert not Path(view_index).is_symlink()
    assert source_index.read_bytes() == b"valid index for a different sample"
    assert "-@ 3" in capture.call_args.args[0]
    binding.close()


@pytest.mark.parametrize("input_kind", ["regular", "symlink"])
def test_atomic_input_replacement_after_preflight_cannot_change_the_view_bytes(tmp_path: Path, input_kind: str) -> None:
    """A later consumer must see the bytes whose run-local index was built."""
    alignment, _ = _indexed_cram(tmp_path / "input")
    original_bytes = alignment.read_bytes()
    input_path = alignment
    if input_kind == "symlink":
        input_path = tmp_path / "input-link.cram"
        input_path.symlink_to(alignment)
    output = tmp_path / "output"

    with patch(
        "vntyper.scripts.alignment_preflight.capture_command",
        side_effect=_write_successful_index,
    ):
        view_path, index_path, binding = build_alignment_view(
            str(input_path), str(output), "sample", "cram", _config(), threads=1
        )

    replacement = tmp_path / "replacement.cram"
    replacement.write_bytes(b"different alignment at the same input pathname")
    replacement.replace(alignment)

    assert alignment.read_bytes() == b"different alignment at the same input pathname"
    consumer_bytes = subprocess.check_output(
        [
            sys.executable,
            "-c",
            "from pathlib import Path; import sys; sys.stdout.buffer.write(Path(sys.argv[1]).read_bytes())",
            view_path,
        ]
    )
    assert consumer_bytes == original_bytes
    assert Path(index_path).read_bytes() == b"current index"
    assert binding.is_open is True
    binding.close()
    assert binding.is_open is False
    with pytest.raises(FileNotFoundError):
        Path(view_path).read_bytes()


def test_an_input_that_cannot_be_descriptor_bound_refuses_before_index_work(tmp_path: Path) -> None:
    """An unavailable immutable pathname binding must fail before preflight commands."""
    alignment, _ = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    real_open = os.open

    def reject_alignment_open(
        path: str | os.PathLike[str],
        flags: int,
        mode: int = 0o777,
        *,
        dir_fd: int | None = None,
    ) -> int:
        if os.fspath(path) == str(alignment):
            raise OSError("descriptor binding unavailable")
        return real_open(path, flags, mode, dir_fd=dir_fd)

    with (
        patch("os.open", side_effect=reject_alignment_open),
        patch(
            "vntyper.scripts.alignment_preflight.capture_command",
            side_effect=_write_successful_index,
        ) as capture,
        pytest.raises(RuntimeError, match="stable run-local alignment binding"),
    ):
        build_alignment_view(str(alignment), str(output), "sample", "cram", _config(), threads=1)

    capture.assert_not_called()
    assert alignment.read_bytes() == b"patient alignment"
    assert not os.path.lexists(output / "sample.cram")


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
    real_unlink = os.unlink

    def reject_blind_view_unlink(path: str | os.PathLike[str]) -> None:
        if Path(path) == view:
            raise AssertionError("blind view unlink")
        real_unlink(path)

    with (
        patch(
            "vntyper.scripts.alignment_preflight.capture_command",
            side_effect=_write_successful_index,
        ),
        patch("vntyper.scripts.alignment_preflight.os.unlink", side_effect=reject_blind_view_unlink),
    ):
        view_path, _, binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=1
        )

    assert Path(view_path).is_symlink()
    assert Path(view_path).samefile(alignment)
    binding.close()


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

    with patch("vntyper.scripts.preflight_command_io.subprocess.run") as run:
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

    with patch("vntyper.scripts.preflight_command_io.subprocess.run", return_value=completed):
        success, output = capture_command("command", str(log_file))

    assert (success, output) == (True, "new log\n")
    assert victim.read_bytes() == b"outside victim"
    assert log_file.is_file()
    assert not log_file.is_symlink()
    assert log_file.read_text() == "new log\n"


def test_capture_command_turns_an_oserror_into_a_safe_logged_diagnostic(tmp_path: Path) -> None:
    """An OS-level spawn failure is a parseable failure, not an escaped exception."""
    log_file = tmp_path / "command.log"

    with patch("vntyper.scripts.preflight_command_io.subprocess.run", side_effect=OSError("cannot spawn")):
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


def test_an_existing_cram_csi_is_ignored_and_a_trusted_crai_is_built(tmp_path: Path) -> None:
    """A source CSI cannot prove that it belongs to this CRAM."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    output.mkdir()
    alignment = input_dir / "sample.cram"
    source_csi = input_dir / "sample.cram.csi"
    alignment.write_bytes(b"CRAM\x02")
    source_csi.write_bytes(b"CSI\x01")

    with patch(
        "vntyper.scripts.alignment_preflight.capture_command",
        side_effect=_write_successful_index,
    ) as capture:
        _, index_path, binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=1
        )

    assert Path(index_path) == output / "sample.cram.crai"
    assert Path(index_path).read_bytes() == b"current index"
    assert not Path(index_path).is_symlink()
    assert source_csi.read_bytes() == b"CSI\x01"
    capture.assert_called_once()
    binding.close()


def test_a_stale_cram_csi_cannot_override_the_selected_crai(tmp_path: Path) -> None:
    """Selecting CRAI removes a stale co-located CSI before downstream htslib calls."""
    alignment, source_crai = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    stale_source = tmp_path / "stale.cram.csi"
    stale_source.write_bytes(b"stale CSI")
    stale_view_csi = output / "sample.cram.csi"
    stale_view_csi.symlink_to(stale_source)

    with patch(
        "vntyper.scripts.alignment_preflight.capture_command",
        side_effect=_write_successful_index,
    ):
        _, index_path, binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=1
        )

    assert Path(index_path).read_bytes() == b"current index"
    assert source_crai.read_bytes() == b"patient index"
    assert not os.path.lexists(stale_view_csi)
    binding.close()


def test_an_existing_index_is_protected_while_the_view_index_is_rebuilt(tmp_path: Path) -> None:
    """Rebuilding a trusted view index never mutates the supplied source index."""
    alignment, source_index = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"

    with patch(
        "vntyper.scripts.alignment_preflight.capture_command",
        side_effect=_write_successful_index,
    ) as capture:
        view_path, index_path, binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=4
        )

    assert Path(view_path).samefile(alignment)
    assert Path(index_path).read_bytes() == b"current index"
    assert not Path(index_path).is_symlink()
    assert source_index.read_bytes() == b"patient index"
    capture.assert_called_once()
    binding.close()


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
        _, index_path, binding = build_alignment_view(
            str(alignment), str(output), "sample", "bam", _config(), threads=3
        )

    assert Path(index_path) == output / "sample.bam.bai"
    assert Path(index_path).read_bytes() == b"current index"
    assert not (input_dir / "sample.bam.bai").exists()
    assert not (input_dir / "sample.bai").exists()
    binding.close()


def test_a_locally_built_index_is_safely_rebuilt_on_a_same_input_rerun(tmp_path: Path) -> None:
    """A reserved regular index from one run must not make a sequential rerun fail."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    input_dir.mkdir()
    alignment = input_dir / "sample.cram"
    alignment.write_bytes(b"patient alignment")
    build_number = 0

    def build_index(
        command: str,
        log_file: str,
        cwd: str | None = None,
        *,
        protected_paths: tuple[str | Path, ...] = (),
    ) -> tuple[bool, str]:
        nonlocal build_number
        del log_file, cwd, protected_paths
        build_number += 1
        arguments = shlex.split(command)
        Path(arguments[arguments.index("-o") + 1]).write_bytes(f"index {build_number}".encode())
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_index):
        first_view, first_index, first_binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=2
        )
        first_binding.close()
        assert not os.path.lexists(first_view)
        second_view, second_index, second_binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=2
        )

    assert first_view == second_view
    assert first_index == second_index
    assert Path(second_view).samefile(alignment)
    assert Path(second_index).is_file()
    assert not Path(second_index).is_symlink()
    assert Path(second_index).read_bytes() == b"index 2"
    assert alignment.read_bytes() == b"patient alignment"
    assert not (input_dir / "sample.cram.crai").exists()
    assert not (input_dir / "sample.crai").exists()
    assert Path(second_view).read_bytes() == b"patient alignment"
    second_binding.close()


def test_a_local_index_is_safely_replaced_when_a_different_input_reuses_the_output(tmp_path: Path) -> None:
    """Sequential output reuse must update both view and generated index without patient-tree writes."""
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"
    output = tmp_path / "output"
    first_dir.mkdir()
    second_dir.mkdir()
    first_alignment = first_dir / "sample.cram"
    second_alignment = second_dir / "sample.cram"
    first_alignment.write_bytes(b"first patient alignment")
    second_alignment.write_bytes(b"second patient alignment")
    build_number = 0

    def build_index(
        command: str,
        log_file: str,
        cwd: str | None = None,
        *,
        protected_paths: tuple[str | Path, ...] = (),
    ) -> tuple[bool, str]:
        nonlocal build_number
        del log_file, cwd, protected_paths
        build_number += 1
        arguments = shlex.split(command)
        Path(arguments[arguments.index("-o") + 1]).write_bytes(f"index {build_number}".encode())
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_index):
        _, _, first_binding = build_alignment_view(
            str(first_alignment), str(output), "sample", "cram", _config(), threads=2
        )
        first_binding.close()
        view_path, index_path, second_binding = build_alignment_view(
            str(second_alignment), str(output), "sample", "cram", _config(), threads=2
        )

    assert Path(view_path).samefile(second_alignment)
    assert Path(index_path).is_file()
    assert not Path(index_path).is_symlink()
    assert Path(index_path).read_bytes() == b"index 2"
    assert first_alignment.read_bytes() == b"first patient alignment"
    assert second_alignment.read_bytes() == b"second patient alignment"
    assert not tuple(first_dir.glob("*.crai"))
    assert not tuple(second_dir.glob("*.crai"))
    assert Path(view_path).read_bytes() == b"second patient alignment"
    second_binding.close()


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


def test_a_view_already_pointing_at_the_same_input_gets_a_fresh_index(tmp_path: Path) -> None:
    """A correct view link remains usable while its untrusted index is replaced."""
    alignment, source_index = _indexed_cram(tmp_path / "input")
    output = tmp_path / "output"
    output.mkdir()
    view = output / "sample.cram"
    view.symlink_to(os.path.relpath(alignment, output))
    (output / "sample.cram.crai").symlink_to(os.path.relpath(source_index, output))

    with patch(
        "vntyper.scripts.alignment_preflight.capture_command",
        side_effect=_write_successful_index,
    ):
        view_path, index_path, binding = build_alignment_view(
            str(alignment), str(output), "sample", "cram", _config(), threads=2
        )

    assert Path(view_path).samefile(alignment)
    assert Path(index_path).read_bytes() == b"current index"
    assert source_index.read_bytes() == b"patient index"
    binding.close()


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

    with patch("vntyper.scripts.preflight_command_io.subprocess.run", return_value=completed):
        result = capture_command("command", str(log_file), cwd=str(tmp_path))

    assert result == (expected_success, output)
    assert log_file.read_text() == output
