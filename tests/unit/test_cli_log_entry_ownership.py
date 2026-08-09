"""Pipeline logging must establish exclusive ownership before opening its file."""

from __future__ import annotations

import os
from pathlib import Path
from unittest import mock

import pytest

from vntyper import cli

pytestmark = pytest.mark.unit


def _config(reference: Path) -> dict[str, object]:
    """Build the configuration surface the CLI reads before dispatch."""
    return {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19"},
        "reference_data": {"bwa_reference_hg19": str(reference)},
    }


@pytest.mark.parametrize("output_route", ["direct", "symlinked"])
@pytest.mark.parametrize("entry_kind", ["symlink", "hardlink", "fifo"])
def test_default_pipeline_log_rejects_an_unowned_unsafe_entry_before_any_setup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    output_route: str,
    entry_kind: str,
) -> None:
    """A planted log entry is unsafe even when it is unrelated to every input."""
    actual_output = tmp_path / "actual-output"
    actual_output.mkdir()
    output = actual_output
    if output_route == "symlinked":
        output = tmp_path / "output-link"
        output.symlink_to(actual_output, target_is_directory=True)
    log_file = actual_output / "pipeline.log"
    external = tmp_path / "unrelated-external.log"
    if entry_kind == "fifo":
        os.mkfifo(log_file)
        entry_inode = os.lstat(log_file).st_ino
        protected_inode = None
        protected_bytes = None
    else:
        external.write_bytes(b"unrelated-external-bytes")
        if entry_kind == "symlink":
            log_file.symlink_to(external)
        else:
            log_file.hardlink_to(external)
        entry_inode = os.lstat(log_file).st_ino
        protected_inode = external.stat().st_ino
        protected_bytes = external.read_bytes()
    fastq = tmp_path / "reads.fastq.gz"
    fastq.write_bytes(b"operator-fastq")
    setup = mock.Mock()
    handler = mock.Mock()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: _config(tmp_path / "reference.fa"))
    monkeypatch.setattr(cli, "setup_logging", setup)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handler)

    with mock.patch.object(Path, "mkdir", autospec=True) as mkdir, pytest.raises(SystemExit) as raised:
        cli.main(["pipeline", "-o", str(output), "--fastq1", str(fastq)])

    assert raised.value.code == 1
    assert os.lstat(log_file).st_ino == entry_inode
    if protected_bytes is not None:
        assert protected_inode is not None
        assert external.stat().st_ino == protected_inode
        assert external.read_bytes() == protected_bytes
    mkdir.assert_not_called()
    setup.assert_not_called()
    handler.assert_not_called()


@pytest.mark.parametrize("output_route", ["direct", "symlinked"])
def test_default_pipeline_log_accepts_a_single_link_regular_rerun(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    output_route: str,
) -> None:
    """A regular rerun log remains valid through direct and symlinked output roots."""
    actual_output = tmp_path / "actual-output"
    actual_output.mkdir()
    log_file = actual_output / "pipeline.log"
    log_file.write_bytes(b"previous-run-log")
    original_inode = log_file.stat().st_ino
    output = actual_output
    if output_route == "symlinked":
        output = tmp_path / "output-link"
        output.symlink_to(actual_output, target_is_directory=True)
    fastq = tmp_path / "reads.fastq.gz"
    fastq.write_bytes(b"operator-fastq")
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: _config(tmp_path / "reference.fa"))
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(["pipeline", "-o", str(output), "--fastq1", str(fastq)])

    assert events == ["setup", "handler"]
    assert log_file.stat().st_ino == original_inode
    assert log_file.read_bytes() == b"previous-run-log"
