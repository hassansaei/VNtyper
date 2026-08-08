"""Captured-command I/O tests for alignment preflight."""

from __future__ import annotations

import signal
import subprocess
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from vntyper.scripts import alignment_preflight, preflight_command_io

pytestmark = pytest.mark.unit


def test_alignment_preflight_reexports_the_extracted_capture_command() -> None:
    """Existing imports keep working after command execution is extracted."""
    assert alignment_preflight.capture_command is preflight_command_io.capture_command


def test_capture_command_timeout_kills_and_reaps_the_shell_process_group(tmp_path: Path) -> None:
    """A timed-out probe cannot leave samtools or a pipeline child running."""
    process = Mock(pid=4321, returncode=-signal.SIGKILL)
    process.communicate.side_effect = [
        subprocess.TimeoutExpired("probe", 0.5, output="partial output"),
        ("partial output", None),
    ]
    log = tmp_path / "probe.log"

    with (
        patch("vntyper.scripts.preflight_command_io.subprocess.Popen", return_value=process) as popen,
        patch("vntyper.scripts.preflight_command_io.os.killpg") as killpg,
    ):
        success, output = preflight_command_io.capture_command("probe", str(log), timeout_seconds=0.5)

    assert success is False
    assert "timed out after 0.5 seconds" in output
    assert log.read_text() == output
    assert popen.call_args.kwargs["start_new_session"] is True
    killpg.assert_called_once_with(4321, signal.SIGKILL)
    assert process.communicate.call_count == 2


def test_capture_command_timeout_reaps_when_process_group_kill_reports_an_os_error(tmp_path: Path) -> None:
    """A failed group-kill attempt must not skip the mandatory child reap."""
    process = Mock(pid=4321, returncode=-signal.SIGKILL)
    process.communicate.side_effect = [
        subprocess.TimeoutExpired("probe", 0.5),
        ("terminated output", None),
    ]
    log = tmp_path / "probe.log"

    with (
        patch("vntyper.scripts.preflight_command_io.subprocess.Popen", return_value=process),
        patch("vntyper.scripts.preflight_command_io.os.killpg", side_effect=PermissionError("not permitted")),
    ):
        success, output = preflight_command_io.capture_command("probe", str(log), timeout_seconds=0.5)

    assert success is False
    assert "timed out after 0.5 seconds" in output
    assert process.communicate.call_count == 2
    process.kill.assert_called_once_with()
