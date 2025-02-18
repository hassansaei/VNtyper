#!/usr/bin/env python3
# tests/unit/test_utils.py

"""
Unit tests for utility functions.
Includes testing for command execution and BAM file validation.
"""

import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
from vntyper.scripts.utils import run_command, validate_bam_file


def test_run_command_success(tmp_path):
    """
    Test successful execution of a shell command.
    """
    log_file = tmp_path / "cmd.log"
    cmd = "echo 'Hello test_run_command'"
    result = run_command(cmd, str(log_file))
    assert result is True, "Expected run_command to succeed."
    assert "Hello test_run_command" in log_file.read_text()


@patch("subprocess.Popen")
def test_run_command_failure(mock_popen, tmp_path):
    """
    Test failure scenario for a shell command execution.
    """
    log_file = tmp_path / "fail.log"
    process_mock = MagicMock()
    process_mock.stdout = [b"Simulated error\n"]
    process_mock.wait.return_value = 1
    process_mock.returncode = 1
    mock_popen.return_value = process_mock

    cmd = "bad_command"
    ret = run_command(cmd, str(log_file))
    assert not ret, "Expected run_command to fail."
    assert "Simulated error" in log_file.read_text()


def test_validate_bam_file_success(tmp_path, test_config):
    """
    Test validation of a valid BAM file.
    We'll create a local file and mock out run_command so samtools quickcheck passes.
    """
    bam_file = tmp_path / "temp_valid.bam"
    bam_file.touch()

    # Mock out run_command to simulate a passing samtools quickcheck.
    with patch("vntyper.scripts.utils.run_command", return_value=True):
        validate_bam_file(str(bam_file))  # Should not raise ValueError.


def test_validate_bam_file_nonexistent():
    """
    Test behavior when a nonexistent BAM file is validated.
    """
    with pytest.raises(ValueError) as exc:
        validate_bam_file("nonexistent.bam")
    assert "does not exist" in str(
        exc.value
    ), "Expected ValueError for nonexistent file."
