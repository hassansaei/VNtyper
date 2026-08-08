"""Validation of the preflight-error artifact before it crosses into Redis."""

from __future__ import annotations

import json
import os
import socket
import threading
import time
from pathlib import Path

import pytest
from app.job_failures import read_preflight_failure

pytestmark = pytest.mark.unit


def _artifact(message: str = "Unable to resolve CRAM reference: contig=chr1, M5=digest.") -> dict:
    """Return the complete on-disk preflight failure contract."""
    return {
        "code": "reference_unresolved",
        "message": message,
        "candidates": [["cli", "full.fa", "probe exited non-zero"]],
    }


def test_reader_returns_only_code_and_message_from_the_exact_artifact(tmp_path: Path) -> None:
    """Candidate detail remains on disk; only curated fields cross into Redis."""
    (tmp_path / "preflight_error.json").write_text(json.dumps(_artifact()), encoding="utf-8")

    assert read_preflight_failure(tmp_path) == {
        "code": "reference_unresolved",
        "message": "Unable to resolve CRAM reference: contig=chr1, M5=digest.",
    }


@pytest.mark.parametrize(
    "payload",
    [
        {"code": "reference_unresolved", "message": "missing candidates"},
        _artifact("Cannot decode /opt/vntyper/private/sample.cram"),
        _artifact(r"Cannot decode C:\worker\private\sample.cram"),
    ],
)
def test_reader_rejects_malformed_or_path_bearing_artifacts(tmp_path: Path, payload: dict) -> None:
    """An unauthenticated endpoint never trusts an incomplete or path-bearing file."""
    (tmp_path / "preflight_error.json").write_text(json.dumps(payload), encoding="utf-8")

    assert read_preflight_failure(tmp_path) is None


def test_reader_does_not_follow_a_preflight_error_symlink(tmp_path: Path) -> None:
    """A pipeline-created filename cannot be redirected to an arbitrary worker file."""
    outside = tmp_path / "private.json"
    outside.write_text(json.dumps(_artifact()), encoding="utf-8")
    output = tmp_path / "output"
    output.mkdir()
    (output / "preflight_error.json").symlink_to(outside)

    assert read_preflight_failure(output) is None


@pytest.mark.parametrize("entry_kind", ["fifo", "socket", "directory", "hardlink", "symlink"])
def test_reader_rejects_every_non_single_link_regular_artifact_without_blocking(
    tmp_path: Path, entry_kind: str
) -> None:
    """Special files and aliases cannot block the worker or expose another inode."""
    output = tmp_path / "output"
    output.mkdir()
    artifact = output / "preflight_error.json"
    protected = tmp_path / "protected.json"
    protected.write_text(json.dumps(_artifact()), encoding="utf-8")
    unblocker: threading.Thread | None = None
    bound_socket: socket.socket | None = None
    if entry_kind == "fifo":
        os.mkfifo(artifact)

        def _release_blocking_reader() -> None:
            time.sleep(0.25)
            descriptor = os.open(artifact, os.O_RDWR | os.O_NONBLOCK)
            os.close(descriptor)

        unblocker = threading.Thread(target=_release_blocking_reader)
        unblocker.start()
    elif entry_kind == "socket":
        bound_socket = socket.socket(socket.AF_UNIX)
        bound_socket.bind(str(artifact))
    elif entry_kind == "directory":
        artifact.mkdir()
    elif entry_kind == "hardlink":
        os.link(protected, artifact)
    else:
        artifact.symlink_to(protected)

    started = time.monotonic()
    try:
        result = read_preflight_failure(output)
        elapsed = time.monotonic() - started
    finally:
        if bound_socket is not None:
            bound_socket.close()
        if unblocker is not None:
            unblocker.join(timeout=1)

    assert result is None
    assert elapsed < 0.15


def test_reader_rejects_a_symlinked_output_directory(tmp_path: Path) -> None:
    """Directory traversal is refused before opening the fixed artifact name."""
    protected = tmp_path / "protected"
    protected.mkdir()
    (protected / "preflight_error.json").write_text(json.dumps(_artifact()), encoding="utf-8")
    output = tmp_path / "output"
    output.symlink_to(protected, target_is_directory=True)

    assert read_preflight_failure(output) is None
