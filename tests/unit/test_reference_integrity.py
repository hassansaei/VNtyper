"""Reference integrity decisions: retry only stale reused bytes, and only once."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from vntyper.scripts.reference_integrity import fetch_verified_asset, verify_existing_asset

pytestmark = pytest.mark.unit

URL = "https://example.invalid/reference.fa"
GOOD = b"pinned reference bytes"
DIGEST = hashlib.sha256(GOOD).hexdigest()


def test_digest_correct_preexisting_bytes_are_reused_unchanged(tmp_path: Path) -> None:
    target = tmp_path / "reference.fa"
    target.write_bytes(GOOD)
    before = target.read_bytes()
    attempts = 0

    def download(url: str, destination: Path) -> bool:
        nonlocal attempts
        attempts += 1
        return False

    assert fetch_verified_asset(target.name, URL, target, DIGEST, download) == DIGEST
    assert attempts == 1
    assert target.read_bytes() == before


def test_a_stale_preexisting_asset_is_replaced_exactly_once(tmp_path: Path) -> None:
    target = tmp_path / "reference.fa"
    target.write_bytes(b"stale reused bytes")
    attempts = 0

    def download(url: str, destination: Path) -> bool:
        nonlocal attempts
        attempts += 1
        if destination.exists():
            return False
        destination.write_bytes(GOOD)
        return True

    assert fetch_verified_asset(target.name, URL, target, DIGEST, download) == DIGEST
    assert attempts == 2
    assert target.read_bytes() == GOOD


def test_a_fresh_checksum_mismatch_fails_without_retry_or_retry_promise(tmp_path: Path) -> None:
    target = tmp_path / "reference.fa"
    attempts = 0

    def download(url: str, destination: Path) -> bool:
        nonlocal attempts
        attempts += 1
        destination.write_bytes(b"wrong upstream bytes")
        return True

    with pytest.raises(ValueError) as excinfo:
        fetch_verified_asset(target.name, URL, target, DIGEST, download)

    assert attempts == 1
    assert not target.exists()
    assert "freshly downloaded bytes do not match" in str(excinfo.value)
    assert "retry" not in str(excinfo.value).lower()


def test_a_failed_initial_download_removes_partial_bytes(tmp_path: Path) -> None:
    target = tmp_path / "reference.fa"

    def download(url: str, destination: Path) -> bool:
        destination.write_bytes(b"partial initial download")
        raise RuntimeError("connection reset")

    with pytest.raises(RuntimeError, match="connection reset"):
        fetch_verified_asset(target.name, URL, target, DIGEST, download)

    assert not target.exists()


def test_a_failed_replacement_download_removes_partial_bytes(tmp_path: Path) -> None:
    """A transport failure during the sole repair must not leave bytes for expansion."""
    target = tmp_path / "reference.fa"
    target.write_bytes(b"stale reused bytes")
    attempts = 0

    def download(url: str, destination: Path) -> bool:
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            return False
        destination.write_bytes(b"partial replacement")
        raise RuntimeError("connection reset")

    with pytest.raises(RuntimeError, match="connection reset"):
        fetch_verified_asset(target.name, URL, target, DIGEST, download)

    assert attempts == 2
    assert not target.exists(), "failed replacement bytes must not survive for a consumer to use"


def test_a_replacement_checksum_mismatch_is_removed_without_a_second_retry(tmp_path: Path) -> None:
    target = tmp_path / "reference.fa"
    target.write_bytes(b"stale reused bytes")
    attempts = 0

    def download(url: str, destination: Path) -> bool:
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            return False
        destination.write_bytes(b"wrong replacement bytes")
        return True

    with pytest.raises(ValueError) as excinfo:
        fetch_verified_asset(target.name, URL, target, DIGEST, download)

    assert attempts == 2
    assert not target.exists()
    assert "re-downloaded bytes do not match" in str(excinfo.value)


def test_a_downloader_that_refuses_the_replacement_fails_closed(tmp_path: Path) -> None:
    target = tmp_path / "reference.fa"
    target.write_bytes(b"stale reused bytes")
    attempts = 0

    def download(url: str, destination: Path) -> bool:
        nonlocal attempts
        attempts += 1
        if attempts == 2:
            destination.write_bytes(b"unverified race bytes")
        return False

    with pytest.raises(RuntimeError, match="no transfer occurred"):
        fetch_verified_asset(target.name, URL, target, DIGEST, download)

    assert attempts == 2
    assert not target.exists()


def test_failed_cleanup_does_not_claim_the_replacement_made_no_transfer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    target = tmp_path / "reference.fa"
    target.write_bytes(b"stale reused bytes")
    attempts = 0
    original_unlink = Path.unlink

    def download(url: str, destination: Path) -> bool:
        nonlocal attempts
        attempts += 1
        if attempts == 2:
            destination.write_bytes(b"unverified race bytes")
        return False

    def fail_for_unverified_replacement(path: Path, missing_ok: bool = False) -> None:
        if path == target and path.exists() and path.read_bytes() == b"unverified race bytes":
            raise OSError("read-only filesystem")
        original_unlink(path, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", fail_for_unverified_replacement)

    with pytest.raises(OSError, match="read-only filesystem"):
        fetch_verified_asset(target.name, URL, target, DIGEST, download)

    assert target.read_bytes() == b"unverified race bytes"
    assert "no transfer occurred" not in caplog.text


def test_verify_existing_asset_removes_mismatching_bytes(tmp_path: Path) -> None:
    target = tmp_path / "reference.fa"
    target.write_bytes(b"wrong")

    with pytest.raises(ValueError, match="removed mismatched reference.fa"):
        verify_existing_asset(target.name, target, DIGEST)

    assert not target.exists()
