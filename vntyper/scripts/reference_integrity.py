"""Digest-gated reuse and one-shot repair of downloaded reference assets."""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path

from vntyper.scripts.reference_bundle import verify_sha256

logger = logging.getLogger(__name__)

DownloadFile = Callable[[str, Path], bool]


def _download_or_discard(url: str, target_path: Path, download: DownloadFile) -> bool:
    """Run an atomic downloader defensively, discarding bytes after a failed transfer."""
    try:
        return download(url, target_path)
    except RuntimeError:
        target_path.unlink(missing_ok=True)
        raise


def verify_existing_asset(name: str, target_path: Path, expected_sha256: str) -> str:
    """Verify one asset and remove bytes that do not match their pinned digest.

    Args:
        name: Operator-facing asset name.
        target_path: File to verify.
        expected_sha256: Pinned SHA-256 digest.

    Returns:
        str: ``expected_sha256`` after successful verification.

    Raises:
        ValueError: If the file does not match. The mismatching bytes are removed first.
    """
    try:
        verify_sha256(target_path, expected_sha256)
    except ValueError as mismatch:
        target_path.unlink(missing_ok=True)
        message = f"{mismatch}; removed mismatched {name}"
        logger.error(message)
        raise ValueError(message) from mismatch
    return expected_sha256


def fetch_verified_asset(
    name: str,
    url: str,
    target_path: Path,
    expected_sha256: str,
    download: DownloadFile,
) -> str:
    """Download or reuse an asset, repairing one stale pre-existing copy in-run.

    Digest configuration must be resolved before this function is called. That keeps
    missing or contradictory pins distinct from byte mismatches: only a checksum
    mismatch on bytes reused from before this call is retried. A fresh mismatch fails
    immediately, and a failed replacement leaves no unverified bytes behind.

    Args:
        name: Operator-facing asset name.
        url: Source URL passed to ``download``.
        target_path: Final local path for the asset.
        expected_sha256: Already-resolved pinned SHA-256 digest.
        download: Atomic downloader returning True for a fresh transfer and False for reuse.

    Returns:
        str: ``expected_sha256`` after successful verification.

    Raises:
        RuntimeError: If either download fails.
        ValueError: If freshly downloaded or replacement bytes miss the pinned digest.
    """
    fetched = _download_or_discard(url, target_path, download)
    try:
        return verify_existing_asset(name, target_path, expected_sha256)
    except ValueError as mismatch:
        if fetched:
            message = f"{mismatch}; freshly downloaded bytes do not match the pinned digest"
            logger.error(message)
            raise ValueError(message) from mismatch

        logger.warning(f"Pre-existing {name} is stale; downloading it afresh")
        refetched = _download_or_discard(url, target_path, download)
        if not refetched:
            target_path.unlink(missing_ok=True)
            message = f"Failed to re-download stale {name}: downloader reported that no transfer occurred"
            logger.error(message)
            raise RuntimeError(message) from mismatch
        try:
            return verify_existing_asset(name, target_path, expected_sha256)
        except ValueError as refetched_mismatch:
            message = f"{refetched_mismatch}; re-downloaded bytes do not match the pinned digest"
            logger.error(message)
            raise ValueError(message) from refetched_mismatch
