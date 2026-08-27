"""Atomic publication of derived adVNTR result tables."""

from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def invalidate_advntr_artifact(destination: str | Path) -> None:
    """Remove an adVNTR artifact before beginning a new attempt.

    Args:
        destination: Raw or derived path whose bytes belong only to the previous attempt.

    Raises:
        RuntimeError: If an existing result cannot be removed.
    """
    artifact_path = Path(destination)
    try:
        artifact_path.unlink(missing_ok=True)
    except OSError as error:
        message = f"Could not invalidate stale adVNTR artifact {artifact_path}: {error}"
        logger.error(message)
        raise RuntimeError(message) from error


def publish_advntr_result(frame: pd.DataFrame, destination: str | Path) -> None:
    """Atomically publish one complete adVNTR result table.

    Args:
        frame: Final result rows in their public column order.
        destination: Public TSV path to replace only after the candidate is complete.

    Raises:
        RuntimeError: If the candidate cannot be written, flushed, or installed.
    """
    result_path = Path(destination)
    temporary_path: Path | None = None
    published = False
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=result_path.parent,
            prefix=f".{result_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as candidate:
            temporary_path = Path(candidate.name)
            os.fchmod(candidate.fileno(), 0o644)
            frame.to_csv(candidate, sep="\t", index=False)
            candidate.flush()
            os.fsync(candidate.fileno())
        os.replace(temporary_path, result_path)
        published = True
    except Exception as error:
        message = f"Failed to publish adVNTR result {result_path}: {error}"
        logger.error(message)
        raise RuntimeError(message) from error
    finally:
        if temporary_path is not None and not published:
            try:
                temporary_path.unlink()
            except OSError as cleanup_error:
                cleanup_message = (
                    f"Failed to remove incomplete adVNTR result candidate {temporary_path}: {cleanup_error}"
                )
                logger.error(cleanup_message)
