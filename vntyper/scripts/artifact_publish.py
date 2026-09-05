"""Atomic installation and cleanup helpers for generated pipeline artifacts (#314).

Public BAM and BAI files produced by subprocesses are written to temporary, deterministic
partial paths (``<final>.partial``) and atomically moved into their permanent destination
via ``os.replace`` only upon verified success. This prevents interrupted or killed
processes from leaving corrupt, truncated, or header-only BAM files under public names.
"""

from __future__ import annotations

import contextlib
import logging
import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger(__name__)

PARTIAL_SUFFIX = ".partial"


def partial_path(final: str | Path) -> Path:
    """Return the deterministic sibling partial path for a target file.

    Args:
        final: The target final destination path.

    Returns:
        Path: Sibling path ending with ``.partial``.
    """
    path = Path(final)
    return path.with_name(f"{path.name}{PARTIAL_SUFFIX}")


def publish_partial(partial: str | Path, final: str | Path) -> Path:
    """Atomically install a completed partial file to its final destination.

    Uses ``os.replace`` to ensure atomic installation across POSIX filesystems within
    the same directory.

    Args:
        partial: Path to the partial artifact.
        final: Path to the permanent destination.

    Returns:
        Path: The resolved destination path.

    Raises:
        FileNotFoundError: If the partial file does not exist.
        OSError: If replacing fails.
    """
    partial_path_obj = Path(partial)
    final_path_obj = Path(final)
    if not partial_path_obj.exists():
        raise FileNotFoundError(f"Partial artifact not found: {partial_path_obj}")
    os.replace(partial_path_obj, final_path_obj)
    logger.debug("Published %s -> %s", partial_path_obj, final_path_obj)
    return final_path_obj


def discard_partial(partial: str | Path | None) -> None:
    """Remove a partial artifact if it exists, suppressing missing file errors.

    Args:
        partial: Path to the partial artifact, or None.
    """
    if partial is None:
        return
    path = Path(partial)
    with contextlib.suppress(FileNotFoundError, OSError):
        if path.is_file() or path.is_symlink():
            path.unlink()
            logger.debug("Discarded partial artifact: %s", path)


@contextmanager
def partial_output(final: str | Path, *, check_non_empty: bool = False) -> Iterator[Path]:
    """Context manager yielding a temporary partial path, publishing on success.

    Before yielding, any existing stale partial artifact at the target path is
    removed. On normal block exit, the partial file is checked (and optionally
    verified to be non-empty) and published to ``final`` via :func:`publish_partial`.
    If an exception occurs within the block, the partial file is unlinked via
    :func:`discard_partial` and the exception is re-raised.

    Args:
        final: Destination file to create atomically.
        check_non_empty: If True, raise OSError if the partial file has size 0.

    Yields:
        Path: The deterministic partial path to write to.

    Raises:
        FileNotFoundError: If the block completes without creating the partial file.
        OSError: If ``check_non_empty`` is True and the partial file has size 0.
    """
    target = Path(final)
    partial = partial_path(target)
    discard_partial(partial)
    try:
        yield partial
        if not partial.exists():
            raise FileNotFoundError(f"Expected partial artifact not created: {partial}")
        if check_non_empty and partial.stat().st_size == 0:
            raise OSError(f"Generated artifact is empty: {partial}")
        publish_partial(partial, target)
    except BaseException:
        discard_partial(partial)
        raise
