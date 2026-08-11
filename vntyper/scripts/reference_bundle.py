"""Fetch, verify and install reference bundles without ever leaving a partial tree.

`install_references.py` orchestrates; this module owns the parts that must not be got
wrong: digest verification against a value committed in this repository, archive
extraction that cannot write outside its root, and an activation that is atomic from
the point of view of the next run.
"""

from __future__ import annotations

import hashlib
import inspect
import logging
import shutil
import tarfile
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_CHUNK = 1024 * 1024

# `extractall(filter=...)` landed in 3.12 and was backported to 3.10.12 and 3.11.4.
# `requires-python` is ">=3.10", so it cannot be passed unconditionally.
_EXTRACT_KWARGS: dict[str, Any] = (
    {"filter": "data"} if "filter" in inspect.signature(tarfile.TarFile.extractall).parameters else {}
)


def sha256_of(path: Path) -> str:
    """Return the hex SHA-256 of a file, read in chunks.

    Args:
        path: File to digest.

    Returns:
        str: Lowercase hex digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str) -> None:
    """Fail closed unless a file matches the digest committed in this repository.

    The release's own `SHA256SUMS` is co-hosted with the assets it describes, so it
    corroborates but cannot be the root of trust. The expected value comes from
    `install_references_config.json`, which is a base-image content-hash input - so
    changing any reference byte necessarily changes the base tag too.

    Args:
        path: File to check.
        expected: Lowercase hex SHA-256.

    Raises:
        ValueError: If the digests differ.
    """
    actual = sha256_of(path)
    if actual != expected:
        message = f"checksum mismatch for {path.name}: expected {expected}, got {actual}"
        logger.error(message)
        raise ValueError(message)
    logger.info(f"  ✓ verified {path.name}")


def _is_within(root: Path, candidate: Path) -> bool:
    """True when `candidate` resolves inside `root`."""
    try:
        candidate.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def safe_extract(archive: Path, destination: Path) -> None:
    """Extract a tar archive, rejecting any member that could write outside the root.

    The explicit member loop below is the guarantee, not `tarfile`'s `filter=`
    argument: `requires-python` is `>=3.10` and `filter=` only exists from 3.10.12,
    3.11.4 and 3.12, so passing it unconditionally would raise `TypeError` on an
    early 3.10. It is applied as defence in depth where available.

    Args:
        archive: `.tar.gz` to unpack.
        destination: Directory to unpack into; created if absent.

    Raises:
        ValueError: On an absolute path, a `..` component, a device or FIFO member,
            or a link resolving outside the destination. Per AGENTS.md the convention is
            `logger.error(message)` then `raise`, with no custom exception class.
    """
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            name = Path(member.name)
            if name.is_absolute():
                message = f"{archive.name}: absolute path in member '{member.name}'"
                logger.error(message)
                raise ValueError(message)
            if ".." in name.parts:
                message = f"{archive.name}: member '{member.name}' escapes the archive root"
                logger.error(message)
                raise ValueError(message)
            if not (member.isfile() or member.isdir() or member.issym() or member.islnk()):
                message = f"{archive.name}: member '{member.name}' is not a regular file or link"
                logger.error(message)
                raise ValueError(message)
            # Symlink targets are relative to the member's own directory; hard-link
            # targets are relative to the archive root. Resolving both the same way
            # misjudges one of them, and on a Python without `filter=` this loop is the
            # only protection there is.
            if member.issym():
                linked = destination / name.parent / member.linkname
            elif member.islnk():
                linked = destination / member.linkname
            else:
                linked = None
            if linked is not None and not _is_within(destination, linked):
                message = f"{archive.name}: link '{member.name}' escapes the archive root"
                logger.error(message)
                raise ValueError(message)
        tar.extractall(path=destination, **_EXTRACT_KWARGS)


@contextmanager
def staged_install(target: Path, *, seed_from_existing: bool = True) -> Iterator[Path]:
    """Build a reference tree beside its destination and activate it only on success.

    The staging directory is a sibling of `target`, so activation is a rename on the
    same filesystem. On any exception the staging directory is removed and any
    pre-existing installation is restored, or - if it cannot be restored - preserved
    under a named path and reported. The previous tree is deleted only after activation
    has succeeded, so no failure path can leave the caller with nothing.

    `seed_from_existing` copies the current tree into the staging directory first, so an
    install of one assembly does not erase another. Without it,
    `install-references --references hg38` after `--references hg19` would silently
    delete hg19, and installing into the tracked `reference/` directory would delete the
    retained `README.md` and `pseudonymize*` files.

    Args:
        target: Final directory.
        seed_from_existing: Copy the existing tree into staging before yielding, so the
            install merges rather than replaces. False only for a deliberate clean install.

    Yields:
        Path: The staging directory to populate.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}.staging.", dir=target.parent))
    previous: Path | None = None
    if seed_from_existing and target.exists():
        shutil.copytree(target, staging, dirs_exist_ok=True, symlinks=True)
    try:
        yield staging
        if target.exists():
            previous = Path(tempfile.mkdtemp(prefix=f".{target.name}.previous.", dir=target.parent))
            previous.rmdir()
            target.rename(previous)
        staging.rename(target)
        staging = None  # type: ignore[assignment]
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        if previous is not None:
            # Restore only into a vacant slot, and NEVER delete `previous` on this path:
            # if something else recreated `target`, or the rename back fails, the previous
            # installation is the only surviving copy. Leaving it on disk under a named
            # path is recoverable; deleting it is not.
            if not target.exists():
                try:
                    previous.rename(target)
                    previous = None
                except OSError as restore_error:
                    logger.error(f"could not restore {target} from {previous}: {restore_error}")
            if previous is not None:
                logger.error(f"previous reference tree preserved at {previous}; restore it by hand")
        raise
    else:
        # Activation is confirmed; only now is the old tree disposable.
        if previous is not None:
            shutil.rmtree(previous, ignore_errors=True)
