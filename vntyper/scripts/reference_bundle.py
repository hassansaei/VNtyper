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
import stat
import tarfile
import tempfile
import zipfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from typing import Any

logger = logging.getLogger(__name__)

_CHUNK = 1024 * 1024

# `extractall(filter=...)` and `tarfile.FilterError` landed together in 3.12 and were
# backported to 3.10.12 and 3.11.4. `requires-python` is ">=3.10", so neither can be
# referenced unconditionally: an early 3.10 has no `filter=` keyword to pass and no
# `FilterError` class to catch. An empty tuple in an `except` clause never matches, which
# is exactly right on an interpreter where the filter cannot run in the first place.
_HAS_EXTRACTION_FILTER = "filter" in inspect.signature(tarfile.TarFile.extractall).parameters
_EXTRACT_KWARGS: dict[str, Any] = {"filter": "data"} if _HAS_EXTRACTION_FILTER else {}
_FILTER_ERRORS: tuple[type[BaseException], ...] = (
    (tarfile.FilterError,) if _HAS_EXTRACTION_FILTER else ()  # type: ignore[attr-defined,unused-ignore]
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


def _resolved_destination_root(archive: Path, destination: Path) -> Path:
    """Create and resolve an extraction root without following a root symlink."""
    if destination.is_symlink():
        message = f"{archive.name}: destination {destination} is a symlink"
        logger.error(message)
        raise ValueError(message)
    destination.mkdir(parents=True, exist_ok=True)
    return destination.resolve()


def safe_extract(archive: Path, destination: Path) -> list[str]:
    """Extract a tar archive, rejecting any member that could write outside the root.

    **Links are rejected outright rather than resolved**, and that is the whole defence
    against a class of attack that per-member target resolution does not stop. A link's
    target is judged *before* extraction, but extracting an earlier member can change
    what a later one means. Three members in order - `pivot` symlinked to `.`, `escape`
    symlinked to `pivot/..`, then a regular file `escape/owned` - each pass a lexical
    pre-check (`destination/pivot/..` normalises to `destination`, and `escape/owned`
    has neither a `..` component nor an absolute path), and the regular member is then
    written *through* `escape` into the destination's parent.

    Rejection is safe here rather than merely conservative: VNtyper's bundles are built
    by `scripts/bundle_release.py` from a fixed list of regular files and never contain
    a link, and `safe_extract` is only ever pointed at those bundles, after their digest
    has been verified against a value committed in this repository.

    A link **already present in `destination`, rather than in the archive**, is a
    separate hazard the checks above do not cover: they judge what the archive
    contains, not what the destination already has. `staged_install` seeds its staging
    directory from the existing tree with `symlinks=True`, so an `alignment` symlink
    pointing outside the reference root survives into a fresh install unchanged, and an
    ordinary, fully-verified member such as `alignment/chr1.hg19.fa` would be written
    straight through it. Every member's destination path is therefore resolved and
    confirmed to stay inside `destination` before extraction runs, which holds
    regardless of what was already on disk - and, unlike the link-target judgement
    above, is safe to do ahead of extraction: nothing in the loop above lets an archive
    member create a symlink, so no member extracted during *this* call can change what
    a later one in the same archive resolves to.

    Absolute paths and `..` components are rejected by the same loop. `filter="data"` is
    applied as defence in depth **where the interpreter has it** - `requires-python` is
    `>=3.10` and `filter=` only exists from 3.10.12, 3.11.4 and 3.12, so passing it
    unconditionally would raise `TypeError` on an early 3.10. That filter signals refusal
    with `tarfile.FilterError`, which is not a `ValueError`, so it is translated below:
    the contract this function documents has to hold on every interpreter, not only the
    ones without the filter.

    Args:
        archive: `.tar.gz` to unpack.
        destination: Directory to unpack into; created if absent.

    Returns:
        list[str]: The extracted regular-file member names, for provenance recording.

    Raises:
        ValueError: If `destination` is a symlink, on an absolute path, a `..`
            component, a symbolic or hard link, a device or FIFO member, a member whose
            destination resolves outside `destination` through a symlink already present
            there, or anything `tarfile`'s own `data` filter refuses. Per AGENTS.md the
            convention is `logger.error(message)` then `raise`, with no custom exception
            class.
    """
    destination_root = _resolved_destination_root(archive, destination)
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        for member in members:
            name = Path(member.name)
            if name.is_absolute():
                message = f"{archive.name}: absolute path in member '{member.name}'"
                logger.error(message)
                raise ValueError(message)
            if ".." in name.parts:
                message = f"{archive.name}: member '{member.name}' escapes the archive root"
                logger.error(message)
                raise ValueError(message)
            if member.issym() or member.islnk():
                message = (
                    f"{archive.name}: member '{member.name}' is a link, and a reference bundle ships "
                    "regular files only. A link's meaning depends on what has already been extracted, "
                    "so it cannot be validated ahead of extraction."
                )
                logger.error(message)
                raise ValueError(message)
            if not (member.isfile() or member.isdir()):
                message = f"{archive.name}: member '{member.name}' is not a regular file or directory"
                logger.error(message)
                raise ValueError(message)
            resolved_target = (destination / name).resolve()
            if not resolved_target.is_relative_to(destination_root):
                message = (
                    f"{archive.name}: member '{member.name}' resolves to {resolved_target}, outside "
                    f"destination root {destination_root} - a symlink already present in the destination "
                    "redirected it there"
                )
                logger.error(message)
                raise ValueError(message)
        try:
            tar.extractall(path=destination, **_EXTRACT_KWARGS)
        except _FILTER_ERRORS as refused:
            message = f"{archive.name}: tarfile's 'data' extraction filter refused the archive: {refused}"
            logger.error(message)
            raise ValueError(message) from refused
        return [member.name for member in members if member.isfile()]


def safe_extract_zip(archive: Path, destination: Path) -> list[str]:
    """Extract a ZIP archive without allowing a member to leave its root.

    The standard library sanitises absolute and parent components during ZIP
    extraction, but it follows directory symlinks that already exist in the destination.
    Staged reference installs copy such links from the prior tree, so every resolved
    member target is checked before extraction and hostile names are rejected rather than
    silently rewritten. Link entries are also refused because reference archives ship
    regular files and directories only.

    Args:
        archive: `.zip` to unpack.
        destination: Directory to unpack into; created if absent.

    Returns:
        list[str]: The extracted regular-file member names, for provenance recording.

    Raises:
        ValueError: If `destination` is a symlink, or a member is absolute, contains a
            parent component, declares a link, or resolves outside `destination` through
            a pre-existing symlink.
    """
    destination_root = _resolved_destination_root(archive, destination)
    with zipfile.ZipFile(archive, "r") as bundle:
        members = bundle.infolist()
        for member in members:
            name = PurePosixPath(member.filename)
            if name.is_absolute():
                message = f"{archive.name}: absolute path in member '{member.filename}'"
                logger.error(message)
                raise ValueError(message)
            if ".." in name.parts:
                message = f"{archive.name}: member '{member.filename}' escapes the archive root"
                logger.error(message)
                raise ValueError(message)
            if stat.S_ISLNK(member.external_attr >> 16):
                message = (
                    f"{archive.name}: member '{member.filename}' is a link, and a reference archive "
                    "ships regular files only"
                )
                logger.error(message)
                raise ValueError(message)
            resolved_target = (destination / member.filename).resolve()
            if not resolved_target.is_relative_to(destination_root):
                message = (
                    f"{archive.name}: member '{member.filename}' resolves to {resolved_target}, outside "
                    f"destination root {destination_root} - a symlink already present in the destination "
                    "redirected it there"
                )
                logger.error(message)
                raise ValueError(message)
        bundle.extractall(path=destination)
        return [member.filename for member in members if not member.is_dir()]


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
    # A plain flag rather than blanking `staging`: once the rename has happened that path
    # *is* the live installation, and the cleanup below must not be able to delete it.
    activated = False
    try:
        # Seeding only reads `target`, but a mid-copy failure (disk full, permission
        # error, an interrupted copy of a large reference tree) must still hit the
        # `except` below - otherwise every failed retry leaks another staging directory.
        if seed_from_existing and target.exists():
            shutil.copytree(target, staging, dirs_exist_ok=True, symlinks=True)
        yield staging
        if target.exists():
            previous = Path(tempfile.mkdtemp(prefix=f".{target.name}.previous.", dir=target.parent))
            previous.rmdir()
            target.rename(previous)
        staging.rename(target)
        activated = True
    except BaseException:
        if not activated:
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
            # Only report a preserved path that actually exists: if `target.rename(previous)`
            # itself is what failed, `previous.rmdir()` already vacated the reserved name
            # and nothing was ever moved there - `target` still holds it safely, and
            # reporting the vacant `previous` path would send an operator hunting for
            # data that was never displaced.
            if previous is not None and previous.exists():
                logger.error(f"previous reference tree preserved at {previous}; restore it by hand")
        raise
    else:
        # Activation is confirmed; only now is the old tree disposable.
        if previous is not None:
            shutil.rmtree(previous, ignore_errors=True)
