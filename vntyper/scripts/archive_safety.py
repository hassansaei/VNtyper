"""Fail-closed, atomic archive creation for VNtyper result trees."""

import logging
import os
import shutil
import stat
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)

_ARCHIVE_SUFFIXES = {"zip": ".zip", "gztar": ".tar.gz"}


def create_safe_archive(base_name: str | Path, archive_format: str, root_dir: str | Path) -> str:
    """Create an archive only when every source entry is owned regular data.

    The complete tree is checked with ``lstat`` before the archiver can open a
    member. Symbolic links, hard-linked files and special filesystem entries
    are rejected because each can refer to bytes that the result directory
    does not exclusively own. The archive is built under a temporary sibling
    directory and installed with ``os.replace`` only after creation succeeds.

    Args:
        base_name: Destination path without its format suffix.
        archive_format: ``zip`` or shutil's ``gztar`` name.
        root_dir: Result directory whose contents should be archived.

    Returns:
        str: The installed archive path.

    Raises:
        ValueError: If the format, root, or any entry is unsafe.
        OSError: If validation or archive creation cannot access the filesystem.
    """
    if archive_format not in _ARCHIVE_SUFFIXES:
        msg = f"Unsupported archive format: {archive_format}"
        logger.error(msg)
        raise ValueError(msg)

    root = Path(root_dir)
    root_metadata = root.lstat()
    if stat.S_ISLNK(root_metadata.st_mode):
        msg = "Refusing to archive unsafe symbolic link used as result root."
        logger.error(msg)
        raise ValueError(msg)
    if not stat.S_ISDIR(root_metadata.st_mode):
        msg = f"Archive root is not a directory: {root}"
        logger.error(msg)
        raise ValueError(msg)

    for directory, directory_names, file_names in os.walk(root, topdown=True, followlinks=False):
        for entry_name in (*directory_names, *file_names):
            entry = Path(directory, entry_name)
            metadata = entry.lstat()
            relative_entry = entry.relative_to(root)
            if stat.S_ISLNK(metadata.st_mode):
                msg = f"Refusing to archive unsafe symbolic link '{relative_entry}'."
                logger.error(msg)
                raise ValueError(msg)
            if stat.S_ISDIR(metadata.st_mode):
                continue
            if stat.S_ISREG(metadata.st_mode):
                if metadata.st_nlink > 1:
                    msg = f"Refusing to archive unsafe hard-linked file '{relative_entry}'."
                    logger.error(msg)
                    raise ValueError(msg)
                continue
            msg = f"Refusing to archive unsupported filesystem entry '{relative_entry}'."
            logger.error(msg)
            raise ValueError(msg)

    archive_base = Path(base_name)
    archive_path = Path(f"{archive_base}{_ARCHIVE_SUFFIXES[archive_format]}")
    temporary_dir = Path(tempfile.mkdtemp(prefix=f".{archive_base.name}.archive-", dir=str(archive_base.parent)))
    try:
        temporary_archive = shutil.make_archive(
            base_name=str(temporary_dir / "payload"),
            format=archive_format,
            root_dir=str(root),
            base_dir=".",
        )
        os.replace(temporary_archive, archive_path)
    finally:
        shutil.rmtree(temporary_dir, ignore_errors=True)

    return str(archive_path)
