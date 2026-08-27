"""Bounded, atomic transport for reference assets."""

from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path
from urllib.request import urlopen

logger = logging.getLogger(__name__)

#: Socket-level ceiling for one reference download. Reference genomes are large, but
#: a stalled mirror must not hang an installation indefinitely.
DOWNLOAD_TIMEOUT_SECONDS = 300

#: Bound each response read so multi-gigabyte references are never buffered in memory.
DOWNLOAD_CHUNK_SIZE = 1024 * 1024


def download_file(url: str, dest_path: Path) -> bool:
    """Download one reference asset through a partial file.

    The transfer lands only after its response reaches EOF. Existing destinations are
    deliberately reused because callers verify their committed digest separately.

    Args:
        url: URL to download.
        dest_path: Final local path for the asset.

    Returns:
        True when the asset was fetched; False when ``dest_path`` already existed.

    Raises:
        RuntimeError: If opening, streaming, or landing the download fails.
    """
    if dest_path.exists():
        logger.info(f"File already exists at {dest_path}. Skipping download.")
        return False

    logger.info(f"Downloading from {url} to {dest_path}...")
    partial_path: Path | None = None
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=dest_path.parent,
            prefix=f".{dest_path.name}.",
            suffix=".partial",
            delete=False,
        ) as sink:
            partial_path = Path(sink.name)
            os.fchmod(sink.fileno(), 0o644)
            with urlopen(url, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response:
                while chunk := response.read(DOWNLOAD_CHUNK_SIZE):
                    sink.write(chunk)
        os.replace(partial_path, dest_path)
    except Exception as error:
        try:
            if partial_path is not None:
                partial_path.unlink(missing_ok=True)
        finally:
            message = f"Failed to download {url}: {error}"
            logger.error(message)
            raise RuntimeError(message) from error

    logger.info(f"Successfully downloaded {dest_path.name}")
    return True
