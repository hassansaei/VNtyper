"""Atomic output installation for the mutation harness."""

import os
import tempfile
from pathlib import Path


def atomic_write_text(destination: Path, content: str) -> None:
    """Install complete UTF-8 text at ``destination`` atomically.

    Args:
        destination: Final path to replace.
        content: Complete report content to install.
    """
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        dir=destination.parent,
    )
    temporary = Path(temporary_name)
    try:
        stream = os.fdopen(descriptor, "w", encoding="utf-8")
        descriptor = -1
        with stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        temporary.unlink(missing_ok=True)
