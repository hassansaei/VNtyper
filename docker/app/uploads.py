"""Constraints applied to uploads at the boundary: where they land, and how big.

The filename attached to an upload is chosen by the client and is not sanitised
by the web framework. It decides where the file is written and travels onwards
into the analysis pipeline, so it is constrained here, at the boundary, rather
than relying on any later layer to cope with a name it was never designed for.

Two rules apply, and both must hold:

1. The name must match a strict allowlist -- a plain, single-segment filename
   made of ASCII letters, digits, dot, dash, underscore and plus, ending in one
   of the extensions the caller expects. The extension is matched without regard
   to case, because sequencers and LIMS exports routinely upper-case it and
   refusing `SAMPLE.BAM` buys nothing.
2. The resulting path must resolve to a direct child of the job directory.

Rule 2 is unreachable past rule 1 today; it is kept as an independent check so
that widening the allowlist later cannot silently widen where files land.

The allowlist stays ASCII-only. A non-ASCII name raises normalisation and
encoding questions -- which of several byte sequences is "the" name, and which
of them the pipeline's tooling will reproduce -- for no clinical benefit, so it
is refused, and the refusal says so explicitly.

Names are never repaired, and never rewritten. A name that fails is refused, so
the caller can answer the request with an error instead of storing the upload
under a name the client did not ask for; a name that passes is used exactly as
it was sent, case included.

Size is constrained on the same terms. `save_upload_bounded` counts the bytes it
copies and stops the moment the running total passes the caller's ceiling, so
the limit holds on what actually arrives rather than on what the client said was
coming. A refused copy removes what it had already written.
"""

import logging
import os
import re
from collections.abc import Sequence
from contextlib import suppress
from functools import cache
from typing import BinaryIO

logger = logging.getLogger(__name__)

# Alignment inputs accepted by the job endpoint.
ALIGNMENT_EXTENSIONS: tuple[str, ...] = ("bam", "cram")

# Index files accepted alongside them. `tasks.run_vntyper_job` looks for
# `<alignment>.bai`, and the input-format docs tell users to supply `.bam.bai`
# or `.bai`, so both index extensions have to be accepted here.
INDEX_EXTENSIONS: tuple[str, ...] = ("bai", "crai")

# POSIX filesystems cap a single name at 255 bytes. Checked explicitly rather
# than folded into the pattern, where the extension length would skew the bound.
# The allowlist below is ASCII-only, so for any name that gets that far the
# character count and the byte count are the same number.
MAX_FILENAME_LENGTH = 255

# Uploads are copied a piece at a time so the memory one request costs does not
# grow with the size of the file it carries, and so the running total can be
# checked between pieces rather than only after everything has been written.
UPLOAD_CHUNK_SIZE = 1024 * 1024

# A leading alphanumeric rules out dotfiles, `.`, `..` and leading-dash names.
# The body excludes every path separator and every character with a meaning to a
# shell, so the name stays a single, inert path segment. `+` is in the set
# because it is neither of those and is how tumour/normal pairs are
# conventionally named; the characters that are absent are absent deliberately.
_SAFE_STEM = r"[A-Za-z0-9][A-Za-z0-9._+-]*"

# What the allowlist accepts, in the words a caller whose name was refused needs.
SAFE_NAME_DESCRIPTION = (
    "Use a plain filename of ASCII letters, digits, dot, dash, underscore or plus, starting with a letter or a digit"
)


@cache
def _name_pattern(allowed_extensions: tuple[str, ...]) -> re.Pattern[str]:
    """Compile the allowlist pattern for a set of extensions.

    Args:
        allowed_extensions: Lowercase extensions, without the leading dot.

    Returns:
        re.Pattern[str]: A pattern anchored at both ends, matching the extension
            without regard to case. `\\Z` is used rather than `$`, because `$`
            also matches immediately before a trailing newline. `re.ASCII`
            accompanies `re.IGNORECASE` and is load-bearing rather than
            decorative: a Unicode-aware case-insensitive match folds four
            non-ASCII letters onto ASCII ones, so `[A-Za-z]` would start
            accepting them and the allowlist would quietly stop being
            ASCII-only.
    """
    alternation = "|".join(re.escape(extension) for extension in allowed_extensions)
    return re.compile(rf"{_SAFE_STEM}\.(?:{alternation})\Z", re.IGNORECASE | re.ASCII)


def _describe(allowed_extensions: Sequence[str]) -> str:
    """Render the accepted extensions for an error message.

    Args:
        allowed_extensions: Extensions, without the leading dot.

    Returns:
        str: A human-readable list such as `.bam or .cram`.
    """
    labels = ["." + extension for extension in allowed_extensions]
    if len(labels) < 2:
        return "".join(labels)
    return f"{', '.join(labels[:-1])} or {labels[-1]}"


def safe_upload_path(
    job_input_dir: str,
    filename: str | None,
    allowed_extensions: Sequence[str] = ALIGNMENT_EXTENSIONS,
) -> str:
    """Build a destination path for an uploaded file.

    Args:
        job_input_dir: The per-job input directory.
        filename: The client-supplied filename. Untrusted.
        allowed_extensions: Extensions this upload slot accepts, without the
            leading dot. Defaults to the alignment extensions.

    Returns:
        str: A path guaranteed to sit directly inside `job_input_dir`, using the
            supplied filename unchanged.

    Raises:
        ValueError: If the filename is absent, over-long, does not match the
            allowlist, or would not resolve inside `job_input_dir`.
    """
    extensions = tuple(allowed_extensions)
    candidate = filename or ""

    if len(candidate) > MAX_FILENAME_LENGTH or not _name_pattern(extensions).match(candidate):
        msg = (
            f"Uploaded filename is not an acceptable {_describe(extensions)} name. "
            f"{SAFE_NAME_DESCRIPTION}, ending in the extension in any case, "
            f"and at most {MAX_FILENAME_LENGTH} characters long"
        )
        logger.error(msg)
        raise ValueError(msg)

    destination = os.path.join(job_input_dir, candidate)
    if os.path.dirname(os.path.realpath(destination)) != os.path.realpath(job_input_dir):
        msg = "Uploaded filename does not resolve inside the job directory"
        logger.error(msg)
        raise ValueError(msg)

    return destination


def save_upload_bounded(
    source: BinaryIO,
    destination: str,
    max_bytes: int,
    chunk_size: int = UPLOAD_CHUNK_SIZE,
) -> int:
    """Copy an upload to disk, writing no more than `max_bytes` of it.

    The size a client declares for its request is advisory, so the ceiling is
    applied to the bytes that are actually read. The copy stops as soon as the
    running total passes `max_bytes`, and the partly written file is removed
    before the error propagates, so a refused upload costs the volume nothing.

    Args:
        source: The upload's byte stream, read `chunk_size` bytes at a time.
        destination: Path to write to. Overwritten if it already exists.
        max_bytes: The largest number of bytes accepted. A source of exactly
            this size is written in full; one byte more is refused.
        chunk_size: Bytes requested per read. Bounds the memory used.

    Returns:
        int: The number of bytes written.

    Raises:
        ValueError: If the source holds more than `max_bytes` bytes.
        OSError: If the destination cannot be written.
    """
    written = 0
    try:
        with open(destination, "wb") as handle:
            while True:
                chunk = source.read(chunk_size)
                if not chunk:
                    break
                written += len(chunk)
                if written > max_bytes:
                    msg = f"Upload exceeds the maximum accepted size of {max_bytes} bytes"
                    logger.error(msg)
                    raise ValueError(msg)
                handle.write(chunk)
    except Exception:
        # Includes the refusal above: whatever reached the volume is reclaimed
        # before the caller ever sees the error.
        with suppress(OSError):
            os.remove(destination)
        raise

    return written
