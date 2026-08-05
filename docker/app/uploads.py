"""Safe construction of upload destination paths.

The filename attached to an upload is chosen by the client and is not sanitised
by the web framework. It decides where the file is written and travels onwards
into the analysis pipeline, so it is constrained here, at the boundary, rather
than relying on any later layer to cope with a name it was never designed for.

Two rules apply, and both must hold:

1. The name must match a strict allowlist -- a plain, single-segment filename
   made of ASCII letters, digits, dot, dash and underscore, ending in one of the
   extensions the caller expects.
2. The resulting path must resolve to a direct child of the job directory.

Rule 2 is unreachable past rule 1 today; it is kept as an independent check so
that widening the allowlist later cannot silently widen where files land.

Names are never repaired. A name that fails is refused, so the caller can answer
the request with an error instead of storing the upload under a name the client
did not ask for.
"""

import logging
import os
import re
from collections.abc import Sequence
from functools import cache

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

# A leading alphanumeric rules out dotfiles, `.`, `..` and leading-dash names.
# The body excludes every path separator and every character with a meaning to a
# shell, so the name stays a single, inert path segment.
_SAFE_STEM = r"[A-Za-z0-9][A-Za-z0-9._-]*"


@cache
def _name_pattern(allowed_extensions: tuple[str, ...]) -> re.Pattern[str]:
    """Compile the allowlist pattern for a set of extensions.

    Args:
        allowed_extensions: Lowercase extensions, without the leading dot.

    Returns:
        re.Pattern[str]: A pattern anchored at both ends. `\\Z` is used rather
            than `$`, because `$` also matches immediately before a trailing
            newline.
    """
    alternation = "|".join(re.escape(extension) for extension in allowed_extensions)
    return re.compile(rf"{_SAFE_STEM}\.(?:{alternation})\Z")


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
        msg = f"Uploaded filename is not an acceptable {_describe(extensions)} name"
        logger.error(msg)
        raise ValueError(msg)

    destination = os.path.join(job_input_dir, candidate)
    if os.path.dirname(os.path.realpath(destination)) != os.path.realpath(job_input_dir):
        msg = "Uploaded filename does not resolve inside the job directory"
        logger.error(msg)
        raise ValueError(msg)

    return destination
