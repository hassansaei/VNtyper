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
that widening the allowlist later cannot silently widen where files land. The
realpath comparison alone cannot see a symlinked job directory (both sides
resolve to the link's target), so `save_upload_bounded` additionally opens the
parent as a real directory and creates the file relative to that descriptor.

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

import hashlib
import logging
import os
import re
from collections.abc import Sequence
from contextlib import suppress
from dataclasses import dataclass
from functools import cache
from typing import Protocol

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

# The service and its Docker image run on Linux, where these flags are
# available. Using the attributes directly is deliberate: an unsupported host
# must fail closed at import rather than silently dropping the symlink checks.
_PARENT_DIRECTORY_FLAGS = os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW | os.O_DIRECTORY
_UPLOAD_FILE_FLAGS = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW

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


class ByteSource(Protocol):
    """A stream `save_upload_bounded` can copy from.

    Stated as the one operation the copy performs, rather than as `BinaryIO`.
    `BinaryIO` additionally demands `seek`, `tell`, `write`, `fileno`, `closed`,
    iteration and the context-manager methods -- none of which is used here, and
    all of which a perfectly serviceable stand-in may lack. Starlette's
    `UploadFile.file` satisfies this, and so does any file object, `io.BytesIO`,
    or a wrapper that counts what was asked of it.

    `size` is required rather than defaulted: this module always says how much it
    wants, and asking for no less than it uses keeps the protocol satisfied by
    sources whose `read` takes a mandatory argument.
    """

    def read(self, size: int, /) -> bytes:
        """Return up to `size` bytes from the stream.

        Args:
            size: The number of bytes requested.

        Returns:
            bytes: Up to `size` bytes, empty once the stream is exhausted.
        """
        ...  # pragma: no cover - structural declaration, never executed


@dataclass(frozen=True)
class UploadReceipt:
    """Bytes written and the identity of the created regular file."""

    bytes_written: int
    identity: tuple[int, int]
    sha256: str


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
    source: ByteSource,
    destination: str,
    max_bytes: int,
    chunk_size: int = UPLOAD_CHUNK_SIZE,
    *,
    directory_descriptor: int | None = None,
    capture_identity: bool = False,
) -> int | UploadReceipt:
    """Copy an upload to disk, writing no more than `max_bytes` of it.

    The size a client declares for its request is advisory, so the ceiling is
    applied to the bytes that are actually read. The copy stops as soon as the
    running total passes `max_bytes`, and the partly written file is removed
    before the error propagates, so a refused upload costs the volume nothing.

    The destination is opened relative to a descriptor of its parent directory,
    with `O_NOFOLLOW` on both opens and `O_EXCL` on the file itself: the parent
    must be a real directory (not a symlink swapped in after the job directory
    was created), and the destination must be a new file this call creates.
    Only the final parent component is constrained, so a mounted or symlinked
    root above the per-job directory remains supported.

    Args:
        source: The upload's byte stream, read `chunk_size` bytes at a time.
        destination: Path to write to. Must not already exist.
        max_bytes: The largest number of bytes accepted. A source of exactly
            this size is written in full; one byte more is refused.
        chunk_size: Bytes requested per read. Bounds the memory used.
        directory_descriptor: Optional already-bound destination directory.
            The copy duplicates it, so the caller retains ownership.
        capture_identity: Return the created file's device and inode with the
            byte count for an identity-bound worker handoff.

    Returns:
        The number of bytes written, or an upload receipt when identity capture
        was requested.

    Raises:
        RuntimeError: If the parent is not a usable real directory, or the
            destination cannot be created as a fresh file.
        ValueError: If the source holds more than `max_bytes` bytes.
        OSError: If wrapping, writing or closing the created file fails.
    """
    parent = os.path.dirname(destination) or "."
    name = os.path.basename(destination)
    if directory_descriptor is None:
        try:
            owned_directory_descriptor = os.open(parent, _PARENT_DIRECTORY_FLAGS)
        except OSError as error:
            msg = f"Upload refused: its job directory is not a usable real directory: {error}"
            logger.error(msg)
            raise RuntimeError(msg) from error
    else:
        owned_directory_descriptor = os.dup(directory_descriptor)

    written = 0
    digest = hashlib.sha256()
    destination_created = False
    copy_completed = False
    identity: tuple[int, int] | None = None
    try:
        try:
            try:
                descriptor = os.open(name, _UPLOAD_FILE_FLAGS, 0o644, dir_fd=owned_directory_descriptor)
                destination_created = True
                metadata = os.fstat(descriptor)
                identity = (metadata.st_dev, metadata.st_ino)
            except OSError as error:
                msg = f"Upload refused: its destination could not be created as a new file: {error}"
                logger.error(msg)
                raise RuntimeError(msg) from error

            try:
                handle = os.fdopen(descriptor, "wb")
            except BaseException:
                # `fdopen` takes ownership only when it succeeds. If wrapping
                # fails, preserve that error even if closing the raw descriptor
                # reports a second failure.
                try:
                    os.close(descriptor)
                except OSError as close_error:
                    logger.warning(f"Closing an unwrapped upload descriptor also failed: {close_error}")
                raise

            try:
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
                    digest.update(chunk)
            except BaseException:
                # A flush/close error is secondary once reading, writing or the
                # size guard has already failed. Report it without replacing
                # the exception that explains why this upload was refused.
                try:
                    handle.close()
                except OSError as close_error:
                    logger.warning(f"Closing a failed upload also failed: {close_error}")
                raise
            else:
                # With no active copy error, a flush/close failure is itself
                # the operation failure and must reach the caller.
                handle.close()

            copy_completed = True
        finally:
            if destination_created and not copy_completed:
                # Includes size, read, write, fdopen and close failures.
                # Reclaim through the descriptor that received the bytes, so a
                # path swap cannot redirect cleanup to another location.
                with suppress(OSError):
                    os.unlink(name, dir_fd=owned_directory_descriptor)
    except BaseException:
        # As with the file descriptor, a parent-close failure is secondary to
        # an active refusal or copy error and must not replace it.
        try:
            os.close(owned_directory_descriptor)
        except OSError as close_error:
            logger.warning(f"Closing the failed upload's parent directory also failed: {close_error}")
        raise
    else:
        # When the copy itself succeeded, a parent-close failure is the only
        # failure and remains visible to the caller.
        os.close(owned_directory_descriptor)

    if capture_identity:
        assert identity is not None
        return UploadReceipt(written, identity, digest.hexdigest())
    return written
