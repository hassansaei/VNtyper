"""Pure local-versus-remote policy for CRAM reference URI values."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

logger = logging.getLogger(__name__)

_HEADER_URI_SCHEME = re.compile(r"^(?P<scheme>[A-Za-z][A-Za-z0-9+.-]*):", re.IGNORECASE)
_REF_PATH_REMOTE_URI = re.compile(r"(?:^|:)(?P<scheme>[A-Za-z][A-Za-z0-9+.-]*)://", re.IGNORECASE)


@dataclass(frozen=True)
class RemoteHeaderReference:
    """Path-free identity of one remote reference declared by an SQ line."""

    contig: str
    scheme: str


def header_reference_scheme(value: str) -> str | None:
    """Return a non-file scheme anchored at the start of one SAM ``UR`` value.

    Args:
        value: One complete ``UR`` tag value.

    Returns:
        The normalized remote scheme, or ``None`` for local, relative, and file
        URI values.
    """
    match = _HEADER_URI_SCHEME.match(value)
    if match is None:
        return None
    scheme = match.group("scheme").lower()
    return None if scheme == "file" else scheme


def ref_path_remote_scheme(value: str) -> str | None:
    """Return a remote URL scheme embedded in a colon-separated ``REF_PATH``.

    Args:
        value: A complete htslib ``REF_PATH`` search list.

    Returns:
        The normalized remote scheme, or ``None`` for local entries and file
        URLs.
    """
    for match in _REF_PATH_REMOTE_URI.finditer(value):
        scheme = match.group("scheme").lower()
        if scheme != "file":
            return scheme
    return None


def allow_ambient_reference_resolution(config: dict) -> bool:
    """Read the ambient-reference waiver without applying truthiness coercion.

    Args:
        config: Pipeline configuration; a missing setting uses the shipped
            default ``False``.

    Returns:
        The configured boolean.

    Raises:
        ValueError: If the setting is not an actual boolean.
    """
    configured = config.get("cram", {}).get("allow_ambient_reference_resolution", False)
    if not isinstance(configured, bool):
        message = "cram.allow_ambient_reference_resolution must be true or false"
        logger.error(message)
        raise ValueError(message)
    return configured


def first_remote_header_reference(header: str) -> RemoteHeaderReference | None:
    """Find the first SQ line whose UR tag declares a remote reference.

    Args:
        header: SAM header text already read from a CRAM input.

    Returns:
        The contig and remote scheme without the URI path, or ``None`` when all
        SQ reference values are local.
    """
    for line in header.splitlines():
        fields = line.split("\t")
        if not fields or fields[0] != "@SQ":
            continue
        contig = "unknown"
        for field in fields[1:]:
            key, separator, value = field.partition(":")
            if separator and key == "SN":
                contig = value
                break
        for field in fields[1:]:
            key, separator, value = field.partition(":")
            if separator and key == "UR":
                scheme = header_reference_scheme(value)
                if scheme is not None:
                    return RemoteHeaderReference(contig=contig, scheme=scheme)
    return None


def enforce_header_reference_policy(header: str, *, allow_ambient: bool) -> None:
    """Reject a remote SQ URI unless the operator explicitly allows ambient access.

    Args:
        header: SAM header text already read from a CRAM input.
        allow_ambient: Whether the existing ambient-resolution waiver is enabled.

    Raises:
        ValueError: If default local-only mode encounters a remote SQ URI.
    """
    if allow_ambient:
        return
    remote = first_remote_header_reference(header)
    if remote is None:
        return
    message = (
        "Remote CRAM header reference is disabled by policy: "
        f"contig={remote.contig}, scheme={remote.scheme}. Replace the @SQ UR with a local path, relative path, or "
        "file-scheme reference, or set cram.allow_ambient_reference_resolution=true to accept network access."
    )
    logger.error(message)
    raise ValueError(message)
