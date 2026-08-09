"""Pure local-versus-remote policy for CRAM reference URI values."""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

logger = logging.getLogger(__name__)

_HEADER_URI_SCHEME = re.compile(r"^(?P<scheme>[A-Za-z][A-Za-z0-9+.-]*):", re.IGNORECASE)
_REF_PATH_REMOTE_URI = re.compile(r"(?:^|:)(?P<scheme>[A-Za-z][A-Za-z0-9+.-]*)://", re.IGNORECASE)
_INVALID_LOCAL_HEADER_REFERENCE = "Invalid local CRAM header reference URI"


@dataclass(frozen=True)
class RemoteHeaderReference:
    """Path-free identity of one remote reference declared by an SQ line."""

    contig: str
    scheme: str


@dataclass(frozen=True)
class LocalHeaderReference:
    """One local SQ reference path associated with its contig and digest."""

    contig: str
    m5: str | None
    path: str


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


def remote_ref_path_suffix(value: str) -> str | None:
    """Return a remote REF_PATH suffix without any local search entry."""
    match = _REF_PATH_REMOTE_URI.search(value)
    if match is None:
        return None
    remote = value[match.start("scheme") :]
    scheme_separator = remote.find("://")
    if re.search(r":(?=/)", remote[scheme_separator + 3 :]):
        raise ValueError("Ambient REF_PATH cannot retain a local entry after a remote URL")
    return remote


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
                    return RemoteHeaderReference(contig=quote(contig, safe="._:-"), scheme=scheme)
    return None


def local_header_references(header: str, input_alignment: str | Path) -> tuple[LocalHeaderReference, ...]:
    """Resolve local SQ ``UR`` values with their contig and digest identity.

    Args:
        header: SAM header text already read from the owned CRAM boundary.
        input_alignment: Original operator CRAM path used as the base for relative UR values.

    Returns:
        Local reference records in first-occurrence order. Remote schemes are
        omitted for the separate ambient-reference policy.

    Raises:
        ValueError: If a file URI is ambiguous or does not name a local path.
    """
    base_directory = Path(os.path.abspath(input_alignment)).parent
    result: list[LocalHeaderReference] = []
    seen: set[tuple[str, str | None, str]] = set()
    for line in header.splitlines():
        fields = line.split("\t")
        if not fields or fields[0] != "@SQ":
            continue
        tags = {key: value for field in fields[1:] for key, separator, value in (field.partition(":"),) if separator}
        contig = tags.get("SN", "unknown")
        m5 = tags.get("M5")
        for field in fields[1:]:
            key, separator, value = field.partition(":")
            if not separator or key != "UR" or header_reference_scheme(value) is not None:
                continue
            scheme_match = _HEADER_URI_SCHEME.match(value)
            if scheme_match is not None:
                parsed = urlsplit(value)
                if (
                    parsed.scheme.lower() != "file"
                    or parsed.netloc.lower() not in {"", "localhost"}
                    or parsed.query
                    or parsed.fragment
                ):
                    raise ValueError(_INVALID_LOCAL_HEADER_REFERENCE)
                try:
                    candidate_value = unquote(parsed.path, errors="strict")
                except UnicodeDecodeError as error:
                    raise ValueError(_INVALID_LOCAL_HEADER_REFERENCE) from error
            else:
                candidate_value = value
            if not candidate_value or "\x00" in candidate_value:
                raise ValueError(_INVALID_LOCAL_HEADER_REFERENCE)
            candidate = Path(candidate_value)
            absolute = str(candidate if candidate.is_absolute() else base_directory / candidate)
            normalized = os.path.abspath(os.path.normpath(absolute))
            identity = (contig, m5, normalized)
            if identity not in seen:
                seen.add(identity)
                result.append(LocalHeaderReference(contig=contig, m5=m5, path=normalized))
    return tuple(result)


def local_header_reference_paths(header: str, input_alignment: str | Path) -> tuple[str, ...]:
    """Return deduplicated local SQ reference paths in header order.

    Args:
        header: SAM header text already read from a CRAM input.
        input_alignment: Original operator CRAM path used as the base for relative UR values.

    Returns:
        Absolute local filesystem candidates in first-occurrence order.
    """
    result: list[str] = []
    seen: set[str] = set()
    for reference in local_header_references(header, input_alignment):
        if reference.path not in seen:
            seen.add(reference.path)
            result.append(reference.path)
    return tuple(result)


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
