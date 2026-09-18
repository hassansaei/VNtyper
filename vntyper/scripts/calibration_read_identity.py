"""Pure primary-read fingerprints for duplicate audits, never patient identity."""

from __future__ import annotations

import hashlib
import logging
import re
from collections.abc import Iterable, Mapping, Set
from dataclasses import dataclass
from typing import NoReturn

from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_BASES = frozenset("ACGTRYSWKMBDHVN")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_QUERY_OPERATIONS = frozenset({0, 1, 4, 7, 8})


@dataclass(frozen=True)
class PrimaryReadRecord:
    """Alignment identity fields, excluding headers and auxiliary tags.

    Sequence and qualities have alignment orientation. Unmapped FASTQ adapters
    supply forward orientation and flags that explicitly retain mate identity.
    None sequence/quality and hard clipping are audit defects, not an assertion
    that the source read can be reconstructed.
    """

    name: str
    sequence: str | None
    qualities: tuple[int, ...] | None
    mate: int
    flags: int
    mapping_quality: int
    contig: str | None
    position_zero_based: int
    cigar: tuple[tuple[int, int], ...]
    mate_contig: str | None
    mate_position_zero_based: int
    template_length: int


@dataclass(frozen=True)
class ReadIdentityTokens:
    """Digest-only record tokens and explicit reconstruction limitations."""

    alignment: str
    named_sequence: str
    unnamed_sequence: str
    reliable: bool
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class MultisetFingerprint:
    """A domain-separated sorted-multiset digest preserving every occurrence."""

    sha256: str
    record_count: int


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _integer(value: object, field: str, *, minimum: int = -(2**53 - 1), maximum: int = 2**53 - 1) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        _fail(f"read identity {field} must be an integer within its supported range")


def _validate_record(record: PrimaryReadRecord) -> None:
    if not isinstance(record, PrimaryReadRecord):
        _fail("read identity requires a PrimaryReadRecord")
    if not isinstance(record.name, str) or not record.name or any(char.isspace() for char in record.name):
        _fail("read identity name must be non-empty and contain no whitespace")
    _integer(record.flags, "flags", minimum=0, maximum=0xFFFF)
    _integer(record.mapping_quality, "mapping quality", minimum=0, maximum=255)
    _integer(record.mate, "mate", minimum=0, maximum=2)
    mate_flags = record.flags & (0x40 | 0x80)
    expected_mate_flags = {0: 0, 1: 0x40, 2: 0x80}[record.mate]
    if mate_flags != expected_mate_flags or bool(record.flags & 1) != bool(record.mate):
        _fail("read identity mate disagrees with alignment flags")
    for contig, position in (
        (record.contig, record.position_zero_based),
        (record.mate_contig, record.mate_position_zero_based),
    ):
        _integer(position, "position", minimum=-1)
        if contig is None:
            if position != -1:
                _fail("read identity absent contig requires absent position")
        elif not isinstance(contig, str) or not contig or position < 0:
            _fail("read identity mapped coordinate requires contig and position")
    _integer(record.template_length, "template length")
    if record.sequence is not None and (
        not isinstance(record.sequence, str) or not record.sequence or not set(record.sequence.upper()) <= _BASES
    ):
        _fail("read identity sequence must contain supported non-empty DNA symbols")
    if record.qualities is not None:
        if (
            not isinstance(record.qualities, tuple)
            or record.sequence is None
            or len(record.qualities) != len(record.sequence)
        ):
            _fail("read identity qualities must be an immutable tuple matching the sequence")
        for quality in record.qualities:
            _integer(quality, "base quality", minimum=0, maximum=93)
    if not isinstance(record.cigar, tuple):
        _fail("read identity CIGAR must be an immutable tuple")
    query_length = 0
    for operation in record.cigar:
        if not isinstance(operation, tuple) or len(operation) != 2:
            _fail("read identity CIGAR operations must be immutable pairs")
        code, length = operation
        _integer(code, "CIGAR code", minimum=0, maximum=8)
        _integer(length, "CIGAR length", minimum=1)
        if code in _QUERY_OPERATIONS:
            query_length += length
    if record.cigar and record.sequence is not None and query_length != len(record.sequence):
        _fail("read identity CIGAR query length differs from sequence length")
    if not record.flags & 4 and (not record.cigar or record.contig is None):
        _fail("read identity mapped record requires CIGAR and coordinates")


def canonical_sequence_name(name: str, mate: int) -> str:
    """Normalize only an explicit terminal mate suffix for sequence identity v1.

    The raw name remains part of alignment identity. A terminal ``/1`` or
    ``/2`` is removed only when it agrees with the explicit mate ordinal. No
    other prefix, suffix, or description transformation is permitted.

    Args:
        name: Original FASTQ header name or BAM/CRAM QNAME.
        mate: Explicit mate ordinal: zero for unpaired, one or two for paired.

    Returns:
        The name used by named sequence identity.

    Raises:
        ValueError: If the name, mate, or terminal mate suffix is invalid.
    """
    if not isinstance(name, str) or not name or any(char.isspace() for char in name):
        _fail("read identity name must be non-empty and contain no whitespace")
    _integer(mate, "mate", minimum=0, maximum=2)
    if mate == 0:
        return name
    matching_suffix = f"/{mate}"
    conflicting_suffix = "/2" if mate == 1 else "/1"
    if name.endswith(conflicting_suffix):
        _fail("read identity terminal mate suffix conflicts with the explicit mate ordinal")
    if name.endswith(matching_suffix):
        normalized = name[: -len(matching_suffix)]
        if not normalized:
            _fail("read identity terminal mate suffix requires a non-empty base name")
        return normalized
    return name


def read_identity_tokens(record: PrimaryReadRecord) -> ReadIdentityTokens | None:
    """Produce canonical record tokens for an external-sort duplicate audit.

    Args:
        record: Decoded alignment-oriented record fields; no aux tags or headers.

    Returns:
        Alignment, named-read and unnamed-read digests, or None for secondary
        and supplementary alignments. Logical matches are evidence duplicates
        for adjudication; they never establish that two people are identical.

    Raises:
        ValueError: If record fields are inconsistent or unsupported.
    """
    _validate_record(record)
    if record.flags & (0x100 | 0x800):
        return None
    reasons: list[str] = []
    sequence = record.sequence.upper() if record.sequence is not None else None
    qualities = record.qualities
    if sequence is None:
        reasons.append("missing_sequence")
    if qualities is None:
        reasons.append("missing_qualities")
    if any(code == 5 for code, _length in record.cigar):
        reasons.append("hard_clipped_sequence")
    if record.flags & 0x10:
        sequence = sequence.translate(_COMPLEMENT)[::-1] if sequence is not None else None
        qualities = tuple(reversed(qualities)) if qualities is not None else None
    content = {
        "sequence": sequence,
        "qualities": list(qualities) if qualities is not None else None,
        "mate": record.mate,
    }
    named_content = {**content, "name": canonical_sequence_name(record.name, record.mate)}
    alignment = {
        **content,
        "name": record.name,
        "flags": record.flags,
        "mapping_quality": record.mapping_quality,
        "contig": record.contig,
        "position_zero_based": record.position_zero_based,
        "cigar": [list(item) for item in record.cigar],
        "mate_contig": record.mate_contig,
        "mate_position_zero_based": record.mate_position_zero_based,
        "template_length": record.template_length,
    }
    return ReadIdentityTokens(
        canonical_sha256({"schema_version": "alignment-identity-v1", "record": alignment}),
        canonical_sha256({"schema_version": "named-read-identity-v1", "record": named_content}),
        canonical_sha256({"schema_version": "unnamed-read-identity-v1", "record": content}),
        not reasons,
        tuple(reasons),
    )


def digest_sorted_tokens(tokens: Iterable[str]) -> MultisetFingerprint:
    """Hash sorted fixed-width record tokens in bounded memory, keeping repeats.

    Args:
        tokens: Non-empty ascending stream from a private external sort. Equal
            adjacent values remain separate records; never pass a set here.

    Returns:
        The domain-separated digest and occurrence count. Sorting is a separate
        I/O responsibility; this function verifies order without buffering reads.

    Raises:
        ValueError: If tokens are empty, malformed or not sorted.
    """
    if not isinstance(tokens, Iterable) or isinstance(tokens, (str, bytes, Mapping, Set)):
        _fail("read multiset requires an ordered stream preserving duplicate tokens")
    digest = hashlib.sha256(b"vntyper-read-multiset-v1\n")
    previous = ""
    count = 0
    for token in tokens:
        if not isinstance(token, str) or _SHA256.fullmatch(token) is None:
            _fail("read multiset token must be a lowercase SHA256")
        if token < previous:
            _fail("read multiset tokens must be sorted")
        digest.update(token.encode("ascii") + b"\n")
        previous = token
        count += 1
    if count == 0:
        _fail("read multiset must contain at least one primary record")
    return MultisetFingerprint(digest.hexdigest(), count)
