"""Bounded external-sort fingerprints of logical reads, with private temporaries."""

from __future__ import annotations

import heapq
import logging
import os
import tempfile
from collections.abc import Iterable
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.calibration_read_identity import (
    PrimaryReadRecord,
    digest_sorted_tokens,
    read_identity_tokens,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LogicalReadFingerprint:
    """Order-independent evidence identities, counts, and reconstruction defects."""

    alignment_sha256: str
    named_sequence_sha256: str
    unnamed_sequence_sha256: str
    primary_record_count: int
    sequence_identity_reliable: bool
    reasons: tuple[str, ...]


def _resource_bound(value: int, name: str, minimum: int, maximum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        message = f"read fingerprint {name} must be an integer in [{minimum}, {maximum}]"
        logger.error(message)
        raise ValueError(message)


def _write_chunk(values: list[str], directory: Path) -> Path:
    descriptor, filename = tempfile.mkstemp(prefix="tokens-", dir=directory)
    with os.fdopen(descriptor, "w", encoding="ascii", newline="\n") as handle:
        handle.writelines(value + "\n" for value in sorted(values))
    return Path(filename)


def _merge_chunks(paths: list[Path], directory: Path) -> Path:
    descriptor, filename = tempfile.mkstemp(prefix="merged-", dir=directory)
    with os.fdopen(descriptor, "w", encoding="ascii", newline="\n") as output, ExitStack() as stack:
        handles = [stack.enter_context(path.open(encoding="ascii")) for path in paths]
        output.writelines(heapq.merge(*handles))
    for path in paths:
        path.unlink()
    return Path(filename)


def _reduce_chunks(paths: list[Path], directory: Path, fan_in: int) -> Path:
    while len(paths) > 1:
        reduced: list[Path] = []
        for start in range(0, len(paths), fan_in):
            group = paths[start : start + fan_in]
            reduced.append(group[0] if len(group) == 1 else _merge_chunks(group, directory))
        paths = reduced
    return paths[0]


def fingerprint_read_records(
    records: Iterable[PrimaryReadRecord],
    *,
    temporary_parent: Path | None = None,
    chunk_records: int = 100_000,
    merge_fan_in: int = 32,
) -> LogicalReadFingerprint:
    """Fingerprint a decoded record stream without retaining all reads in memory.

    Args:
        records: Alignment records owned/closed by the calling input adapter.
        temporary_parent: Parent for a newly created private temporary directory.
        chunk_records: Maximum primary records held per sorting chunk, at most a million.
        merge_fan_in: Maximum chunk files opened by one merge, between two and 64.

    Returns:
        Three multiset digests and explicit evidence limitations. Temporary files
        contain only record digests, have mode 0600 within a mode-0700 directory,
        and are removed on success or failure. Biological identity and truth
        adjudication remain separate from these logical-read comparisons.

    Raises:
        ValueError: If resource bounds, records, or primary evidence are invalid.
        OSError: If private sorting or input iteration encounters filesystem errors.
    """
    _resource_bound(chunk_records, "chunk size", 1, 1_000_000)
    _resource_bound(merge_fan_in, "merge fan-in", 2, 64)
    count = 0
    reasons: set[str] = set()
    buffers: list[list[str]] = [[], [], []]
    chunks: list[list[Path]] = [[], [], []]
    with tempfile.TemporaryDirectory(prefix="vntyper-read-audit-", dir=temporary_parent) as temporary:
        directory = Path(temporary)
        for record in records:
            tokens = read_identity_tokens(record)
            if tokens is None:
                continue
            count += 1
            reasons.update(tokens.reasons)
            for buffer, token in zip(
                buffers, (tokens.alignment, tokens.named_sequence, tokens.unnamed_sequence), strict=True
            ):
                buffer.append(token)
            if len(buffers[0]) == chunk_records:
                for buffer, stream_chunks in zip(buffers, chunks, strict=True):
                    stream_chunks.append(_write_chunk(buffer, directory))
                    buffer.clear()
        if count == 0:
            message = "read fingerprint requires at least one primary record"
            logger.error(message)
            raise ValueError(message)
        for buffer, stream_chunks in zip(buffers, chunks, strict=True):
            if buffer:
                stream_chunks.append(_write_chunk(buffer, directory))
                buffer.clear()
        digests: list[str] = []
        for stream_chunks in chunks:
            merged = _reduce_chunks(stream_chunks, directory, merge_fan_in)
            with merged.open(encoding="ascii") as handle:
                fingerprint = digest_sorted_tokens(line.removesuffix("\n") for line in handle)
            if fingerprint.record_count != count:
                message = "read fingerprint sorting changed the primary record count"
                logger.error(message)
                raise ValueError(message)
            digests.append(fingerprint.sha256)
        return LogicalReadFingerprint(digests[0], digests[1], digests[2], count, not reasons, tuple(sorted(reasons)))
