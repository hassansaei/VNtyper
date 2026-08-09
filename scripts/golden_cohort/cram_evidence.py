"""Fail-closed A-178-2 decisions over recorded CRAM read-set evidence."""

from __future__ import annotations

import re
from typing import Any

MAX_GUARD_COUNT_DIGITS = 20
GUARD_MESSAGE_PATTERN = (
    rf"idxstats reports ([0-9]{{1,{MAX_GUARD_COUNT_DIGITS}}}) placed-unmapped reads; using stream scan"
)
PLACED_UNMAPPED_GUARD_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(GUARD_MESSAGE_PATTERN),
    re.compile(rf"ValueError: {GUARD_MESSAGE_PATTERN}"),
    re.compile(
        rf"[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}},[0-9]{{3}} - "
        rf"vntyper\.scripts\.idxstats_parsing - ERROR - {GUARD_MESSAGE_PATTERN}"
    ),
)


def parse_placed_unmapped_guard_count(stderr: str) -> int | None:
    """Return the unique placed-unmapped count from the stable guard diagnostic.

    Args:
        stderr: Complete pipeline stderr, including propagated logger output.

    Returns:
        The diagnostic count when every exact occurrence agrees, otherwise ``None``.
    """
    counts: set[int] = set()
    for line in stderr.splitlines():
        for pattern in PLACED_UNMAPPED_GUARD_PATTERNS:
            match = pattern.fullmatch(line)
            if match is not None:
                counts.add(int(match.group(1)))
                break
    if len(counts) != 1:
        return None
    return counts.pop()


def validate_cram_evidence(case: dict[str, Any], record: dict[str, Any]) -> list[str]:
    """Validate one candidate CRAM result against its declared evidence policy.

    Args:
        case: Runtime matrix case with an optional ``cram_evidence_expectation``.
        record: Result record containing raw-indexed and produced-unmapped evidence.

    Returns:
        Human-readable problems. An empty list is the only passing decision.
    """
    expectation = case.get("cram_evidence_expectation")
    if expectation is None:
        return []

    problems: list[str] = []
    expected_raw = expectation.get("raw_indexed_read_set")
    actual_raw = record.get("raw_indexed_read_set")
    if actual_raw != expected_raw:
        problems.append(f"A-178-2 raw indexed evidence differs: expected {expected_raw}, got {actual_raw}")

    scan = case.get("effective_unmapped_scan")
    actual_stream = record.get("unmapped_read_set")
    expected_stream = expectation.get("stream_read_set")
    if expectation.get("indexed_authorized") is True:
        if scan not in {"indexed", "stream"}:
            problems.append(f"A-178-2 evidence has unsupported effective scan mode {scan!r}")
            return problems
        if actual_stream != expected_stream:
            problems.append(f"A-178-2 stream evidence differs: expected {expected_stream}, got {actual_stream}")
        actual_loss = record.get("raw_indexed_loss")
        if actual_loss is None:
            problems.append("A-178-2 authorized evidence did not record the raw indexed loss")
        elif isinstance(expected_raw, dict) and isinstance(expected_stream, dict):
            raw_count = expected_raw.get("count")
            stream_count = expected_stream.get("count")
            if isinstance(raw_count, int) and isinstance(stream_count, int):
                expected_loss = stream_count - raw_count
                if actual_loss != expected_loss:
                    problems.append(f"A-178-2 raw indexed loss differs: expected {expected_loss}, got {actual_loss}")
        return problems

    if scan == "indexed":
        expected_guard_count = expectation.get("placed_unmapped_guard_count")
        actual_guard_count = record.get("placed_unmapped_guard_count")
        if type(expected_guard_count) is not int or expected_guard_count <= 0:
            problems.append(f"A-178-2 indexed guard expectation is missing or malformed: {expected_guard_count}")
        elif actual_guard_count != expected_guard_count:
            problems.append(
                f"A-178-2 indexed guard count differs: expected {expected_guard_count}, got {actual_guard_count}"
            )
        if actual_stream is not None:
            problems.append("A-178-2 indexed rejection produced an unmapped BAM instead of failing before work")
        return problems

    if scan != "stream":
        problems.append(f"A-178-2 evidence has unsupported effective scan mode {scan!r}")
        return problems

    if actual_stream != expected_stream:
        problems.append(f"A-178-2 stream evidence differs: expected {expected_stream}, got {actual_stream}")
    actual_loss = record.get("raw_indexed_loss")
    if actual_loss is None:
        problems.append("A-178-2 stream evidence did not record the raw indexed loss")
    elif isinstance(expected_raw, dict) and isinstance(expected_stream, dict):
        raw_count = expected_raw.get("count")
        stream_count = expected_stream.get("count")
        if isinstance(raw_count, int) and isinstance(stream_count, int):
            expected_loss = stream_count - raw_count
            if actual_loss != expected_loss:
                problems.append(f"A-178-2 raw indexed loss differs: expected {expected_loss}, got {actual_loss}")
    return problems
