"""Fail-closed A-178-2 decisions over recorded CRAM read-set evidence."""

from __future__ import annotations

from typing import Any


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
    if scan == "indexed":
        if actual_stream is not None:
            problems.append("A-178-2 indexed rejection produced an unmapped BAM instead of failing before work")
        return problems

    if scan != "stream":
        problems.append(f"A-178-2 evidence has unsupported effective scan mode {scan!r}")
        return problems

    expected_stream = expectation.get("stream_read_set")
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
