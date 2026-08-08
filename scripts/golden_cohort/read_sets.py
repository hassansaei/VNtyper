"""Pure read-set evidence for golden-cohort CRAM equivalence checks."""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Iterable
from typing import TypedDict

logger = logging.getLogger(__name__)


class ReadSetEvidence(TypedDict):
    """Stable evidence that two extraction modes emitted the same BAM records."""

    count: int
    sorted_read_name_sha256: str


def _parse_count(count_output: str) -> int:
    """Parse the exact non-negative integer emitted by ``samtools view -c``."""
    rendered_count = count_output.strip()
    if not rendered_count.isdecimal():
        msg = f"samtools view -c returned an invalid non-negative count: {count_output!r}"
        logger.error(msg)
        raise ValueError(msg)
    return int(rendered_count)


def read_name_from_sam_record(record: str, line_number: int) -> str:
    """Extract one non-empty QNAME from a headerless SAM record.

    Args:
        record: One line emitted by ``samtools view``.
        line_number: One-based position used in a fail-closed diagnostic.

    Returns:
        The record's QNAME.

    Raises:
        ValueError: If the record has no non-empty first SAM field.
    """
    read_name, separator, _remaining_fields = record.partition("\t")
    if not separator or not read_name:
        msg = f"samtools view emitted a malformed SAM record at line {line_number}: {record.rstrip()!r}"
        logger.error(msg)
        raise ValueError(msg)
    return read_name


def summarize_sorted_read_names(count_output: str, sorted_read_names: Iterable[str]) -> ReadSetEvidence:
    """Hash an already-sorted QNAME stream without retaining it in memory.

    Args:
        count_output: Standard output from ``samtools view -c``.
        sorted_read_names: QNAMEs sorted under the C locale, one per line.

    Returns:
        The validated record count and incremental SHA-256 digest.

    Raises:
        ValueError: If the count is invalid, a QNAME is empty, or totals disagree.
    """
    count = _parse_count(count_output)
    digest = hashlib.sha256()
    digested_count = 0
    for line_number, line in enumerate(sorted_read_names, start=1):
        read_name = line.removesuffix("\n")
        if not read_name or "\t" in read_name or "\n" in read_name or "\r" in read_name:
            msg = f"sorted read-name stream contains an invalid QNAME at line {line_number}: {line.rstrip()!r}"
            logger.error(msg)
            raise ValueError(msg)
        digest.update(f"{read_name}\n".encode())
        digested_count += 1

    if digested_count != count:
        msg = f"samtools view -c reported {count} records but emitted {digested_count} records for the digest"
        logger.error(msg)
        raise ValueError(msg)
    return {"count": count, "sorted_read_name_sha256": digest.hexdigest()}


def summarize_unmapped_read_set(count_output: str, records_output: str) -> ReadSetEvidence:
    """Validate in-memory SAM text for focused pure tests and small callers.

    The golden-cohort runner uses :func:`summarize_sorted_read_names` over an external-sort
    stream; this convenience seam keeps record parsing independently testable.

    Args:
        count_output: Standard output from ``samtools view -c``.
        records_output: Headerless SAM records from ``samtools view``.

    Returns:
        The validated record count and sorted-QNAME SHA-256 digest.

    Raises:
        ValueError: If either output is malformed or the two commands disagree.
    """

    read_names: list[str] = []
    for line_number, record in enumerate(records_output.splitlines(), start=1):
        read_names.append(read_name_from_sam_record(record, line_number))
    return summarize_sorted_read_names(count_output, (f"{read_name}\n" for read_name in sorted(read_names)))
