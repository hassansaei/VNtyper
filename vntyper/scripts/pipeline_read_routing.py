"""Count and route conversion FASTQs before downstream pipeline stages."""

from __future__ import annotations

import gzip
import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any, TextIO

from vntyper.scripts.read_layout import classify_layout, route_fastqs

logger = logging.getLogger(__name__)

_FASTQ_KEYS = ("r1", "r2", "other", "single")


def _open_fastq(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r", encoding="utf-8")


def count_fastq_records(path: str | Path, *, lines_per_record: int) -> int:
    """Count complete records in one plain or gzipped FASTQ.

    A zero-byte file returns immediately. Every non-zero file is read so its
    record count can be reported and paired outputs can be checked for parity.

    Args:
        path: FASTQ path produced by alignment conversion.
        lines_per_record: Configured number of lines in one FASTQ record.

    Returns:
        The number of complete FASTQ records.

    Raises:
        ValueError: If ``lines_per_record`` is invalid, the file cannot be
            inspected or it contains an incomplete record.
    """
    if isinstance(lines_per_record, bool) or not isinstance(lines_per_record, int) or lines_per_record <= 0:
        msg = "lines_per_record must be a positive integer."
        logger.error(msg)
        raise ValueError(msg)

    fastq_path = Path(path)
    try:
        if fastq_path.stat().st_size == 0:
            return 0
        with _open_fastq(fastq_path) as handle:
            line_count = sum(1 for _ in handle)
    except (OSError, EOFError, UnicodeError) as exc:
        msg = f"Could not count FASTQ records in {fastq_path}: {exc}"
        logger.error(msg)
        raise ValueError(msg) from exc

    if line_count % lines_per_record != 0:
        msg = (
            f"FASTQ {fastq_path} has {line_count} lines, which is not divisible by "
            f"the configured {lines_per_record} lines per record."
        )
        logger.error(msg)
        raise ValueError(msg)
    return line_count // lines_per_record


def _configured_lines_per_record(config: Mapping[str, Any]) -> int:
    utils_config = config.get("utils", {})
    if not isinstance(utils_config, Mapping):
        msg = "Configuration key utils.fastq_validation_lines must be a positive integer."
        logger.error(msg)
        raise ValueError(msg)
    value = utils_config.get("fastq_validation_lines", 4)
    if isinstance(value, bool) or not isinstance(value, int) or value != 4:
        msg = "Configuration key utils.fastq_validation_lines must be 4 for FASTQ records."
        logger.error(msg)
        raise ValueError(msg)
    return value


def route_converted_fastqs(
    produced: tuple[str, str, str, str],
    config: Mapping[str, Any],
) -> tuple[str, str | None]:
    """Count and route all four FASTQs emitted by alignment conversion.

    Args:
        produced: Paths in R1, R2, other and singleton order.
        config: Pipeline configuration. Replacement configurations may omit
            ``utils.fastq_validation_lines`` and retain the shipped default of
            four lines per record.

    Returns:
        The first routed FASTQ and the optional paired FASTQ.

    Raises:
        ValueError: If counting fails, the layout is mixed or empty, or any
            produced FASTQ would be stranded.
    """
    paths = dict(zip(_FASTQ_KEYS, produced, strict=True))
    lines_per_record = _configured_lines_per_record(config)
    counts = {key: count_fastq_records(path, lines_per_record=lines_per_record) for key, path in paths.items()}
    layout = classify_layout(**counts)
    consumed, stranded = route_fastqs(layout, paths, counts)
    details = ", ".join(f"{paths[key]}: {counts[key]} records" for key in _FASTQ_KEYS)

    if layout in {"mixed", "empty"} or stranded or len(consumed) not in {1, 2}:
        msg = f"FASTQ layout '{layout}' cannot be consumed without dropping reads. Produced FASTQs: {details}"
        logger.error(msg)
        raise ValueError(msg)

    logger.info(f"Detected {layout} FASTQ layout. Produced FASTQs: {details}")
    return consumed[0], consumed[1] if len(consumed) == 2 else None
