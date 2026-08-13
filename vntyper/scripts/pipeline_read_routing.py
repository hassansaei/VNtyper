"""Count and route conversion FASTQs before downstream pipeline stages."""

from __future__ import annotations

import codecs
import gzip
import json
import logging
import zlib
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import IO, Any

from vntyper.scripts.read_layout import classify_layout, route_fastqs

logger = logging.getLogger(__name__)

_FASTQ_KEYS = ("r1", "r2", "other", "single")

#: Read size for :func:`_count_lines`. Module-level so a test can shrink it and place a
#: multi-byte character across a chunk boundary; nothing else should change it.
_COUNT_CHUNK_BYTES = 1 << 20


def _open_fastq(path: Path) -> gzip.GzipFile | IO[bytes]:
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def _count_lines(path: Path) -> int:
    """Count lines without decoding the file into text, validating its encoding.

    Classifying the read layout needs the number of newlines and nothing else, so
    materialising 3.85M ``str`` objects to reach it was pure cost -- measured 1044 ms
    against 723 ms for this. ``zcat | wc -l`` is 710 ms on the same input, so gzip
    decompression is the floor and there is nothing further to win here.

    The text reader validated UTF-8 as a side effect, and this restores that
    deliberately rather than dropping it: a non-UTF-8 FASTQ still fails. The decoder is
    **incremental** because strict-decoding each chunk on its own rejects valid files --
    a multi-byte character straddling a read boundary raises "unexpected end of data" --
    while ``final=True`` still catches a genuinely truncated final character.

    Args:
        path: FASTQ file, gzipped or plain.

    Returns:
        int: The number of lines, counting a final line with no terminator.

    Raises:
        UnicodeDecodeError: If the bytes are not valid UTF-8.
        OSError: If the file cannot be read.
        zlib.error: If a gzipped file cannot be decompressed.
    """
    decoder = codecs.getincrementaldecoder("utf-8")("strict")
    lines = 0
    last_byte = b"\n"
    with _open_fastq(path) as handle:
        while True:
            chunk = handle.read(_COUNT_CHUNK_BYTES)
            if not chunk:
                break
            lines += chunk.count(b"\n")
            last_byte = chunk[-1:]
            decoder.decode(chunk)
    decoder.decode(b"", final=True)
    if last_byte != b"\n":
        lines += 1
    return lines


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
        line_count = _count_lines(fastq_path)
    except (OSError, EOFError, UnicodeError, zlib.error) as exc:
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
) -> tuple[str, ...]:
    """Count and route all four FASTQs emitted by alignment conversion.

    Args:
        produced: Paths in R1, R2, other and singleton order.
        config: Pipeline configuration. Replacement configurations may omit
            ``utils.fastq_validation_lines`` and retain the shipped default of
            four lines per record.

    Returns:
        Every non-empty routed FASTQ in R1/R2/other/single order.

    Raises:
        ValueError: If counting fails, the layout is invalid or empty, or any
            produced FASTQ would be stranded.
    """
    paths = dict(zip(_FASTQ_KEYS, produced, strict=True))
    lines_per_record = _configured_lines_per_record(config)
    # zlib releases the GIL while decompressing, so the four counts genuinely overlap:
    # measured 723 ms serial against 463 ms here. The futures are resolved in
    # `_FASTQ_KEYS` order rather than as they complete, which keeps the *first* corrupt
    # file reported the same one the serial comprehension reported -- otherwise which
    # path the ValueError names becomes a scheduling race.
    with ThreadPoolExecutor(max_workers=len(paths)) as pool:
        futures = {
            key: pool.submit(count_fastq_records, path, lines_per_record=lines_per_record)
            for key, path in paths.items()
        }
        counts = {key: future.result() for key, future in futures.items()}
    layout = classify_layout(**counts)
    consumed, stranded = route_fastqs(layout, paths, counts)
    details = ", ".join(f"{paths[key]}: {counts[key]} records" for key in _FASTQ_KEYS)

    if layout == "invalid":
        msg = f"FASTQ layout 'invalid': mate outputs are inconsistent. Produced FASTQs: {details}"
        logger.error(msg)
        raise ValueError(msg)

    if layout == "empty" or stranded or not 1 <= len(consumed) <= 4:
        msg = f"FASTQ layout '{layout}' cannot be consumed without dropping reads. Produced FASTQs: {details}"
        logger.error(msg)
        raise ValueError(msg)

    record = {
        "counts": counts,
        "layout": layout,
        "selected": [Path(path).name for path in consumed],
    }
    logger.info(f"READ_SET_ROUTING {json.dumps(record, sort_keys=True, separators=(',', ':'))}")
    return consumed
