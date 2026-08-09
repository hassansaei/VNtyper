"""Captured-command evidence for lossless indexed unmapped-read recovery."""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path

from vntyper.scripts.command_builders import (
    build_samtools_idxstats_command,
    build_samtools_unmapped_indexed_count_command,
)
from vntyper.scripts.idxstats_parsing import choose_scan, parse_idxstats

logger = logging.getLogger(__name__)

CaptureCommand = Callable[[str, str], tuple[bool, str]]


def select_unmapped_scan(
    view_path: str,
    config: dict,
    threads: int,
    output_dir: str,
    output_name: str,
    *,
    file_format: str,
    capture: CaptureCommand,
) -> str:
    """Choose a lossless scan from idxstats and exact indexed-consumer evidence.

    Args:
        view_path: Run-local indexed alignment path.
        config: Pipeline configuration. Missing keys use shipped defaults.
        threads: Thread count for samtools.
        output_dir: Run output directory for command logs.
        output_name: Base name for command logs.
        file_format: Alignment format whose scan policy is being selected.
        capture: Captured-command runner supplied by alignment preflight.

    Returns:
        ``"indexed"`` only when both evidence sources prove it complete; otherwise
        ``"stream"``.

    Raises:
        ValueError: If scan mode is invalid or forced indexed mode would lose reads.
    """
    samtools_path = config.get("tools", {}).get("samtools", "samtools")
    configured = config.get(file_format, {}).get("unmapped_scan", "auto")
    idxstats_command = build_samtools_idxstats_command(
        samtools_path=samtools_path,
        in_bam=view_path,
        threads=threads,
    )
    idxstats_log = str(Path(output_dir) / f"{output_name}_idxstats.log")
    exit_ok, output = capture(idxstats_command, idxstats_log)
    indexed_count_exit_ok = False
    indexed_count_output: str | None = None
    counts = parse_idxstats(output) if exit_ok else None
    if configured in {"auto", "indexed"} and counts is not None and counts[0] == 0:
        indexed_count_command = build_samtools_unmapped_indexed_count_command(
            samtools_path=samtools_path,
            in_bam=view_path,
            threads=threads,
        )
        indexed_count_log = str(Path(output_dir) / f"{output_name}_indexed_unmapped_count.log")
        indexed_count_exit_ok, indexed_count_output = capture(indexed_count_command, indexed_count_log)
    scan, reason = choose_scan(
        configured,
        output,
        exit_ok,
        indexed_count_text=indexed_count_output,
        indexed_count_exit_ok=indexed_count_exit_ok,
    )
    logger.info(f"Selected {scan} unmapped-read scan: {reason}")
    return scan
