"""Focused planning and execution for converting an alignment view to FASTQ."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from vntyper.scripts.command_builders import build_bam_to_fastq_command

logger = logging.getLogger(__name__)


class CommandRunner(Protocol):
    """Callable contract used to execute a shell command."""

    def __call__(
        self,
        command: str,
        log_file: str,
        critical: bool = False,
        cwd: str | Path | None = None,
    ) -> bool:
        """Execute a command and report whether it succeeded."""
        ...


@dataclass(frozen=True)
class AlignmentFastqConversionPaths:
    """All paths owned by one alignment-to-FASTQ conversion."""

    fastq_r1: Path
    fastq_r2: Path
    fastq_other: Path
    fastq_single: Path
    sort_tmp_prefix: Path
    log_file: Path


def plan_alignment_fastq_conversion(
    *,
    output: str | Path,
    output_name: str,
) -> AlignmentFastqConversionPaths:
    """Derive conversion outputs and the name-sort spill prefix.

    Args:
        output: Run output directory.
        output_name: Base name for conversion artifacts.

    Returns:
        All conversion paths, including a sort prefix inside ``output``.
    """
    output_path = Path(output)
    return AlignmentFastqConversionPaths(
        fastq_r1=output_path / f"{output_name}_R1.fastq.gz",
        fastq_r2=output_path / f"{output_name}_R2.fastq.gz",
        fastq_other=output_path / f"{output_name}_other.fastq.gz",
        fastq_single=output_path / f"{output_name}_single.fastq.gz",
        sort_tmp_prefix=output_path / f"{output_name}_sort_tmp",
        log_file=output_path / f"{output_name}_sort_fastq.log",
    )


def run_alignment_fastq_conversion(
    *,
    paths: AlignmentFastqConversionPaths,
    final_bam: str | Path,
    samtools_path: str,
    threads: int,
    command_runner: CommandRunner,
) -> tuple[str, str, str, str]:
    """Name-sort an alignment and write the four FASTQ routing outputs.

    Args:
        paths: Planned conversion output paths.
        final_bam: Alignment view to convert.
        samtools_path: samtools invocation from config.
        threads: Thread count for sort and FASTQ conversion.
        command_runner: Injected command executor from the calling stage.

    Returns:
        Paths to R1, R2, other, and singleton FASTQs, in that order.

    Raises:
        RuntimeError: If the conversion command reports failure.
    """
    command = build_bam_to_fastq_command(
        samtools_path=samtools_path,
        in_bam=final_bam,
        threads=threads,
        fastq_r1=paths.fastq_r1,
        fastq_r2=paths.fastq_r2,
        fastq_other=paths.fastq_other,
        fastq_single=paths.fastq_single,
        sort_tmp_prefix=paths.sort_tmp_prefix,
    )
    logger.info(f"Executing BAM to FASTQ conversion with command: {command}")

    if not command_runner(str(command), str(paths.log_file), critical=True):
        message = "BAM to FASTQ conversion failed."
        logger.error(message)
        raise RuntimeError(message)
    logger.info("BAM to FASTQ conversion completed.")

    return (
        str(paths.fastq_r1),
        str(paths.fastq_r2),
        str(paths.fastq_other),
        str(paths.fastq_single),
    )
