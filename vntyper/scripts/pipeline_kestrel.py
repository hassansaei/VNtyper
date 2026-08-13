"""Run and summarize the pipeline's Kestrel stage."""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast

from vntyper.scripts.summary import record_step
from vntyper.scripts.summary_steps import STEP_KESTREL

logger = logging.getLogger(__name__)


def run_kestrel_stage(
    *,
    fastq_files: tuple[str, ...],
    dirs: Mapping[str, str],
    config: Mapping[str, Any],
    sample_name: str,
    log_level: int,
    cwd: str,
    summary: dict[str, Any],
    summary_file_path: str,
    runner: Callable[..., None],
    threads: int = 4,
) -> None:
    """Forward one complete routed read set and record success atomically.

    Args:
        fastq_files: Non-empty immutable FASTQ paths in routing order.
        dirs: Pipeline stage directories.
        config: Pipeline configuration.
        sample_name: Sample name embedded by Kestrel.
        log_level: Numeric Python logging level.
        cwd: Working directory captured at pipeline entry.
        summary: Pipeline summary payload.
        summary_file_path: Incremental JSON summary destination.
        runner: Patchable Kestrel runner owned by the calling pipeline module.
        threads: The run's total thread budget, forwarded to the runner for the
            KAnalyze counting step. Defaults to 4, matching ``config.json``'s
            ``default_values.threads``, so a caller that has not been updated still
            binds.

    Raises:
        ValueError: If the routed tuple is empty, malformed, or contains duplicates.
        Exception: Any runner or summary failure, without retaining a success step.
    """
    if (
        not isinstance(fastq_files, tuple)
        or not fastq_files
        or any(not isinstance(path, str) or not path for path in fastq_files)
    ):
        raise ValueError("Kestrel requires a non-empty tuple of FASTQ paths.")
    resolved = tuple(Path(path).resolve(strict=False) for path in fastq_files)
    if len(set(resolved)) != len(resolved):
        raise ValueError("Kestrel FASTQ inputs contain duplicate paths.")

    kestrel_dir = Path(dirs["kestrel"])
    tools = cast(Mapping[str, Any], config["tools"])
    reference_data = cast(Mapping[str, Any], config["reference_data"])
    start = datetime.now(timezone.utc).replace(tzinfo=None)
    runner(
        vcf_path=kestrel_dir / "output.vcf",
        output_dir=kestrel_dir,
        fastq_files=fastq_files,
        reference_vntr=reference_data["muc1_reference_vntr"],
        kestrel_path=tools["kestrel"],
        config=config,
        sample_name=sample_name,
        log_level=log_level,
        cwd=cwd,
        threads=threads,
    )
    end = datetime.now(timezone.utc).replace(tzinfo=None)

    steps = summary["steps"]
    previous_length = len(steps)
    try:
        record_step(
            summary,
            STEP_KESTREL,
            str(kestrel_dir / "kestrel_result.tsv"),
            "tsv",
            "run_kestrel(...)",
            start,
            end,
            write_summary_path=summary_file_path,
        )
    except Exception:
        del steps[previous_length:]
        raise
