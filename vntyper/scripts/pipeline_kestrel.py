"""Run and summarize the pipeline's Kestrel stage."""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast

from vntyper.scripts.run_configuration import (
    resolve_compatibility_component,
    resolve_compatibility_runtime_component,
)
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
    resolved_component: Mapping[str, object] | None = None,
    runtime_component: Mapping[str, object] | None = None,
    custom_context_active: bool = False,
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
        resolved_component: Immutable Kestrel decision component for this run.
        runtime_component: Immutable excluded Kestrel runtime component.
        custom_context_active: Whether an explicit custom profile owns this run.

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

    # Which counting path produced this result, recorded so a result can always be
    # attributed to the code that produced it. It is a **run-level** property, so it goes
    # beside the other run-level provenance rather than into a step record: `record_step`
    # derives `parsed_result` from the result file and cannot carry it, the step's
    # `command` string is compared by the golden-cohort gate as `pipeline_step_records`,
    # and trap 5 forbids inventing a new step name.
    #
    from vntyper.scripts.kestrel_counting import DEFAULT_KANALYZE_PATH

    explicit_context = resolved_component is not None or runtime_component is not None or custom_context_active
    decision = resolve_compatibility_component(
        "kestrel",
        resolved_component,
        custom_context_active=custom_context_active,
    )
    runtime = resolve_compatibility_runtime_component("kestrel", runtime_component)

    # Both conditions, not just the setting. `plan_kestrel_invocations` splits only when
    # `split_counting` is true **and** a kanalyze path is configured, so a replacement
    # config with `"kanalyze": ""` counts internally while `split_counting` is still
    # true. Reading only the setting would record "split" for a run that did not split,
    # which is worse than recording nothing: the field exists to attribute a result to
    # the code that produced it.
    settings = cast(Mapping[str, Any], runtime.get("kestrel_settings", {}))
    kanalyze_path = config.get("tools", {}).get("kanalyze", DEFAULT_KANALYZE_PATH)
    splitting = bool(settings.get("split_counting", True)) and bool(kanalyze_path)
    summary["kestrel_counting_mode"] = "split" if splitting else "internal"

    kestrel_dir = Path(dirs["kestrel"])
    tools = cast(Mapping[str, Any], config["tools"])
    reference_data = cast(Mapping[str, Any], config["reference_data"])
    start = datetime.now(timezone.utc).replace(tzinfo=None)
    runner_kwargs: dict[str, object] = {
        "vcf_path": kestrel_dir / "output.vcf",
        "output_dir": kestrel_dir,
        "fastq_files": fastq_files,
        "reference_vntr": reference_data["muc1_reference_vntr"],
        "kestrel_path": tools["kestrel"],
        "config": config,
        "sample_name": sample_name,
        "log_level": log_level,
        "cwd": cwd,
        "threads": threads,
    }
    if explicit_context:
        runner_kwargs.update(
            resolved_component=decision,
            runtime_component=runtime,
            custom_context_active=custom_context_active,
        )
    runner(**runner_kwargs)
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
