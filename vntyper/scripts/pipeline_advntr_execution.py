"""Focused native genotype invocation from one resolved pipeline configuration."""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

from vntyper.modules.advntr.advntr_decision_config import project_advntr_settings
from vntyper.modules.advntr.advntr_genotyping import resolve_advntr_threads
from vntyper.scripts.pipeline_advntr_run_context import AdvntrRunContext
from vntyper.scripts.pipeline_caller_native import caller_native_policy_argv, prepare_caller_native_execution
from vntyper.scripts.pipeline_research_advntr import research_advntr_policy_argv
from vntyper.scripts.run_configuration import RunConfiguration, cast_mapping

logger = logging.getLogger(__name__)


def execute_advntr_genotype(
    *,
    configuration: RunConfiguration,
    native_context: AdvntrRunContext,
    config: Mapping[str, Any],
    alignment: Path,
    output: str | Path,
    cwd: str,
    threads: int,
    additional_commands: str | None,
    background: Path | None,
    invoke: Callable[..., int],
) -> int:
    """Execute the original stage API with explicit verified calibrated arguments.

    An approved bundle supplies its preflight-verified arguments; a research decision
    profile supplies the arguments rendered from its own calibrated legacy values.

    Args:
        configuration: One resolved pipeline decision configuration.
        native_context: Run-owned model snapshot and observed executable prefix.
        config: Pipeline nondecision configuration.
        alignment: Prepared native alignment.
        output: Native stage output directory.
        cwd: Pipeline working directory.
        threads: Pipeline thread count.
        additional_commands: Optional legacy native extension override.
        background: Exact approved background snapshot when calibrated.
        invoke: Native stage process boundary, injectable for orchestration tests.

    Returns:
        Native stage status, unchanged for the caller to handle.

    Raises:
        ValueError: If calibrated runtime assets or background changed before execution,
            or a research profile selects exact adVNTR mode.
    """
    runtime = configuration.advntr_runtime
    if additional_commands is not None:
        runtime = {
            **runtime,
            "settings": {**cast_mapping(runtime["settings"]), "additional_commands": additional_commands},
        }
    kwargs: dict[str, Any] = {}
    calibration = configuration.caller_calibration
    if calibration is not None and calibration.bundle.advntr_policy is not None:
        observed = prepare_caller_native_execution(calibration, native_context)
        kwargs["calibrated_policy_arguments"] = caller_native_policy_argv(calibration, observed, background)
    elif calibration is None and "calibrated_calling" in configuration.advntr:
        # A research profile carries calibrated adVNTR values but no bundle, so its
        # explicit policy is rendered from the profile with the stage's own -t value.
        settings = project_advntr_settings(configuration.advntr, runtime).command_mapping()
        research = research_advntr_policy_argv(configuration, resolve_advntr_threads(settings, threads))
        if research is not None:
            kwargs["calibrated_policy_arguments"] = research
    return invoke(
        native_context.model_snapshot,
        alignment,
        output,
        "output",
        config={**config, "tools": dict(native_context.tools)},
        cwd=cwd,
        pipeline_threads=threads,
        resolved_component=configuration.advntr,
        runtime_component=runtime,
        custom_context_active=configuration.decision_profile.source == "explicit-cli",
        advntr_version=native_context.version,
        **kwargs,
    )
