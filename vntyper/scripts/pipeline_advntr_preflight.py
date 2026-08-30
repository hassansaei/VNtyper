"""Pure adVNTR enablement and model-reference planning for the pipeline."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, NoReturn

from vntyper.scripts.reference_resolution import resolve_from_mapping

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AdvntrPreflight:
    """Resolved adVNTR state needed before production I/O begins."""

    enabled: bool
    reference: str | None


def select_advntr_reference(config: dict[str, Any], reference_assembly: str) -> str | None:
    """Return the adVNTR database configured for an assembly's coordinate system.

    Args:
        config: Complete pipeline configuration.
        reference_assembly: Supported assembly label.

    Returns:
        Database path, or ``None`` when its configuration key has no path.
    """
    resolved = resolve_from_mapping("advntr", reference_assembly, config.get("reference_data", {}))
    return resolved.value if resolved is not None else None


def _refuse_invalid_preflight(message: str, *, log_error: bool) -> NoReturn:
    if log_error:
        logger.error(message)
    raise ValueError(message)


def _plan_advntr_preflight(
    config: dict[str, Any],
    extra_modules: list[str],
    module_args: dict[str, dict[str, Any]],
    reference_assembly: str,
    *,
    log_errors: bool,
) -> AdvntrPreflight:
    if "advntr" not in extra_modules:
        return AdvntrPreflight(enabled=False, reference=None)

    reference = module_args.get("advntr", {}).get("advntr_reference")
    if reference is not None and not isinstance(reference, str):
        _refuse_invalid_preflight(f"Invalid advntr_reference: {reference}", log_error=log_errors)
    if not reference:
        reference = select_advntr_reference(config, reference_assembly)
    elif reference in {"hg19", "hg38"}:
        reference = config.get("reference_data", {}).get(f"advntr_reference_vntr_{reference}")
    else:
        _refuse_invalid_preflight(f"Invalid advntr_reference: {reference}", log_error=log_errors)

    if not reference:
        _refuse_invalid_preflight("adVNTR reference path not found in configuration.", log_error=log_errors)
    if not isinstance(reference, str) or not reference.strip():
        _refuse_invalid_preflight("adVNTR reference path must be a non-empty string.", log_error=log_errors)

    return AdvntrPreflight(enabled=True, reference=reference)


def plan_advntr_preflight(
    config: dict[str, Any],
    extra_modules: list[str],
    module_args: dict[str, dict[str, Any]],
    reference_assembly: str,
) -> AdvntrPreflight:
    """Resolve the optional adVNTR model path without touching the filesystem.

    Args:
        config: Complete pipeline configuration.
        extra_modules: Requested optional pipeline modules.
        module_args: Per-module configuration overrides.
        reference_assembly: Assembly selected for the pipeline run.

    Returns:
        Typed enablement and resolved model-reference plan.

    Raises:
        ValueError: If an override or resolved reference is invalid.
    """
    return _plan_advntr_preflight(
        config,
        extra_modules,
        module_args,
        reference_assembly,
        log_errors=True,
    )


def plan_valid_advntr_preflight(
    config: dict[str, Any],
    extra_modules: list[str],
    module_args: dict[str, dict[str, Any]],
    reference_assembly: str,
) -> AdvntrPreflight | None:
    """Resolve valid adVNTR state without taking over strict validation.

    The CLI log guard runs before logging setup. Invalid override or configuration
    state belongs to the later strict pipeline planner, so this variant returns None
    without logging rather than moving that failure earlier.

    Args:
        config: Complete pipeline configuration.
        extra_modules: Requested optional pipeline modules.
        module_args: Per-module configuration overrides.
        reference_assembly: Assembly selected for the pipeline run.

    Returns:
        A valid preflight plan, or None when strict validation must decide later.
    """
    try:
        return _plan_advntr_preflight(
            config,
            extra_modules,
            module_args,
            reference_assembly,
            log_errors=False,
        )
    except (AttributeError, TypeError, ValueError):
        return None
