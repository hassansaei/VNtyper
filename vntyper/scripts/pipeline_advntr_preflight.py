"""Pure adVNTR enablement and model-reference planning for the pipeline."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

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
    if "advntr" not in extra_modules:
        return AdvntrPreflight(enabled=False, reference=None)

    reference = module_args.get("advntr", {}).get("advntr_reference")
    if not reference:
        reference = select_advntr_reference(config, reference_assembly)
    elif reference in {"hg19", "hg38"}:
        reference = config.get("reference_data", {}).get(f"advntr_reference_vntr_{reference}")
    else:
        msg = f"Invalid advntr_reference: {reference}"
        logger.error(msg)
        raise ValueError(msg)

    if not reference:
        msg = "adVNTR reference path not found in configuration."
        logger.error(msg)
        raise ValueError(msg)
    if not isinstance(reference, str) or not reference.strip():
        msg = "adVNTR reference path must be a non-empty string."
        logger.error(msg)
        raise ValueError(msg)

    return AdvntrPreflight(enabled=True, reference=reference)
