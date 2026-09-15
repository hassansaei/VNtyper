"""Frozen per-run decision components resolved once at the CLI boundary."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal

from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    resolve_decision_profile,
    resolve_research_decision_profile,
)

if TYPE_CHECKING:
    from vntyper.scripts.pipeline_caller_configuration import CallerPipelineConfiguration

StageName = Literal["kestrel", "advntr", "shark", "nomenclature", "cross_match", "dominance"]

_KESTREL_RUNTIME_PATH = Path(__file__).with_name("kestrel_config.json")
_ADVNTR_RUNTIME_PATH = Path(__file__).parents[1] / "modules" / "advntr" / "advntr_config.json"
_SHARK_RUNTIME_PATH = Path(__file__).parents[1] / "modules" / "shark" / "shark_config.json"


def _freeze(value: object) -> object:
    """Recursively freeze JSON-compatible mappings and arrays."""
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(child) for child in value)
    return value


@dataclass(frozen=True)
class RunConfiguration:
    """One immutable resolved profile and its stage-specific components."""

    decision_profile: ResolvedDecisionProfile
    kestrel: Mapping[str, object]
    advntr: Mapping[str, object]
    shark: Mapping[str, object]
    nomenclature: Mapping[str, object]
    cross_match: Mapping[str, object]
    dominance: Mapping[str, object]
    kestrel_runtime: Mapping[str, object]
    advntr_runtime: Mapping[str, object]
    shark_runtime: Mapping[str, object]
    caller_calibration: CallerPipelineConfiguration | None = None


def _load_runtime_components() -> dict[str, Mapping[str, object]]:
    """Load excluded runtime and presentation values without decision leaves."""
    kestrel_sidecar = load_strict_json_object(_KESTREL_RUNTIME_PATH.read_bytes())
    advntr_sidecar = load_strict_json_object(_ADVNTR_RUNTIME_PATH.read_bytes())
    shark_sidecar = load_strict_json_object(_SHARK_RUNTIME_PATH.read_bytes())

    advntr_settings = cast_mapping(advntr_sidecar["advntr_settings"])
    return {
        "kestrel": {
            "kestrel_settings": kestrel_sidecar["kestrel_settings"],
            "subthreshold_note": kestrel_sidecar["subthreshold_note"],
        },
        "advntr": {
            "settings": {
                "additional_commands": advntr_settings["additional_commands"],
                "threads": advntr_settings["threads"],
            }
        },
        "shark": shark_sidecar,
    }


def resolve_run_configuration(
    path: str | Path | None = None,
    *,
    research_profile: str | Path | None = None,
    calibration_bundle: str | Path | None = None,
    calibration_context: str | Path | None = None,
) -> RunConfiguration:
    """Resolve and recursively freeze all decision components once.

    Exactly one decision source is permitted: the packaged default, one explicit
    profile, one derived caller research profile, or one approved portable bundle.
    A research profile carries derived cutoffs and no deployment approval, so it is
    admitted through its own argument rather than by relaxing either of the others.

    Args:
        path: Explicit complete decision profile, or None for the package default.
        research_profile: Derived caller research profile, exclusive with the others.
        calibration_bundle: Approved portable caller bundle, exclusive with path.
        calibration_context: Paired explicit applicability context.

    Returns:
        Frozen run configuration.

    Raises:
        ValueError: If more than one decision source is supplied, or if the bundle and
            its context are not paired.
    """
    if (calibration_bundle is None) != (calibration_context is None):
        raise ValueError("calibration bundle and context must be paired together")
    if path is not None and calibration_bundle is not None:
        raise ValueError("decision profile and calibration bundle are exclusive")
    if research_profile is not None and (path is not None or calibration_bundle is not None):
        raise ValueError("research decision profile is exclusive with an explicit profile or calibration bundle")
    caller_calibration = None
    if calibration_bundle is not None and calibration_context is not None:
        from vntyper.scripts.pipeline_caller_configuration import resolve_caller_pipeline_configuration

        caller_calibration = resolve_caller_pipeline_configuration(Path(calibration_bundle), Path(calibration_context))
    if caller_calibration is not None:
        profile = caller_calibration.bundle.profile
    elif research_profile is not None:
        profile = resolve_research_decision_profile(research_profile)
    else:
        profile = resolve_decision_profile(path)
    frozen = {name: _freeze(component) for name, component in profile.components.items()}
    runtime = {name: _freeze(component) for name, component in _load_runtime_components().items()}
    return RunConfiguration(
        decision_profile=profile,
        caller_calibration=caller_calibration,
        kestrel=cast_mapping(frozen["kestrel"]),
        advntr=cast_mapping(frozen["advntr"]),
        shark=cast_mapping(frozen["shark"]),
        nomenclature=cast_mapping(frozen["nomenclature"]),
        cross_match=cast_mapping(frozen["cross_match"]),
        dominance=cast_mapping(frozen["dominance"]),
        kestrel_runtime=cast_mapping(runtime["kestrel"]),
        advntr_runtime=cast_mapping(runtime["advntr"]),
        shark_runtime=cast_mapping(runtime["shark"]),
    )


def resolve_compatibility_component(
    stage: StageName,
    resolved_component: Mapping[str, object] | None,
    *,
    custom_context_active: bool,
) -> Mapping[str, object]:
    """Resolve a stage component without silently mixing custom and packaged policy.

    Args:
        stage: Decision-profile component name.
        resolved_component: Explicit component already resolved for the run.
        custom_context_active: Whether the caller is operating inside a custom run.

    Returns:
        The explicit component, or the packaged component for a legacy direct caller.

    Raises:
        ValueError: If a custom context omits its explicit component.
    """
    if resolved_component is not None:
        return resolved_component
    if custom_context_active:
        display = {"advntr": "adVNTR", "cross_match": "cross-match"}.get(stage, stage.capitalize())
        raise ValueError(f"custom {display} run context requires an explicit resolved component")
    return cast_mapping(getattr(resolve_run_configuration(), stage))


def resolve_compatibility_runtime_component(
    stage: Literal["kestrel", "advntr", "shark"],
    runtime_component: Mapping[str, object] | None,
) -> Mapping[str, object]:
    """Return an explicit frozen runtime component or the packaged compatibility value."""
    if runtime_component is not None:
        return runtime_component
    return cast_mapping(getattr(resolve_run_configuration(), f"{stage}_runtime"))


def cast_mapping(value: object) -> Mapping[str, object]:
    """Narrow a frozen top-level component to its required mapping type."""
    if not isinstance(value, Mapping):
        raise ValueError("decision profile component must be a mapping")
    return value
