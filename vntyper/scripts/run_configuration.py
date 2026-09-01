"""Frozen per-run decision components resolved once at the CLI boundary."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Literal

from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.decision_profile import ResolvedDecisionProfile, resolve_decision_profile

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


def resolve_run_configuration(path: str | Path | None = None) -> RunConfiguration:
    """Resolve and recursively freeze all decision components once.

    Args:
        path: Explicit complete decision profile, or None for the package default.

    Returns:
        Frozen run configuration.
    """
    profile = resolve_decision_profile(path)
    if profile.profile_kind == "generated":
        packaged = resolve_decision_profile()
        if profile.components["dominance"] != packaged.components["dominance"]:
            raise ValueError(
                "generated dominance policy is not active until its PR-D consumer is installed; use neutral values"
            )
    frozen = {name: _freeze(component) for name, component in profile.components.items()}
    runtime = {name: _freeze(component) for name, component in _load_runtime_components().items()}
    return RunConfiguration(
        decision_profile=profile,
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
