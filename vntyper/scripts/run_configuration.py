"""Frozen per-run decision components resolved once at the CLI boundary."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

from vntyper.scripts.decision_profile import ResolvedDecisionProfile, resolve_decision_profile


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


def resolve_run_configuration(path: str | Path | None = None) -> RunConfiguration:
    """Resolve and recursively freeze all decision components once.

    Args:
        path: Explicit complete decision profile, or None for the package default.

    Returns:
        Frozen run configuration.
    """
    profile = resolve_decision_profile(path)
    frozen = {name: _freeze(component) for name, component in profile.components.items()}
    return RunConfiguration(
        decision_profile=profile,
        kestrel=cast_mapping(frozen["kestrel"]),
        advntr=cast_mapping(frozen["advntr"]),
        shark=cast_mapping(frozen["shark"]),
        nomenclature=cast_mapping(frozen["nomenclature"]),
        cross_match=cast_mapping(frozen["cross_match"]),
        dominance=cast_mapping(frozen["dominance"]),
    )


def cast_mapping(value: object) -> Mapping[str, object]:
    """Narrow a frozen top-level component to its required mapping type."""
    if not isinstance(value, Mapping):
        raise ValueError("decision profile component must be a mapping")
    return value
