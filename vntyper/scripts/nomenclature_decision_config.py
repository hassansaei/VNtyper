"""Typed projection of one resolved nomenclature decision component."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from vntyper.scripts.identity_reconciliation import IdentityReconciliationPolicy

_COMPONENT_FIELDS = {
    "advntr",
    "canonical_unit",
    "identity_reconciliation",
    "known_variants",
    "motifs",
    "sources",
    "thresholds",
    "unit_length",
}


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    return value


def _positive_integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _string_mapping(value: object, label: str) -> Mapping[str, str]:
    raw = _mapping(value, label)
    if any(
        not isinstance(key, str) or not key or not isinstance(child, str) or not child for key, child in raw.items()
    ):
        raise ValueError(f"{label} must map non-empty strings to non-empty strings")
    return MappingProxyType(dict(raw))  # type: ignore[arg-type]


@dataclass(frozen=True)
class NomenclatureDecisionConfig:
    """Every nomenclature decision value resolved for one run."""

    canonical_unit: str
    unit_length: int
    motifs: Mapping[str, str]
    mappable_repeat_units: Mapping[int, str]
    rotation_offset: int
    caller_of: Mapping[str, str]
    known_variants: Mapping[str, str]
    bam_flank: int
    bam_thin_haplotype_record_support: int
    identity_reconciliation: IdentityReconciliationPolicy


def decision_config_from_component(component: Mapping[str, object]) -> NomenclatureDecisionConfig:
    """Validate and type one already resolved nomenclature component.

    Args:
        component: Complete nomenclature projection from a resolved profile.

    Returns:
        Immutable typed nomenclature configuration.

    Raises:
        ValueError: If the component is incomplete or malformed.
    """
    if set(component) != _COMPONENT_FIELDS:
        raise ValueError(
            f"nomenclature component fields differ: expected {sorted(_COMPONENT_FIELDS)}, got {sorted(component)}"
        )
    canonical_unit = component["canonical_unit"]
    if not isinstance(canonical_unit, str) or not canonical_unit:
        raise ValueError("nomenclature canonical_unit must be a non-empty string")
    unit_length = _positive_integer(component["unit_length"], "nomenclature unit_length")
    if len(canonical_unit) != unit_length:
        raise ValueError("nomenclature canonical_unit length must equal unit_length")

    motifs = _string_mapping(component["motifs"], "nomenclature motifs")
    if any(len(sequence) != unit_length for sequence in motifs.values()):
        raise ValueError("every nomenclature motif must equal unit_length")

    advntr = _mapping(component["advntr"], "nomenclature advntr")
    if set(advntr) != {"mappable_repeat_units", "rotation_offset"}:
        raise ValueError("nomenclature advntr fields differ")
    raw_mappable = _string_mapping(advntr["mappable_repeat_units"], "nomenclature mappable repeat units")
    try:
        mappable = MappingProxyType({int(key): value for key, value in raw_mappable.items()})
    except ValueError as error:
        raise ValueError("nomenclature mappable repeat-unit keys must be integers") from error
    rotation_offset = _positive_integer(advntr["rotation_offset"], "nomenclature rotation_offset")

    sources = _mapping(component["sources"], "nomenclature sources")
    if set(sources) != {"caller_of"}:
        raise ValueError("nomenclature sources fields differ")
    caller_of = _string_mapping(sources["caller_of"], "nomenclature caller_of")
    known_variants = _string_mapping(component["known_variants"], "nomenclature known_variants")

    thresholds = _mapping(component["thresholds"], "nomenclature thresholds")
    expected_thresholds = {"bam_flank", "bam_thin_haplotype_record_support"}
    if set(thresholds) != expected_thresholds:
        raise ValueError("nomenclature threshold fields differ")

    identity = _mapping(component["identity_reconciliation"], "nomenclature identity_reconciliation")
    return NomenclatureDecisionConfig(
        canonical_unit=canonical_unit,
        unit_length=unit_length,
        motifs=motifs,
        mappable_repeat_units=mappable,
        rotation_offset=rotation_offset,
        caller_of=caller_of,
        known_variants=known_variants,
        bam_flank=_positive_integer(thresholds["bam_flank"], "nomenclature bam_flank"),
        bam_thin_haplotype_record_support=_positive_integer(
            thresholds["bam_thin_haplotype_record_support"],
            "nomenclature bam_thin_haplotype_record_support",
        ),
        identity_reconciliation=IdentityReconciliationPolicy.from_component(identity, caller_of=caller_of),
    )
