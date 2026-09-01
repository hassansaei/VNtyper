"""Focused projection of typed identity reconciliation onto legacy nomenclature."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from vntyper.scripts.identity_reconciliation import (
    IdentityReconciliationObservation,
    IdentityReconciliationPolicy,
    reconcile_identity_observations,
    select_compatibility_presentation_index,
)
from vntyper.scripts.molecular_identity import MolecularIdentity
from vntyper.scripts.nomenclature_evidence import (
    FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
    low_support_flag_for_source,
)

if TYPE_CHECKING:
    from vntyper.scripts.nomenclature import Nomenclature
    from vntyper.scripts.nomenclature_decision_config import NomenclatureDecisionConfig


@dataclass(frozen=True)
class IdentityAwareNomenclatureResult:
    """Presentation call paired with the typed reconciler's selected identity."""

    call: Nomenclature
    selected_identity: MolecularIdentity | None


def reconcile_identity_aware_nomenclature(
    *calls: Nomenclature,
    support: int | None = None,
    supports: Mapping[str, int | None] | None = None,
    identity_observations: tuple[IdentityReconciliationObservation, ...],
    identity_policy: IdentityReconciliationPolicy,
    decision_config: NomenclatureDecisionConfig | None = None,
) -> IdentityAwareNomenclatureResult:
    """Evaluate typed identity evidence once and return presentation plus identity.

    Args:
        *calls: Legacy presentation calls in their existing order.
        support: Legacy scalar support for compatibility-only callers.
        supports: Legacy support keyed by source.
        identity_observations: Typed observations explicitly bound to calls.
        identity_policy: Source-unit identity reconciliation thresholds.
        decision_config: Explicit resolved nomenclature values for the run.

    Returns:
        One legacy presentation call paired with the typed selected identity.

    Raises:
        ValueError: If typed observations, bindings, or policy are inconsistent.
    """
    from vntyper.scripts.nomenclature import (
        FLAG_CALLER_DISAGREEMENT,
        FLAG_KNOWN_VARIANT,
        FLAG_MOTIF_CONTEXT_DIVERGES,
        FLAG_REPRESENTATION_ONLY,
        KNOWN_VARIANTS,
        Nomenclature,
        _reconcile_legacy,
        _undetermined,
    )

    if not isinstance(identity_policy, IdentityReconciliationPolicy):
        raise ValueError("Typed identity observations require an explicit identity policy")
    if not calls:
        return IdentityAwareNomenclatureResult(_undetermined("unknown", 0, "reconciled", ()), None)
    implicit_binding = len(identity_observations) == len(calls) and all(
        observation.presentation_call_index is None for observation in identity_observations
    )
    call_indices = tuple(
        index if implicit_binding else observation.presentation_call_index
        for index, observation in enumerate(identity_observations)
    )
    bound_indices = tuple(index for index in call_indices if index is not None)
    if tuple(sorted(bound_indices)) != tuple(range(len(calls))):
        raise ValueError("Typed reconciliation requires exactly one identity observation per call")
    if not implicit_binding and any(
        observation.presentation_call_index is None and not observation.is_unbound_advntr_row
        for observation in identity_observations
    ):
        raise ValueError("An unbound observation must be a complete adVNTR row")
    if any(
        observation.display_name != call.name
        or observation.source != call.source
        or observation.event != call.event
        or observation.net_length != call.net_length
        for observation, call_index in zip(identity_observations, call_indices, strict=True)
        if call_index is not None
        for call in (calls[call_index],)
    ):
        raise ValueError("Typed identity observations must match their presentation calls")
    admissible = tuple(
        index
        for index in range(len(calls))
        if any(
            call_index == index and observation.disposition.value == identity_policy.admissible_disposition
            for observation, call_index in zip(identity_observations, call_indices, strict=True)
        )
    )
    decision_indices = admissible or tuple(range(len(calls)))
    decision_calls = tuple(calls[index] for index in decision_indices)
    relative_selection = select_compatibility_presentation_index(decision_calls)
    if relative_selection is None:  # pragma: no cover - calls was established above
        return IdentityAwareNomenclatureResult(_undetermined("unknown", 0, "reconciled", ()), None)
    compatibility_call_index = decision_indices[relative_selection]
    compatibility_observation_index = next(
        index for index, call_index in enumerate(call_indices) if call_index == compatibility_call_index
    )
    compatibility_presentation = _reconcile_legacy(
        *decision_calls,
        support=support,
        supports=supports,
        decision_config=decision_config,
    )
    result = reconcile_identity_observations(
        identity_observations,
        identity_policy,
        presentation_selected_observation_index=compatibility_observation_index,
    )
    selected_observation_index = result.presentation_selected_observation_index
    if selected_observation_index is None:  # pragma: no cover - calls established a compatibility selection
        raise ValueError("Typed reconciliation did not select a presentation observation")
    selected_call_index = call_indices[selected_observation_index]
    if selected_call_index is None:
        raise ValueError("Typed reconciliation selected an unbound row observation for presentation")
    presentation = (
        compatibility_presentation if selected_call_index == compatibility_call_index else calls[selected_call_index]
    )
    flags = set().union(*(call.flags for call in decision_calls))
    if result.caller_disagreement or (result.identity is None and FLAG_CALLER_DISAGREEMENT in presentation.flags):
        flags.add(FLAG_CALLER_DISAGREEMENT)
    if any(identity_observations[index].translation.context_diverges for index in result.backing_observation_indices):
        flags.add(FLAG_MOTIF_CONTEXT_DIVERGES)
    for source in result.low_support_sources:
        flags.add(low_support_flag_for_source(source))
    if result.identity is None:
        flags.update(compatibility_presentation.flags)
    if FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT in presentation.flags:
        flags.add(FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT)
    known_variants = decision_config.known_variants if decision_config is not None else KNOWN_VARIANTS
    if presentation.name is not None:
        flags.add(FLAG_KNOWN_VARIANT if presentation.name in known_variants else FLAG_REPRESENTATION_ONLY)
    if result.event_disagreement or presentation.name is None:
        call = _undetermined("unknown", presentation.net_length, "reconciled", tuple(sorted(flags)))
    else:
        tier = "B" if result.tier == "abstained" else result.tier
        call = Nomenclature(
            name=presentation.name,
            event=presentation.event,
            unit=presentation.unit,
            tier=tier,
            flags=tuple(sorted(flags)),
            ambiguity=presentation.ambiguity,
            repeat_form=presentation.repeat_form,
            net_length=presentation.net_length,
            source="reconciled" if len(result.backing_sources) > 1 else presentation.source,
        )
    return IdentityAwareNomenclatureResult(call, result.identity)
