"""Truth-bound caller observation arms built from replayed cutoff-grid endpoints.

Replay evidence states what a caller decided under one policy; the cohort manifest
states what is true. This module is the only place the two meet, and its whole job is
to keep apart three outcomes that a single Boolean cannot express: a positive call, a
complete negative call, and evidence that was never assessable. A nonempty candidate
population in which every candidate was filtered out is a genuine negative. An empty
capture, or an adVNTR locus whose audit failed, is unavailable, and is never coerced
into a negative just because a negative is convenient for a denominator.

Tier-A identity is deliberately absent from every arm built here. Production assigns
that tier downstream of the capture boundary, during nomenclature reconciliation, and
publishes it as ``Nomenclature_Tier`` on the final caller table
(:func:`vntyper.scripts.calibration_run_projection.build_shipped_projection` is what
:func:`vntyper.scripts.calibration_cohort_callers.read_native_observation` reads it
from). No replayed endpoint carries that column, so no replayed identity may claim the
tier; the alternative would be inventing a second tier rule that production does not
use.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_caller_metrics import CallerObservation, validate_caller_observations
from vntyper.scripts.calibration_cohort_manifest import CohortSample
from vntyper.scripts.calibration_cutoff_advntr import (
    AdvntrCutoffGridResult,
    AdvntrCutoffPolicyResult,
    AdvntrCutoffSample,
)
from vntyper.scripts.calibration_cutoff_kestrel import KestrelGridObservation, KestrelGridReplay

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "calibration-cutoff-observations-v1"

#: Replayed evidence cannot establish the production identity tier; see the module docstring.
_NO_TIER_A: tuple[str, ...] = ()


@dataclass(frozen=True)
class EndpointTallies:
    """Per-policy distribution of outcomes that are not binary calls.

    Attributes:
        confidence: Policy to assigned confidence label to count, over called endpoints
            only; an endpoint with no selected candidate has no confidence to report.
        flags: Policy to flag text to count, over called endpoints only.
        dispositions: Policy to replay disposition to count, over the whole roster.
    """

    confidence: Mapping[str, Mapping[str, int]]
    flags: Mapping[str, Mapping[str, int]]
    dispositions: Mapping[str, Mapping[str, int]]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _roster(samples: Sequence[CohortSample]) -> dict[str, CohortSample]:
    if not isinstance(samples, (tuple, list)) or not samples:
        _fail("cutoff observations require a nonempty declared cohort roster")
    roster: dict[str, CohortSample] = {}
    for sample in samples:
        if not isinstance(sample, CohortSample):
            _fail("cutoff observations require CohortSample roster members")
        if sample.sample_id in roster:
            _fail(f"cutoff observations roster repeats sample {sample.sample_id}")
        roster[sample.sample_id] = sample
    return dict(sorted(roster.items()))


def _require_same_roster(observed: Sequence[str], roster: Sequence[str], label: str) -> None:
    missing = sorted(set(roster) - set(observed))
    unexpected = sorted(set(observed) - set(roster))
    if missing or unexpected:
        _fail(f"cutoff observations {label} roster mismatch; missing {missing}, unexpected {unexpected}")


def _truth(sample: CohortSample) -> tuple[bool | None, tuple[str, ...] | None]:
    """Return the declared truth label and the identity denominator it supports."""
    if sample.truth_variant is not None and sample.genotype is not True:
        _fail(
            f"cutoff observations sample {sample.sample_id} declares truth_variant "
            f"{sample.truth_variant!r} with genotype {sample.genotype!r}"
        )
    if sample.truth_variant is not None:
        return sample.genotype, (sample.truth_variant,)
    # An established negative has the empty identity tuple; anything else is unavailable.
    return sample.genotype, () if sample.genotype is False else None


def _kestrel_call(row: KestrelGridObservation) -> bool | None:
    if not isinstance(row, KestrelGridObservation):
        _fail("cutoff observations require KestrelGridObservation replay endpoints")
    if row.disposition == "called":
        if row.called_positive is not True:
            _fail("cutoff observations Kestrel called disposition must carry a positive endpoint")
        return True
    if row.disposition == "no-call":
        if row.called_positive is not False:
            _fail("cutoff observations Kestrel complete filtered no-call must be a negative endpoint")
        return False
    if row.disposition == "unassessable-no-candidates":
        # The replay already decided whether a native Negative placeholder proves the run
        # executed and found nothing; that decision is evidence, so it is honoured here.
        if row.called_positive is True:
            _fail("cutoff observations Kestrel empty capture cannot be a positive endpoint")
        return row.called_positive
    _fail(f"cutoff observations Kestrel disposition {row.disposition!r} is unsupported")


def _advntr_call(row: AdvntrCutoffSample) -> bool | None:
    if not isinstance(row, AdvntrCutoffSample):
        _fail("cutoff observations require AdvntrCutoffSample replay outcomes")
    if not row.assessable:
        if row.called_positive is not None:
            _fail("cutoff observations unassessable adVNTR evidence cannot carry a call")
        return None
    if row.called_positive is None:
        _fail("cutoff observations assessable adVNTR evidence requires a Boolean call")
    return row.called_positive


def kestrel_observation_arms(
    replay: KestrelGridReplay,
    samples: Sequence[CohortSample],
) -> dict[str, tuple[CallerObservation, ...]]:
    """One arm per replayed policy, carrying detection and exact-identity endpoints.

    Args:
        replay: Policy-major Kestrel grid replay over the declared capture roster.
        samples: Declared cohort rows, one per independent group.

    Returns:
        Policy ID to the complete arm, every arm covering the roster exactly once in the
        same deterministic group order and agreeing on key, group, truth and truth
        identity. Tier-A identities are always empty; see the module docstring.

    Raises:
        ValueError: If the replay is not a :class:`KestrelGridReplay`, the roster is
            empty, duplicated or malformed, the replay and the roster disagree about
            which samples exist, or an endpoint's disposition contradicts its call.
    """
    roster = _roster(samples)
    if not isinstance(replay, KestrelGridReplay):
        _fail("cutoff observations require a KestrelGridReplay")
    _require_same_roster(tuple(replay.sample_keys), tuple(roster), "Kestrel replay")
    arms: dict[str, tuple[CallerObservation, ...]] = {}
    for policy_id in replay.policy_ids:
        endpoints = replay.observations.get(policy_id)
        if endpoints is None:
            _fail(f"cutoff observations Kestrel replay has no endpoints for policy {policy_id}")
        _require_same_roster(tuple(endpoints), tuple(roster), f"Kestrel policy {policy_id}")
        rows = []
        for key, sample in roster.items():
            endpoint = endpoints[key]
            called = _kestrel_call(endpoint)
            identity = endpoint.canonical_identity
            variants = (identity,) if called is True and identity is not None else ()
            truth_positive, truth_variants = _truth(sample)
            rows.append(
                CallerObservation(key, sample.group_id, truth_positive, truth_variants, called, variants, _NO_TIER_A)
            )
        arms[policy_id] = validate_caller_observations(tuple(rows))
    return arms


def advntr_observation_arms(
    result: AdvntrCutoffGridResult,
    samples: Sequence[CohortSample],
) -> dict[str, tuple[CallerObservation, ...]]:
    """One arm per replayed adVNTR policy, retaining audit failures as unavailable.

    The upstream locus results are opaque and carry no canonical molecular identity, so
    every called adVNTR endpoint has an empty called-identity tuple: the call is
    observed, the identity is simply not established by this evidence.

    Args:
        result: Native adVNTR cutoff grid result over the declared capture roster.
        samples: Declared cohort rows, one per independent group.

    Returns:
        Policy ID to the complete arm, in the same deterministic group order as every
        other arm built from the same roster.

    Raises:
        ValueError: If the result is not an :class:`AdvntrCutoffGridResult`, the roster
            is empty, duplicated or malformed, a policy repeats a sample or a policy ID,
            the roster and the replayed samples disagree, or an outcome's assessability
            contradicts its call.
    """
    roster = _roster(samples)
    if not isinstance(result, AdvntrCutoffGridResult) or not result.policies:
        _fail("cutoff observations require a nonempty AdvntrCutoffGridResult")
    arms: dict[str, tuple[CallerObservation, ...]] = {}
    for policy in result.policies:
        if not isinstance(policy, AdvntrCutoffPolicyResult) or policy.policy_id in arms:
            _fail("cutoff observations require unique AdvntrCutoffPolicyResult policies")
        keys = tuple(row.key for row in policy.samples)
        if len(set(keys)) != len(keys):
            _fail(f"cutoff observations adVNTR policy {policy.policy_id} repeats a sample")
        _require_same_roster(keys, tuple(roster), f"adVNTR policy {policy.policy_id}")
        outcomes = {row.key: row for row in policy.samples}
        rows = []
        for key, sample in roster.items():
            truth_positive, truth_variants = _truth(sample)
            rows.append(
                CallerObservation(
                    key,
                    sample.group_id,
                    truth_positive,
                    truth_variants,
                    _advntr_call(outcomes[key]),
                    (),
                    _NO_TIER_A,
                )
            )
        arms[policy.policy_id] = validate_caller_observations(tuple(rows))
    return arms


def _validated_arms(
    arms: Mapping[str, Sequence[CallerObservation]], label: str
) -> dict[str, tuple[CallerObservation, ...]]:
    if not isinstance(arms, Mapping) or not arms:
        _fail(f"cutoff observations {label} arms require a nonempty mapping")
    validated: dict[str, tuple[CallerObservation, ...]] = {}
    for name, rows in arms.items():
        if not isinstance(name, str) or not name or name.strip() != name:
            _fail(f"cutoff observations {label} policy IDs must be nonempty trimmed strings")
        validated[name] = validate_caller_observations(tuple(rows))
    return validated


def _shared(rows: Sequence[CallerObservation]) -> tuple[tuple[str, str, bool | None, tuple[str, ...] | None], ...]:
    return tuple((row.key, row.group_key, row.truth_positive, row.truth_variants) for row in rows)


def _require_agreement(arms: Sequence[tuple[CallerObservation, ...]]) -> tuple[CallerObservation, ...]:
    expected = _shared(arms[0])
    if any(_shared(rows) != expected for rows in arms):
        _fail("cutoff observations arms require one identical specimen, group and truth roster")
    return arms[0]


def _union_call(left: bool | None, right: bool | None) -> bool | None:
    if left is True or right is True:
        return True
    if left is None or right is None:
        return None
    return False


def _union_row(one: CallerObservation, other: CallerObservation) -> CallerObservation:
    called = _union_call(one.called_positive, other.called_positive)
    variants = tuple(sorted(set(one.called_variants) | set(other.called_variants))) if called is True else ()
    tier = tuple(value for value in sorted(set(one.tier_a_variants) | set(other.tier_a_variants)) if value in variants)
    return CallerObservation(one.key, one.group_key, one.truth_positive, one.truth_variants, called, variants, tier)


def union_observation_arms(
    kestrel: Mapping[str, Sequence[CallerObservation]],
    advntr: Mapping[str, Sequence[CallerObservation]],
) -> dict[str, tuple[CallerObservation, ...]]:
    """Combine both callers' arms so a sample is called when either caller calls it.

    Unknown propagates rather than disappearing: ``None`` with ``True`` is ``True``,
    because one caller having called the sample settles it, while ``None`` with
    ``False`` stays ``None``, because the silent caller might have called it.

    Args:
        kestrel: Kestrel arms over the shared roster, keyed by policy ID.
        advntr: adVNTR arms over the same roster, keyed by policy ID.

    Returns:
        Every Kestrel/adVNTR pairing, keyed ``"<kestrel_id>+<advntr_id>"`` in
        deterministic sorted order, carrying the union of the components' identities.

    Raises:
        ValueError: If either inventory is empty or malformed, the arms do not share one
            identical specimen, group and truth roster, or two pairings would produce
            the same joined arm ID.
    """
    left = _validated_arms(kestrel, "Kestrel")
    right = _validated_arms(advntr, "adVNTR")
    _require_agreement([*left.values(), *right.values()])
    arms: dict[str, tuple[CallerObservation, ...]] = {}
    for kestrel_id in sorted(left):
        for advntr_id in sorted(right):
            arm_id = f"{kestrel_id}+{advntr_id}"
            if arm_id in arms:
                _fail(f"cutoff observations union arm IDs must be unique; {arm_id} is ambiguous")
            rows = tuple(_union_row(*pair) for pair in zip(left[kestrel_id], right[advntr_id], strict=True))
            arms[arm_id] = validate_caller_observations(rows)
    return arms


def _counts(values: Mapping[str, int]) -> Mapping[str, int]:
    return MappingProxyType(dict(sorted(values.items())))


def endpoint_tallies(replay: KestrelGridReplay) -> EndpointTallies:
    """Count the non-binary Kestrel endpoints each replayed policy produced.

    Args:
        replay: Policy-major Kestrel grid replay.

    Returns:
        Immutable per-policy confidence, flag and disposition distributions. Confidence
        and flags cover called endpoints only; dispositions cover the whole roster.

    Raises:
        ValueError: If the replay is not a :class:`KestrelGridReplay`, a policy has no
            endpoints, or an endpoint's disposition contradicts its call.
    """
    if not isinstance(replay, KestrelGridReplay):
        _fail("cutoff observations require a KestrelGridReplay")
    confidence: dict[str, Mapping[str, int]] = {}
    flags: dict[str, Mapping[str, int]] = {}
    dispositions: dict[str, Mapping[str, int]] = {}
    for policy_id in replay.policy_ids:
        endpoints = replay.observations.get(policy_id)
        if endpoints is None:
            _fail(f"cutoff observations Kestrel replay has no endpoints for policy {policy_id}")
        by_confidence: dict[str, int] = {}
        by_flag: dict[str, int] = {}
        by_disposition: dict[str, int] = {}
        for key in replay.sample_keys:
            row = endpoints[key]
            _kestrel_call(row)
            by_disposition[row.disposition] = by_disposition.get(row.disposition, 0) + 1
            if row.confidence is not None:
                by_confidence[row.confidence] = by_confidence.get(row.confidence, 0) + 1
            if row.flag is not None:
                by_flag[row.flag] = by_flag.get(row.flag, 0) + 1
        confidence[policy_id] = _counts(by_confidence)
        flags[policy_id] = _counts(by_flag)
        dispositions[policy_id] = _counts(by_disposition)
    return EndpointTallies(MappingProxyType(confidence), MappingProxyType(flags), MappingProxyType(dispositions))


def _observation_document(row: CallerObservation) -> dict[str, object]:
    return {
        "key": row.key,
        "group_key": row.group_key,
        "truth_positive": row.truth_positive,
        "truth_variants": None if row.truth_variants is None else list(row.truth_variants),
        "called_positive": row.called_positive,
        "called_variants": list(row.called_variants),
        "tier_a_variants": list(row.tier_a_variants),
    }


def observation_arms_document(arms: Mapping[str, Sequence[CallerObservation]]) -> dict[str, object]:
    """Project validated arms as canonical ``calibration-cutoff-observations-v1`` JSON.

    Args:
        arms: Complete arms over one shared roster, keyed by policy ID.

    Returns:
        Fresh JSON-compatible content: the shared roster keys once, then every arm's
        observations in the same order. ``truth_variants`` stays null when the identity
        is unavailable, which is not the empty list an established negative carries.

    Raises:
        ValueError: If the inventory is empty or malformed, an observation is invalid,
            or the arms disagree about key, group, truth or truth identity.
    """
    validated = _validated_arms(arms, "document")
    reference = _require_agreement(list(validated.values()))
    return {
        "schema_version": SCHEMA_VERSION,
        "policy_ids": sorted(validated),
        "keys": [row.key for row in reference],
        "arms": [
            {"policy_id": name, "observations": [_observation_document(row) for row in validated[name]]}
            for name in sorted(validated)
        ],
    }


__all__ = [
    "SCHEMA_VERSION",
    "EndpointTallies",
    "advntr_observation_arms",
    "endpoint_tallies",
    "kestrel_observation_arms",
    "observation_arms_document",
    "union_observation_arms",
]
