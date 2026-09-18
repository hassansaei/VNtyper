"""Bounded component cutoff grids with genuine Kestrel inclusion relaxations."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import product
from types import MappingProxyType

from vntyper.scripts.calibration_caller_policy import (
    CallerPolicyScalar,
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)

_FLOOR = "/components/kestrel/confidence_assignment/reporting_floor"
_LOW = "/components/kestrel/confidence_assignment/depth_score_thresholds/low"
_GG = "/components/kestrel/alt_filtering/gg_depth_score_threshold"
_MODE = "/components/advntr/calibrated_calling/mode"
_CUTOFF = "/components/advntr/calibrated_calling/cutoff"
_SUPPORT = "/components/advntr/calibrated_calling/minimum_read_support"
_DEFAULT_K = {
    "reporting_floor": [0.001, 0.002, 0.003, 0.00469, 0.00515, 0.006, 0.008, 0.01, 0.015, 0.02, 0.03],
    "gg_depth_score_threshold": [0.002, 0.00469, 0.008, 0.015, 0.03],
}
_DEFAULT_A = {"cutoff": [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05], "minimum_read_support": [2, 3, 5, 8]}


@dataclass(frozen=True)
class CutoffCandidate:
    """A full validated policy and its changed component values."""

    candidate_id: str
    policy: CallerPolicyValues
    parameters: Mapping[str, CallerPolicyScalar]


@dataclass(frozen=True)
class CutoffGrid:
    """Independent component inventories; their Cartesian product needs no new replay.

    Attributes:
        kestrel: Kestrel candidates, baseline first.
        advntr: adVNTR candidates, baseline first, empty when the arm is skipped or
            the baseline does not require adVNTR.
        advntr_skipped_reason: Why a required adVNTR arm is empty, or ``None`` when
            nothing was skipped. An empty arm with no reason means the baseline
            simply does not require adVNTR.
    """

    kestrel: tuple[CutoffCandidate, ...]
    advntr: tuple[CutoffCandidate, ...]
    advntr_skipped_reason: str | None


def _values(value: object, *, integer: bool, strict: bool, cap: int) -> tuple[int | float, ...]:
    if not isinstance(value, list) or not value or len(value) > cap:
        raise ValueError("cutoff grid values must be nonempty bounded lists")
    parsed: list[int | float] = []
    for item in value:
        if type(item) not in {int, float} or (isinstance(item, float) and not math.isfinite(item)):
            raise ValueError("cutoff grid values must be finite numbers, never booleans")
        if integer:
            if type(item) is not int or item < 1:
                raise ValueError("cutoff minimum read support must be a positive integer")
        elif (strict and not 0 < item < 1) or (not strict and not 0 <= item <= 1):
            raise ValueError("cutoff fractions lie in [0,1], or strictly (0,1) for adVNTR")
        parsed.append(item if integer else float(item))
    return tuple(sorted(set(parsed)))


def _component(value: object, names: tuple[str, str], cap: int) -> tuple[tuple[int | float, ...], ...]:
    if not isinstance(value, dict) or set(value) != set(names):
        raise ValueError("cutoff component fields differ from the closed schema")
    return tuple(
        _values(value[name], integer=name == "minimum_read_support", strict=name == "cutoff", cap=cap) for name in names
    )


def _candidates(
    baseline: CallerPolicyValues,
    changes: list[dict[str, CallerPolicyScalar]],
    prefix: str,
) -> tuple[CutoffCandidate, ...]:
    entries = [CutoffCandidate("baseline", baseline, MappingProxyType({}))]
    seen = {baseline.sha256}
    raw = caller_policy_values_document(baseline)
    for change in changes:
        values = dict(baseline.values)
        values.update(change)
        policy = decode_caller_policy_values({**raw, "values": values})
        if policy.sha256 in seen:
            continue
        seen.add(policy.sha256)
        parameters = MappingProxyType({key: value for key, value in change.items() if baseline.values[key] != value})
        entries.append(CutoffCandidate(f"{prefix}-{len(entries):04d}", policy, parameters))
    return tuple(entries)


def build_cutoff_grid(
    baseline: CallerPolicyValues,
    grid_document: object | None = None,
    *,
    max_candidates: int = 4096,
) -> CutoffGrid:
    """Build separately replayable Kestrel and adVNTR cutoff candidates.

    Args:
        baseline: Complete typed policy from the actual capture; always retained.
        grid_document: Optional closed calibration-cutoff-grid-v1 JSON object.
        max_candidates: Positive bound on the combined component Cartesian product.

    Returns:
        Component policies with stable IDs and immutable changed-value metadata.
        Kestrel low confidence boundary moves down with a lower reporting floor;
        the high boundary and other confidence partitions remain at baseline.
        adVNTR candidates keep the baseline calibrated-calling mode. Only a legacy
        baseline is replayable, so a non-legacy baseline yields an empty adVNTR arm
        and a stated ``advntr_skipped_reason`` instead of unusable candidates. The
        requested adVNTR values are still validated and counted against the cap.

    Raises:
        ValueError: For invalid policies, grid values/schema or an excessive product.
    """
    if type(max_candidates) is not int or max_candidates < 1:
        raise ValueError("maximum cutoff candidates must be a positive integer")
    caller_policy_values_document(baseline)
    has_advntr = "advntr" in baseline.required_callers
    if grid_document is None:
        grid_document = {
            "schema_version": "calibration-cutoff-grid-v1",
            "kestrel": _DEFAULT_K,
            **({"advntr": _DEFAULT_A} if has_advntr else {}),
        }
    expected = {"schema_version", "kestrel"} | ({"advntr"} if has_advntr else set())
    if not isinstance(grid_document, dict) or set(grid_document) != expected:
        raise ValueError("cutoff grid fields differ from the baseline caller inventory")
    if grid_document["schema_version"] != "calibration-cutoff-grid-v1":
        raise ValueError("cutoff grid schema_version is unsupported")
    floors, gg_values = _component(
        grid_document["kestrel"], ("reporting_floor", "gg_depth_score_threshold"), max_candidates
    )
    ad_values = (
        _component(grid_document["advntr"], ("cutoff", "minimum_read_support"), max_candidates) if has_advntr else ()
    )
    # Check lower bounds before materializing products, then exact deduplicated sizes.
    k_size = len(floors) * len(gg_values)
    a_size = len(ad_values[0]) * len(ad_values[1]) if has_advntr else 1
    if k_size * a_size > max_candidates:
        raise ValueError("combined cutoff candidate product exceeds the configured cap")
    baseline_low = baseline.values[_LOW]
    assert isinstance(baseline_low, (int, float))
    kestrel = _candidates(
        baseline,
        [{_FLOOR: floor, _LOW: min(baseline_low, floor), _GG: gg} for floor, gg in product(floors, gg_values)],
        "kestrel",
    )
    # The adVNTR replay refuses any candidate whose mode differs from the captured
    # baseline, so candidates inherit the baseline mode rather than forcing legacy.
    baseline_mode = baseline.values[_MODE] if has_advntr else None
    requested = (
        _candidates(
            baseline,
            [
                {_MODE: baseline_mode, _CUTOFF: float(cutoff), _SUPPORT: support}
                for cutoff, support in product(*ad_values)
            ],
            "advntr",
        )
        if has_advntr
        else ()
    )
    # The cap bounds the requested product, whether or not the adVNTR arm is replayable.
    if len(kestrel) * max(1, len(requested)) > max_candidates:
        raise ValueError("combined cutoff candidate product including baseline exceeds the configured cap")
    skipped = (
        None
        if not has_advntr or baseline_mode == "legacy"
        else (
            f"adVNTR baseline calibrated_calling mode is {baseline_mode!r}, not 'legacy'; the adVNTR replay refuses "
            "candidates whose mode differs from the captured baseline, so no adVNTR candidate was emitted"
        )
    )
    return CutoffGrid(kestrel, () if skipped is not None else requested, skipped)
