"""adVNTR legacy cutoff axes: breakpoints observed in native replay, and the checks binding them.

adVNTR's legacy frameshift caller scores a candidate only when its read support reaches
``minimum_read_support`` and calls it when the p-value is strictly below ``cutoff``
(adVNTR ``frameshift_decisions.passes_support`` / ``legacy_call``). A sample is called when
any assessable locus returns a call. So under the policy ``(c, s)``::

    called(sample) = any visit with read_support >= s and pvalue < c

Holding the support at baseline, only each sample's smallest eligible p-value ``q`` matters
and the partition changes only where ``c`` crosses a ``q``; the smallest float calling a
sample with ``q`` is ``nextafter(q, +inf)``. Holding the cutoff at baseline, only each
sample's largest support among visits below the cutoff matters. That rule is used here only
to *choose* which values to test: every tested value is replayed natively, and the
replay-consistency check proves the rule matched the native replay for every tested
candidate. Completeness of the enumeration is a separate, by-construction property that
holds only for an uncapped axis (spec §14.2).

The statistics come from a native replay at a permissive projection of each axis, not from
the capture file: the capture only scores visits its own baseline support admitted.
:func:`derive_advntr_axes` runs those probes in one grid of their own, separate from the
candidate grid, and derives every axis once from all samples and once per outer fold from
that fold's training samples (``calibration_cutoff_folds``). Because the probe and the
candidate grid are two native evaluations, :func:`require_same_evidence` proves afterwards
that both replayed the same capture bytes, records and tool.
"""

from __future__ import annotations

import hashlib
import logging
import math
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from types import MappingProxyType
from typing import Final, NoReturn, cast

from vntyper.modules.advntr.advntr_calibration_policy import advntr_canonical_sha256
from vntyper.modules.advntr.advntr_replay import replay_locus_document
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_observations import decode_advntr_baseline_calls
from vntyper.scripts.calibration_caller_policy import (
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_cutoff_advntr import AdvntrCutoffGridResult, evaluate_advntr_cutoff_grid
from vntyper.scripts.calibration_cutoff_axes import (
    ADVNTR_CUTOFF,
    ADVNTR_MIN_SUPPORT,
    AxisBreakpoints,
    axis_candidates,
    observed_axis,
)
from vntyper.scripts.calibration_cutoff_folds import fold_axis
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate
from vntyper.scripts.calibration_cutoff_observations import union_observation_arms
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)

_PREFIX: Final[str] = "/components/advntr/calibrated_calling/"
_MODE: Final[str] = f"{_PREFIX}mode"
_POINTERS: Final[Mapping[str, str]] = MappingProxyType(
    {ADVNTR_CUTOFF: f"{_PREFIX}cutoff", ADVNTR_MIN_SUPPORT: f"{_PREFIX}minimum_read_support"}
)
#: Loosening-first probe ladders; the policy decoder decides which rung is admissible.
PROBE_LADDERS: Final[Mapping[str, tuple[float | int, ...]]] = MappingProxyType(
    {ADVNTR_CUTOFF: (math.nextafter(1.0, 0.0), 0.5, 0.1, 0.05), ADVNTR_MIN_SUPPORT: (1, 2)}
)
#: Decision dispositions that carry a legacy p-value; every other disposition is skipped.
_SCORED: Final[frozenset[str]] = frozenset({"called", "cutoff"})
#: Why a run stops when its probe and candidate replays cannot be shown to share evidence.
EVIDENCE_UNBOUND: Final[str] = "cutoff optimize adVNTR probe and candidate replays are not bound to the same evidence"
#: How many mismatched samples a parity failure names.
_NAMED_MISMATCHES: Final[int] = 5


@dataclass(frozen=True)
class AdvntrVisit:
    """One scored adVNTR decision visit: the read support and the exact legacy p-value."""

    read_support: int
    pvalue: Fraction


#: Per-sample scored visits; ``None`` marks an unassessable sample.
SampleVisits = Mapping[str, tuple[AdvntrVisit, ...] | None]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _axis_pointer(axis: str) -> str:
    if axis not in _POINTERS:
        _fail(f"adVNTR cutoff axis must be one of {sorted(_POINTERS)}")
    return _POINTERS[axis]


def _require_legacy(policy: CallerPolicyValues, role: str = "baseline") -> None:
    if "advntr" not in policy.required_callers:
        _fail(f"adVNTR cutoff axes require a {role} policy that includes adVNTR")
    mode = policy.values[_MODE]
    if mode != "legacy":
        _fail(f"adVNTR cutoff axes require the legacy calibrated_calling mode; the {role} mode is {mode!r}")


def _baseline_fraction(baseline: CallerPolicyValues, axis: str) -> Fraction:
    return Fraction(cast(int | float, baseline.values[_POINTERS[axis]]))


def _with_value(baseline: CallerPolicyValues, pointer: str, value: float | int) -> CallerPolicyValues:
    document = caller_policy_values_document(baseline)
    values = dict(baseline.values)
    values[pointer] = value
    return decode_caller_policy_values({**document, "values": values})


def advntr_probe_policy(axis: str, baseline: CallerPolicyValues) -> CallerPolicyValues:
    """The most permissive admissible projection of one adVNTR axis.

    Only the axis pointer moves; every other value, ``mode`` included, keeps its baseline
    value. The rungs of :data:`PROBE_LADDERS` are tried loosening-first and the first one the
    policy decoder accepts is returned.

    Args:
        axis: :data:`ADVNTR_CUTOFF` or :data:`ADVNTR_MIN_SUPPORT`.
        baseline: The complete legacy baseline policy.

    Returns:
        The decoded probe policy.

    Raises:
        ValueError: For an unknown axis, a non-legacy or adVNTR-less baseline, or when the
            decoder refuses every rung of the probe ladder.
    """
    pointer = _axis_pointer(axis)
    _require_legacy(baseline)
    for rung in PROBE_LADDERS[axis]:
        try:
            return _with_value(baseline, pointer, rung)
        except ValueError:
            continue
    _fail(f"adVNTR cutoff axis {axis} has no admissible permissive projection")


def _visit(value: object) -> AdvntrVisit | None:
    if not isinstance(value, Mapping):
        _fail("adVNTR replay decision visit must be an object")
    if value.get("disposition") not in _SCORED:
        return None
    statistic, plan = value.get("statistic"), value.get("plan")
    if not isinstance(statistic, Mapping) or not isinstance(plan, Mapping):
        _fail("adVNTR replay decision visit statistic and plan must be objects")
    support, pvalue = plan.get("read_support"), statistic.get("pvalue")
    if isinstance(support, bool) or not isinstance(support, int) or support < 0:
        _fail("adVNTR replay decision visit read support must be a non-negative integer")
    if pvalue is None or isinstance(pvalue, bool) or not isinstance(pvalue, (int, float)):
        _fail("adVNTR replay scored visit must carry a numeric p-value")
    if not math.isfinite(pvalue) or not 0 <= pvalue <= 1:
        _fail("adVNTR replay scored visit p-value must be finite and lie in [0, 1]")
    return AdvntrVisit(support, Fraction(pvalue))


def probe_visits(result: AdvntrCutoffGridResult, policy_id: str) -> dict[str, tuple[AdvntrVisit, ...] | None]:
    """Every scored visit of one replayed policy, per sample.

    Only the dispositions ``called`` and ``cutoff`` carry a legacy p-value; visits with any
    other disposition (``legacy-nonfinite``, ``insufficient-read-support``,
    ``outside-boundary``) are skipped, and a scored visit must be well-formed (spec §14.5).

    The grid records no policy values, so this cannot check the replayed policy's mode:
    callers must replay a legacy-mode policy (:func:`advntr_probe_policy` keeps the
    legacy baseline mode), because an exact-mode p-value has different semantics.
    :func:`check_replay_consistency` checks the mode of every candidate it binds.

    Args:
        result: A native replay grid that contains ``policy_id``.
        policy_id: The replayed probe policy to read.

    Returns:
        Sample key to its scored visits across all loci, or ``None`` for an unassessable
        sample.

    Raises:
        ValueError: If the policy is absent from the grid, or a scored visit is malformed
            (including a ``called``/``cutoff`` visit with a null statistic or p-value).
    """
    matches = [row for row in result.policies if row.policy_id == policy_id]
    if len(matches) != 1:
        _fail(f"adVNTR probe policy {policy_id} is absent from the replay grid")
    visits: dict[str, tuple[AdvntrVisit, ...] | None] = {}
    for sample in matches[0].samples:
        if not sample.assessable:
            visits[sample.key] = None
            continue
        found: list[AdvntrVisit] = []
        for locus in sample.loci:
            raw = replay_locus_document(locus).get("decision_visits")
            if not isinstance(raw, list):
                _fail("adVNTR replay decision visits must be a list")
            found.extend(visit for visit in (_visit(item) for item in raw) if visit is not None)
        visits[sample.key] = tuple(found)
    return visits


def sample_statistics(axis: str, visits: SampleVisits, baseline: CallerPolicyValues) -> dict[str, Fraction]:
    """Each assessable sample's decisive statistic on one axis.

    On the cutoff axis it is the smallest p-value among visits with support at least the
    baseline support; on the support axis, the largest support among visits whose p-value
    is strictly below the baseline cutoff.

    Args:
        axis: :data:`ADVNTR_CUTOFF` or :data:`ADVNTR_MIN_SUPPORT`.
        visits: Per-sample scored visits from :func:`probe_visits`.
        baseline: The complete legacy baseline policy holding the other axis fixed.

    Returns:
        Sample key to its statistic; unassessable samples and samples without an eligible
        visit are absent.

    Raises:
        ValueError: For an unknown axis or a non-legacy or adVNTR-less baseline.
    """
    _axis_pointer(axis)
    _require_legacy(baseline)
    cutoff = _baseline_fraction(baseline, ADVNTR_CUTOFF)
    support = _baseline_fraction(baseline, ADVNTR_MIN_SUPPORT)
    statistics: dict[str, Fraction] = {}
    for key in sorted(visits):
        items = visits[key]
        if items is None:
            continue
        if axis == ADVNTR_CUTOFF:
            eligible = [visit.pvalue for visit in items if visit.read_support >= support]
            if eligible:
                statistics[key] = min(eligible)
        else:
            supports = [Fraction(visit.read_support) for visit in items if visit.pvalue < cutoff]
            if supports:
                statistics[key] = max(supports)
    return statistics


def _is_float_value(value: Fraction) -> bool:
    """Breakpoints and sentinels round through float, so a statistic must be one exactly."""
    try:
        return Fraction(float(value)) == value
    except OverflowError:
        return False


def _next_up(value: Fraction) -> Fraction:
    return Fraction(math.nextafter(float(value), math.inf))


def _cutoff_sentinel(observed: set[Fraction], anchor: Fraction) -> Fraction | None:
    """The cutoff calling no rejectable sample: min(q), or the smallest positive float when min(q) is 0."""
    if not observed:
        return None
    low = min(observed)
    sentinel = low if low > 0 else Fraction(math.nextafter(0.0, 1.0))
    return None if anchor <= sentinel else sentinel


def derive_advntr_axis(
    axis: str, statistics: Mapping[str, Fraction], *, baseline: CallerPolicyValues, max_values: int | None
) -> AxisBreakpoints:
    """Breakpoints of one adVNTR axis from per-sample decisive statistics.

    Cutoff axis: ``nextafter(q, +inf)`` for every distinct ``q``, plus the sentinel
    ``min(q)`` (the smallest positive float when ``min(q)`` is 0) unless the baseline
    already rejects every sample. Support axis: every distinct observed support, plus the
    sentinel ``max + 1`` unless the baseline already exceeds it. The baseline value is
    always kept, and the decoder screens every value.

    Args:
        axis: :data:`ADVNTR_CUTOFF` or :data:`ADVNTR_MIN_SUPPORT`.
        statistics: Per-sample exact statistics from :func:`sample_statistics`.
        baseline: The complete legacy baseline policy.
        max_values: Optional cap applied by rank-space subsampling.

    Returns:
        The observed-breakpoint axis.

    Raises:
        ValueError: For an unknown axis, a non-legacy or adVNTR-less baseline, a statistic
            that is not an exact float value, or an invalid cap.
    """
    _axis_pointer(axis)
    _require_legacy(baseline)
    observed = set(statistics.values())
    if any(not isinstance(value, Fraction) for value in observed):
        _fail("adVNTR cutoff axis statistics must be exact Fractions")
    if any(not _is_float_value(value) for value in observed):
        _fail("adVNTR cutoff axis statistics must be exactly representable as floats")
    anchor = _baseline_fraction(baseline, axis)
    sentinel: Fraction | None
    if axis == ADVNTR_CUTOFF:
        candidates = {_next_up(value) for value in observed}
        sentinel = _cutoff_sentinel(observed, anchor)
    else:
        candidates = set(observed)
        top = max(observed) if observed else None
        sentinel = None if top is None or anchor > top else top + 1
    return observed_axis(
        axis,
        sorted(candidates),
        baseline=baseline,
        observed_count=len(observed),
        sentinel=sentinel,
        max_values=max_values,
    )


def unrejectable_samples(axis: str, statistics: Mapping[str, Fraction]) -> int:
    """Samples no admissible cutoff can reject.

    Args:
        axis: :data:`ADVNTR_CUTOFF` or :data:`ADVNTR_MIN_SUPPORT`.
        statistics: Per-sample statistics from :func:`sample_statistics`.

    Returns:
        The count of p-values of exactly zero on the cutoff axis (a cutoff must exceed 0);
        always 0 on the support axis.

    Raises:
        ValueError: For an unknown axis.
    """
    _axis_pointer(axis)
    return sum(1 for value in statistics.values() if value == 0) if axis == ADVNTR_CUTOFF else 0


def predicted_call(visits: tuple[AdvntrVisit, ...] | None, cutoff: float, support: int) -> bool | None:
    """The legacy sample call the module rule predicts.

    Args:
        visits: One sample's scored visits, or ``None`` for an unassessable sample.
        cutoff: The legacy cutoff; a visit calls when its p-value is strictly below it.
        support: The minimum read support; a visit is eligible when its support reaches it.

    Returns:
        Whether any eligible visit calls, or ``None`` for an unassessable sample.
    """
    if visits is None:
        return None
    threshold = Fraction(cutoff)
    return any(visit.read_support >= support and visit.pvalue < threshold for visit in visits)


def _moved_axis(candidate: CutoffCandidate) -> str:
    moved = [axis for axis, pointer in _POINTERS.items() if pointer in candidate.parameters]
    return ", ".join(moved) if moved else "baseline"


def check_replay_consistency(
    result: AdvntrCutoffGridResult, candidates: Sequence[CutoffCandidate], visits: SampleVisits
) -> int:
    """Prove every candidate's native calls equal the calls the module rule predicts.

    An unassessable sample is consistent only when the replay also leaves it uncalled
    (``None``); it never counts as a call either way.

    Args:
        result: The native replay grid holding every candidate's calls.
        candidates: The tested candidates of one axis; each must be a legacy adVNTR policy.
        visits: The per-sample scored visits the axis was derived from.

    Returns:
        The number of candidates checked.

    Raises:
        ValueError: If a candidate policy is not legacy adVNTR, was not replayed, or was
            replayed over a different roster, or naming the axis, the candidate, its tested
            cutoff and support, and the number of disagreeing samples.
    """
    native = {row.policy_id: {sample.key: sample.called_positive for sample in row.samples} for row in result.policies}
    for candidate in candidates:
        _require_legacy(candidate.policy, f"candidate {candidate.candidate_id}")
        calls = native.get(candidate.candidate_id)
        if calls is None:
            _fail(f"adVNTR replay consistency: candidate {candidate.candidate_id} was not replayed")
        if set(calls) != set(visits):
            _fail(f"adVNTR replay consistency: candidate {candidate.candidate_id} replayed a different roster")
        cutoff = cast(float, candidate.policy.values[_POINTERS[ADVNTR_CUTOFF]])
        support = cast(int, candidate.policy.values[_POINTERS[ADVNTR_MIN_SUPPORT]])
        wrong = [key for key in sorted(visits) if predicted_call(visits[key], cutoff, support) != calls[key]]
        if wrong:
            _fail(
                f"adVNTR replay consistency failed on axis {_moved_axis(candidate)} for candidate "
                f"{candidate.candidate_id} (cutoff={cutoff!r}, support={support}): the native replay disagrees "
                f"with the legacy decision rule on {len(wrong)} samples, so the replayed operating points "
                "cannot be trusted"
            )
    return len(candidates)


def _capture_records(raw: bytes) -> tuple[tuple[int, str], ...]:
    """The capture's ``(VNTR identifier, canonical record digest)`` pairs, in record order."""
    records = []
    for line in raw.splitlines():
        document = load_strict_json_object(line)
        locus = document.get("locus")
        vntr_id = locus.get("vntr_id") if isinstance(locus, Mapping) else None
        if isinstance(vntr_id, bool) or not isinstance(vntr_id, int) or vntr_id <= 0:
            _fail("adVNTR capture record lacks a positive VNTR identifier")
        records.append((vntr_id, advntr_canonical_sha256(document)))
    return tuple(records)


def advntr_baseline_parity(
    result: AdvntrCutoffGridResult, anchor_id: str, capture_paths: Mapping[str, Path]
) -> dict[str, object]:
    """Require the anchor replay to reproduce each capture's own native baseline calls.

    Each capture's recorded decisions are decoded independently of the replay; a sample
    whose capture failed its calibration audit is unassessable (``None``) and must be
    unassessable in the replay too. Each capture record is also bound to the replay: its
    canonical digest must equal the ``capture_record_sha256`` the anchor replayed for that
    sample and locus, so the parity is proven against the very bytes that were replayed.

    Args:
        result: The native replay grid holding the anchor policy.
        anchor_id: The replayed policy ID of the baseline anchor.
        capture_paths: Sample key to its capture file; the roster must equal the replay's.

    Returns:
        ``{"proven": True, "sample_count": n, "mismatches": []}``.

    Raises:
        ValueError: If the anchor was not replayed, the rosters or a sample's loci differ,
            a capture is malformed or is not the record the anchor replayed, or any
            sample's replayed call differs from its capture (naming the first few samples).
    """
    rows = [row for row in result.policies if row.policy_id == anchor_id]
    if len(rows) != 1:
        _fail(f"adVNTR baseline parity: anchor {anchor_id} was not replayed")
    samples = {sample.key: sample for sample in rows[0].samples}
    if set(samples) != set(capture_paths):
        _fail("adVNTR baseline parity: the replay and the captures disagree about the roster")
    mismatches = []
    for key in sorted(capture_paths):
        raw = read_regular_path(capture_paths[key])
        records = _capture_records(raw)
        vntr_ids = tuple(vntr_id for vntr_id, _ in records)
        loci = {locus.vntr_id: locus for locus in samples[key].loci}
        if sorted(vntr_ids) != sorted(locus.vntr_id for locus in samples[key].loci):
            _fail(f"adVNTR baseline parity: the capture for {key} covers different loci than its replay")
        for vntr_id, digest in records:
            if replay_locus_document(loci[vntr_id]).get("capture_record_sha256") != digest:
                _fail(
                    f"adVNTR baseline parity: the capture record for {key} at VNTR {vntr_id} is not the record "
                    "the anchor replayed"
                )
        calls, assessable = decode_advntr_baseline_calls(raw, vntr_ids)
        native = bool(calls) if assessable else None
        if native != samples[key].called_positive:
            mismatches.append(key)
    if mismatches:
        _fail(
            "adVNTR baseline parity failed: the anchor replay disagrees with the native capture decisions "
            f"for {len(mismatches)} samples (first: {', '.join(mismatches[:_NAMED_MISMATCHES])})"
        )
    return {"proven": True, "sample_count": len(samples), "mismatches": []}


def capture_snapshot(capture_paths: Mapping[str, Path]) -> dict[str, str]:
    """The exact content digest of every capture, taken before any adVNTR replay.

    Args:
        capture_paths: Sample key to its capture file.

    Returns:
        Sample key to the SHA-256 of its capture bytes.
    """
    return {key: hashlib.sha256(read_regular_path(path)).hexdigest() for key, path in sorted(capture_paths.items())}


def _record_digests(result: AdvntrCutoffGridResult) -> dict[tuple[str, int], str]:
    """Every replayed ``(sample, locus)`` record digest; the policies of one grid must agree."""
    digests: dict[tuple[str, int], str] = {}
    for row in result.policies:
        for sample in row.samples:
            for locus in sample.loci:
                digest = replay_locus_document(locus).get("capture_record_sha256")
                if not isinstance(digest, str) or digests.setdefault((sample.key, locus.vntr_id), digest) != digest:
                    _fail(EVIDENCE_UNBOUND)
    return digests


def require_same_evidence(
    probe: AdvntrCutoffGridResult,
    main: AdvntrCutoffGridResult,
    snapshot: Mapping[str, str],
    capture_paths: Mapping[str, Path],
) -> None:
    """Prove the probe grid and the candidate grid replayed the same evidence with the same tool.

    The axis statistics come from the probe grid and the scored arms from the candidate grid,
    so the replay-consistency check binds the two only if both saw identical inputs.

    Args:
        probe: The probe grid the axes were derived from.
        main: The candidate grid the arms are scored from.
        snapshot: :func:`capture_snapshot` taken before the probe grid ran.
        capture_paths: The same sample-to-capture mapping, re-hashed now.

    Raises:
        ValueError: If the tool identities, capture policies or per-locus record digests
            differ, or a capture's bytes changed since the snapshot.
    """
    if (
        probe.capabilities != main.capabilities
        or probe.capture_policy_sha256 != main.capture_policy_sha256
        or _record_digests(probe) != _record_digests(main)
        or capture_snapshot(capture_paths) != dict(snapshot)
    ):
        _fail(EVIDENCE_UNBOUND)


@dataclass(frozen=True)
class AdvntrAxisSearch:
    """The adVNTR axes of one optimize run, derived from their probe replays.

    Attributes:
        derived: Each merged axis (full data plus every fold) with its candidates.
        fold_values: Axis name to outer fold to the values that fold's training samples derived.
        visits: Axis name to the probe's per-sample scored visits.
        unrejectable: Axis name to the count of samples no admissible value can reject.
        probe: The probe grid (baseline plus one permissive projection per axis).
        probe_seconds: Wall time of the probe grid.
    """

    derived: tuple[tuple[AxisBreakpoints, tuple[CutoffCandidate, ...]], ...]
    fold_values: Mapping[str, Mapping[int, frozenset[int | float]]]
    visits: Mapping[str, SampleVisits]
    unrejectable: Mapping[str, int]
    probe: AdvntrCutoffGridResult
    probe_seconds: float


def derive_advntr_axes(
    names: Sequence[str],
    baseline: CallerPolicyValues,
    capture_paths: Mapping[str, Path],
    *,
    executable_path: Path,
    output: Path,
    max_values: int | None,
    assignments: Mapping[str, int],
) -> AdvntrAxisSearch:
    """Replay the probes once and derive each adVNTR axis, full-data and fold-local.

    Args:
        names: The requested adVNTR axes, in order.
        baseline: The complete legacy baseline policy.
        capture_paths: Sample key to its adVNTR capture.
        executable_path: The pinned installed adVNTR executable.
        output: New private directory for the probe grid.
        max_values: Optional per-inventory breakpoint cap.
        assignments: Sample key to outer fold (``outer_fold_assignments``); empty when
            there are too few groups to cross-validate.

    Returns:
        The derived axes, their fold inventories, the probe visits, the unrejectable
        sample counts and the probe grid.

    Raises:
        ValueError: For an unknown axis, a non-legacy or adVNTR-less baseline or probe, or
            any failure of the probe replay or the derivation.
    """
    for name in names:
        _axis_pointer(name)
    probes = {f"probe-{name}": advntr_probe_policy(name, baseline) for name in names}
    for name in names:
        _require_legacy(probes[f"probe-{name}"], f"probe {name}")
    started = time.monotonic()
    probe = evaluate_advntr_cutoff_grid(
        capture_paths,
        {"baseline": baseline, **probes},
        baseline_policy_id="baseline",
        executable_path=executable_path,
        output=output,
    )
    seconds = time.monotonic() - started
    derived: list[tuple[AxisBreakpoints, tuple[CutoffCandidate, ...]]] = []
    fold_values: dict[str, Mapping[int, frozenset[int | float]]] = {}
    visits: dict[str, SampleVisits] = {}
    unrejectable: dict[str, int] = {}
    for name in names:
        visits[name] = probe_visits(probe, f"probe-{name}")
        statistics = sample_statistics(name, visits[name], baseline)

        def derive(stats: Mapping[str, Fraction], axis: str = name) -> AxisBreakpoints:
            return derive_advntr_axis(axis, stats, baseline=baseline, max_values=max_values)

        axis, fold_values[name] = fold_axis(derive, statistics, assignments)
        derived.append((axis, axis_candidates(baseline, axis)))
        unrejectable[name] = unrejectable_samples(name, statistics)
    return AdvntrAxisSearch(tuple(derived), fold_values, visits, unrejectable, probe, seconds)


def combine_arms(
    kestrel: Mapping[str, tuple[CallerObservation, ...]], advntr: Mapping[str, tuple[CallerObservation, ...]]
) -> dict[str, tuple[CallerObservation, ...]]:
    """Pair each policy's Kestrel arm with the same policy's adVNTR arm, and only that one.

    The full Cartesian product of the two inventories would pair a Kestrel policy with an
    adVNTR policy that was never replayed beside it, so only the diagonal is kept.

    Args:
        kestrel: Policy ID to its Kestrel arm.
        advntr: Policy ID to its adVNTR arm; it must cover every Kestrel policy.

    Returns:
        Policy ID to the either-caller union of its two arms.

    Raises:
        ValueError: If a Kestrel policy has no adVNTR arm, or the union refuses the pair.
    """
    combined: dict[str, tuple[CallerObservation, ...]] = {}
    for policy_id, rows in kestrel.items():
        partner = advntr.get(policy_id)
        if partner is None:
            _fail(f"cutoff optimize adVNTR replay produced no arm for policy {policy_id}")
        combined[policy_id] = union_observation_arms({policy_id: rows}, {policy_id: partner})[
            f"{policy_id}+{policy_id}"
        ]
    return combined
