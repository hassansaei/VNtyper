"""Derive caller cutoffs from labelled cohort evidence and export a research profile.

This is the run behind ``vntyper calibrate optimize``. It reads a cohort with truth,
replays the captured evidence over the decision breakpoints the cohort itself produces,
scores every resulting operating point, selects one under an explicitly declared
objective, and writes the selected policy out in the exact form
``vntyper pipeline --research-decision-profile`` accepts. What the published record says
lives in ``calibration_cutoff_document``; the bytes that reach disk live in
``calibration_cutoff_report``.

Three properties are load-bearing and are checked rather than assumed.

*Baseline parity.* The baseline is taken from the captures, because that is the policy
that actually produced the evidence, and the candidate that reproduces it must replay to
the same per-sample outcome. A curve computed on evidence that cannot reproduce the
shipped result describes nothing, so a parity failure aborts the run instead of being
recorded as a caveat.

*The breakpoints come from the data.* A threshold compared with ``>=`` or ``<=`` changes a
decision only at values the data actually take, so the complete set of distinct outcomes
is obtained by testing exactly the observed values. They are discovered by replaying each
capture once at a permissive projection of the axis -- permissive so that no row is lost
before its value can be observed -- and the projection itself is screened through the
policy decoder rather than hard-coded, so a tightened bound moves the probe instead of
producing an inadmissible policy.

*The exported profile round-trips.* The written profile is read back through the runtime
resolver and projected into a caller policy again; a profile that does not recover the
selected policy exactly fails the run rather than being published.

Leakage control happens before anything is scored: only the first-seen representative of
each declared ``group_id`` is retained, so biological duplicates cannot be counted twice
or straddle a training split. Every axis is derived once from all samples and once per
outer fold from that fold's training samples alone (``calibration_cutoff_folds``), and a
fold may select only the baseline and its own training-derived candidates, so no held-out
value can create, remove or displace a threshold its own fold uses. Kestrel captures are
required for every caller selection, because they carry the complete baseline policy of
both callers. The output directory holds cohort results and is written ``0700`` with
``0600`` files throughout.

*adVNTR axes.* ``--caller advntr`` searches the legacy cutoff and read-support axes and
scores adVNTR-only arms; ``--caller both`` searches the union of the Kestrel and adVNTR
axes, one axis at a time, each candidate holding the other caller at baseline, and scores
the either-caller union. adVNTR breakpoints come from a separate probe grid: the baseline
and one permissive projection per axis replayed natively under ``advntr-probe``. Every
candidate is then replayed natively again, and four checks bind the result: the probe and
candidate grids replayed the same capture bytes, records and tool; every candidate's native
calls equal the calls the legacy rule predicts from the probe statistics (replay
consistency); the adVNTR anchor reproduces each capture's own recorded baseline calls
(adVNTR parity); and the candidate grid executed exactly one adVNTR replay per distinct
adVNTR policy -- one in all when only Kestrel axes are searched.
"""

from __future__ import annotations

import hashlib
import logging
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from types import MappingProxyType
from typing import Any, Final, NoReturn, cast

from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_caller_profile import build_caller_generated_profile
from vntyper.scripts.calibration_cohort_manifest import CohortSample, read_cohort_manifest
from vntyper.scripts.calibration_cohort_metrics import group_folds
from vntyper.scripts.calibration_cutoff_advntr import (
    AdvntrCutoffGridResult,
    evaluate_advntr_cutoff_grid,
    read_capture_policy,
    require_one_execution_per_signature,
)
from vntyper.scripts.calibration_cutoff_advntr_axes import (
    AdvntrAxisSearch,
    advntr_baseline_parity,
    capture_snapshot,
    check_replay_consistency,
    combine_arms,
    derive_advntr_axes,
    require_same_evidence,
)
from vntyper.scripts.calibration_cutoff_axes import (
    ACTIVE_REGION,
    ADVNTR_CUTOFF,
    ADVNTR_MIN_SUPPORT,
    ALT_DEPTH_BAND,
    DEPTH_FLOOR_LINKED,
    DEPTH_SCORE_HIGH,
    GG_GATE_INDEPENDENT,
    AxisBreakpoints,
    axis_caller,
    axis_candidates,
    axis_comparison,
    axis_document,
    declared_axis,
    derive_axis,
    eligible_statistic_values,
)
from vntyper.scripts.calibration_cutoff_document import (
    BASELINE_ID,
    PROFILE_NAME,
    USAGE_HINT,
    CutoffReportInputs,
    DerivedAxis,
    build_cutoff_report_document,
)
from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms, outer_fold_assignments
from vntyper.scripts.calibration_cutoff_folds import fold_axis, fold_candidate_inventories
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate
from vntyper.scripts.calibration_cutoff_inputs import primary_samples, read_cutoff_captures
from vntyper.scripts.calibration_cutoff_kestrel import KestrelGridReplay, replay_kestrel_grid
from vntyper.scripts.calibration_cutoff_observations import (
    advntr_observation_arms,
    kestrel_observation_arms,
)
from vntyper.scripts.calibration_cutoff_report import write_cutoff_reports
from vntyper.scripts.calibration_cutoff_selection import SearchSpec
from vntyper.scripts.calibration_kestrel_capture import KestrelCapture, decode_kestrel_capture
from vntyper.scripts.calibration_kestrel_replay import kestrel_replay_prefilter_frame, replay_kestrel_capture
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.decision_profile import resolve_research_decision_profile
from vntyper.scripts.pipeline_research_advntr import (
    RESEARCH_LEGACY_ONLY,
    project_caller_policy,
    research_capture_differences,
    research_policy_argv,
)
from vntyper.version import __version__

logger = logging.getLogger(__name__)

GENERATOR_VERSION: Final[str] = f"vntyper-calibrate-optimize/{__version__}"
_CALLERS: Final[frozenset[str]] = frozenset({"kestrel", "advntr", "both"})
_ADVNTR_MODE: Final[str] = "/components/advntr/calibrated_calling/mode"

#: Per axis: the exact production comparator, and a ladder of extreme values probed in
#: loosening-first order. The first rung the policy decoder accepts is the permissive
#: projection used to observe breakpoints. Nothing here asserts that a rung is admissible
#: -- the decoder is asked -- so a tightened policy bound moves the probe to the next rung
#: instead of producing a policy that cannot be replayed.
AXIS_PROBE: Final[Mapping[str, tuple[str, tuple[float, ...]]]] = MappingProxyType(
    {
        DEPTH_FLOOR_LINKED: (">=", (0.0, 1e-9, 1e-6, 1e-4, 1e-3)),
        GG_GATE_INDEPENDENT: (">=", (0.0, 1e-9, 1e-6, 1e-4, 1e-3)),
        DEPTH_SCORE_HIGH: ("<=", (1.0, 0.5, 0.1, 0.01)),
        ALT_DEPTH_BAND: ("<=", (98, 64, 32, 16, 8, 4, 2, 1, 0)),
        ACTIVE_REGION: ("<=", (1_000_000, 100_000, 10_000, 1_000, 200, 0)),
    }
)

#: Every axis optimize can search: the Kestrel probe axes, then the two adVNTR legacy axes.
_ALL_AXES: Final[tuple[str, ...]] = (*AXIS_PROBE, ADVNTR_CUTOFF, ADVNTR_MIN_SUPPORT)

#: Axis name to its production comparator, which is all the report layer needs.
AXIS_COMPARISON: Final[Mapping[str, str]] = MappingProxyType({name: axis_comparison(name) for name in _ALL_AXES})

#: The axes searched when ``--axis`` is not given, per caller selection.
_DEFAULT_AXES: Final[Mapping[str, tuple[str, ...]]] = MappingProxyType(
    {"kestrel": (DEPTH_FLOOR_LINKED,), "advntr": (ADVNTR_CUTOFF,), "both": (DEPTH_FLOOR_LINKED, ADVNTR_CUTOFF)}
)


@dataclass(frozen=True)
class _Request:
    """One fully validated optimize invocation."""

    manifest: Path
    captures: Path
    spec: SearchSpec
    caller: str
    axes: tuple[str, ...]
    max_breakpoints: int | None
    folds: int
    seed: int
    workers: int
    advntr_executable: Path | None


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _write_private(path: Path, raw: bytes) -> None:
    """Write one artifact byte-exactly with owner-only permissions."""
    path.write_bytes(raw)
    path.chmod(0o600)


def _digest(path: Path) -> str:
    """The exact content digest of one declared local input."""
    return hashlib.sha256(read_regular_path(path)).hexdigest()


def _validate_arguments(args: object, output: Path) -> _Request:
    """Reject every unusable option before any cohort byte is read.

    Args:
        args: Parsed namespace carrying manifest/captures Paths, the declared objective
            and its optional rate floors, the caller selection, the repeated axis names,
            and the fold/seed/worker/breakpoint knobs.
        output: Empty staged private directory supplied by the CLI atomic adapter.

    Returns:
        The validated request.

    Raises:
        ValueError: For a malformed path, an unsupported caller, axis or objective, an axis
            of a caller the selection does not search, an adVNTR run without an
            executable, or a non-empty staging directory.
    """
    manifest = getattr(args, "manifest", None)
    captures = getattr(args, "captures", None)
    executable = getattr(args, "advntr_executable", None)
    if not isinstance(manifest, Path) or not isinstance(captures, Path):
        _fail("cutoff optimize manifest and captures must be Paths")
    if executable is not None and not isinstance(executable, Path):
        _fail("cutoff optimize manifest and captures must be Paths")
    caller = getattr(args, "caller", "kestrel")
    if caller not in _CALLERS:
        _fail("cutoff optimize caller must be kestrel, advntr, or both")
    requested = getattr(args, "axes", None)
    if requested is None:
        # ``--axis`` is repeatable, so argparse leaves it unset rather than defaulted.
        requested = list(_DEFAULT_AXES[caller])
    if isinstance(requested, str) or not isinstance(requested, Sequence) or not requested:
        _fail("cutoff optimize requires at least one axis")
    axes = tuple(requested)
    unsupported = [name for name in axes if name not in _ALL_AXES]
    if unsupported or len(set(axes)) != len(axes):
        _fail(f"cutoff axis name must be one of {sorted(_ALL_AXES)}, each declared at most once")
    advntr_axes = sorted(name for name in axes if axis_caller(name) == "advntr")
    kestrel_axes = sorted(name for name in axes if axis_caller(name) == "kestrel")
    if caller == "kestrel" and advntr_axes:
        _fail(
            f"cutoff optimize --caller kestrel cannot search the adVNTR axis {advntr_axes}; use --caller advntr or both"
        )
    if caller == "advntr" and kestrel_axes:
        _fail(
            f"cutoff optimize --caller advntr cannot search the Kestrel axis {kestrel_axes}; use --caller kestrel or both"
        )
    objective = getattr(args, "objective", None)
    if not isinstance(objective, str):
        _fail("unsupported cutoff search objective")
    spec = SearchSpec(
        objective=objective,
        min_sensitivity=getattr(args, "min_sensitivity", None),
        min_specificity=getattr(args, "min_specificity", None),
    )
    workers = getattr(args, "workers", 1)
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        _fail("cutoff optimize workers must be a positive integer")
    maximum = getattr(args, "max_breakpoints", None)
    if maximum is not None and (type(maximum) is not int or maximum < 3):
        _fail("cutoff optimize max_breakpoints must be an integer of at least 3")
    folds = getattr(args, "folds", 5)
    seed = getattr(args, "seed", 20260915)
    group_folds({}, folds=folds, seed=seed)
    if caller in {"advntr", "both"} and executable is None:
        _fail("cutoff optimize with adVNTR requires --advntr-executable")
    if not output.is_dir() or output.is_symlink() or any(output.iterdir()):
        _fail("cutoff optimize output must be an empty staged directory")
    return _Request(manifest, captures, spec, caller, axes, maximum, folds, seed, workers, executable)


def _decode_capture(path: Path, key: str) -> KestrelCapture:
    """Decode one complete capture, naming the sample whose evidence is unusable."""
    try:
        raw = read_regular_path(path)
    except (OSError, ValueError):
        _fail(f"cutoff optimize capture for {key} is missing, unsafe, or unreadable")
    try:
        return decode_kestrel_capture(load_strict_json_object(raw))
    except ValueError:
        _fail(f"cutoff optimize capture for {key} is not strict JSON or not a complete capture")


def _capture_baseline(captures: Mapping[str, KestrelCapture]) -> CallerPolicyValues:
    """The one policy every capture was produced under; disagreement is a defect."""
    policies = {capture.baseline_policy.sha256: capture.baseline_policy for capture in captures.values()}
    if len(policies) != 1:
        _fail("cutoff optimize captures declare more than one baseline policy")
    return next(iter(policies.values()))


def _permissive_axis(name: str, baseline: CallerPolicyValues) -> tuple[CutoffCandidate, str]:
    """Discover the most permissive admissible projection of one axis, and its statistic.

    Args:
        name: One of the supported axis names.
        baseline: The shipped complete policy the projection is built on.

    Returns:
        The complete permissive policy as a candidate, and the column the axis reads.

    Raises:
        ValueError: If the policy decoder refuses every probed extreme of the axis.
    """
    ladder = AXIS_PROBE[name][1]
    probe = declared_axis(name, list(ladder), baseline=baseline)
    for value, reason in probe.rejected:
        logger.info("cutoff axis %s probe rung %r refused: %s", name, value, reason)
    accepted = [value for value in ladder if value in probe.values]
    if not accepted:
        _fail(f"cutoff optimize found no admissible permissive projection for axis {name}")
    candidates = axis_candidates(baseline, probe)
    chosen = next(candidate for candidate, value in zip(candidates, probe.values, strict=True) if value == accepted[0])
    return chosen, str(axis_document(probe)["statistic"])


def _observed_breakpoints(
    captures: Mapping[str, KestrelCapture], policy: CallerPolicyValues, statistic: str
) -> dict[str, tuple[Fraction, ...]]:
    """Replay every capture once at the permissive projection and collect its values."""
    observed: dict[str, tuple[Fraction, ...]] = {}
    for key in sorted(captures):
        capture = captures[key]
        result = replay_kestrel_capture(capture, policy, capture_policy_sha256=capture.provenance.capture_policy_sha256)
        frame = kestrel_replay_prefilter_frame(result)
        # An empty capture has no candidate rows at all, so it contributes no breakpoint;
        # its prefilter frame does not even carry the gate columns to be filtered on.
        observed[key] = () if frame.empty else eligible_statistic_values(frame, statistic)
    return observed


def _derive_axes(
    request: _Request,
    baseline: CallerPolicyValues,
    captures: Mapping[str, KestrelCapture],
    assignments: Mapping[str, int],
) -> tuple[tuple[DerivedAxis, ...], dict[str, Mapping[int, frozenset[int | float]]]]:
    """Turn each requested Kestrel axis into observed breakpoints and complete replayable policies.

    Returns:
        The derived Kestrel axes, each merged with its fold-local values, and per axis the
        values each outer fold's training samples derived.
    """
    derived: list[DerivedAxis] = []
    fold_values: dict[str, Mapping[int, frozenset[int | float]]] = {}
    for name in request.axes:
        if axis_caller(name) != "kestrel":
            continue
        permissive, statistic = _permissive_axis(name, baseline)
        observed = _observed_breakpoints(captures, permissive.policy, statistic)

        def derive(values: Mapping[str, tuple[Fraction, ...]], axis: str = name) -> AxisBreakpoints:
            return derive_axis(axis, values, baseline=baseline, max_values=request.max_breakpoints)

        axis, fold_values[name] = fold_axis(derive, observed, assignments)
        derived.append((axis, axis_candidates(baseline, axis)))
    return tuple(derived), fold_values


def _require_research_capture(advntr_paths: Mapping[str, Path], baseline: CallerPolicyValues) -> None:
    """Refuse adVNTR evidence whose capture settings an exported research profile cannot reproduce.

    The research runtime renders adVNTR's capture settings from fixed constants
    (``RESEARCH_CAPTURE_PARAMETERS``) and only the caller values from the profile, so a
    cutoff derived under any other capture setting would be applied under different semantics.
    """
    differing = research_capture_differences(read_capture_policy(advntr_paths), baseline)
    if differing:
        _fail(
            "cutoff optimize adVNTR captures were produced under capture settings the research runtime cannot "
            f"reproduce: {', '.join(differing)}; an exported profile would run adVNTR under different semantics"
        )


def _anchor(candidates: Sequence[CutoffCandidate]) -> CutoffCandidate:
    """The candidate that reproduces the shipped policy exactly."""
    anchors = [candidate for candidate in candidates if not candidate.parameters]
    if len(anchors) != 1:
        _fail("cutoff optimize axis must contain exactly one candidate reproducing the baseline")
    return anchors[0]


def _endpoint(replay: KestrelGridReplay, policy_id: str, key: str) -> tuple[object, ...]:
    """The comparable part of one replayed endpoint: what the caller actually decided."""
    row = replay.observations[policy_id][key]
    return (row.disposition, row.called_positive, row.confidence, row.flag, row.canonical_identity)


def _prove_baseline_parity(replay: KestrelGridReplay, anchors: Mapping[str, str]) -> dict[str, Any]:
    """Require the anchor candidates to reproduce the capture baseline sample by sample.

    Args:
        replay: The policy-major grid replay, which already refused any native Kestrel
            result that disagrees with its own capture replay.
        anchors: Axis name to the candidate id that reproduces the shipped policy.

    Returns:
        The parity record published in the report.

    Raises:
        ValueError: If any anchor disagrees with the baseline arm for any sample.
    """
    mismatches = [
        (axis, anchor, key)
        for axis, anchor in sorted(anchors.items())
        for key in replay.sample_keys
        if _endpoint(replay, BASELINE_ID, key) != _endpoint(replay, anchor, key)
    ]
    if mismatches:
        _fail(
            "cutoff optimize baseline parity failed: the candidate reproducing the shipped values disagrees "
            f"with the captured baseline for {len(mismatches)} sample/axis pairs"
        )
    modes = list(replay.baseline_parity.values())
    return {
        "proven": True,
        "anchor_candidate_ids": dict(sorted(anchors.items())),
        "native_exact_count": modes.count("native-exact"),
        "capture_replay_authoritative_count": modes.count("capture-replay-authoritative"),
        "mismatches": [],
    }


def _project_policy(components: Mapping[str, object], policy: CallerPolicyValues) -> CallerPolicyValues:
    """Read the resolved runtime components back into a complete caller policy."""
    return project_caller_policy(components, tuple(policy.values), policy.required_callers)


def _export_profile(output: Path, request: _Request, selected: CallerPolicyValues) -> dict[str, Any]:
    """Write the selected policy as a research profile and prove it round-trips.

    Args:
        output: Staged private output directory.
        request: The validated request, supplying the seed and the input digests.
        selected: The policy the objective selected.

    Returns:
        The profile record published in the report, including the round-trip boolean.

    Raises:
        ValueError: If the selected policy carries a non-legacy adVNTR mode (refused before
            anything is written), if the written profile does not resolve back to the
            selected policy, or if its adVNTR arguments cannot be rendered.
    """
    mode = selected.values.get(_ADVNTR_MODE) if "advntr" in selected.required_callers else "legacy"
    if mode != "legacy":
        _fail(f"cutoff optimize cannot export a research profile with adVNTR mode {mode!r}: {RESEARCH_LEGACY_ONLY}")
    profile = build_caller_generated_profile(
        selected,
        dataset_manifest_hash=_digest(request.manifest),
        partition_manifest_hash=_digest(request.captures),
        seed=request.seed,
        generator_version=GENERATOR_VERSION,
    )
    path = output / PROFILE_NAME
    _write_private(path, profile.canonical_bytes)
    resolved = resolve_research_decision_profile(path)
    if _project_policy(resolved.components, selected) != selected:
        _fail("cutoff optimize research profile did not round-trip to the selected policy")
    if "advntr" in selected.required_callers:
        # Render the profile exactly as the research runtime will; a nominal thread count suffices.
        research_policy_argv(
            cast(Mapping[str, object], resolved.components["advntr"]),
            cast(Mapping[str, object], resolved.components["kestrel"]),
            1,
        )
    return {
        "status": "available",
        "path": PROFILE_NAME,
        "profile_id": resolved.profile_id,
        "sha256": hashlib.sha256(profile.canonical_bytes).hexdigest(),
        "round_trip_sha256": resolved.digest,
        "round_trip_matches_selected_policy": True,
    }


def _advntr_arms(
    request: _Request,
    output: Path,
    advntr_paths: Mapping[str, Path],
    policies: Mapping[str, CallerPolicyValues],
    anchor_id: str,
    primary: Sequence[CohortSample],
    search: AdvntrAxisSearch | None,
    snapshot: Mapping[str, str],
) -> tuple[dict[str, tuple[CallerObservation, ...]], AdvntrCutoffGridResult, dict[str, Any]]:
    """Replay every candidate natively and prove the replay is the one the axes describe.

    Returns:
        The adVNTR arms (with the baseline arm), the candidate grid, and the adVNTR parity,
        replay-consistency and main-grid wall-time records for the report.
    """
    assert request.advntr_executable is not None
    started = time.monotonic()
    result = evaluate_advntr_cutoff_grid(
        advntr_paths,
        policies,
        baseline_policy_id=anchor_id,
        executable_path=request.advntr_executable,
        output=output / "advntr",
    )
    seconds = time.monotonic() - started
    require_one_execution_per_signature(result, policies)
    checked = 0
    if search is not None:
        require_same_evidence(search.probe, result, snapshot, advntr_paths)
        for axis, candidates in search.derived:
            checked += check_replay_consistency(result, candidates, search.visits[axis.axis])
    parity = advntr_baseline_parity(result, anchor_id, advntr_paths)
    arms = advntr_observation_arms(result, primary)
    arms[BASELINE_ID] = arms[anchor_id]
    consistency = None if search is None else {"checked_candidates": checked, "mismatches": []}
    return arms, result, {"advntr_parity": parity, "replay_consistency": consistency, "main_seconds": seconds}


def run_cutoff_optimization(args: object, output: Path) -> bool:
    """Derive, score, select and export caller cutoffs from labelled cohort evidence.

    Args:
        args: Parsed namespace; see :func:`_validate_arguments` for the fields read.
        output: Empty private staging directory owned by the CLI atomic adapter.

    Returns:
        Whether the declared objective selected a cutoff. An objective whose constraints
        no candidate satisfies still writes the complete report and returns ``False``.

    Raises:
        ValueError: For malformed arguments or evidence, for captures that disagree about
            the baseline policy, for a Kestrel or adVNTR baseline parity failure, for adVNTR
            probe and candidate replays that do not share their evidence, disagree with the
            legacy rule, or execute an unexpected number of adVNTR policies, or for a
            research profile that does not resolve back to the selected policy.
    """
    request = _validate_arguments(args, output)
    samples = read_cohort_manifest(request.manifest)
    primary = primary_samples(samples)
    declared = read_cutoff_captures(
        request.captures, samples, caller="kestrel" if request.caller == "kestrel" else "both"
    )
    keys = tuple(sample.sample_id for sample in primary)
    if not set(keys) <= set(declared.kestrel):
        _fail("cutoff optimize capture manifest must declare every primary cohort sample")
    capture_paths = {key: declared.kestrel[key] for key in keys}
    native_paths = {key: declared.native_kestrel[key] for key in keys}
    captures = {key: _decode_capture(path, key) for key, path in sorted(capture_paths.items())}
    baseline = _capture_baseline(captures)
    # The same grouped, truth-stratified folds evaluate_cutoff_arms computes; every axis is
    # derived per fold from training samples before anything is scored.
    assignments = outer_fold_assignments(
        [(sample.sample_id, sample.group_id, sample.genotype) for sample in primary],
        folds=request.folds,
        seed=request.seed,
    )
    derived, fold_values = _derive_axes(request, baseline, captures, assignments)
    advntr_names = tuple(name for name in request.axes if axis_caller(name) == "advntr")
    advntr_paths = {key: declared.advntr[key] for key in keys} if request.caller != "kestrel" else {}
    if advntr_paths:
        _require_research_capture(advntr_paths, baseline)
    snapshot = capture_snapshot(advntr_paths) if advntr_names else {}
    search: AdvntrAxisSearch | None = None
    if advntr_names:
        assert request.advntr_executable is not None
        search = derive_advntr_axes(
            advntr_names,
            baseline,
            advntr_paths,
            executable_path=request.advntr_executable,
            output=output / "advntr-probe",
            max_values=request.max_breakpoints,
            assignments=assignments,
        )
        derived = (*derived, *search.derived)
        fold_values = {**fold_values, **search.fold_values}
    order = {name: index for index, name in enumerate(request.axes)}
    derived = tuple(sorted(derived, key=lambda entry: order[entry[0].axis]))
    inventories = fold_candidate_inventories(derived, fold_values)
    policies = {candidate.candidate_id: candidate.policy for _, candidates in derived for candidate in candidates}
    anchors = {axis.axis: _anchor(candidates) for axis, candidates in derived}
    try:
        replay = replay_kestrel_grid(capture_paths, policies, workers=request.workers, native_paths=native_paths)
    except ValueError as error:
        _fail(f"cutoff optimize Kestrel grid replay failed, so baseline parity is unproven: {error}")
    parity = _prove_baseline_parity(replay, {axis: anchor.candidate_id for axis, anchor in anchors.items()})
    arms: dict[str, tuple[CallerObservation, ...]] = kestrel_observation_arms(replay, primary)
    advntr_result: AdvntrCutoffGridResult | None = None
    checks: dict[str, Any] = {"advntr_parity": None, "replay_consistency": None, "main_seconds": None}
    if request.caller != "kestrel":
        # Under --caller advntr the Kestrel arms are replayed but not scored: they still
        # bind the baseline through Kestrel parity and provenance.
        advntr, advntr_result, checks = _advntr_arms(
            request, output, advntr_paths, policies, anchors[request.axes[0]].candidate_id, primary, search, snapshot
        )
        arms = advntr if request.caller == "advntr" else combine_arms(arms, advntr)
    if {row.key for row in arms[BASELINE_ID]} != set(captures):
        _fail("cutoff optimize arm rows are not keyed by the capture sample identifiers")
    evaluation = evaluate_cutoff_arms(
        arms,
        baseline_id=BASELINE_ID,
        spec=request.spec,
        folds=request.folds,
        seed=request.seed,
        fold_inventories=inventories,
        # Candidate IDs number the merged inventory, which held-out values help build, so a
        # tie is broken by policy content instead (Codex final review, X1).
        tie_keys={policy_id: policy.sha256 for policy_id, policy in policies.items()},
    )
    # The inventories were built on the folds computed above; the evaluation's own folds must
    # be the same allocation, not just the same fold numbers, or admissibility would misfire.
    if {row["key"]: row["fold"] for row in evaluation["rows"] if row["fold"] is not None} != assignments:
        _fail("cutoff optimize fold inventories were derived on different outer folds than the evaluation")
    document = build_cutoff_report_document(
        CutoffReportInputs(
            objective=request.spec,
            caller=request.caller,
            folds=request.folds,
            seed=request.seed,
            workers=request.workers,
            max_breakpoints=request.max_breakpoints,
            samples=samples,
            primary=primary,
            derived=derived,
            anchors=anchors,
            arms=arms,
            replay=replay,
            parity=parity,
            evaluation=evaluation,
            advntr_result=advntr_result,
            comparisons=AXIS_COMPARISON,
            cohort_manifest_sha256=_digest(request.manifest),
            capture_manifest_sha256=_digest(request.captures),
            generator_version=GENERATOR_VERSION,
            advntr_parity=checks["advntr_parity"],
            replay_consistency=checks["replay_consistency"],
            advntr_search=search,
            advntr_main_seconds=checks["main_seconds"],
        )
    )
    selected = evaluation["final_selection"]["policy_id"]
    if selected is not None:
        document["profile"] = _export_profile(output, request, policies.get(selected, baseline))
    write_cutoff_reports(output, document)
    return selected is not None


__all__ = [
    "AXIS_COMPARISON",
    "AXIS_PROBE",
    "BASELINE_ID",
    "GENERATOR_VERSION",
    "PROFILE_NAME",
    "USAGE_HINT",
    "run_cutoff_optimization",
]
