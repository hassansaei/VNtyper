"""Static caller comparison reports from complete frozen populations."""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

from jinja2 import Environment, FileSystemLoader, StrictUndefined, select_autoescape

import vntyper
from vntyper.scripts.calibration_caller_curves import (
    CallerCurvePoint,
    CallerCurves,
    CallerOperatingPoint,
    build_caller_curves,
)
from vntyper.scripts.calibration_caller_metrics import CallerMetrics, CallerObservation, calculate_caller_metrics
from vntyper.scripts.calibration_caller_roster import (
    CallerEligibleRoster,
    bind_caller_observations,
    caller_eligible_roster_document,
)

logger = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_PHASES = {"policy-selection", "validation", "locked-heldout", "development-assessment"}


@dataclass(frozen=True)
class CallerReportCandidate:
    """Fixed production outcomes and reasons supplied by the acceptance engine."""

    candidate_id: str
    observations: tuple[CallerObservation, ...]
    rejection_reasons: tuple[str, ...]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object) -> None:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        _fail("caller report identities require lowercase SHA256 digests")


def _strings(values: object) -> None:
    if (
        not isinstance(values, tuple)
        or any(not isinstance(value, str) or not value or value.strip() != value for value in values)
        or len(values) != len(set(values))
    ):
        _fail("caller report labels must be unique non-empty strings in an immutable tuple")


def _metric_rows(metrics: CallerMetrics) -> list[tuple[str, str]]:
    result = [
        ("Eligible independent groups", str(metrics.eligible_count)),
        ("Known / unknown truth", f"{metrics.known_truth_count} / {metrics.unknown_truth_count}"),
        (
            "TP / FP / TN / FN",
            f"{metrics.true_positives} / {metrics.false_positives} / {metrics.true_negatives} / {metrics.false_negatives}",
        ),
        ("Positive / negative truth no-calls", f"{metrics.positive_no_calls} / {metrics.negative_no_calls}"),
    ]
    for label, rate in (
        ("Sensitivity", metrics.sensitivity),
        ("Specificity", metrics.specificity),
        ("False-positive rate", metrics.false_positive_rate),
        ("Precision (study prevalence)", metrics.precision),
        ("No-call rate", metrics.no_call_rate),
        ("Assessability", metrics.assessability),
        ("Exact variant-set recovery", metrics.exact_variant_recovery),
    ):
        value = (
            "undefined (0/0; empty metric denominator)"
            if rate.estimate is None or rate.lower is None or rate.upper is None
            else (
                f"{rate.events}/{rate.total} ({100 * float(rate.estimate):.2f}%; "
                f"central 95% CI {100 * float(rate.lower):.2f}–{100 * float(rate.upper):.2f}%)"
            )
        )
        result.append((label, value))
    upper = metrics.fpr_one_sided_upper
    result.extend(
        (
            ("FPR one-sided 95% upper bound", "undefined" if upper is None else f"{100 * float(upper):.4f}%"),
            (
                "Study prevalence among known truth",
                "undefined" if metrics.study_prevalence is None else str(metrics.study_prevalence),
            ),
            ("Positive truth missing variant identity", str(metrics.positive_truth_missing_identity)),
            ("Identity-assessable groups", str(metrics.identity_assessable_count)),
            ("Wrong identity groups", str(metrics.wrong_identity_groups)),
            ("Wrong tier-A identity groups", str(metrics.wrong_tier_a_identity_groups)),
        )
    )
    return result


def _plots(
    curves: tuple[CallerCurves, ...], candidates: dict[str, CallerReportCandidate], phase: str
) -> list[dict[str, object]]:
    if not isinstance(curves, tuple):
        _fail("caller report curves must be an immutable tuple")
    if curves and phase != "policy-selection":
        _fail("cutoff-search curves are only permitted on policy-selection evidence")
    result: list[dict[str, object]] = []
    for curve in curves:
        if not isinstance(curve, CallerCurves) or not isinstance(curve.points, tuple) or not curve.points:
            _fail("caller report requires typed curve points")
        if any(not isinstance(point, CallerCurvePoint) for point in curve.points):
            _fail("caller report requires typed curve points")
        if any(point.candidate_id not in candidates for point in curve.points):
            _fail("caller report curve references an absent candidate")
        rebuilt = build_caller_curves(
            tuple(
                CallerOperatingPoint(
                    point.candidate_id,
                    point.threshold,
                    curve.fixed_policy_sha256,
                    candidates[point.candidate_id].observations,
                )
                for point in curve.points
            ),
            comparison=curve.comparison,
            phase=phase,
        )
        if rebuilt != curve:
            _fail("caller report curve differs from its measured candidate outcomes")
        for kind in ("ROC", "Precision–recall"):
            points = []
            for point in curve.points:
                x_value = point.false_positive_rate if kind == "ROC" else point.sensitivity
                y_value = point.sensitivity if kind == "ROC" else point.precision
                if y_value is None:
                    continue
                points.append(
                    {
                        "x": 40 + 320 * float(x_value),
                        "y": 340 - 320 * float(y_value),
                        "label": f"cutoff {curve.comparison} {point.threshold}; candidate {point.candidate_id}",
                    }
                )
            result.append(
                {
                    "kind": kind,
                    "policy": curve.fixed_policy_sha256,
                    "points": points,
                    "line": " ".join(f"{point['x']},{point['y']}" for point in points),
                    "no_calls": curve.no_call_count,
                    "eligible": curve.eligible_count,
                    "x_label": "False-positive rate" if kind == "ROC" else "Recall",
                    "y_label": "Sensitivity" if kind == "ROC" else "Precision",
                }
            )
    return result


def render_caller_comparison(
    candidates: Sequence[CallerReportCandidate],
    roster: CallerEligibleRoster,
    *,
    phase: str,
    baseline_id: str,
    selected_id: str | None,
    protocol_sha256: str,
    evidence_sha256: str,
    required_strata: tuple[str, ...],
    curves: tuple[CallerCurves, ...] = (),
) -> str:
    """Render escaped offline HTML, recomputing all rates from complete outcomes.

    Args:
        candidates: Fixed baseline and candidate observations with gate reasons.
        roster: Population and stratum memberships frozen before performance.
        phase: Policy-selection, validation, locked-heldout or development-assessment.
        baseline_id: Exact baseline candidate identity.
        selected_id: Selected/frozen candidate, or None after failed selection.
        protocol_sha256: Frozen protocol identity, supplied by the study engine.
        evidence_sha256: Authorized phase evidence identity.
        required_strata: Every predeclared mandatory population, including sparse ones.
        curves: Verified comparable scalar-family operating curves, selection only.

    Returns:
        Standalone HTML with counts, intervals, mandatory strata, rejection reasons
        and optional SVG plots. Sample keys and raw reads are never rendered.
        Rendering is descriptive; it does not authorize exposure or promotion.

    Raises:
        ValueError: If rosters, truth, phases, identities or curves are inconsistent.
    """
    if not isinstance(phase, str) or phase not in _PHASES:
        _fail("unsupported caller report phase")
    for value in (baseline_id, protocol_sha256, evidence_sha256):
        _digest(value)
    if selected_id is not None:
        _digest(selected_id)
    _strings(required_strata)
    if not required_strata:
        _fail("caller reports require the frozen mandatory strata")
    caller_eligible_roster_document(roster)
    if not isinstance(candidates, (tuple, list)) or not candidates:
        _fail("caller reports require non-empty candidate observations")
    by_id: dict[str, CallerReportCandidate] = {}
    truth: dict[str, tuple[bool | None, tuple[str, ...] | None]] | None = None
    entries: list[dict[str, object]] = []
    for candidate in candidates:
        if not isinstance(candidate, CallerReportCandidate) or not isinstance(candidate.observations, tuple):
            _fail("caller reports require immutable CallerReportCandidate values")
        _digest(candidate.candidate_id)
        _strings(candidate.rejection_reasons)
        if candidate.candidate_id in by_id:
            _fail("caller report candidate identities must be unique")
        rows = bind_caller_observations(candidate.observations, roster)
        observed_truth = {row.key: (row.truth_positive, row.truth_variants) for row in rows}
        if truth is not None and truth != observed_truth:
            _fail("caller report baseline and candidate truth differ")
        truth = observed_truth
        by_id[candidate.candidate_id] = candidate
        populations: list[dict[str, object]] = [
            {"name": "Pooled", "rows": _metric_rows(calculate_caller_metrics(rows))}
        ]
        for name in required_strata:
            keys = {member.key for member in roster.members if name in member.strata}
            subset = tuple(row for row in rows if row.key in keys)
            populations.append(
                {"name": name, "rows": _metric_rows(calculate_caller_metrics(subset)) if subset else None}
            )
        label = (
            "Baseline"
            if candidate.candidate_id == baseline_id
            else (
                ("Selected candidate" if phase == "policy-selection" else "Frozen candidate")
                if candidate.candidate_id == selected_id
                else "Rejected"
                if candidate.rejection_reasons
                else "Unselected"
            )
        )
        entries.append(
            {
                "id": candidate.candidate_id,
                "label": label,
                "reasons": candidate.rejection_reasons,
                "populations": populations,
            }
        )
    if baseline_id not in by_id or (selected_id is not None and selected_id not in by_id):
        _fail("caller report baseline and selected identities must name supplied candidates")
    if selected_id == baseline_id:
        _fail("caller report selected candidate must differ from baseline")
    if phase != "policy-selection":
        if selected_id is None:
            _fail("non-selection caller reports require a frozen candidate")
        if set(by_id) != {baseline_id, selected_id}:
            _fail("non-selection reports require exactly the frozen baseline and candidate")
    elif selected_id is not None and by_id[selected_id].rejection_reasons:
        _fail("a rejected policy-selection candidate cannot be selected")
    environment = Environment(
        loader=FileSystemLoader(str(Path(vntyper.__file__).resolve().parent / "templates")),
        autoescape=select_autoescape(["html"]),
        undefined=StrictUndefined,
    )
    return environment.get_template("calibration_caller_report.html").render(
        phase=phase,
        selected_id=selected_id,
        protocol=protocol_sha256,
        evidence=evidence_sha256,
        roster=roster.sha256,
        entries=sorted(entries, key=lambda entry: str(entry["id"])),
        plots=_plots(curves, by_id, phase),
    )
