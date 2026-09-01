"""Derived static tables and HTML for calibration artifact directories."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path

from vntyper.scripts.calibration_artifact_io import write_json
from vntyper.scripts.calibration_objective import CandidateEvaluation
from vntyper.scripts.calibration_report import decode_calibration_report, write_calibration_report
from vntyper.scripts.calibration_scalar_replay import replay_scalar_dominance
from vntyper.scripts.calibration_statistics import clopper_pearson_interval, deterministic_curves, joint_surface
from vntyper.scripts.calibration_workflow import ExtractedEvidence
from vntyper.scripts.decision_profile import ResolvedDecisionProfile
from vntyper.version import __version__


def write_evaluation_artifacts(
    output: Path,
    *,
    phase: str,
    profile: ResolvedDecisionProfile,
    evidence: ExtractedEvidence,
    evaluation: CandidateEvaluation,
    accessed_roles: tuple[str, ...],
) -> None:
    """Write complete deterministic metrics tables and static HTML."""
    if phase not in {"fitted", "validation", "held-out"}:
        raise ValueError(f"unsupported calibration reporting phase: {phase}")
    if not isinstance(profile, ResolvedDecisionProfile) or not isinstance(evidence, ExtractedEvidence):
        raise ValueError("calibration reporting requires resolved profile and evidence")
    if not isinstance(evaluation, CandidateEvaluation):
        raise ValueError("calibration reporting requires CandidateEvaluation")
    analysis = _analyse(profile, evidence)
    protocol = evidence.study.protocol
    intervals = _interval_rows(evaluation, analysis["mutated"], analysis["detected"])
    curves = _curve_rows(analysis["scores"])
    surface = joint_surface(protocol.candidate_grid)
    surface_rows = tuple(" ".join(f"{name}={_text(value)}" for name, value in row) for row in surface)
    report_raw = {
        "schema_version": "calibration-report-v1",
        "phase": phase,
        "profile_sha256": profile.digest,
        "protocol_sha256": protocol.sha256,
        "evidence_sha256": evidence.dataset_sha256,
        "objective": protocol.objective,
        "tier_metrics": analysis["tier_metrics"],
        "abstentions": analysis["abstentions"],
        "provenance": _provenance(evidence, accessed_roles),
        "statistics": {
            "intervals": list(intervals),
            "roc_rows": [row["roc"] for row in curves],
            "pr_rows": [row["pr"] for row in curves],
            "joint_surface_rows": list(surface_rows),
        },
        "limitations": [
            "Local append-only custody guards are not proof of independent custody; closure requires a named external custodian.",
            "Small-n strata and absent boundary cells widen uncertainty and remain visible rather than pooled away.",
            "Reporting an interval is not a clinical safety claim.",
        ],
    }
    report = decode_calibration_report(report_raw)
    write_json(output / "report_data.json", report_raw)
    write_calibration_report(output / "report.html", report)
    write_json(
        output / "grid.json",
        {
            "schema_version": "calibration-grid-v1",
            "objective": protocol.objective,
            "candidate_grid": {
                name: [_text(value) for value in values] for name, values in protocol.candidate_grid.items()
            },
        },
    )
    write_json(
        output / "intervals.json",
        {"schema_version": "calibration-intervals-v1", "rows": list(intervals)},
    )
    _write_lines(output / "roc.tsv", "roc_row", (row["roc"] for row in curves))
    _write_lines(output / "pr.tsv", "pr_row", (row["pr"] for row in curves))
    _write_lines(output / "joint_surface.tsv", "candidate", surface_rows)
    _write_abstentions(output / "abstentions.tsv", analysis["abstentions"])


def _analyse(profile: ResolvedDecisionProfile, evidence: ExtractedEvidence) -> dict[str, object]:
    component = profile.components["dominance"]
    if not isinstance(component, Mapping):
        raise ValueError("calibration reporting dominance component must be an object")
    labels = {row.manifest_key: row for row in evidence.labels.rows}
    baseline = evidence.baseline["expected"]
    if not isinstance(baseline, Mapping):
        raise ValueError("calibration reporting baseline expected projection must be an object")
    rows = baseline["rows"]
    if not isinstance(rows, tuple) or len(rows) != len(evidence.features.rows):
        raise ValueError("calibration reporting rows must align with features")
    tiers: dict[str, dict[str, object]] = {
        tier: {"tier": tier, "displayed": 0, "exact": 0, "wrong": 0} for tier in ("A", "B", "C", "D")
    }
    reasons: Counter[str] = Counter()
    scores: list[tuple[Fraction, bool]] = []
    mutated = detected = applicable = 0
    for feature, baseline_row in zip(evidence.features.rows, rows, strict=True):
        label = labels[feature.manifest_key]
        if not isinstance(baseline_row, Mapping):
            raise ValueError("calibration reporting baseline row must be an object")
        decision = replay_scalar_dominance(feature.features, component)
        applicable += decision.applicable
        if label.truth_status == "mutated":
            mutated += 1
            detected += decision.selected_identity is not None
        if decision.abstention_reason is not None:
            reasons[decision.abstention_reason] += 1
        tier = baseline_row.get("tier")
        if decision.selected_identity is not None and isinstance(tier, str):
            counts = tiers.setdefault(tier, {"tier": tier, "displayed": 0, "exact": 0, "wrong": 0})
            _increment(counts, "displayed")
            is_exact = (
                decision.selected_identity == label.expected_identity
                and baseline_row.get("name") == label.expected_display_name
            )
            _increment(counts, "exact" if is_exact else "wrong")
        score_value = feature.features.get("haplotype_record_share", 0)
        if isinstance(score_value, bool) or not isinstance(score_value, (int, float)):
            raise ValueError("calibration reporting record share must be numeric")
        scores.append((Fraction(str(score_value)), label.truth_status == "mutated"))
    split = "+".join(sorted({member.role for member in evidence.study.partitions.members if member.key in labels}))
    abstentions = [
        {"split": split, "reason": reason, "count": count, "rate": str(Fraction(count, applicable))}
        for reason, count in sorted(reasons.items())
    ]
    return {
        "tier_metrics": list(tiers.values()),
        "abstentions": abstentions,
        "scores": tuple(scores),
        "mutated": mutated,
        "detected": detected,
    }


def _interval_rows(evaluation: CandidateEvaluation, mutated: object, detected: object) -> tuple[str, ...]:
    if isinstance(mutated, bool) or not isinstance(mutated, int) or mutated <= 0:
        raise ValueError("calibration reporting requires mutated members")
    if isinstance(detected, bool) or not isinstance(detected, int):
        raise ValueError("calibration reporting detected count must be an integer")
    exact = clopper_pearson_interval(detected, mutated)
    return (
        f"two-sided 95% binary-detection interval {exact.lower}–{exact.upper}",
        f"one-sided paired detection lower bound {evaluation.detection_lower_bound}",
        f"one-sided paired macro exact-recovery lower bound {evaluation.macro_exact_lower_bound}",
        f"Holm-adjusted family p-value {evaluation.holm_adjusted_p_value}",
    )


def _curve_rows(scores: object) -> tuple[dict[str, str], ...]:
    if not isinstance(scores, tuple):
        raise ValueError("calibration reporting scores must be a tuple")
    labels = {label for _, label in scores}
    if labels != {False, True}:
        return ()
    return tuple(
        {
            "roc": (
                f"threshold={point.threshold} TP={point.true_positives} FP={point.false_positives} "
                f"TN={point.true_negatives} FN={point.false_negatives} TPR={point.true_positive_rate} "
                f"FPR={point.false_positive_rate}"
            ),
            "pr": f"threshold={point.threshold} precision={point.precision} recall={point.recall}",
        }
        for point in deterministic_curves(scores)
    )


def _provenance(evidence: ExtractedEvidence, accessed_roles: tuple[str, ...]) -> dict[str, list[str]]:
    features = tuple(row.features for row in evidence.features.rows)
    labels = evidence.labels.rows

    def values(name: str) -> list[str]:
        return sorted({str(row[name]) for row in features if row.get(name) not in {None, ""}})

    assays = values("assay_class")
    depths = values("active_region_depth")
    mutation_classes = sorted({row.mutation_class for row in labels})
    mutated = sum(row.truth_status == "mutated" for row in labels)
    controls = len(labels) - mutated
    boundary = [
        f"{name}: {', '.join(_text(value) for value in values)}"
        for name, values in evidence.study.protocol.candidate_grid.items()
    ]
    return {
        "software_versions": values("tool_version") or [f"VNtyper {__version__}"],
        "reference_versions": values("reference_version") or ["not recorded in calibration evidence"],
        "sample_composition": [f"total={len(labels)}", f"mutated={mutated}", f"controls={controls}"],
        "assays": assays or ["not recorded in calibration evidence"],
        "depths": depths or ["not recorded in calibration evidence"],
        "read_lengths": ["not recorded in calibration evidence"],
        "independent_array_size": ["not recorded independently in this evidence bundle"],
        "mutation_classes": mutation_classes,
        "manifest_hashes": [
            f"study={evidence.study_sha256}",
            f"partition={evidence.study.partitions.sha256}",
            f"dataset={evidence.dataset_sha256}",
        ],
        "access_attempts": [f"{role}: one authorized operation access" for role in accessed_roles],
        "boundary_coverage": boundary,
        "seeds": [f"bootstrap and candidate-family seed={evidence.study.protocol.seed}"],
    }


def _write_lines(path: Path, heading: str, rows) -> None:
    path.write_text(heading + "\n" + "".join(f"{row}\n" for row in rows), encoding="utf-8")


def _write_abstentions(path: Path, rows: object) -> None:
    if not isinstance(rows, list):
        raise ValueError("calibration reporting abstentions must be a list")
    text = "split\treason\tcount\trate\n" + "".join(
        f"{row['split']}\t{row['reason']}\t{row['count']}\t{row['rate']}\n" for row in rows
    )
    path.write_text(text, encoding="utf-8")


def _text(value: object) -> str:
    return str(value)


def _increment(counts: dict[str, object], field: str) -> None:
    value = counts[field]
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"calibration reporting {field} count must be an integer")
    counts[field] = value + 1
