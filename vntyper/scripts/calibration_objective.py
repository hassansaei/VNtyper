"""Literal safety-first calibration objective and admissibility rules."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction

from vntyper.scripts.calibration_contract import CalibrationProtocol, CandidateMetrics
from vntyper.scripts.canonical_json import canonical_sha256

_DOMINANCE_FIELDS = {
    "enabled",
    "minimum_record_count_margin",
    "minimum_record_share",
    "minimum_record_share_margin",
    "xd_veto",
    "abstain_on_inadmissible_advntr",
}
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class CandidateEvaluation:
    """Metrics plus protocol-fixed admissibility evidence for one candidate."""

    metrics: CandidateMetrics
    detection_lower_bound: Fraction
    macro_exact_lower_bound: Fraction
    stratum_counts: tuple[int, ...]
    holm_adjusted_p_value: Fraction

    def __post_init__(self) -> None:
        """Validate exact paired bounds, stratum counts, and multiplicity evidence."""
        if not isinstance(self.metrics, CandidateMetrics):
            raise ValueError("candidate evaluation metrics must be CandidateMetrics")
        for value, label in (
            (self.detection_lower_bound, "detection lower bound"),
            (self.macro_exact_lower_bound, "macro exact lower bound"),
            (self.holm_adjusted_p_value, "Holm-adjusted p-value"),
        ):
            if not isinstance(value, Fraction):
                raise ValueError(f"candidate {label} must be an exact Fraction")
        if not isinstance(self.stratum_counts, tuple) or not self.stratum_counts:
            raise ValueError("candidate stratum counts must be a non-empty tuple")
        if any(isinstance(count, bool) or not isinstance(count, int) or count < 0 for count in self.stratum_counts):
            raise ValueError("candidate stratum counts must be non-negative integers")
        if not 0 <= self.holm_adjusted_p_value <= 1:
            raise ValueError("candidate Holm-adjusted p-value must be between zero and one")
        if self.metrics.wrong_displayed_names_all_tiers < self.metrics.wrong_tier_a_displayed_names:
            raise ValueError("all-tier wrong displayed names must include wrong Tier-A displayed names")


@dataclass(frozen=True)
class OutcomeObservation:
    """One truth-bound candidate outcome used to calculate objective metrics."""

    key: str
    assay_class: str
    mutation_class: str
    expected_identity: str | None
    expected_display_name: str | None
    selected_identity: str | None
    displayed_name: str | None
    tier: str | None
    abstained: bool
    applicable: bool
    baseline_applicable: bool
    baseline_tier: str | None

    def __post_init__(self) -> None:
        """Validate the distinct mutated/control and selected/abstained states."""
        for value, label in (
            (self.key, "outcome key"),
            (self.assay_class, "assay class"),
            (self.mutation_class, "mutation class"),
        ):
            if not isinstance(value, str) or not value:
                raise ValueError(f"calibration {label} must be a non-empty string")
        for boolean_value, label in (
            (self.abstained, "abstained"),
            (self.applicable, "applicable"),
            (self.baseline_applicable, "baseline applicable"),
        ):
            if not isinstance(boolean_value, bool):
                raise ValueError(f"calibration outcome {label} must be Boolean")
        if (self.expected_identity is None) != (self.expected_display_name is None):
            raise ValueError("calibration expected identity and displayed name must be jointly present or absent")
        if (self.selected_identity is None) != (self.displayed_name is None):
            raise ValueError("calibration selected identity and displayed name must be jointly present or absent")
        if self.abstained and self.selected_identity is not None:
            raise ValueError("calibration abstention cannot carry a selected identity")
        if self.selected_identity is not None and self.tier is None:
            raise ValueError("calibration selected identity requires a fixed reconciliation tier")


@dataclass(frozen=True)
class ObjectiveSummary:
    """Candidate metrics plus counts for the predeclared strata."""

    metrics: CandidateMetrics
    stratum_counts: tuple[int, ...]


def calculate_metrics(
    observations: Sequence[OutcomeObservation],
    *,
    profile_sha256: str,
    free_parameter_count: int,
    required_strata: Sequence[str],
) -> ObjectiveSummary:
    """Calculate exact safety metrics without dropping abstentions from denominators."""
    if not isinstance(observations, Sequence) or not observations:
        raise ValueError("calibration objective requires outcome observations")
    if any(not isinstance(row, OutcomeObservation) for row in observations):
        raise ValueError("calibration objective rows must be OutcomeObservation values")
    if not isinstance(profile_sha256, str) or len(profile_sha256) != 64:
        raise ValueError("calibration objective profile SHA-256 must contain 64 characters")
    if _SHA256_PATTERN.fullmatch(profile_sha256) is None:
        raise ValueError("calibration objective profile SHA-256 must be lowercase hexadecimal")
    if isinstance(free_parameter_count, bool) or not isinstance(free_parameter_count, int) or free_parameter_count < 0:
        raise ValueError("calibration objective free parameter count must be a non-negative integer")
    if not isinstance(required_strata, Sequence) or not required_strata:
        raise ValueError("calibration objective requires predeclared strata")
    strata = tuple(required_strata)
    if any(not isinstance(stratum, str) or not stratum for stratum in strata) or len(strata) != len(set(strata)):
        raise ValueError("calibration objective strata must be unique non-empty strings")

    mutated = tuple(row for row in observations if row.expected_identity is not None)
    controls = tuple(row for row in observations if row.expected_identity is None)
    if not mutated:
        raise ValueError("calibration objective requires mutated truth members")
    observed_strata = {f"{row.assay_class}:{row.mutation_class}" for row in mutated}
    undeclared_strata = observed_strata - set(strata)
    if undeclared_strata:
        raise ValueError(f"calibration objective observations use undeclared strata: {sorted(undeclared_strata)}")
    exact_by_stratum: dict[str, tuple[int, int]] = {}
    for stratum in strata:
        members = tuple(row for row in mutated if f"{row.assay_class}:{row.mutation_class}" == stratum)
        exact_by_stratum[stratum] = (
            sum(row.selected_identity == row.expected_identity for row in members),
            len(members),
        )
    counts = tuple(total for _, total in exact_by_stratum.values())
    macro_exact = sum(
        (Fraction(exact, total) if total else Fraction(0) for exact, total in exact_by_stratum.values()),
        start=Fraction(0),
    ) / len(strata)
    detected = sum(row.selected_identity is not None for row in mutated)
    wrong_names = tuple(
        row for row in mutated if row.displayed_name is not None and row.displayed_name != row.expected_display_name
    )
    applicable_denominator = tuple(row for row in observations if row.baseline_applicable)
    abstentions = sum(row.abstained for row in applicable_denominator)
    tier_a_was_required = any(row.baseline_tier == "A" for row in mutated)
    tier_a_reachable = not tier_a_was_required or any(
        row.tier == "A" and row.selected_identity is not None for row in mutated
    )
    wrong_tier_a = sum(row.tier == "A" for row in wrong_names)
    control_findings = sum(row.selected_identity is not None for row in controls)
    detection_sensitivity = Fraction(detected, len(mutated))
    abstention_fraction = Fraction(abstentions, len(applicable_denominator)) if applicable_denominator else Fraction(0)
    applicability_matches = all(row.applicable == row.baseline_applicable for row in observations)
    payload = {
        "candidate_profile_sha256": profile_sha256,
        "wrong_tier_a_displayed_names": wrong_tier_a,
        "control_findings": control_findings,
        "wrong_displayed_names_all_tiers": len(wrong_names),
        "macro_exact_recovery": str(macro_exact),
        "binary_detection_sensitivity": str(detection_sensitivity),
        "free_parameter_count": free_parameter_count,
        "abstention_fraction": str(abstention_fraction),
        "tier_a_reachable": tier_a_reachable,
        "applicability_matches": applicability_matches,
    }
    metrics = CandidateMetrics(
        profile_sha256,
        wrong_tier_a,
        control_findings,
        len(wrong_names),
        macro_exact,
        detection_sensitivity,
        free_parameter_count,
        abstention_fraction,
        tier_a_reachable,
        applicability_matches,
        canonical_sha256(payload),
    )
    return ObjectiveSummary(metrics, counts)


def lexicographic_safety_key(metrics: CandidateMetrics) -> tuple[object, ...]:
    """Return the literal approved safety-first objective tuple."""
    if not isinstance(metrics, CandidateMetrics):
        raise ValueError("lexicographic objective metrics must be CandidateMetrics")
    return (
        metrics.wrong_tier_a_displayed_names,
        metrics.control_findings,
        metrics.wrong_displayed_names_all_tiers,
        -metrics.macro_exact_recovery,
        -metrics.binary_detection_sensitivity,
        metrics.free_parameter_count,
        metrics.candidate_profile_sha256,
    )


def count_free_parameters(component: Mapping[str, object]) -> int:
    """Count non-neutral generated dominance leaves."""
    if not isinstance(component, Mapping) or set(component) != _DOMINANCE_FIELDS:
        raise ValueError("generated dominance component fields differ from the closed contract")
    enabled = component["enabled"]
    count_margin = component["minimum_record_count_margin"]
    share = component["minimum_record_share"]
    share_margin = component["minimum_record_share_margin"]
    xd_veto = component["xd_veto"]
    advntr_veto = component["abstain_on_inadmissible_advntr"]
    if not isinstance(enabled, bool) or not isinstance(advntr_veto, bool):
        raise ValueError("generated dominance Boolean leaves must be Boolean")
    if isinstance(count_margin, bool) or not isinstance(count_margin, int) or count_margin < 0:
        raise ValueError("generated dominance count margin must be a non-negative integer")
    for value, label in ((share, "record share"), (share_margin, "record-share margin")):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ValueError(f"generated dominance {label} must be between zero and one")
    if xd_veto not in {"disabled", "missingness", "concentration", "discordance"}:
        raise ValueError("generated dominance XD veto is unsupported")
    return sum(
        (
            int(enabled),
            int(count_margin != 0),
            int(share != 0),
            int(share_margin != 0),
            int(xd_veto != "disabled"),
            int(advntr_veto),
        )
    )


def select_candidate(
    evaluations: Sequence[CandidateEvaluation],
    protocol: CalibrationProtocol,
) -> CandidateEvaluation | None:
    """Select the lexicographically best candidate after all hard constraints."""
    if not isinstance(protocol, CalibrationProtocol):
        raise ValueError("candidate selection protocol must be a CalibrationProtocol")
    if not isinstance(evaluations, Sequence) or any(
        not isinstance(evaluation, CandidateEvaluation) for evaluation in evaluations
    ):
        raise ValueError("candidate evaluations must contain CandidateEvaluation values")
    admissible = tuple(evaluation for evaluation in evaluations if _is_admissible(evaluation, protocol))
    if not admissible:
        return None
    return min(admissible, key=lambda evaluation: lexicographic_safety_key(evaluation.metrics))


def _is_admissible(evaluation: CandidateEvaluation, protocol: CalibrationProtocol) -> bool:
    metrics = evaluation.metrics
    return (
        metrics.wrong_tier_a_displayed_names == 0
        and metrics.control_findings == 0
        and metrics.tier_a_reachable
        and metrics.applicability_matches
        and metrics.abstention_fraction <= protocol.maximum_abstention_fraction
        and metrics.free_parameter_count <= protocol.maximum_free_parameters
        and len(evaluation.stratum_counts) == len(protocol.required_strata)
        and all(count >= protocol.minimum_stratum_count for count in evaluation.stratum_counts)
        and evaluation.detection_lower_bound >= 0
        and evaluation.macro_exact_lower_bound >= 0
        and evaluation.holm_adjusted_p_value <= Fraction(1, 20)
    )
