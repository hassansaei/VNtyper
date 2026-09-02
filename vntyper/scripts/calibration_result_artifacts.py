"""Pure closed encoders for calibration metrics and result attestations."""

from __future__ import annotations

import logging

from vntyper.scripts.calibration_contract import CandidateMetrics
from vntyper.scripts.calibration_objective import CandidateEvaluation
from vntyper.scripts.calibration_workflow import ExtractedEvidence
from vntyper.scripts.decision_profile import ResolvedDecisionProfile

logger = logging.getLogger(__name__)


def encode_metrics(metrics: CandidateMetrics) -> dict[str, object]:
    """Encode exact candidate metrics into the strict version-1 artifact."""
    return {
        "schema_version": "calibration-metrics-v1",
        "candidate_profile_sha256": metrics.candidate_profile_sha256,
        "wrong_tier_a_displayed_names": metrics.wrong_tier_a_displayed_names,
        "control_findings": metrics.control_findings,
        "wrong_displayed_names_all_tiers": metrics.wrong_displayed_names_all_tiers,
        "macro_exact_recovery": str(metrics.macro_exact_recovery),
        "binary_detection_sensitivity": str(metrics.binary_detection_sensitivity),
        "free_parameter_count": metrics.free_parameter_count,
        "abstention_fraction": str(metrics.abstention_fraction),
        "tier_a_reachable": metrics.tier_a_reachable,
        "applicability_matches": metrics.applicability_matches,
    }


def encode_evaluation(evaluation: CandidateEvaluation) -> dict[str, object]:
    """Encode one candidate evaluation and exact statistical bounds."""
    return {
        "schema_version": "calibration-candidate-evaluation-v1",
        "metrics_sha256": evaluation.metrics.sha256,
        "detection_lower_bound": str(evaluation.detection_lower_bound),
        "macro_exact_lower_bound": str(evaluation.macro_exact_lower_bound),
        "stratum_counts": list(evaluation.stratum_counts),
        "holm_adjusted_p_value": str(evaluation.holm_adjusted_p_value),
    }


def encode_attestation(
    role: str,
    profile: ResolvedDecisionProfile,
    evidence: ExtractedEvidence,
    metrics_sha256: str,
    passed: bool,
) -> dict[str, object]:
    """Encode a profile-, protocol-, evidence-, and metrics-bound result."""
    return encode_attestation_hashes(
        role,
        profile.digest,
        evidence.study.protocol.sha256,
        evidence.dataset_sha256,
        metrics_sha256,
        passed,
    )


def encode_attestation_hashes(
    role: str,
    profile_sha256: str,
    protocol_sha256: str,
    evidence_sha256: str,
    metrics_sha256: str,
    passed: bool,
) -> dict[str, object]:
    """Encode an attestation from already validated exact hashes."""
    return {
        "schema_version": "calibration-attestation-v1",
        "role": role,
        "status": "passed" if passed else "failed",
        "profile_sha256": profile_sha256,
        "protocol_sha256": protocol_sha256,
        "evidence_sha256": evidence_sha256,
        "metrics_sha256": metrics_sha256,
    }
