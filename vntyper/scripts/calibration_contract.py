"""Strict immutable contracts shared by the calibration workflow."""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.canonical_json import canonical_sha256

CalibrationRole = Literal["training", "policy-selection", "validation", "locked-heldout"]
EvidenceProvenance = Literal["development", "external-custodian"]

_ROLES = frozenset({"training", "policy-selection", "validation", "locked-heldout"})
_PROVENANCE = frozenset({"development", "external-custodian"})
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_PROTOCOL_FIELDS = {
    "objective",
    "bootstrap_iterations",
    "bootstrap_interval",
    "multiplicity_method",
    "seed",
    "maximum_free_parameters",
    "minimum_stratum_count",
    "maximum_abstention_fraction",
    "assay_classes",
    "mutation_classes",
    "candidate_grid",
}
_GRID_FIELDS = {
    "minimum_record_count_margin",
    "minimum_record_share",
    "minimum_record_share_margin",
    "xd_veto",
}
_DOMINANCE_FIELDS = {
    "enabled",
    "minimum_record_count_margin",
    "minimum_record_share",
    "minimum_record_share_margin",
    "xd_veto",
    "abstain_on_inadmissible_advntr",
}


@dataclass(frozen=True)
class CandidateRule:
    """One member of the protocol-frozen finite dominance grid."""

    minimum_record_count_margin: int
    minimum_record_share: Fraction
    minimum_record_share_margin: Fraction
    xd_veto: str

    def as_component(self) -> Mapping[str, object]:
        """Project this candidate into a generated dominance component."""
        return MappingProxyType(
            {
                "enabled": True,
                "minimum_record_count_margin": self.minimum_record_count_margin,
                "minimum_record_share": float(self.minimum_record_share),
                "minimum_record_share_margin": float(self.minimum_record_share_margin),
                "xd_veto": self.xd_veto,
                "abstain_on_inadmissible_advntr": False,
            }
        )


@dataclass(frozen=True)
class CalibrationProtocol:
    """Validated version-1 study protocol and complete candidate family."""

    objective: str
    bootstrap_iterations: int
    bootstrap_interval: str
    multiplicity_method: str
    seed: int
    maximum_free_parameters: int
    minimum_stratum_count: int
    maximum_abstention_fraction: Fraction
    assay_classes: tuple[str, ...]
    mutation_classes: tuple[str, ...]
    required_strata: tuple[str, ...]
    candidate_grid: Mapping[str, tuple[object, ...]]
    candidates: tuple[CandidateRule, ...]
    sha256: str


@dataclass(frozen=True)
class CalibrationPolicy:
    """Hash-bound generated dominance policy artifact."""

    profile_sha256: str
    dominance: Mapping[str, object]
    sha256: str


@dataclass(frozen=True)
class EvidenceManifest:
    """Hashes binding one role-specific evidence artifact."""

    role: CalibrationRole
    provenance: EvidenceProvenance
    protocol_sha256: str
    partition_manifest_sha256: str
    features_sha256: str
    labels_sha256: str
    baseline_sha256: str
    sha256: str


@dataclass(frozen=True)
class CandidateMetrics:
    """Exact objective inputs for one candidate profile."""

    candidate_profile_sha256: str
    wrong_tier_a_displayed_names: int
    control_findings: int
    wrong_displayed_names_all_tiers: int
    macro_exact_recovery: Fraction
    binary_detection_sensitivity: Fraction
    free_parameter_count: int
    abstention_fraction: Fraction
    tier_a_reachable: bool
    applicability_matches: bool
    sha256: str


@dataclass(frozen=True)
class CalibrationAttestation:
    """Role- and hash-bound validation or held-out result."""

    role: CalibrationRole
    status: str
    profile_sha256: str
    protocol_sha256: str
    evidence_sha256: str
    metrics_sha256: str
    sha256: str


def validate_evidence_bindings(
    evidence: EvidenceManifest,
    *,
    protocol_sha256: str,
    partition_manifest_sha256: str,
    features_sha256: str,
    labels_sha256: str,
    baseline_sha256: str,
) -> None:
    """Require every evidence hash to match the artifact opened by the caller."""
    if not isinstance(evidence, EvidenceManifest):
        raise ValueError("calibration evidence must be an EvidenceManifest")
    expected = (
        _digest(protocol_sha256, "opened protocol"),
        _digest(partition_manifest_sha256, "opened partition manifest"),
        _digest(features_sha256, "opened features"),
        _digest(labels_sha256, "opened labels"),
        _digest(baseline_sha256, "opened baseline"),
    )
    actual = (
        evidence.protocol_sha256,
        evidence.partition_manifest_sha256,
        evidence.features_sha256,
        evidence.labels_sha256,
        evidence.baseline_sha256,
    )
    if actual != expected:
        raise ValueError("calibration evidence hashes do not match the opened artifacts")


def validate_attestation_bindings(
    attestation: CalibrationAttestation,
    *,
    role: CalibrationRole,
    profile_sha256: str,
    protocol_sha256: str,
    evidence_sha256: str,
    metrics_sha256: str,
) -> None:
    """Require an attestation to bind the exact role and evaluated artifacts."""
    if not isinstance(attestation, CalibrationAttestation):
        raise ValueError("calibration attestation must be a CalibrationAttestation")
    expected = (
        role,
        _digest(profile_sha256, "evaluated profile"),
        _digest(protocol_sha256, "evaluated protocol"),
        _digest(evidence_sha256, "evaluated evidence"),
        _digest(metrics_sha256, "evaluated metrics"),
    )
    actual = (
        attestation.role,
        attestation.profile_sha256,
        attestation.protocol_sha256,
        attestation.evidence_sha256,
        attestation.metrics_sha256,
    )
    if actual != expected:
        raise ValueError("calibration attestation role or hashes do not match the evaluated artifacts")


def decode_protocol(value: object) -> CalibrationProtocol:
    """Decode the strict finite version-1 calibration protocol."""
    if isinstance(value, Mapping):
        for field, label in (("assay_classes", "assay classes"), ("mutation_classes", "mutation classes")):
            if field not in value:
                raise ValueError(f"calibration protocol requires predeclared {label}")
    raw = _exact_object(value, _PROTOCOL_FIELDS, "calibration protocol")
    if raw["objective"] != "lexicographic-safety-v1":
        raise ValueError("calibration protocol objective must be lexicographic-safety-v1")
    if raw["bootstrap_iterations"] != 10_000:
        raise ValueError("calibration protocol bootstrap iterations must be 10000")
    if raw["bootstrap_interval"] != "percentile":
        raise ValueError("calibration protocol bootstrap interval must be percentile")
    if raw["multiplicity_method"] != "holm":
        raise ValueError("calibration protocol multiplicity method must be holm")
    seed = _nonnegative_integer(raw["seed"], "calibration protocol seed")
    maximum_parameters = _nonnegative_integer(
        raw["maximum_free_parameters"], "calibration protocol maximum free parameters"
    )
    minimum_stratum_count = _positive_integer(
        raw["minimum_stratum_count"], "calibration protocol minimum stratum count"
    )
    maximum_abstention = _unit_fraction(
        raw["maximum_abstention_fraction"], "calibration protocol maximum abstention fraction"
    )
    assay_classes = _sorted_unique_strings(raw["assay_classes"], "calibration protocol assay classes")
    mutation_classes = _sorted_unique_strings(raw["mutation_classes"], "calibration protocol mutation classes")
    required_strata = tuple(f"{assay}:{mutation}" for assay, mutation in product(assay_classes, mutation_classes))
    grid = _decode_grid(raw["candidate_grid"])
    count_margins = cast(tuple[int, ...], grid["minimum_record_count_margin"])
    shares = cast(tuple[Fraction, ...], grid["minimum_record_share"])
    share_margins = cast(tuple[Fraction, ...], grid["minimum_record_share_margin"])
    xd_vetoes = cast(tuple[str, ...], grid["xd_veto"])
    candidates = tuple(
        CandidateRule(count_margin, share, share_margin, xd_veto)
        for count_margin, share, share_margin, xd_veto in product(
            count_margins,
            shares,
            share_margins,
            xd_vetoes,
        )
    )
    return CalibrationProtocol(
        objective="lexicographic-safety-v1",
        bootstrap_iterations=10_000,
        bootstrap_interval="percentile",
        multiplicity_method="holm",
        seed=seed,
        maximum_free_parameters=maximum_parameters,
        minimum_stratum_count=minimum_stratum_count,
        maximum_abstention_fraction=maximum_abstention,
        assay_classes=assay_classes,
        mutation_classes=mutation_classes,
        required_strata=required_strata,
        candidate_grid=MappingProxyType(grid),
        candidates=candidates,
        sha256=canonical_sha256(raw),
    )


def decode_policy(value: object) -> CalibrationPolicy:
    """Decode one strict generated-policy artifact."""
    raw = _exact_object(value, {"schema_version", "profile_sha256", "dominance"}, "calibration policy")
    _version(raw, "calibration-policy-v1", "calibration policy")
    profile_sha256 = _digest(raw["profile_sha256"], "calibration policy profile")
    dominance = _decode_dominance(raw["dominance"])
    return CalibrationPolicy(profile_sha256, MappingProxyType(dominance), canonical_sha256(raw))


def decode_evidence_manifest(value: object) -> EvidenceManifest:
    """Decode a strict role-bound evidence manifest."""
    raw = _exact_object(
        value,
        {
            "schema_version",
            "role",
            "provenance",
            "protocol_sha256",
            "partition_manifest_sha256",
            "features_sha256",
            "labels_sha256",
            "baseline_sha256",
        },
        "calibration evidence manifest",
    )
    _version(raw, "calibration-evidence-v1", "calibration evidence manifest")
    role = _role(raw["role"])
    provenance = _provenance(raw["provenance"])
    if role == "locked-heldout" and provenance != "external-custodian":
        raise ValueError("locked held-out evidence requires external custodian provenance")
    return EvidenceManifest(
        role,
        provenance,
        _digest(raw["protocol_sha256"], "evidence protocol"),
        _digest(raw["partition_manifest_sha256"], "evidence partition manifest"),
        _digest(raw["features_sha256"], "evidence features"),
        _digest(raw["labels_sha256"], "evidence labels"),
        _digest(raw["baseline_sha256"], "evidence baseline"),
        canonical_sha256(raw),
    )


def decode_metrics(value: object) -> CandidateMetrics:
    """Decode exact candidate metrics without floating-point objective rates."""
    raw = _exact_object(
        value,
        {
            "schema_version",
            "candidate_profile_sha256",
            "wrong_tier_a_displayed_names",
            "control_findings",
            "wrong_displayed_names_all_tiers",
            "macro_exact_recovery",
            "binary_detection_sensitivity",
            "free_parameter_count",
            "abstention_fraction",
            "tier_a_reachable",
            "applicability_matches",
        },
        "calibration metrics",
    )
    _version(raw, "calibration-metrics-v1", "calibration metrics")
    tier_a_reachable = _boolean(raw["tier_a_reachable"], "metrics Tier A reachability")
    applicability_matches = _boolean(raw["applicability_matches"], "metrics applicability match")
    return CandidateMetrics(
        _digest(raw["candidate_profile_sha256"], "metrics candidate profile"),
        _nonnegative_integer(raw["wrong_tier_a_displayed_names"], "wrong Tier-A displayed names"),
        _nonnegative_integer(raw["control_findings"], "control findings"),
        _nonnegative_integer(raw["wrong_displayed_names_all_tiers"], "wrong displayed names"),
        _fraction(raw["macro_exact_recovery"], "macro exact recovery"),
        _fraction(raw["binary_detection_sensitivity"], "binary detection sensitivity"),
        _nonnegative_integer(raw["free_parameter_count"], "free parameter count"),
        _unit_fraction(raw["abstention_fraction"], "abstention fraction"),
        tier_a_reachable,
        applicability_matches,
        canonical_sha256(raw),
    )


def decode_attestation(value: object) -> CalibrationAttestation:
    """Decode a validation or held-out hash-bound attestation."""
    raw = _exact_object(
        value,
        {
            "schema_version",
            "role",
            "status",
            "profile_sha256",
            "protocol_sha256",
            "evidence_sha256",
            "metrics_sha256",
        },
        "calibration attestation",
    )
    _version(raw, "calibration-attestation-v1", "calibration attestation")
    role = _role(raw["role"])
    if role not in {"validation", "locked-heldout"}:
        raise ValueError("calibration attestation role must be validation or locked-heldout")
    status = raw["status"]
    if status not in {"passed", "failed"}:
        raise ValueError("calibration attestation status must be passed or failed")
    return CalibrationAttestation(
        role,
        cast(str, status),
        _digest(raw["profile_sha256"], "attestation profile"),
        _digest(raw["protocol_sha256"], "attestation protocol"),
        _digest(raw["evidence_sha256"], "attestation evidence"),
        _digest(raw["metrics_sha256"], "attestation metrics"),
        canonical_sha256(raw),
    )


def _decode_grid(value: object) -> dict[str, tuple[object, ...]]:
    raw = _exact_object(value, _GRID_FIELDS, "calibration candidate grid")
    count_margins = _unique_list(raw["minimum_record_count_margin"], "candidate count margins")
    shares = _unique_list(raw["minimum_record_share"], "candidate record shares")
    share_margins = _unique_list(raw["minimum_record_share_margin"], "candidate share margins")
    xd_vetoes = _unique_list(raw["xd_veto"], "candidate XD vetoes")
    return {
        "minimum_record_count_margin": tuple(
            _nonnegative_integer(item, "candidate count margin") for item in count_margins
        ),
        "minimum_record_share": tuple(_unit_fraction(item, "candidate record share") for item in shares),
        "minimum_record_share_margin": tuple(
            _unit_fraction(item, "candidate record-share margin") for item in share_margins
        ),
        "xd_veto": tuple(_calibration_xd_veto(item) for item in xd_vetoes),
    }


def _decode_dominance(value: object) -> dict[str, object]:
    raw = _exact_object(value, _DOMINANCE_FIELDS, "calibration dominance policy")
    return {
        "enabled": _boolean(raw["enabled"], "dominance enabled"),
        "minimum_record_count_margin": _nonnegative_integer(
            raw["minimum_record_count_margin"], "dominance count margin"
        ),
        "minimum_record_share": float(_unit_fraction(raw["minimum_record_share"], "dominance record share")),
        "minimum_record_share_margin": float(
            _unit_fraction(raw["minimum_record_share_margin"], "dominance record-share margin")
        ),
        "xd_veto": _xd_veto(raw["xd_veto"]),
        "abstain_on_inadmissible_advntr": _boolean(raw["abstain_on_inadmissible_advntr"], "dominance adVNTR veto"),
    }


def _exact_object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    return value


def _version(value: Mapping[str, object], expected: str, label: str) -> None:
    if value["schema_version"] != expected:
        raise ValueError(f"{label} schema version must be {expected}")


def _unique_list(value: object, label: str) -> list[object]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{label} must be a non-empty list")
    if len({repr(item) for item in value}) != len(value):
        raise ValueError(f"{label} must not contain duplicate values")
    return value


def _sorted_unique_strings(value: object, label: str) -> tuple[str, ...]:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or value != sorted(value)
        or len(value) != len(set(value))
    ):
        raise ValueError(f"{label} must be sorted unique non-empty strings")
    return tuple(value)


def _nonnegative_integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{label} must be a non-negative integer")
    return value


def _positive_integer(value: object, label: str) -> int:
    result = _nonnegative_integer(value, label)
    if result == 0:
        raise ValueError(f"{label} must be positive")
    return result


def _fraction(value: object, label: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{label} must be an exact fraction")
    try:
        parsed = Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError) as error:
        raise ValueError(f"{label} must be an exact fraction") from error
    return parsed


def _unit_fraction(value: object, label: str) -> Fraction:
    parsed = _fraction(str(value) if isinstance(value, float) else value, label)
    if not 0 <= parsed <= 1:
        raise ValueError(f"{label} must be between zero and one")
    return parsed


def _xd_veto(value: object) -> str:
    if value not in {"disabled", "missingness", "concentration", "discordance"}:
        raise ValueError(f"candidate dominance XD veto is unsupported: {value!r}")
    return cast(str, value)


def _calibration_xd_veto(value: object) -> str:
    parsed = _xd_veto(value)
    if parsed not in {"disabled", "missingness"}:
        raise ValueError(f"candidate XD veto cannot be replayed from calibration-v1 scalar evidence: {parsed}")
    return parsed


def _boolean(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be Boolean")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} SHA-256 must be 64 lowercase hexadecimal characters")
    return value


def _role(value: object) -> CalibrationRole:
    if value not in _ROLES:
        raise ValueError(f"unsupported calibration role: {value!r}")
    return cast(CalibrationRole, value)


def _provenance(value: object) -> EvidenceProvenance:
    if value not in _PROVENANCE:
        raise ValueError(f"unsupported calibration evidence provenance: {value!r}")
    return cast(EvidenceProvenance, value)
