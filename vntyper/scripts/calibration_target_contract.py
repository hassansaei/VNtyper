"""Closed target studies with outcome-independent baseline plans and four roles.

A length baseline plan declares a training-only recipe. Its digest is deliberately
separate from the fitted training-mean artifact used by a candidate's baseline
binding. Neither a declaration nor its exposure-ledger identity grants access.
"""

from __future__ import annotations

import logging
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_caller_policy import (
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_caller_protocol import CallerProtocol, caller_protocol_document, decode_caller_protocol
from vntyper.scripts.calibration_candidate import (
    CandidateApplicability,
    CandidateProducer,
    candidate_applicability_document,
    candidate_producer_document,
    decode_candidate_applicability,
    decode_candidate_producer,
)
from vntyper.scripts.calibration_exposure import require_digest
from vntyper.scripts.calibration_length_protocol import LengthProtocol, decode_length_protocol, length_protocol_document
from vntyper.scripts.calibration_manifest import PartitionManifest, decode_partition_manifest
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)
Target = Literal["callers", "length"]
_FIELDS = {"schema_version", "target", "protocol", "partitions", "baseline", "applicability", "exposure_ledger_id"}
_ROLES = {"training", "policy-selection", "validation", "locked-heldout"}


@dataclass(frozen=True)
class CallerBaselinePlan:
    """Full baseline policy and separately frozen capture/reference/tool assets."""

    policy: CallerPolicyValues
    assets_sha256: str
    producer: CandidateProducer
    sha256: str


@dataclass(frozen=True)
class LengthBaselinePlan:
    """Training-only baseline recipe and fixed measurement/numerical policy."""

    annotation_sha256: str
    counting_policy_sha256: str
    maximum_condition_number: float
    producer: CandidateProducer
    sha256: str


@dataclass(frozen=True)
class TargetStudy:
    """One immutable finite search declaration; never scientific approval."""

    target: Target
    protocol: CallerProtocol | LengthProtocol
    partitions: PartitionManifest
    baseline: CallerBaselinePlan | LengthBaselinePlan
    applicability: CandidateApplicability
    exposure_ledger_id: str
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"target study {label} fields differ")
    return value


def _caller_baseline(value: object) -> CallerBaselinePlan:
    raw = _object(value, {"schema_version", "policy", "assets_sha256", "producer"}, "caller baseline")
    if raw["schema_version"] != "calibration-caller-baseline-plan-v1":
        _fail("target study caller baseline schema is unsupported")
    return CallerBaselinePlan(
        decode_caller_policy_values(raw["policy"]),
        require_digest(raw["assets_sha256"], "baseline assets"),
        decode_candidate_producer(raw["producer"]),
        canonical_sha256(raw),
    )


def _length_baseline(value: object) -> LengthBaselinePlan:
    raw = _object(
        value,
        {
            "schema_version",
            "baseline_kind",
            "annotation_sha256",
            "counting_policy_sha256",
            "maximum_condition_number",
            "producer",
        },
        "length baseline",
    )
    if raw["schema_version"] != "calibration-length-baseline-plan-v1" or raw["baseline_kind"] != "training-mean-v1":
        _fail("target study length baseline must declare the training-only mean recipe")
    condition = raw["maximum_condition_number"]
    if (
        isinstance(condition, bool)
        or not isinstance(condition, (int, float))
        or not 1 < condition <= sys.float_info.max
    ):
        _fail("target study maximum condition number must be finite and greater than one")
    return LengthBaselinePlan(
        require_digest(raw["annotation_sha256"], "baseline annotation"),
        require_digest(raw["counting_policy_sha256"], "baseline counting policy"),
        float(condition),
        decode_candidate_producer(raw["producer"]),
        canonical_sha256(raw),
    )


def decode_target_study(value: object) -> TargetStudy:
    """Decode a target-specific search, baseline plan, and role declaration.

    Args:
        value: Closed calibration-study-v2 object declared before outcome access.

    Returns:
        Immutable target contract with independent study and baseline-plan digests.

    Raises:
        ValueError: For schema drift, incompatible targets, missing roles, leakage,
            undeclared assay classes, or inconsistent policy/applicability bindings.
    """
    raw = _object(value, _FIELDS, "declaration")
    if raw["schema_version"] != "calibration-study-v2":
        _fail("target study schema must be calibration-study-v2")
    target_value = raw["target"]
    if not isinstance(target_value, str) or target_value not in {"callers", "length"}:
        _fail("target study target must be callers or length")
    target = cast(Target, target_value)
    applicability = decode_candidate_applicability(raw["applicability"], target=target)
    baseline: CallerBaselinePlan | LengthBaselinePlan
    protocol: CallerProtocol | LengthProtocol
    if target == "callers":
        baseline = _caller_baseline(raw["baseline"])
        protocol = decode_caller_protocol(raw["protocol"], baseline_policy=baseline.policy)
        if baseline.policy.required_callers != applicability.required_callers:
            _fail("target study caller composition differs from applicability")
    else:
        baseline = _length_baseline(raw["baseline"])
        protocol = decode_length_protocol(raw["protocol"])
        if baseline.counting_policy_sha256 != applicability.counting_policy_sha256:
            _fail("target study length counting policy differs from applicability")
    partitions = decode_partition_manifest(raw["partitions"])
    if {member.role for member in partitions.members} != _ROLES:
        _fail("target study requires all four promotion roles")
    for member in partitions.members:
        if member.assay_class not in applicability.assay_classes:
            _fail("target study partition assay is outside declared applicability")
        if member.role == "locked-heldout" and member.provenance != "external-custodian":
            _fail("target study locked membership requires external-custodian provenance")
    return TargetStudy(
        target,
        protocol,
        partitions,
        baseline,
        applicability,
        require_digest(raw["exposure_ledger_id"], "ledger identity"),
        canonical_sha256(raw),
    )


def _baseline_document(baseline: CallerBaselinePlan | LengthBaselinePlan) -> dict[str, object]:
    if isinstance(baseline, CallerBaselinePlan):
        document: dict[str, object] = {
            "schema_version": "calibration-caller-baseline-plan-v1",
            "policy": caller_policy_values_document(baseline.policy),
            "assets_sha256": baseline.assets_sha256,
            "producer": candidate_producer_document(baseline.producer),
        }
        if _caller_baseline(document) != baseline:
            _fail("target study caller baseline identity differs")
        return document
    if not isinstance(baseline, LengthBaselinePlan):
        _fail("target study baseline must be typed")
    document = {
        "schema_version": "calibration-length-baseline-plan-v1",
        "baseline_kind": "training-mean-v1",
        "annotation_sha256": baseline.annotation_sha256,
        "counting_policy_sha256": baseline.counting_policy_sha256,
        "maximum_condition_number": baseline.maximum_condition_number,
        "producer": candidate_producer_document(baseline.producer),
    }
    if _length_baseline(document) != baseline:
        _fail("target study length baseline identity differs")
    return document


def target_study_document(study: TargetStudy) -> dict[str, object]:
    """Project an immutable study after revalidating all content and digests.

    Args:
        study: Previously decoded target study.

    Returns:
        Fresh canonical-content mapping with no paths or outcome data added.

    Raises:
        ValueError: If typed content, target, nested identities or hashes changed.
    """
    if not isinstance(study, TargetStudy):
        _fail("target study projection requires a typed declaration")
    partitions = {
        "schema_version": "calibration-partitions-v1",
        "members": [
            {
                "key": member.key,
                "role": member.role,
                "provenance": member.provenance,
                "assay_class": member.assay_class,
                "groups": {key: list(values) for key, values in member.groups.items()},
            }
            for member in study.partitions.members
        ],
    }
    if decode_partition_manifest(partitions) != study.partitions:
        _fail("target study partition identity differs")
    protocol = (
        caller_protocol_document(study.protocol)
        if isinstance(study.protocol, CallerProtocol)
        else length_protocol_document(study.protocol)
    )
    document: dict[str, object] = {
        "schema_version": "calibration-study-v2",
        "target": study.target,
        "protocol": protocol,
        "partitions": partitions,
        "baseline": _baseline_document(study.baseline),
        "applicability": candidate_applicability_document(study.applicability, target=study.target),
        "exposure_ledger_id": study.exposure_ledger_id,
    }
    if decode_target_study(document) != study:
        _fail("target study identity differs from its content")
    return document
