"""Closed local feature/truth evidence for length training and evaluation."""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import NoReturn, cast

from vntyper.scripts.calibration_length import LengthTrainingRow
from vntyper.scripts.calibration_length_evaluation import EvaluationPhase, LengthEvaluationRow
from vntyper.scripts.calibration_length_metrics import LengthEligibleRoster, length_eligible_roster_document
from vntyper.scripts.calibration_length_protocol import LengthProtocol, length_protocol_document
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_features import encode_length_features
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION

logger = logging.getLogger(__name__)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_PHASES = {"policy-selection", "validation", "locked-heldout", "development-assessment"}
_ROW_FIELDS = {
    "key",
    "group_key",
    "features_sha256",
    "truth_boundary_definition",
    "total_truth_repeat_count",
    "evidence_domain",
}
_TRAINING_FIELDS = {
    "schema_version",
    "role",
    "study_sha256",
    "partition_sha256",
    "training_roster_sha256",
    "run_manifest_sha256",
    "rows",
    "sha256",
}
_ROLE_FIELDS = {
    "schema_version",
    "phase",
    "study_sha256",
    "partition_sha256",
    "protocol_sha256",
    "eligible_roster_sha256",
    "run_manifest_sha256",
    "rows",
    "sha256",
}


@dataclass(frozen=True)
class LengthTrainingEvidence:
    """Exact local training rows and their input/run commitments."""

    study_sha256: str
    partition_sha256: str
    training_roster_sha256: str
    run_manifest_sha256: str
    rows: tuple[LengthTrainingRow, ...]
    sha256: str


@dataclass(frozen=True)
class LengthRoleEvidence:
    """Exact local evaluation rows and their input/run commitments."""

    phase: EvaluationPhase
    study_sha256: str
    partition_sha256: str
    protocol_sha256: str
    eligible_roster_sha256: str
    run_manifest_sha256: str
    rows: tuple[LengthEvaluationRow, ...]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"length evidence {label} fields differ from the closed contract")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"length evidence {label} must be a lowercase SHA256 digest")
    return value


def _row_document(row: LengthTrainingRow | LengthEvaluationRow) -> dict[str, object]:
    if not isinstance(row.key, str) or not row.key or row.key != row.key.strip():
        _fail("length evidence row key must be non-empty trimmed text")
    if not isinstance(row.group_key, str) or not row.group_key or row.group_key != row.group_key.strip():
        _fail("length evidence row group key must be non-empty trimmed text")
    if row.truth_boundary_definition != TARGET_BOUNDARY_DEFINITION:
        _fail("length evidence truth boundary definition is incompatible")
    if not isinstance(row.evidence_domain, str) or row.evidence_domain not in {"synthetic", "external"}:
        _fail("length evidence domain must be synthetic or external")
    if isinstance(row.total_truth_repeat_count, bool) or not isinstance(row.total_truth_repeat_count, (int, float)):
        _fail("length evidence truth must be a positive integral repeat count")
    try:
        truth = float(row.total_truth_repeat_count)
    except OverflowError:
        _fail("length evidence truth must be a positive integral repeat count")
    if not math.isfinite(truth) or truth <= 0 or not truth.is_integer():
        _fail("length evidence truth must be a positive integral repeat count")
    encode_length_features(row.features)
    if row.features.manifest_key != row.key:
        _fail("length evidence feature manifest key differs from its row key")
    return {
        "key": row.key,
        "group_key": row.group_key,
        "features_sha256": row.features.sha256,
        "truth_boundary_definition": row.truth_boundary_definition,
        "total_truth_repeat_count": row.total_truth_repeat_count,
        "evidence_domain": row.evidence_domain,
    }


def _checked_rows(
    rows: Sequence[LengthTrainingRow] | Sequence[LengthEvaluationRow],
    roster: LengthEligibleRoster,
    expected_type: type[LengthTrainingRow] | type[LengthEvaluationRow],
) -> tuple[LengthTrainingRow | LengthEvaluationRow, ...]:
    length_eligible_roster_document(roster)
    if not isinstance(rows, (tuple, list)) or not rows:
        _fail("length evidence rows must be a non-empty typed sequence")
    if any(not isinstance(row, expected_type) for row in rows):
        _fail("length evidence rows contain the wrong typed row kind")
    checked = tuple(rows)
    if tuple(row.group_key for row in checked) != tuple(sorted(row.group_key for row in checked)):
        _fail("length evidence typed rows must be sorted by group key")
    documents = tuple(_row_document(row) for row in checked)
    if len({row["key"] for row in documents}) != len(documents) or len({row["group_key"] for row in documents}) != len(
        documents
    ):
        _fail("length evidence rows must have unique keys and groups")
    expected = {(member.key, member.group_key) for member in roster.members}
    if {(row.key, row.group_key) for row in checked} != expected:
        _fail("length evidence rows must match the frozen roster exactly")
    return checked


def _training_payload(evidence: LengthTrainingEvidence) -> dict[str, object]:
    return {
        "schema_version": "calibration-length-training-evidence-v1",
        "role": "training",
        "study_sha256": evidence.study_sha256,
        "partition_sha256": evidence.partition_sha256,
        "training_roster_sha256": evidence.training_roster_sha256,
        "run_manifest_sha256": evidence.run_manifest_sha256,
        "rows": [_row_document(row) for row in evidence.rows],
    }


def bind_length_training_evidence(
    rows: Sequence[LengthTrainingRow],
    roster: LengthEligibleRoster,
    *,
    study_sha256: str,
    partition_sha256: str,
    run_manifest_sha256: str,
) -> LengthTrainingEvidence:
    """Bind exact opened training rows before fitting.

    Args:
        rows: Complete ordered training rows from the authorized adapter.
        roster: Exact separately frozen training roster.
        study_sha256: Exact target study identity.
        partition_sha256: Exact private training partition identity.
        run_manifest_sha256: Exact measurement run manifest identity.

    Returns:
        Immutable training evidence whose digest metadata must bind.

    Raises:
        ValueError: If identities, rows or roster differ.
    """
    checked = _checked_rows(rows, roster, LengthTrainingRow)
    training_rows = cast(tuple[LengthTrainingRow, ...], checked)
    if any(row.role != "training" for row in training_rows):
        _fail("length training evidence rows must have the training role")
    evidence = LengthTrainingEvidence(
        _digest(study_sha256, "study digest"),
        _digest(partition_sha256, "partition digest"),
        roster.sha256,
        _digest(run_manifest_sha256, "run manifest digest"),
        training_rows,
        "",
    )
    return replace(evidence, sha256=canonical_sha256(_training_payload(evidence)))


def length_training_evidence_document(evidence: LengthTrainingEvidence) -> dict[str, object]:
    """Project immutable training evidence without feature bodies or paths.

    Args:
        evidence: Previously bound typed training evidence.

    Returns:
        Closed canonical JSON-compatible local document.

    Raises:
        ValueError: If typed content or digest differs.
    """
    if not isinstance(evidence, LengthTrainingEvidence) or not isinstance(evidence.rows, tuple):
        _fail("length training evidence projection requires immutable typed evidence")
    for name in ("study_sha256", "partition_sha256", "training_roster_sha256", "run_manifest_sha256"):
        _digest(getattr(evidence, name), name)
    if any(not isinstance(row, LengthTrainingRow) or row.role != "training" for row in evidence.rows):
        _fail("length training evidence projection requires training rows")
    _internal_rows(evidence.rows)
    payload = _training_payload(evidence)
    if canonical_sha256(payload) != evidence.sha256:
        _fail("length training evidence digest differs from its canonical content")
    return {**payload, "sha256": evidence.sha256}


def decode_length_training_evidence(
    value: object, *, rows: Sequence[LengthTrainingRow], roster: LengthEligibleRoster
) -> LengthTrainingEvidence:
    """Decode training evidence by re-deriving it from opened rows.

    Args:
        value: Parsed closed training evidence document.
        rows: Already opened typed training rows.
        roster: Exact separately frozen training roster.

    Returns:
        Immutable evidence retaining the opened rows.

    Raises:
        ValueError: If serialized and opened evidence differ.
    """
    raw = _object(value, _TRAINING_FIELDS, "training root")
    if raw["schema_version"] != "calibration-length-training-evidence-v1" or raw["role"] != "training":
        _fail("length training evidence schema or role is unsupported")
    _raw_rows(raw["rows"])
    evidence = bind_length_training_evidence(
        rows,
        roster,
        study_sha256=_digest(raw["study_sha256"], "study digest"),
        partition_sha256=_digest(raw["partition_sha256"], "partition digest"),
        run_manifest_sha256=_digest(raw["run_manifest_sha256"], "run manifest digest"),
    )
    if dict(raw) != length_training_evidence_document(evidence):
        _fail("length training evidence differs from the supplied typed rows")
    return evidence


def _role_payload(evidence: LengthRoleEvidence) -> dict[str, object]:
    return {
        "schema_version": "calibration-length-role-evidence-v1",
        "phase": evidence.phase,
        "study_sha256": evidence.study_sha256,
        "partition_sha256": evidence.partition_sha256,
        "protocol_sha256": evidence.protocol_sha256,
        "eligible_roster_sha256": evidence.eligible_roster_sha256,
        "run_manifest_sha256": evidence.run_manifest_sha256,
        "rows": [_row_document(row) for row in evidence.rows],
    }


def bind_length_role_evidence(
    rows: Sequence[LengthEvaluationRow],
    roster: LengthEligibleRoster,
    protocol: LengthProtocol,
    *,
    study_sha256: str,
    partition_sha256: str,
    run_manifest_sha256: str,
) -> LengthRoleEvidence:
    """Bind exact opened feature/truth rows to one evaluation role.

    Args:
        rows: Complete ordered evaluation rows from the authorized adapter.
        roster: Exact role-specific eligible roster.
        protocol: Frozen role-specific evaluation protocol.
        study_sha256: Study or development-assessment declaration identity.
        partition_sha256: Exact private input partition identity.
        run_manifest_sha256: Exact measurement run manifest identity.

    Returns:
        Immutable local role evidence with feature-content commitments.

    Raises:
        ValueError: If identities, rows, protocol or roster differ.
    """
    length_protocol_document(protocol)
    if roster.sha256 != protocol.eligible_roster_sha256:
        _fail("length role evidence roster differs from the protocol")
    checked = cast(tuple[LengthEvaluationRow, ...], _checked_rows(rows, roster, LengthEvaluationRow))
    phases = {row.phase for row in checked}
    if len(phases) != 1 or next(iter(phases)) not in _PHASES:
        _fail("length role evidence rows must use one supported phase")
    evidence = LengthRoleEvidence(
        cast(EvaluationPhase, next(iter(phases))),
        _digest(study_sha256, "study digest"),
        _digest(partition_sha256, "partition digest"),
        protocol.sha256,
        roster.sha256,
        _digest(run_manifest_sha256, "run manifest digest"),
        checked,
        "",
    )
    return replace(evidence, sha256=canonical_sha256(_role_payload(evidence)))


def length_role_evidence_document(evidence: LengthRoleEvidence) -> dict[str, object]:
    """Project local role evidence without feature bodies or source paths.

    Args:
        evidence: Previously bound typed role evidence.

    Returns:
        Closed canonical JSON-compatible local document.

    Raises:
        ValueError: If typed content is mutable or its digest is stale.
    """
    if not isinstance(evidence, LengthRoleEvidence) or not isinstance(evidence.rows, tuple):
        _fail("length role evidence projection requires immutable typed evidence")
    for name in (
        "study_sha256",
        "partition_sha256",
        "protocol_sha256",
        "eligible_roster_sha256",
        "run_manifest_sha256",
    ):
        _digest(getattr(evidence, name), name)
    if (
        not isinstance(evidence.phase, str)
        or evidence.phase not in _PHASES
        or any(not isinstance(row, LengthEvaluationRow) or row.phase != evidence.phase for row in evidence.rows)
    ):
        _fail("length role evidence projection phase or row kind differs")
    _internal_rows(evidence.rows)
    payload = _role_payload(evidence)
    if canonical_sha256(payload) != evidence.sha256:
        _fail("length role evidence digest differs from its canonical content")
    return {**payload, "sha256": evidence.sha256}


def _raw_rows(value: object) -> None:
    if not isinstance(value, list) or any(not isinstance(row, Mapping) or set(row) != _ROW_FIELDS for row in value):
        _fail("length evidence row fields differ from the closed contract")


def _internal_rows(rows: tuple[LengthTrainingRow, ...] | tuple[LengthEvaluationRow, ...]) -> None:
    for row in rows:
        _row_document(row)
    group_keys = tuple(row.group_key for row in rows)
    if group_keys != tuple(sorted(group_keys)):
        _fail("length evidence projected rows must be sorted by group key")
    if len({row.key for row in rows}) != len(rows) or len(set(group_keys)) != len(rows):
        _fail("length evidence projected rows must have unique keys and groups")


def decode_length_role_evidence(
    value: object,
    *,
    rows: Sequence[LengthEvaluationRow],
    roster: LengthEligibleRoster,
    protocol: LengthProtocol,
) -> LengthRoleEvidence:
    """Decode role evidence by re-deriving it from opened feature rows.

    Args:
        value: Parsed closed role evidence document.
        rows: Already opened typed rows supplied by the controller.
        roster: Exact role-specific eligible roster.
        protocol: Exact role-specific protocol.

    Returns:
        Immutable evidence retaining the opened typed rows.

    Raises:
        ValueError: If serialized and opened evidence differ in any field.
    """
    raw = _object(value, _ROLE_FIELDS, "role root")
    if raw["schema_version"] != "calibration-length-role-evidence-v1":
        _fail("length role evidence schema is unsupported")
    _raw_rows(raw["rows"])
    phase = raw["phase"]
    if not isinstance(phase, str) or phase not in _PHASES:
        _fail("length role evidence phase is unsupported")
    evidence = bind_length_role_evidence(
        rows,
        roster,
        protocol,
        study_sha256=_digest(raw["study_sha256"], "study digest"),
        partition_sha256=_digest(raw["partition_sha256"], "partition digest"),
        run_manifest_sha256=_digest(raw["run_manifest_sha256"], "run manifest digest"),
    )
    if dict(raw) != length_role_evidence_document(evidence):
        _fail("length role evidence differs from the supplied typed rows")
    return evidence
