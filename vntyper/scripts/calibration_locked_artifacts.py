"""Locked-role declarations and externally supplied payload decoding."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from types import MappingProxyType
from typing import cast

from vntyper.scripts.calibration_artifact_io import freeze_json
from vntyper.scripts.calibration_features import decode_feature_artifact, decode_label_artifact
from vntyper.scripts.calibration_manifest import StudyDeclaration, decode_study_declaration
from vntyper.scripts.calibration_run_extraction import decode_run_hashes
from vntyper.scripts.calibration_workflow import ExtractedEvidence
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)


def build_locked_declaration_documents(
    study: StudyDeclaration,
    keys: tuple[str, ...],
    runs: Mapping[str, object],
) -> tuple[dict[str, object], dict[str, object]]:
    """Build value-free member and run commitments for a future custodian.

    Args:
        study: Validated complete study declaration.
        keys: Locked-role member keys in partition order.
        runs: Complete immutable run declarations supplied to ordinary extraction.

    Returns:
        Member-declaration and run-commitment documents.

    Raises:
        ValueError: If inputs or locked keys do not match the study.
    """
    if not isinstance(study, StudyDeclaration) or not isinstance(keys, tuple) or not isinstance(runs, Mapping):
        raise ValueError("locked calibration declaration requires a study, key tuple, and run mapping")
    expected = tuple(member.key for member in study.partitions.members if member.role == "locked-heldout")
    if keys != expected or any(key not in runs for key in keys):
        raise ValueError("locked calibration declaration keys differ from the partition or run commitments")
    members = {member.key: member for member in study.partitions.members}
    declaration: dict[str, object] = {
        "schema_version": "calibration-locked-member-declaration-v1",
        "study_sha256": study.sha256,
        "protocol_sha256": study.protocol.sha256,
        "partition_manifest_sha256": study.partitions.sha256,
        "members": [
            {
                "key": key,
                "role": "locked-heldout",
                "assay_class": members[key].assay_class,
                "groups": {name: list(values) for name, values in members[key].groups.items()},
            }
            for key in keys
        ],
    }
    commitments: dict[str, object] = {
        "schema_version": "calibration-locked-run-commitments-v1",
        "study_sha256": study.sha256,
        "runs": {key: runs[key] for key in keys},
    }
    return declaration, commitments


def decode_locked_payload(raw: bytes) -> ExtractedEvidence:
    """Decode custodian-opened locked bytes into immutable evaluation evidence.

    Args:
        raw: Exact locked payload bytes opened after the durable precommit.

    Returns:
        Immutable role-specific evidence and exact dataset digest.

    Raises:
        ValueError: If the payload schema or nested artifacts are invalid.
    """
    payload = load_strict_json_object(raw)
    expected = {"schema_version", "study", "features", "labels", "baseline", "run_hashes"}
    if set(payload) != expected or payload["schema_version"] != "calibration-locked-payload-v1":
        raise ValueError("locked calibration payload fields or schema version differ")
    study = decode_study_declaration(payload["study"])
    features = decode_feature_artifact(payload["features"])
    labels = decode_label_artifact(payload["labels"])
    hashes = decode_run_hashes(_mapping(payload["run_hashes"], "run hashes"))
    baseline = dict(_mapping(payload["baseline"], "baseline"))
    dataset = canonical_sha256(
        {
            "study_sha256": study.sha256,
            "features_sha256": features.sha256,
            "labels_sha256": labels.sha256,
            "baseline_sha256": canonical_sha256(baseline),
            "run_artifact_sha256": hashes,
        }
    )
    return ExtractedEvidence(
        study,
        features,
        labels,
        cast(Mapping[str, object], freeze_json(baseline)),
        MappingProxyType({key: MappingProxyType(dict(value)) for key, value in hashes.items()}),
        study.sha256,
        dataset,
        dataset,
    )


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"locked calibration {label} must be an object")
    return value
