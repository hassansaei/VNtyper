"""Pure normalization and role projection of predeclared calibration runs."""

from __future__ import annotations

import logging
from collections.abc import Mapping

from vntyper.scripts.calibration_manifest import StudyDeclaration
from vntyper.scripts.calibration_run_extraction import decode_run_artifact_declaration

logger = logging.getLogger(__name__)

ROLES = ("training", "policy-selection", "validation", "locked-heldout")


def normalize_run_commitments(study: StudyDeclaration, runs: Mapping[str, object]) -> dict[str, object]:
    """Validate every predeclared run without opening any run-root artifact."""
    if not isinstance(study, StudyDeclaration) or not isinstance(runs, Mapping):
        raise ValueError("calibration run commitments require a study and run mapping")
    expected = tuple(member.key for member in study.partitions.members)
    if set(runs) != set(expected):
        raise ValueError("calibration run commitments must match the complete partition manifest")
    normalized: dict[str, object] = {}
    for key in expected:
        declaration = decode_run_artifact_declaration(runs[key])
        if declaration.member_key != key:
            raise ValueError("calibration run commitment member key differs from its partition key")
        normalized[key] = {
            "member_key": key,
            "root": str(declaration.root),
            "artifacts": dict(declaration.artifact_sha256),
        }
    return normalized


def role_run_commitments(study: StudyDeclaration, commitments: Mapping[str, object], role: str) -> dict[str, object]:
    """Project normalized complete declarations for exactly one study role."""
    if role not in ROLES:
        raise ValueError(f"unsupported calibration commitment role: {role}")
    keys = tuple(member.key for member in study.partitions.members if member.role == role)
    if any(key not in commitments for key in keys):
        raise ValueError("calibration role run commitments are incomplete")
    return {key: commitments[key] for key in keys}


def role_run_artifacts(study: StudyDeclaration, commitments: Mapping[str, object], role: str) -> dict[str, object]:
    """Project only exact artifact hashes from one role's normalized declarations."""
    projected = role_run_commitments(study, commitments, role)
    artifacts: dict[str, object] = {}
    for key, raw in projected.items():
        if not isinstance(raw, Mapping) or not isinstance(raw.get("artifacts"), Mapping):
            raise ValueError("calibration normalized run commitment is malformed")
        artifacts[key] = dict(raw["artifacts"])
    return artifacts
