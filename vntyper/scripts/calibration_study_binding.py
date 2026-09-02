"""Pure profile-to-study commitments for generated calibration candidates."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from vntyper.scripts.calibration_manifest import StudyDeclaration
from vntyper.scripts.calibration_run_commitments import ROLES, role_run_artifacts, role_run_commitments
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_MANIFEST_ROLES = ("policy-selection", "training", "validation")
_FIELDS = {
    "schema_version",
    "study_sha256",
    "protocol_sha256",
    "partition_manifest_sha256",
    "objective",
    "seed",
    "role_manifests_sha256",
    "run_commitments_sha256",
    "role_run_commitments_sha256",
    "role_run_artifacts_sha256",
    "role_evidence_sha256",
    "dataset_manifest_sha256",
}


@dataclass(frozen=True)
class StudyBinding:
    """Closed hashes binding one generated candidate to its fitting study."""

    study_sha256: str
    protocol_sha256: str
    partition_manifest_sha256: str
    objective: str
    seed: int
    role_manifests_sha256: Mapping[str, str]
    run_commitments_sha256: str
    role_run_commitments_sha256: Mapping[str, str]
    role_run_artifacts_sha256: Mapping[str, str]
    role_evidence_sha256: Mapping[str, str]
    dataset_manifest_sha256: str
    document: Mapping[str, object]


def build_study_binding(
    study: StudyDeclaration,
    role_manifests: Mapping[str, object],
    run_commitments: Mapping[str, object],
) -> StudyBinding:
    """Build the exact profile dataset commitment for a fitting study.

    Args:
        study: Validated complete study declaration.
        role_manifests: Training and policy-selection role manifests.
        run_commitments: Exact immutable run-artifact hash mappings used to fit.

    Returns:
        A decoded immutable study binding.

    Raises:
        ValueError: If any required role or commitment is missing or malformed.
    """
    if not isinstance(study, StudyDeclaration):
        raise ValueError("calibration study binding requires a StudyDeclaration")
    manifests = _role_hashes(role_manifests)
    if not isinstance(run_commitments, Mapping) or not run_commitments:
        raise ValueError("calibration study binding requires non-empty run commitments")
    run_sha256 = canonical_sha256(run_commitments)
    role_commitments = {role: canonical_sha256(role_run_commitments(study, run_commitments, role)) for role in ROLES}
    role_artifacts = {role: canonical_sha256(role_run_artifacts(study, run_commitments, role)) for role in ROLES}
    role_evidence = {
        role: _role_evidence_hash(role_manifests[role], role_run_artifacts(study, run_commitments, role))
        for role in _MANIFEST_ROLES
    }
    identity = {
        "study_sha256": study.sha256,
        "protocol_sha256": study.protocol.sha256,
        "partition_manifest_sha256": study.partitions.sha256,
        "objective": study.protocol.objective,
        "seed": study.protocol.seed,
        "role_manifests_sha256": manifests,
        "run_commitments_sha256": run_sha256,
        "role_run_commitments_sha256": role_commitments,
        "role_run_artifacts_sha256": role_artifacts,
        "role_evidence_sha256": role_evidence,
    }
    document = {
        "schema_version": "calibration-study-binding-v1",
        **identity,
        "dataset_manifest_sha256": canonical_sha256(identity),
    }
    return decode_study_binding(document)


def decode_study_binding(value: object) -> StudyBinding:
    """Decode a strict study/profile binding document.

    Args:
        value: Candidate JSON object.

    Returns:
        The immutable validated binding.

    Raises:
        ValueError: If fields, types, hashes, or the derived dataset digest differ.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"calibration study binding fields differ: expected {sorted(_FIELDS)}, got {actual}")
    if value["schema_version"] != "calibration-study-binding-v1":
        raise ValueError("calibration study binding schema version is unsupported")
    study_sha256 = _digest(value["study_sha256"], "study")
    protocol_sha256 = _digest(value["protocol_sha256"], "protocol")
    partition_sha256 = _digest(value["partition_manifest_sha256"], "partition manifest")
    objective = value["objective"]
    if objective != "lexicographic-safety-v1":
        raise ValueError("calibration study binding objective must be lexicographic-safety-v1")
    seed = value["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("calibration study binding seed must be a non-negative integer")
    manifests = _decoded_role_hashes(value["role_manifests_sha256"])
    run_sha256 = _digest(value["run_commitments_sha256"], "run commitments")
    role_commitments = _decoded_hashes(value["role_run_commitments_sha256"], ROLES, "role run commitments")
    role_artifacts = _decoded_hashes(value["role_run_artifacts_sha256"], ROLES, "role run artifacts")
    role_evidence = _decoded_hashes(value["role_evidence_sha256"], _MANIFEST_ROLES, "role evidence")
    dataset_sha256 = _digest(value["dataset_manifest_sha256"], "dataset manifest")
    identity = {
        "study_sha256": study_sha256,
        "protocol_sha256": protocol_sha256,
        "partition_manifest_sha256": partition_sha256,
        "objective": objective,
        "seed": seed,
        "role_manifests_sha256": dict(manifests),
        "run_commitments_sha256": run_sha256,
        "role_run_commitments_sha256": dict(role_commitments),
        "role_run_artifacts_sha256": dict(role_artifacts),
        "role_evidence_sha256": dict(role_evidence),
    }
    if canonical_sha256(identity) != dataset_sha256:
        raise ValueError("calibration study binding dataset manifest digest differs")
    document = MappingProxyType(
        {"schema_version": "calibration-study-binding-v1", **identity, "dataset_manifest_sha256": dataset_sha256}
    )
    return StudyBinding(
        study_sha256,
        protocol_sha256,
        partition_sha256,
        objective,
        seed,
        manifests,
        run_sha256,
        role_commitments,
        role_artifacts,
        role_evidence,
        dataset_sha256,
        document,
    )


def validate_study_binding(
    binding: StudyBinding,
    study: StudyDeclaration,
    role_manifests: Mapping[str, object] | None = None,
    run_commitments: Mapping[str, object] | None = None,
) -> None:
    """Verify a binding against its opened study and optional fitting artifacts.

    Args:
        binding: Decoded binding to verify.
        study: Opened complete study declaration.
        role_manifests: Opened fitting role manifests, when their payload is authorized.
        run_commitments: Opened fitting run hashes, when their payload is authorized.

    Raises:
        ValueError: If any supplied study or artifact binding differs.
    """
    if not isinstance(binding, StudyBinding) or not isinstance(study, StudyDeclaration):
        raise ValueError("calibration study binding validation requires decoded values")
    if (
        binding.study_sha256,
        binding.protocol_sha256,
        binding.partition_manifest_sha256,
        binding.objective,
        binding.seed,
    ) != (
        study.sha256,
        study.protocol.sha256,
        study.partitions.sha256,
        study.protocol.objective,
        study.protocol.seed,
    ):
        raise ValueError("calibration study binding differs from the opened study or protocol")
    if role_manifests is not None:
        if not isinstance(role_manifests, Mapping) or not set(role_manifests) <= set(_MANIFEST_ROLES):
            raise ValueError("calibration study binding role manifest commitments differ")
        observed_manifests = {role: canonical_sha256(value) for role, value in role_manifests.items()}
        expected_manifests = {role: binding.role_manifests_sha256[role] for role in role_manifests}
        if observed_manifests != expected_manifests:
            raise ValueError("calibration study binding role manifest commitments differ")
    if run_commitments is not None and (
        not isinstance(run_commitments, Mapping) or canonical_sha256(run_commitments) != binding.run_commitments_sha256
    ):
        raise ValueError("calibration study binding run commitments differ")
    if run_commitments is not None:
        expected_commitments = {
            role: canonical_sha256(role_run_commitments(study, run_commitments, role)) for role in ROLES
        }
        expected_artifacts = {
            role: canonical_sha256(role_run_artifacts(study, run_commitments, role)) for role in ROLES
        }
        if expected_commitments != dict(binding.role_run_commitments_sha256) or expected_artifacts != dict(
            binding.role_run_artifacts_sha256
        ):
            raise ValueError("calibration study binding role run commitments differ")


def validate_role_run_artifacts(binding: StudyBinding, role: str, artifacts: Mapping[str, object]) -> None:
    """Require opened role artifact hashes to equal the predeclared study commitments."""
    if not isinstance(binding, StudyBinding) or role not in ROLES or not isinstance(artifacts, Mapping):
        raise ValueError("calibration role run-artifact validation inputs are invalid")
    normalized = {key: dict(value) if isinstance(value, Mapping) else value for key, value in artifacts.items()}
    if canonical_sha256(normalized) != binding.role_run_artifacts_sha256[role]:
        raise ValueError(f"calibration {role} run artifacts differ from the predeclared study commitment")


def _role_hashes(value: Mapping[str, object]) -> dict[str, str]:
    if not isinstance(value, Mapping) or set(value) != set(_MANIFEST_ROLES):
        raise ValueError(f"calibration study binding requires role manifests {_MANIFEST_ROLES}")
    return {role: canonical_sha256(value[role]) for role in _MANIFEST_ROLES}


def _decoded_role_hashes(value: object) -> Mapping[str, str]:
    return _decoded_hashes(value, _MANIFEST_ROLES, "role manifest")


def _role_evidence_hash(manifest: object, artifacts: Mapping[str, object]) -> str:
    if not isinstance(manifest, Mapping):
        raise ValueError("calibration role manifest must be an object")
    required = {"study_sha256", "features_sha256", "labels_sha256", "baseline_sha256"}
    if not required <= set(manifest):
        raise ValueError("calibration role manifest is incomplete")
    return canonical_sha256(
        {
            "study_sha256": manifest["study_sha256"],
            "features_sha256": manifest["features_sha256"],
            "labels_sha256": manifest["labels_sha256"],
            "baseline_sha256": manifest["baseline_sha256"],
            "run_artifact_sha256": artifacts,
        }
    )


def _decoded_hashes(value: object, keys: tuple[str, ...], label: str) -> Mapping[str, str]:
    if not isinstance(value, Mapping) or set(value) != set(keys):
        raise ValueError(f"calibration study binding {label} hashes must contain {keys}")
    return MappingProxyType({key: _digest(value[key], f"{key} {label}") for key in keys})


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"calibration {label} hash must be lowercase SHA-256")
    return value
