"""Trusted curator projection from normalized intake to sealed target truth."""

from __future__ import annotations

import hashlib
import logging
import os
from collections.abc import Mapping
from pathlib import Path
from typing import NoReturn

from vntyper.scripts.calibration_caller_observations import decode_caller_truth
from vntyper.scripts.calibration_identity import IdentityAudit, build_partition_members, identity_audit_document
from vntyper.scripts.calibration_intake_contract import IntakeDeclaration, TruthRecord, encode_intake
from vntyper.scripts.calibration_length_metrics import decode_length_eligible_roster
from vntyper.scripts.calibration_manifest import connected_leakage_groups
from vntyper.scripts.calibration_role_source import RoleSource, decode_role_source, role_source_document
from vntyper.scripts.calibration_target_contract import TargetStudy, target_study_document
from vntyper.scripts.calibration_target_runs import TargetRunAsset, TargetRuns, target_runs_document
from vntyper.scripts.calibration_truth import ConversionRegistry, converted_length_target, length_target
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION

logger = logging.getLogger(__name__)
_ROLES = {"training", "policy-selection", "validation"}
_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _private_json(path: Path, value: object) -> bytes:
    raw = canonical_json_bytes(value)
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    os.chmod(path.parent, 0o700)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | _CLOEXEC | _NOFOLLOW, 0o600)
    try:
        os.fchmod(descriptor, 0o600)
        view = memoryview(raw)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("target curation write made no progress")
            view = view[written:]
    finally:
        os.close(descriptor)
    return raw


def _validate_staging(
    staging: Path,
    role: str,
    declaration: IntakeDeclaration,
    audit: IdentityAudit,
    study: TargetStudy,
    runs: TargetRuns,
) -> None:
    if not staging.is_dir() or staging.is_symlink():
        _fail("target curation staging must be an existing regular directory")
    names = {path.name for path in staging.iterdir()}
    if names - {"roles", "curation"}:
        _fail("target curation staging inventory contains unrelated files")
    roles_root = staging / "roles"
    curation_root = staging / "curation"
    existing_roles: set[str] = set()
    if roles_root.exists():
        if not roles_root.is_dir() or roles_root.is_symlink():
            _fail("target curation role staging is not a regular directory")
        existing_roles = {path.name for path in roles_root.iterdir()}
        if existing_roles - _ROLES or role in existing_roles:
            _fail("target curation role was already prepared or is unknown")
        for existing in existing_roles:
            root = roles_root / existing
            if (
                root.is_symlink()
                or not root.is_dir()
                or {path.name for path in root.iterdir()}
                != {
                    "source.json",
                    "truth.json",
                }
            ):
                _fail("target curation existing role inventory differs")
    if curation_root.exists():
        if not curation_root.is_dir() or curation_root.is_symlink():
            _fail("target curation audit staging is not a regular directory")
        if {path.name for path in curation_root.iterdir()} != {f"{existing}.json" for existing in existing_roles}:
            _fail("target curation audit inventory differs from prepared roles")
    elif existing_roles:
        _fail("target curation audit inventory differs from prepared roles")
    specimens = {specimen.key: specimen for specimen in declaration.specimens}
    artifacts = {artifact.key: artifact for artifact in declaration.artifacts}
    for existing in existing_roles:
        try:
            curation = load_strict_json_object((curation_root / f"{existing}.json").read_bytes())
            source = decode_role_source(
                load_strict_json_object((roles_root / existing / "source.json").read_bytes()),
                study=study,
                runs=runs,
                expected_role=existing,
            )
        except (OSError, ValueError) as error:
            raise ValueError("target curation existing role content is invalid") from error
        fields = {
            "schema_version": "calibration-target-curation-audit-v1",
            "target": study.target,
            "role": existing,
            "intake_sha256": declaration.sha256,
            "identity_audit_sha256": audit.sha256,
            "study_sha256": study.sha256,
            "run_manifest_sha256": runs.sha256,
            "source_sha256": source.sha256,
        }
        previous = {
            member.key: specimens[artifacts[member.key].specimen_key].previously_examined
            for member in study.partitions.members
            if member.role == existing
        }
        if (
            set(curation) != {*fields, "previously_examined"}
            or any(curation[name] != value for name, value in fields.items())
            or curation["previously_examined"] != previous
        ):
            _fail("target curation existing role differs from the current trust context")


def _members(
    declaration: IntakeDeclaration, audit: IdentityAudit, study: TargetStudy, role: str
) -> tuple[tuple[str, ...], dict[str, str]]:
    encode_intake(declaration)
    identity_audit_document(declaration, audit)
    target_study_document(study)
    projected = {member.key: member for member in build_partition_members(declaration, audit)}
    role_members = tuple(member for member in study.partitions.members if member.role == role)
    if not role_members or any(
        member.key not in projected or projected[member.key] != member for member in role_members
    ):
        _fail("target curation role membership differs from the audited intake")
    specimen_by_artifact = {artifact.key: artifact.specimen_key for artifact in declaration.artifacts}
    return tuple(member.key for member in role_members), specimen_by_artifact


def _truth_by_specimen(declaration: IntakeDeclaration) -> dict[str, tuple[TruthRecord, ...]]:
    grouped: dict[str, list[TruthRecord]] = {}
    for row in declaration.truth:
        grouped.setdefault(row.specimen_key, []).append(row)
    return {key: tuple(values) for key, values in grouped.items()}


def _caller_value(rows: tuple[TruthRecord, ...]) -> tuple[str, list[str] | None] | None:
    values: set[tuple[str, tuple[str, ...] | None]] = set()
    for row in rows:
        if row.status != "confirmed":
            continue
        variants: tuple[str, ...] | None = row.variants
        if row.genotype in {"positive", "unknown"} and not variants:
            variants = None
        if row.genotype == "negative" and variants:
            _fail("target curation negative truth cannot carry variant identities")
        if row.genotype == "unknown" and variants:
            _fail("target curation unknown truth cannot carry variant identities")
        values.add((row.genotype, variants))
    if not values:
        return None
    if len(values) != 1:
        _fail("target curation caller truth differs across confirmed records")
    genotype, variants = values.pop()
    return genotype, None if variants is None else list(variants)


def _length_value(rows: tuple[TruthRecord, ...], registry: ConversionRegistry | None) -> int | None:
    values: set[int] = set()
    for row in rows:
        if row.status != "confirmed" or row.length is None or row.length.measurement != "exact":
            continue
        length = row.length
        if (
            length.unit == "repeat-count"
            and length.conversion_id is None
            and length.boundary_definition == TARGET_BOUNDARY_DEFINITION
        ):
            value = length_target(row)
        else:
            if registry is None:
                _fail("target curation length truth requires an explicit conversion registry")
            value = converted_length_target(row, registry, target_boundary_definition=TARGET_BOUNDARY_DEFINITION)
        if value is None or not value.is_integer():
            _fail("target curation length truth must produce an exact integral total")
        values.add(int(value))
    if not values:
        return None
    if len(values) != 1:
        _fail("target curation length truth differs across confirmed records")
    return values.pop()


def _identities(
    keys: tuple[str, ...], declaration: IntakeDeclaration, audit: IdentityAudit, runs: TargetRuns
) -> dict[str, list[dict[str, str]]]:
    artifacts = {artifact.key: artifact for artifact in declaration.artifacts}
    specimens = {specimen.key: specimen for specimen in declaration.specimens}
    inputs = {run.manifest_key: run.input_sha256 for run in runs.runs}
    result = {}
    for key in keys:
        artifact = artifacts[key]
        specimen = specimens[artifact.specimen_key]
        logical = audit.fingerprints[key].logical
        if key not in inputs:
            _fail("target curation role member is absent from the run manifest")
        tokens = {
            ("specimen", canonical_sha256(specimen.individual_key or specimen.key)),
            ("named-readset", logical.named_sequence_sha256),
            ("unnamed-readset", logical.unnamed_sequence_sha256),
            ("physical-readset", inputs[key]),
        }
        if specimen.family_key is not None:
            tokens.add(("family", canonical_sha256(specimen.family_key)))
        result[key] = [{"namespace": namespace, "sha256": digest} for namespace, digest in sorted(tokens)]
    return result


def _source_documents(
    declaration: IntakeDeclaration,
    audit: IdentityAudit,
    study: TargetStudy,
    runs: TargetRuns,
    role: str,
    truth_path: Path,
    evidence_domains: Mapping[str, str],
    strata_by_key: Mapping[str, tuple[str, ...]],
    registry: ConversionRegistry | None,
) -> tuple[RoleSource, dict[str, object], dict[str, object]]:
    role_keys, specimen_by_artifact = _members(declaration, audit, study, role)
    truth_rows = _truth_by_specimen(declaration)
    representatives = set(audit.primary_artifact_by_group.values())
    included: list[str] = []
    excluded: list[dict[str, str]] = []
    values: dict[str, object] = {}
    for key in role_keys:
        specimen = specimen_by_artifact[key]
        if specimen in audit.quarantined_specimens:
            excluded.append({"key": key, "reason": "identity-quarantined"})
            continue
        if key not in representatives:
            excluded.append({"key": key, "reason": "nonrepresentative"})
            continue
        value = (
            _caller_value(truth_rows.get(specimen, ()))
            if study.target == "callers"
            else _length_value(truth_rows.get(specimen, ()), registry)
        )
        if value is None:
            excluded.append({"key": key, "reason": "truth-unavailable"})
            continue
        included.append(key)
        values[key] = value
    if not included:
        _fail("target curation role has no eligible truth-bearing representative")
    if not isinstance(evidence_domains, Mapping) or set(evidence_domains) != set(included):
        _fail("target curation evidence domains must match the eligible roster exactly")
    if not isinstance(strata_by_key, Mapping) or set(strata_by_key) != set(included):
        _fail("target curation strata must match the eligible roster exactly")
    if any(not isinstance(strata_by_key[key], tuple) for key in included):
        _fail("target curation strata values must be immutable tuples")
    groups = connected_leakage_groups(study.partitions)
    roster = [
        {"key": key, "group_key": groups[key], "strata": list(strata_by_key[key])}
        for key in sorted(included, key=lambda value: groups[value])
    ]
    if study.target == "length":
        decode_length_eligible_roster(roster)
        truth: dict[str, object] = {
            "schema_version": "calibration-length-truth-v1",
            "boundary_definition": TARGET_BOUNDARY_DEFINITION,
            "rows": [{"key": key, "total_repeat_count": values[key]} for key in sorted(included)],
        }
    else:
        truth = {
            "schema_version": "calibration-caller-truth-v1",
            "rows": [
                {"key": key, "genotype": values[key][0], "variants": values[key][1]}  # type: ignore[index]
                for key in sorted(included)
            ],
        }
        decode_caller_truth(truth, tuple(sorted(included)))
    truth_raw = canonical_json_bytes(truth)
    asset = TargetRunAsset(truth_path, hashlib.sha256(truth_raw).hexdigest(), len(truth_raw))
    source_raw: dict[str, object] = {
        "schema_version": "calibration-role-source-v2",
        "study_sha256": study.sha256,
        "role": role,
        "roster": roster,
        "excluded": sorted(excluded, key=lambda row: row["key"]),
        "evidence_domains": {key: evidence_domains[key] for key in sorted(included)},
        "identities_by_key": _identities(role_keys, declaration, audit, runs),
        "truth_asset": {"path": str(asset.path), "sha256": asset.sha256, "size_bytes": asset.size_bytes},
        "run_manifest_sha256": runs.sha256,
    }
    source = decode_role_source(source_raw, study=study, runs=runs, expected_role=role)
    return source, truth, source_raw


def prepare_target_role_source(
    declaration: IntakeDeclaration,
    audit: IdentityAudit,
    study: TargetStudy,
    runs: TargetRuns,
    staging: Path,
    *,
    installed_root: Path,
    role: str,
    evidence_domains: Mapping[str, str],
    strata_by_key: Mapping[str, tuple[str, ...]],
    conversion_registry: ConversionRegistry | None = None,
) -> RoleSource:
    """Split one authorized nonlocked intake role into sealed target truth.

    Args:
        declaration: Already-opened normalized intake at the trusted curator stage.
        audit: Recomputed identity audit bound to that intake.
        study: Frozen target study declared before target outcomes.
        runs: Exact target run commitments.
        staging: Empty private staging tree used for output writes.
        installed_root: Final absolute tree path encoded in the truth commitment.
        role: Training, policy-selection, or validation role to prepare.
        evidence_domains: Explicit synthetic/external domain for each eligible key.
        strata_by_key: Explicit overlapping strata for each eligible key.
        conversion_registry: Required only when length truth needs conversion.

    Returns:
        Strict role source whose truth asset points into ``installed_root``.

    Raises:
        ValueError: If trust bindings, eligibility, truth, conversion, or paths differ.
    """
    if role not in _ROLES:
        _fail("target curation cannot prepare locked or unknown roles from normalized intake")
    if (
        not isinstance(staging, Path)
        or not isinstance(installed_root, Path)
        or not installed_root.is_absolute()
        or installed_root != Path(os.path.normpath(installed_root))
        or installed_root.name in {"", ".", ".."}
    ):
        _fail("target curation requires staging and a normalized absolute installed root")
    _validate_staging(staging, role, declaration, audit, study, runs)
    target_runs_document(runs)
    if study.target != runs.target:
        _fail("target curation study and runs target differ")
    truth_path = installed_root / "roles" / role / "truth.json"
    source, truth, source_raw = _source_documents(
        declaration,
        audit,
        study,
        runs,
        role,
        truth_path,
        evidence_domains,
        strata_by_key,
        conversion_registry,
    )
    specimens = {specimen.key: specimen for specimen in declaration.specimens}
    artifacts = {artifact.key: artifact for artifact in declaration.artifacts}
    previous = {
        key: specimens[artifacts[key].specimen_key].previously_examined
        for key in sorted(member.key for member in study.partitions.members if member.role == role)
    }
    curation = {
        "schema_version": "calibration-target-curation-audit-v1",
        "target": study.target,
        "role": role,
        "intake_sha256": declaration.sha256,
        "identity_audit_sha256": audit.sha256,
        "study_sha256": study.sha256,
        "run_manifest_sha256": runs.sha256,
        "previously_examined": previous,
        "source_sha256": source.sha256,
    }
    _private_json(staging / "roles" / role / "truth.json", truth)
    _private_json(staging / "roles" / role / "source.json", source_raw)
    _private_json(staging / "curation" / f"{role}.json", curation)
    role_source_document(source, study=study, runs=runs)
    return source
