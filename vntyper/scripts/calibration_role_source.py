"""Outcome-free metadata for separately sealed target calibration roles.

The trusted intake operator supplies complete identity tokens and freezes truth
eligibility before caller or length outcomes are observed. This contract checks
those declarations against the study and raw input commitments; it does not
infer laboratory identity or read a combined truth table.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_caller_roster import (
    CallerEligibleRoster,
    caller_eligible_roster_document,
    decode_caller_eligible_roster,
)
from vntyper.scripts.calibration_exposure import decode_exposure_identities, require_digest
from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    decode_length_eligible_roster,
    length_eligible_roster_document,
)
from vntyper.scripts.calibration_manifest import connected_leakage_groups
from vntyper.scripts.calibration_target_contract import TargetStudy, target_study_document
from vntyper.scripts.calibration_target_runs import (
    TargetRunAsset,
    TargetRuns,
    decode_target_asset,
    target_runs_document,
)
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)
_FIELDS = {
    "schema_version",
    "study_sha256",
    "role",
    "roster",
    "excluded",
    "evidence_domains",
    "identities_by_key",
    "truth_asset",
    "run_manifest_sha256",
}
_ROLES = {"training", "policy-selection", "validation", "locked-heldout"}
_EXCLUSIONS = {"truth-unavailable", "identity-quarantined", "nonrepresentative"}


@dataclass(frozen=True)
class RoleSource:
    """Role metadata and sealed truth commitment, with no opened outcomes."""

    study_sha256: str
    role: str
    roster: CallerEligibleRoster | LengthEligibleRoster
    excluded: tuple[tuple[str, str], ...]
    evidence_domains: Mapping[str, str]
    identities_by_key: Mapping[str, tuple[tuple[str, str], ...]]
    truth_asset: TargetRunAsset
    run_manifest_sha256: str
    sha256: str

    @property
    def keys(self) -> tuple[str, ...]:
        """Return eligible artifact keys in canonical artifact order."""
        return tuple(sorted(member.key for member in self.roster.members))

    @property
    def identities(self) -> tuple[tuple[str, str], ...]:
        """Return all identities exposed when this role's truth asset is read."""
        return tuple(sorted({token for values in self.identities_by_key.values() for token in values}))


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"calibration role source {label} fields differ")
    return value


def _exclusions(value: object) -> tuple[tuple[str, str], ...]:
    if not isinstance(value, list):
        _fail("calibration role source exclusions must be a list")
    rows = []
    for row in value:
        raw = _object(row, {"key", "reason"}, "exclusion")
        key, reason = raw["key"], raw["reason"]
        if not isinstance(key, str) or not key or not isinstance(reason, str) or reason not in _EXCLUSIONS:
            _fail("calibration role source exclusion is not a pre-outcome eligibility reason")
        rows.append((key, reason))
    if rows != sorted(rows) or len({key for key, _ in rows}) != len(rows):
        _fail("calibration role source exclusions must be unique and sorted")
    return tuple(rows)


def _identities(value: object, all_keys: set[str], runs: TargetRuns) -> Mapping[str, tuple[tuple[str, str], ...]]:
    raw = _object(value, all_keys, "audited identities")
    decoded = {key: decode_exposure_identities(raw[key]) for key in sorted(raw)}
    for run in runs.runs:
        if run.manifest_key in all_keys and ("physical-readset", run.input_sha256) not in decoded[run.manifest_key]:
            _fail("calibration role source physical read identity differs from its committed run")
    return MappingProxyType(decoded)


def _document(source: RoleSource) -> dict[str, object]:
    roster = (
        caller_eligible_roster_document(source.roster)
        if isinstance(source.roster, CallerEligibleRoster)
        else length_eligible_roster_document(source.roster)
    )
    return {
        "schema_version": "calibration-role-source-v2",
        "study_sha256": source.study_sha256,
        "role": source.role,
        "roster": roster,
        "excluded": [{"key": key, "reason": reason} for key, reason in source.excluded],
        "evidence_domains": dict(source.evidence_domains),
        "identities_by_key": {
            key: [{"namespace": namespace, "sha256": digest} for namespace, digest in tokens]
            for key, tokens in source.identities_by_key.items()
        },
        "truth_asset": {
            "path": str(source.truth_asset.path),
            "sha256": source.truth_asset.sha256,
            "size_bytes": source.truth_asset.size_bytes,
        },
        "run_manifest_sha256": source.run_manifest_sha256,
    }


def decode_role_source(value: object, *, study: TargetStudy, runs: TargetRuns, expected_role: str) -> RoleSource:
    """Check metadata against one study role before reading any source outcomes.

    Args:
        value: Closed role-source metadata, supplied by the intake operator.
        study: Frozen complete study with all four role memberships.
        runs: Exact local run commitments, which may span several roles.
        expected_role: Role authorized by the caller, never inferred from the file.

    Returns:
        Immutable eligible roster, exclusions, identity tokens and truth commitment.

    Raises:
        ValueError: For changed membership, group, target, input or other bindings.
    """
    target_study_document(study)
    target_runs_document(runs)
    raw = _object(value, _FIELDS, "root")
    if raw["schema_version"] != "calibration-role-source-v2":
        _fail("calibration role source schema is unsupported")
    if expected_role not in _ROLES or raw["role"] != expected_role:
        _fail("calibration role source is outside the authorized role")
    if raw["study_sha256"] != study.sha256 or raw["run_manifest_sha256"] != runs.sha256 or study.target != runs.target:
        _fail("calibration role source study, runs or target binding differs")
    roster = (
        decode_caller_eligible_roster(raw["roster"])
        if study.target == "callers"
        else decode_length_eligible_roster(raw["roster"])
    )
    excluded = _exclusions(raw["excluded"])
    included_keys = {member.key for member in roster.members}
    excluded_keys = {key for key, _ in excluded}
    expected_keys = {member.key for member in study.partitions.members if member.role == expected_role}
    if included_keys & excluded_keys or included_keys | excluded_keys != expected_keys:
        _fail("calibration role source must account for every declared role member exactly once")
    groups = connected_leakage_groups(study.partitions)
    if any(groups[member.key] != member.group_key for member in roster.members):
        _fail("calibration role source group identity differs from its study partition")
    domains = _object(raw["evidence_domains"], included_keys, "evidence domains")
    domain_values = {}
    for key, domain in domains.items():
        if not isinstance(domain, str) or domain not in {"synthetic", "external"}:
            _fail("calibration role source evidence domain is unsupported")
        domain_values[key] = domain
    identities = _identities(raw["identities_by_key"], expected_keys, runs)
    source = RoleSource(
        require_digest(raw["study_sha256"], "role study"),
        expected_role,
        roster,
        excluded,
        MappingProxyType(domain_values),
        identities,
        decode_target_asset(raw["truth_asset"]),
        require_digest(raw["run_manifest_sha256"], "role runs"),
        "",
    )
    return RoleSource(
        source.study_sha256,
        source.role,
        source.roster,
        source.excluded,
        source.evidence_domains,
        source.identities_by_key,
        source.truth_asset,
        source.run_manifest_sha256,
        canonical_sha256(_document(source)),
    )


def role_source_document(source: RoleSource, *, study: TargetStudy, runs: TargetRuns) -> dict[str, object]:
    """Project verified metadata after rechecking its context and identity.

    Args:
        source: Previously decoded role metadata.
        study: Exact original study.
        runs: Exact original run manifest.

    Returns:
        Fresh local JSON-compatible metadata, including its sealed truth path.

    Raises:
        ValueError: If the typed content or its contextual binding changed.
    """
    if not isinstance(source, RoleSource):
        _fail("calibration role source projection requires typed metadata")
    document = _document(source)
    checked = decode_role_source(document, study=study, runs=runs, expected_role=source.role)
    if checked != source:
        _fail("calibration role source typed content differs from its identity")
    return document
