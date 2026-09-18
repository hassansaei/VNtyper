"""Outcome-free metadata for previously examined target assessment evidence."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import NoReturn, cast

from vntyper.scripts.calibration_caller_roster import (
    CallerEligibleRoster,
    caller_eligible_roster_document,
    decode_caller_eligible_roster,
)
from vntyper.scripts.calibration_candidate import CandidateEnvelope, candidate_document
from vntyper.scripts.calibration_exposure import decode_exposure_identities, require_digest
from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    decode_length_eligible_roster,
    length_eligible_roster_document,
)
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
    "target",
    "evidence_role",
    "promotion_eligible",
    "candidate_sha256",
    "study_sha256",
    "run_manifest_sha256",
    "partition_sha256",
    "roster",
    "identities_by_key",
    "evidence_domains",
    "previously_examined",
    "truth_asset",
}
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))


@dataclass(frozen=True)
class DevelopmentSource:
    """Sealed truth and audited identities for one nonpromotable assessment."""

    target: str
    evidence_role: str
    promotion_eligible: bool
    candidate_sha256: str
    study_sha256: str
    run_manifest_sha256: str
    partition_sha256: str
    roster: CallerEligibleRoster | LengthEligibleRoster
    identities_by_key: Mapping[str, tuple[tuple[str, str], ...]]
    evidence_domains: Mapping[str, str]
    previously_examined: Mapping[str, bool]
    truth_asset: TargetRunAsset
    sha256: str

    @property
    def keys(self) -> tuple[str, ...]:
        """Return eligible artifact keys in canonical order."""
        return tuple(sorted(member.key for member in self.roster.members))

    @property
    def identities(self) -> tuple[tuple[str, str], ...]:
        """Return the complete deduplicated exposure identity set."""
        return tuple(sorted({token for values in self.identities_by_key.values() for token in values}))


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"calibration development source {label} fields differ")
    return value


def _roster_document(source: DevelopmentSource) -> list[dict[str, object]]:
    if isinstance(source.roster, CallerEligibleRoster):
        return caller_eligible_roster_document(source.roster)
    if isinstance(source.roster, LengthEligibleRoster):
        return length_eligible_roster_document(source.roster)
    _fail("calibration development source roster is not target-specific")


def _document(source: DevelopmentSource) -> dict[str, object]:
    return {
        "schema_version": "calibration-development-source-v1",
        "target": source.target,
        "evidence_role": source.evidence_role,
        "promotion_eligible": source.promotion_eligible,
        "candidate_sha256": source.candidate_sha256,
        "study_sha256": source.study_sha256,
        "run_manifest_sha256": source.run_manifest_sha256,
        "partition_sha256": source.partition_sha256,
        "roster": _roster_document(source),
        "identities_by_key": {
            key: [{"namespace": namespace, "sha256": digest} for namespace, digest in values]
            for key, values in source.identities_by_key.items()
        },
        "evidence_domains": dict(source.evidence_domains),
        "previously_examined": dict(source.previously_examined),
        "truth_asset": {
            "path": str(source.truth_asset.path),
            "sha256": source.truth_asset.sha256,
            "size_bytes": source.truth_asset.size_bytes,
        },
    }


def decode_development_source(value: object, *, candidate: CandidateEnvelope, runs: TargetRuns) -> DevelopmentSource:
    """Decode assessment metadata without opening truth or result assets.

    Args:
        value: Closed source metadata prepared before assessment.
        candidate: Exact research-only candidate being assessed.
        runs: Exact target run commitments for the development population.

    Returns:
        Immutable source metadata and its aggregate-free content identity.

    Raises:
        ValueError: If membership, prior-exposure status, identities, or bindings differ.
    """
    candidate_document(candidate)
    target_runs_document(runs)
    raw = _object(value, _FIELDS, "root")
    if raw["schema_version"] != "calibration-development-source-v1":
        _fail("calibration development source schema is unsupported")
    if raw["evidence_role"] != "development-assessment":
        _fail("calibration development source requires the development-assessment role")
    if raw["promotion_eligible"] is not False:
        _fail("calibration development source must remain ineligible for promotion")
    target = raw["target"]
    if target != candidate.target or target != runs.target:
        _fail("calibration development source target differs from candidate or runs")
    if raw["candidate_sha256"] != candidate.sha256:
        _fail("calibration development source candidate binding differs")
    if raw["study_sha256"] != candidate.study_sha256:
        _fail("calibration development source study binding differs")
    if raw["run_manifest_sha256"] != runs.sha256:
        _fail("calibration development source run binding differs")
    roster = (
        decode_caller_eligible_roster(raw["roster"])
        if target == "callers"
        else decode_length_eligible_roster(raw["roster"])
    )
    keys = {member.key for member in roster.members}
    if not keys:
        _fail("calibration development source roster must not be empty")
    prior = _object(raw["previously_examined"], keys, "previously examined members")
    if any(value is not True for value in prior.values()):
        _fail("calibration development source requires every member to be previously examined")
    domains = _object(raw["evidence_domains"], keys, "evidence domains")
    if any(value not in {"synthetic", "external"} for value in domains.values()):
        _fail("calibration development source evidence domain is unsupported")
    raw_identities = _object(raw["identities_by_key"], keys, "audited identities")
    identities = {key: decode_exposure_identities(raw_identities[key]) for key in sorted(keys)}
    selected_runs = [run for run in runs.runs if run.manifest_key in keys]
    if {run.manifest_key for run in selected_runs} != keys:
        _fail("calibration development source runs do not cover the exact roster")
    for run in selected_runs:
        if ("physical-readset", run.input_sha256) not in identities[run.manifest_key]:
            _fail("calibration development source physical read identity differs from its run")
    source = DevelopmentSource(
        target,
        "development-assessment",
        False,
        require_digest(raw["candidate_sha256"], "development candidate"),
        require_digest(raw["study_sha256"], "development study"),
        require_digest(raw["run_manifest_sha256"], "development runs"),
        require_digest(raw["partition_sha256"], "development partition"),
        roster,
        MappingProxyType(identities),
        MappingProxyType({key: cast(str, domains[key]) for key in sorted(keys)}),
        MappingProxyType(dict.fromkeys(sorted(keys), True)),
        decode_target_asset(raw["truth_asset"]),
        "",
    )
    return DevelopmentSource(
        source.target,
        source.evidence_role,
        source.promotion_eligible,
        source.candidate_sha256,
        source.study_sha256,
        source.run_manifest_sha256,
        source.partition_sha256,
        source.roster,
        source.identities_by_key,
        source.evidence_domains,
        source.previously_examined,
        source.truth_asset,
        canonical_sha256(_document(source)),
    )


def development_source_document(
    source: DevelopmentSource, *, candidate: CandidateEnvelope, runs: TargetRuns
) -> dict[str, object]:
    """Project and contextually revalidate development source metadata.

    Args:
        source: Previously decoded immutable metadata.
        candidate: Exact bound research candidate.
        runs: Exact bound development run commitments.

    Returns:
        Fresh JSON-compatible closed metadata document.

    Raises:
        ValueError: If typed content, immutability, context, or digest differs.
    """
    if (
        not isinstance(source, DevelopmentSource)
        or not isinstance(source.identities_by_key, _MAPPING_PROXY_TYPE)
        or not isinstance(source.evidence_domains, _MAPPING_PROXY_TYPE)
        or not isinstance(source.previously_examined, _MAPPING_PROXY_TYPE)
    ):
        _fail("calibration development source requires immutable typed content")
    document = _document(source)
    checked = decode_development_source(document, candidate=candidate, runs=runs)
    if checked != source:
        _fail("calibration development source typed content differs from its identity")
    return document
