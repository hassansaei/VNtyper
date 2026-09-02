"""Lineage-bound finalization of one externally held calibration evaluation."""

from __future__ import annotations

import hashlib
from collections.abc import Callable, Mapping
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.calibration_artifact_io import load_object, thaw_json, write_checksums, write_json
from vntyper.scripts.calibration_contract import CandidateMetrics, decode_metrics
from vntyper.scripts.calibration_custodian_import import (
    CustodianAuthority,
    CustodianImport,
    load_custodian_import_header,
    validate_custodian_authority_bindings,
)
from vntyper.scripts.calibration_custody import ensure_candidate_retired
from vntyper.scripts.calibration_locked_artifacts import decode_locked_payload
from vntyper.scripts.calibration_objective import CandidateEvaluation
from vntyper.scripts.calibration_reporting import write_evaluation_artifacts
from vntyper.scripts.calibration_result_artifacts import encode_attestation_hashes, encode_metrics
from vntyper.scripts.calibration_study_binding import validate_role_run_artifacts, validate_study_binding
from vntyper.scripts.calibration_workflow import ExtractedEvidence, evaluate_locked_candidate
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.decision_profile import ResolvedDecisionProfile, resolve_decision_profile


@dataclass(frozen=True)
class _CompletedLockedResult:
    evidence: ExtractedEvidence
    evaluation: CandidateEvaluation
    metrics_document: Mapping[str, object]
    passed: bool


def evaluate_artifact_bundle(
    profile_path: Path,
    evidence_path: Path,
    output: Path,
    *,
    evaluator: Callable[[ResolvedDecisionProfile, ExtractedEvidence], CandidateEvaluation],
    metrics_decoder: Callable[[object], CandidateMetrics] = decode_metrics,
    selector: Callable[..., CandidateEvaluation | None],
) -> bool:
    """Evaluate one profile through a pinned, exact custodian import."""
    profile = resolve_decision_profile(profile_path)
    with load_custodian_import_header(evidence_path) as imported:
        metadata = _profile_metadata(profile)
        _validate_header(profile, metadata, imported)
        authority = imported.authority
        custody = evidence_path.parent / f".{evidence_path.name}.calibration-custody"

        def locked_evaluator(raw: bytes) -> _CompletedLockedResult:
            evidence = decode_locked_payload(raw)
            _validate_locked_evidence(profile, metadata, imported, raw, evidence)
            evaluation = evaluator(profile, evidence)
            metrics_document = encode_metrics(evaluation.metrics)
            passed = selector((evaluation,), evidence.study.protocol) is not None
            return _CompletedLockedResult(evidence, evaluation, metrics_document, passed)

        locked = evaluate_locked_candidate(
            profile,
            imported.read_locked_payload,
            authority.protocol_sha256,
            authority.locked_payload_sha256,
            custody,
            evaluator=locked_evaluator,
        )
        try:
            result = locked.result
            if not isinstance(result, _CompletedLockedResult):
                raise ValueError("locked calibration evaluator returned an invalid result")
            metrics = metrics_decoder(result.metrics_document)
            if not result.passed:
                retirement = ensure_candidate_retired(
                    custody,
                    profile.digest,
                    authority.locked_payload_sha256,
                    "completed-failed-locked-heldout-evaluation",
                )
            else:
                retirement = None
            write_evaluation_artifacts(
                output,
                phase="held-out",
                profile=profile,
                evidence=result.evidence,
                evaluation=result.evaluation,
                accessed_roles=("locked-heldout",),
            )
            output.mkdir(parents=True, exist_ok=True)
            write_json(output / "metrics.json", dict(result.metrics_document))
            write_json(
                output / "attestation.json",
                encode_attestation_hashes(
                    "locked-heldout",
                    profile.digest,
                    authority.protocol_sha256,
                    result.evidence.dataset_sha256,
                    metrics.sha256,
                    result.passed,
                ),
            )
            _write_custody_artifacts(output, authority, locked.receipt.path, retirement)
            write_checksums(output)
            return result.passed
        except BaseException as error:
            with suppress(BaseException):
                ensure_candidate_retired(
                    custody,
                    profile.digest,
                    authority.locked_payload_sha256,
                    f"exceptional-locked-finalization:{type(error).__name__}",
                )
            raise


def _validate_header(
    profile: ResolvedDecisionProfile, metadata: Mapping[str, object], imported: CustodianImport
) -> None:
    authority = imported.authority
    validation = imported.validation
    binding = imported.study_binding
    expected = (
        validation.status,
        validation.profile_sha256,
        validation.protocol_sha256,
        validation.study_sha256,
        validation.partition_manifest_sha256,
        validation.profile_dataset_sha256,
        validation.study_binding_sha256,
        validation.run_commitments_sha256,
        validation.role_run_commitments_sha256,
        validation.role_run_artifacts_sha256,
        validation.evidence_sha256,
    )
    observed = (
        "passed",
        profile.digest,
        binding.protocol_sha256,
        binding.study_sha256,
        binding.partition_manifest_sha256,
        binding.dataset_manifest_sha256,
        canonical_sha256(thaw_json(binding.document)),
        binding.run_commitments_sha256,
        binding.role_run_commitments_sha256["validation"],
        binding.role_run_artifacts_sha256["validation"],
        binding.role_evidence_sha256["validation"],
    )
    if expected != observed:
        raise ValueError("custodian import validation lineage differs from the installed study binding")
    if (
        authority.profile_sha256 != profile.digest
        or authority.profile_dataset_sha256 != binding.dataset_manifest_sha256
        or metadata.get("dataset_manifest_hash") != binding.dataset_manifest_sha256
        or metadata.get("partition_manifest_hash") != binding.partition_manifest_sha256
        or authority.validation_evidence_sha256 != validation.evidence_sha256
        or authority.run_commitments_sha256 != binding.run_commitments_sha256
        or authority.validation_role_run_commitments_sha256 != binding.role_run_commitments_sha256["validation"]
        or authority.validation_role_run_artifacts_sha256 != binding.role_run_artifacts_sha256["validation"]
        or authority.locked_role_run_commitments_sha256 != binding.role_run_commitments_sha256["locked-heldout"]
        or authority.locked_role_run_artifacts_sha256 != binding.role_run_artifacts_sha256["locked-heldout"]
    ):
        raise ValueError("custodian authority, profile, or study-binding lineage differs")


def _validate_locked_evidence(
    profile: ResolvedDecisionProfile,
    metadata: Mapping[str, object],
    imported: CustodianImport,
    raw: bytes,
    evidence: ExtractedEvidence,
) -> None:
    authority = imported.authority
    binding = imported.study_binding
    validation = imported.validation
    validate_study_binding(binding, evidence.study)
    locked_keys = tuple(member.key for member in evidence.study.partitions.members if member.role == "locked-heldout")
    if (
        tuple(row.manifest_key for row in evidence.features.rows) != locked_keys
        or tuple(row.manifest_key for row in evidence.labels.rows) != locked_keys
        or tuple(sorted(evidence.run_hashes)) != locked_keys
    ):
        raise ValueError("custodian import locked artifacts do not match the locked partition members")
    validate_role_run_artifacts(binding, "locked-heldout", evidence.run_hashes)
    run_hashes_sha256 = canonical_sha256({key: dict(value) for key, value in sorted(evidence.run_hashes.items())})
    validate_custodian_authority_bindings(
        authority,
        study_sha256=evidence.study.sha256,
        protocol_sha256=evidence.study.protocol.sha256,
        partition_manifest_sha256=evidence.study.partitions.sha256,
        profile_sha256=profile.digest,
        profile_dataset_sha256=str(metadata["dataset_manifest_hash"]),
        locked_payload_sha256=hashlib.sha256(raw).hexdigest(),
        locked_dataset_sha256=evidence.dataset_sha256,
        locked_run_hashes_sha256=run_hashes_sha256,
        validation_attestation_sha256=validation.sha256,
        validation_evidence_sha256=validation.evidence_sha256,
        study_binding_sha256=canonical_sha256(thaw_json(binding.document)),
        run_commitments_sha256=binding.run_commitments_sha256,
        validation_role_run_commitments_sha256=binding.role_run_commitments_sha256["validation"],
        validation_role_run_artifacts_sha256=binding.role_run_artifacts_sha256["validation"],
        locked_role_run_commitments_sha256=binding.role_run_commitments_sha256["locked-heldout"],
        locked_role_run_artifacts_sha256=binding.role_run_artifacts_sha256["locked-heldout"],
    )
    if (
        metadata.get("objective") != evidence.study.protocol.objective
        or metadata.get("seed") != evidence.study.protocol.seed
    ):
        raise ValueError("generated profile objective or seed differs from the custodian study")


def _write_custody_artifacts(
    output: Path, authority: CustodianAuthority, receipt_path: Path, retirement: Path | None
) -> None:
    write_json(
        output / "custody_limitations.json",
        {
            "schema_version": "calibration-custody-limitations-v1",
            "local_verification_proves_external_independence": False,
            "statement": (
                "This local verifier checks a supplied named external-custodian attestation and exact bindings; "
                "without a separately configured cryptographic trust root it does not prove external independence."
            ),
            "custodian_name": authority.custodian_name,
            "attestation_id": authority.attestation_id,
        },
    )
    write_json(
        output / "access.json", {"schema_version": "calibration-access-v1", "accessed_roles": ["locked-heldout"]}
    )
    write_json(output / "custody_receipt.json", load_object(receipt_path, "custody receipt"))
    if retirement is not None:
        write_json(output / "retirement.json", load_object(retirement, "calibration retirement"))


def _profile_metadata(profile: ResolvedDecisionProfile) -> Mapping[str, object]:
    metadata = profile.document.get("generated_metadata")
    if not isinstance(metadata, Mapping):
        raise ValueError("custodian evaluation requires generated profile metadata")
    return metadata
