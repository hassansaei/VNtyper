"""One-use target confirmation around fixed scientific evaluators."""

from __future__ import annotations

import logging
from argparse import Namespace
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

from vntyper.scripts.calibration_artifact_io import write_checksums, write_json
from vntyper.scripts.calibration_candidate import CandidateEnvelope
from vntyper.scripts.calibration_exposure import ExposureReceipt, exposure_receipt_document
from vntyper.scripts.calibration_exposure_io import record_exposure
from vntyper.scripts.calibration_role_source import RoleSource, decode_role_source
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.calibration_target_access import TargetConfirmation, preflight_target_confirmation
from vntyper.scripts.calibration_target_asset_io import read_target_asset
from vntyper.scripts.calibration_target_attestation import (
    TargetLockedAttestation,
    TargetValidationAttestation,
    decode_target_locked_attestation,
    decode_target_validation_attestation,
    encode_target_locked_attestation,
    encode_target_validation_attestation,
    target_locked_attestation_document,
    target_validation_attestation_document,
)
from vntyper.scripts.calibration_target_authority import decode_target_custodian_authority
from vntyper.scripts.calibration_target_completion import encode_target_completion
from vntyper.scripts.calibration_target_contract import TargetStudy
from vntyper.scripts.calibration_target_custody import claim_target_confirmation, preflight_target_custody
from vntyper.scripts.calibration_target_runs import TargetRuns, decode_target_runs
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ConfirmationOutcome:
    """Actual fixed-model metrics, report, and scientific acceptance."""

    metrics: Mapping[str, object]
    report: str
    passed: bool


@dataclass(frozen=True)
class ConfirmationBackend:
    """Verified candidate and fixed evaluator; no fit or selection operation."""

    candidate: CandidateEnvelope
    study: TargetStudy
    evaluate: Callable[[TargetRuns, RoleSource, ExposureReceipt, bytes | None], ConfirmationOutcome]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(path: Path) -> dict[str, object]:
    return load_strict_json_object(read_regular_path(path))


def load_confirmation_backend(target: str, profile_dir: Path) -> ConfirmationBackend:
    """Load a verified research payload and bind its fixed target evaluator.

    Args:
        target: Exactly length or callers.
        profile_dir: Original fit output, including training-only baseline evidence.

    Returns:
        Verified candidate, study, and an evaluator that cannot select a replacement.

    Raises:
        ValueError: If the target, payload, or profile provenance is invalid.
    """
    if target == "length":
        from vntyper.scripts.calibration_length_controller import evaluate_fixed_length_source
        from vntyper.scripts.calibration_length_evaluation import length_evaluation_document
        from vntyper.scripts.calibration_length_profile import load_length_research_profile

        length_profile = load_length_research_profile(profile_dir)

        def evaluate_length(runs, source, receipt, opened_truth):
            result, report = evaluate_fixed_length_source(
                length_profile, runs, source, receipt, opened_truth=opened_truth
            )
            selected = [
                item for item in result.candidates if item.candidate_id == length_profile.selected_protocol_candidate_id
            ]
            if result.phase != source.role or len(selected) != 1:
                _fail("length confirmation did not evaluate its exact fixed candidate and role")
            acceptance = selected[0].acceptance
            return ConfirmationOutcome(
                length_evaluation_document(result), report, acceptance is not None and acceptance.status == "passed"
            )

        return ConfirmationBackend(length_profile.candidate, length_profile.study, evaluate_length)
    if target == "callers":
        from vntyper.scripts.calibration_caller_artifacts import caller_evaluation_document
        from vntyper.scripts.calibration_caller_controller import (
            evaluate_fixed_caller_source,
            load_caller_research_profile,
        )

        caller_profile = load_caller_research_profile(profile_dir)

        def evaluate_callers(runs, source, receipt, opened_truth):
            result, report = evaluate_fixed_caller_source(
                caller_profile, runs, source, receipt, opened_truth=opened_truth
            )
            selected = [
                item for item in result.candidates if item.candidate_id == caller_profile.selected_protocol_candidate_id
            ]
            if result.phase != source.role or len(selected) != 1:
                _fail("caller confirmation did not evaluate its exact fixed candidate and role")
            return ConfirmationOutcome(
                caller_evaluation_document(result), report, selected[0].acceptance.status == "passed"
            )

        return ConfirmationBackend(caller_profile.candidate, caller_profile.study, evaluate_callers)
    _fail("confirmation target must be callers or length")


def _attestation(
    confirmation: TargetConfirmation,
    receipt: ExposureReceipt,
    outcome: ConfirmationOutcome,
    consumption_sha256: str | None,
) -> TargetValidationAttestation | TargetLockedAttestation:
    candidate, study = confirmation.candidate, confirmation.study
    common = {
        "target": candidate.target,
        "status": "passed" if outcome.passed else "failed",
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "study_sha256": study.sha256,
        "protocol_sha256": study.protocol.sha256,
        "partition_sha256": study.partitions.sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "evidence_sha256": confirmation.evidence_sha256,
        "run_manifest_sha256": confirmation.run_manifest_sha256,
        "metrics_sha256": canonical_sha256(outcome.metrics),
        "exposure_ledger_id": study.exposure_ledger_id,
        "exposure_receipt_sha256": receipt.sha256,
    }
    if confirmation.role == "validation":
        return decode_target_validation_attestation(encode_target_validation_attestation(**common))
    if confirmation.validation is None or confirmation.authority is None or consumption_sha256 is None:
        _fail("locked confirmation outcome is missing its prior authorization or consumption")
    return decode_target_locked_attestation(
        encode_target_locked_attestation(
            **common,
            validation_attestation_sha256=confirmation.validation.sha256,
            custodian_authority_sha256=confirmation.authority.sha256,
            custody_consumption_receipt_sha256=consumption_sha256,
        )
    )


def _completion(confirmation: TargetConfirmation, locked: TargetLockedAttestation) -> dict[str, object]:
    validation, authority = confirmation.validation, confirmation.authority
    if validation is None or authority is None or locked.status != "passed":
        _fail("target completion requires passed validation, locked outcome, and external authority")
    return encode_target_completion(
        target=locked.target,
        candidate_sha256=locked.candidate_sha256,
        candidate_id=locked.candidate_id,
        study_sha256=locked.study_sha256,
        protocol_sha256=locked.protocol_sha256,
        partition_sha256=locked.partition_sha256,
        baseline_sha256=locked.baseline_sha256,
        validation_attestation_sha256=validation.sha256,
        locked_heldout_attestation_sha256=locked.sha256,
        custodian_authority_sha256=authority.sha256,
        locked_heldout_evidence_sha256=locked.evidence_sha256,
        locked_heldout_run_manifest_sha256=locked.run_manifest_sha256,
        exposure_ledger_id=locked.exposure_ledger_id,
        validation_exposure_receipt_sha256=validation.exposure_receipt_sha256,
        locked_heldout_exposure_receipt_sha256=locked.exposure_receipt_sha256,
        custody_consumption_receipt_sha256=locked.custody_consumption_receipt_sha256,
    )


def confirm_calibration_bundle(args: Namespace, output: Path, *, role: str) -> bool:
    """Evaluate one authorized role and durably finish its one-use custody claim.

    Args:
        args: Explicit target, research profile, role evidence, ledger, and custody;
            locked evaluation additionally requires validation and authority files.
        output: Empty atomic staging directory owned by the CLI.
        role: Validation or locked-heldout, fixed by the requested operation.

    Returns:
        Actual fixed candidate scientific acceptance. A failed result is still complete.

    Raises:
        ValueError: For incompatible metadata, prior exposure, custody, or integrity.
        OSError: For failed I/O; started or consumed claims remain unavailable for retry.
    """
    if role not in {"validation", "locked-heldout"}:
        _fail("confirmation role must be validation or locked-heldout")
    backend = load_confirmation_backend(args.target, args.profile)
    runs = decode_target_runs(_object(args.evidence / "runs.json"))
    source = decode_role_source(
        _object(args.evidence / "source.json"), study=backend.study, runs=runs, expected_role=role
    )
    validation = authority = None
    if role == "locked-heldout":
        validation = decode_target_validation_attestation(_object(args.validation))
        authority = decode_target_custodian_authority(_object(args.authority))
        if authority.locked_payload_sha256 != source.truth_asset.sha256:
            _fail("locked source truth commitment differs from the external authority")
    confirmation = preflight_target_confirmation(
        candidate=backend.candidate,
        study=backend.study,
        role=role,
        evidence_sha256=source.sha256,
        run_manifest_sha256=runs.sha256,
        validation=validation,
        authority=authority,
    )
    forbidden = (args.profile, args.evidence, output)
    preflight_target_custody(args.custody, confirmation, forbidden_roots=forbidden)
    receipt = record_exposure(
        args.exposure_ledger,
        expected_ledger_id=backend.study.exposure_ledger_id,
        target=args.target,
        role=role,
        study_sha256=backend.study.sha256,
        partition_sha256=backend.study.partitions.sha256,
        evidence_sha256=source.sha256,
        identities=[{"namespace": name, "sha256": digest} for name, digest in source.identities],
        forbidden_roots=forbidden,
    )
    with claim_target_confirmation(args.custody, confirmation, receipt, forbidden_roots=forbidden) as claim:
        opened = None
        if role == "locked-heldout":
            opened = claim.open_locked_payload(lambda: read_target_asset(source.truth_asset))
        outcome = backend.evaluate(runs, source, receipt, None if opened is None else opened.payload)
        if not isinstance(outcome, ConfirmationOutcome) or type(outcome.passed) is not bool or not outcome.report:
            _fail("confirmation evaluator returned an incomplete scientific outcome")
        attestation = _attestation(confirmation, receipt, outcome, None if opened is None else opened.receipt_sha256)
        write_json(output / "metrics.json", outcome.metrics)
        (output / "report.html").write_text(outcome.report, encoding="utf-8")
        write_json(output / "exposure-receipt.json", exposure_receipt_document(receipt))
        if isinstance(attestation, TargetValidationAttestation):
            write_json(output / "validation-attestation.json", target_validation_attestation_document(attestation))
        else:
            write_json(output / "locked-attestation.json", target_locked_attestation_document(attestation))
            if opened is not None:
                write_json(output / "consumption-receipt.json", dict(opened.receipt))
            if outcome.passed:
                write_json(output / "completion.json", _completion(confirmation, attestation))
        write_checksums(output)
        claim.finish(attestation)
    return outcome.passed
