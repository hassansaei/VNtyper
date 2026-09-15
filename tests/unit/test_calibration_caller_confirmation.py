"""Invented complete caller evidence exercises real fit and one-use confirmation."""

from argparse import Namespace
from pathlib import Path

import pytest

from tests.builders import kestrel_config
from tests.unit.test_calibration_caller_protocol import protocol_document
from tests.unit.test_calibration_kestrel_replay import _GG, _candidate, _capture, _raw
from tests.unit.test_calibration_manifest import _manifest, _member
from tests.unit.test_calibration_target_contract import study_document
from tests.unit.test_calibration_target_runs import run_document
from vntyper.scripts.calibration_artifact_io import load_object, write_json
from vntyper.scripts.calibration_caller_policy import caller_policy_values_document
from vntyper.scripts.calibration_candidate import candidate_producer_document
from vntyper.scripts.calibration_exposure_io import initialize_exposure_ledger
from vntyper.scripts.calibration_kestrel_capture import kestrel_capture_document
from vntyper.scripts.calibration_manifest import connected_leakage_groups
from vntyper.scripts.calibration_target_contract import decode_target_study
from vntyper.scripts.calibration_target_custody import initialize_target_custody
from vntyper.scripts.calibration_target_runs import decode_target_runs
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit

IDENTITY = "MUC1-X-60-coding-v1|60|59|-|C"
NEGATIVE_TSV = (
    "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\t"
    "Estimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n"
    "None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n"
)


def _asset(path):
    import hashlib

    raw = path.read_bytes()
    return {"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def caller_fit_fixture(tmp_path, *, failed_role=None):
    """Two weak positives and one empty control per independent invented role."""
    root = tmp_path / "evidence"
    root.mkdir()
    ledger = tmp_path / "exposure.jsonl"
    ledger_id = initialize_exposure_ledger(ledger)
    config = kestrel_config(
        **{
            "confidence_assignment.reporting_floor": 0.003,
            "confidence_assignment.depth_score_thresholds.low": 0.003,
        }
    )
    positive = _capture(_raw(depth_alt=4, depth_region=1000), config)
    negative = _capture(_raw().iloc[0:0], config)
    baseline = positive.baseline_policy
    candidate = _candidate(baseline, **{_GG: 0.003})
    raw_protocol = protocol_document(baseline, [(candidate, 1)])
    raw_protocol["required_strata"] = ["all"]
    raw_protocol["acceptance"].update(
        minimum_positive_groups=2,
        minimum_negative_groups=1,
        maximum_fpr_upper=1.0,
        minimum_sensitivity_delta_lower=-1.0,
    )
    roles = {
        f"{prefix}-{kind}": role
        for prefix, role in (("t", "training"), ("s", "policy-selection"), ("v", "validation"), ("l", "locked-heldout"))
        for kind in ("p1", "p2", "n")
    }
    members = [_member(key, role, assay_class="capture") for key, role in sorted(roles.items())]
    for member in members:
        if member["role"] == "locked-heldout":
            member["provenance"] = "external-custodian"
    raw_study = study_document("callers")
    raw_study.update(protocol=raw_protocol, partitions=_manifest(*members), exposure_ledger_id=ledger_id)
    raw_study["applicability"]["required_callers"] = ["kestrel"]
    raw_study["baseline"]["policy"] = caller_policy_values_document(baseline)
    study = decode_target_study(raw_study)
    groups = connected_leakage_groups(study.partitions)
    write_json(root / "study.json", raw_study)
    native = root / "baseline.tsv"
    native.write_text(NEGATIVE_TSV)
    raw_runs = {"schema_version": "calibration-runs-v2", "target": "callers", "runs": []}
    for key in sorted(roles):
        capture = negative if key.endswith("-n") else positive
        capture_path = root / (key + "-capture.json")
        write_json(capture_path, kestrel_capture_document(capture))
        for policy in sorted((baseline, candidate), key=lambda item: item.sha256):
            row = run_document(target="callers")["runs"][0]
            row.update(
                manifest_key=key,
                policy_sha256=policy.sha256,
                input_sha256=canonical_sha256(key),
                producer_sha256=canonical_sha256(candidate_producer_document(study.baseline.producer)),
                baseline_assets_sha256=study.baseline.assets_sha256,
                capture_policy_sha256=capture.provenance.capture_policy_sha256,
                execution_kind="baseline-rerun" if policy == baseline else "scalar-replay",
                assets={"kestrel_capture": _asset(capture_path), "kestrel_result": _asset(native)},
            )
            raw_runs["runs"].append(row)
    runs = decode_target_runs(raw_runs)
    write_json(root / "runs.json", raw_runs)
    for role in sorted(set(roles.values())):
        keys = sorted(key for key in roles if roles[key] == role)
        truth_path = root / "roles" / role / "truth.json"
        truth_rows = [
            {
                "key": key,
                "genotype": "negative" if key.endswith("-n") else "positive",
                "variants": [] if key.endswith("-n") else [IDENTITY],
            }
            for key in keys
        ]
        if failed_role == role:
            # Independent truth has no carriers: confirmation must be insufficient, not a pass.
            truth_rows = [{"key": key, "genotype": "negative", "variants": []} for key in keys]
        write_json(truth_path, {"schema_version": "calibration-caller-truth-v1", "rows": truth_rows})
        source = {
            "schema_version": "calibration-role-source-v2",
            "study_sha256": study.sha256,
            "role": role,
            "roster": [
                {"key": key, "group_key": groups[key], "strata": ["all"]}
                for key in sorted(keys, key=lambda item: groups[item])
            ],
            "excluded": [],
            "evidence_domains": dict.fromkeys(keys, "synthetic"),
            "identities_by_key": {
                key: [
                    {"namespace": "physical-readset", "sha256": canonical_sha256(key)},
                    {"namespace": "specimen", "sha256": canonical_sha256("specimen-" + key)},
                ]
                for key in keys
            },
            "truth_asset": _asset(truth_path),
            "run_manifest_sha256": runs.sha256,
        }
        write_json(root / "roles" / role / "source.json", source)
        write_json(root / "roles" / role / "runs.json", raw_runs)
    return Namespace(evidence=root, exposure_ledger=ledger, objective="caller-safety-v1"), study


@pytest.mark.parametrize("passed", [True, False])
def test_real_caller_fit_and_fixed_validation(tmp_path, passed):
    from vntyper.scripts.calibration_caller_controller import fit_caller_bundle, load_caller_research_profile
    from vntyper.scripts.calibration_confirmation_controller import confirm_calibration_bundle

    fit_args, study = caller_fit_fixture(tmp_path, failed_role=None if passed else "validation")
    profile = tmp_path / "profile"
    profile.mkdir()
    assert fit_caller_bundle(fit_args, profile) is True
    research = load_caller_research_profile(profile)
    selected = study.protocol.candidates[0].candidate_id
    assert research.selected_protocol_candidate_id == selected
    custody = tmp_path / "custody"
    initialize_target_custody(custody, exposure_ledger_id=study.exposure_ledger_id)
    args = Namespace(
        target="callers",
        profile=profile,
        evidence=fit_args.evidence / "roles" / "validation",
        exposure_ledger=fit_args.exposure_ledger,
        custody=custody,
    )
    output = tmp_path / "validation-output"
    output.mkdir()
    assert confirm_calibration_bundle(args, output, role="validation") is passed
    attestation = load_object(output / "validation-attestation.json", "test")
    assert attestation["status"] == ("passed" if passed else "failed")
    metrics = load_object(output / "metrics.json", "test")
    assert metrics["phase"] == "validation"
    assert len(metrics["candidates"]) == 1
    assert metrics["candidates"][0]["candidate_id"] == selected
    assert metrics["selection"]["status"] == "not-applicable"
    with pytest.raises(ValueError):
        confirm_calibration_bundle(args, output, role="validation")


@pytest.mark.parametrize("passed", [True, False])
def test_real_locked_caller_truth_is_consumed_once_before_evaluation(tmp_path, monkeypatch, passed):
    from vntyper.scripts import calibration_caller_controller as callers
    from vntyper.scripts import calibration_confirmation_controller as confirmation
    from vntyper.scripts.calibration_target_authority import encode_target_custodian_authority

    fit_args, study = caller_fit_fixture(tmp_path, failed_role=None if passed else "locked-heldout")
    profile = tmp_path / "profile"
    profile.mkdir()
    assert callers.fit_caller_bundle(fit_args, profile)
    custody = tmp_path / "custody"
    initialize_target_custody(custody, exposure_ledger_id=study.exposure_ledger_id)
    args = Namespace(
        target="callers",
        profile=profile,
        evidence=fit_args.evidence / "roles" / "validation",
        exposure_ledger=fit_args.exposure_ledger,
        custody=custody,
    )
    validation_output = tmp_path / "validation-output"
    validation_output.mkdir()
    assert confirmation.confirm_calibration_bundle(args, validation_output, role="validation")
    args.validation = validation_output / "validation-attestation.json"
    validation = load_object(args.validation, "test")
    args.evidence = fit_args.evidence / "roles" / "locked-heldout"
    source = load_object(args.evidence / "source.json", "test")
    args.authority = tmp_path / "authority.json"
    common = {
        key: validation[key]
        for key in (
            "target",
            "candidate_sha256",
            "candidate_id",
            "study_sha256",
            "protocol_sha256",
            "partition_sha256",
            "baseline_sha256",
            "exposure_ledger_id",
        )
    }
    write_json(
        args.authority,
        encode_target_custodian_authority(
            **common,
            validation_attestation_sha256=canonical_sha256(validation),
            validation_evidence_sha256=validation["evidence_sha256"],
            validation_run_manifest_sha256=validation["run_manifest_sha256"],
            validation_exposure_receipt_sha256=validation["exposure_receipt_sha256"],
            locked_heldout_evidence_sha256=canonical_sha256(source),
            locked_heldout_run_manifest_sha256=source["run_manifest_sha256"],
            locked_payload_sha256=source["truth_asset"]["sha256"],
            custodian_name="invented independent operator",
            attestation_id="invented approval",
        ),
    )
    truth_path = Path(source["truth_asset"]["path"])
    reads = []
    original_read = confirmation.read_target_asset
    original_json = callers.read_target_json

    def consume(asset):
        assert asset.path == truth_path
        assert list(custody.glob("*.locked-consumption.json"))
        reads.append(asset.path)
        return original_read(asset)

    def forbid_reopen(asset):
        assert asset.path != truth_path, "locked truth reopened by fixed evaluator"
        return original_json(asset)

    monkeypatch.setattr(confirmation, "read_target_asset", consume)
    monkeypatch.setattr(callers, "read_target_json", forbid_reopen)
    output = tmp_path / "locked-output"
    output.mkdir()
    assert confirmation.confirm_calibration_bundle(args, output, role="locked-heldout") is passed
    assert reads == [truth_path]
    locked = load_object(output / "locked-attestation.json", "test")
    assert locked["status"] == ("passed" if passed else "failed")
    assert (output / "completion.json").exists() is passed
    assert (output / "consumption-receipt.json").is_file()
    if passed:
        from vntyper.scripts.calibration_caller_bundle import load_caller_model_bundle
        from vntyper.scripts.calibration_export import export_calibration_bundle

        exported = tmp_path / "approved"
        exported.mkdir()
        args.evaluation = output / "locked-attestation.json"
        args.completion = output / "completion.json"
        assert export_calibration_bundle(args, exported)
        bundle = load_caller_model_bundle(exported)
        assert bundle.candidate.sha256 == callers.load_caller_research_profile(profile).candidate.sha256
        assert bundle.advntr_policy is None
        assert bundle.profile.canonical_bytes == (profile / "payload" / "decision-profile.json").read_bytes()
    with pytest.raises(ValueError):
        confirmation.confirm_calibration_bundle(args, output, role="locked-heldout")
    assert reads == [truth_path]
