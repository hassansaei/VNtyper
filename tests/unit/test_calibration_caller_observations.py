from __future__ import annotations

import hashlib
import json
from copy import deepcopy

import pytest

from tests.builders import kestrel_config
from tests.unit.test_calibration_kestrel_replay import _capture, _raw
from vntyper.scripts.calibration_caller_observations import (
    caller_truth_document,
    decode_advntr_baseline_calls,
    decode_advntr_replay_calls,
    decode_caller_truth,
    native_caller_observation,
    replayed_caller_observation,
)
from vntyper.scripts.calibration_caller_roster import EligibleCallerMember

pytestmark = pytest.mark.unit


def _truth() -> dict[str, object]:
    return {
        "schema_version": "calibration-caller-truth-v1",
        "rows": [
            {"key": "case-known", "genotype": "positive", "variants": ["identity-a"]},
            {"key": "case-missing", "genotype": "positive", "variants": None},
            {"key": "control", "genotype": "negative", "variants": []},
            {"key": "unknown", "genotype": "unknown", "variants": None},
        ],
    }


def test_caller_truth_preserves_known_missing_negative_and_unknown_identity() -> None:
    truth = decode_caller_truth(_truth(), ("case-known", "case-missing", "control", "unknown"))

    assert truth.rows[0].variants == ("identity-a",)
    assert truth.rows[1].variants is None
    assert truth.rows[2].variants == ()
    assert truth.rows[3].genotype == "unknown"
    assert caller_truth_document(truth) == _truth()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda raw: raw.update(extra=True), "fields"),
        (lambda raw: raw["rows"].append(deepcopy(raw["rows"][0])), "duplicate"),
        (lambda raw: raw["rows"][0].update(genotype=True), "genotype"),
        (lambda raw: raw["rows"][2].update(variants=None), "negative"),
        (lambda raw: raw["rows"][3].update(variants=[]), "unknown"),
    ],
)
def test_caller_truth_rejects_closed_schema_and_semantic_hazards(mutation, message: str) -> None:
    raw = _truth()
    mutation(raw)

    with pytest.raises(ValueError, match=message):
        decode_caller_truth(raw, ("case-known", "case-missing", "control", "unknown"))


def test_native_negative_tsv_uses_capture_disposition_without_guessing() -> None:
    truth = decode_caller_truth(_truth(), ("case-known", "case-missing", "control", "unknown"))
    negative = (
        b"Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\t"
        b"Estimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n"
        b"None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n"
    )

    observation = native_caller_observation(
        EligibleCallerMember("control", "group-control", ("all",)),
        truth,
        kestrel_tsv=negative,
        advntr_tsv=None,
        disposition="zero-candidate",
        source_evidence_sha256="a" * 64,
    )

    assert observation.observation.called_positive is False
    assert observation.observation.called_variants == ()
    assert observation.disposition == "zero-candidate"


def test_native_tsv_rejects_malformed_rows_and_call_disposition_mismatch() -> None:
    truth = decode_caller_truth(_truth(), ("case-known", "case-missing", "control", "unknown"))
    member = EligibleCallerMember("control", "group-control", ("all",))
    malformed = b"Motif\tConfidence\nNone\tNegative\textra\n"

    with pytest.raises(ValueError, match="malformed TSV"):
        native_caller_observation(
            member,
            truth,
            kestrel_tsv=malformed,
            advntr_tsv=None,
            disposition="zero-candidate",
            source_evidence_sha256="a" * 64,
        )
    positive = b"Motifs\tPOS\tREF\tALT\nX-X\t67\tG\tGG\n"
    with pytest.raises(ValueError, match="projection|identity"):
        native_caller_observation(
            member,
            truth,
            kestrel_tsv=positive,
            advntr_tsv=None,
            disposition="no-call",
            source_evidence_sha256="a" * 64,
        )


def test_replayed_observation_uses_complete_kestrel_capture_and_retains_no_call() -> None:
    truth = decode_caller_truth(_truth(), ("case-known", "case-missing", "control", "unknown"))
    config = kestrel_config()
    capture = _capture(_raw(depth_alt=7, depth_region=500), config)
    member = EligibleCallerMember("case-known", "group-case", ("all",))

    called = replayed_caller_observation(
        member,
        truth,
        kestrel_capture=capture,
        policy=capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
        advntr_calls=(),
        advntr_assessable=True,
        source_evidence_sha256="b" * 64,
    )
    strict = dict(capture.baseline_policy.values)
    strict["/components/kestrel/confidence_assignment/reporting_floor"] = 1
    from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values

    no_call_policy = decode_caller_policy_values(
        {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": ["kestrel"],
            "values": strict,
        }
    )
    no_call = replayed_caller_observation(
        member,
        truth,
        kestrel_capture=capture,
        policy=no_call_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
        advntr_calls=(),
        advntr_assessable=True,
        source_evidence_sha256="c" * 64,
    )

    assert called.observation.called_positive is True
    assert called.observation.called_variants
    assert no_call.observation.called_positive is None
    assert no_call.disposition == "no-call"


def test_replayed_observation_distinguishes_zero_candidates_from_audit_failure() -> None:
    truth = decode_caller_truth(_truth(), ("case-known", "case-missing", "control", "unknown"))
    capture = _capture(_raw().iloc[0:0], kestrel_config())
    member = EligibleCallerMember("control", "group-control", ("all",))

    row = replayed_caller_observation(
        member,
        truth,
        kestrel_capture=capture,
        policy=capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
        advntr_calls=(),
        advntr_assessable=True,
        source_evidence_sha256="d" * 64,
    )

    assert row.observation.called_positive is False
    assert row.disposition == "zero-candidate"

    failed = replayed_caller_observation(
        member,
        truth,
        kestrel_capture=capture,
        policy=capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
        advntr_calls=(),
        advntr_assessable=False,
        source_evidence_sha256="e" * 64,
    )
    assert failed.observation.called_positive is None
    assert failed.disposition == "unsupported"


def _json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode() + b"\n"


def _advntr_record() -> dict[str, object]:
    return {
        "schema_version": "advntr-frameshift-capture-v2",
        "completion": "completed-vntr",
        "producer": {"package_version": "2.4.0", "build_id": "a" * 64, "source_revision": "b" * 40},
        "assets": {},
        "model_locus": {},
        "loaded_background": None,
        "capture_policy": {},
        "caller_policy": {},
        "locus": {"vntr_id": 17},
        "unit_geometry": {},
        "reference_order": [],
        "flank_boundaries": {},
        "decision_visits": [
            {
                "plan": {"state": "I22_2_G_LEN1", "read_support": 8},
                "mean_coverage": 20.0,
                "statistic": {"called": True, "pvalue": 0.01},
            }
        ],
        "warnings": [],
        "occurrences": [],
        "spans": [],
        "evidence_rows": [],
        "candidate_traversal": {},
    }


def test_advntr_baseline_calls_preserve_native_decision_receipts_and_audit() -> None:
    record = _advntr_record()
    calls, assessable = decode_advntr_baseline_calls(_json(record), (17,))
    assert calls == ({"state": "I22_2_G_LEN1", "read_support": 8, "mean_coverage": 20.0, "pvalue": 0.01},)
    assert assessable is True

    record["warnings"] = [
        {"origin": "calibration-audit", "ordinal": None, "state": "invented", "disposition": "attribution-outside-trials"}
    ]
    assert decode_advntr_baseline_calls(_json(record), (17,))[1] is False


def test_advntr_replay_output_binds_manifest_policy_capture_and_complete_roster() -> None:
    record_raw = _json(_advntr_record())
    record_line = record_raw.rstrip(b"\n")
    record_sha = hashlib.sha256(record_line).hexdigest()
    manifest = {
        "schema_version": "advntr-frameshift-replay-manifest-v1",
        "captures": [{"key": "case-known", "filename": "capture.jsonl", "sha256": hashlib.sha256(record_raw).hexdigest(), "vntr_ids": [17]}],
    }
    policy = {"schema_version": "advntr-frameshift-replay-policy-v1", "capture_policy": {}, "caller_policy": {}}
    locus = {
        "schema_version": "advntr-frameshift-replay-result-v1",
        "vntr_id": 17,
        "capture_record_sha256": record_sha,
        "policy_sha256": hashlib.sha256(_json(policy).rstrip(b"\n")).hexdigest(),
        "capture_producer": _advntr_record()["producer"],
        "capture_assets": {},
        "loaded_background_sha256": None,
        "baseline_parity": True,
        "decision_visits": [],
        "calls": [{"state": "I22_2_G_LEN1", "read_support": 9, "mean_coverage": 20.0, "pvalue": 0.005}],
        "warnings": [],
        "capture_audit": {"attribution_outside_trials": [], "calibrated_policy_domain_errors": []},
    }
    manifest_raw, policy_raw = _json(manifest), _json(policy)
    output = {
        "schema_version": "advntr-frameshift-replay-output-v1",
        "manifest_file_sha256": hashlib.sha256(manifest_raw).hexdigest(),
        "manifest_sha256": hashlib.sha256(manifest_raw.rstrip(b"\n")).hexdigest(),
        "policy_file_sha256": hashlib.sha256(policy_raw).hexdigest(),
        "policy_sha256": hashlib.sha256(policy_raw.rstrip(b"\n")).hexdigest(),
        "background_file_sha256": None,
        "replay_producer": _advntr_record()["producer"],
        "results": [{"key": "case-known", "capture_sha256": hashlib.sha256(record_raw).hexdigest(), "vntrs": [{"vntr_id": 17, "result": locus}]}],
    }

    calls, assessable = decode_advntr_replay_calls(
        _json(output), manifest_raw=manifest_raw, policy_raw=policy_raw, capture_raw=record_raw,
        expected_key="case-known", expected_vntr_ids=(17,), expected_background_sha256=None,
    )
    assert calls[0]["read_support"] == 9
    assert assessable is True

    output["manifest_file_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="binding"):
        decode_advntr_replay_calls(
            _json(output), manifest_raw=manifest_raw, policy_raw=policy_raw, capture_raw=record_raw,
            expected_key="case-known", expected_vntr_ids=(17,), expected_background_sha256=None,
        )
