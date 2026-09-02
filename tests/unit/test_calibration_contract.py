"""Strict calibration protocol and artifact contracts."""

from copy import deepcopy
from typing import cast

import pytest

from vntyper.scripts.calibration_contract import (
    decode_attestation,
    decode_evidence_manifest,
    decode_metrics,
    decode_policy,
    decode_protocol,
    validate_attestation_bindings,
    validate_evidence_bindings,
)

pytestmark = pytest.mark.unit


def synthetic_protocol() -> dict[str, object]:
    """Return the finite test-only protocol declared by the implementation plan."""
    return {
        "objective": "lexicographic-safety-v1",
        "bootstrap_iterations": 10000,
        "bootstrap_interval": "percentile",
        "multiplicity_method": "holm",
        "seed": 295,
        "maximum_free_parameters": 4,
        "minimum_stratum_count": 2,
        "maximum_abstention_fraction": 0.25,
        "assay_classes": ["capture-short-read"],
        "mutation_classes": ["duplication"],
        "candidate_grid": {
            "minimum_record_count_margin": [1, 2],
            "minimum_record_share": [0.5, 0.75],
            "minimum_record_share_margin": [0.0, 0.25],
            "xd_veto": ["disabled", "missingness"],
        },
    }


def test_protocol_decodes_exact_finite_grid_and_hash() -> None:
    protocol = decode_protocol(synthetic_protocol())

    assert protocol.objective == "lexicographic-safety-v1"
    assert protocol.bootstrap_iterations == 10_000
    assert protocol.maximum_abstention_fraction.numerator == 1
    assert protocol.maximum_abstention_fraction.denominator == 4
    assert protocol.assay_classes == ("capture-short-read",)
    assert protocol.mutation_classes == ("duplication",)
    assert protocol.required_strata == ("capture-short-read:duplication",)
    assert len(protocol.candidates) == 16
    assert len(protocol.sha256) == 64
    assert protocol.candidates[0].minimum_record_count_margin == 1
    assert protocol.candidates[-1].xd_veto == "missingness"


@pytest.mark.parametrize("field", ["assay_classes", "mutation_classes"])
def test_protocol_requires_predeclared_study_dimensions(field: str) -> None:
    raw = synthetic_protocol()
    del raw[field]

    with pytest.raises(ValueError, match=field.replace("_", " ")):
        decode_protocol(raw)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("assay_classes", []),
        ("assay_classes", ("capture-short-read",)),
        ("assay_classes", [""]),
        ("assay_classes", [1]),
        ("assay_classes", ["genome-short-read", "capture-short-read"]),
        ("assay_classes", ["capture-short-read", "capture-short-read"]),
        ("mutation_classes", []),
        ("mutation_classes", ("duplication",)),
        ("mutation_classes", [""]),
        ("mutation_classes", [1]),
        ("mutation_classes", ["duplication", "deletion"]),
        ("mutation_classes", ["duplication", "duplication"]),
    ],
)
def test_protocol_study_dimensions_are_sorted_unique_and_non_empty(field: str, value: object) -> None:
    raw = synthetic_protocol()
    raw[field] = value

    with pytest.raises(ValueError, match="sorted unique non-empty strings"):
        decode_protocol(raw)


def test_protocol_digest_and_cross_product_bind_every_predeclared_stratum() -> None:
    raw = synthetic_protocol()
    raw["assay_classes"] = ["capture-short-read", "genome-short-read"]
    raw["mutation_classes"] = ["deletion", "duplication"]

    protocol = decode_protocol(raw)

    assert protocol.required_strata == (
        "capture-short-read:deletion",
        "capture-short-read:duplication",
        "genome-short-read:deletion",
        "genome-short-read:duplication",
    )
    changed = deepcopy(raw)
    changed["mutation_classes"] = ["duplication"]
    assert decode_protocol(changed).sha256 != protocol.sha256


@pytest.mark.parametrize(
    ("assay_classes", "mutation_classes"),
    [
        (["capture:short-read"], ["duplication"]),
        (["capture-short-read"], ["dup:lication"]),
        (["a", "a:b"], ["b:c", "c"]),
    ],
)
def test_protocol_rejects_ambiguous_stratum_delimiters(
    assay_classes: list[str],
    mutation_classes: list[str],
) -> None:
    raw = synthetic_protocol()
    raw["assay_classes"] = assay_classes
    raw["mutation_classes"] = mutation_classes

    with pytest.raises(ValueError, match="must not contain ':'"):
        decode_protocol(raw)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("objective",), "f1"),
        (("bootstrap_iterations",), 9999),
        (("bootstrap_interval",), "basic"),
        (("multiplicity_method",), "none"),
        (("seed",), -1),
        (("maximum_free_parameters",), -1),
        (("minimum_stratum_count",), 0),
        (("maximum_abstention_fraction",), 1.1),
        (("candidate_grid", "minimum_record_share"), []),
        (("candidate_grid", "xd_veto"), ["winner"]),
        (("candidate_grid", "xd_veto"), ["concentration"]),
        (("candidate_grid", "xd_veto"), ["discordance"]),
    ],
)
def test_protocol_rejects_values_outside_closed_v1_contract(path: tuple[str, ...], value: object) -> None:
    raw = synthetic_protocol()
    target = raw
    for key in path[:-1]:
        target = target[key]  # type: ignore[assignment,index]
    target[path[-1]] = value

    with pytest.raises(ValueError, match="protocol|candidate"):
        decode_protocol(raw)


def test_protocol_rejects_unknown_or_missing_fields_and_duplicate_grid_values() -> None:
    unknown = synthetic_protocol()
    unknown["extra"] = True
    missing = synthetic_protocol()
    del missing["seed"]
    duplicate = synthetic_protocol()
    duplicate["candidate_grid"]["minimum_record_share"] = [0.5, 0.5]  # type: ignore[index]

    for raw in (unknown, missing, duplicate):
        with pytest.raises(ValueError):
            decode_protocol(raw)


def test_policy_evidence_metrics_and_attestation_are_strict_and_hash_bound() -> None:
    digest = "a" * 64
    policy = {
        "schema_version": "calibration-policy-v1",
        "profile_sha256": digest,
        "dominance": {
            "enabled": True,
            "minimum_record_count_margin": 1,
            "minimum_record_share": 0.5,
            "minimum_record_share_margin": 0.0,
            "xd_veto": "disabled",
            "abstain_on_inadmissible_advntr": False,
        },
    }
    evidence = {
        "schema_version": "calibration-evidence-v1",
        "role": "policy-selection",
        "provenance": "development",
        "protocol_sha256": digest,
        "partition_manifest_sha256": "b" * 64,
        "features_sha256": "c" * 64,
        "labels_sha256": "d" * 64,
        "baseline_sha256": "e" * 64,
    }
    metrics = {
        "schema_version": "calibration-metrics-v1",
        "candidate_profile_sha256": digest,
        "wrong_tier_a_displayed_names": 0,
        "control_findings": 0,
        "wrong_displayed_names_all_tiers": 1,
        "macro_exact_recovery": "3/4",
        "binary_detection_sensitivity": "4/5",
        "free_parameter_count": 1,
        "abstention_fraction": "1/10",
        "tier_a_reachable": True,
        "applicability_matches": True,
    }
    attestation = {
        "schema_version": "calibration-attestation-v1",
        "role": "validation",
        "status": "passed",
        "profile_sha256": digest,
        "protocol_sha256": "b" * 64,
        "evidence_sha256": "c" * 64,
        "metrics_sha256": "d" * 64,
    }

    assert decode_policy(policy).profile_sha256 == digest
    assert decode_evidence_manifest(evidence).role == "policy-selection"
    assert decode_metrics(metrics).macro_exact_recovery.numerator == 3
    assert decode_attestation(attestation).status == "passed"

    decoder_cases = (
        (decode_policy, policy),
        (decode_evidence_manifest, evidence),
        (decode_metrics, metrics),
        (decode_attestation, attestation),
    )
    for decoder, raw in decoder_cases:
        changed = deepcopy(cast(dict[str, object], raw))
        changed["unexpected"] = True
        with pytest.raises(ValueError):
            decoder(changed)


def test_development_evidence_cannot_be_relabeled_as_locked_heldout() -> None:
    evidence = {
        "schema_version": "calibration-evidence-v1",
        "role": "locked-heldout",
        "provenance": "development",
        "protocol_sha256": "a" * 64,
        "partition_manifest_sha256": "b" * 64,
        "features_sha256": "c" * 64,
        "labels_sha256": "d" * 64,
        "baseline_sha256": "e" * 64,
    }

    with pytest.raises(ValueError, match="held-out|custodian"):
        decode_evidence_manifest(evidence)


def test_evidence_and_attestation_bind_exact_opened_hashes_and_roles() -> None:
    evidence = decode_evidence_manifest(
        {
            "schema_version": "calibration-evidence-v1",
            "role": "validation",
            "provenance": "development",
            "protocol_sha256": "a" * 64,
            "partition_manifest_sha256": "b" * 64,
            "features_sha256": "c" * 64,
            "labels_sha256": "d" * 64,
            "baseline_sha256": "e" * 64,
        }
    )
    attestation = decode_attestation(
        {
            "schema_version": "calibration-attestation-v1",
            "role": "validation",
            "status": "passed",
            "profile_sha256": "f" * 64,
            "protocol_sha256": "a" * 64,
            "evidence_sha256": evidence.sha256,
            "metrics_sha256": "1" * 64,
        }
    )

    validate_evidence_bindings(
        evidence,
        protocol_sha256="a" * 64,
        partition_manifest_sha256="b" * 64,
        features_sha256="c" * 64,
        labels_sha256="d" * 64,
        baseline_sha256="e" * 64,
    )
    validate_attestation_bindings(
        attestation,
        role="validation",
        profile_sha256="f" * 64,
        protocol_sha256="a" * 64,
        evidence_sha256=evidence.sha256,
        metrics_sha256="1" * 64,
    )
    with pytest.raises(ValueError, match="hash"):
        validate_evidence_bindings(
            evidence,
            protocol_sha256="0" * 64,
            partition_manifest_sha256="b" * 64,
            features_sha256="c" * 64,
            labels_sha256="d" * 64,
            baseline_sha256="e" * 64,
        )
    with pytest.raises(ValueError, match="role|hash"):
        validate_attestation_bindings(
            attestation,
            role="locked-heldout",
            profile_sha256="f" * 64,
            protocol_sha256="a" * 64,
            evidence_sha256=evidence.sha256,
            metrics_sha256="1" * 64,
        )
