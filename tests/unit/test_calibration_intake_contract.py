"""Strict normalized calibration intake contracts."""

from __future__ import annotations

import copy
from collections.abc import Callable
from dataclasses import FrozenInstanceError
from hashlib import sha256
from typing import Any

import pytest

from vntyper.scripts.calibration_intake_contract import (
    canonical_intake_bytes,
    decode_intake,
    encode_intake,
)
from vntyper.scripts.calibration_truth import length_target
from vntyper.scripts.canonical_json import load_strict_json_object

pytestmark = pytest.mark.unit

_DIGEST = "a" * 64
_GROUPS = {
    "individual-family": ["individual-family:sample-001"],
    "simulated-pair": ["simulated-pair:sample-001"],
    "backbone-seed-lineage": ["backbone-seed-lineage:sample-001"],
    "replicate-rerun": ["replicate-rerun:sample-001"],
    "depth-series-source": ["depth-series-source:sample-001"],
    "batch": ["batch:sample-001"],
    "repeat-context": ["repeat-context:sample-001"],
}


def synthetic_intake() -> dict[str, Any]:
    """Return one complete synthetic development intake document."""
    return {
        "schema_version": "calibration-intake-v1",
        "specimens": [
            {
                "key": "sample-001",
                "individual_key": "individual-001",
                "family_key": None,
                "identity_status": "confirmed",
                "previously_examined": False,
            }
        ],
        "artifacts": [
            {
                "key": "artifact-001",
                "specimen_key": "sample-001",
                "path": "/synthetic/artifact-001.bam",
                "format": "BAM",
                "mate_path": None,
                "assembly": "synthetic-assembly-v1",
                "assay_class": "synthetic-short-read",
                "input_scope": "full",
                "preprocessing_id": "synthetic-preprocessing-v1",
                "replicate_group": "replicate-001",
                "expected_sha256": _DIGEST,
            }
        ],
        "aliases": [{"alias": "source-sample-001", "specimen_key": "sample-001", "evidence": "explicit-user"}],
        "truth": [
            {
                "specimen_key": "sample-001",
                "genotype": "positive",
                "variants": ["synthetic-variant-001"],
                "length": None,
                "method": "synthetic-method-v1",
                "source_digest": _DIGEST,
                "source_row": 1,
                "status": "confirmed",
            }
        ],
        "assignments": [
            {
                "specimen_key": "sample-001",
                "role": "training",
                "provenance": "development",
                "groups": copy.deepcopy(_GROUPS),
            }
        ],
    }


def test_intake_decodes_immutable_typed_rows_and_preserves_genotype_only_truth() -> None:
    declaration = decode_intake(synthetic_intake())

    assert declaration.specimens[0].key == "sample-001"
    assert declaration.artifacts[0].format == "BAM"
    assert declaration.aliases[0].alias == "source-sample-001"
    assert declaration.truth[0].genotype == "positive"
    assert declaration.truth[0].length is None
    assert declaration.assignments[0].role == "training"
    assert length_target(declaration.truth[0]) is None
    assert len(declaration.sha256) == 64
    with pytest.raises(FrozenInstanceError):
        declaration.specimens[0].key = "changed"  # type: ignore[misc]
    with pytest.raises(TypeError):
        declaration.assignments[0].groups["batch"] = ("changed",)  # type: ignore[index]


def test_intake_encoding_is_canonical_and_independent_of_input_row_order() -> None:
    raw = synthetic_intake()
    second = copy.deepcopy(raw["specimens"][0])
    second["key"] = "sample-002"
    second["individual_key"] = "individual-002"
    raw["specimens"] = [second, raw["specimens"][0]]
    second_artifact = copy.deepcopy(raw["artifacts"][0])
    second_artifact["key"] = "artifact-002"
    second_artifact["specimen_key"] = "sample-002"
    raw["artifacts"] = [second_artifact, raw["artifacts"][0]]
    second_assignment = copy.deepcopy(raw["assignments"][0])
    second_assignment["specimen_key"] = "sample-002"
    second_assignment["groups"] = {name: [f"{name}:sample-002"] for name in _GROUPS}
    raw["assignments"] = [second_assignment, raw["assignments"][0]]

    declaration = decode_intake(raw)
    encoded = encode_intake(declaration)
    specimens = encoded["specimens"]

    assert isinstance(specimens, list)
    assert [row["key"] for row in specimens] == ["sample-001", "sample-002"]
    assert canonical_intake_bytes(declaration).endswith(b"\n")
    assert declaration.sha256 == sha256(canonical_intake_bytes(declaration)).hexdigest()
    assert decode_intake(encoded) == declaration


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda raw: raw.update(schema_version="calibration-intake-v2"), "schema version"),
        (lambda raw: raw.update(specimens=[]), "specimens must be a non-empty"),
        (lambda raw: raw.update(assignments=[]), "exactly one assignment"),
        (lambda raw: raw["specimens"][0].update(identity_status="maybe"), "identity status"),
        (lambda raw: raw["specimens"][0].update(previously_examined=1), "previously examined"),
        (lambda raw: raw["artifacts"][0].update(format="FASTQ"), "artifact format"),
        (lambda raw: raw["artifacts"][0].update(input_scope="unknown"), "input scope"),
        (lambda raw: raw["artifacts"][0].update(expected_sha256="A" * 64), "expected sha256"),
        (lambda raw: raw["aliases"][0].update(evidence="filename"), "alias evidence"),
        (lambda raw: raw["truth"][0].update(genotype="affected"), "truth genotype"),
        (lambda raw: raw["truth"][0].update(status="withheld"), "truth status"),
        (lambda raw: raw["truth"][0].update(source_digest="not-a-digest"), "source digest"),
        (lambda raw: raw["assignments"][0].update(role="test"), "assignment role"),
        (lambda raw: raw["assignments"][0].update(provenance="self-reported"), "assignment provenance"),
        (lambda raw: raw["assignments"][0]["groups"].update(batch=[]), "group batch.*non-empty"),
    ],
)
def test_intake_rejects_invalid_closed_vocabularies_and_required_values(
    mutate: Callable[[dict[str, Any]], None], message: str
) -> None:
    raw = synthetic_intake()
    mutate(raw)

    with pytest.raises(ValueError, match=message):
        decode_intake(raw)


@pytest.mark.parametrize(
    ("collection", "extra_field"),
    [
        (None, "unexpected"),
        ("specimens", "display_name"),
        ("artifacts", "resolved_path"),
        ("aliases", "confidence"),
        ("truth", "eligible"),
        ("assignments", "target"),
    ],
)
def test_intake_rejects_unknown_root_and_row_fields(collection: str | None, extra_field: str) -> None:
    raw = synthetic_intake()
    if collection is None:
        raw[extra_field] = "forbidden"
    else:
        rows = raw[collection]
        assert isinstance(rows, list)
        row = rows[0]
        assert isinstance(row, dict)
        row[extra_field] = "forbidden"

    with pytest.raises(ValueError, match="fields differ"):
        decode_intake(raw)


def test_strict_json_parser_rejects_duplicate_intake_fields_before_typed_decode() -> None:
    raw = '{"schema_version":"calibration-intake-v1","specimens":[],"specimens":[]}'

    with pytest.raises(ValueError, match="duplicate JSON object key: specimens"):
        decode_intake(load_strict_json_object(raw))


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda raw: raw["specimens"][0].update(key=""), "specimen key"),
        (lambda raw: raw["artifacts"][0].update(key=""), "artifact key"),
        (lambda raw: raw["aliases"][0].update(alias=""), "alias"),
        (lambda raw: raw["assignments"][0].update(specimen_key=""), "specimen key"),
        (lambda raw: raw["truth"][0].update(source_row=True), "source row"),
        (lambda raw: raw["truth"][0].update(source_row=float("nan")), "source row"),
    ],
)
def test_intake_rejects_empty_keys_bool_counts_and_nonfinite_numbers(
    mutate: Callable[[dict[str, Any]], None], message: str
) -> None:
    raw = synthetic_intake()
    mutate(raw)

    with pytest.raises(ValueError, match=message):
        decode_intake(raw)


@pytest.mark.parametrize("collection", ["specimens", "artifacts", "assignments"])
def test_intake_rejects_duplicate_primary_row_keys(collection: str) -> None:
    raw = synthetic_intake()
    rows = raw[collection]
    assert isinstance(rows, list)
    rows.append(copy.deepcopy(rows[0]))

    with pytest.raises(ValueError, match="unique"):
        decode_intake(raw)


def test_intake_allows_independent_truth_sources_but_rejects_a_duplicate_source_row() -> None:
    independent = synthetic_intake()
    second = copy.deepcopy(independent["truth"][0])
    second["source_digest"] = "b" * 64
    independent["truth"].append(second)

    assert len(decode_intake(independent).truth) == 2

    duplicate = synthetic_intake()
    duplicate["truth"].append(copy.deepcopy(duplicate["truth"][0]))
    with pytest.raises(ValueError, match="truth source rows.*unique"):
        decode_intake(duplicate)


@pytest.mark.parametrize(
    ("collection", "label"),
    [("artifacts", "artifact"), ("aliases", "alias"), ("truth", "truth")],
)
def test_intake_rejects_rows_joined_to_an_unknown_specimen(collection: str, label: str) -> None:
    raw = synthetic_intake()
    raw[collection][0]["specimen_key"] = "sample-999"

    with pytest.raises(ValueError, match=f"{label}.*unknown specimen"):
        decode_intake(raw)


def test_aliases_are_explicit_and_cannot_conflict_or_shadow_a_specimen_key() -> None:
    raw = synthetic_intake()
    conflicting = copy.deepcopy(raw["specimens"][0])
    conflicting["key"] = "sample-002"
    conflicting["individual_key"] = "individual-002"
    raw["specimens"].append(conflicting)
    raw["aliases"] = [
        {"alias": "source-id", "specimen_key": "sample-001", "evidence": "explicit-user"},
        {"alias": "source-id", "specimen_key": "sample-002", "evidence": "source-crosswalk"},
    ]

    with pytest.raises(ValueError, match="alias.*unique|conflicting alias"):
        decode_intake(raw)

    raw = synthetic_intake()
    raw["aliases"] = [{"alias": "sample-001", "specimen_key": "sample-001", "evidence": "explicit-user"}]
    with pytest.raises(ValueError, match="alias.*specimen key|ambiguous alias"):
        decode_intake(raw)


def test_fastq_pair_requires_an_explicit_mate_and_single_file_formats_forbid_one() -> None:
    pair = synthetic_intake()
    pair["artifacts"][0].update(format="FASTQ_PAIR", mate_path=None)  # type: ignore[index]
    bam = synthetic_intake()
    bam["artifacts"][0].update(mate_path="/synthetic/mate.fastq.gz")  # type: ignore[index]

    for raw in (pair, bam):
        with pytest.raises(ValueError, match="mate"):
            decode_intake(raw)


def test_locked_membership_only_rows_are_allowed_without_artifact_or_truth() -> None:
    raw = synthetic_intake()
    raw["artifacts"] = []
    raw["truth"] = []
    raw["assignments"][0].update(role="locked-heldout", provenance="external-custodian")

    declaration = decode_intake(raw)

    assert declaration.artifacts == ()
    assert declaration.truth == ()
    assert declaration.assignments[0].role == "locked-heldout"


def test_locked_membership_requires_external_custodian_provenance() -> None:
    raw = synthetic_intake()
    raw["artifacts"] = []
    raw["truth"] = []
    raw["assignments"][0]["role"] = "locked-heldout"

    with pytest.raises(ValueError, match="external custodian"):
        decode_intake(raw)


def test_nonlocked_membership_requires_an_artifact_but_truth_remains_optional() -> None:
    no_artifact = synthetic_intake()
    no_artifact["artifacts"] = []
    no_truth = synthetic_intake()
    no_truth["truth"] = []

    with pytest.raises(ValueError, match="requires an input artifact"):
        decode_intake(no_artifact)
    assert decode_intake(no_truth).truth == ()


def test_missing_truth_is_preserved_with_null_length() -> None:
    raw = synthetic_intake()
    raw["truth"][0].update(genotype="unknown", variants=[], status="missing", length=None)

    truth = decode_intake(raw).truth[0]

    assert truth.status == "missing"
    assert truth.length is None


@pytest.mark.parametrize("bearing_collection", ["artifacts", "truth"])
def test_locked_membership_rejects_path_or_truth_bearing_rows(bearing_collection: str) -> None:
    raw = synthetic_intake()
    raw["assignments"][0].update(role="locked-heldout", provenance="external-custodian")
    other_collection = "truth" if bearing_collection == "artifacts" else "artifacts"
    raw[other_collection] = []

    with pytest.raises(ValueError, match="locked.*artifact|locked.*truth"):
        decode_intake(raw)


@pytest.mark.parametrize("role", ["validation", "locked-heldout"])
def test_unresolved_identity_cannot_enter_confirmatory_roles(role: str) -> None:
    raw = synthetic_intake()
    raw["specimens"][0]["identity_status"] = "unresolved"
    raw["assignments"][0]["role"] = role
    if role == "locked-heldout":
        raw["assignments"][0]["provenance"] = "external-custodian"
        raw["artifacts"] = []
        raw["truth"] = []

    with pytest.raises(ValueError, match="unresolved.*validation|unresolved.*locked"):
        decode_intake(raw)


@pytest.mark.parametrize("role", ["validation", "locked-heldout"])
def test_previously_examined_specimens_cannot_enter_confirmatory_roles(role: str) -> None:
    raw = synthetic_intake()
    raw["specimens"][0]["previously_examined"] = True
    raw["assignments"][0]["role"] = role
    if role == "locked-heldout":
        raw["assignments"][0]["provenance"] = "external-custodian"
        raw["artifacts"] = []
        raw["truth"] = []

    with pytest.raises(ValueError, match="previously examined.*validation|previously examined.*locked"):
        decode_intake(raw)
