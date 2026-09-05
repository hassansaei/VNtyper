"""Decision-profile snapshot and summary-schema provenance contracts."""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)
from vntyper.scripts.profile_provenance import (
    LEGACY_DECISION_PROFILE_REVISION,
    DecisionProfileProvenance,
    profile_summary_fields,
    resolve_summary_profile,
    snapshot_decision_profile,
    verify_profile_snapshot,
)
from vntyper.scripts.summary import start_summary
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_KESTREL

pytestmark = pytest.mark.unit

PROFILE_FIELDS = {
    "decision_profile_id",
    "decision_profile_revision",
    "decision_profile_kind",
    "decision_profile_source",
    "decision_profile_sha256",
    "decision_profile_snapshot",
}

GOOD_IDENTITY_ROW = {
    "Motif": "5",
    "Confidence": "High_Precision",
    "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|C",
    "Molecular_Identity_Status": "unique",
    "Equivalent_Representation_Count": 1,
    "Identity_Hypothesis_Count": 1,
}


def _schema_three_summary(profile, *, rows: list[dict[str, object]] | None = None) -> dict[str, object]:
    return {
        "schema_version": 3,
        "decision_policy": "legacy-selection-v1",
        "advntr_evidence_digest": None,
        **profile_summary_fields(profile),
        "steps": [] if rows is None else [{"step": STEP_KESTREL, "parsed_result": {"comments": [], "data": rows}}],
    }


def _snapshot_packaged(tmp_path: Path):
    profile = load_packaged_decision_profile()
    snapshot = tmp_path / "provenance" / "decision_profile.json"
    snapshot_decision_profile(profile, snapshot)
    return profile, snapshot


def _record_advntr_evidence_digest(summary: dict[str, object], profile: ResolvedDecisionProfile) -> None:
    advntr = profile.components["advntr"]
    assert isinstance(advntr, Mapping)
    artifact_evidence = advntr["artifact_evidence"]
    assert isinstance(artifact_evidence, Mapping)
    summary["advntr_evidence_digest"] = artifact_evidence["digest"]


def test_profile_summary_fields_and_snapshot_are_exact_and_path_free(tmp_path: Path) -> None:
    profile, snapshot = _snapshot_packaged(tmp_path)

    assert profile_summary_fields(profile) == {
        "decision_profile_id": profile.profile_id,
        "decision_profile_revision": profile.profile_revision,
        "decision_profile_kind": "packaged",
        "decision_profile_source": "package",
        "decision_profile_sha256": profile.digest,
        "decision_profile_snapshot": "provenance/decision_profile.json",
    }
    assert snapshot.read_bytes() == profile.canonical_bytes
    assert not Path(profile_summary_fields(profile)["decision_profile_snapshot"]).is_absolute()


def test_profile_projection_and_snapshot_reject_unresolved_objects(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="resolved profile"):
        profile_summary_fields(object())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="resolved profile"):
        snapshot_decision_profile(object(), tmp_path / "profile.json")  # type: ignore[arg-type]


def test_current_summary_is_schema_three_and_records_the_resolved_profile() -> None:
    profile = load_packaged_decision_profile()

    current = start_summary(decision_profile=profile)

    assert current["schema_version"] == 3
    assert current["decision_policy"] == "legacy-selection-v1"
    assert current["advntr_evidence_digest"] is None
    assert {key: current[key] for key in PROFILE_FIELDS} == profile_summary_fields(profile)


def test_schema_three_snapshot_verifies_to_the_recorded_profile(tmp_path: Path) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)

    recorded = resolve_summary_profile(_schema_three_summary(profile, rows=[GOOD_IDENTITY_ROW]), tmp_path)

    assert recorded.profile_id == profile.profile_id
    assert recorded.profile_revision == profile.profile_revision
    assert recorded.profile_kind == profile.profile_kind
    assert recorded.source == profile.source
    assert recorded.sha256 == profile.digest
    assert recorded.snapshot_path == "provenance/decision_profile.json"
    assert recorded.revision == profile.profile_revision


@pytest.mark.parametrize("schema_version", [None, 1, 2])
def test_legacy_summary_never_reads_or_infers_a_current_profile(tmp_path: Path, schema_version: int | None) -> None:
    summary: dict[str, object] = {"steps": []}
    if schema_version is not None:
        summary["schema_version"] = schema_version

    recorded = resolve_summary_profile(summary, tmp_path / "does-not-exist")

    assert recorded.profile_id is None
    assert recorded.sha256 is None
    assert recorded.revision == LEGACY_DECISION_PROFILE_REVISION


@pytest.mark.parametrize("schema_version", [1, 2])
def test_legacy_schema_rejects_schema_three_profile_fields_as_a_downgrade(tmp_path: Path, schema_version: int) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = {"schema_version": schema_version, "steps": [], **profile_summary_fields(profile)}

    with pytest.raises(ValueError, match="legacy.*decision profile|schema 3.*field"):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize("missing", ["advntr_evidence_digest", *sorted(PROFILE_FIELDS)])
def test_schema_three_requires_every_profile_and_pr_b_evidence_field(tmp_path: Path, missing: str) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    del summary[missing]

    with pytest.raises(ValueError, match="schema 3|required"):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("decision_profile_snapshot", "/operator/private/profile.json", "relative snapshot"),
        ("decision_profile_sha256", "0" * 64, "digest mismatch"),
        ("decision_profile_kind", "explicit-custom", "kind|source|metadata"),
        ("decision_profile_source", "/operator/private/profile.json", "source"),
    ],
)
def test_schema_three_rejects_leaks_mismatches_and_invalid_relationships(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    summary[field] = value

    with pytest.raises(ValueError, match=message):
        resolve_summary_profile(summary, tmp_path)


def test_schema_three_rejects_noncanonical_snapshot_bytes(tmp_path: Path) -> None:
    profile, snapshot = _snapshot_packaged(tmp_path)
    snapshot.write_text(json.dumps(json.loads(profile.canonical_bytes), indent=2), encoding="utf-8")

    with pytest.raises(ValueError, match="canonical"):
        resolve_summary_profile(_schema_three_summary(profile), tmp_path)


def test_package_source_requires_the_checked_in_packaged_hash(tmp_path: Path) -> None:
    profile, snapshot = _snapshot_packaged(tmp_path)
    document = json.loads(profile.canonical_bytes)
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    inventory["/components/kestrel/duplicate_flagging/flag_name"]["value"] = "Forged_Package_Default"
    forged = canonical_json_bytes(document)
    snapshot.write_bytes(forged)
    summary = _schema_three_summary(profile)
    summary["decision_profile_sha256"] = hashlib.sha256(forged).hexdigest()

    with pytest.raises(ValueError, match="checked-in packaged profile"):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("decision_profile_id", "", "non-empty"),
        ("decision_profile_kind", "foreign", "kind is unsupported"),
        ("decision_profile_source", "explicit-cli", "explicit-cli source"),
        ("decision_profile_sha256", "0" * 63, "SHA-256"),
        ("advntr_evidence_digest", "x" * 64, "adVNTR evidence digest"),
        ("decision_policy", None, "decision_policy"),
        ("decision_policy", "", "decision_policy"),
    ],
)
def test_schema_three_rejects_invalid_top_level_values(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    summary[field] = value

    with pytest.raises(ValueError, match=message):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize("summary", [[], {"schema_version": True}, {"schema_version": 4}])
def test_summary_profile_rejects_nonobjects_and_unknown_schema_versions(summary: object, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="summary|schema version"):
        resolve_summary_profile(summary, tmp_path)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("Equivalent_Representation_Count", True, "non-negative integer"),
        ("Equivalent_Representation_Count", "01", "non-negative integer"),
        ("Equivalent_Representation_Count", -1, "non-negative integer"),
        ("Molecular_Identity", 1, "must be a string"),
        ("Molecular_Identity_Status", "legacy-selected-among-multiple", "multiple hypotheses"),
    ],
)
def test_schema_three_rejects_additional_identity_type_and_count_contradictions(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    row = copy.deepcopy(GOOD_IDENTITY_ROW)
    row[field] = value

    with pytest.raises(ValueError, match=message):
        resolve_summary_profile(_schema_three_summary(profile, rows=[row]), tmp_path)


def test_schema_three_accepts_tsv_string_counts_and_consistent_unresolved_row(tmp_path: Path) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    resolved = copy.deepcopy(GOOD_IDENTITY_ROW)
    resolved["Equivalent_Representation_Count"] = "1"
    resolved["Identity_Hypothesis_Count"] = "1"
    unresolved = {
        **GOOD_IDENTITY_ROW,
        "Molecular_Identity": "",
        "Molecular_Identity_Status": "unresolved",
        "Equivalent_Representation_Count": "0",
        "Identity_Hypothesis_Count": "1",
    }

    assert resolve_summary_profile(_schema_three_summary(profile, rows=[resolved, unresolved]), tmp_path).sha256 == (
        profile.digest
    )


def test_schema_three_rejects_identity_fields_on_a_negative_row(tmp_path: Path) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    negative = {
        "Motif": "None",
        "POS": "None",
        "Confidence": "Negative",
        "Molecular_Identity": "",
        "Molecular_Identity_Status": "unresolved",
        "Equivalent_Representation_Count": 0,
        "Identity_Hypothesis_Count": 0,
    }

    with pytest.raises(ValueError, match="negative caller rows"):
        resolve_summary_profile(_schema_three_summary(profile, rows=[negative]), tmp_path)


@pytest.mark.parametrize(
    ("steps", "message"),
    [
        ("not-an-array", "steps must be an array"),
        ([1], "step must be an object"),
        ([{"step": STEP_KESTREL, "parsed_result": {"data": {}}}], "step data must be an array"),
        ([{"step": STEP_ADVNTR, "parsed_result": {"data": [1]}}], "result row must be an object"),
    ],
)
def test_schema_three_rejects_malformed_caller_step_containers(
    tmp_path: Path,
    steps: object,
    message: str,
) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    summary["steps"] = steps
    if isinstance(steps, list) and any(isinstance(step, dict) and step.get("step") == STEP_ADVNTR for step in steps):
        _record_advntr_evidence_digest(summary, profile)

    with pytest.raises(ValueError, match=message):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize("parsed_result", [None, [], "malformed"])
def test_schema_three_rejects_nonmapping_caller_results(tmp_path: Path, parsed_result: object) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    summary["steps"] = [{"step": STEP_KESTREL, "parsed_result": parsed_result}]

    with pytest.raises(ValueError, match="parsed_result must be an object"):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize("disposition", [None, "excluded-artifact"])
def test_schema_three_positive_advntr_row_requires_governed_evidence_disposition(
    tmp_path: Path, disposition: str | None
) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    row = {
        "VID": "25561",
        **GOOD_IDENTITY_ROW,
    }
    if disposition is not None:
        row["Evidence_Disposition"] = disposition
    summary = _schema_three_summary(profile)
    _record_advntr_evidence_digest(summary, profile)
    summary["steps"] = [{"step": STEP_ADVNTR, "parsed_result": {"data": [row]}}]

    with pytest.raises(ValueError, match="Evidence_Disposition"):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize(
    ("state", "disposition"),
    [
        ("I10_2_A_LEN1", "admissible"),
        ("NOT_GOVERNED", "identity-insufficient"),
    ],
)
def test_schema_three_binds_advntr_disposition_to_the_governed_variant_state(
    tmp_path: Path,
    state: str,
    disposition: str,
) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    _record_advntr_evidence_digest(summary, profile)
    summary["steps"] = [
        {
            "step": STEP_ADVNTR,
            "parsed_result": {
                "data": [
                    {
                        "VID": "25561",
                        "Variant": state,
                        "Evidence_Disposition": disposition,
                        **GOOD_IDENTITY_ROW,
                    }
                ]
            },
        }
    ]

    with pytest.raises(ValueError, match="governed Variant state"):
        resolve_summary_profile(summary, tmp_path)


def test_schema_three_advntr_step_requires_a_non_null_evidence_digest(tmp_path: Path) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    summary["steps"] = [{"step": STEP_ADVNTR, "parsed_result": {"data": []}}]

    with pytest.raises(ValueError, match="non-null evidence digest"):
        resolve_summary_profile(summary, tmp_path)


def test_schema_three_binds_recorded_evidence_digest_to_the_profile_component(tmp_path: Path) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    summary["advntr_evidence_digest"] = "0" * 64

    with pytest.raises(ValueError, match="profile.*adVNTR evidence digest|evidence digest.*profile"):
        resolve_summary_profile(summary, tmp_path)


def test_snapshot_verification_rejects_legacy_provenance_and_metadata_mismatch(tmp_path: Path) -> None:
    profile, snapshot = _snapshot_packaged(tmp_path)
    legacy = DecisionProfileProvenance(None, None, None, None, None, None)
    with pytest.raises(ValueError, match="recorded schema-3 provenance"):
        verify_profile_snapshot(legacy, snapshot)

    summary = _schema_three_summary(profile)
    summary["decision_profile_id"] = "different-recorded-id"
    with pytest.raises(ValueError, match="metadata differs"):
        resolve_summary_profile(summary, tmp_path)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda row: row.pop("Molecular_Identity_Status"), "identity fields"),
        (lambda row: row.__setitem__("Identity_Hypothesis_Count", 2), "unique"),
        (lambda row: row.__setitem__("Equivalent_Representation_Count", 0), "resolved"),
        (lambda row: row.__setitem__("Molecular_Identity_Status", "unsupported"), "status"),
        (
            lambda row: row.__setitem__("Molecular_Identity", "MUC1-X-60-coding-v2|60|59|-|C"),
            "Unsupported molecular identity reference",
        ),
    ],
)
def test_schema_three_rejects_invalid_positive_identity_rows(tmp_path: Path, mutate, message: str) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    row = copy.deepcopy(GOOD_IDENTITY_ROW)
    mutate(row)

    with pytest.raises(ValueError, match=message):
        resolve_summary_profile(_schema_three_summary(profile, rows=[row]), tmp_path)


def test_schema_three_accepts_unchanged_negative_caller_rows(tmp_path: Path) -> None:
    profile, _snapshot = _snapshot_packaged(tmp_path)
    summary = _schema_three_summary(profile)
    _record_advntr_evidence_digest(summary, profile)
    summary["steps"] = [
        {
            "step": STEP_KESTREL,
            "parsed_result": {"data": [{"Motif": "None", "POS": "None", "Confidence": "Negative"}]},
        },
        {"step": STEP_ADVNTR, "parsed_result": {"data": [{"VID": "Negative"}]}},
    ]

    assert resolve_summary_profile(summary, tmp_path).sha256 == profile.digest


def test_explicit_profile_summary_records_source_not_operator_path(tmp_path: Path) -> None:
    packaged = load_packaged_decision_profile()
    document = json.loads(packaged.canonical_bytes)
    document["profile_id"] = "unit-test-explicit"
    document["profile_revision"] = "test-1"
    document["profile_kind"] = "explicit-custom"
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    inventory["/components/kestrel/duplicate_flagging/flag_name"]["value"] = "Unit_Test_Duplicate"
    explicit = parse_decision_profile(json.dumps(document), packaged_document=packaged.document)
    snapshot_decision_profile(explicit, tmp_path / "provenance" / "decision_profile.json")

    summary = _schema_three_summary(explicit)
    recorded = resolve_summary_profile(summary, tmp_path)

    assert recorded.source == "explicit-cli"
    assert recorded.profile_kind == "explicit-custom"
    assert recorded.snapshot_path == "provenance/decision_profile.json"
    assert str(tmp_path) not in json.dumps(summary)


def test_failed_atomic_snapshot_leaves_no_destination_or_temporary_file(tmp_path: Path, monkeypatch) -> None:
    from vntyper.scripts import profile_provenance

    profile = load_packaged_decision_profile()
    destination = tmp_path / "provenance" / "decision_profile.json"

    def fail_replace(_source: object, _destination: object) -> None:
        raise OSError("injected install failure")

    monkeypatch.setattr(profile_provenance.os, "replace", fail_replace)

    with pytest.raises(RuntimeError, match="Failed to snapshot decision profile"):
        snapshot_decision_profile(profile, destination)

    assert not destination.exists()
    assert list(destination.parent.iterdir()) == []


def test_verify_profile_snapshot_accepts_historical_revision_1_packaged_profile(tmp_path: Path) -> None:
    from vntyper.scripts.canonical_json import load_strict_json_object
    from vntyper.scripts.profile_provenance import HISTORICAL_REVISION_1_PACKAGED_SHA256

    packaged = load_packaged_decision_profile()
    rev1 = load_strict_json_object(packaged.canonical_bytes)
    rev1["profile_revision"] = "1"
    del rev1["inventory"]["/components/kestrel/confidence_assignment/reporting_floor"]

    snapshot = tmp_path / "provenance" / "decision_profile.json"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    raw = canonical_json_bytes(rev1)
    snapshot.write_bytes(raw)

    provenance = DecisionProfileProvenance(
        profile_id="vntyper-packaged-default",
        profile_revision="1",
        profile_kind="packaged",
        source="package",
        sha256=HISTORICAL_REVISION_1_PACKAGED_SHA256,
        snapshot_path=str(snapshot),
    )

    verified = verify_profile_snapshot(provenance, snapshot)
    assert verified == provenance

    # Tampered byte verification against HISTORICAL_REVISION_1_PACKAGED_SHA256
    tampered = dict(rev1)
    tampered["inventory"] = dict(rev1["inventory"])
    tampered["inventory"]["/components/kestrel/duplicate_flagging/flag_name"] = {
        "class": "explicit-custom",
        "value": "Changed_Name",
    }
    tampered_bytes = canonical_json_bytes(tampered)
    tampered_digest = hashlib.sha256(tampered_bytes).hexdigest()
    snapshot.write_bytes(tampered_bytes)
    tampered_prov = DecisionProfileProvenance(
        profile_id="vntyper-packaged-default",
        profile_revision="1",
        profile_kind="packaged",
        source="package",
        sha256=tampered_digest,
        snapshot_path=str(snapshot),
    )
    with pytest.raises(
        ValueError, match="recorded package-source snapshot does not match the checked-in packaged profile"
    ):
        verify_profile_snapshot(tampered_prov, snapshot)
