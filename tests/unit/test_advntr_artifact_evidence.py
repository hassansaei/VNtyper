"""Canonical governed evidence for recurrent adVNTR State strings (#267)."""

from __future__ import annotations

import json
from collections import Counter
from copy import deepcopy
from dataclasses import FrozenInstanceError
from pathlib import Path

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping
from vntyper.modules.advntr import artifact_evidence as artifact_evidence_module
from vntyper.modules.advntr.artifact_evidence import (
    ASSERTION,
    ArtifactEvidence,
    evidence_disposition_for_state,
    load_artifact_evidence,
    load_packaged_artifact_evidence,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.flagging import add_flags
from vntyper.scripts.molecular_identity import EvidenceDisposition

pytestmark = pytest.mark.unit

EXPECTED_STATES = (
    "I10_2_A_LEN1",
    "D8_2&D9_2&I9_2_A_LEN9",
    "D2_2&I2_2_C_LEN5",
    "I39_2_A_LEN4",
    "I52_2_A_LEN7",
    "I45_2_A_LEN4",
    "D45_2&I45_2_A_LEN2",
    "D14_2&I14_2_G_LEN14",
    "D58_2&D59_2",
    "I60_2_A_LEN10",
    "I14_2_G_LEN16",
    "I18_2_T_LEN1",
    "I21_2_G_LEN4",
    "D29_2&I29_2_A_LEN2",
    "D8_2&I8_2_A_LEN20",
    "D20_2&D21_2",
    "D21_2&D22_2",
    "I14_2_A_LEN1",
    "I11_2_G_LEN1",
    "I26_7_A_LEN25",
    "D17_2&D18_2&D19_2&D20_2&D21_2",
    "I14_2_C_LEN4",
    "I23_6_G_LEN1",
    "I21_2_T_LEN1",
)
CONFIRMED_STATE = "D58_2&D59_2"
EXPECTED_ARTIFACT_DIGEST = "8bb68bd5fba539feee6feb240f113aaa24fc65b5a1e55776c58cea83db5654b0"
MODULE_DIR = Path(advntr_genotyping.__file__).parent
ARTIFACT_PATH = MODULE_DIR / "advntr_artifact_evidence.json"
DIGEST_PATH = MODULE_DIR / "advntr_artifact_evidence.sha256"


def _flagging_config() -> dict:
    return json.loads((MODULE_DIR / "advntr_config.json").read_text(encoding="utf-8"))["flagging_rules"]


def _one_edit_near_match(state: str) -> str:
    replacement = "2" if state[-1] != "2" else "3"
    return state[:-1] + replacement


@pytest.fixture(scope="module")
def evidence() -> ArtifactEvidence:
    """Load the packaged governed evidence once for immutable-value checks."""
    return load_packaged_artifact_evidence()


class TestStrictCanonicalJson:
    def test_duplicate_object_keys_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="duplicate JSON object key: state"):
            load_strict_json_object(b'{"state":"a","state":"b"}')

    @pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity"])
    def test_nonfinite_json_constants_are_rejected(self, constant: str) -> None:
        with pytest.raises(ValueError, match="non-finite JSON constant"):
            load_strict_json_object(f'{{"value":{constant}}}'.encode())

    def test_top_level_json_must_be_an_object(self) -> None:
        with pytest.raises(ValueError, match="top-level JSON value must be an object"):
            load_strict_json_object(b"[]")

    def test_rfc8785_known_answer_bytes_and_digest_are_literal(self) -> None:
        value = {"b": 1, "a": "é"}

        assert canonical_json_bytes(value) == b'{"a":"\xc3\xa9","b":1}\n'
        assert canonical_sha256(value) == "2f8b31e48ae7ba3230b3c8499e02f919d89bdfebec9f9bf3872e4dee8458bf14"


class TestPackagedArtifact:
    def test_packaged_bytes_are_canonical_with_one_terminal_newline(self) -> None:
        raw = ARTIFACT_PATH.read_bytes()
        decoded = load_strict_json_object(raw)

        assert raw == canonical_json_bytes(decoded)
        assert raw.endswith(b"\n") and not raw.endswith(b"\n\n")

    def test_sidecar_matches_the_literal_canonical_digest(self) -> None:
        expected = DIGEST_PATH.read_text(encoding="ascii").strip()

        assert expected == EXPECTED_ARTIFACT_DIGEST
        assert expected == canonical_sha256(load_strict_json_object(ARTIFACT_PATH.read_bytes()))

    def test_distribution_manifest_names_both_governed_files(self) -> None:
        pyproject = Path("pyproject.toml").read_text(encoding="utf-8")

        assert '"modules/advntr/advntr_artifact_evidence.json"' in pyproject
        assert '"modules/advntr/advntr_artifact_evidence.sha256"' in pyproject


class TestTypedEvidence:
    def test_exact_governed_scope_and_limited_assertion(self, evidence: ArtifactEvidence) -> None:
        assert evidence.schema_version == 1
        assert evidence.cohort_name == "renome"
        assert evidence.advntr_version_upper_bound_exclusive == "2.0.4"
        assert evidence.exact_advntr_version is None
        assert (
            evidence.assertion
            == ASSERTION
            == ("A carried-forward recurrent adVNTR State is insufficient for molecular identity.")
        )

    def test_exactly_24_distinct_active_states_are_loaded(self, evidence: ArtifactEvidence) -> None:
        assert evidence.active_states == frozenset(EXPECTED_STATES)
        assert len(evidence.entries) == len(evidence.active_states) == 24
        assert all(entry.active for entry in evidence.entries)

    def test_one_entry_is_confirmed_and_23_are_pending(self, evidence: ArtifactEvidence) -> None:
        statuses = Counter(entry.status for entry in evidence.entries)

        assert statuses == {"confirmed_artifact": 1, "pending_renome_revalidation": 23}
        assert {entry.state for entry in evidence.entries if entry.status == "confirmed_artifact"} == {CONFIRMED_STATE}

    def test_unknown_cohort_measurements_remain_json_null(self, evidence: ArtifactEvidence) -> None:
        assert all(entry.count is None for entry in evidence.entries)
        assert all(entry.denominator is None for entry in evidence.entries)
        assert all(entry.frequency is None for entry in evidence.entries)

    def test_loaded_types_are_frozen(self, evidence: ArtifactEvidence) -> None:
        with pytest.raises(FrozenInstanceError):
            evidence.cohort_name = "invented"  # type: ignore[misc]
        with pytest.raises(FrozenInstanceError):
            evidence.entries[0].state = "invented"  # type: ignore[misc]

    def test_all_active_states_are_identity_insufficient_regardless_of_status(self, evidence: ArtifactEvidence) -> None:
        assert {entry.status for entry in evidence.entries} == {
            "confirmed_artifact",
            "pending_renome_revalidation",
        }
        assert {evidence_disposition_for_state(entry.state, evidence) for entry in evidence.entries} == {
            EvidenceDisposition("identity-insufficient")
        }

    def test_unlisted_state_is_admissible(self, evidence: ArtifactEvidence) -> None:
        assert evidence_disposition_for_state("I22_2_G_LEN1", evidence) == EvidenceDisposition("admissible")

    @pytest.mark.parametrize("state", ["", None, 1])
    def test_disposition_rejects_invalid_state(self, evidence: ArtifactEvidence, state: object) -> None:
        with pytest.raises(ValueError, match="requires a non-empty State string"):
            evidence_disposition_for_state(state, evidence)  # type: ignore[arg-type]

    def test_disposition_rejects_unverified_evidence(self) -> None:
        with pytest.raises(ValueError, match="requires verified ArtifactEvidence"):
            evidence_disposition_for_state("I22_2_G_LEN1", object())  # type: ignore[arg-type]

    def test_file_loader_rejects_a_wrong_expected_digest(self, tmp_path: Path) -> None:
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(ARTIFACT_PATH.read_bytes())

        with pytest.raises(ValueError, match="artifact evidence digest mismatch"):
            load_artifact_evidence(artifact, expected_digest="0" * 64)

    def test_file_loader_accepts_the_literal_expected_digest(self, tmp_path: Path) -> None:
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(ARTIFACT_PATH.read_bytes())

        loaded = load_artifact_evidence(artifact, expected_digest=EXPECTED_ARTIFACT_DIGEST)

        assert loaded.digest == EXPECTED_ARTIFACT_DIGEST

    def test_file_loader_rejects_noncanonical_bytes(self, tmp_path: Path) -> None:
        artifact = tmp_path / "evidence.json"
        artifact.write_text('{"schema_version": 1}\n', encoding="utf-8")

        with pytest.raises(ValueError, match="must use canonical RFC 8785 bytes"):
            load_artifact_evidence(artifact)

    @pytest.mark.parametrize(
        ("field", "invalid", "message"),
        [
            ("schema_version", True, "schema_version must be 1"),
            ("schema_version", 2, "schema_version must be 1"),
            ("assertion", "broader claim", "assertion differs"),
            ("cohort_name", "other", "cohort_name must be renome"),
            ("advntr_version_upper_bound_exclusive", "2.0.5", "version upper bound must be 2.0.4"),
            ("exact_advntr_version", "", "exact_advntr_version must be a non-empty string or null"),
            ("exact_advntr_version", 204, "exact_advntr_version must be a non-empty string or null"),
            ("entries", [], "entries must be a non-empty list"),
            ("entries", {}, "entries must be a non-empty list"),
        ],
    )
    def test_document_contract_rejects_invalid_root_values(
        self, tmp_path: Path, field: str, invalid: object, message: str
    ) -> None:
        document = load_strict_json_object(ARTIFACT_PATH.read_bytes())
        document[field] = invalid
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(canonical_json_bytes(document))

        with pytest.raises(ValueError, match=message):
            load_artifact_evidence(artifact)

    @pytest.mark.parametrize("edit", ["missing", "extra"])
    def test_document_contract_requires_exact_root_fields(self, tmp_path: Path, edit: str) -> None:
        document = load_strict_json_object(ARTIFACT_PATH.read_bytes())
        if edit == "missing":
            del document["assertion"]
        else:
            document["invented"] = None
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(canonical_json_bytes(document))

        with pytest.raises(ValueError, match="artifact evidence root fields differ"):
            load_artifact_evidence(artifact)

    @pytest.mark.parametrize(
        ("field", "invalid", "message"),
        [
            ("state", "", "state must be a non-empty string"),
            ("state", 1, "state must be a non-empty string"),
            ("status", "invented", "status is unsupported"),
            ("status", [], "status is unsupported"),
            ("active", 1, "active must be a Boolean"),
            ("count", -1, "count must be a non-negative integer or null"),
            ("count", True, "count must be a non-negative integer or null"),
            ("denominator", 1.5, "denominator must be a non-negative integer or null"),
            ("frequency", -0.1, "frequency must be between 0 and 1 or null"),
            ("frequency", 1.1, "frequency must be between 0 and 1 or null"),
            ("frequency", True, "frequency must be between 0 and 1 or null"),
        ],
    )
    def test_document_contract_rejects_invalid_entry_values(
        self, tmp_path: Path, field: str, invalid: object, message: str
    ) -> None:
        document = load_strict_json_object(ARTIFACT_PATH.read_bytes())
        document["entries"][0][field] = invalid
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(canonical_json_bytes(document))

        with pytest.raises(ValueError, match=message):
            load_artifact_evidence(artifact)

    def test_document_contract_rejects_nonobject_entry(self, tmp_path: Path) -> None:
        document = load_strict_json_object(ARTIFACT_PATH.read_bytes())
        document["entries"][0] = "not an object"
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(canonical_json_bytes(document))

        with pytest.raises(ValueError, match=r"entries\[0\] must be an object"):
            load_artifact_evidence(artifact)

    @pytest.mark.parametrize("edit", ["missing", "extra"])
    def test_document_contract_requires_exact_entry_fields(self, tmp_path: Path, edit: str) -> None:
        document = load_strict_json_object(ARTIFACT_PATH.read_bytes())
        if edit == "missing":
            del document["entries"][0]["status"]
        else:
            document["entries"][0]["invented"] = None
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(canonical_json_bytes(document))

        with pytest.raises(ValueError, match=r"entries\[0\] fields differ"):
            load_artifact_evidence(artifact)

    def test_document_contract_rejects_duplicate_states(self, tmp_path: Path) -> None:
        document = load_strict_json_object(ARTIFACT_PATH.read_bytes())
        duplicate = deepcopy(document["entries"][0])
        document["entries"].append(duplicate)
        artifact = tmp_path / "evidence.json"
        artifact.write_bytes(canonical_json_bytes(document))

        with pytest.raises(ValueError, match="State strings must be unique"):
            load_artifact_evidence(artifact)

    def test_packaged_loader_rejects_malformed_digest_sidecar(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        digest = tmp_path / "evidence.sha256"
        digest.write_text("NOT-A-DIGEST\n", encoding="ascii")
        monkeypatch.setattr(artifact_evidence_module, "_PACKAGED_DIGEST", digest)

        with pytest.raises(ValueError, match="64 lowercase hexadecimal characters"):
            load_packaged_artifact_evidence()


class TestIndependentReachabilityOracle:
    def test_literal_oracle_equals_the_production_polymorphic_rule_inventory(self) -> None:
        configured = _flagging_config()["Polymorphic_Call"]
        production_states = configured["all"][0]["right"]["literal"]

        assert set(production_states) == set(EXPECTED_STATES)
        assert len(production_states) == len(EXPECTED_STATES)

    def test_every_literal_state_fires_exactly_the_polymorphic_rule(self) -> None:
        configured = _flagging_config()["Polymorphic_Call"]

        for state in EXPECTED_STATES:
            row = pd.DataFrame([{"Variant": state, "RU": "2", "NumberOfSupportingReads": 50}])
            flagged = add_flags(row, {"Polymorphic_Call": configured})
            assert flagged.iloc[0]["Flag"] == "Polymorphic_Call", state

    def test_one_edit_near_matches_do_not_fire(self) -> None:
        configured = _flagging_config()["Polymorphic_Call"]
        near_matches = tuple(_one_edit_near_match(state) for state in EXPECTED_STATES)

        assert len(set(near_matches)) == len(EXPECTED_STATES)
        assert not set(near_matches).intersection(EXPECTED_STATES)
        for state in near_matches:
            row = pd.DataFrame([{"Variant": state, "RU": "2", "NumberOfSupportingReads": 50}])
            flagged = add_flags(row, {"Polymorphic_Call": configured})
            assert flagged.iloc[0]["Flag"] == "Not flagged", state
