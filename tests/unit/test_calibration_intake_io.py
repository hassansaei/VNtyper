"""Atomic local intake bundle production from strict synthetic declarations."""

import hashlib
import json
import stat
from copy import deepcopy
from dataclasses import replace
from importlib import import_module
from pathlib import Path
from unittest.mock import call, patch

import pytest

from tests.unit.test_calibration_intake_contract import synthetic_intake
from vntyper.scripts.calibration_identity import ArtifactFingerprint
from vntyper.scripts.calibration_read_fingerprints import LogicalReadFingerprint
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

pytestmark = pytest.mark.unit


def _digest(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _fingerprint(key: str) -> ArtifactFingerprint:
    return ArtifactFingerprint(
        key,
        _digest(key + ":bytes"),
        None,
        LogicalReadFingerprint(
            _digest(key + ":alignment"),
            _digest(key + ":named"),
            _digest(key + ":unnamed"),
            2,
            True,
            (),
        ),
    )


def _write_declaration(tmp_path: Path, raw: dict | None = None) -> tuple[Path, dict]:
    value = deepcopy(synthetic_intake() if raw is None else raw)
    for artifact in value["artifacts"]:
        artifact["path"] = str(tmp_path / f"{artifact['key']}.bam")
        artifact["expected_sha256"] = None
        if artifact["mate_path"] is not None:
            artifact["mate_path"] = str(tmp_path / f"{artifact['key']}-mate.fastq")
    path = tmp_path / "intake.json"
    path.write_bytes(canonical_json_bytes(value))
    return path, value


def test_prepare_intake_bundle_writes_exact_private_canonical_artifacts(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path, _ = _write_declaration(tmp_path)
    output = tmp_path / "bundle"
    fingerprint = _fingerprint("artifact-001")

    with patch.object(module, "fingerprint_input_artifact", return_value=fingerprint) as reader:
        result = module.prepare_intake_bundle(
            declaration_path,
            output,
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
            temporary_parent=tmp_path,
        )

    assert reader.call_args_list == [
        call(result.declaration.artifacts[0], reference_path=None, reference_sha256=None, temporary_parent=tmp_path)
    ]
    assert result.output == output
    assert result.audit.fingerprints == {"artifact-001": fingerprint}
    assert tuple(member.key for member in result.partition_members) == ("artifact-001",)
    assert {path.name for path in output.iterdir()} == {
        "normalized.json",
        "dedup_audit.json",
        "partitions.json",
        "provenance.json",
    }
    assert all(stat.S_IMODE(path.stat().st_mode) == 0o600 for path in output.iterdir())
    assert (
        canonical_json_bytes(load_strict_json_object((output / "normalized.json").read_bytes()))
        == (output / "normalized.json").read_bytes()
    )
    assert canonical_sha256(load_strict_json_object((output / "normalized.json").read_bytes())) == (
        result.declaration.sha256
    )
    provenance = load_strict_json_object((output / "provenance.json").read_bytes())
    assert provenance["schema_version"] == "calibration-intake-provenance-v2"
    assert provenance["cram_references"] == []
    assert provenance["intake_sha256"] == result.declaration.sha256
    assert provenance["identity_audit_sha256"] == result.audit.sha256
    assert provenance["partition_sha256"] == canonical_sha256(
        load_strict_json_object((output / "partitions.json").read_bytes())
    )
    assert result.provenance_sha256 == canonical_sha256(provenance)
    assert [item["path"] for item in provenance["files"]] == [
        "dedup_audit.json",
        "normalized.json",
        "partitions.json",
    ]
    for item in provenance["files"]:
        payload = (output / item["path"]).read_bytes()
        assert item == {"path": item["path"], "size_bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


def test_cram_references_are_exactly_assembly_bound_and_artifacts_are_read_once_in_key_order(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    raw = synthetic_intake()
    second = deepcopy(raw["artifacts"][0])
    second.update(key="artifact-000", path="/replaced", format="CRAM", assembly="assembly-b")
    raw["artifacts"].append(second)
    declaration_path, _ = _write_declaration(tmp_path, raw)
    reference = tmp_path / "reference.fa"
    reference.write_text(">synthetic\nACGT\n", encoding="utf-8")
    reference_sha256 = hashlib.sha256(reference.read_bytes()).hexdigest()
    fingerprints = {key: _fingerprint(key) for key in ("artifact-000", "artifact-001")}
    observed_keys: list[str] = []

    def read_cram(artifact, **_kwargs):
        observed_keys.append(artifact.key)
        return read_module.ArtifactReadEvidence(fingerprints[artifact.key], reference_sha256)

    def read_other(artifact, **_kwargs):
        observed_keys.append(artifact.key)
        return fingerprints[artifact.key]

    read_module = import_module("vntyper.scripts.calibration_read_io")
    with (
        patch.object(module, "fingerprint_input_artifact_evidence", side_effect=read_cram) as cram_reader,
        patch.object(module, "fingerprint_input_artifact", side_effect=read_other) as reader,
    ):
        result = module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={"assembly-b": module.PinnedCramReference(reference.resolve(), reference_sha256)},
        )

    assert cram_reader.call_args_list[0].args[0].key == "artifact-000"
    assert cram_reader.call_args_list[0].kwargs == {
        "reference_path": reference.resolve(),
        "reference_sha256": reference_sha256,
        "temporary_parent": None,
    }
    assert reader.call_args_list[0].args[0].key == "artifact-001"
    assert reader.call_args_list[0].kwargs == {
        "reference_path": None,
        "reference_sha256": None,
        "temporary_parent": None,
    }
    assert observed_keys == ["artifact-000", "artifact-001"]
    assert set(result.audit.fingerprints) == {"artifact-000", "artifact-001"}


def test_verified_cram_reference_path_and_digest_are_bound_into_provenance(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    read_module = import_module("vntyper.scripts.calibration_read_io")
    raw = synthetic_intake()
    raw["artifacts"][0].update(format="CRAM", assembly="assembly-a")
    declaration_path, _ = _write_declaration(tmp_path, raw)
    reference_a = module.PinnedCramReference((tmp_path / "reference-a.fa").resolve(), "a" * 64)
    reference_b = module.PinnedCramReference((tmp_path / "reference-b.fa").resolve(), "b" * 64)

    def build(output: Path, reference):
        evidence = read_module.ArtifactReadEvidence(_fingerprint("artifact-001"), reference.sha256)
        with patch.object(module, "fingerprint_input_artifact_evidence", return_value=evidence):
            result = module.prepare_intake_bundle(
                declaration_path,
                output,
                preprocessing_priority=("synthetic-preprocessing-v1",),
                cram_references={"assembly-a": reference},
            )
        return result, load_strict_json_object((output / "provenance.json").read_bytes())

    first, first_document = build(tmp_path / "first", reference_a)
    second, second_document = build(tmp_path / "second", reference_b)

    assert first_document["schema_version"] == "calibration-intake-provenance-v2"
    assert first_document["cram_references"] == [
        {"assembly": "assembly-a", "path": str(reference_a.path), "sha256": reference_a.sha256}
    ]
    assert second_document["cram_references"] == [
        {"assembly": "assembly-a", "path": str(reference_b.path), "sha256": reference_b.sha256}
    ]
    assert first.provenance_sha256 != second.provenance_sha256


def test_cram_reference_observed_digest_must_match_the_declared_binding(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    read_module = import_module("vntyper.scripts.calibration_read_io")
    raw = synthetic_intake()
    raw["artifacts"][0].update(format="CRAM", assembly="assembly-a")
    declaration_path, _ = _write_declaration(tmp_path, raw)
    reference = module.PinnedCramReference((tmp_path / "reference.fa").resolve(), "a" * 64)
    evidence = read_module.ArtifactReadEvidence(_fingerprint("artifact-001"), "b" * 64)

    with (
        patch.object(module, "fingerprint_input_artifact_evidence", return_value=evidence),
        pytest.raises(ValueError, match="observed CRAM reference"),
    ):
        module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={"assembly-a": reference},
        )
    assert not (tmp_path / "bundle").exists()


def test_lazy_cram_reference_loader_runs_only_after_atomic_output_preflight(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    output = tmp_path / "bundle"
    output.mkdir()
    called = False

    def load_references():
        nonlocal called
        called = True
        return {}

    with pytest.raises(ValueError, match="already exists"):
        module.prepare_intake_bundle(
            tmp_path / "missing.json",
            output,
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_reference_loader=load_references,
        )
    assert called is False


@pytest.mark.parametrize("sources", ["neither", "both", "invalid-loader", "invalid-result"])
def test_cram_reference_source_boundary_is_closed(sources: str, tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path, _ = _write_declaration(tmp_path)
    kwargs: dict[str, object] = {}
    message = "exactly one"
    if sources == "both":
        kwargs = {"cram_references": {}, "cram_reference_loader": dict}
    elif sources == "invalid-loader":
        kwargs = {"cram_reference_loader": []}
        message = "callable"
    elif sources == "invalid-result":
        kwargs = {"cram_reference_loader": list}
        message = "return a mapping"

    with pytest.raises(ValueError, match=message):
        module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            **kwargs,  # type: ignore[arg-type]
        )
    assert not (tmp_path / "bundle").exists()


def test_locked_membership_only_is_retained_without_an_artifact_read(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    raw = synthetic_intake()
    locked_specimen = deepcopy(raw["specimens"][0])
    locked_specimen.update(key="locked-member", individual_key="locked-individual")
    locked_assignment = deepcopy(raw["assignments"][0])
    locked_assignment.update(
        specimen_key="locked-member",
        role="locked-heldout",
        provenance="external-custodian",
        groups={name: [f"{name}:locked"] for name in locked_assignment["groups"]},
    )
    raw["specimens"].append(locked_specimen)
    raw["assignments"].append(locked_assignment)
    declaration_path, _ = _write_declaration(tmp_path, raw)

    with patch.object(module, "fingerprint_input_artifact", return_value=_fingerprint("artifact-001")) as reader:
        result = module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )

    reader.assert_called_once()
    assert {row.key for row in result.declaration.specimens} == {"locked-member", "sample-001"}
    assert tuple(member.key for member in result.partition_members) == ("artifact-001",)


def test_all_locked_declaration_has_no_local_artifacts_and_reads_nothing(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    raw = synthetic_intake()
    raw["artifacts"] = []
    raw["aliases"] = []
    raw["truth"] = []
    raw["assignments"][0].update(role="locked-heldout", provenance="external-custodian")
    declaration_path, _ = _write_declaration(tmp_path, raw)

    with (
        patch.object(module, "fingerprint_input_artifact") as reader,
        pytest.raises(ValueError, match="no_local_artifacts"),
    ):
        module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    reader.assert_not_called()
    assert not (tmp_path / "bundle").exists()


@pytest.mark.parametrize("invalid", ["duplicate", "unknown", "relative"])
def test_declaration_is_strictly_decoded_and_paths_validated_before_artifact_reads(
    tmp_path: Path, invalid: str
) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path, raw = _write_declaration(tmp_path)
    if invalid == "duplicate":
        declaration_path.write_text('{"schema_version":"calibration-intake-v1","schema_version":"x"}', encoding="utf-8")
    elif invalid == "unknown":
        raw["extra"] = True
        declaration_path.write_text(json.dumps(raw), encoding="utf-8")
    else:
        raw["artifacts"][0]["path"] = "relative.bam"
        declaration_path.write_bytes(canonical_json_bytes(raw))

    with patch.object(module, "fingerprint_input_artifact") as reader, pytest.raises(ValueError):
        module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    reader.assert_not_called()


@pytest.mark.parametrize("mode", ["missing", "extra", "invalid_pin"])
def test_cram_reference_registry_is_closed_and_validated_before_artifact_reads(tmp_path: Path, mode: str) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    raw = synthetic_intake()
    raw["artifacts"][0].update(format="CRAM", assembly="assembly-a")
    declaration_path, _ = _write_declaration(tmp_path, raw)
    reference = tmp_path / "reference.fa"
    reference.write_text(">synthetic\nAC\n", encoding="utf-8")
    if mode == "missing":
        references = {}
    elif mode == "extra":
        references = {
            "assembly-a": module.PinnedCramReference(reference.resolve(), _digest("reference")),
            "unused": module.PinnedCramReference(reference.resolve(), _digest("unused")),
        }
    else:
        references = {
            "assembly-a": replace(module.PinnedCramReference(reference.resolve(), _digest("reference")), sha256=True)
        }

    with patch.object(module, "fingerprint_input_artifact") as reader, pytest.raises(ValueError, match="reference"):
        module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references=references,
        )
    reader.assert_not_called()


def test_source_change_or_fingerprint_failure_leaves_no_partial_bundle(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path, _ = _write_declaration(tmp_path)
    output = tmp_path / "bundle"

    def change_source(_artifact, **_kwargs):
        declaration_path.write_bytes(b"{}\n")
        return _fingerprint("artifact-001")

    with (
        patch.object(module, "fingerprint_input_artifact", side_effect=change_source),
        pytest.raises(RuntimeError, match="changed"),
    ):
        module.prepare_intake_bundle(
            declaration_path,
            output,
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    assert not output.exists()
    assert not tuple(tmp_path.glob(".bundle.*"))

    declaration_path, _ = _write_declaration(tmp_path)
    with (
        patch.object(module, "fingerprint_input_artifact", side_effect=RuntimeError("synthetic read failure")),
        pytest.raises(RuntimeError, match="synthetic read failure"),
    ):
        module.prepare_intake_bundle(
            declaration_path,
            output,
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    assert not output.exists()
    assert not tuple(tmp_path.glob(".bundle.*"))


def test_existing_output_is_rejected_without_reading_declaration_or_artifacts(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path = tmp_path / "missing.json"
    output = tmp_path / "bundle"
    output.mkdir()

    with (
        patch.object(module, "fingerprint_input_artifact") as reader,
        pytest.raises(ValueError, match="already exists"),
    ):
        module.prepare_intake_bundle(
            declaration_path,
            output,
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    reader.assert_not_called()


def test_declaration_symlink_is_rejected_without_artifact_reads(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path, _ = _write_declaration(tmp_path)
    linked = tmp_path / "linked-intake.json"
    linked.symlink_to(declaration_path)

    with patch.object(module, "fingerprint_input_artifact") as reader, pytest.raises(ValueError, match="symlink"):
        module.prepare_intake_bundle(
            linked,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    reader.assert_not_called()


def test_declaration_open_requests_nonblocking_file_type_admission(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path, _ = _write_declaration(tmp_path)
    descriptor = module.os.open(declaration_path, module.os.O_RDONLY)

    with patch.object(module.os, "open", return_value=descriptor) as opener:
        snapshot = module._open_source(declaration_path)
    module.os.close(snapshot.descriptor)

    assert opener.call_args.args[1] & module.os.O_NONBLOCK


def test_fifo_declaration_without_a_writer_is_rejected_without_hanging_or_artifact_reads(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path = tmp_path / "intake.fifo"
    module.os.mkfifo(declaration_path)

    with patch.object(module, "fingerprint_input_artifact") as reader, pytest.raises(ValueError, match="regular"):
        module.prepare_intake_bundle(
            declaration_path,
            tmp_path / "bundle",
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    reader.assert_not_called()


def test_unexpected_staging_file_fails_exact_manifest_and_is_cleaned(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_intake_io")
    declaration_path, _ = _write_declaration(tmp_path)
    output = tmp_path / "bundle"
    real_write = module._write_private_json

    def contaminate(path: Path, value: object) -> bytes:
        raw = real_write(path, value)
        if path.name == "normalized.json":
            (path.parent / "unexpected.json").write_bytes(b"{}\n")
        return raw

    with (
        patch.object(module, "fingerprint_input_artifact", return_value=_fingerprint("artifact-001")),
        patch.object(module, "_write_private_json", side_effect=contaminate),
        pytest.raises(ValueError, match="inventory"),
    ):
        module.prepare_intake_bundle(
            declaration_path,
            output,
            preprocessing_priority=("synthetic-preprocessing-v1",),
            cram_references={},
        )
    assert not output.exists()
    assert not tuple(tmp_path.glob(".bundle.*"))
