"""Generated read bundles are deterministic, private and atomically published."""

import hashlib
import json
import os
from importlib import import_module
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.unit.test_calibration_sim_protocol import decode, simulation_case, simulation_protocol

pytestmark = pytest.mark.unit


def _module():
    return import_module("calibration_sim.generation")


def test_generate_installs_complete_deterministic_bundle_with_independent_origins(tmp_path: Path):
    protocol = decode(simulation_protocol())
    first = _module().generate_simulation_bundle(protocol, tmp_path / "first")
    second = _module().generate_simulation_bundle(protocol, tmp_path / "second")
    assert first.manifest_sha256 == second.manifest_sha256
    manifest_path = first.output / "manifest.json"
    assert hashlib.sha256(manifest_path.read_bytes()).hexdigest() == first.manifest_sha256
    manifest = json.loads(manifest_path.read_bytes())
    assert manifest["protocol_sha256"] == protocol.sha256
    assert manifest["pair_count"] == 3
    assert manifest["generated_bases"] == 48
    assert manifest["independent_group_count"] == 1
    assert manifest["evidence_status"] == "generated-inputs-only"
    assert manifest["generator"]["python_version"]
    assert set(manifest["generator"]["source_sha256"]) == {
        "haplotypes.py",
        "reads.py",
        "protocol.py",
        "generation.py",
    }
    expected = {
        "protocol.json",
        "case-a/haplotypes.fa",
        "case-a/truth.json",
        "case-a/origins.tsv",
        "case-a/reads_R1.fastq",
        "case-a/reads_R2.fastq",
    }
    assert {row["path"] for row in manifest["files"]} == expected
    assert [row["path"] for row in manifest["files"]] == sorted(expected)
    for entry in manifest["files"]:
        path = first.output / entry["path"]
        assert path.stat().st_mode & 0o777 == 0o600
        data = path.read_bytes()
        assert len(data) == entry["size"]
        assert hashlib.sha256(data).hexdigest() == entry["sha256"]
        assert data == (second.output / entry["path"]).read_bytes()
    assert first.output.stat().st_mode & 0o777 == 0o700
    assert (first.output / "case-a").stat().st_mode & 0o777 == 0o700
    truth = json.loads((first.output / "case-a/truth.json").read_bytes())
    assert truth["allele_repeat_counts"] == [2, 2]
    assert truth["total_repeat_count"] == 4
    assert truth["caller_positive"] is False
    assert truth["truth_variant_ids"] == []
    assert truth["primary"] is True
    reads = (first.output / "case-a/reads_R1.fastq").read_text().splitlines()
    assert len(reads) == 12
    assert reads[0] == "@synthetic-read-000000000000/1"
    assert reads[2] == "+"
    assert len(reads[1]) == len(reads[3]) == 8
    origins = (first.output / "case-a/origins.tsv").read_text().splitlines()
    assert len(origins) == 4
    assert origins[0] == "read\tallele_index\tfragment_start\tfragment_end\tread1_substitutions\tread2_substitutions"


def test_existing_output_is_refused_before_read_generation(tmp_path: Path):
    output = tmp_path / "existing"
    output.mkdir()
    with patch.object(_module(), "generate_read_pairs") as generator:
        with pytest.raises(ValueError, match="already exists"):
            _module().generate_simulation_bundle(decode(simulation_protocol()), output)
        generator.assert_not_called()
    assert list(output.iterdir()) == []


def test_failure_mid_generation_removes_partial_output(tmp_path: Path):
    output = tmp_path / "failed"
    with (
        patch.object(_module(), "generate_read_pairs", side_effect=RuntimeError("synthetic generator failure")),
        pytest.raises(RuntimeError, match="synthetic generator failure"),
    ):
        _module().generate_simulation_bundle(decode(simulation_protocol()), output)
    assert not output.exists()
    assert list(tmp_path.iterdir()) == []


def test_zero_coverage_still_has_empty_read_files_and_explicit_positive_truth(tmp_path: Path):
    case = simulation_case(caller_positive=True, truth_variant_ids=["insertion-1"])
    case["reads"]["pair_count"] = 0
    result = _module().generate_simulation_bundle(decode(simulation_protocol(case)), tmp_path / "zero")
    assert (result.output / "case-a/reads_R1.fastq").read_bytes() == b""
    assert (result.output / "case-a/reads_R2.fastq").read_bytes() == b""
    truth = json.loads((result.output / "case-a/truth.json").read_bytes())
    assert truth["caller_positive"] is True
    assert truth["truth_variant_ids"] == ["insertion-1"]


def test_changed_read_parameters_change_bundle_identity(tmp_path: Path):
    raw = simulation_protocol()
    first = _module().generate_simulation_bundle(decode(raw), tmp_path / "first")
    raw["cases"][0]["reads"]["substitution_rate"] = 1
    second = _module().generate_simulation_bundle(decode(raw), tmp_path / "second")
    assert first.manifest_sha256 != second.manifest_sha256
    assert (first.output / "case-a/reads_R1.fastq").read_bytes() != (
        second.output / "case-a/reads_R1.fastq"
    ).read_bytes()


def test_generated_manifest_rejects_missing_extra_or_corrupt_files(tmp_path: Path):
    result = _module().generate_simulation_bundle(decode(simulation_protocol()), tmp_path / "bundle")
    _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)
    target = result.output / "case-a/reads_R1.fastq"
    original = target.read_bytes()
    target.write_bytes(original + b"changed")
    with pytest.raises(ValueError):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)
    target.write_bytes(original)
    (result.output / "extra").write_bytes(b"extra")
    with pytest.raises(ValueError, match="inventory"):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)
    (result.output / "extra").unlink()
    target.unlink()
    with pytest.raises(ValueError):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)


def test_verification_uses_externally_bound_manifest_digest_and_rejects_symlinks(tmp_path: Path):
    result = _module().generate_simulation_bundle(decode(simulation_protocol()), tmp_path / "bundle")
    with pytest.raises(ValueError, match="manifest"):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256="0" * 64)
    target = result.output / "case-a/truth.json"
    external = tmp_path / "truth.json"
    target.rename(external)
    target.symlink_to(external)
    with pytest.raises(ValueError):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)


def test_extra_empty_directory_is_not_omitted_from_inventory_validation(tmp_path: Path):
    result = _module().generate_simulation_bundle(decode(simulation_protocol()), tmp_path / "bundle")
    (result.output / "undeclared").mkdir()
    with pytest.raises(ValueError, match="inventory"):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)


@pytest.mark.parametrize(
    "field,value",
    [
        ("generator", {}),
        ("pair_count", True),
        ("independent_group_count", True),
        ("generator", {"python_implementation": "CPython", "python_version": "3.12.13", "source_sha256": {}}),
    ],
)
def test_manifest_shape_is_checked_even_when_digest_is_supplied_by_caller(tmp_path: Path, field, value):
    raw = simulation_protocol()
    raw["cases"][0]["reads"]["pair_count"] = 1
    result = _module().generate_simulation_bundle(decode(raw), tmp_path / "bundle")
    path = result.output / "manifest.json"
    manifest = json.loads(path.read_bytes())
    manifest[field] = value
    data = json.dumps(manifest).encode()
    path.write_bytes(data)
    with pytest.raises(ValueError, match="manifest|generator"):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=hashlib.sha256(data).hexdigest())


def test_manifest_duplicate_keys_are_invalid_even_with_an_exact_digest(tmp_path: Path):
    result = _module().generate_simulation_bundle(decode(simulation_protocol()), tmp_path / "bundle")
    path = result.output / "manifest.json"
    data = path.read_bytes().replace(b"{", b'{"pair_count":3,', 1)
    path.write_bytes(data)
    with pytest.raises(ValueError):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=hashlib.sha256(data).hexdigest())


def test_file_mutation_during_streamed_inventory_hash_is_rejected(tmp_path: Path):
    module = _module()
    payload = tmp_path / "reads.fastq"
    payload.write_bytes(b"original")
    real_read = os.read
    changed = False

    def mutate(descriptor, count):
        nonlocal changed
        data = real_read(descriptor, count)
        if data and not changed:
            changed = True
            payload.write_bytes(b"changed and extended")
        return data

    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with patch.object(module.os, "read", side_effect=mutate), pytest.raises(ValueError, match="changed"):
            module._hash_child(parent, payload.name)
    finally:
        os.close(parent)
    assert changed


def test_stream_wrapper_failure_closes_owned_file_descriptor(tmp_path: Path):
    before = set(os.listdir("/proc/self/fd"))
    with patch.object(_module().os, "fdopen", side_effect=OSError("wrapper failed")), pytest.raises(OSError):
        _module()._private_file(tmp_path / "payload")
    assert set(os.listdir("/proc/self/fd")) == before


def test_symlink_case_directory_and_nonregular_payload_are_refused(tmp_path: Path):
    result = _module().generate_simulation_bundle(decode(simulation_protocol()), tmp_path / "bundle")
    folder = result.output / "case-a"
    outside = tmp_path / "case-a"
    folder.rename(outside)
    folder.symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)
    folder.unlink()
    outside.rename(folder)
    target = folder / "reads_R1.fastq"
    target.unlink()
    os.mkfifo(target)
    with pytest.raises(ValueError, match="nonregular"):
        _module().verify_generated_bundle(result.output, expected_manifest_sha256=result.manifest_sha256)


def test_linked_root_is_rejected_before_opening_manifest(tmp_path: Path):
    result = _module().generate_simulation_bundle(decode(simulation_protocol()), tmp_path / "bundle")
    linked = tmp_path / "link"
    linked.symlink_to(result.output, target_is_directory=True)
    with patch.object(_module(), "read_regular_path") as reader, pytest.raises(ValueError, match="root"):
        _module().verify_generated_bundle(linked, expected_manifest_sha256=result.manifest_sha256)
    reader.assert_not_called()
