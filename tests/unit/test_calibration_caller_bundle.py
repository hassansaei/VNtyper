"""Closed portable caller bundles admit caller-v2 profiles only after approval."""

import hashlib
from importlib import import_module

import pytest

from tests.unit.test_calibration_export import _dual_setup, _setup
from vntyper.scripts.calibration_export import export_calibration_bundle
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

pytestmark = pytest.mark.unit


def loader():
    return import_module("vntyper.scripts.calibration_caller_bundle")


def bundle(root, mode="kestrel"):
    args = _setup(root, "callers")[0] if mode == "kestrel" else _dual_setup(root, mode)
    export_calibration_bundle(args, args.output)
    return args.output


@pytest.mark.parametrize("mode", ["kestrel", "legacy", "exact"])
def test_load_exported_inventory_preserves_resolved_policy_and_assets(mode, tmp_path):
    root = bundle(tmp_path, mode)
    result = loader().load_caller_model_bundle(root)
    assert result.candidate.target == "callers"
    assert result.approval.candidate_sha256 == result.candidate.sha256
    assert result.payload.sha256 == result.candidate.payload_sha256
    assert result.profile.profile_kind == "generated"
    assert result.caller_policy.required_callers == (("kestrel",) if mode == "kestrel" else ("advntr", "kestrel"))
    assert (result.advntr_policy is None) is (mode == "kestrel")
    assert (result.background_bytes is None) is (mode != "exact")
    if mode == "exact":
        assert result.advntr_policy.background_sha256 == hashlib.sha256(result.background_bytes).hexdigest()
    assert result.sha256 == hashlib.sha256((root / "checksums.json").read_bytes()).hexdigest()
    assert not hasattr(result, "source_path")


@pytest.mark.parametrize("mode", ["kestrel", "legacy", "exact"])
@pytest.mark.parametrize("change", ["missing", "extra", "symlink", "tampered", "checksum-inventory", "checksum-hash"])
def test_runtime_loader_refuses_partial_open_or_changed_artifacts(tmp_path, mode, change):
    root = bundle(tmp_path, mode)
    path = root / "candidate.json"
    if change == "missing":
        path.unlink()
    elif change == "extra":
        (root / "training.json").write_text("{}")
    elif change == "symlink":
        real = tmp_path / "real-candidate.json"
        path.rename(real)
        path.symlink_to(real)
    elif change == "tampered":
        path.write_bytes(path.read_bytes() + b" ")
    else:
        checksums_path = root / "checksums.json"
        document = load_strict_json_object(checksums_path.read_bytes())
        if change == "checksum-inventory":
            document["files"]["evidence.json"] = "a" * 64
        else:
            document["files"]["candidate.json"] = "not-a-hash"
        checksums_path.write_bytes(canonical_json_bytes(document))
    with pytest.raises(ValueError):
        loader().load_caller_model_bundle(root)


@pytest.mark.parametrize(
    "name", ["candidate.json", "payload-manifest.json", "portable-approval.json", "checksums.json"]
)
def test_forged_outer_checksums_do_not_relax_canonical_documents(tmp_path, name):
    root = bundle(tmp_path)
    path = root / name
    path.write_bytes(b" " + path.read_bytes())
    if name != "checksums.json":
        checksums_path = root / "checksums.json"
        checksums = load_strict_json_object(checksums_path.read_bytes())
        checksums["files"][name] = hashlib.sha256(path.read_bytes()).hexdigest()
        checksums_path.write_bytes(canonical_json_bytes(checksums))
    with pytest.raises(ValueError, match="canonical"):
        loader().load_caller_model_bundle(root)


def test_research_profile_and_symlinked_bundle_cannot_activate(tmp_path):
    args, _ = _setup(tmp_path, "callers")
    with pytest.raises(ValueError):
        loader().load_caller_model_bundle(args.profile)
    export_calibration_bundle(args, args.output)
    alias = tmp_path / "alias"
    alias.symlink_to(args.output, target_is_directory=True)
    with pytest.raises(ValueError):
        loader().load_caller_model_bundle(alias)


def test_outer_checksum_cannot_authorize_a_different_candidate(tmp_path):
    root = bundle(tmp_path)
    candidate = load_strict_json_object((root / "candidate.json").read_bytes())
    candidate["payload_sha256"] = "0" * 64
    (root / "candidate.json").write_bytes(canonical_json_bytes(candidate))
    checksums = load_strict_json_object((root / "checksums.json").read_bytes())
    checksums["files"]["candidate.json"] = hashlib.sha256((root / "candidate.json").read_bytes()).hexdigest()
    (root / "checksums.json").write_bytes(canonical_json_bytes(checksums))
    with pytest.raises(ValueError):
        loader().load_caller_model_bundle(root)


def test_runtime_loader_requires_path_and_existing_directory(tmp_path):
    with pytest.raises(ValueError, match="Path"):
        loader().load_caller_model_bundle(str(tmp_path))
    with pytest.raises(ValueError, match="missing"):
        loader().load_caller_model_bundle(tmp_path / "missing")


@pytest.mark.parametrize("change", ["schema", "fields", "files"])
def test_closed_checksum_schema(tmp_path, change):
    root = bundle(tmp_path)
    path = root / "checksums.json"
    checksums = load_strict_json_object(path.read_bytes())
    if change == "schema":
        checksums["schema_version"] = "other"
    elif change == "fields":
        checksums["extra"] = True
    else:
        checksums["files"] = []
    path.write_bytes(canonical_json_bytes(checksums))
    with pytest.raises(ValueError, match="checksum"):
        loader().load_caller_model_bundle(root)


@pytest.mark.parametrize("mode", ["kestrel", "legacy", "exact"])
def test_bundle_can_snapshot_exact_frozen_bytes_without_reopening_source(tmp_path, mode):
    root = bundle(tmp_path, mode)
    decoded = loader().load_caller_model_bundle(root)
    expected = {path.name: path.read_bytes() for path in root.iterdir()}
    for path in root.iterdir():
        path.unlink()
    assert loader().caller_model_bundle_files(decoded) == expected


def test_bundle_snapshot_rejects_replaced_digest_or_background(tmp_path):
    from dataclasses import replace

    decoded = loader().load_caller_model_bundle(bundle(tmp_path, "exact"))
    for changed in (replace(decoded, sha256="0" * 64), replace(decoded, background_bytes=b"changed")):
        with pytest.raises(ValueError):
            loader().caller_model_bundle_files(changed)
