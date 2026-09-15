"""Portable export verifies opened payload bytes and exports no research outcomes."""

from argparse import Namespace
from importlib import import_module
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.unit import test_calibration_portable_projection as approval_fixtures
from tests.unit.test_calibration_caller_profile import _build
from tests.unit.test_length_model_bundle import _write_bundle
from vntyper.scripts.calibration_candidate import candidate_document, decode_candidate
from vntyper.scripts.calibration_target_attestation import (
    target_locked_attestation_document,
    target_validation_attestation_document,
)
from vntyper.scripts.calibration_target_authority import target_custodian_authority_document
from vntyper.scripts.calibration_target_completion import target_completion_document
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.length_model_bundle import load_length_model_bundle

pytestmark = pytest.mark.unit


def exporter():
    return import_module("vntyper.scripts.calibration_export")


def _evidence(root, candidate):
    with patch.object(approval_fixtures, "_candidate", return_value=candidate):
        _, validation, locked, authority, completion = approval_fixtures._artifacts(target=candidate.target)
    docs = {
        "validation": target_validation_attestation_document(validation),
        "evaluation": target_locked_attestation_document(locked),
        "authority": target_custodian_authority_document(authority),
        "completion": target_completion_document(completion),
    }
    for name, doc in docs.items():
        (root / f"{name}.json").write_bytes(canonical_json_bytes(doc))
    return {name: root / f"{name}.json" for name in docs}


def _setup(root: Path, target="length"):
    profile = root / "research"
    profile.mkdir()
    payload = profile / "payload"
    payload.mkdir()
    if target == "length":
        docs = _write_bundle(root / "template")
        candidate = decode_candidate(docs["candidate.json"])
        manifest = docs["payload-manifest.json"]
        for name in ("length-model.json", "length-annotation.json"):
            (payload / name).write_bytes(canonical_json_bytes(docs[name]))
    else:
        generated = _build()
        (payload / "decision-profile.json").write_bytes(generated.canonical_bytes)
        descriptor = {
            "schema_version": "caller-bundle-v2",
            "required_callers": ["kestrel"],
            "components": {
                "decision-profile.json": generated.digest,
                "advntr-policy.json": None,
                "background.json": None,
            },
        }
        (payload / "caller-bundle.json").write_bytes(canonical_json_bytes(descriptor))
        import hashlib

        manifest = [
            {"path": p.name, "size_bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(payload.iterdir())
        ]
        raw = candidate_document(approval_fixtures._candidate(target="callers"))
        raw["payload_sha256"] = canonical_sha256(manifest)
        raw["partition_sha256"] = "b" * 64
        raw["candidate_id"] = canonical_sha256({k: v for k, v in raw.items() if k != "candidate_id"})
        candidate = decode_candidate(raw)
    (profile / "candidate.json").write_bytes(canonical_json_bytes(candidate_document(candidate)))
    (profile / "payload-manifest.json").write_bytes(canonical_json_bytes(manifest))
    # A research profile is deliberately richer than its portable projection.
    (profile / "study.json").write_text('{"private_marker":"never-export-this"}')
    (profile / "report.html").write_text("never-export-this")
    (profile / "checksums.json").write_text('{"files":{"anything":"untrusted-claim"}}')
    args = Namespace(target=target, profile=profile, output=root / "portable", **_evidence(root, candidate))
    args.output.mkdir()
    return args, candidate


@pytest.mark.parametrize("target", ["length", "callers"])
def test_export_opens_exact_payload_and_writes_only_approved_portable_files(tmp_path, target):
    args, candidate = _setup(tmp_path, target)
    assert exporter().export_calibration_bundle(args, args.output) is True
    payload_names = {p.name for p in (args.profile / "payload").iterdir()}
    assert {p.name for p in args.output.iterdir()} == payload_names | {
        "candidate.json",
        "payload-manifest.json",
        "portable-approval.json",
        "checksums.json",
    }
    assert all(b"never-export-this" not in p.read_bytes() for p in args.output.iterdir())
    approval = load_strict_json_object((args.output / "portable-approval.json").read_bytes())
    assert approval["candidate_sha256"] == candidate.sha256
    assert "custodian_name" not in approval and "metrics" not in approval
    for name in payload_names:
        assert (args.output / name).read_bytes() == (args.profile / "payload" / name).read_bytes()
    if target == "length":
        assert load_length_model_bundle(args.output).candidate == candidate
    else:
        from vntyper.scripts.decision_profile import load_packaged_decision_profile, parse_decision_profile

        with pytest.raises(ValueError):
            parse_decision_profile(
                (args.output / "decision-profile.json").read_bytes(),
                packaged_document=load_packaged_decision_profile().document,
            )


@pytest.mark.parametrize("field", ["validation", "evaluation"])
def test_unpassed_evidence_cannot_export(tmp_path, field):
    args, _ = _setup(tmp_path)
    path = getattr(args, field)
    doc = load_strict_json_object(path.read_bytes())
    doc["status"] = "failed"
    path.write_bytes(canonical_json_bytes(doc))
    with pytest.raises(ValueError, match="passed"):
        exporter().export_calibration_bundle(args, args.output)
    assert not list(args.output.iterdir())


@pytest.mark.parametrize("field", ["validation", "evaluation", "authority", "completion"])
def test_wrong_candidate_lineage_fails_before_writes(tmp_path, field):
    args, _ = _setup(tmp_path)
    path = getattr(args, field)
    doc = load_strict_json_object(path.read_bytes())
    doc["candidate_sha256"] = "0" * 64
    path.write_bytes(canonical_json_bytes(doc))
    with pytest.raises(ValueError):
        exporter().export_calibration_bundle(args, args.output)
    assert not list(args.output.iterdir())


@pytest.mark.parametrize(
    "change",
    ["payload-byte", "manifest-byte", "extra-payload", "missing-payload", "wrong-target", "unsupported-target"],
)
def test_current_payload_bytes_win_over_fake_research_checksums(tmp_path, change):
    args, _ = _setup(tmp_path)
    if change == "payload-byte":
        p = args.profile / "payload" / "length-model.json"
        p.write_bytes(p.read_bytes() + b" ")
    elif change == "manifest-byte":
        p = args.profile / "payload-manifest.json"
        p.write_bytes(p.read_bytes().replace(b'"size_bytes":', b'"size_bytes":1', 1))
    elif change == "extra-payload":
        (args.profile / "payload" / "evidence.json").write_text("{}")
    elif change == "missing-payload":
        (args.profile / "payload" / "length-model.json").unlink()
    else:
        args.target = "callers" if change == "wrong-target" else "dominance"
    with pytest.raises(ValueError):
        exporter().export_calibration_bundle(args, args.output)
    assert not list(args.output.iterdir())


@pytest.mark.parametrize("where", ["profile", "payload", "candidate.json", "payload/length-model.json", "authority"])
def test_export_rejects_symlinked_boundary_assets(tmp_path, where):
    args, _ = _setup(tmp_path)
    path = args.profile if where == "profile" else args.profile / where
    if where == "authority":
        path = args.authority
    if where == "payload":
        path = args.profile / "payload"
    moved = path.with_name(path.name + "-real")
    path.rename(moved)
    path.symlink_to(moved, target_is_directory=moved.is_dir())
    with pytest.raises(ValueError):
        exporter().export_calibration_bundle(args, args.output)
    assert not list(args.output.iterdir())


def test_caller_descriptor_composition_cannot_disagree_with_candidate(tmp_path):
    args, _ = _setup(tmp_path, "callers")
    path = args.profile / "payload" / "caller-bundle.json"
    raw = load_strict_json_object(path.read_bytes())
    raw["required_callers"] = ["advntr", "kestrel"]
    path.write_bytes(canonical_json_bytes(raw))
    with pytest.raises(ValueError):
        exporter().export_calibration_bundle(args, args.output)


def test_nonempty_staging_directory_is_not_overwritten(tmp_path):
    args, _ = _setup(tmp_path)
    protected = args.output / "candidate.json"
    protected.write_text("protected")
    with pytest.raises(ValueError):
        exporter().export_calibration_bundle(args, args.output)
    assert protected.read_text() == "protected"


def _rebind_payload(args):
    import hashlib

    manifest = [
        {"path": p.name, "size_bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
        for p in sorted((args.profile / "payload").iterdir())
    ]
    raw = load_strict_json_object((args.profile / "candidate.json").read_bytes())
    raw["payload_sha256"] = canonical_sha256(manifest)
    raw["candidate_id"] = canonical_sha256({k: v for k, v in raw.items() if k != "candidate_id"})
    candidate = decode_candidate(raw)
    (args.profile / "candidate.json").write_bytes(canonical_json_bytes(raw))
    (args.profile / "payload-manifest.json").write_bytes(canonical_json_bytes(manifest))
    _evidence(args.profile.parent, candidate)


@pytest.mark.parametrize("change", ["length-binding", "caller-partition", "packaged-profile", "descriptor-hash"])
def test_fully_rebound_bytes_still_require_target_payload_semantics_and_atomic_rollback(tmp_path, change):
    from vntyper.scripts.calibration_atomic_io import atomic_output
    from vntyper.scripts.decision_profile import load_packaged_decision_profile

    target = "length" if change == "length-binding" else "callers"
    args, _ = _setup(tmp_path, target)
    if change == "length-binding":
        path = args.profile / "payload" / "length-model.json"
        raw = load_strict_json_object(path.read_bytes())
        raw["study_sha256"] = "0" * 64
        path.write_bytes(canonical_json_bytes(raw))
    else:
        path = args.profile / "payload" / "decision-profile.json"
        raw = load_strict_json_object(path.read_bytes())
        if change == "caller-partition":
            raw["generated_metadata"]["partition_manifest_hash"] = "0" * 64
            path.write_bytes(canonical_json_bytes(raw))
        elif change == "packaged-profile":
            path.write_bytes(load_packaged_decision_profile().canonical_bytes)
        descriptor_path = args.profile / "payload" / "caller-bundle.json"
        descriptor = load_strict_json_object(descriptor_path.read_bytes())
        import hashlib

        descriptor["components"]["decision-profile.json"] = (
            "0" * 64 if change == "descriptor-hash" else hashlib.sha256(path.read_bytes()).hexdigest()
        )
        descriptor_path.write_bytes(canonical_json_bytes(descriptor))
    _rebind_payload(args)
    args.output.rmdir()
    with pytest.raises(ValueError):
        atomic_output(args.output, lambda staging: exporter().export_calibration_bundle(args, staging))
    assert not args.output.exists()


@pytest.mark.parametrize("name", ["candidate.json", "payload-manifest.json"])
def test_noncanonical_headers_fail_instead_of_being_silently_rehashed(tmp_path, name):
    args, _ = _setup(tmp_path)
    path = args.profile / name
    path.write_bytes(b" " + path.read_bytes())
    with pytest.raises(ValueError, match="canonical"):
        exporter().export_calibration_bundle(args, args.output)


@pytest.mark.parametrize("name", ["candidate.json", "payload-manifest.json"])
def test_duplicate_header_fields_fail_even_before_manifest_validation(tmp_path, name):
    args, _ = _setup(tmp_path)
    path = args.profile / name
    raw = path.read_bytes()
    path.write_bytes(
        raw.replace(b'"sha256":', b'"sha256":"duplicate","sha256":', 1)
        if name.startswith("payload")
        else raw.replace(b'"target":', b'"target":"length","target":', 1)
    )
    with pytest.raises(ValueError, match="duplicate"):
        exporter().export_calibration_bundle(args, args.output)


def test_profile_root_is_pinned_across_header_and_payload_opens(tmp_path, monkeypatch):
    args, candidate = _setup(tmp_path)
    module = exporter()
    original = module.SecureDirectoryReader.read_files
    moved = tmp_path / "moved-research"

    def swap(reader, names):
        raw = original(reader, names)
        if reader.path == args.profile and "candidate.json" in names:
            args.profile.rename(moved)
            args.profile.mkdir()
            (args.profile / "payload").mkdir()
            (args.profile / "payload" / "injected.txt").write_text("not-the-approved-payload")
        return raw

    monkeypatch.setattr(module.SecureDirectoryReader, "read_files", swap)
    assert module.export_calibration_bundle(args, args.output)
    assert (
        load_strict_json_object((args.output / "candidate.json").read_bytes())["candidate_id"] == candidate.candidate_id
    )
    assert not (args.output / "injected.txt").exists()


@pytest.mark.parametrize("field", ["profile", "validation", "evaluation", "authority", "completion"])
def test_missing_required_paths_are_explicit_failures(tmp_path, field):
    args, _ = _setup(tmp_path)
    setattr(args, field, None)
    with pytest.raises(ValueError, match="explicit path"):
        exporter().export_calibration_bundle(args, args.output)


def test_missing_headers_and_fifo_payload_fail_without_blocking(tmp_path):
    import os

    args, _ = _setup(tmp_path)
    header = args.profile / "candidate.json"
    saved = header.read_bytes()
    header.unlink()
    with pytest.raises(ValueError, match="headers"):
        exporter().export_calibration_bundle(args, args.output)
    header.write_bytes(saved)
    path = args.profile / "payload" / "length-model.json"
    path.unlink()
    os.mkfifo(path)
    with pytest.raises(ValueError, match="regular"):
        exporter().export_calibration_bundle(args, args.output)


def test_nofollow_capability_is_required(tmp_path, monkeypatch):
    args, _ = _setup(tmp_path)
    monkeypatch.setattr(exporter().os, "O_NOFOLLOW", 0)
    with pytest.raises(ValueError, match="no-follow"):
        exporter().export_calibration_bundle(args, args.output)


def _dual_setup(root, mode="exact"):
    import hashlib

    from tests.unit.test_calibration_advntr_runtime_policy import _caller
    from tests.unit.test_calibration_portable_background import background
    from vntyper.scripts.calibration_advntr_runtime_policy import (
        advntr_runtime_policy_document,
        build_advntr_runtime_policy,
    )
    from vntyper.scripts.calibration_caller_profile import build_caller_generated_profile
    from vntyper.scripts.calibration_portable_background import project_portable_background

    args, _ = _setup(root, "callers")
    caller = _caller(mode=mode)
    payload = args.profile / "payload"
    generated = build_caller_generated_profile(
        caller, dataset_manifest_hash="a" * 64, partition_manifest_hash="b" * 64, seed=17, generator_version="synthetic"
    )
    (payload / "decision-profile.json").write_bytes(generated.canonical_bytes)
    background_sha = None
    if mode == "exact":
        background_bytes = canonical_json_bytes(project_portable_background(background()))
        (payload / "background.json").write_bytes(background_bytes)
        background_sha = hashlib.sha256(background_bytes).hexdigest()
    sidecar = build_advntr_runtime_policy(
        caller,
        model_sha256="a" * 64,
        background_sha256=background_sha,
        capture_policy_sha256="c" * 64,
        advntr_revision="d" * 40,
    )
    sidecar_bytes = canonical_json_bytes(advntr_runtime_policy_document(sidecar))
    (payload / "advntr-policy.json").write_bytes(sidecar_bytes)
    descriptor = {
        "schema_version": "caller-bundle-v2",
        "required_callers": ["advntr", "kestrel"],
        "components": {
            "decision-profile.json": generated.digest,
            "advntr-policy.json": hashlib.sha256(sidecar_bytes).hexdigest(),
            "background.json": background_sha,
        },
    }
    (payload / "caller-bundle.json").write_bytes(canonical_json_bytes(descriptor))
    raw = load_strict_json_object((args.profile / "candidate.json").read_bytes())
    raw["applicability"]["required_callers"] = ["advntr", "kestrel"]
    raw["candidate_id"] = canonical_sha256({k: v for k, v in raw.items() if k != "candidate_id"})
    (args.profile / "candidate.json").write_bytes(canonical_json_bytes(raw))
    _rebind_payload(args)
    return args


@pytest.mark.parametrize("mode", ["legacy", "exact"])
def test_dual_caller_export_binds_sidecar_values_and_conditional_background(tmp_path, mode):
    args = _dual_setup(tmp_path, mode)
    assert exporter().export_calibration_bundle(args, args.output)
    sidecar = load_strict_json_object((args.output / "advntr-policy.json").read_bytes())
    assert sidecar["mode"] == mode
    assert (args.output / "background.json").exists() is (mode == "exact")
    assert {p.name for p in args.output.iterdir()} == {
        "candidate.json",
        "payload-manifest.json",
        "portable-approval.json",
        "checksums.json",
        "caller-bundle.json",
        "decision-profile.json",
        "advntr-policy.json",
    } | ({"background.json"} if mode == "exact" else set())


@pytest.mark.parametrize("change", ["sidecar-values", "sidecar-background", "background-narrative", "background-rate"])
def test_rebound_dual_caller_payload_still_checks_semantics_and_portable_provenance(tmp_path, change):
    import hashlib

    args = _dual_setup(tmp_path)
    payload = args.profile / "payload"
    sidecar = load_strict_json_object((payload / "advntr-policy.json").read_bytes())
    background = load_strict_json_object((payload / "background.json").read_bytes())
    if change == "sidecar-values":
        sidecar["minimum_read_support"] = 100
    elif change == "sidecar-background":
        sidecar["background_sha256"] = "0" * 64
    elif change == "background-narrative":
        background["provenance"] = "invented private narrative that must never leave research"
    else:
        background["default_probability"] = 0
    background_bytes = canonical_json_bytes(background)
    (payload / "background.json").write_bytes(background_bytes)
    if change != "sidecar-background":
        sidecar["background_sha256"] = hashlib.sha256(background_bytes).hexdigest()
    sidecar_bytes = canonical_json_bytes(sidecar)
    (payload / "advntr-policy.json").write_bytes(sidecar_bytes)
    descriptor = load_strict_json_object((payload / "caller-bundle.json").read_bytes())
    descriptor["components"]["background.json"] = hashlib.sha256(background_bytes).hexdigest()
    descriptor["components"]["advntr-policy.json"] = hashlib.sha256(sidecar_bytes).hexdigest()
    (payload / "caller-bundle.json").write_bytes(canonical_json_bytes(descriptor))
    _rebind_payload(args)
    with pytest.raises(ValueError):
        exporter().export_calibration_bundle(args, args.output)
    assert not list(args.output.iterdir())


def test_even_a_rebound_length_manifest_cannot_add_research_files(tmp_path):
    args, _ = _setup(tmp_path)
    (args.profile / "payload" / "training.json").write_text("{}")
    _rebind_payload(args)
    with pytest.raises(ValueError, match="exactly model and annotation"):
        exporter().export_calibration_bundle(args, args.output)
    assert not list(args.output.iterdir())


def test_legacy_generated_decision_profile_is_not_a_caller_v2_payload(tmp_path):
    import hashlib

    from vntyper.scripts.calibration_profiles import build_generated_profile
    from vntyper.scripts.decision_profile import load_packaged_decision_profile

    args, _ = _setup(tmp_path, "callers")
    packaged = load_packaged_decision_profile()
    generated = build_generated_profile(
        packaged.components["dominance"],
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=17,
        objective="lexicographic-safety-v1",
        generator_version="invented-test",
    )
    payload = args.profile / "payload"
    (payload / "decision-profile.json").write_bytes(generated.canonical_bytes)
    descriptor = load_strict_json_object((payload / "caller-bundle.json").read_bytes())
    descriptor["components"]["decision-profile.json"] = hashlib.sha256(generated.canonical_bytes).hexdigest()
    (payload / "caller-bundle.json").write_bytes(canonical_json_bytes(descriptor))
    _rebind_payload(args)
    with pytest.raises(ValueError, match="schema-v2"):
        exporter().export_calibration_bundle(args, args.output)
    assert not list(args.output.iterdir())
