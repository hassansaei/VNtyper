"""Strict no-follow loading of approved portable length model bundles."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256

pytestmark = pytest.mark.unit

_FILES = {
    "candidate.json",
    "payload-manifest.json",
    "portable-approval.json",
    "length-model.json",
    "length-annotation.json",
}


def _annotation() -> dict[str, object]:
    return {
        "schema_version": "length-annotation-v1",
        "assembly": "synthetic-build-v1",
        "contig": "synthetic-contig",
        "accepted_contig_aliases": ["alias-contig"],
        "reference_fasta_sha256": "a" * 64,
        "coordinate_system": "zero-based-half-open",
        "boundary_definition": "complete-core-plus-invariant-units-v1",
        "repeat_unit_bp": 2,
        "regions": {
            "CORE": [{"start": 10, "end": 14}],
            "INVARIANT": [{"start": 16, "end": 18}],
            "ARRAY": {"start": 10, "end": 19},
            "LEFT_FLANK": {"start": 7, "end": 10},
            "RIGHT_FLANK": {"start": 19, "end": 22},
        },
        "array_boundary_geometry": {"array_only_bp": 3, "target_only_bp": 0},
        "target_boundary_conversion_sha256": None,
        "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
        "annotation_provenance": "generated-synthetic-test-contract",
        "annotation_version": "synthetic-1",
    }


def _applicability() -> dict[str, object]:
    return {
        "domain": "synthetic",
        "assemblies": ["synthetic-build-v1"],
        "assay_classes": ["synthetic-short-read"],
        "input_scopes": ["regional"],
        "preprocessing_ids": ["synthetic-preprocessing-v1"],
        "aligner_name": "synthetic-aligner",
        "aligner_version": "1.0",
        "aligner_arguments_sha256": "c" * 64,
        "primary_secondary_marking": "primary-only",
        "counting_policy_sha256": "b" * 64,
    }


def _producer() -> dict[str, object]:
    return {
        "name": "synthetic-length-fit",
        "version": "1.0",
        "source_revision": "f" * 40,
        "tool_versions": {"numpy": "1.26"},
        "feature_schema_sha256": "0" * 64,
    }


def _model(annotation_sha256: str) -> dict[str, object]:
    return {
        "schema_version": "length-model-v1",
        "target": {
            "name": "total_diploid_repeat_count",
            "boundary_definition": "complete-core-plus-invariant-units-v1",
        },
        "unit": "repeat_units",
        "model_kind": "affine-A",
        "feature_order": ["A"],
        "intercept": 10.0,
        "coefficients": [50.0],
        "annotation_sha256": annotation_sha256,
        "counting_policy_sha256": "b" * 64,
        "applicability": _applicability(),
        "qc": {
            "minimum_denominator_mean_depth": 10.0,
            "minimum_denominator_covered_fraction": 0.9,
            "minimum_denominator_supporting_fragments": 100,
            "fragment_evidence_kind": "read-pair-identity-qc-proxy",
        },
        "feature_bounds": {"A": {"minimum": 0.5, "maximum": 3.0}},
        "study_sha256": "d" * 64,
        "training_evidence_sha256": "e" * 64,
        "producer": _producer(),
    }


def _write_bundle(root: Path) -> dict[str, dict[str, object] | list[dict[str, object]]]:
    root.mkdir(exist_ok=True)
    annotation = _annotation()
    annotation_raw = canonical_json_bytes(annotation)
    model = _model(hashlib.sha256(annotation_raw).hexdigest())
    model_raw = canonical_json_bytes(model)
    manifest = [
        {
            "path": "length-annotation.json",
            "size_bytes": len(annotation_raw),
            "sha256": hashlib.sha256(annotation_raw).hexdigest(),
        },
        {
            "path": "length-model.json",
            "size_bytes": len(model_raw),
            "sha256": hashlib.sha256(model_raw).hexdigest(),
        },
    ]
    candidate = {
        "schema_version": "calibration-candidate-v2",
        "target": "length",
        "study_sha256": model["study_sha256"],
        "baseline_sha256": "5" * 64,
        "partition_sha256": "4" * 64,
        "training_evidence_sha256": model["training_evidence_sha256"],
        "selection_evidence_sha256": "6" * 64,
        "payload_sha256": canonical_sha256(manifest),
        "applicability": _applicability(),
        "producer": _producer(),
        "status": "research-only",
    }
    candidate["candidate_id"] = canonical_sha256(candidate)
    candidate_sha256 = canonical_sha256(candidate)
    approval = {
        "schema_version": "calibration-portable-approval-v2",
        "target": "length",
        "disposition": "passed-validation-and-locked-heldout",
        "candidate_sha256": candidate_sha256,
        "candidate_id": candidate["candidate_id"],
        "payload_sha256": candidate["payload_sha256"],
        "study_sha256": candidate["study_sha256"],
        "protocol_sha256": "7" * 64,
        "partition_sha256": candidate["partition_sha256"],
        "baseline_sha256": candidate["baseline_sha256"],
        "exposure_ledger_id": "8" * 64,
        "applicability_sha256": canonical_sha256(_applicability()),
        "producer_sha256": canonical_sha256(_producer()),
        "validation_attestation_sha256": "9" * 64,
        "locked_heldout_attestation_sha256": "a" * 64,
        "custodian_authority_sha256": "b" * 64,
        "custody_completion_sha256": "c" * 64,
    }
    documents: dict[str, dict[str, object] | list[dict[str, object]]] = {
        "candidate.json": candidate,
        "payload-manifest.json": manifest,
        "portable-approval.json": approval,
        "length-model.json": model,
        "length-annotation.json": annotation,
    }
    for name, document in documents.items():
        (root / name).write_bytes(canonical_json_bytes(document))
    checksums = {
        "schema_version": "calibration-checksums-v1",
        "files": {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in sorted(_FILES)},
    }
    (root / "checksums.json").write_bytes(canonical_json_bytes(checksums))
    return documents


def _rewrite_checksums(root: Path) -> None:
    checksums = {
        "schema_version": "calibration-checksums-v1",
        "files": {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in sorted(_FILES)},
    }
    (root / "checksums.json").write_bytes(canonical_json_bytes(checksums))


def _resign_bundle(root: Path, documents: dict[str, dict[str, object] | list[dict[str, object]]]) -> None:
    payload_rows: list[dict[str, object]] = []
    for name in ("length-annotation.json", "length-model.json"):
        raw = canonical_json_bytes(documents[name])
        (root / name).write_bytes(raw)
        payload_rows.append({"path": name, "size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()})
    documents["payload-manifest.json"] = payload_rows
    candidate = documents["candidate.json"]
    approval = documents["portable-approval.json"]
    assert isinstance(candidate, dict) and isinstance(approval, dict)
    candidate["payload_sha256"] = canonical_sha256(payload_rows)
    candidate["candidate_id"] = canonical_sha256(
        {key: value for key, value in candidate.items() if key != "candidate_id"}
    )
    approval["payload_sha256"] = candidate["payload_sha256"]
    approval["candidate_id"] = candidate["candidate_id"]
    approval["candidate_sha256"] = canonical_sha256(candidate)
    for name in ("payload-manifest.json", "candidate.json", "portable-approval.json"):
        (root / name).write_bytes(canonical_json_bytes(documents[name]))
    _rewrite_checksums(root)


def test_runtime_length_bundle_loads_exact_approved_model_and_annotation(tmp_path: Path) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    _write_bundle(tmp_path)
    bundle = load_length_model_bundle(tmp_path)

    assert bundle.model.annotation_sha256 == bundle.annotation.sha256
    assert bundle.candidate.payload_sha256 == bundle.payload.sha256
    assert bundle.approval.candidate_sha256 == bundle.candidate.sha256
    assert bundle.sha256 == canonical_sha256(
        {
            "schema_version": "calibration-checksums-v1",
            "files": {name: hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() for name in sorted(_FILES)},
        }
    )
    assert not hasattr(bundle, "source_path")


@pytest.mark.parametrize("missing", sorted(_FILES | {"checksums.json"}))
def test_runtime_bundle_requires_exact_six_file_inventory(tmp_path: Path, missing: str) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    _write_bundle(tmp_path)
    (tmp_path / missing).unlink()
    with pytest.raises(ValueError, match="inventory"):
        load_length_model_bundle(tmp_path)


def test_runtime_bundle_rejects_extra_file_and_symlinked_entry(tmp_path: Path) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    _write_bundle(tmp_path)
    (tmp_path / "extra.json").write_text("{}")
    with pytest.raises(ValueError, match="inventory"):
        load_length_model_bundle(tmp_path)
    (tmp_path / "extra.json").unlink()
    (tmp_path / "length-model.json").unlink()
    (tmp_path / "length-model.json").symlink_to("candidate.json")
    with pytest.raises(ValueError, match="regular|symlink"):
        load_length_model_bundle(tmp_path)


@pytest.mark.parametrize("name", sorted(_FILES))
def test_runtime_bundle_rejects_each_changed_file_before_semantic_use(tmp_path: Path, name: str) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    _write_bundle(tmp_path)
    with (tmp_path / name).open("ab") as stream:
        stream.write(b" ")
    with pytest.raises(ValueError, match="checksum"):
        load_length_model_bundle(tmp_path)


def test_runtime_bundle_rejects_noncanonical_json_even_with_updated_checksum(tmp_path: Path) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    _write_bundle(tmp_path)
    raw = (tmp_path / "portable-approval.json").read_bytes()
    (tmp_path / "portable-approval.json").write_bytes(raw.rstrip() + b"  \n")
    _rewrite_checksums(tmp_path)
    with pytest.raises(ValueError, match="canonical"):
        load_length_model_bundle(tmp_path)


@pytest.mark.parametrize("mutation", ["noncanonical", "malformed", "nonfinite"])
def test_runtime_bundle_rejects_noncanonical_or_malformed_payload_manifest(tmp_path: Path, mutation: str) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    _write_bundle(tmp_path)
    raw = (tmp_path / "payload-manifest.json").read_bytes()
    if mutation == "noncanonical":
        changed = raw.rstrip() + b"  \n"
    elif mutation == "nonfinite":
        changed = b"[NaN]\n"
    else:
        changed = b"{\n"
    (tmp_path / "payload-manifest.json").write_bytes(changed)
    _rewrite_checksums(tmp_path)
    with pytest.raises(ValueError, match="payload manifest"):
        load_length_model_bundle(tmp_path)


@pytest.mark.parametrize("mutation", ["schema", "files-type", "inventory", "digest"])
def test_runtime_bundle_rejects_malformed_checksum_contract(tmp_path: Path, mutation: str) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    _write_bundle(tmp_path)
    checksums: dict[str, object] = {
        "schema_version": "calibration-checksums-v1",
        "files": {name: hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() for name in sorted(_FILES)},
    }
    if mutation == "schema":
        checksums["schema_version"] = "calibration-checksums-v2"
    elif mutation == "files-type":
        checksums["files"] = []
    elif mutation == "inventory":
        assert isinstance(checksums["files"], dict)
        checksums["files"].pop("candidate.json")
    else:
        assert isinstance(checksums["files"], dict)
        checksums["files"]["candidate.json"] = "F" * 64
    (tmp_path / "checksums.json").write_bytes(canonical_json_bytes(checksums))
    with pytest.raises(ValueError, match="checksum"):
        load_length_model_bundle(tmp_path)


def test_runtime_bundle_rejects_research_model_only_payload(tmp_path: Path) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    documents = _write_bundle(tmp_path)
    model_raw = (tmp_path / "length-model.json").read_bytes()
    manifest = [
        {
            "path": "length-model.json",
            "size_bytes": len(model_raw),
            "sha256": hashlib.sha256(model_raw).hexdigest(),
        }
    ]
    (tmp_path / "payload-manifest.json").write_bytes(canonical_json_bytes(manifest))
    candidate = documents["candidate.json"]
    assert isinstance(candidate, dict)
    candidate["payload_sha256"] = canonical_sha256(manifest)
    candidate["candidate_id"] = canonical_sha256(
        {key: value for key, value in candidate.items() if key != "candidate_id"}
    )
    (tmp_path / "candidate.json").write_bytes(canonical_json_bytes(candidate))
    approval = documents["portable-approval.json"]
    assert isinstance(approval, dict)
    approval["payload_sha256"] = candidate["payload_sha256"]
    approval["candidate_id"] = candidate["candidate_id"]
    approval["candidate_sha256"] = canonical_sha256(candidate)
    (tmp_path / "portable-approval.json").write_bytes(canonical_json_bytes(approval))
    _rewrite_checksums(tmp_path)

    with pytest.raises(ValueError, match="exactly.*annotation.*model|payload"):
        load_length_model_bundle(tmp_path)


@pytest.mark.parametrize("mismatch", ["annotation", "study", "training", "applicability", "producer"])
def test_runtime_bundle_rejects_semantic_model_candidate_mismatches(tmp_path: Path, mismatch: str) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    documents = _write_bundle(tmp_path)
    model = documents["length-model.json"]
    assert isinstance(model, dict)
    if mismatch == "annotation":
        model["annotation_sha256"] = "f" * 64
    elif mismatch == "study":
        model["study_sha256"] = "f" * 64
    elif mismatch == "training":
        model["training_evidence_sha256"] = "f" * 64
    elif mismatch == "applicability":
        model["applicability"] = {**_applicability(), "assay_classes": ["other"]}
    else:
        model["producer"] = {**_producer(), "version": "2.0"}
    _resign_bundle(tmp_path, documents)

    with pytest.raises(ValueError, match="model|annotation"):
        load_length_model_bundle(tmp_path)


@pytest.mark.parametrize("mismatch", ["boundary", "assembly"])
def test_runtime_bundle_rejects_annotation_outside_model_target_or_build(tmp_path: Path, mismatch: str) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    documents = _write_bundle(tmp_path)
    annotation = documents["length-annotation.json"]
    model = documents["length-model.json"]
    assert isinstance(annotation, dict) and isinstance(model, dict)
    if mismatch == "boundary":
        annotation["boundary_definition"] = "different-target-boundary"
    else:
        annotation["assembly"] = "different-build"
    model["annotation_sha256"] = canonical_sha256(annotation)
    _resign_bundle(tmp_path, documents)
    with pytest.raises(ValueError, match="boundary|assembly"):
        load_length_model_bundle(tmp_path)


def test_runtime_bundle_rejects_wrong_type_and_direct_bundle_path_symlink(tmp_path: Path) -> None:
    from vntyper.scripts.length_model_bundle import load_length_model_bundle

    with pytest.raises(ValueError):
        load_length_model_bundle("bundle")  # type: ignore[arg-type]
    bundle = tmp_path / "bundle"
    _write_bundle(bundle)
    linked = tmp_path / "linked"
    linked.symlink_to(bundle, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        load_length_model_bundle(linked)
