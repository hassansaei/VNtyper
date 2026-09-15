"""Strict pipeline configuration for optional VNTR length measurement."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import cast

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256
from vntyper.scripts.length_annotation import decode_length_annotation
from vntyper.scripts.length_feature_provenance import decode_length_feature_context
from vntyper.scripts.length_model_bundle import LengthModelBundle, load_length_model_bundle

pytestmark = pytest.mark.unit


def _annotation_raw() -> dict[str, object]:
    return {
        "schema_version": "length-annotation-v1",
        "assembly": "synthetic-build-v1",
        "contig": "synthetic-contig",
        "accepted_contig_aliases": ["alias-contig"],
        "reference_fasta_sha256": "2" * 64,
        "coordinate_system": "zero-based-half-open",
        "boundary_definition": "complete-core-plus-invariant-units-v1",
        "repeat_unit_bp": 2,
        "regions": {
            "CORE": [{"start": 2, "end": 4}],
            "INVARIANT": [{"start": 6, "end": 8}],
            "ARRAY": {"start": 2, "end": 9},
            "LEFT_FLANK": {"start": 0, "end": 2},
            "RIGHT_FLANK": {"start": 9, "end": 11},
        },
        "array_boundary_geometry": {"array_only_bp": 3, "target_only_bp": 0},
        "target_boundary_conversion_sha256": None,
        "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
        "annotation_provenance": "generated-synthetic-test-contract",
        "annotation_version": "synthetic-1",
    }


def _context_raw(annotation_sha256: str) -> dict[str, object]:
    return {
        "schema_version": "length-feature-measurement-context-v1",
        "manifest_key": "synthetic-member-1",
        "input_sha256": "1" * 64,
        "assembly": "synthetic-build-v1",
        "assay_class": "synthetic-short-read",
        "input_scope": "regional",
        "original_contig": "synthetic-contig",
        "reference_fasta_sha256": "2" * 64,
        "annotation_sha256": annotation_sha256,
        "aligner": {
            "name": "synthetic-aligner",
            "version": "1.0",
            "arguments_sha256": "4" * 64,
            "primary_secondary_marking": "primary-only",
        },
        "fragment_reader": {
            "name": "pysam",
            "version": "0.22.1",
            "htslib_version": "1.18",
            "alignment_semantics": "explicit-filtered-aligned-pairs-v1",
        },
        "preprocessing_id": "synthetic-preprocessing-v1",
        "counting_policy": {
            "policy_id": "primary-mapq0-baseq0-overlap-count-v1",
            "samtools_revision": "samtools=1.20;htslib=1.20",
            "minimum_mapping_quality": 0,
            "minimum_base_quality": 0,
            "excluded_alignment_flags": ["UNMAP", "SECONDARY", "QCFAIL", "DUP"],
            "supplementary_alignment_policy": "included-unless-excluded-by-another-flag",
            "overlap_policy": "count-overlapping-mates-independently",
            "base_counting_policy": "one-per-aligned-covered-base",
            "zero_coverage_policy": "emit-zero-for-every-queried-position",
            "queried_intervals": [{"start": 0, "end": 11}],
        },
    }


def _sidecar_raw() -> dict[str, object]:
    annotation = decode_length_annotation(_annotation_raw())
    return {
        "schema_version": "length-pipeline-context-v1",
        "evidence_domain": "synthetic",
        "measurement_context": _context_raw(annotation.sha256),
    }


def _write_applicable_bundle(root: Path) -> tuple[LengthModelBundle, dict[str, object]]:
    from tests.unit.test_length_model_bundle import _resign_bundle, _write_bundle

    documents = _write_bundle(root)
    annotation_raw = cast(dict[str, object], documents["length-annotation.json"])
    annotation = decode_length_annotation(annotation_raw)
    context_raw = _context_raw(annotation.sha256)
    context_raw.update(
        {
            "reference_fasta_sha256": annotation.reference_fasta_sha256,
            "aligner": {
                **cast(dict[str, object], context_raw["aligner"]),
                "arguments_sha256": "c" * 64,
            },
            "counting_policy": {
                **cast(dict[str, object], context_raw["counting_policy"]),
                "queried_intervals": [{"start": 7, "end": 22}],
            },
        }
    )
    context = decode_length_feature_context(context_raw)
    model = cast(dict[str, object], documents["length-model.json"])
    candidate = cast(dict[str, object], documents["candidate.json"])
    applicability = {
        **cast(dict[str, object], model["applicability"]),
        "counting_policy_sha256": context.counting_policy_sha256,
    }
    model["counting_policy_sha256"] = context.counting_policy_sha256
    model["applicability"] = applicability
    candidate["applicability"] = applicability
    approval = cast(dict[str, object], documents["portable-approval.json"])
    approval["applicability_sha256"] = canonical_sha256(applicability)
    _resign_bundle(root, documents)
    return load_length_model_bundle(root), {
        "schema_version": "length-pipeline-context-v1",
        "evidence_domain": "synthetic",
        "measurement_context": context_raw,
    }


def test_pipeline_context_is_closed_immutable_and_hash_bound() -> None:
    from vntyper.scripts.pipeline_length import decode_length_pipeline_context, encode_length_pipeline_context

    context = decode_length_pipeline_context(_sidecar_raw())

    assert encode_length_pipeline_context(context) == _sidecar_raw()
    assert context.sha256 == canonical_sha256(_sidecar_raw())
    assert context.evidence_domain == "synthetic"
    with pytest.raises(ValueError, match="digest"):
        encode_length_pipeline_context(replace(context, sha256="f" * 64))


@pytest.mark.parametrize(
    "raw",
    [
        {**_sidecar_raw(), "extra": True},
        {**_sidecar_raw(), "evidence_domain": "unknown"},
        {**_sidecar_raw(), "evidence_domain": True},
    ],
)
def test_pipeline_context_rejects_open_or_invalid_content(raw: dict[str, object]) -> None:
    from vntyper.scripts.pipeline_length import decode_length_pipeline_context

    with pytest.raises(ValueError):
        decode_length_pipeline_context(raw)


def test_measurement_only_configuration_binds_annotation_context_and_policy() -> None:
    from vntyper.scripts.pipeline_length import (
        build_length_pipeline_configuration,
        encode_length_pipeline_configuration,
    )

    annotation = decode_length_annotation(_annotation_raw())
    context = decode_length_feature_context(_context_raw(annotation.sha256))
    sidecar = __import__(
        "vntyper.scripts.pipeline_length", fromlist=["decode_length_pipeline_context"]
    ).decode_length_pipeline_context(_sidecar_raw())

    configuration = build_length_pipeline_configuration(annotation=annotation, pipeline_context=sidecar, bundle=None)
    encoded = encode_length_pipeline_configuration(configuration)

    assert configuration.measurement_enabled is True
    assert configuration.annotation == annotation
    assert configuration.measurement_context == context
    assert configuration.model is None
    assert encoded == {
        "schema_version": "length-pipeline-configuration-v1",
        "measurement_enabled": True,
        "annotation_sha256": annotation.sha256,
        "model_sha256": None,
        "model_bundle_sha256": None,
        "portable_approval_sha256": None,
        "candidate_id": None,
        "context_sha256": sidecar.sha256,
        "counting_policy_sha256": context.counting_policy_sha256,
        "evidence_domain": "synthetic",
    }
    assert configuration.sha256 == canonical_sha256(encoded)


def test_measurement_configuration_rejects_context_geometry_or_observable_applicability_drift() -> None:
    from vntyper.scripts.pipeline_length import build_length_pipeline_configuration, decode_length_pipeline_context

    annotation = decode_length_annotation(_annotation_raw())
    raw = _sidecar_raw()
    measurement = dict(cast(dict[str, object], raw["measurement_context"]))
    measurement["assembly"] = "different-build"
    raw["measurement_context"] = measurement
    context = decode_length_pipeline_context(raw)

    with pytest.raises(ValueError, match="assembly"):
        build_length_pipeline_configuration(annotation=annotation, pipeline_context=context, bundle=None)


def test_resolver_enforces_disabled_and_measurement_only_path_matrix(tmp_path: Path) -> None:
    from vntyper.scripts.pipeline_length import resolve_length_pipeline_configuration

    disabled = resolve_length_pipeline_configuration(
        measurement_enabled=False, model_path=None, annotation_path=None, context_path=None
    )
    assert disabled.measurement_enabled is False
    assert disabled.annotation is None
    assert disabled.sha256 == canonical_sha256(
        {
            "schema_version": "length-pipeline-configuration-v1",
            "measurement_enabled": False,
            "annotation_sha256": None,
            "model_sha256": None,
            "model_bundle_sha256": None,
            "portable_approval_sha256": None,
            "candidate_id": None,
            "context_sha256": None,
            "counting_policy_sha256": None,
            "evidence_domain": None,
        }
    )
    from vntyper.scripts.pipeline_length import encode_length_pipeline_configuration

    assert encode_length_pipeline_configuration(disabled)["measurement_enabled"] is False

    annotation = decode_length_annotation(_annotation_raw())
    annotation_path = tmp_path / "annotation.json"
    context_path = tmp_path / "context.json"
    annotation_path.write_bytes(canonical_json_bytes(_annotation_raw()))
    context_path.write_bytes(canonical_json_bytes(_sidecar_raw()))
    enabled = resolve_length_pipeline_configuration(
        measurement_enabled=True, model_path=None, annotation_path=annotation_path, context_path=context_path
    )
    assert enabled.annotation == annotation

    for model_path, annotation_value, context_value in (
        (None, None, context_path),
        (None, annotation_path, None),
        (None, annotation_path, context_path),
    ):
        with pytest.raises(ValueError):
            resolve_length_pipeline_configuration(
                measurement_enabled=False,
                model_path=model_path,
                annotation_path=annotation_value,
                context_path=context_value,
            )

    with pytest.raises(ValueError, match="length-context"):
        resolve_length_pipeline_configuration(
            measurement_enabled=True, model_path=None, annotation_path=annotation_path, context_path=None
        )
    with pytest.raises(ValueError, match="length-annotation"):
        resolve_length_pipeline_configuration(
            measurement_enabled=True, model_path=None, annotation_path=None, context_path=context_path
        )


def test_resolver_rejects_noncanonical_and_symlinked_sidecars(tmp_path: Path) -> None:
    from vntyper.scripts.pipeline_length import resolve_length_pipeline_configuration

    annotation_path = tmp_path / "annotation.json"
    annotation_path.write_bytes(canonical_json_bytes(_annotation_raw()))
    context_path = tmp_path / "context.json"
    context_path.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError):
        resolve_length_pipeline_configuration(
            measurement_enabled=True, model_path=None, annotation_path=annotation_path, context_path=context_path
        )
    context_path.unlink()
    target = tmp_path / "real-context.json"
    target.write_bytes(canonical_json_bytes(_sidecar_raw()))
    context_path.symlink_to(target)
    with pytest.raises(ValueError, match="symlink|unreadable"):
        resolve_length_pipeline_configuration(
            measurement_enabled=True, model_path=None, annotation_path=annotation_path, context_path=context_path
        )


def test_model_bundle_implies_measurement_and_supplies_annotation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from vntyper.scripts import pipeline_length

    annotation = decode_length_annotation(_annotation_raw())
    context_path = tmp_path / "context.json"
    context_path.write_bytes(canonical_json_bytes(_sidecar_raw()))
    bundle_path = tmp_path / "approved-bundle"
    bundle_path.mkdir()
    sentinel_bundle = SimpleNamespace(annotation=annotation)
    sentinel_configuration = object()
    observed: dict[str, object] = {}

    def fake_load(path: Path) -> object:
        observed["path"] = path
        return sentinel_bundle

    def fake_build(*, annotation: object, pipeline_context: object, bundle: object) -> object:
        observed.update(annotation=annotation, pipeline_context=pipeline_context, bundle=bundle)
        return sentinel_configuration

    monkeypatch.setattr(pipeline_length, "load_length_model_bundle", fake_load)
    monkeypatch.setattr(pipeline_length, "build_length_pipeline_configuration", fake_build)

    result = pipeline_length.resolve_length_pipeline_configuration(
        measurement_enabled=False,
        model_path=bundle_path,
        annotation_path=None,
        context_path=context_path,
    )

    assert result is sentinel_configuration
    assert observed["path"] == bundle_path
    assert observed["annotation"] == annotation
    assert observed["bundle"] is sentinel_bundle


def test_model_bundle_rejects_separate_annotation_before_loading_files(tmp_path: Path) -> None:
    from vntyper.scripts.pipeline_length import resolve_length_pipeline_configuration

    with pytest.raises(ValueError, match="supplies its annotation"):
        resolve_length_pipeline_configuration(
            measurement_enabled=False,
            model_path=tmp_path / "bundle",
            annotation_path=tmp_path / "annotation.json",
            context_path=tmp_path / "context.json",
        )


def test_real_approved_bundle_resolves_model_and_revalidates_configuration(tmp_path: Path) -> None:
    from vntyper.scripts.pipeline_length import (
        encode_length_pipeline_configuration,
        resolve_length_pipeline_configuration,
    )

    bundle, sidecar = _write_applicable_bundle(tmp_path / "bundle")
    context_path = tmp_path / "context.json"
    context_path.write_bytes(canonical_json_bytes(sidecar))

    configuration = resolve_length_pipeline_configuration(
        measurement_enabled=False,
        model_path=tmp_path / "bundle",
        annotation_path=None,
        context_path=context_path,
    )
    encoded = encode_length_pipeline_configuration(configuration)

    assert configuration.measurement_enabled is True
    assert configuration.model == bundle.model
    assert configuration.model_bundle_sha256 == bundle.sha256
    assert configuration.portable_approval_sha256 == bundle.approval.sha256
    assert encoded["model_sha256"] == bundle.model.sha256


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("evidence_domain", "external", "applicability"),
        ("assay_class", "different-assay", "applicability"),
        ("preprocessing_id", "different-preprocessing", "applicability"),
    ],
)
def test_model_mode_rejects_observable_applicability_drift(
    tmp_path: Path, field: str, value: str, message: str
) -> None:
    from vntyper.scripts.pipeline_length import build_length_pipeline_configuration, decode_length_pipeline_context

    bundle, sidecar = _write_applicable_bundle(tmp_path / "bundle")
    if field == "evidence_domain":
        sidecar[field] = value
    else:
        nested = dict(cast(dict[str, object], sidecar["measurement_context"]))
        nested[field] = value
        sidecar["measurement_context"] = nested
    context = decode_length_pipeline_context(sidecar)

    with pytest.raises(ValueError, match=message):
        build_length_pipeline_configuration(
            annotation=bundle.annotation,
            pipeline_context=context,
            bundle=bundle,
        )


def test_configuration_encoder_rejects_directly_replaced_state(tmp_path: Path) -> None:
    from vntyper.scripts.pipeline_length import (
        build_length_pipeline_configuration,
        decode_length_pipeline_context,
        encode_length_pipeline_configuration,
    )

    annotation = decode_length_annotation(_annotation_raw())
    configuration = build_length_pipeline_configuration(
        annotation=annotation,
        pipeline_context=decode_length_pipeline_context(_sidecar_raw()),
        bundle=None,
    )

    with pytest.raises(ValueError, match="approval"):
        encode_length_pipeline_configuration(replace(configuration, candidate_id="invented"))
    with pytest.raises(ValueError, match="digest"):
        encode_length_pipeline_configuration(replace(configuration, sha256="f" * 64))
    with pytest.raises(ValueError, match="typed"):
        encode_length_pipeline_configuration(True)  # type: ignore[arg-type]


def test_resolver_rejects_nonboolean_flag_and_nonpath_values() -> None:
    from vntyper.scripts.pipeline_length import resolve_length_pipeline_configuration

    with pytest.raises(ValueError, match="boolean"):
        resolve_length_pipeline_configuration(
            measurement_enabled=cast(bool, 1),
            model_path=None,
            annotation_path=None,
            context_path=None,
        )
    with pytest.raises(ValueError, match="model path"):
        resolve_length_pipeline_configuration(
            measurement_enabled=True,
            model_path="bundle",  # type: ignore[arg-type]
            annotation_path=None,
            context_path=Path("context.json"),
        )
