"""A payload digest binds every file and a unique, safe layout."""

from copy import deepcopy
from dataclasses import FrozenInstanceError
from importlib import import_module

import pytest

from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def manifest_rows():
    return [
        {"path": "background.json", "size_bytes": 11, "sha256": "a" * 64},
        {"path": "decision-profile.json", "size_bytes": 23, "sha256": "b" * 64},
    ]


def test_manifest_binds_all_files_and_is_immutable():
    payload = import_module("vntyper.scripts.calibration_payload")
    raw = manifest_rows()
    manifest = payload.decode_payload_manifest(raw)

    assert manifest.sha256 == canonical_sha256(raw)
    assert tuple(item.path for item in manifest.files) == ("background.json", "decision-profile.json")
    assert manifest.files[0].size_bytes == 11
    raw[0]["sha256"] = "c" * 64
    assert manifest.files[0].sha256 == "a" * 64
    assert payload.payload_manifest_document(manifest) == manifest_rows()
    with pytest.raises(FrozenInstanceError):
        manifest.files[0].size_bytes = 100


@pytest.mark.parametrize(
    "path",
    [
        "",
        "/model.json",
        "../model",
        "a/../model",
        "./model",
        "a//model",
        "a/",
        "a\\model",
        "C:model",
        "a\x00b",
        "a\nb",
        "payload-manifest.json",
    ],
)
def test_manifest_rejects_unsafe_noncanonical_or_self_referential_paths(path):
    payload = import_module("vntyper.scripts.calibration_payload")
    raw = [{"path": path, "size_bytes": 1, "sha256": "a" * 64}]
    with pytest.raises(ValueError, match="path"):
        payload.decode_payload_manifest(raw)


@pytest.mark.parametrize("value", [None, {}, (), [], "manifest"])
def test_manifest_requires_a_nonempty_json_array(value):
    payload = import_module("vntyper.scripts.calibration_payload")
    with pytest.raises(ValueError, match="non-empty list"):
        payload.decode_payload_manifest(value)


@pytest.mark.parametrize(
    "field,value",
    [
        ("size_bytes", True),
        ("size_bytes", -1),
        ("size_bytes", 1.5),
        ("sha256", "A" * 64),
        ("sha256", "g" * 64),
        ("sha256", None),
        ("path", 3),
    ],
)
def test_manifest_rejects_wrong_field_types_and_invalid_digests(field, value):
    payload = import_module("vntyper.scripts.calibration_payload")
    raw = manifest_rows()
    raw[0][field] = value
    with pytest.raises(ValueError, match=field):
        payload.decode_payload_manifest(raw)


@pytest.mark.parametrize(
    "rows",
    [[{"path": "a", "size_bytes": 1}], [None], [{"path": "a", "size_bytes": 1, "sha256": "a" * 64, "extra": True}]],
)
def test_manifest_rejects_missing_or_unknown_fields(rows):
    payload = import_module("vntyper.scripts.calibration_payload")
    with pytest.raises(ValueError, match="fields"):
        payload.decode_payload_manifest(rows)


@pytest.mark.parametrize("paths", [["b", "a"], ["a", "a"], ["a", "a/b"], ["a/b", "a/b/c"]])
def test_manifest_layout_has_unique_sorted_files_and_no_file_parent(paths):
    payload = import_module("vntyper.scripts.calibration_payload")
    rows = [{"path": path, "size_bytes": 0, "sha256": "a" * 64} for path in paths]
    with pytest.raises(ValueError, match="layout"):
        payload.decode_payload_manifest(rows)


def test_nested_payloads_and_empty_files_are_valid():
    payload = import_module("vntyper.scripts.calibration_payload")
    rows = [{"path": "assets/background.json", "size_bytes": 0, "sha256": "a" * 64}]
    assert payload.payload_manifest_document(payload.decode_payload_manifest(rows)) == rows


def test_observations_require_exact_files_sizes_and_digests():
    payload = import_module("vntyper.scripts.calibration_payload")
    manifest = payload.decode_payload_manifest(manifest_rows())
    observed = {"background.json": (11, "a" * 64), "decision-profile.json": (23, "b" * 64)}
    assert payload.validate_payload_observations(manifest, observed) is None

    changed = deepcopy(observed)
    changed["background.json"] = (11, "c" * 64)
    with pytest.raises(ValueError, match="differs"):
        payload.validate_payload_observations(manifest, changed)
    changed["background.json"] = (12, "a" * 64)
    with pytest.raises(ValueError, match="differs"):
        payload.validate_payload_observations(manifest, changed)
    with pytest.raises(ValueError, match="files"):
        payload.validate_payload_observations(manifest, {"background.json": (11, "a" * 64)})
    with pytest.raises(ValueError, match="files"):
        payload.validate_payload_observations(manifest, {**observed, "extra": (0, "d" * 64)})


@pytest.mark.parametrize(
    "observed",
    [
        None,
        [],
        {"background.json": (True, "a" * 64), "decision-profile.json": (23, "b" * 64)},
        {"background.json": [11, "a" * 64], "decision-profile.json": (23, "b" * 64)},
    ],
)
def test_observations_do_not_coerce_or_accept_boolean_size(observed):
    payload = import_module("vntyper.scripts.calibration_payload")
    manifest = payload.decode_payload_manifest(manifest_rows())
    with pytest.raises(ValueError):
        payload.validate_payload_observations(manifest, observed)


def test_public_projection_and_validation_refuse_undecoded_manifest():
    payload = import_module("vntyper.scripts.calibration_payload")
    with pytest.raises(ValueError, match="PayloadManifest"):
        payload.payload_manifest_document(manifest_rows())
    with pytest.raises(ValueError, match="PayloadManifest"):
        payload.validate_payload_observations(manifest_rows(), {})
