"""A payload digest binds every file and a unique, safe layout."""

from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
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


def test_manifest_public_boundaries_recompute_the_canonical_digest():
    payload = import_module("vntyper.scripts.calibration_payload")
    manifest = payload.decode_payload_manifest(manifest_rows())
    forged = replace(manifest, sha256="f" * 64)

    with pytest.raises(ValueError, match="canonical"):
        payload.payload_manifest_document(forged)
    with pytest.raises(ValueError, match="canonical"):
        payload.validate_payload_observations(forged, {})


def caller_bundle_document(callers=None, *, advntr_digest="b" * 64, background_digest=None):
    if callers is None:
        callers = ["advntr", "kestrel"]
    return {
        "schema_version": "caller-bundle-v2",
        "required_callers": callers,
        "components": {
            "decision-profile.json": "a" * 64,
            "advntr-policy.json": advntr_digest,
            "background.json": background_digest,
        },
    }


def test_caller_bundle_descriptor_is_closed_hash_bound_and_immutable():
    payload = import_module("vntyper.scripts.calibration_payload")
    raw = caller_bundle_document(background_digest="c" * 64)
    expected = deepcopy(raw)
    descriptor = payload.decode_caller_bundle_descriptor(raw)

    assert descriptor.sha256 == canonical_sha256(expected)
    assert descriptor.required_callers == ("advntr", "kestrel")
    raw["components"]["decision-profile.json"] = "d" * 64
    assert payload.caller_bundle_descriptor_document(descriptor) == expected
    with pytest.raises(FrozenInstanceError):
        descriptor.background_sha256 = None


def test_caller_descriptor_entry_does_not_change_the_r4_manifest_list_hash():
    payload = import_module("vntyper.scripts.calibration_payload")
    descriptor = payload.decode_caller_bundle_descriptor(caller_bundle_document())
    rows = [
        {"path": "caller-bundle.json", "size_bytes": 301, "sha256": descriptor.sha256},
        {"path": "decision-profile.json", "size_bytes": 23, "sha256": "a" * 64},
    ]

    assert payload.decode_payload_manifest(rows).sha256 == canonical_sha256(rows)


@pytest.mark.parametrize(
    "change",
    [
        "root-extra",
        "root-missing",
        "component-extra",
        "component-missing",
        "schema",
        "unsorted-callers",
        "missing-kestrel",
        "unknown-caller",
        "missing-advntr-policy",
        "unexpected-advntr-policy",
        "bad-component-digest",
    ],
)
def test_caller_bundle_descriptor_rejects_open_or_inconsistent_content(change):
    payload = import_module("vntyper.scripts.calibration_payload")
    raw = caller_bundle_document()
    if change == "root-extra":
        raw["extra"] = True
    elif change == "root-missing":
        del raw["components"]
    elif change == "component-extra":
        raw["components"]["extra.json"] = None
    elif change == "component-missing":
        del raw["components"]["background.json"]
    elif change == "schema":
        raw["schema_version"] = "caller-bundle-v1"
    elif change == "unsorted-callers":
        raw["required_callers"] = ["kestrel", "advntr"]
    elif change == "missing-kestrel":
        raw["required_callers"] = ["advntr"]
    elif change == "unknown-caller":
        raw["required_callers"] = ["kestrel", "other"]
    elif change == "missing-advntr-policy":
        raw["components"]["advntr-policy.json"] = None
    elif change == "unexpected-advntr-policy":
        raw = caller_bundle_document(["kestrel"], advntr_digest="b" * 64)
    else:
        raw["components"]["decision-profile.json"] = "bad"

    with pytest.raises(ValueError):
        payload.decode_caller_bundle_descriptor(raw)


def test_kestrel_only_bundle_omits_advntr_policy_and_descriptor_boundary_revalidates():
    payload = import_module("vntyper.scripts.calibration_payload")
    raw = caller_bundle_document(["kestrel"], advntr_digest=None)
    descriptor = payload.decode_caller_bundle_descriptor(raw)
    assert payload.caller_bundle_descriptor_document(descriptor) == raw

    forged = replace(descriptor, sha256="f" * 64)
    with pytest.raises(ValueError, match="canonical"):
        payload.caller_bundle_descriptor_document(forged)
