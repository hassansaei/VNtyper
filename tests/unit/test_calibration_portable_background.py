"""Portable background projection redacts only narrative provenance before binding."""

from copy import deepcopy
from importlib import import_module

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

pytestmark = pytest.mark.unit


def module():
    return import_module("vntyper.scripts.calibration_portable_background")


def background():
    return {
        "schema": "advntr.frameshift.background",
        "version": 1,
        "provenance": "invented source label; 12 training groups; development notes",
        "default_probability": 0.03125,
        "states": {"D1_1": 0.015625, "I2_1_T_LEN1&D3_1": 0.0625},
    }


def test_projection_preserves_every_model_parameter_and_removes_free_text_before_hashing():
    raw = background()
    original = deepcopy(raw)
    projected = module().project_portable_background(raw)
    assert raw == original
    assert projected == {**raw, "provenance": "VNtyper calibrated background; portable runtime parameters only"}
    assert projected["states"] is not raw["states"]
    assert canonical_sha256(projected) != canonical_sha256(raw)
    assert b"invented source" not in canonical_json_bytes(projected)
    assert module().validate_portable_background(load_strict_json_object(canonical_json_bytes(projected))) == projected
    assert module().project_portable_background(projected) == projected
    projected["states"]["D1_1"] = 0.5
    assert raw == original


def test_unprojected_free_provenance_is_rejected_at_export_boundary():
    with pytest.raises(ValueError, match="provenance"):
        module().validate_portable_background(background())


@pytest.mark.parametrize(
    "change",
    [
        "extra",
        "missing",
        "schema",
        "version",
        "provenance",
        "states",
        "state-name",
        "state-probability",
        "default-probability",
    ],
)
def test_malformed_source_is_not_silently_repaired(change):
    raw = background()
    if change == "extra":
        raw["training_metrics"] = {}
    elif change == "missing":
        del raw["states"]
    elif change == "schema":
        raw["schema"] = "other"
    elif change == "version":
        raw["version"] = True
    elif change == "provenance":
        raw["provenance"] = " "
    elif change == "states":
        raw["states"] = []
    elif change == "state-name":
        raw["states"] = {" D1_1": 0.1}
    elif change == "state-probability":
        raw["states"] = {"D1_1": True}
    else:
        raw["default_probability"] = float("nan")
    with pytest.raises(ValueError):
        module().project_portable_background(raw)


@pytest.mark.parametrize("probability", [0, 1, -0.1, float("inf"), "0.1"])
def test_background_probabilities_are_strictly_between_zero_and_one(probability):
    raw = background()
    raw["default_probability"] = probability
    with pytest.raises(ValueError):
        module().project_portable_background(raw)
