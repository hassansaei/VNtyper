from __future__ import annotations

from dataclasses import replace

import pytest

from vntyper.scripts.calibration_advntr_runtime_policy import (
    advntr_runtime_policy_document,
    build_advntr_runtime_policy,
    decode_advntr_runtime_policy,
    validate_advntr_runtime_policy,
)
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

pytestmark = pytest.mark.unit


def _caller(*, mode: str = "exact"):
    return decode_caller_policy_values(
        {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": ["advntr", "kestrel"],
            "values": {
                "/components/advntr/calibrated_calling/adapter_filter": False,
                "/components/advntr/calibrated_calling/cutoff": 0.001,
                "/components/advntr/calibrated_calling/minimum_read_match_ratio": 0.9,
                "/components/advntr/calibrated_calling/minimum_read_support": 3,
                "/components/advntr/calibrated_calling/mode": mode,
                "/components/advntr/calibrated_calling/prune_reverse": False,
                "/components/advntr/calibrated_calling/rare_unit_fraction": None,
                "/components/kestrel/alt_filtering/gg_depth_score_threshold": 0.00469,
                "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": 20,
                "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": 100,
                "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": 21,
                "/components/kestrel/confidence_assignment/depth_score_thresholds/high": 0.00515,
                "/components/kestrel/confidence_assignment/depth_score_thresholds/low": 0.00469,
                "/components/kestrel/confidence_assignment/reporting_floor": 0.003,
                "/components/kestrel/confidence_assignment/var_active_region_threshold": 200,
            },
        }
    )


def test_runtime_policy_binds_selected_values_and_native_assets() -> None:
    caller = _caller()
    runtime = build_advntr_runtime_policy(
        caller,
        model_sha256="a" * 64,
        background_sha256="b" * 64,
        capture_policy_sha256="c" * 64,
        advntr_revision="d" * 40,
    )

    document = advntr_runtime_policy_document(runtime)

    assert document == {
        "schema_version": "vntyper-advntr-calibrated-calling-v1",
        "mode": "exact",
        "cutoff": caller.values["/components/advntr/calibrated_calling/cutoff"],
        "minimum_read_support": caller.values["/components/advntr/calibrated_calling/minimum_read_support"],
        "rare_unit_fraction": caller.values["/components/advntr/calibrated_calling/rare_unit_fraction"],
        "adapter_filter": caller.values["/components/advntr/calibrated_calling/adapter_filter"],
        "minimum_read_match_ratio": caller.values["/components/advntr/calibrated_calling/minimum_read_match_ratio"],
        "prune_reverse": caller.values["/components/advntr/calibrated_calling/prune_reverse"],
        "background_sha256": "b" * 64,
        "model_sha256": "a" * 64,
        "capture_policy_sha256": "c" * 64,
        "advntr_revision": "d" * 40,
    }
    assert decode_advntr_runtime_policy(document) == runtime
    assert validate_advntr_runtime_policy(runtime, caller, background_raw_sha256="b" * 64) is runtime


def test_legacy_runtime_policy_forbids_background() -> None:
    caller = _caller(mode="legacy")
    runtime = build_advntr_runtime_policy(
        caller,
        model_sha256="a" * 64,
        background_sha256=None,
        capture_policy_sha256="c" * 64,
        advntr_revision="d" * 40,
    )

    assert validate_advntr_runtime_policy(runtime, caller, background_raw_sha256=None) is runtime
    with pytest.raises(ValueError, match="legacy.*background"):
        build_advntr_runtime_policy(
            caller,
            model_sha256="a" * 64,
            background_sha256="b" * 64,
            capture_policy_sha256="c" * 64,
            advntr_revision="d" * 40,
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("cutoff", True, "cutoff"),
        ("minimum_read_support", True, "minimum_read_support"),
        ("minimum_read_match_ratio", float("inf"), "minimum_read_match_ratio"),
        ("extra", 1, "fields"),
    ],
)
def test_runtime_policy_rejects_unknown_bool_and_nonfinite_fields(field: str, value: object, message: str) -> None:
    document = advntr_runtime_policy_document(
        build_advntr_runtime_policy(
            _caller(),
            model_sha256="a" * 64,
            background_sha256="b" * 64,
            capture_policy_sha256="c" * 64,
            advntr_revision="d" * 40,
        )
    )
    document[field] = value

    with pytest.raises(ValueError, match=message):
        decode_advntr_runtime_policy(document)


def test_runtime_policy_validation_rejects_forged_typed_content_and_raw_background() -> None:
    caller = _caller()
    runtime = build_advntr_runtime_policy(
        caller,
        model_sha256="a" * 64,
        background_sha256="b" * 64,
        capture_policy_sha256="c" * 64,
        advntr_revision="d" * 40,
    )

    with pytest.raises(ValueError, match="selected caller policy"):
        validate_advntr_runtime_policy(runtime, _caller(mode="legacy"), background_raw_sha256="b" * 64)
    with pytest.raises(ValueError, match="background bytes"):
        validate_advntr_runtime_policy(runtime, caller, background_raw_sha256="e" * 64)
    with pytest.raises(ValueError, match="canonical content"):
        advntr_runtime_policy_document(replace(runtime, model_sha256="f" * 64))


def test_integral_float_boundaries_survive_rfc8785_json_roundtrip() -> None:
    raw = advntr_runtime_policy_document(
        build_advntr_runtime_policy(
            _caller(),
            model_sha256="a" * 64,
            background_sha256="b" * 64,
            capture_policy_sha256="c" * 64,
            advntr_revision="d" * 40,
        )
    )
    raw["cutoff"] = 1.0
    raw["minimum_read_match_ratio"] = 1.0
    raw["rare_unit_fraction"] = 1.0

    decoded = decode_advntr_runtime_policy(load_strict_json_object(canonical_json_bytes(raw)))

    assert decoded.cutoff == 1.0
    assert decoded.minimum_read_match_ratio == 1.0
    assert decoded.rare_unit_fraction == 1.0
