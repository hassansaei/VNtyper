"""Training-only adVNTR background staging and provenance tests."""

from dataclasses import replace
from pathlib import Path
from types import MappingProxyType
from unittest.mock import patch

import pytest

from tests.unit.test_calibration_caller_controller import _caller_study
from tests.unit.test_calibration_portable_background import background
from vntyper.modules.advntr.advntr_background import BackgroundFitResult
from vntyper.modules.advntr.advntr_calibration_policy import AdvntrCapabilities, AdvntrToolPin
from vntyper.scripts.calibration_caller_observations import CallerTruth, CallerTruthRow
from vntyper.scripts.calibration_caller_roster import EligibleCallerMember
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

pytestmark = pytest.mark.unit


def test_staging_includes_known_controls_and_cases_but_explicitly_excludes_unknown(tmp_path: Path) -> None:
    module = __import__("vntyper.scripts.calibration_caller_background_training", fromlist=["x"])
    members = tuple(
        EligibleCallerMember(key, f"group-{index}", ("all",))
        for index, key in enumerate(("negative-key", "positive-key", "unknown-key"))
    )
    source = type(
        "Source",
        (),
        {"keys": tuple(member.key for member in members), "roster": type("R", (), {"members": members})(),
         "truth_asset": object()},
    )()
    rows = (
        CallerTruthRow("negative-key", "negative", ()),
        CallerTruthRow("positive-key", "positive", ("variant-a",)),
        CallerTruthRow("unknown-key", "unknown", None),
    )
    truth = CallerTruth(rows, MappingProxyType({row.key: row for row in rows}), "1" * 64)
    runs = tuple(
        type(
            "Run", (),
            {"manifest_key": member.key, "vntr_ids": (17,), "capture_policy_sha256": "c" * 64,
             "assets": {"advntr_capture": object(), "advntr_model": type("A", (), {"sha256": "d" * 64})()}},
        )()
        for member in members
    )
    record = {
        "producer": {"package_version": "2.4.1", "build_id": "b" * 64, "source_revision": "a" * 40},
        "assets": {"model_sha256": "d" * 64}, "capture_policy": {},
    }
    with (
        patch.object(module, "decode_caller_truth", return_value=truth),
        patch.object(module, "read_target_json", return_value={}),
        patch.object(module, "select_target_runs", return_value=runs),
        patch.object(module, "read_target_asset", side_effect=[b"negative", b"positive"]) as read_capture,
        patch.object(module, "_capture_records", return_value=(record,)),
        patch.object(module, "advntr_canonical_sha256", return_value="c" * 64),
    ):
        labels_path, diagnostic_path, negative, positive, unknown, pin = module._stage_training(
            tmp_path, _caller_study(), source, object(), type("P", (), {"baseline_policy_sha256": "2" * 64})()
        )

    labels = load_strict_json_object(labels_path.read_bytes())["samples"]
    assert [row["truth"] for row in labels] == [False, True]
    assert all(row["array_length"] is None for row in labels)
    assert "negative-key" not in labels_path.read_text()
    assert negative == ("negative-key",)
    assert positive == ("positive-key",)
    assert unknown == ("unknown-key",)
    assert read_capture.call_count == 2
    assert pin.source_revision == "a" * 40
    assert load_strict_json_object(diagnostic_path.read_bytes())["minimum_read_support"] == 3


def test_native_fit_is_portable_projected_and_bound_to_training_receipt(tmp_path: Path) -> None:
    module = __import__("vntyper.scripts.calibration_caller_background_training", fromlist=["x"])
    study = _caller_study()
    protocol = study.protocol
    output = (tmp_path / "background fit").resolve()
    runs = type("Runs", (), {"sha256": "6" * 64})()
    source = type("Source", (), {"sha256": "7" * 64})()
    receipt = type("Receipt", (), {"sha256": "8" * 64})()

    def native_fit(plan, _pin):
        plan.output_directory.mkdir()
        raw = canonical_json_bytes(background())
        (plan.output_directory / "caller-training.background.json").write_bytes(raw)
        capabilities = AdvntrCapabilities("2.4.1", "b" * 64, "a" * 40, (), (), (), (), "c" * 64)
        return BackgroundFitResult(capabilities, "9" * 64, "a" * 64, MappingProxyType({"native": "b" * 64}))

    with (
        patch.object(module, "_require_context", return_value=protocol),
        patch.object(
            module,
            "_stage_training",
            return_value=(
                output / "capture/labels.json", output / "capture/policy.json",
                ("negative",), ("positive",), ("unknown",), AdvntrToolPin("2.4.1", "b" * 64, "a" * 40),
            ),
        ),
        patch.object(module, "fit_background", side_effect=native_fit),
    ):
        result = module.fit_caller_training_background(
            study, runs, source, receipt, argv_prefix=("/path with spaces/python",), output=output
        )

    assert b"invented source label" not in result.background_bytes
    assert result.negative_keys == ("negative",)
    assert result.diagnostic_positive_keys == ("positive",)
    assert result.excluded_unknown_keys == ("unknown",)
    assert module.training_background_document(result)["training_evidence_sha256"] == result.training_evidence_sha256
    with pytest.raises(ValueError, match="bytes"):
        module.training_background_document(replace(result, background_bytes=b"changed"))
    assert (output / "training-background.json").is_file()


def test_failed_native_fit_cleans_private_staging(tmp_path: Path) -> None:
    module = __import__("vntyper.scripts.calibration_caller_background_training", fromlist=["x"])
    study = _caller_study()
    output = (tmp_path / "failed").resolve()
    with (
        patch.object(module, "_require_context", return_value=study.protocol),
        patch.object(module, "_stage_training", side_effect=RuntimeError("failed fit")),
        pytest.raises(RuntimeError, match="failed fit"),
    ):
        module.fit_caller_training_background(
            study, type("Runs", (), {})(), type("Source", (), {})(), type("Receipt", (), {})(),
            argv_prefix=("advntr",), output=output,
        )
    assert not output.exists()
