"""Training-only adVNTR background staging and provenance tests."""

from dataclasses import replace
from pathlib import Path
from types import MappingProxyType, SimpleNamespace
from unittest.mock import patch

import pytest

from tests.unit.test_calibration_caller_controller import _caller_study
from tests.unit.test_calibration_caller_protocol import changed_policy, protocol_document
from tests.unit.test_calibration_portable_background import background
from vntyper.modules.advntr.advntr_background import BackgroundFitResult
from vntyper.modules.advntr.advntr_calibration_policy import AdvntrCapabilities, AdvntrToolPin
from vntyper.scripts.calibration_caller_observations import CallerTruth, CallerTruthRow
from vntyper.scripts.calibration_caller_protocol import decode_caller_protocol
from vntyper.scripts.calibration_caller_roster import EligibleCallerMember
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

pytestmark = pytest.mark.unit


def _context(protocol=None):
    study = _caller_study()
    study = SimpleNamespace(
        target=study.target,
        baseline=study.baseline,
        protocol=study.protocol if protocol is None else protocol,
        sha256=study.sha256,
        partitions=study.partitions,
        exposure_ledger_id=study.exposure_ledger_id,
    )
    runs = SimpleNamespace(sha256="6" * 64)
    identities = (("physical-readset", "1" * 64), ("specimen", "2" * 64))
    source = SimpleNamespace(
        role="training",
        study_sha256=study.sha256,
        run_manifest_sha256=runs.sha256,
        identities=identities,
    )
    receipt = SimpleNamespace(
        target="callers",
        role="training",
        study_sha256=study.sha256,
        partition_sha256=study.partitions.sha256,
        evidence_sha256=source.sha256 if hasattr(source, "sha256") else "7" * 64,
        membership_sha256=canonical_sha256([{"namespace": name, "sha256": digest} for name, digest in identities]),
        exposure_ledger_id=study.exposure_ledger_id,
    )
    source.sha256 = receipt.evidence_sha256
    return study, runs, source, receipt


def _receipt(module):
    portable = canonical_json_bytes(background())
    native = BackgroundFitResult(
        AdvntrCapabilities("2.4.1", "b" * 64, "a" * 40, (), (), (), (), "c" * 64),
        "9" * 64,
        "a" * 64,
        MappingProxyType({"native.json": "b" * 64}),
    )
    result = module._result(
        SimpleNamespace(sha256="1" * 64),
        SimpleNamespace(sha256="2" * 64),
        SimpleNamespace(sha256="3" * 64),
        SimpleNamespace(sha256="4" * 64),
        portable,
        native,
        ("negative",),
        ("positive",),
        ("unknown",),
    )
    return portable, result, module.training_background_document(result)


@pytest.mark.parametrize("change", ["source", "receipt", "membership"])
def test_context_refuses_forged_training_source_and_exposure_bindings(change: str) -> None:
    module = __import__("vntyper.scripts.calibration_caller_background_training", fromlist=["x"])
    study, runs, source, receipt = _context()
    if change == "source":
        source.role = "policy-selection"
    elif change == "receipt":
        receipt.evidence_sha256 = "8" * 64
    else:
        receipt.membership_sha256 = "8" * 64
    with (
        patch.object(module, "target_study_document"),
        patch.object(module, "target_runs_document"),
        patch.object(module, "role_source_document"),
        patch.object(module, "exposure_receipt_document"),
        pytest.raises(ValueError, match="training source|receipt membership"),
    ):
        module._require_context(study, runs, source, receipt)


def test_context_accepts_exact_baseline_when_all_candidates_are_legacy() -> None:
    module = __import__("vntyper.scripts.calibration_caller_background_training", fromlist=["x"])
    baseline = _caller_study().protocol.baseline_policy
    legacy = changed_policy(**{"/components/advntr/calibrated_calling/mode": "legacy"})
    protocol = decode_caller_protocol(protocol_document(baseline, [(legacy, 1)]), baseline_policy=baseline)
    study, runs, source, receipt = _context(protocol)
    with (
        patch.object(module, "target_study_document"),
        patch.object(module, "target_runs_document"),
        patch.object(module, "role_source_document"),
        patch.object(module, "exposure_receipt_document"),
    ):
        assert module._require_context(study, runs, source, receipt) == protocol


@pytest.mark.parametrize("change", ["overlap", "empty-artifacts", "artifact-hash", "training-hash", "receipt-hash"])
def test_training_background_decoder_refuses_forged_membership_and_hashes(change: str) -> None:
    module = __import__("vntyper.scripts.calibration_caller_background_training", fromlist=["x"])
    portable, _result, document = _receipt(module)
    if change == "overlap":
        document["diagnostic_positive_keys"] = ["negative"]
    elif change == "empty-artifacts":
        document["fitter_artifact_sha256"] = {}
    elif change == "artifact-hash":
        document["fitter_artifact_sha256"] = {"native.json": "invalid"}
    elif change == "training-hash":
        document["training_evidence_sha256"] = "0" * 64
    else:
        document["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        module.decode_training_background_document(document, portable)


def test_staging_includes_known_controls_and_cases_but_explicitly_excludes_unknown(tmp_path: Path) -> None:
    module = __import__("vntyper.scripts.calibration_caller_background_training", fromlist=["x"])
    members = tuple(
        EligibleCallerMember(key, f"group-{index}", ("all",))
        for index, key in enumerate(("negative-key", "positive-key", "unknown-key"))
    )
    source = type(
        "Source",
        (),
        {
            "keys": tuple(member.key for member in members),
            "roster": type("R", (), {"members": members})(),
            "truth_asset": object(),
        },
    )()
    rows = (
        CallerTruthRow("negative-key", "negative", ()),
        CallerTruthRow("positive-key", "positive", ("variant-a",)),
        CallerTruthRow("unknown-key", "unknown", None),
    )
    truth = CallerTruth(rows, MappingProxyType({row.key: row for row in rows}), "1" * 64)
    runs = tuple(
        type(
            "Run",
            (),
            {
                "manifest_key": member.key,
                "vntr_ids": (17,),
                "capture_policy_sha256": "c" * 64,
                "assets": {"advntr_capture": object(), "advntr_model": type("A", (), {"sha256": "d" * 64})()},
            },
        )()
        for member in members
    )
    record = {
        "producer": {"package_version": "2.4.1", "build_id": "b" * 64, "source_revision": "a" * 40},
        "assets": {"model_sha256": "d" * 64},
        "capture_policy": {},
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
                output / "capture/labels.json",
                output / "capture/policy.json",
                ("negative",),
                ("positive",),
                ("unknown",),
                AdvntrToolPin("2.4.1", "b" * 64, "a" * 40),
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
            study,
            type("Runs", (), {})(),
            type("Source", (), {})(),
            type("Receipt", (), {})(),
            argv_prefix=("advntr",),
            output=output,
        )
    assert not output.exists()
