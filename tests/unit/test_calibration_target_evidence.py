"""Target-v2 evidence extraction leaves sealed outcomes unopened."""

from argparse import Namespace
from importlib import import_module
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.unit.test_calibration_length_controller import fit_fixture

pytestmark = pytest.mark.unit


def test_length_target_extract_verifies_sources_without_opening_truth_and_is_fit_usable(tmp_path):
    module = import_module("vntyper.scripts.calibration_target_evidence")
    controller = import_module("vntyper.scripts.calibration_length_controller")
    args, study = fit_fixture(tmp_path)
    extracted = tmp_path / "extracted"
    extracted.mkdir()
    assert module.extract_target_evidence(
        args.evidence / "study.json",
        args.evidence / "runs.json",
        args.evidence,
        extracted,
        expected_target="length",
        length_annotation_path=args.evidence / "annotation.json",
    )
    assert {path.name for path in extracted.iterdir()} == {
        "annotation.json",
        "checksums.json",
        "roles",
        "runs.json",
        "study.json",
    }
    for role in ("training", "policy-selection", "validation", "locked-heldout"):
        assert {path.name for path in (extracted / "roles" / role).iterdir()} == {"source.json"}
    output = tmp_path / "fit"
    output.mkdir()
    fit_args = Namespace(evidence=extracted, exposure_ledger=args.exposure_ledger, objective="length-total-v1")
    assert controller.fit_length_bundle(fit_args, output)
    assert controller.load_length_research_profile(output).study == study


def test_target_extract_does_not_parse_or_hash_sealed_truth(tmp_path):
    module = import_module("vntyper.scripts.calibration_target_evidence")
    args, _ = fit_fixture(tmp_path)
    for path in (args.evidence / "roles").glob("*/truth.json"):
        path.write_bytes(b"not-json-and-not-the-committed-bytes")
    extracted = tmp_path / "extracted"
    extracted.mkdir()
    with patch(
        "vntyper.scripts.calibration_target_asset_io.read_target_asset", side_effect=AssertionError("opened outcome")
    ) as reader:
        assert module.extract_target_evidence(
            args.evidence / "study.json",
            args.evidence / "runs.json",
            args.evidence,
            extracted,
            expected_target="length",
            length_annotation_path=args.evidence / "annotation.json",
        )
    reader.assert_not_called()


def test_target_extract_accepts_relative_cli_paths_with_absolute_truth_commitments(tmp_path, monkeypatch):
    module = import_module("vntyper.scripts.calibration_target_evidence")
    args, _ = fit_fixture(tmp_path)
    output = tmp_path / "extracted"
    output.mkdir()
    monkeypatch.chdir(tmp_path)

    assert module.extract_target_evidence(
        Path("evidence/study.json"),
        Path("evidence/runs.json"),
        Path("evidence"),
        output,
        expected_target="length",
        length_annotation_path=Path("evidence/annotation.json"),
    )
    assert (output / "checksums.json").is_file()


def test_target_extract_relative_path_normalization_does_not_follow_role_symlinks(tmp_path, monkeypatch):
    module = import_module("vntyper.scripts.calibration_target_evidence")
    args, _ = fit_fixture(tmp_path)
    role = args.evidence / "roles" / "training"
    real_role = tmp_path / "training-source"
    role.rename(real_role)
    role.symlink_to(real_role, target_is_directory=True)
    output = tmp_path / "extracted"
    output.mkdir()
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValueError, match="role source inventory"):
        module.extract_target_evidence(
            Path("evidence/study.json"),
            Path("evidence/runs.json"),
            Path("evidence"),
            output,
            expected_target="length",
            length_annotation_path=Path("evidence/annotation.json"),
        )


@pytest.mark.parametrize(
    "expected_target,annotation,match",
    [
        ("callers", True, "target"),
        ("length", False, "annotation"),
        ("unsupported", True, "target"),
    ],
)
def test_target_extract_rejects_cli_target_drift_or_missing_length_annotation(
    tmp_path, expected_target, annotation, match
):
    module = import_module("vntyper.scripts.calibration_target_evidence")
    args, _ = fit_fixture(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    annotation_path = args.evidence / "annotation.json" if annotation else None
    with pytest.raises(ValueError, match=match):
        module.extract_target_evidence(
            args.evidence / "study.json",
            args.evidence / "runs.json",
            args.evidence,
            output,
            expected_target=expected_target,
            length_annotation_path=annotation_path,
        )
    assert not tuple(output.iterdir())
