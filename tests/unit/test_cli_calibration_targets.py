"""Target CLI dispatch preserves dominance and refuses mismatched objectives."""

import logging
from pathlib import Path

import pytest

from vntyper.scripts import cli_calibrate
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


def test_exact_caller_training_executable_is_an_explicit_path():
    arguments = [
        "calibrate",
        "fit",
        "--target",
        "callers",
        "--evidence",
        "evidence",
        "--objective",
        "caller-safety-v1",
        "--output",
        "output",
    ]
    assert build_parser().parse_args(arguments).advntr_executable is None
    assert build_parser().parse_args(
        arguments + ["--advntr-executable", "/opt/native/advntr"]
    ).advntr_executable == Path("/opt/native/advntr")


@pytest.mark.parametrize("target,objective", [("length", "length-total-v1"), ("callers", "caller-safety-v1")])
def test_target_fit_dispatch_is_explicit_and_atomic(target, objective, tmp_path, monkeypatch):
    parser = build_parser()
    output = tmp_path / "result"
    args = parser.parse_args(
        [
            "calibrate",
            "fit",
            "--target",
            target,
            "--evidence",
            str(tmp_path / "evidence"),
            "--objective",
            objective,
            "--exposure-ledger",
            str(tmp_path / "ledger"),
            "--output",
            str(output),
        ]
    )

    def produce(observed, staging):
        assert observed.target == target
        assert observed.objective == objective
        assert staging != output
        (staging / "result.json").write_text("{}\n")
        return True

    monkeypatch.setitem(cli_calibrate.TARGET_OPERATIONS, (target, "fit"), produce)
    monkeypatch.setitem(cli_calibrate.OPERATIONS, "fit", lambda *_: pytest.fail("dominance fitter invoked"))
    cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert (output / "result.json").read_text() == "{}\n"


@pytest.mark.parametrize("target,objective", [("length", "caller-safety-v1"), ("dominance", "length-total-v1")])
def test_wrong_objective_is_usage_error_before_output(target, objective, tmp_path, monkeypatch):
    parser = build_parser()
    args = parser.parse_args(
        [
            "calibrate",
            "fit",
            "--target",
            target,
            "--evidence",
            "evidence",
            "--objective",
            objective,
            "--output",
            str(tmp_path / "result"),
        ]
    )
    monkeypatch.setattr(cli_calibrate, "_atomic_output", lambda *_: pytest.fail("output opened"))
    with pytest.raises(SystemExit) as error:
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert error.value.code == 2


def test_v2_requires_external_exposure_ledger_before_output(tmp_path, monkeypatch):
    parser = build_parser()
    args = parser.parse_args(
        [
            "calibrate",
            "fit",
            "--target",
            "length",
            "--evidence",
            "evidence",
            "--objective",
            "length-total-v1",
            "--output",
            str(tmp_path / "result"),
        ]
    )
    monkeypatch.setattr(cli_calibrate, "_atomic_output", lambda *_: pytest.fail("output opened"))
    with pytest.raises(SystemExit) as error:
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert error.value.code == 2


def test_assess_parser_requires_fixed_profile_and_declared_inputs():
    args = build_parser().parse_args(
        [
            "calibrate",
            "assess",
            "--target",
            "length",
            "--profile",
            "profile",
            "--intake",
            "intake",
            "--runs",
            "runs.json",
            "--exposure-ledger",
            "ledger",
            "--output",
            "out",
        ]
    )
    assert args.calibration_operation == "assess"
    assert args.profile == Path("profile")
    assert args.intake == Path("intake")
    assert args.runs == Path("runs.json")
    assert not hasattr(args, "objective")


def test_pipeline_caller_bundle_and_context_paths_are_explicit():
    args = build_parser().parse_args(
        ["pipeline", "--bam", "reads.bam", "--calibration-bundle", "approved", "--calibration-context", "context.json"]
    )
    assert args.calibration_bundle == Path("approved")
    assert args.calibration_context == Path("context.json")


@pytest.mark.parametrize("target", ["length", "callers"])
def test_export_dispatch_requires_all_approval_artifacts_but_no_new_exposure(target, tmp_path, monkeypatch):
    parser = build_parser()
    args = parser.parse_args(
        [
            "calibrate",
            "export",
            "--target",
            target,
            "--profile",
            "profile",
            "--validation",
            "validation.json",
            "--evaluation",
            "locked.json",
            "--authority",
            "authority.json",
            "--completion",
            "completion.json",
            "--output",
            str(tmp_path / "portable"),
        ]
    )

    def produce(observed, staging):
        assert observed.validation == Path("validation.json")
        assert observed.completion == Path("completion.json")
        (staging / "approved.json").write_text("{}")
        return True

    monkeypatch.setitem(cli_calibrate.TARGET_OPERATIONS, (target, "export"), produce)
    cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert (args.output / "approved.json").read_text() == "{}"


def test_locked_target_evaluation_requires_prior_validation_and_authority_before_output(tmp_path, monkeypatch):
    parser = build_parser()
    args = parser.parse_args(
        [
            "calibrate",
            "evaluate",
            "--target",
            "length",
            "--profile",
            "profile",
            "--evidence",
            "evidence",
            "--custody",
            "custody",
            "--exposure-ledger",
            "ledger",
            "--output",
            str(tmp_path / "output"),
        ]
    )
    monkeypatch.setattr(cli_calibrate, "_atomic_output", lambda *_: pytest.fail("output opened"))
    with pytest.raises(SystemExit) as error:
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert error.value.code == 2


@pytest.mark.parametrize("target", ["length", "callers"])
def test_target_extract_uses_sealed_sources_without_legacy_truth_arguments(target, tmp_path, monkeypatch):
    parser = build_parser()
    argv = [
        "calibrate",
        "extract",
        "--target",
        target,
        "--study",
        "study.json",
        "--sources",
        "sealed",
        "--runs",
        "runs.json",
        "--output",
        str(tmp_path / "evidence"),
    ]
    if target == "length":
        argv += ["--length-annotation", "annotation.json"]
    args = parser.parse_args(argv)

    def extract(observed, staging):
        assert observed.study == Path("study.json")
        assert observed.sources == Path("sealed")
        assert observed.truth is None
        (staging / "study.json").write_text("{}")
        return True

    monkeypatch.setitem(cli_calibrate.TARGET_OPERATIONS, (target, "extract"), extract)
    cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert (args.output / "study.json").exists()


def test_dominance_extract_still_requires_truth_and_partitions_before_output(tmp_path, monkeypatch):
    parser = build_parser()
    args = parser.parse_args(["calibrate", "extract", "--runs", "runs.json", "--output", str(tmp_path / "output")])
    monkeypatch.setattr(cli_calibrate, "_atomic_output", lambda *_: pytest.fail("output opened"))
    with pytest.raises(SystemExit) as error:
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert error.value.code == 2


def test_real_length_extract_fit_and_fixed_validation_through_cli(tmp_path):
    from tests.unit.test_calibration_length_controller import fit_fixture
    from vntyper.scripts.calibration_artifact_io import load_object, write_json
    from vntyper.scripts.calibration_target_custody import initialize_target_custody

    fixture, study = fit_fixture(tmp_path)
    parser = build_parser()
    evidence = tmp_path / "extracted"
    profile = tmp_path / "candidate"

    def run(tokens):
        args = parser.parse_args(["calibrate", *tokens])
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)

    run(
        [
            "extract",
            "--target",
            "length",
            "--study",
            str(fixture.evidence / "study.json"),
            "--sources",
            str(fixture.evidence),
            "--runs",
            str(fixture.evidence / "runs.json"),
            "--length-annotation",
            str(fixture.evidence / "annotation.json"),
            "--output",
            str(evidence),
        ]
    )
    run(
        [
            "fit",
            "--target",
            "length",
            "--objective",
            "length-total-v1",
            "--evidence",
            str(evidence),
            "--exposure-ledger",
            str(fixture.exposure_ledger),
            "--output",
            str(profile),
        ]
    )
    model = load_object(profile / "payload" / "length-model.json", "test model")
    assert model["coefficients"] == pytest.approx([50])
    assert model["intercept"] == pytest.approx(10)
    confirmation = tmp_path / "confirmation"
    confirmation.mkdir()
    write_json(confirmation / "runs.json", load_object(evidence / "runs.json", "test runs"))
    write_json(
        confirmation / "source.json", load_object(evidence / "roles" / "validation" / "source.json", "test source")
    )
    custody = tmp_path / "custody"
    initialize_target_custody(custody, exposure_ledger_id=study.exposure_ledger_id)
    output = tmp_path / "validation-result"
    run(
        [
            "validate",
            "--target",
            "length",
            "--profile",
            str(profile),
            "--evidence",
            str(confirmation),
            "--exposure-ledger",
            str(fixture.exposure_ledger),
            "--custody",
            str(custody),
            "--output",
            str(output),
        ]
    )
    assert load_object(output / "validation-attestation.json", "test outcome")["status"] == "passed"
    assert load_object(profile / "candidate.json", "test candidate")["status"] == "research-only"
    assert (output / "report.html").is_file()


@pytest.mark.parametrize(
    "name,module_name,function_name",
    [
        ("_fit_length", "calibration_length_controller", "fit_length_bundle"),
        ("_assess_length", "calibration_length_controller", "assess_length_bundle"),
        ("_fit_callers", "calibration_caller_controller", "fit_caller_bundle"),
        ("_assess_callers", "calibration_caller_controller", "assess_caller_bundle"),
        ("_export_target", "calibration_export", "export_calibration_bundle"),
    ],
)
def test_target_adapters_preserve_failed_scientific_outcomes(name, module_name, function_name, tmp_path, monkeypatch):
    from argparse import Namespace
    from importlib import import_module

    provider = import_module("vntyper.scripts." + module_name)
    args = Namespace(target="length")

    def operation(observed, output):
        assert observed is args and output == tmp_path
        return False

    monkeypatch.setattr(provider, function_name, operation)
    assert getattr(cli_calibrate, name)(args, tmp_path) is False


@pytest.mark.parametrize("name,role", [("_validate_target", "validation"), ("_evaluate_target", "locked-heldout")])
def test_confirmation_adapters_fix_the_role(name, role, tmp_path, monkeypatch):
    from argparse import Namespace

    from vntyper.scripts import calibration_confirmation_controller

    args = Namespace(target="length")

    def operation(observed, output, *, role):
        assert observed is args and output == tmp_path
        assert role == expected_role
        return False

    expected_role = role
    monkeypatch.setattr(calibration_confirmation_controller, "confirm_calibration_bundle", operation)
    assert getattr(cli_calibrate, name)(args, tmp_path) is False


@pytest.mark.parametrize(
    "extra",
    [
        ["--target", "length", "--sources", "sources"],
        ["--target", "length", "--sources", "sources", "--study", "study"],
        ["--target", "callers", "--sources", "sources", "--study", "study", "--length-annotation", "annotation"],
        [
            "--target",
            "length",
            "--sources",
            "sources",
            "--study",
            "study",
            "--truth",
            "truth",
            "--length-annotation",
            "annotation",
        ],
        ["--truth", "truth", "--partitions", "partitions", "--study", "study"],
    ],
)
def test_extract_incompatible_target_options_fail_before_output(extra, tmp_path, monkeypatch):
    parser = build_parser()
    args = parser.parse_args(["calibrate", "extract", "--runs", "runs", "--output", str(tmp_path / "output"), *extra])
    monkeypatch.setattr(cli_calibrate, "_atomic_output", lambda *_: pytest.fail("output opened"))
    with pytest.raises(SystemExit) as error:
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert error.value.code == 2


def test_validation_cannot_silently_ignore_locked_authority(tmp_path, monkeypatch):
    parser = build_parser()
    args = parser.parse_args(
        [
            "calibrate",
            "validate",
            "--target",
            "length",
            "--profile",
            "profile",
            "--evidence",
            "evidence",
            "--authority",
            "authority",
            "--custody",
            "custody",
            "--exposure-ledger",
            "ledger",
            "--output",
            str(tmp_path / "out"),
        ]
    )
    monkeypatch.setattr(cli_calibrate, "_atomic_output", lambda *_: pytest.fail("output opened"))
    with pytest.raises(SystemExit) as error:
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    assert error.value.code == 2
