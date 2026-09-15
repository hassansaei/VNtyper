"""Public entry points for standard length estimation and development calibration."""

import logging
from pathlib import Path

import pytest

from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("flags", "enabled"),
    [([], None), (["--estimate-vntr-length"], True), (["--no-estimate-vntr-length"], False)],
)
def test_standard_length_toggle_has_configurable_default(flags: list[str], enabled: bool | None) -> None:
    args = build_parser().parse_args(["pipeline", "--bam", "input.bam", *flags])
    assert args.estimate_vntr_length is enabled
    assert args.standard_length_model is None


def test_locally_fitted_research_model_is_explicitly_separate_from_approved_bundle() -> None:
    args = build_parser().parse_args(["pipeline", "--bam", "input.bam", "--standard-length-model", "model.json"])
    assert args.standard_length_model == Path("model.json")
    assert args.length_model is None


@pytest.mark.parametrize("target", ["auto", "length", "callers", "both"])
def test_cohort_command_accepts_optional_length_truth_and_reference(target: str) -> None:
    args = build_parser().parse_args(
        ["calibrate", "cohort", "--manifest", "samples.tsv", "--target", target, "--output", "out"]
    )
    assert args.manifest == Path("samples.tsv")
    assert args.reference is None
    assert args.target == target
    assert args.folds == 5
    assert args.seed == 20260915
    assert args.count_convention == "source-reported"


def test_cohort_command_preserves_explicit_reference_and_native_policy_inputs() -> None:
    args = build_parser().parse_args(
        [
            "calibrate",
            "cohort",
            "--manifest",
            "samples.tsv",
            "--output",
            "out",
            "--reference",
            "reference.fa",
            "--caller-runs",
            "runs",
            "--caller-policies",
            "policies.json",
            "--folds",
            "3",
            "--seed",
            "17",
            "--count-convention",
            "complete",
        ]
    )
    assert args.reference == Path("reference.fa")
    assert args.caller_runs == Path("runs")
    assert args.caller_policies == Path("policies.json")
    assert (args.folds, args.seed, args.count_convention) == (3, 17, "complete")


def test_standard_model_cannot_be_selected_as_pipeline_log(tmp_path: Path) -> None:
    from vntyper.scripts.cli_logging_safety import validate_pipeline_log_destination

    model = tmp_path / "model.json"
    model.write_text('{"model": "invented"}')
    args = build_parser().parse_args(
        ["pipeline", "--bam", str(tmp_path / "reads" / "input.bam"), "--standard-length-model", str(model)]
    )
    with pytest.raises(ValueError, match="operator"):
        validate_pipeline_log_destination(model, args, {})
    assert model.read_text() == '{"model": "invented"}'


@pytest.mark.parametrize("successful", [True, False])
def test_cohort_dispatch_installs_complete_scientific_outcome(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, successful: bool
) -> None:
    from vntyper.scripts import cli_calibrate

    output = tmp_path / "out"
    parser = build_parser()
    args = parser.parse_args(["calibrate", "cohort", "--manifest", "samples.tsv", "--output", str(output)])
    observed = []

    def run_cohort(observed_args, staging):
        observed.append(observed_args)
        (staging / "report.html").write_text("Complete scientific result")
        return successful

    monkeypatch.setattr("vntyper.scripts.calibration_cohort.run_cohort_calibration", run_cohort)
    if successful:
        cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
    else:
        with pytest.raises(SystemExit) as error:
            cli_calibrate.handle_calibrate(args, {}, parser, logging.INFO, None)
        assert error.value.code == 1
    assert observed == [args]
    assert (output / "report.html").read_text() == "Complete scientific result"
