"""Closed four-command calibration CLI and atomic output installation."""

import logging
from pathlib import Path

import pytest

from vntyper import cli
from vntyper.scripts import cli_calibrate
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


_COMMANDS = {
    "extract": ["--truth", "truth.json", "--partitions", "partitions.json", "--runs", "runs.json"],
    "fit": ["--evidence", "evidence", "--objective", "lexicographic-safety-v1"],
    "validate": ["--profile", "profile.json", "--evidence", "validation"],
    "evaluate": ["--profile", "profile.json", "--evidence", "locked-heldout"],
}


@pytest.mark.parametrize(("operation", "options"), sorted(_COMMANDS.items()))
def test_exact_four_calibration_operations_parse(operation: str, options: list[str]) -> None:
    args = build_parser().parse_args(["calibrate", operation, *options, "--output", "out"])

    assert args.command == "calibrate"
    assert args.calibration_operation == operation
    assert args.output == Path("out")


def test_fit_objective_is_mandatory_and_closed() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit) as missing:
        parser.parse_args(["calibrate", "fit", "--evidence", "evidence", "--output", "out"])
    with pytest.raises(SystemExit) as unknown:
        parser.parse_args(["calibrate", "fit", "--evidence", "evidence", "--objective", "f1", "--output", "out"])

    assert missing.value.code == 2
    assert unknown.value.code == 2


@pytest.mark.parametrize("operation", sorted(_COMMANDS))
def test_every_operation_atomically_installs_only_complete_output(
    operation: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "result"

    def produce(_args, staging: Path) -> bool:
        assert staging.parent == output.parent
        (staging / "complete.json").write_text("{}\n", encoding="utf-8")
        return True

    monkeypatch.setitem(cli_calibrate.OPERATIONS, operation, produce)
    args = build_parser().parse_args(
        ["calibrate", operation, *_path_options(operation, tmp_path), "--output", str(output)]
    )
    cli_calibrate.handle_calibrate(args, {}, build_parser(), logging.INFO, None)

    assert tuple(path.name for path in output.iterdir()) == ("complete.json",)
    assert not tuple(tmp_path.glob(".result.*"))


def test_failed_operation_leaves_no_partial_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    output = tmp_path / "candidate"

    def fail(_args, staging: Path) -> None:
        (staging / "partial.json").write_text("partial\n", encoding="utf-8")
        raise ValueError("snapshotted objective differs")

    monkeypatch.setitem(cli_calibrate.OPERATIONS, "fit", fail)
    args = build_parser().parse_args(["calibrate", "fit", *_path_options("fit", tmp_path), "--output", str(output)])

    with pytest.raises(ValueError, match="objective differs"):
        cli_calibrate.handle_calibrate(args, {}, build_parser(), logging.INFO, None)
    assert not output.exists()
    assert not tuple(tmp_path.glob(".candidate.*"))


@pytest.mark.parametrize("interruption", [KeyboardInterrupt(), SystemExit(9)])
def test_atomic_output_cleans_exact_staging_on_base_exception(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, interruption: BaseException
) -> None:
    output = tmp_path / "candidate"
    neighbour = tmp_path / ".candidate.keep"
    neighbour.mkdir()

    def interrupt(_args, staging: Path) -> bool:
        (staging / "partial.json").write_text("partial\n", encoding="utf-8")
        raise interruption

    monkeypatch.setitem(cli_calibrate.OPERATIONS, "fit", interrupt)
    args = build_parser().parse_args(["calibrate", "fit", *_path_options("fit", tmp_path), "--output", str(output)])

    with pytest.raises(type(interruption)):
        cli_calibrate.handle_calibrate(args, {}, build_parser(), logging.INFO, None)

    assert neighbour.is_dir()
    assert not output.exists()
    assert set(tmp_path.glob(".candidate.*")) == {neighbour}


def test_completed_failed_operation_is_installed_before_cli_exit_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "validation"

    def completed_failure(_args, staging: Path) -> bool:
        (staging / "attestation.json").write_text('{"status":"failed"}\n', encoding="utf-8")
        (staging / "retirement.json").write_text('{"reason":"failed-validation"}\n', encoding="utf-8")
        return False

    monkeypatch.setitem(cli_calibrate.OPERATIONS, "validate", completed_failure)
    args = build_parser().parse_args(
        ["calibrate", "validate", *_path_options("validate", tmp_path), "--output", str(output)]
    )

    with pytest.raises(SystemExit) as failure:
        cli_calibrate.handle_calibrate(args, {}, build_parser(), logging.INFO, None)

    assert failure.value.code == 1
    assert (output / "attestation.json").is_file()
    assert (output / "retirement.json").is_file()


def test_main_maps_completed_failed_operation_to_one_after_install(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "validation"

    def completed_failure(_args, staging: Path) -> bool:
        (staging / "attestation.json").write_text('{"status":"failed"}\n', encoding="utf-8")
        (staging / "retirement.json").write_text('{"reason":"failed-validation"}\n', encoding="utf-8")
        return False

    monkeypatch.setitem(cli_calibrate.OPERATIONS, "validate", completed_failure)
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)

    with pytest.raises(SystemExit) as failure:
        cli.main(
            [
                "calibrate",
                "validate",
                "--profile",
                str(tmp_path / "profile.json"),
                "--evidence",
                str(tmp_path / "evidence"),
                "--output",
                str(output),
            ]
        )

    assert failure.value.code == 1
    assert (output / "attestation.json").is_file()
    assert (output / "retirement.json").is_file()


def test_objective_mismatch_exits_one_without_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    output = tmp_path / "candidate"

    def fail(_args, _staging: Path) -> None:
        raise ValueError("fit objective differs from snapshotted protocol")

    monkeypatch.setitem(cli_calibrate.OPERATIONS, "fit", fail)
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)

    with pytest.raises(SystemExit) as failure:
        cli.main(
            [
                "calibrate",
                "fit",
                "--evidence",
                str(tmp_path / "evidence"),
                "--objective",
                "lexicographic-safety-v1",
                "--output",
                str(output),
            ]
        )

    assert failure.value.code == 1
    assert not output.exists()


def _path_options(operation: str, root: Path) -> list[str]:
    if operation == "extract":
        return [
            "--truth",
            str(root / "truth.json"),
            "--partitions",
            str(root / "partitions.json"),
            "--runs",
            str(root / "runs.json"),
        ]
    if operation == "fit":
        return ["--evidence", str(root / "evidence"), "--objective", "lexicographic-safety-v1"]
    return ["--profile", str(root / "profile.json"), "--evidence", str(root / operation)]
