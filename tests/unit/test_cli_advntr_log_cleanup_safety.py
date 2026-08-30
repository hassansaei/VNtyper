"""The CLI must protect its active log from adVNTR preflight cleanup."""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

from vntyper import cli

pytestmark = pytest.mark.unit


def _config() -> dict[str, object]:
    """Return the configuration values resolved before pipeline dispatch."""
    return {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19", "archive_format": "zip"},
        "reference_data": {},
    }


@pytest.mark.parametrize(
    ("case", "module_spelling"),
    [
        ("exact-output", "advntr"),
        ("redundant-spelling", " ADvNTR "),
        ("output-alias", "shark,AdVnTr"),
        ("selected-archive", "ADVNTR"),
    ],
)
def test_advntr_cleanup_log_collision_is_rejected_before_logging_or_dispatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    case: str,
    module_spelling: str,
) -> None:
    """Every accepted path/module spelling reaches the same pre-open refusal."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    bam = input_dir / "sample.bam"
    bam.write_bytes(b"operator input")
    actual_output = tmp_path / "output"
    output_argument = actual_output
    archive_arguments: list[str] = []

    if case == "exact-output":
        log_file = actual_output / "pipeline_summary.json"
    elif case == "redundant-spelling":
        (actual_output / "existing").mkdir(parents=True)
        log_file = actual_output / "existing" / ".." / "summary_report.html"
    elif case == "output-alias":
        actual_output.mkdir()
        output_argument = tmp_path / "output-alias"
        output_argument.symlink_to(actual_output, target_is_directory=True)
        log_file = actual_output / "pipeline_summary.csv"
    else:
        log_file = Path(f"{actual_output}.tar.gz")
        archive_arguments = ["--archive-results", "--archive-format", "tar.gz"]

    setup = mock.Mock()
    handler = mock.Mock()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: _config())
    monkeypatch.setattr(cli, "setup_logging", setup)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handler)

    with pytest.raises(SystemExit) as raised:
        cli.main(
            [
                "--log-file",
                str(log_file),
                "pipeline",
                "-o",
                str(output_argument),
                "--bam",
                str(bam),
                "--extra-modules",
                module_spelling,
                *archive_arguments,
            ]
        )

    assert raised.value.code == 1
    setup.assert_not_called()
    handler.assert_not_called()
    assert not log_file.exists()


def test_ordinary_advntr_pipeline_log_keeps_setup_before_dispatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The normal pipeline.log is outside the exact cleanup destination set."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    bam = input_dir / "sample.bam"
    bam.write_bytes(b"operator input")
    output = tmp_path / "output"
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: _config())
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(
        [
            "--log-file",
            str(output / "pipeline.log"),
            "pipeline",
            "-o",
            str(output),
            "--bam",
            str(bam),
            "--extra-modules",
            "advntr",
            "--archive-results",
        ]
    )

    assert events == ["setup", "handler"]


def test_malformed_archive_policy_remains_owned_by_pipeline_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The pre-open log guard must not take over unrelated archive validation."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    bam = input_dir / "sample.bam"
    bam.write_bytes(b"operator input")
    output = tmp_path / "output"
    config = _config()
    defaults = config["default_values"]
    assert isinstance(defaults, dict)
    defaults["archive_format"] = "rar"
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(
        [
            "--log-file",
            str(output / "pipeline.log"),
            "pipeline",
            "-o",
            str(output),
            "--bam",
            str(bam),
            "--extra-modules",
            "advntr",
            "--archive-results",
        ]
    )

    assert events == ["setup", "handler"]
