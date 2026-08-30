"""The CLI must protect its active log from adVNTR preflight cleanup."""

from __future__ import annotations

import argparse
from pathlib import Path
from unittest import mock

import pytest

from vntyper import cli
from vntyper.scripts.cli_logging_safety import validate_pipeline_log_destination

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
        ("model-snapshot", "advntr"),
        ("model-snapshot-output-alias", "AdVnTr"),
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
    elif case in {"output-alias", "model-snapshot-output-alias"}:
        actual_output.mkdir()
        output_argument = tmp_path / "output-alias"
        output_argument.symlink_to(actual_output, target_is_directory=True)
        log_file = (
            actual_output / "pipeline_summary.csv"
            if case == "output-alias"
            else actual_output / "advntr" / "advntr_model.db"
        )
    elif case == "selected-archive":
        log_file = Path(f"{actual_output}.tar.gz")
        archive_arguments = ["--archive-results", "--archive-format", "tar.gz"]
    else:
        log_file = actual_output / "advntr" / "advntr_model.db"

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


@pytest.mark.parametrize("model_route", ["exact", "symlinked-source", "named-override"])
def test_selected_operator_advntr_model_is_rejected_before_logging_or_dispatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    model_route: str,
) -> None:
    """The pre-open guard resolves the same model source as direct preflight."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    bam = input_dir / "sample.bam"
    bam.write_bytes(b"operator input")
    hg19_model = tmp_path / "models" / "hg19.db"
    hg38_model = tmp_path / "models" / "hg38.db"
    hg19_model.parent.mkdir()
    hg19_model.write_bytes(b"hg19 operator model")
    hg38_model.write_bytes(b"hg38 operator model")
    configured_hg19 = hg19_model
    selected_model = hg19_model
    module_spelling = "shark, AdVnTr" if model_route == "exact" else "advntr"
    if model_route == "symlinked-source":
        configured_hg19 = tmp_path / "models" / "hg19-link.db"
        configured_hg19.symlink_to(hg19_model)
    elif model_route == "named-override":
        selected_model = hg38_model

    config = _config()
    config["reference_data"] = {
        "advntr_reference_vntr_hg19": str(configured_hg19),
        "advntr_reference_vntr_hg38": str(hg38_model),
    }
    setup = mock.Mock()
    handler = mock.Mock()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", setup)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handler)

    parser = cli.build_parser()
    cli_arguments = [
        "--log-file",
        str(selected_model),
        "pipeline",
        "-o",
        str(tmp_path / "output"),
        "--bam",
        str(bam),
        "--extra-modules",
        module_spelling,
    ]
    if model_route != "exact":
        cli_arguments.extend(("--reference-assembly", "hg19"))
    parsed = parser.parse_args(cli_arguments)
    if model_route == "named-override":
        parsed.advntr_reference = "hg38"
    monkeypatch.setattr(parser, "parse_args", lambda _argv=None: parsed)
    monkeypatch.setattr(cli, "build_parser", lambda: parser)

    with pytest.raises(SystemExit) as raised:
        cli.main([])

    assert raised.value.code == 1
    setup.assert_not_called()
    handler.assert_not_called()
    assert selected_model.read_bytes().endswith(b"operator model")


@pytest.mark.parametrize(
    ("override", "configured_model"),
    [
        ("/literal/model.db", "/configured/model.db"),
        (["hg19"], "/configured/model.db"),
        (None, ["/malformed/model.db"]),
    ],
)
def test_invalid_advntr_model_selection_remains_owned_by_later_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    override: object,
    configured_model: object,
) -> None:
    """The pre-open guard must not raise or log the planner's configuration refusal."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    bam = input_dir / "sample.bam"
    bam.write_bytes(b"operator input")
    config = {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19"},
        "reference_data": {"advntr_reference_vntr_hg19": configured_model},
    }
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    parser = cli.build_parser()
    parsed = parser.parse_args(
        [
            "--log-file",
            str(tmp_path / "safe.log"),
            "pipeline",
            "--bam",
            str(bam),
            "--extra-modules",
            "advntr",
        ]
    )
    parsed.advntr_reference = override
    monkeypatch.setattr(parser, "parse_args", lambda _argv=None: parsed)
    monkeypatch.setattr(cli, "build_parser", lambda: parser)

    cli.main([])

    assert events == ["setup", "handler"]
    assert not [
        record
        for record in caplog.records
        if record.name == "vntyper.scripts.pipeline_advntr_preflight" and record.levelno >= 40
    ]


def test_non_advntr_log_may_equal_the_unused_configured_model(tmp_path: Path) -> None:
    """A model that this run will not read is not operator-owned by the adVNTR guard."""
    model = tmp_path / "model.db"
    model.write_bytes(b"unused model")
    args = argparse.Namespace(
        fastq1=None,
        fastq2=None,
        bam=None,
        cram=None,
        bed_file=None,
        reference_fasta=None,
        reference_assembly="hg19",
        output_dir=str(tmp_path / "output"),
        archive_results=False,
        archive_format="zip",
        extra_modules=["shark"],
    )
    config = {"reference_data": {"advntr_reference_vntr_hg19": str(model)}}

    validate_pipeline_log_destination(model, args, config)


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
