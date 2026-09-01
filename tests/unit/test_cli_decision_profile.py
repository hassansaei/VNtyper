"""CLI boundary tests for resolving a profile before any run artifact."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from vntyper import cli
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


def test_decision_profile_option_exists_only_on_pipeline() -> None:
    parser = build_parser()
    args = parser.parse_args(["pipeline", "--bam", "input.bam", "--decision-profile", "complete.json"])

    assert args.decision_profile == Path("complete.json")
    with pytest.raises(SystemExit) as error:
        parser.parse_args(["report", "--output-dir", "out", "--decision-profile", "complete.json"])
    assert error.value.code == 2


def test_invalid_profile_fails_before_output_or_log_creation(tmp_path: Path) -> None:
    output = tmp_path / "must-not-exist"
    bad = tmp_path / "bad.json"
    bad.write_text("{}", encoding="utf-8")

    with pytest.raises(SystemExit) as error:
        cli.main(
            [
                "pipeline",
                "--bam",
                "input.bam",
                "--output-dir",
                str(output),
                "--decision-profile",
                str(bad),
            ]
        )

    assert error.value.code == 1
    assert not output.exists()


def test_report_and_cohort_never_resolve_current_profile(monkeypatch: pytest.MonkeyPatch) -> None:
    def forbidden(_path: object = None) -> None:
        raise AssertionError("non-pipeline command resolved a current profile")

    monkeypatch.setattr(cli, "resolve_run_configuration", forbidden)
    parser = build_parser()
    assert not hasattr(parser.parse_args(["report", "--output-dir", "out"]), "decision_profile")
    assert not hasattr(
        parser.parse_args(["cohort", "--input-dir", "in", "--output-dir", "out"]),
        "decision_profile",
    )


def test_pipeline_resolves_profile_exactly_once_before_dispatch(monkeypatch: pytest.MonkeyPatch) -> None:
    order: list[str] = []
    sentinel = object()

    def resolve(path: object = None) -> object:
        assert path is None
        order.append("resolve")
        return sentinel

    def handle(
        args: argparse.Namespace,
        config: dict[str, object],
        parser: argparse.ArgumentParser,
        log_level_value: int,
        log_file_str: str | None,
    ) -> None:
        del config, parser, log_level_value, log_file_str
        assert args.run_configuration is sentinel
        order.append("handler")

    monkeypatch.setattr(cli, "resolve_run_configuration", resolve)
    monkeypatch.setattr(cli, "load_config", lambda _path=None: {"cli_defaults": {"log_file": None}})
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: None)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handle)

    cli.main(["pipeline", "--bam", "input.bam"])

    assert order == ["resolve", "handler"]
