"""``handle_pipeline`` refuses to start when a Kestrel file is missing (#338).

The other handler tests stub the preflight so their synthetic configs stay about what
they test. This module is where it runs for real, against a config whose files do or do
not exist, and asserts that ``run_pipeline`` is never reached when one is missing.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts import cli_handlers
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


def _config(root: Path) -> dict:
    """A config whose Kestrel inputs are absolute paths under ``root``."""
    return {
        "tools": {
            "kestrel": str(root / "kestrel.jar"),
            "kanalyze": str(root / "kanalyze.jar"),
            "java_path": "java",
            "samtools": "samtools",
        },
        "reference_data": {
            "muc1_reference_vntr": str(root / "reference" / "motifs.fa"),
            "muc1_motifs_rev_com": str(root / "reference" / "rev_com.fa"),
        },
    }


def _install_all(config: dict) -> None:
    for section in ("tools", "reference_data"):
        for key in ("kestrel", "kanalyze", "muc1_reference_vntr", "muc1_motifs_rev_com"):
            value = config[section].get(key)
            if value:
                Path(value).parent.mkdir(parents=True, exist_ok=True)
                Path(value).write_text("x", encoding="utf-8")


def _handle(tmp_path: Path, config: dict) -> mock.MagicMock:
    parser = build_parser()
    args = parser.parse_args(["pipeline", "-o", str(tmp_path / "out"), "--cram", "in.cram", "--fast-mode"])
    with (
        mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as runner,
        mock.patch.object(cli_handlers, "_resolve_bwa_reference", return_value=None),
    ):
        cli_handlers.handle_pipeline(
            args, config=config, parser=parser, log_level_value=logging.INFO, log_file_str=None
        )
    return runner


def test_a_missing_motif_reference_stops_the_run_before_the_pipeline(tmp_path):
    config = _config(tmp_path)
    _install_all(config)
    Path(config["reference_data"]["muc1_reference_vntr"]).unlink()

    parser = build_parser()
    args = parser.parse_args(["pipeline", "-o", str(tmp_path / "out"), "--cram", "in.cram", "--fast-mode"])
    with (
        mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as runner,
        mock.patch.object(cli_handlers, "_resolve_bwa_reference", return_value=None),
        pytest.raises(ValueError, match="muc1_reference_vntr") as excinfo,
    ):
        cli_handlers.handle_pipeline(
            args, config=config, parser=parser, log_level_value=logging.INFO, log_file_str=None
        )

    runner.assert_not_called()
    assert f"current working directory ({os.getcwd()})" in str(excinfo.value)


def test_a_complete_install_reaches_the_pipeline(tmp_path):
    config = _config(tmp_path)
    _install_all(config)
    runner = _handle(tmp_path, config)
    assert runner.call_count == 1
