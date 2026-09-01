"""Resolved adVNTR settings reach command and pathogenic-frame decisions."""

from __future__ import annotations

import shlex
from pathlib import Path

import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def test_explicit_settings_control_the_command_without_module_state(tmp_path: Path, monkeypatch) -> None:
    database = tmp_path / "model.db"
    alignment = tmp_path / "sample.bam"
    output = tmp_path / "out"
    database.touch()
    alignment.touch()
    output.mkdir()
    commands: list[str] = []

    def record_command(command, *_args, **_kwargs) -> bool:
        commands.append(command)
        return True

    monkeypatch.setattr(advntr, "run_command", record_command)
    run = resolve_run_configuration()
    component = dict(run.advntr)
    component["settings"] = {
        "frameshift_multiplier": 3,
        "max_frameshift": 100,
        "output_format": "tsv",
        "vid": 4242,
    }
    runtime = {"settings": {"additional_commands": "-aln", "threads": 7}}

    result = advntr.run_advntr(
        str(database),
        str(alignment),
        str(output),
        "sample",
        {"tools": {"advntr": "advntr"}},
        resolved_component=component,
        runtime_component=runtime,
        custom_context_active=True,
    )

    assert result == 0
    words = shlex.split(commands[0])
    assert words[words.index("-vid") + 1] == "4242"
    assert words[words.index("-t") + 1] == "7"
    assert words[words.index("-o") + 1].endswith("sample_adVNTR.tsv")
    assert "-aln" in words


def test_custom_command_context_cannot_fall_back_to_packaged_decisions(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="custom adVNTR run context requires an explicit resolved component"):
        advntr.run_advntr(
            str(tmp_path / "model.db"),
            str(tmp_path / "sample.bam"),
            str(tmp_path / "out"),
            "sample",
            {"tools": {"advntr": "advntr"}},
            custom_context_active=True,
        )


def test_explicit_frameshift_settings_control_the_accepted_series() -> None:
    settings = {"frameshift_multiplier": 5, "max_frameshift": 3}

    result = advntr.accepted_frame_magnitudes(advntr.INSERTION_FRAME_OFFSET, settings=settings)

    assert result.tolist() == ["1", "6", "11"]
