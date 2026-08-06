# tests/unit/test_advntr_command.py

"""
Unit tests for :func:`vntyper.modules.advntr.advntr_genotyping.run_advntr`.

``run_advntr`` builds a single shell string and hands it to ``run_command`` (which runs it
under ``shell=True``). Nothing downstream re-derives that string, so the tests here pin it
literally, along with the artefact name it tells adVNTR to write -- ``pipeline.py``
reconstructs that name independently and reads the file back (wave-2 contract C4).

They also cover every early return, because each one is the difference between "adVNTR did
not run" and a confusing failure several stages later.
"""

import logging
import subprocess as sp
from pathlib import Path

import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr

pytestmark = pytest.mark.unit

MAIN_CONFIG = {"tools": {"advntr": "mamba run -n envadvntr advntr"}}


@pytest.fixture
def inputs(tmp_path: Path):
    """A database file, a sorted BAM and an output directory that all exist."""
    db_file = tmp_path / "vntr_db.db"
    db_file.write_text("db")
    sorted_bam = tmp_path / "sample.sorted.bam"
    sorted_bam.write_text("bam")
    output = tmp_path / "advntr"
    output.mkdir()
    return db_file, sorted_bam, output


@pytest.fixture
def captured_command(monkeypatch):
    """Replace ``run_command`` with a spy and return the list of calls it recorded."""
    calls: list[dict] = []

    def spy(command, log_file, critical=False, cwd=None):
        calls.append({"command": command, "log_file": log_file, "critical": critical, "cwd": cwd})
        return True

    monkeypatch.setattr(advntr, "run_command", spy)
    return calls


class TestTheCommandLine:
    def test_the_command_is_built_exactly_as_advntr_expects(self, inputs, captured_command):
        db_file, sorted_bam, output = inputs

        result = advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert result == 0
        assert captured_command[0]["command"] == (
            f"mamba run -n envadvntr advntr genotype -fs -vid 25561 "
            f"--alignment_file {sorted_bam} -o {output}/output_adVNTR.vcf "
            f"-m {db_file} --working_directory {output} -t 1 -aln"
        )

    def test_the_artefact_name_is_built_from_output_name(self, inputs, captured_command):
        """Contract C4: ``pipeline.py`` rebuilds this path rather than consuming a return value."""
        db_file, sorted_bam, output = inputs

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "sample_42", MAIN_CONFIG)

        assert f"-o {output}/sample_42_adVNTR.vcf" in captured_command[0]["command"]

    def test_the_log_file_sits_beside_the_result(self, inputs, captured_command):
        db_file, sorted_bam, output = inputs

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "sample_42", MAIN_CONFIG)

        assert captured_command[0]["log_file"] == str(output / "sample_42_advntr.log")

    def test_the_command_is_critical_and_carries_the_working_directory(self, inputs, captured_command, tmp_path):
        db_file, sorted_bam, output = inputs

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG, cwd=str(tmp_path))

        assert captured_command[0]["critical"] is True
        assert captured_command[0]["cwd"] == str(tmp_path)


class TestSettingsComeFromTheDerivedGlobal:
    """H1 again: ``run_advntr`` reads ``advntr_settings``, never its own ``config`` argument."""

    @pytest.mark.parametrize(
        ("output_format", "expected_extension"),
        [("vcf", ".vcf"), ("tsv", ".tsv"), ("anything-else", ".tsv")],
    )
    def test_the_extension_follows_the_configured_output_format(
        self, inputs, captured_command, monkeypatch, output_format, expected_extension
    ):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(advntr, "advntr_settings", {"output_format": output_format})

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert f"-o {output}/output_adVNTR{expected_extension}" in captured_command[0]["command"]

    def test_threads_vid_and_extra_arguments_come_from_the_settings_global(self, inputs, captured_command, monkeypatch):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"threads": 12, "vid": 999, "additional_commands": "--haploid", "output_format": "tsv"},
        )

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        command = captured_command[0]["command"]
        assert "-vid 999" in command
        assert "-t 12" in command
        assert command.endswith("--haploid")

    def test_patching_the_raw_config_global_does_not_change_the_command(self, inputs, captured_command, monkeypatch):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(advntr, "advntr_config", {"advntr_settings": {"vid": 999, "threads": 12}})

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert "-vid 25561" in captured_command[0]["command"]

    def test_the_defaults_apply_when_the_settings_global_is_empty(self, inputs, captured_command, monkeypatch):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(advntr, "advntr_settings", {})

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        command = captured_command[0]["command"]
        assert "-vid 25561" in command
        assert "-t 8" in command
        assert command.endswith("-aln")
        assert "output_adVNTR.tsv" in command, "the default output_format is tsv"


class TestInputValidation:
    def test_a_missing_database_file_returns_one_without_running_anything(self, tmp_path, captured_command, caplog):
        sorted_bam = tmp_path / "sample.bam"
        sorted_bam.write_text("bam")

        with caplog.at_level(logging.CRITICAL, logger=advntr.logger.name):
            result = advntr.run_advntr(str(tmp_path / "missing.db"), str(sorted_bam), str(tmp_path), "o", MAIN_CONFIG)

        assert result == 1
        assert captured_command == []
        assert [r for r in caplog.records if r.levelno >= logging.CRITICAL]

    def test_a_missing_bam_returns_one_without_running_anything(self, tmp_path, captured_command, caplog):
        db_file = tmp_path / "db.db"
        db_file.write_text("db")

        with caplog.at_level(logging.CRITICAL, logger=advntr.logger.name):
            result = advntr.run_advntr(str(db_file), str(tmp_path / "missing.bam"), str(tmp_path), "o", MAIN_CONFIG)

        assert result == 1
        assert captured_command == []

    def test_a_missing_output_directory_is_created(self, tmp_path, captured_command, caplog):
        db_file = tmp_path / "db.db"
        db_file.write_text("db")
        sorted_bam = tmp_path / "sample.bam"
        sorted_bam.write_text("bam")
        output = tmp_path / "not_yet_there"

        with caplog.at_level(logging.WARNING, logger=advntr.logger.name):
            result = advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert result == 0
        assert output.is_dir()
        assert [r for r in caplog.records if r.levelno == logging.WARNING]

    def test_an_uncreatable_output_directory_returns_one(self, tmp_path, monkeypatch, captured_command, caplog):
        db_file = tmp_path / "db.db"
        db_file.write_text("db")
        sorted_bam = tmp_path / "sample.bam"
        sorted_bam.write_text("bam")

        def boom(*_args, **_kwargs):
            raise PermissionError("read-only file system")

        monkeypatch.setattr(advntr.os, "makedirs", boom)

        with caplog.at_level(logging.CRITICAL, logger=advntr.logger.name):
            result = advntr.run_advntr(str(db_file), str(sorted_bam), str(tmp_path / "nope"), "o", MAIN_CONFIG)

        assert result == 1
        assert captured_command == []

    def test_a_config_without_a_tool_path_raises_rather_than_running_a_bare_name(self, inputs):
        """``config["tools"]["advntr"]`` has no ``.get`` fallback -- trap 2 in AGENTS.md."""
        db_file, sorted_bam, output = inputs

        with pytest.raises(KeyError):
            advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", {"tools": {}})


class TestFailureHandling:
    @pytest.mark.parametrize(
        "outcome",
        [
            pytest.param(False, id="run_command_reports_failure"),
            pytest.param(sp.CalledProcessError(1, "advntr"), id="called_process_error"),
            pytest.param(RuntimeError("critical command failed"), id="unexpected_error"),
        ],
    )
    def test_every_failure_mode_returns_one_and_logs_at_error(self, inputs, monkeypatch, caplog, outcome):
        db_file, sorted_bam, output = inputs

        def failing(*_args, **_kwargs):
            if isinstance(outcome, BaseException):
                raise outcome
            return outcome

        monkeypatch.setattr(advntr, "run_command", failing)

        with caplog.at_level(logging.ERROR, logger=advntr.logger.name):
            result = advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert result == 1
        assert [r for r in caplog.records if r.levelno >= logging.ERROR]


class TestConfigLoading:
    def test_the_default_config_is_the_one_shipped_beside_the_module(self):
        loaded = advntr.load_advntr_config()

        assert set(loaded) == {"advntr_settings", "flagging_rules"}
        assert loaded["advntr_settings"]["vid"] == 25561

    def test_an_explicit_path_is_honoured(self, tmp_path):
        path = tmp_path / "custom.json"
        path.write_text('{"advntr_settings": {"vid": 1}}')

        assert advntr.load_advntr_config(str(path)) == {"advntr_settings": {"vid": 1}}
