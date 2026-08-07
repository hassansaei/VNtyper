# tests/unit/test_shark_filtering.py

"""
Unit tests for :func:`vntyper.modules.shark.shark_filtering.run_shark_filter`.

SHARK is an optional pre-filter: when ``--extra-modules shark`` is given for a FASTQ input,
``pipeline.py`` runs it *before* quality control and then feeds its two output FASTQs into
the rest of the pipeline. Whatever SHARK discards is invisible to Kestrel, so this stage
silently bounds every downstream call.

Two things are pinned here:

* the command string and the two artefact names, since ``pipeline.py`` consumes the
  returned paths and separately records a summary artefact path of its own;
* the fact that ``reference_assembly`` is accepted and then never read -- see
  :class:`TestReferenceAssemblyIsAccceptedAndIgnored`, which is characterisation of a live
  defect, deliberately not a fix.
"""

import logging
import shlex
from pathlib import Path

import pytest

from vntyper.modules.shark import shark_filtering as shark

pytestmark = pytest.mark.unit

MAIN_CONFIG = {"tools": {"shark": "mamba run -n shark_env shark"}}
SHARK_CONFIG = {"shark_settings": {"muc1_region_fasta": "reference/muc1_region_hg19.fa"}}


@pytest.fixture
def captured_command(monkeypatch):
    """Replace ``run_command`` with a spy and return the list of calls it recorded."""
    calls: list[dict] = []

    def spy(command, log_file, critical=False, cwd=None):
        calls.append({"command": command, "log_file": log_file, "critical": critical})
        return True

    monkeypatch.setattr(shark, "run_command", spy)
    return calls


def filter_with(output_dir: Path, **overrides):
    """Call ``run_shark_filter`` with sensible defaults, overriding as needed."""
    kwargs = {
        "fastq_1": "/data/sample_R1.fastq.gz",
        "fastq_2": "/data/sample_R2.fastq.gz",
        "output_dir": str(output_dir),
        "config": SHARK_CONFIG,
        "main_config": MAIN_CONFIG,
        "sample_name": "sample",
        "reference_assembly": "hg19",
        "threads": 4,
    }
    kwargs.update(overrides)
    return shark.run_shark_filter(**kwargs)


class TestTheCommandLine:
    def test_the_command_is_built_exactly_as_shark_expects(self, tmp_path, captured_command):
        filter_with(tmp_path)

        assert captured_command[0]["command"] == (
            "mamba run -n shark_env shark -r reference/muc1_region_hg19.fa "
            "-1 /data/sample_R1.fastq.gz -2 /data/sample_R2.fastq.gz "
            f"-o {tmp_path}/sample_shark_R1.fastq -p {tmp_path}/sample_shark_R2.fastq -t 4"
        )

    def test_the_two_filtered_fastqs_are_named_after_the_sample(self, tmp_path, captured_command):
        first, second = filter_with(tmp_path, sample_name="patient_7")

        assert first == str(tmp_path / "patient_7_shark_R1.fastq")
        assert second == str(tmp_path / "patient_7_shark_R2.fastq")

    def test_the_returned_paths_are_the_ones_named_on_the_command_line(self, tmp_path, captured_command):
        """``pipeline.py`` reassigns ``fastq1``/``fastq2`` to these, so they must agree."""
        first, second = filter_with(tmp_path)

        assert f"-o {first} -p {second}" in captured_command[0]["command"]

    def test_the_log_file_sits_beside_the_output(self, tmp_path, captured_command):
        filter_with(tmp_path, sample_name="patient_7")

        assert captured_command[0]["log_file"] == str(tmp_path / "patient_7_shark.log")

    def test_the_command_is_critical(self, tmp_path, captured_command):
        filter_with(tmp_path)

        assert captured_command[0]["critical"] is True

    def test_threads_are_passed_through(self, tmp_path, captured_command):
        filter_with(tmp_path, threads=16)

        assert captured_command[0]["command"].endswith("-t 16")

    def test_the_tool_name_falls_back_to_a_bare_shark(self, tmp_path, captured_command):
        filter_with(tmp_path, main_config={})

        assert captured_command[0]["command"].startswith("shark -r ")

    def test_the_invocation_is_logged_at_info(self, tmp_path, captured_command, caplog):
        with caplog.at_level(logging.INFO, logger=shark.logger.name):
            filter_with(tmp_path)

        messages = [r.getMessage() for r in caplog.records if r.levelno >= logging.INFO]
        assert any("Running SHARK filtering" in message for message in messages)


class TestConfigResolution:
    """
    Contrast with the adVNTR module (hazard H1): ``run_shark_filter`` really does read its
    ``config`` argument, so a test that patches the module global is the one that proves
    nothing here.
    """

    def test_the_region_fasta_comes_from_the_config_argument(self, tmp_path, captured_command):
        custom = {"shark_settings": {"muc1_region_fasta": "/refs/custom_region.fa"}}

        filter_with(tmp_path, config=custom)

        assert "-r /refs/custom_region.fa " in captured_command[0]["command"]

    def test_patching_the_module_global_has_no_effect(self, tmp_path, captured_command, monkeypatch):
        monkeypatch.setattr(shark, "shark_settings", {"muc1_region_fasta": "/refs/patched.fa"})

        filter_with(tmp_path)

        assert "/refs/patched.fa" not in captured_command[0]["command"]

    @pytest.mark.parametrize(
        "config",
        [
            pytest.param({}, id="no_shark_settings"),
            pytest.param({"shark_settings": {}}, id="empty_shark_settings"),
            pytest.param({"shark_settings": {"muc1_region_fasta": ""}}, id="blank_region_fasta"),
        ],
    )
    def test_a_config_without_a_region_fasta_raises_before_running_anything(self, tmp_path, captured_command, config):
        with pytest.raises(ValueError, match="muc1_region_fasta"):
            filter_with(tmp_path, config=config)

        assert captured_command == []

    def test_load_shark_config_reads_the_file_shipped_beside_the_module(self):
        loaded = shark.load_shark_config()

        assert loaded["shark_settings"]["muc1_region_fasta"].endswith("muc1_region_hg19.fa")

    def test_an_explicit_config_path_is_honoured(self, tmp_path):
        path = tmp_path / "custom.json"
        path.write_text('{"shark_settings": {"muc1_region_fasta": "/refs/x.fa"}}')

        assert shark.load_shark_config(str(path)) == {"shark_settings": {"muc1_region_fasta": "/refs/x.fa"}}


class TestFailureHandling:
    def test_a_failed_run_raises_runtime_error_and_logs_at_error(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(shark, "run_command", lambda *a, **k: False)

        with caplog.at_level(logging.ERROR, logger=shark.logger.name), pytest.raises(RuntimeError, match="SHARK"):
            filter_with(tmp_path)

        assert [r for r in caplog.records if r.levelno >= logging.ERROR]

    def test_a_runtime_error_from_run_command_propagates(self, tmp_path, monkeypatch):
        def boom(*_args, **_kwargs):
            raise RuntimeError("Critical command failed: shark ...")

        monkeypatch.setattr(shark, "run_command", boom)

        with pytest.raises(RuntimeError):
            filter_with(tmp_path)


class TestReferenceAssemblyIsAccceptedAndIgnored:
    """
    **Specification (#187), not a live defect.**

    @hassansaei on #187: "SHARK is sequence-based, not coordinate-based; keep one MUC1
    region FASTA; document that assembly does not change SHARK; optionally warn or
    ignore reference_assembly there instead of building a second FASTA unless
    GRCh37/GRCh38 sequences at MUC1 VNTR differ enough to matter, like the number of
    motif could impact the number of reads being retained by the SHARK model."

    ``run_shark_filter`` takes a ``reference_assembly`` argument and ``pipeline.py`` passes
    the run's actual assembly into it, but the function does not use it to select a
    region: the region FASTA comes solely from ``shark_config.json``, which declares a
    single MUC1 region file. That is intentional, not a bug -- SHARK matches k-mers
    against sequence, not coordinates, so there is nothing for the assembly to select.
    Per the decision above, no hg38 region FASTA is added. The parameter's docstring now
    reads:

        reference_assembly (str): Accepted for API compatibility and **ignored**.
        SHARK matches k-mers against a single MUC1 region FASTA and does not select a
        region by coordinate, so the assembly does not change what it retains. Passing
        anything other than hg19/GRCh37 logs a warning. See issue #187.

    See :class:`TestNonHg19AssemblyLogsAWarning` for that warning.
    """

    def test_hg38_produces_a_byte_identical_command_to_hg19(self, tmp_path, captured_command):
        filter_with(tmp_path, reference_assembly="hg19")
        filter_with(tmp_path, reference_assembly="hg38")

        assert captured_command[0]["command"] == captured_command[1]["command"]

    def test_an_hg38_run_still_filters_against_the_hg19_region(self, tmp_path, captured_command):
        filter_with(tmp_path, reference_assembly="hg38")

        assert "muc1_region_hg19.fa" in captured_command[0]["command"]

    def test_even_a_nonsense_assembly_is_accepted_and_warned_about(self, tmp_path, captured_command, caplog):
        """Specification (#187): a nonsense value is still accepted -- the run proceeds
        regardless -- but no longer passes silently. Any value other than hg19/GRCh37
        now logs a warning naming it.
        """
        with caplog.at_level(logging.WARNING, logger=shark.logger.name):
            filter_with(tmp_path, reference_assembly="not-an-assembly")

        assert captured_command, "the run proceeds regardless"
        messages = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
        assert any("reference_assembly" in m and "not-an-assembly" in m for m in messages), "and now says so"

    def test_the_shipped_config_offers_only_an_hg19_region(self):
        """There is no hg38 alternative to select, so this cannot be fixed in-module."""
        settings = shark.load_shark_config()["shark_settings"]

        assert [key for key in settings if "fasta" in key] == ["muc1_region_fasta"]
        assert "hg38" not in str(settings)


class TestNonHg19AssemblyLogsAWarning:
    """
    Specification (#187). A non-hg19/GRCh37 ``reference_assembly`` does not change which
    region SHARK filters against -- there is only ever the one MUC1 region FASTA -- but
    it now logs a warning naming the value, so a caller who passes hg38 finds out the
    parameter did nothing instead of silently getting hg19-filtered reads (see
    :class:`TestReferenceAssemblyIsAccceptedAndIgnored`).
    """

    @pytest.mark.parametrize("assembly", ["hg38", "GRCh38", "hg38_ensembl"])
    def test_a_non_hg19_assembly_is_warned_about(self, assembly, tmp_path, captured_command, caplog):
        with caplog.at_level(logging.WARNING, logger=shark.logger.name):
            filter_with(tmp_path, reference_assembly=assembly)

        messages = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
        assert any("reference_assembly" in m and assembly in m for m in messages)

    @pytest.mark.parametrize("assembly", ["hg19", "GRCh37", "grch37"])
    def test_hg19_and_grch37_do_not_warn(self, assembly, tmp_path, captured_command, caplog):
        """The warning marks a mismatch, not every call: hg19 and GRCh37 (any case) are silent."""
        with caplog.at_level(logging.WARNING, logger=shark.logger.name):
            filter_with(tmp_path, reference_assembly=assembly)

        messages = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
        assert not [m for m in messages if "reference_assembly" in m]

    def test_the_region_fasta_operand_is_byte_identical_at_hg19_and_hg38(self, tmp_path, captured_command):
        """The substantive claim: assembly does not select a region.

        A warning alone would not prove the parameter is inert; capturing the command
        SHARK is actually given at hg19 and at hg38 and comparing the ``-r`` operand
        does. (``TestReferenceAssemblyIsAccceptedAndIgnored`` already pins the whole
        command as identical; this asserts the specific operand SHARK reads the region
        from, so the claim holds even if unrelated parts of the command ever diverge.)
        """
        filter_with(tmp_path, reference_assembly="hg19")
        filter_with(tmp_path, reference_assembly="hg38")

        def region_operand(command):
            words = shlex.split(command)
            return words[words.index("-r") + 1]

        hg19_operand = region_operand(captured_command[0]["command"])
        hg38_operand = region_operand(captured_command[1]["command"])
        assert hg19_operand == hg38_operand


class TestShellInterpolationIsQuoted:
    """
    ``run_command`` executes with ``shell=True`` (trap 9), which makes quoting the
    caller's job. This caller used to interpolate raw, so a path containing a space
    silently became two arguments and a path containing a ``;`` was arbitrary command
    execution -- and the web service accepts filenames that reach this line.

    ``shark_path`` stays unquoted on purpose: ``config.json`` holds a command *prefix*
    there ("mamba run -n shark_env shark"), which quoting would collapse into one token.
    """

    @pytest.mark.parametrize(
        ("flag", "argument", "hostile"),
        [
            ("-1", "fastq_1", "/data/my samples/R1.fastq.gz"),
            ("-2", "fastq_2", "/data/my samples/R2.fastq.gz"),
            ("-1", "fastq_1", "/data/o'brien/R1.fastq.gz"),
            ("-1", "fastq_1", "/data/x; touch /tmp/pwned/R1.fastq.gz"),
            ("-1", "fastq_1", "/data/$(whoami)/R1.fastq.gz"),
            ("-1", "fastq_1", "/data/`id`/R1.fastq.gz"),
        ],
    )
    def test_a_hostile_path_survives_as_one_shell_word(self, tmp_path, captured_command, flag, argument, hostile):
        filter_with(tmp_path, **{argument: hostile})

        command = captured_command[0]["command"]
        words = shlex.split(command)
        assert words[words.index(flag) + 1] == hostile, f"{hostile!r} did not survive as one word: {command}"

    def test_a_sample_name_with_a_space_survives_in_both_output_paths(self, tmp_path, captured_command):
        filter_with(tmp_path, sample_name="patient 7")

        words = shlex.split(captured_command[0]["command"])
        assert words[words.index("-o") + 1] == str(tmp_path / "patient 7_shark_R1.fastq")
        assert words[words.index("-p") + 1] == str(tmp_path / "patient 7_shark_R2.fastq")

    def test_the_tool_prefix_is_still_three_separate_words(self, tmp_path, captured_command):
        """Quoting `shark_path` would make bash look for one binary called `mamba run ...`."""
        filter_with(tmp_path)

        assert captured_command[0]["command"].startswith("mamba run -n shark_env shark ")
