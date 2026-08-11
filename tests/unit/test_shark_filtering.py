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
* ``reference_assembly`` now *selects* the MUC1 region FASTA -- see
  :class:`TestReferenceAssemblySelectsTheRegion` (#152) -- resolved through
  :func:`vntyper.modules.shark.shark_filtering.select_muc1_region_fasta`, with
  ``main_config["reference_data"]`` (what ``install-references`` writes) taking
  precedence over the shipped ``shark_config.json`` (:class:`TestTheInstalledTreeWins`).
"""

import logging
import shlex
from pathlib import Path

import pytest

from vntyper.modules.shark import shark_filtering as shark

pytestmark = pytest.mark.unit

MAIN_CONFIG = {"tools": {"shark": "mamba run -n shark_env shark"}}
SHARK_CONFIG = {
    "shark_settings": {
        "muc1_region_fasta_hg19": "reference/muc1_region_hg19.fa",
        "muc1_region_fasta_hg38": "reference/muc1_region_hg38.fa",
    }
}


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

        settings = loaded["shark_settings"]
        assert settings["muc1_region_fasta_hg19"].endswith("muc1_region_hg19.fa")
        assert settings["muc1_region_fasta_hg38"].endswith("muc1_region_hg38.fa")

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


class TestReferenceAssemblySelectsTheRegion:
    """#152. SHARK filters on exact k-mers, and 40.6% of the hg38 MUC1 region's
    canonical 17-mers are absent from the hg19 region - so an hg38 run against the
    hg19 reference cannot see them. Measured across seven cohort samples: the hg38
    reference retains 3.2-34.7% more reads. The parameter now selects.
    """

    @pytest.mark.parametrize(
        "assembly,expected",
        [
            ("hg19", "muc1_region_hg19.fa"),
            ("GRCh37", "muc1_region_hg19.fa"),
            ("hg19_ncbi", "muc1_region_hg19.fa"),
            ("hg19_ensembl", "muc1_region_hg19.fa"),
            ("hg38", "muc1_region_hg38.fa"),
            ("GRCh38", "muc1_region_hg38.fa"),
            ("hg38_ncbi", "muc1_region_hg38.fa"),
            ("hg38_ensembl", "muc1_region_hg38.fa"),
        ],
    )
    def test_each_assembly_gets_its_coordinate_systems_region(self, tmp_path, captured_command, assembly, expected):
        filter_with(tmp_path, reference_assembly=assembly)
        assert expected in captured_command[0]["command"]

    def test_no_warning_is_logged_for_a_supported_assembly(self, tmp_path, captured_command, caplog):
        with caplog.at_level(logging.WARNING):
            filter_with(tmp_path, reference_assembly="hg38")
        assert not [r for r in caplog.records if "reference_assembly" in r.getMessage()]


class TestTheInstalledTreeWins:
    """shark_config.json is a separate file from config.json, so install-references
    cannot update it. Without reference_data as the first layer, `install-references
    --output-dir X` would leave SHARK pointing at a CWD-relative path - #163's exact
    bug, reintroduced for SHARK.
    """

    def test_reference_data_overrides_the_shipped_shark_config(self, tmp_path, captured_command):
        main = {**MAIN_CONFIG, "reference_data": {"muc1_region_fasta_hg38": "/custom/refs/muc1_region_hg38.fa"}}
        filter_with(tmp_path, main_config=main, reference_assembly="hg38")
        assert "/custom/refs/muc1_region_hg38.fa" in captured_command[0]["command"]

    def test_a_present_null_in_reference_data_fails_closed(self, tmp_path):
        main = {**MAIN_CONFIG, "reference_data": {"muc1_region_fasta_hg38": None}}
        with pytest.raises(ValueError, match="hg38"):
            filter_with(tmp_path, main_config=main, reference_assembly="hg38")

    def test_a_structurally_legacy_config_still_works(self, tmp_path, captured_command):
        legacy = {"shark_settings": {"muc1_region_fasta": "/legacy/region.fa"}}
        filter_with(tmp_path, config=legacy, reference_assembly="hg19")
        assert "/legacy/region.fa" in captured_command[0]["command"]

    def test_an_incomplete_keyed_config_does_not_masquerade_as_legacy(self, tmp_path):
        """Only hg19 configured, plus a stray flat key: an hg38 run must not take it."""
        partial = {"shark_settings": {"muc1_region_fasta_hg19": "/a.fa", "muc1_region_fasta": "/flat.fa"}}
        with pytest.raises(ValueError, match="hg38"):
            filter_with(tmp_path, config=partial, reference_assembly="hg38")


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
