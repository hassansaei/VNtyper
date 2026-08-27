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
import shlex
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
            f"-m {db_file} --working_directory {output} -t 1"
        )

    def test_the_shipped_default_does_not_request_the_unused_alignment_sidecar(self, inputs, captured_command):
        """The optional ``-aln`` post-processing path can crash after writing a complete VCF."""
        db_file, sorted_bam, output = inputs

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        words = shlex.split(captured_command[0]["command"])
        assert "-aln" not in words
        assert "--aln" not in words

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


class TestShellInterpolationIsQuoted:
    """
    ``run_command`` executes with ``shell=True`` (trap 9), which makes quoting the caller's
    job. Every operand here originates outside ``config["tools"]`` -- the BAM path, the
    output directory, the sample-derived ``output_name`` -- and the web service accepts
    filenames that reach this line, so a space alone splits one argument into two and a
    ``;`` is arbitrary command execution.

    ``advntr_path`` and non-empty ``additional_commands`` stay unquoted on purpose: both hold command
    *fragments* from config (``"mamba run -n envadvntr advntr"``, or an explicit ``"-aln"``), which
    quoting would collapse into a single token.
    """

    #: Single path components, so each can be a real directory under ``tmp_path``.
    HOSTILE = ["my samples", "o'brien", "x; touch pwned", "$(whoami)", "`id`", 'quote"inside']

    @pytest.mark.parametrize("fragment", HOSTILE)
    def test_a_hostile_bam_path_survives_as_one_shell_word(self, tmp_path, captured_command, fragment):
        db_file = tmp_path / "vntr_db.db"
        db_file.write_text("db")
        directory = tmp_path / fragment
        directory.mkdir()
        sorted_bam = directory / "sample.sorted.bam"
        sorted_bam.write_text("bam")
        output = tmp_path / "advntr"
        output.mkdir()

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        words = shlex.split(captured_command[0]["command"])
        assert words[words.index("--alignment_file") + 1] == str(sorted_bam)

    @pytest.mark.parametrize("fragment", HOSTILE)
    def test_a_hostile_output_name_survives_as_one_shell_word(self, inputs, captured_command, fragment):
        db_file, sorted_bam, output = inputs

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), fragment, MAIN_CONFIG)

        words = shlex.split(captured_command[0]["command"])
        assert words[words.index("-o") + 1] == f"{output}/{fragment}_adVNTR.vcf"

    @pytest.mark.parametrize("fragment", HOSTILE)
    def test_a_hostile_output_directory_survives_in_both_places(self, tmp_path, captured_command, fragment):
        db_file = tmp_path / "vntr_db.db"
        db_file.write_text("db")
        sorted_bam = tmp_path / "sample.sorted.bam"
        sorted_bam.write_text("bam")
        output = tmp_path / fragment
        output.mkdir()

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        words = shlex.split(captured_command[0]["command"])
        assert words[words.index("--working_directory") + 1] == str(output)
        assert words[words.index("-o") + 1] == f"{output}/output_adVNTR.vcf"

    def test_the_tool_prefix_and_the_flag_list_stay_separate_words(self, inputs, captured_command, monkeypatch):
        """Quoting either would give bash one token where it needs several."""
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(advntr, "advntr_settings", {**advntr.advntr_settings, "additional_commands": "-aln"})

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        command = captured_command[0]["command"]
        assert command.startswith("mamba run -n envadvntr advntr genotype ")
        assert command.endswith(" -aln")


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
        monkeypatch.setattr(advntr, "advntr_settings", {"output_format": output_format, "threads": 1})

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

    def test_an_empty_settings_global_raises_rather_than_applying_defaults(self, inputs, captured_command, monkeypatch):
        """#247: was ``test_the_defaults_apply_when_the_settings_global_is_empty``.

        Those defaults contradicted the shipped configuration -- ``threads`` defaulted to 8
        while advntr_config.json sets 1, and ``output_format`` defaulted to "tsv" while it sets
        "vcf" -- so dropping either key changed the emitted command with no error at all.
        Configuration is authoritative now, as for the calibration constants in
        confidence_assignment.py:108-111, and a partial configuration fails loudly.
        """
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(advntr, "advntr_settings", {})

        with pytest.raises(KeyError, match="threads"):
            advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert captured_command == [], "no command may be emitted with an unconfigured thread count"

    def test_a_missing_output_format_raises_rather_than_flipping_the_extension(
        self, inputs, captured_command, monkeypatch
    ):
        """The second mismatch, which #247 does not name: the fallback was "tsv" while
        advntr_config.json ships "vcf"."""
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(advntr, "advntr_settings", {"threads": 1})

        with pytest.raises(KeyError, match="output_format"):
            advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert captured_command == []

    def test_a_caller_that_names_no_thread_count_gets_one(self, inputs, captured_command):
        """The ``1`` here is ``run_advntr``'s ``pipeline_threads`` default, not a pin in the
        shipped configuration: advntr_config.json sets ``threads: null``, which inherits.

        This assertion used to be called "the shipped configuration emits one thread" and
        was justified by adVNTR's ``-t`` setting only ``settings.CORES``, which nothing on
        the ``genotype -fs`` short-read path reads (#215). That was true until adVNTR 2.0.0
        moved the Viterbi DP into a ``nogil`` block and threaded the read loop. The
        assertion outlived its reason because the default is 1 at both ends, so it kept
        passing while saying something false about why. What matters now is in
        :class:`TestThreadsInheritThePipelineValue`.
        """
        db_file, sorted_bam, output = inputs

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert "-t 1" in captured_command[0]["command"]


class TestOutputExtensionIsDerivedInOnePlace:
    """#247: pipeline.py repeated advntr_genotyping.py's ``.get("output_format", "tsv")``
    fallback, so the producer of the adVNTR output path and its consumer could drift apart."""

    @pytest.mark.parametrize(
        ("output_format", "expected"),
        [("vcf", ".vcf"), ("tsv", ".tsv"), ("anything-else", ".tsv")],
    )
    def test_the_extension_follows_the_mapping_it_is_given(self, output_format, expected):
        assert advntr.advntr_output_extension({"output_format": output_format}) == expected

    def test_a_missing_key_raises_rather_than_defaulting(self):
        with pytest.raises(KeyError, match="output_format"):
            advntr.advntr_output_extension({})

    def test_it_reads_its_argument_and_not_the_module_global(self, monkeypatch):
        """The mapping is a parameter on purpose. run_advntr reads the import-time
        advntr_settings global while pipeline.py calls load_advntr_config() again and derives
        its own -- two independently loaded states. A helper that read the global would let
        the two disagree while appearing to share one source of truth.
        """
        monkeypatch.setattr(advntr, "advntr_settings", {"output_format": "tsv"})

        assert advntr.advntr_output_extension({"output_format": "vcf"}) == ".vcf"


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


class TestConfigLoading:
    def test_the_default_config_is_the_one_shipped_beside_the_module(self):
        loaded = advntr.load_advntr_config()

        assert set(loaded) == {"advntr_settings", "flagging_rules"}
        assert loaded["advntr_settings"]["vid"] == 25561

    def test_an_explicit_path_is_honoured(self, tmp_path):
        path = tmp_path / "custom.json"
        path.write_text('{"advntr_settings": {"vid": 1}}')

        assert advntr.load_advntr_config(str(path)) == {"advntr_settings": {"vid": 1}}


class TestThreadsInheritThePipelineValue:
    """``"threads": null`` means inherit ``--threads``; an integer overrides it.

    Before adVNTR 2.0.0 this would have been a false affordance: ``-t`` set only
    ``settings.CORES``, which nothing on the ``genotype -fs`` path read (#215, #254).
    2.0.0 moved the Viterbi DP into a ``nogil`` block and threaded the read loop, so the
    flag now has an effect and the pipeline's value is worth passing through.
    """

    def test_null_threads_inherit_the_pipeline_value(self, inputs, captured_command, monkeypatch):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"threads": None, "output_format": "vcf", "vid": 25561, "additional_commands": "-aln"},
        )

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG, pipeline_threads=12)

        assert captured_command[0]["command"].endswith("-t 12 -aln")

    def test_an_explicit_thread_count_overrides_the_pipeline_value(self, inputs, captured_command, monkeypatch):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"threads": 3, "output_format": "vcf", "vid": 25561, "additional_commands": "-aln"},
        )

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG, pipeline_threads=12)

        assert captured_command[0]["command"].endswith("-t 3 -aln")

    def test_a_missing_threads_key_still_raises(self, inputs, monkeypatch):
        """``null`` is not the same as absent; a partial mapping is unsupported input."""
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"output_format": "vcf", "vid": 25561, "additional_commands": "-aln"},
        )

        with pytest.raises(KeyError):
            advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG, pipeline_threads=4)

    def test_the_shipped_config_uses_null_so_the_cli_wins(self):
        """If this becomes an integer again, ``--threads`` silently stops reaching adVNTR."""
        config = advntr.load_advntr_config()

        assert config["advntr_settings"]["threads"] is None

    def test_the_default_is_one_so_callers_that_do_not_pass_it_are_unchanged(
        self, inputs, captured_command, monkeypatch
    ):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"threads": None, "output_format": "vcf", "vid": 25561, "additional_commands": "-aln"},
        )

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert captured_command[0]["command"].endswith("-t 1 -aln")


class TestThreadCountValidation:
    """Bad thread values fail here, naming the key and the file.

    adVNTR's own check (`advntr_commands.py:72-73`) says only "threads cannot be less
    than 1", which names neither where the value came from nor how to fix it. Validating
    before the command is built turns that into an actionable message.
    """

    def test_a_missing_key_names_the_key_and_the_file(self):
        with pytest.raises(KeyError) as excinfo:
            advntr.resolve_advntr_threads({}, 4)

        message = str(excinfo.value)
        assert "threads" in message
        assert "advntr_config.json" in message
        assert "null" in message

    def test_zero_is_rejected(self):
        with pytest.raises(ValueError, match="at least 1"):
            advntr.resolve_advntr_threads({"threads": 0}, 4)

    def test_a_negative_count_is_rejected(self):
        with pytest.raises(ValueError, match="at least 1"):
            advntr.resolve_advntr_threads({"threads": -8}, 4)

    def test_an_inherited_value_is_validated_too(self):
        """`null` delegates to --threads, so a bad pipeline value must not slip through."""
        with pytest.raises(ValueError, match="--threads"):
            advntr.resolve_advntr_threads({"threads": None}, 0)

    def test_a_non_integer_is_rejected_with_its_type(self):
        with pytest.raises(ValueError, match="str"):
            advntr.resolve_advntr_threads({"threads": "4"}, 1)

    def test_a_float_is_rejected(self):
        with pytest.raises(ValueError, match="float"):
            advntr.resolve_advntr_threads({"threads": 2.5}, 1)

    def test_true_is_rejected_rather_than_becoming_one_thread(self):
        """bool is an int subclass; True would otherwise silently mean one thread."""
        with pytest.raises(ValueError, match="bool"):
            advntr.resolve_advntr_threads({"threads": True}, 4)

    def test_an_explicit_count_is_returned(self):
        assert advntr.resolve_advntr_threads({"threads": 6}, 4) == 6

    def test_null_inherits_the_pipeline_value(self):
        assert advntr.resolve_advntr_threads({"threads": None}, 12) == 12


class TestAdditionalCommandsCannotOverrideAManagedOption:
    """``additional_commands`` is interpolated verbatim *after* every option ``run_advntr``
    sets, so any managed option repeated there wins.

    adVNTR declares ``-t``/``--threads`` as one argparse option
    (``advntr/__main__.py:88``) and argparse lets the last occurrence win. Measured against
    a parser built from adVNTR's own declarations, every one of these forms overrides an
    earlier ``-t 12``: ``--threads 3``, ``-t3``, ``--threads=3`` and the abbreviations
    ``--thr``/``--thread`` that ``allow_abbrev`` accepts by default.

    That made ``resolve_advntr_threads`` decorative: it rejects 0, negatives, bools and
    non-integers, and then a config that also says ``additional_commands: "--threads 0"``
    hands adVNTR the very value it refused. The same hole redirects ``-o`` and ``-m``,
    which would write the genotype somewhere ``pipeline.py`` does not read it back from.
    """

    @pytest.mark.parametrize(
        "additional",
        [
            "--threads 3",
            "-t 3",
            "-t3",
            "--threads=3",
            "--thr 3",
            "--thread 3",
        ],
        ids=["long", "short", "short-attached", "long-equals", "abbrev-short", "abbrev-long"],
    )
    def test_every_spelling_of_a_repeated_thread_count_is_refused(
        self, additional, inputs, captured_command, monkeypatch
    ):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"threads": 12, "output_format": "vcf", "additional_commands": additional},
        )

        with pytest.raises(ValueError, match="additional_commands"):
            advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert captured_command == [], "no command may be emitted once a managed option is duplicated"

    @pytest.mark.parametrize(
        "additional",
        ["-o /tmp/elsewhere.vcf", "--outfile /tmp/elsewhere.vcf", "-m /tmp/other.db", "--models /tmp/other.db"],
    )
    def test_redirecting_the_artefact_or_the_model_is_refused(self, additional, inputs, captured_command, monkeypatch):
        """``pipeline.py`` reconstructs the output path independently and reads the file
        back, so an ``-o`` nobody chose is a missing result, not a relocated one."""
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"threads": 1, "output_format": "vcf", "additional_commands": additional},
        )

        with pytest.raises(ValueError, match="additional_commands"):
            advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert captured_command == []

    def test_the_message_names_the_option_and_the_key_that_owns_it(self):
        with pytest.raises(ValueError) as excinfo:
            advntr.resolve_additional_commands({"additional_commands": "--threads 3"})

        message = str(excinfo.value)
        assert "--threads" in message
        assert "'threads'" in message, "the operator needs to be told which key to set instead"

    def test_a_missing_additional_commands_key_adds_no_optional_flags(self):
        """Packaged and missing-key defaults must both avoid adVNTR's optional post-processing."""
        assert advntr.resolve_additional_commands({}) == ""

    @pytest.mark.parametrize("additional", ["-aln", "-aln --haploid", "--fullru", "-u", ""])
    def test_flags_advntr_owns_alone_still_pass_through(self, additional):
        """The guard must not turn ``additional_commands`` into a dead extension point:
        explicit ``-aln`` and ``--haploid``/``--fullru``/``-u`` collide with nothing.
        """
        assert advntr.resolve_additional_commands({"additional_commands": additional}) == additional

    def test_an_explicit_alignment_sidecar_opt_in_reaches_the_command_unchanged(
        self, inputs, captured_command, monkeypatch
    ):
        db_file, sorted_bam, output = inputs
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {**advntr.advntr_settings, "additional_commands": "-aln"},
        )

        advntr.run_advntr(str(db_file), str(sorted_bam), str(output), "output", MAIN_CONFIG)

        assert captured_command[0]["command"].endswith("-aln")

    def test_an_unparseable_value_is_refused_rather_than_interpolated(self):
        """``run_command`` runs the result under ``shell=True``; an unbalanced quote would
        otherwise reach bash as a syntax error several stages later."""
        with pytest.raises(ValueError, match="additional_commands"):
            advntr.resolve_additional_commands({"additional_commands": "-aln 'unclosed"})


class TestTheGuardIsAnAllowListAndNotAPatternMatch:
    """A deny-list has to reimplement argparse's option matching, and the first one here
    did not. Two whole classes of spelling walked past it, both measured against a parser
    built from adVNTR's own declarations and appended to the exact command run_advntr
    emits:

    * **single-dash flag groups.** `consume_optional` decomposes them: a zero-argument
      short option consumes its first character and re-forms `'-' + tail`, so every one of
      `-p`, `-n`, `-e`, `-u` is a free prefix. `-pt3` set threads to 3. `-po /elsewhere.vcf`
      redirected the artefact -- the case the module docstring calls out as the worst,
      because pipeline.py rebuilds that path independently and reads it back.
    * **single-dash abbreviation.** `allow_abbrev` applies to single-dash options too, and
      no `-v` is declared, so `-v 99999` and `-vi 99999` both reach `-vid` and genotype a
      different VNTR.

    The fix inverts the rule: a word starting with `-` must be an option adVNTR declares
    AND one this module does not set. Anything else is refused, including abbreviations,
    attached values and flag groups, because none of them are exact declared spellings.
    That is fail-closed by construction rather than by enumeration, and it costs the same
    transcription the module was already paying.
    """

    @pytest.mark.parametrize(
        "additional",
        ["-pt3", "-nt3", "-et3", "-ut3", "-pt 3", "-po /elsewhere.vcf", "-pm /other.db", "-pa /other.bam"],
    )
    def test_a_flag_group_cannot_smuggle_in_a_managed_option(self, additional):
        with pytest.raises(ValueError, match="additional_commands"):
            advntr.resolve_additional_commands({"additional_commands": additional})

    @pytest.mark.parametrize("additional", ["-v 99999", "-vi 99999"])
    def test_a_single_dash_abbreviation_of_vid_is_refused(self, additional):
        """`-v 25561` is a spelling an operator could plausibly type."""
        with pytest.raises(ValueError, match="additional_commands"):
            advntr.resolve_additional_commands({"additional_commands": additional})

    @pytest.mark.parametrize(
        "additional",
        ["--thr 3", "--thread 3", "--work /tmp", "--model /db", "--out /x.vcf", "--vntr 1"],
    )
    def test_a_long_abbreviation_is_refused(self, additional):
        with pytest.raises(ValueError, match="additional_commands"):
            advntr.resolve_additional_commands({"additional_commands": additional})

    @pytest.mark.parametrize(
        "additional",
        [
            "-aln",
            "--aln",
            "-naive",
            "--naive",
            "-e",
            "-c 30",
            "-u",
            "-p",
            "-n",
            "--haploid",
            "--fullru",
            "--append",
            "--noref_aln",
            "--outfmt vcf",
            "--outfmt=vcf",
            "--min_read_length 100",
            "--vid_file /ids.txt",
            "-aln --haploid --fullru",
            "",
        ],
    )
    def test_every_flag_advntr_owns_alone_still_passes(self, additional):
        """The guard must refuse overrides, not turn additional_commands into a dead
        extension point. These are all declared `genotype` options run_advntr never sets."""
        assert advntr.resolve_additional_commands({"additional_commands": additional}) == additional

    def test_an_option_advntr_does_not_declare_is_refused_as_unknown(self):
        """A typo reaching adVNTR is an opaque argparse error from a subprocess."""
        with pytest.raises(ValueError, match="does not declare"):
            advntr.resolve_additional_commands({"additional_commands": "--haploidd"})

    def test_the_two_refusals_say_different_things(self):
        """A managed option needs "set X instead"; an unknown one needs "no such option"."""
        with pytest.raises(ValueError, match="advntr_settings\\['threads'\\]"):
            advntr.resolve_additional_commands({"additional_commands": "--threads 3"})

    def test_a_negative_number_is_a_value_and_not_an_option(self):
        assert advntr.resolve_additional_commands({"additional_commands": "-c -1.5"}) == "-c -1.5"

    @pytest.mark.parametrize("additional", ["--", "-"])
    def test_the_bare_separators_are_refused(self, additional):
        """Neither is a declared option, and `--` changes how everything after it parses."""
        with pytest.raises(ValueError, match="additional_commands"):
            advntr.resolve_additional_commands({"additional_commands": additional})

    def test_every_managed_option_names_an_owner(self):
        """The refusal path indexes the owners map; a managed option added without one
        turns a clear refusal into a KeyError."""
        assert set(advntr.MANAGED_ADVNTR_OPTIONS) == set(advntr.MANAGED_ADVNTR_OPTION_OWNERS)
