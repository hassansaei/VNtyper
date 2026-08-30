# tests/unit/test_utils.py

"""
Unit tests for utility functions.
Includes testing for command execution and BAM file validation.
"""

import gzip
import json
import logging
import os
import re
from unittest.mock import MagicMock, mock_open, patch

import pandas as pd
import pytest

from vntyper.scripts import utils
from vntyper.scripts.utils import (
    create_output_directories,
    get_tool_version,
    get_tool_versions,
    load_config,
    run_command,
    search,
    setup_logging,
    validate_bam_file,
    validate_fastq_file,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit

UTILS_LOGGER = "vntyper.scripts.utils"


def _completed_process(stdout="", stderr=""):
    """Stand-in for what ``subprocess.run`` returns: only ``.stdout``/``.stderr`` are read."""
    result = MagicMock()
    result.stdout = stdout
    result.stderr = stderr
    return result


def _logged(records, text, *, logger=UTILS_LOGGER, level=logging.ERROR):
    """
    True if some captured record at ``level`` from ``logger`` contains ``text``.

    ``caplog.set_level(level, logger=logger)`` only raises that logger's effective
    level so its records aren't filtered out -- the caplog handler itself is attached
    at the root logger and captures everything that propagates there, from any logger,
    at or above the *root's* level. A bare ``any(text in r.message for r in
    caplog.records)`` therefore passes just as well for a message emitted by a
    different logger (or at a different level) as for the one under test. Pinning both
    ``r.name`` and ``r.levelno`` here closes that gap.
    """
    return any(text in r.message and r.name == logger and r.levelno == level for r in records)


def test_run_command_success(tmp_path):
    """
    Test successful execution of a shell command.
    """
    log_file = tmp_path / "cmd.log"
    cmd = "echo 'Hello test_run_command'"
    result = run_command(cmd, str(log_file))
    assert result is True, "Expected run_command to succeed."
    assert "Hello test_run_command" in log_file.read_text()


@patch("subprocess.Popen")
def test_run_command_failure(mock_popen, tmp_path):
    """
    Test failure scenario for a shell command execution.
    """
    log_file = tmp_path / "fail.log"
    process_mock = MagicMock()
    process_mock.stdout = [b"Simulated error\n"]
    process_mock.wait.return_value = 1
    process_mock.returncode = 1
    mock_popen.return_value = process_mock

    cmd = "bad_command"
    ret = run_command(cmd, str(log_file))
    assert not ret, "Expected run_command to fail."
    assert "Simulated error" in log_file.read_text()


def test_run_command_survives_non_utf8_tool_output(tmp_path):
    """Corrupt tool bytes are visibly replaced instead of aborting the stage."""
    log_file = tmp_path / "cmd.log"

    assert run_command(r"printf '\xff\n'", str(log_file)) is True
    assert "�" in log_file.read_text(encoding="utf-8")


def test_run_command_still_reports_failure_after_non_utf8_output(tmp_path):
    """Corrupt output does not mask the child's non-critical failure status."""
    log_file = tmp_path / "cmd.log"

    assert run_command(r"printf '\xff\n'; exit 3", str(log_file)) is False


def test_run_command_reaps_child_when_stdout_iteration_fails(tmp_path):
    """A broken stdout iterator still closes the pipe and reaps the child."""
    process = MagicMock()
    process.stdout.__iter__.side_effect = OSError("pipe failure")
    with (
        patch.object(utils.subprocess, "Popen", return_value=process),
        pytest.raises(OSError, match="pipe failure"),
    ):
        run_command("true", str(tmp_path / "cmd.log"))

    process.stdout.close.assert_called_once()
    process.wait.assert_called_once()


def test_run_command_reaps_child_when_log_writing_fails(tmp_path):
    """A full command log does not leave its child running or zombified."""
    process = MagicMock()
    process.stdout = MagicMock()
    process.stdout.__iter__.return_value = iter([b"tool output\n"])
    log_file = MagicMock()
    log_file.write.side_effect = OSError("disk full")
    open_file = MagicMock()
    open_file.__enter__.return_value = log_file

    with (
        patch.object(utils.subprocess, "Popen", return_value=process),
        patch("builtins.open", return_value=open_file),
        pytest.raises(OSError, match="disk full"),
    ):
        run_command("true", str(tmp_path / "cmd.log"))

    process.stdout.close.assert_called_once()
    process.wait.assert_called_once()


def test_run_command_waits_when_stdout_close_fails(tmp_path):
    """A pipe-close error cannot prevent reaping the completed child."""
    process = MagicMock()
    process.stdout.__iter__.return_value = iter([])
    process.stdout.close.side_effect = OSError("close failure")
    with (
        patch.object(utils.subprocess, "Popen", return_value=process),
        pytest.raises(OSError, match="close failure"),
    ):
        run_command("true", str(tmp_path / "cmd.log"))

    process.wait.assert_called_once()


def test_validate_bam_file_success(tmp_path, test_config):
    """
    Test validation of a valid BAM file.
    We'll create a local file and mock out run_command so samtools quickcheck passes.
    """
    bam_file = tmp_path / "temp_valid.bam"
    bam_file.touch()

    # Mock out run_command to simulate a passing samtools quickcheck.
    with patch("vntyper.scripts.utils.run_command", return_value=True):
        validate_bam_file(str(bam_file))  # Should not raise ValueError.


def test_validate_bam_file_passes_the_configured_samtools_to_quickcheck(tmp_path):
    """F1: quickcheck honors the configured samtools executable."""
    bam_file = tmp_path / "sample.bam"
    bam_file.write_text("x", encoding="utf-8")

    with patch.object(utils, "run_command", return_value=True) as run:
        validate_bam_file(str(bam_file), samtools_path="/opt/tools/samtools")

    assert run.call_args.args[0] == f"/opt/tools/samtools quickcheck -v {bam_file}"


def test_validate_bam_file_uses_configured_samtools_when_quickcheck_fails(tmp_path):
    """F1: configured samtools reaches quickcheck before its failure is reported."""
    bam_file = tmp_path / "broken.bam"
    bam_file.write_text("x", encoding="utf-8")

    with (
        patch.object(utils, "run_command", return_value=False) as run,
        pytest.raises(ValueError, match="failed quickcheck"),
    ):
        validate_bam_file(str(bam_file), samtools_path="/opt/tools/samtools")

    assert run.call_args.args[0] == f"/opt/tools/samtools quickcheck -v {bam_file}"


def test_validate_bam_file_nonexistent():
    """
    Test behavior when a nonexistent BAM file is validated.
    """
    with pytest.raises(ValueError) as exc:
        validate_bam_file("nonexistent.bam")
    assert "does not exist" in str(exc.value), "Expected ValueError for nonexistent file."


# ---------------------------------------------------------------------------
# setup_logging
# ---------------------------------------------------------------------------


@pytest.fixture
def clean_root_logger():
    """
    ``setup_logging`` mutates the root logger's handlers and level as a side effect.

    Snapshot and restore around each test so these assertions neither leak into later
    tests nor depend on whatever logging state pytest's own machinery (``caplog``,
    ``log_cli``) has already attached to the root logger.
    """
    root_logger = logging.getLogger()
    original_handlers = list(root_logger.handlers)
    original_level = root_logger.level
    try:
        yield root_logger
    finally:
        for handler in list(root_logger.handlers):
            if handler not in original_handlers:
                root_logger.removeHandler(handler)
                handler.close()
        for handler in original_handlers:
            if handler not in root_logger.handlers:
                root_logger.addHandler(handler)
        root_logger.setLevel(original_level)


def test_setup_logging_attaches_file_and_console_handlers(tmp_path, clean_root_logger):
    """utils.py:95-105 -- a log_file gets both a FileHandler and a console StreamHandler,
    both set to the requested level."""
    log_file = tmp_path / "run.log"

    setup_logging(log_level=logging.WARNING, log_file=str(log_file))

    handler_types = {type(h) for h in clean_root_logger.handlers}
    assert logging.FileHandler in handler_types
    assert logging.StreamHandler in handler_types
    assert clean_root_logger.level == logging.WARNING

    file_handlers = [h for h in clean_root_logger.handlers if type(h) is logging.FileHandler]
    assert file_handlers[0].level == logging.WARNING

    logging.getLogger("vntyper.scripts.utils_test_marker").warning("hello marker")
    for handler in clean_root_logger.handlers:
        handler.flush()
    assert "hello marker" in log_file.read_text()


def test_setup_logging_clears_previously_attached_handlers(clean_root_logger):
    """utils.py:89->90 -- hasHandlers()==True: old handlers are cleared, not accumulated."""
    clean_root_logger.addHandler(logging.NullHandler())
    assert clean_root_logger.hasHandlers()

    setup_logging(log_level=logging.INFO, log_file=None)

    assert len(clean_root_logger.handlers) == 1
    assert type(clean_root_logger.handlers[0]) is logging.StreamHandler


def test_setup_logging_with_no_prior_handlers_still_attaches_a_console_handler(clean_root_logger):
    """utils.py:89->92 -- hasHandlers()==False must not skip attaching the console handler,
    and log_file=None (95->102) must not attach a FileHandler."""
    clean_root_logger.handlers = []
    assert not clean_root_logger.hasHandlers()

    setup_logging(log_level=logging.ERROR, log_file=None)

    assert len(clean_root_logger.handlers) == 1
    assert type(clean_root_logger.handlers[0]) is logging.StreamHandler
    assert clean_root_logger.level == logging.ERROR


# ---------------------------------------------------------------------------
# create_output_directories
# ---------------------------------------------------------------------------


def test_create_output_directories_creates_the_five_pipeline_subdirectories(tmp_path):
    """utils.py:118-133 happy path: every directory named in DIRS exists after the call."""
    base = tmp_path / "out"

    dirs = create_output_directories(str(base))

    assert set(dirs) == {
        "base",
        "kestrel",
        "advntr",
        "fastq_bam_processing",
        "alignment_processing",
        "coverage",
    }
    for path in dirs.values():
        assert os.path.isdir(path)


def test_create_output_directories_logs_and_reraises_when_makedirs_fails(tmp_path, caplog):
    """utils.py:134-136 -- a mkdir failure is logged at ERROR and re-raised, not swallowed."""
    base = tmp_path / "out"
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with (
        patch("vntyper.scripts.utils.os.makedirs", side_effect=PermissionError("denied")),
        pytest.raises(PermissionError),
    ):
        create_output_directories(str(base))

    assert _logged(caplog.records, "Failed to create directory")


# ---------------------------------------------------------------------------
# get_tool_version
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "command, version_flag, stdout, expected",
    [
        ("fastp", "--version", "banner\nfastp 0.23.4\n", "0.23.4"),
        ("samtools", "--version", "Using htslib\nsamtools 1.17\n", "1.17"),
        ("bwa", "", "Program: bwa\nVersion: 0.7.17-r1188\n", "0.7.17-r1188"),
        ("/opt/bin/advntr", "", "l0\nl1\nadVNTR 1.4.0: info\n", "1.4.0"),
        ("java -jar kestrel.jar", "-h", "l0\nl1\nkestrel: 1.0.1\n", "1.0.1"),
        (
            "java",
            "--version",
            "openjdk 17.0.2 2022-01-18\nOpenJDK Runtime Environment\n",
            "openjdk 17.0.2 2022-01-18",
        ),
        ("tool-x", "", "whatever\n", "unknown"),
    ],
    ids=["fastp", "samtools", "bwa", "advntr", "kestrel", "java_path", "unrecognized"],
)
def test_get_tool_version_parses_the_expected_tool_format(command, version_flag, stdout, expected):
    """utils.py:160-186 -- each tool-specific parsing branch, happy path."""
    with patch("vntyper.scripts.utils.subprocess.run", return_value=_completed_process(stdout=stdout)):
        assert get_tool_version(command, version_flag) == expected


@pytest.mark.parametrize(
    "command, version_flag, stdout",
    [
        ("fastp", "--version", "nothing useful\n"),
        ("samtools", "--version", "nope\n"),
        ("bwa", "", "Program: bwa\nUsage: bwa command\n"),
        ("advntr", "", "only one line\n"),
        ("java -jar kestrel.jar", "-h", "no info here\n"),
    ],
    ids=["fastp", "samtools", "bwa", "advntr", "kestrel"],
)
def test_get_tool_version_returns_unknown_when_the_marker_is_absent(command, version_flag, stdout):
    """utils.py:161-183 -- the "marker not in output" arm of each tool branch."""
    with patch("vntyper.scripts.utils.subprocess.run", return_value=_completed_process(stdout=stdout)):
        assert get_tool_version(command, version_flag) == "unknown"


def test_get_tool_version_falls_back_to_stderr_when_stdout_is_blank():
    """utils.py:157 -- ``result.stdout.strip() or result.stderr.strip()``: some tools
    print their version banner to stderr rather than stdout."""
    with patch(
        "vntyper.scripts.utils.subprocess.run",
        return_value=_completed_process(stdout="   ", stderr="banner\nsamtools 1.18\n"),
    ):
        assert get_tool_version("samtools", "--version") == "1.18"


def test_get_tool_version_returns_unknown_and_logs_when_the_binary_is_missing(caplog):
    """utils.py:188-190 -- FileNotFoundError is swallowed into 'unknown', not raised."""
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)
    with patch("vntyper.scripts.utils.subprocess.run", side_effect=FileNotFoundError("no such file")):
        assert get_tool_version("fastp", "--version") == "unknown"
    assert _logged(caplog.records, "Command not found")


def test_get_tool_version_returns_unknown_and_logs_on_permission_denied(caplog):
    """utils.py:191-193 -- PermissionError is swallowed into 'unknown', not raised."""
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)
    with patch("vntyper.scripts.utils.subprocess.run", side_effect=PermissionError("denied")):
        assert get_tool_version("fastp", "--version") == "unknown"
    assert _logged(caplog.records, "Permission denied")


def test_get_tool_version_returns_unknown_and_logs_on_malformed_tool_output(caplog):
    """utils.py:194-196 -- a version line with no second token raises IndexError inside
    the parser; it must not escape get_tool_version."""
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)
    with patch(
        "vntyper.scripts.utils.subprocess.run",
        return_value=_completed_process(stdout="only\nfastp\n"),
    ):
        assert get_tool_version("fastp", "--version") == "unknown"
    assert _logged(caplog.records, "Failed to parse version")


def test_tool_version_unexpected_failure_returns_unknown(caplog):
    """utils.py:197-199 -- an unexpected subprocess failure remains observable."""
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)
    with patch("vntyper.scripts.utils.subprocess.run", side_effect=OSError("exec failed")):
        assert get_tool_version("fastp", "--version") == "unknown"
    assert _logged(caplog.records, "Failed to get version for fastp: exec failed")


# ---------------------------------------------------------------------------
# get_tool_versions
# ---------------------------------------------------------------------------


def test_get_tool_versions_maps_each_named_tool_to_its_parsed_version():
    """The config-driven fan-out to get_tool_version.

    The fan-out is now over the intersection of ``config["tools"]`` and the caller's
    declared set rather than over the whole config (#262), so this names all three
    tools to keep measuring exactly what it measured before: that each configured tool
    reaches ``get_tool_version`` with its own command and version flag.
    """
    config = {
        "tools": {
            "fastp": "/usr/bin/fastp",
            "samtools": "/usr/bin/samtools",
            "java_path": "/usr/bin/java",
        }
    }

    def fake_get_tool_version(command, version_flag):
        return f"{command}::{version_flag}"

    with patch("vntyper.scripts.utils.get_tool_version", side_effect=fake_get_tool_version) as mocked:
        versions = get_tool_versions(config, tools_in_use={"fastp", "samtools", "java_path"})

    assert versions == {
        "fastp": "/usr/bin/fastp::",
        "samtools": "/usr/bin/samtools::",
        "java_path": "/usr/bin/java::--version",
    }
    mocked.assert_any_call("/usr/bin/fastp", "")
    mocked.assert_any_call("/usr/bin/java", "--version")


def test_get_tool_versions_with_no_tools_configured_returns_an_empty_dict():
    """config.get("tools", {}) with nothing configured, whatever the caller declares."""
    assert get_tool_versions({}, tools_in_use={"samtools", "kestrel"}) == {}


#: Tail of the real ``java -jar kestrel.jar -h`` output: the last line is what
#: ``get_tool_version`` parses (``output.split("\n")[-1].split(": ")[1]``, utils.py:182).
KESTREL_HELP_OUTPUT = "Usage: kestrel [options]\n  -h[<TOPIC>]  Show help\nkestrel version: 1.0.1"

#: First line of ``java --version``, which is all ``get_tool_version`` keeps (utils.py:185).
JAVA_VERSION_OUTPUT = "openjdk 11.0.23-internal 2024-04-16\nOpenJDK Runtime Environment (build 11.0.23-internal)"


def test_get_tool_versions_passes_the_kestrel_help_flag_exactly_once():
    """
    SPECIFICATION: kestrel is invoked as ``java -jar <jar> -h``, once.

    get_tool_versions used to build ``command`` for kestrel as
    ``f"{java_path} {version_flag}"`` (utils.py:232), which already embeds the
    ``-jar "<path>" -h`` flag, and then call ``get_tool_version(command,
    version_flag)`` with that SAME ``version_flag`` again (utils.py:233).
    get_tool_version concatenates ``shlex.split(command) + shlex.split(version_flag)``
    (utils.py:155), so the argv reaching ``subprocess.run`` was
    ``java -jar k.jar -h -jar k.jar -h``. That happened to exit 0 against the pinned
    kestrel build -- its arg parser stops at the first ``-h`` -- but nothing in
    get_tool_versions depends on kestrel continuing to ignore trailing garbage: a
    stricter parser returns "unknown" (utils.py:183), raises IndexError into
    "unknown" (utils.py:194-196), or yields a usage banner that parses into a
    plausible-looking garbage version (``"Usage: kestrel [options]".split(": ")[1]``).

    Filed as issue-utils-kestrel-duplicate-version-flag.md and fixed on the
    repository owner's instruction: ``version_flag`` is cleared once it has been
    folded into ``command``. Asserting the exact argv, not just its length, is the
    point -- the operand list is the whole defect.
    """
    config = {
        "tools": {
            "kestrel": "/opt/kestrel/kestrel-1.0.1.jar",
            "java_path": "/usr/bin/java",
        }
    }

    with patch(
        "vntyper.scripts.utils.subprocess.run",
        return_value=_completed_process(),
    ) as mocked_run:
        get_tool_versions(config, tools_in_use={"kestrel", "java_path"})

    kestrel_calls = [call for call in mocked_run.call_args_list if "-jar" in call.args[0]]
    assert len(kestrel_calls) == 1, "expected exactly one subprocess.run call built for kestrel"

    assert kestrel_calls[0].args[0] == [
        "/usr/bin/java",
        "-jar",
        "/opt/kestrel/kestrel-1.0.1.jar",
        "-h",
    ]


def test_get_tool_versions_reports_the_kestrel_version_from_its_help_banner():
    """
    SPECIFICATION: removing the duplicated flag must not move the reported version.

    ``get_tool_version`` dispatches its kestrel parse on the *command* string --
    ``if "java" in command and "kestrel" in command`` (utils.py:179) -- so the fix
    had to keep the pre-built ``command`` and drop the second ``version_flag``. The
    other candidate fix (stop pre-building ``command``, pass ``java_path`` alone and
    let get_tool_version append the flag) breaks that dispatch: with the shipped
    config's ``java_path == "java"``, ``"kestrel" in command`` is then False, and the
    ``command.startswith("java")`` branch (utils.py:184-185) takes over and reports
    the *first* line of kestrel's help banner -- ``"Usage: kestrel [options]"`` --
    as the kestrel version. That was confirmed by inducing it here. The argv
    assertion above still passes under that mutation, so this test is the only thing
    holding the reported version in place.
    """
    config = {
        "tools": {
            "kestrel": "vntyper/dependencies/kestrel/kestrel.jar",
            "java_path": "java",
        }
    }

    def fake_run(argv, **kwargs):
        return _completed_process(stdout=KESTREL_HELP_OUTPUT if "-jar" in argv else JAVA_VERSION_OUTPUT)

    with patch("vntyper.scripts.utils.subprocess.run", side_effect=fake_run):
        versions = get_tool_versions(config, tools_in_use={"kestrel", "java_path"})

    assert versions == {
        "kestrel": "1.0.1",
        "java_path": "openjdk 11.0.23-internal 2024-04-16",
    }


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------


def test_search_matches_case_insensitively_by_default():
    """utils.py:250-257 -- the happy path, default case=False."""
    df = pd.DataFrame({"Gene": ["MUC1", "muc1", "TP53"], "Value": [1, 2, 3]})

    result = search("muc1", df)

    assert list(result["Gene"]) == ["MUC1", "muc1"]


def test_search_is_case_sensitive_when_requested():
    """utils.py:250-257 -- case=True narrows the match."""
    df = pd.DataFrame({"Gene": ["MUC1", "muc1", "TP53"], "Value": [1, 2, 3]})

    result = search("MUC1", df, case=True)

    assert list(result["Gene"]) == ["MUC1"]


def test_search_logs_and_reraises_on_an_invalid_regex(caplog):
    """utils.py:258-260 -- a regex compile error is logged at ERROR and re-raised."""
    df = pd.DataFrame({"Gene": ["MUC1", "muc1", "TP53"]})
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with pytest.raises(re.error):
        search("(unclosed", df)

    assert _logged(caplog.records, "Error during regex search")


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------


def test_load_config_reads_a_user_provided_path(tmp_path):
    """utils.py:273-279 -- the explicit-path happy path."""
    config_path = tmp_path / "cfg.json"
    config_path.write_text(json.dumps({"a": 1}))

    assert load_config(str(config_path)) == {"a": 1}


def test_load_config_raises_and_logs_on_malformed_json(tmp_path, caplog):
    """utils.py:280-282 -- a JSONDecodeError from a bad user config is logged and re-raised."""
    config_path = tmp_path / "cfg.json"
    config_path.write_text("{not valid json")
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with pytest.raises(json.JSONDecodeError):
        load_config(str(config_path))

    assert _logged(caplog.records, "Error decoding JSON")


def test_load_config_raises_and_logs_on_an_unreadable_user_config(tmp_path, caplog):
    """utils.py:283-285 -- any other read failure (here: the path is a directory) is
    logged and re-raised, not swallowed into a default."""
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with pytest.raises(OSError):
        load_config(str(tmp_path))

    assert _logged(caplog.records, "Unexpected error loading config file")


def test_load_config_falls_back_to_the_package_default_when_no_path_is_given():
    """utils.py:287-292 -- config_path=None takes the packaged-config branch."""
    with patch(
        "vntyper.scripts.utils.pkg_resources.open_text",
        mock_open(read_data=json.dumps({"default": True})),
    ) as mocked_open:
        assert load_config(None) == {"default": True}

    mocked_open.assert_called_once_with("vntyper", "config.json")


def test_load_config_falls_back_to_the_package_default_when_the_given_path_is_missing(tmp_path):
    """utils.py:273 -- os.path.exists()==False on a provided path also takes the
    default-config branch."""
    missing_path = tmp_path / "does-not-exist.json"

    with patch(
        "vntyper.scripts.utils.pkg_resources.open_text",
        mock_open(read_data=json.dumps({"default": True})),
    ):
        assert load_config(str(missing_path)) == {"default": True}


def test_load_config_exits_when_the_package_default_config_is_unreadable(caplog):
    """utils.py:293-296 -- a broken packaged config.json exits the process rather than
    returning a partial config the rest of the pipeline would silently misread."""
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with (
        patch("vntyper.scripts.utils.pkg_resources.open_text", side_effect=OSError("boom")),
        pytest.raises(SystemExit) as excinfo,
    ):
        load_config(None)

    assert excinfo.value.code == 1
    assert _logged(caplog.records, "Default config file not found")


# ---------------------------------------------------------------------------
# validate_bam_file
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("file_path", [None, ""])
def test_validate_bam_file_rejects_a_falsy_path(file_path):
    """utils.py:316-318 -- no path at all, before any filesystem check."""
    with pytest.raises(ValueError, match="No alignment file provided"):
        validate_bam_file(file_path)


def test_validate_bam_file_rejects_the_wrong_extension(tmp_path):
    """utils.py:325-327 -- an existing file with neither .bam nor .cram."""
    bad_file = tmp_path / "notes.txt"
    bad_file.touch()

    with pytest.raises(ValueError, match="Invalid alignment file extension"):
        validate_bam_file(str(bad_file))


def test_validate_bam_file_accepts_the_cram_extension(tmp_path):
    """utils.py:325 -- .cram, not just .bam, is a valid extension: this module's
    docstring says the function was extended to cover CRAM inputs too."""
    cram_file = tmp_path / "sample.cram"
    cram_file.touch()

    with patch("vntyper.scripts.utils.run_command", return_value=True):
        validate_bam_file(str(cram_file))  # must not raise


@pytest.mark.parametrize(
    "filename",
    ["SAMPLE.CRAM", "SAMPLE.BAM", "sample.Cram", "sample.Bam", "NA12878.CRAM"],
    ids=["upper-cram", "upper-bam", "mixed-cram", "mixed-bam", "upper-cram-realistic"],
)
def test_validate_bam_file_accepts_an_uppercase_alignment_extension(tmp_path, filename):
    """An upper-case extension is an accepted upload name, so it must survive validation.

    Three layers see the same file and they have to agree on what a CRAM is:

    * ``docker/app/uploads.py`` compiles its allowlist with ``re.IGNORECASE``
      (and says so, deliberately: "sequencers and LIMS exports routinely
      upper-case it"), so ``SAMPLE.CRAM`` is a name the endpoint stores;
    * ``docker/app/tasks.py`` lower-cases the suffix before choosing ``--cram``
      or ``--bam``, so the worker routes it correctly, and
      ``tests/unit/web/test_cram_alignment_handoff.py::test_the_flag_is_chosen_case_insensitively``
      pins that;
    * this function then decided the extension with a case-*sensitive*
      ``endswith((".bam", ".cram"))`` and raised ``ValueError``.

    So an upload the endpoint accepted and the worker routed died in the
    validator, and the web test could not see it: it inspects the argument
    vector the worker builds and never runs that command through
    ``validate_bam_file``. This test closes exactly that gap.

    Nothing downstream needs the lower-case name. ``samtools`` identifies BAM
    and CRAM from the file's own signature, not from its extension --
    ``samtools quickcheck -v`` on files named ``SAMPLE.BAM`` and
    ``SAMPLE.CRAM`` was measured returning 0 -- and ``pipeline.py`` branches on
    which CLI flag carried the path, never on the path's suffix.

    Args:
        tmp_path: Pytest temporary directory.
        filename: An alignment filename whose extension is not all lower-case.
    """
    alignment_file = tmp_path / filename
    alignment_file.touch()

    with patch("vntyper.scripts.utils.run_command", return_value=True):
        validate_bam_file(str(alignment_file))  # must not raise


@pytest.mark.parametrize(
    "filename",
    ["notes.TXT", "reads.SAM", "sample.CRAM.gz", "sample.BAMBOO"],
    ids=["txt", "sam", "gzipped-cram", "longer-suffix"],
)
def test_validate_bam_file_still_rejects_a_wrong_extension_whatever_its_case(tmp_path, filename):
    """Case-folding the suffix must widen the check to casing, and to nothing else.

    ``.sam`` is the one that matters: an uncompressed SAM is not something
    ``pipeline.py``'s BAM path can slice, and lower-casing before comparing
    must not be mistaken for "accept anything alignment-shaped".

    Args:
        tmp_path: Pytest temporary directory.
        filename: A filename whose extension is neither ``.bam`` nor ``.cram``.
    """
    bad_file = tmp_path / filename
    bad_file.touch()

    with pytest.raises(ValueError, match="Invalid alignment file extension"):
        validate_bam_file(str(bad_file))


def test_validate_bam_file_passes_cwd_through_to_run_command(tmp_path):
    """
    utils.py:334 -- ``cwd`` is threaded into the quickcheck ``run_command`` call
    unchanged, not dropped or replaced.

    AGENTS.md trap 7: every tool and reference path in ``config.json`` is relative to
    the process working directory, and ``pipeline.py`` pins ``project_root =
    os.getcwd()`` and threads it into this call (and the Java/samtools calls
    elsewhere) specifically so quickcheck resolves paths against the right directory.
    Dropping ``cwd=cwd`` here leaves every other assertion in this file passing while
    quickcheck silently resolves against whatever directory the process happened to
    inherit.
    """
    bam_file = tmp_path / "sample.bam"
    bam_file.touch()

    with patch("vntyper.scripts.utils.run_command", return_value=True) as mocked_run_command:
        validate_bam_file(str(bam_file), cwd="/somewhere")

    _, kwargs = mocked_run_command.call_args
    assert kwargs["cwd"] == "/somewhere"


def test_validate_bam_file_raises_when_run_command_reports_quickcheck_failure(tmp_path):
    """utils.py:336-337 -- run_command()==False is turned into a ValueError naming the file."""
    bam_file = tmp_path / "broken.bam"
    bam_file.touch()

    with (
        patch("vntyper.scripts.utils.run_command", return_value=False),
        pytest.raises(ValueError, match="failed quickcheck"),
    ):
        validate_bam_file(str(bam_file))


def test_validate_bam_file_raises_valueerror_when_the_real_quickcheck_fails(tmp_path, caplog):
    """
    SPECIFICATION: a failed quickcheck raises the documented ValueError, naming the file.

    This is the end-to-end counterpart of the test above: nothing is stubbed between
    validate_bam_file and ``subprocess.Popen``, so it exercises the real
    ``run_command`` and pins the type a caller actually sees.

    validate_bam_file used to call ``run_command(..., critical=True, ...)``, and
    run_command's critical path *raises* ``RuntimeError`` rather than returning False
    (utils.py:66-73, pinned by
    ``test_run_command_contract.py::test_a_critical_failure_raises_runtime_error``).
    The ``if not success:`` branch was therefore unreachable by construction and the
    documented ``Raises: ValueError`` was never true on a real failure. Passing
    ``critical=False`` makes the documented contract real: validate_bam_file aborts
    the run itself, by raising, so run_command does not also need to.

    ValueError, not RuntimeError, is the type this repository uses for "the input the
    user handed us is not usable" -- it is what validate_bam_file's own three earlier
    checks raise, what its sibling validate_fastq_file raises for every failure, what
    ``cli.py:110`` catches to turn a usage error into a clean exit, and what
    ``test_pipeline_guards.py`` already models a validate_bam_file failure as.
    ``RuntimeError`` is this repository's "an execution we depended on failed"
    (``cli_handlers.py:406``); a truncated BAM is a verdict about the input, not an
    infrastructure failure. Filed as
    issue-utils-validate-bam-file-wrong-exception-type.md and fixed on the repository
    owner's instruction, in the same fail-loud direction as issue #185.
    """
    bam_file = tmp_path / "broken.bam"
    bam_file.touch()
    failing_process = MagicMock()
    failing_process.stdout = [b"[E::hts_open] truncated file\n"]
    failing_process.wait.return_value = 1
    failing_process.returncode = 1
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with (
        patch("vntyper.scripts.utils.subprocess.Popen", return_value=failing_process),
        pytest.raises(ValueError) as excinfo,
    ):
        validate_bam_file(str(bam_file))

    assert str(excinfo.value) == f"Alignment file failed quickcheck: {bam_file}"
    assert _logged(caplog.records, f"Alignment file failed quickcheck: {bam_file}")


def test_validate_bam_file_passes_the_log_dir_through_to_run_command(tmp_path):
    """#201: the log is an artefact of the run and belongs with the run's output.

    ``run_command`` opens its log file before it spawns anything, so deriving
    that path from the input alignment made a read-only input mount an unhandled
    ``PermissionError`` raised before quickcheck ever ran.
    """
    bam_file = tmp_path / "in" / "sample.bam"
    bam_file.parent.mkdir()
    bam_file.write_text("dummy")
    out = tmp_path / "out"
    out.mkdir()

    with patch("vntyper.scripts.utils.run_command", return_value=True) as mocked_run_command:
        validate_bam_file(str(bam_file), log_dir=str(out))

    log_file = mocked_run_command.call_args[0][1]
    assert log_file == str(out / "sample.bam.quickcheck.log")
    assert str(bam_file.parent) not in log_file


def test_validate_bam_file_without_a_log_dir_still_logs_beside_the_input(tmp_path):
    """``log_dir=None`` keeps today's behaviour, deliberately.

    That default is the issue author's explicit recommendation -- it keeps #201 a
    contained change rather than a breaking one -- so it is pinned rather than
    left to drift. ``pipeline.py`` is what opts in, at both of its call sites.
    """
    bam_file = tmp_path / "in" / "sample.bam"
    bam_file.parent.mkdir()
    bam_file.write_text("dummy")

    with patch("vntyper.scripts.utils.run_command", return_value=True) as mocked_run_command:
        validate_bam_file(str(bam_file))

    assert mocked_run_command.call_args[0][1] == f"{bam_file}.quickcheck.log"


def test_validate_bam_file_names_the_log_after_the_alignment_not_its_path(tmp_path):
    """Two runs of two same-named alignments must not overwrite one log.

    The log is named for the alignment's *basename* inside ``log_dir``, so
    ``a/sample.bam`` and ``b/sample.bam`` validated into different output
    directories keep separate logs -- and neither reaches back into its own
    input tree.
    """
    out = tmp_path / "out"
    out.mkdir()
    logs = []
    for source in ("a", "b"):
        bam_file = tmp_path / source / "sample.bam"
        bam_file.parent.mkdir()
        bam_file.write_text("dummy")
        with patch("vntyper.scripts.utils.run_command", return_value=True) as mocked_run_command:
            validate_bam_file(str(bam_file), log_dir=str(out))
        logs.append(mocked_run_command.call_args[0][1])

    assert logs == [str(out / "sample.bam.quickcheck.log")] * 2


# ---------------------------------------------------------------------------
# validate_fastq_file
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("file_path", [None, ""])
def test_validate_fastq_file_rejects_a_falsy_path(file_path):
    """utils.py:352-354 -- no path at all, before any filesystem check."""
    with pytest.raises(ValueError, match="No FASTQ file provided"):
        validate_fastq_file(file_path)


def test_validate_fastq_file_rejects_a_nonexistent_file():
    """utils.py:356-358."""
    with pytest.raises(ValueError, match="does not exist"):
        validate_fastq_file("nonexistent.fastq")


def test_validate_fastq_file_rejects_the_wrong_extension(tmp_path):
    """utils.py:360-363."""
    bad_file = tmp_path / "reads.txt"
    bad_file.touch()

    with pytest.raises(ValueError, match="Invalid FASTQ file extension"):
        validate_fastq_file(str(bad_file))


def test_validate_fastq_file_accepts_a_well_formed_plain_fastq(tmp_path):
    """utils.py:365-374 -- happy path, uncompressed."""
    fastq_file = tmp_path / "reads.fastq"
    fastq_file.write_text("@SEQ_ID\nGATTACA\n+\nIIIIIII\n")

    validate_fastq_file(str(fastq_file))  # must not raise


def test_validate_fastq_file_accepts_a_gzipped_fastq(tmp_path):
    """utils.py:367 -- the gzip.open branch of open_func."""
    fastq_file = tmp_path / "reads.fastq.gz"
    with gzip.open(fastq_file, "wt") as handle:
        handle.write("@SEQ_ID\nGATTACA\n+\nIIIIIII\n")

    validate_fastq_file(str(fastq_file))  # must not raise


def test_validate_fastq_file_rejects_an_empty_file(tmp_path):
    """utils.py:369-373 -- the very first readline() already returns ""."""
    fastq_file = tmp_path / "empty.fastq"
    fastq_file.touch()

    with pytest.raises(ValueError, match="incomplete or empty"):
        validate_fastq_file(str(fastq_file))


def test_validate_fastq_file_rejects_a_truncated_read(tmp_path):
    """utils.py:369-373 -- fewer than 4 lines for the first read."""
    fastq_file = tmp_path / "truncated.fastq"
    fastq_file.write_text("@SEQ_ID\nGATTACA\n")

    with pytest.raises(ValueError, match="incomplete or empty"):
        validate_fastq_file(str(fastq_file))


@pytest.mark.parametrize(
    "filename",
    ["reads.FASTQ.GZ", "reads.FASTQ", "reads.Fq.Gz", "reads.FQ"],
    ids=["upper-gz", "upper-plain", "mixed-gz", "upper-fq"],
)
def test_validate_fastq_file_extension_check_stays_case_sensitive(tmp_path, filename):
    """The FASTQ validator does NOT case-fold, and that asymmetry is on purpose.

    ``validate_bam_file`` was made case-insensitive because the upload
    allowlist above it already accepts ``SAMPLE.CRAM``; the obvious follow-on
    is to do the same here. It would be wrong, and this test is what stops it.

    ``fastp`` decides whether to decompress from a case-sensitive ``.gz``
    suffix. Measured with the pinned fastp v0.23.4, on the same 5-read gzipped
    FASTQ written under two names::

        reads.fastq.gz -> exit=0  reads_written=5
        reads.FASTQ.GZ -> exit=0  reads_written=0

    So accepting ``reads.FASTQ.GZ`` here would trade a clear ``ValueError`` at
    the boundary for a silent zero-read run that reaches Kestrel and produces a
    genotype from nothing. There is also no layer above this one that disagrees:
    the web service has no FASTQ endpoint, so unlike the alignment case there is
    no accepted-then-rejected contradiction to resolve.

    This is a decision taken in this fix lane with the measurement attached, not
    a behaviour signed off by the repository owner; it is a candidate for review
    if FASTQ upload is ever added or fastp's detection changes.

    Args:
        tmp_path: Pytest temporary directory.
        filename: A FASTQ filename whose extension is not all lower-case.
    """
    fastq_file = tmp_path / filename
    fastq_file.touch()

    with pytest.raises(ValueError, match="Invalid FASTQ file extension"):
        validate_fastq_file(str(fastq_file))


def test_validate_fastq_file_wraps_unexpected_read_errors(tmp_path, caplog):
    """utils.py:375-377 -- an extension-valid but corrupt .gz surfaces through the generic
    except, logged and re-raised as the original exception type (not swallowed)."""
    fastq_file = tmp_path / "corrupt.fastq.gz"
    fastq_file.write_text("this is not gzip data")
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with pytest.raises(gzip.BadGzipFile):
        validate_fastq_file(str(fastq_file))

    assert _logged(caplog.records, "Error validating FASTQ file")


# ---------------------------------------------------------------------------
# get_tool_versions probes only what the run will invoke (#262)
# ---------------------------------------------------------------------------

#: Every tool ``config.json`` ships, so a test can name any subset as in use.
ALL_TOOLS = {
    "fastp": "fastp",
    "samtools": "samtools",
    "bwa": "bwa",
    "advntr": "mamba run -n envadvntr advntr",
    "shark": "mamba run -n shark_env shark",
    "kanalyze": "vntyper/dependencies/kestrel/kanalyze.jar",
    "kestrel": "vntyper/dependencies/kestrel/kestrel.jar",
    "java_path": "java",
}


def _probed(tools, tools_in_use):
    """Record the commands ``get_tool_versions`` actually shells out for.

    Args:
        tools: The ``config["tools"]`` mapping.
        tools_in_use: The names the caller declares this run will invoke.

    Returns:
        tuple[list[str], dict[str, str]]: The probed commands and the returned versions.
    """
    commands: list[str] = []

    def record(command, version_flag):
        commands.append(command)
        return "1.0"

    with patch("vntyper.scripts.utils.get_tool_version", side_effect=record):
        versions = get_tool_versions({"tools": tools}, tools_in_use=tools_in_use)
    return commands, versions


def test_only_the_named_tools_are_probed():
    """A Kestrel-only BAM run needs neither adVNTR, SHARK, fastp nor BWA.

    ``mamba run -n envadvntr advntr`` costs 315 ms on every run -- adVNTR's own import,
    not mamba's overhead -- and the result is only ever logged.
    """
    commands, versions = _probed(ALL_TOOLS, {"samtools", "kestrel", "java_path"})

    assert not any("advntr" in command or "shark" in command for command in commands)
    assert not any(command.startswith(("fastp", "bwa")) for command in commands)
    assert set(versions) == {"samtools", "kestrel", "java_path"}


def test_kanalyze_is_never_probed_even_when_it_is_in_use():
    """It is a JAR with no version flag of its own; probing it would execute the JAR.

    It ships and is versioned with kestrel, so kestrel's version already covers it.
    """
    commands, versions = _probed(ALL_TOOLS, {"samtools", "kanalyze", "kestrel"})

    assert not any("kanalyze" in command for command in commands)
    assert "kanalyze" not in versions


def test_shark_is_never_probed_even_when_it_is_in_use():
    """get_tool_version has no SHARK branch and returns "unknown" unconditionally.

    Gating a probe that can only ever answer "unknown" would spend 36 ms to learn
    nothing. Dropping it is the honest choice; a real answer needs a parser, not a gate.
    """
    commands, versions = _probed(ALL_TOOLS, {"samtools", "shark"})

    assert not any("shark" in command for command in commands)
    assert "shark" not in versions


def test_advntr_is_probed_when_it_is_in_use():
    """The saving must come from not running adVNTR, not from never reporting it."""
    commands, versions = _probed(ALL_TOOLS, {"samtools", "advntr"})

    assert any("advntr" in command for command in commands)
    assert versions["advntr"] == "1.0"


def test_a_replacement_config_without_kanalyze_still_works():
    """--config-path replaces the whole config, so the key may be absent (trap 2)."""
    tools = {name: command for name, command in ALL_TOOLS.items() if name != "kanalyze"}

    _, versions = _probed(tools, {"samtools", "kanalyze", "kestrel"})

    assert set(versions) == {"samtools", "kestrel"}


def test_a_name_in_use_that_the_config_does_not_declare_is_skipped():
    """The caller names what the run invokes; the config decides what exists."""
    _, versions = _probed({"samtools": "samtools"}, {"samtools", "bwa"})

    assert versions == {"samtools": "1.0"}


def test_the_kestrel_help_flag_is_still_folded_in_exactly_once():
    """The special case must survive the filtering, argv and all."""
    config = {"tools": {"kestrel": "/opt/kestrel/kestrel-1.0.1.jar", "java_path": "/usr/bin/java"}}

    with patch("vntyper.scripts.utils.subprocess.run", return_value=_completed_process()) as mocked_run:
        get_tool_versions(config, tools_in_use={"kestrel", "java_path"})

    kestrel_calls = [call for call in mocked_run.call_args_list if "-jar" in call.args[0]]
    assert len(kestrel_calls) == 1
    assert kestrel_calls[0].args[0] == ["/usr/bin/java", "-jar", "/opt/kestrel/kestrel-1.0.1.jar", "-h"]
