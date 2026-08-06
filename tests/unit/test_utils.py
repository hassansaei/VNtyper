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


def test_get_tool_version_returns_unknown_and_logs_on_unexpected_error(caplog):
    """utils.py:197-199 -- the catch-all: an unbalanced quote in ``command`` makes
    shlex.split raise ValueError before a process is ever started."""
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)
    with patch("vntyper.scripts.utils.subprocess.run") as mock_run:
        assert get_tool_version('fastp "unterminated', "") == "unknown"
    mock_run.assert_not_called()
    assert _logged(caplog.records, "Failed to get version for")


# ---------------------------------------------------------------------------
# get_tool_versions
# ---------------------------------------------------------------------------


def test_get_tool_versions_maps_each_configured_tool_to_its_parsed_version():
    """utils.py:213-235 -- the config-driven fan-out to get_tool_version."""
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
        versions = get_tool_versions(config)

    assert versions == {
        "fastp": "/usr/bin/fastp::",
        "samtools": "/usr/bin/samtools::",
        "java_path": "/usr/bin/java::--version",
    }
    mocked.assert_any_call("/usr/bin/fastp", "")
    mocked.assert_any_call("/usr/bin/java", "--version")


def test_get_tool_versions_with_no_tools_configured_returns_an_empty_dict():
    """utils.py:213-214 -- config.get("tools", {}) with nothing configured."""
    assert get_tool_versions({}) == {}


def test_get_tool_versions_kestrel_command_duplicates_the_help_flag_today():
    """
    CHARACTERISATION of a live defect. Do not "fix" this here.

    get_tool_versions builds ``command`` for kestrel as
    ``f"{java_path} {version_flag}"`` (utils.py:232), which already embeds the
    ``-jar "<path>" -h`` flag, and then calls
    ``get_tool_version(command, version_flag)`` with that SAME ``version_flag``
    again (utils.py:233). get_tool_version concatenates
    ``shlex.split(command) + shlex.split(version_flag)`` (utils.py:155), so the
    argv handed to ``subprocess.run`` duplicates ``-jar <path> -h``: it becomes
    ``java -jar k.jar -h -jar k.jar -h`` instead of ``java -jar k.jar -h``.

    This does not currently crash -- kestrel's arg parser tolerates the extra
    tokens in the cases observed in this repo's fixtures -- but it is not the
    intended invocation, and a stricter arg parser would make version
    detection silently report "unknown" (get_tool_version's blanket
    ``except Exception`` -> "unknown" is exactly the "silently wrong answer"
    shape this codebase keeps producing). Filed as
    issue-utils-kestrel-duplicate-version-flag.md; not fixed here per
    AGENT-RULES section 8 and the orchestrator's do-not-fix instruction.
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
        get_tool_versions(config)

    kestrel_calls = [call for call in mocked_run.call_args_list if "-jar" in call.args[0]]
    assert len(kestrel_calls) == 1, "expected exactly one subprocess.run call built for kestrel"
    argv = kestrel_calls[0].args[0]

    assert argv == [
        "/usr/bin/java",
        "-jar",
        "/opt/kestrel/kestrel-1.0.1.jar",
        "-h",
        "-jar",
        "/opt/kestrel/kestrel-1.0.1.jar",
        "-h",
    ], "the -jar/-h pair should appear exactly once if this is ever fixed"


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


def test_validate_bam_file_quickcheck_failure_raises_runtimeerror_not_valueerror_today(tmp_path):
    """
    CHARACTERISATION of a live defect. Do not "fix" this here.

    validate_bam_file's own docstring promises
    ``Raises: ValueError: If any validation check fails.``, and utils.py:334-337
    reads as though a failed quickcheck raises ValueError. But the call at
    utils.py:334 hardcodes ``critical=True``, and run_command's contract (see
    ``test_run_command_contract.py::test_a_critical_failure_raises_runtime_error``)
    is: a critical failure RAISES RuntimeError, it never returns False. So on a
    real quickcheck failure, ``if not success:`` at utils.py:335 is
    unreachable -- the caller gets ``RuntimeError("Critical command failed: ...")``
    instead of the documented ValueError, and an ``except ValueError`` around a
    ``validate_bam_file()`` call will not catch it.

    ``pipeline.py``'s own call sites (pipeline.py:219,221) sit inside a
    ``try`` whose matching clause is a blanket ``except Exception:``
    (pipeline.py:715), so this is currently benign there; any other caller
    that catches ``ValueError`` specifically is not protected. Filed as
    issue-utils-validate-bam-file-wrong-exception-type.md; not fixed here.
    """
    bam_file = tmp_path / "broken.bam"
    bam_file.touch()
    failing_process = MagicMock()
    failing_process.stdout = [b"[E::hts_open] truncated file\n"]
    failing_process.wait.return_value = 1
    failing_process.returncode = 1

    with (
        patch("vntyper.scripts.utils.subprocess.Popen", return_value=failing_process),
        pytest.raises(RuntimeError, match="Critical command failed"),
    ):
        validate_bam_file(str(bam_file))


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


def test_validate_fastq_file_wraps_unexpected_read_errors(tmp_path, caplog):
    """utils.py:375-377 -- an extension-valid but corrupt .gz surfaces through the generic
    except, logged and re-raised as the original exception type (not swallowed)."""
    fastq_file = tmp_path / "corrupt.fastq.gz"
    fastq_file.write_text("this is not gzip data")
    caplog.set_level(logging.ERROR, logger=UTILS_LOGGER)

    with pytest.raises(gzip.BadGzipFile):
        validate_fastq_file(str(fastq_file))

    assert _logged(caplog.records, "Error validating FASTQ file")
