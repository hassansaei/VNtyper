"""
Contract tests for :func:`vntyper.scripts.utils.run_command`.

``run_command`` is the single choke point through which every external tool in the
pipeline is executed: Kestrel (Java), samtools, bwa, fastp, bcftools and both optional
modules all go through it. Three of the repository traps meet in this one function, and
before this file existed, deleting any of them broke no test at all:

* **Trap 9** - ``shell=True`` with ``executable="/bin/bash"`` is deliberate. The CRAM
  unmapped-read path builds a command containing bash process substitution, which no
  ``shell=False`` argv list can express.
* **Trap 7** - the ``cwd=`` passthrough. ``pipeline.py`` pins
  ``project_root = os.getcwd()`` and threads it into the Java and samtools calls,
  because every tool and reference path in ``config.json`` is relative to the process
  working directory.
* The failure log level. A non-critical failure is the pipeline's way of saying "this
  stage did not work but we are carrying on", and it must be visible at the default
  log level rather than only under ``--log-level DEBUG``.

Every test here mocks ``subprocess.Popen``; nothing in this file starts a real process.
"""

import logging
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from vntyper.scripts.utils import run_command

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit

UTILS_LOGGER = "vntyper.scripts.utils"

# A command that only bash can run: process substitution is a bash-ism and is exactly
# what the CRAM unmapped-read path emits.
BASH_ONLY_COMMAND = 'samtools view -b -h -T ref.fa sample.cram | tee >(samtools view -c > "counts.txt")'


def _fake_process(returncode: int = 0, output: tuple[bytes, ...] = (b"child stdout line\n",)) -> MagicMock:
    """
    Build a stand-in for the object ``subprocess.Popen`` returns.

    Args:
        returncode (int): The exit status the fake child reports.
        output (tuple[bytes, ...]): The lines the fake child writes to stdout.

    Returns:
        MagicMock: An object exposing the ``stdout``, ``wait`` and ``returncode``
        attributes that ``run_command`` uses.
    """
    process = MagicMock(name="Popen")
    process.stdout = list(output)
    process.wait.return_value = returncode
    process.returncode = returncode
    return process


def test_the_command_is_handed_to_bash_as_a_single_shell_string(tmp_path):
    """
    Trap 9: the command must reach bash as one string, via ``shell=True``.

    Rewriting this call to ``shell=False`` with an argv list - the usual hardening
    reflex - silently breaks the CRAM unmapped-read branch, because process
    substitution cannot be expressed as an argv list. Defence in depth for the
    injection surface belongs in the command builders (``shlex.quote``), not here.
    """
    log_file = tmp_path / "cmd.log"

    with patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process()) as popen:
        assert run_command(BASH_ONLY_COMMAND, str(log_file)) is True

    popen.assert_called_once()
    call_args, call_kwargs = popen.call_args
    assert call_args[0] == BASH_ONLY_COMMAND, "the command must be passed through verbatim"
    assert isinstance(call_args[0], str), "shell=True takes a string; an argv list means shell=False crept back in"
    assert call_kwargs["shell"] is True, "trap 9: shell=True is required for CRAM process substitution"
    assert call_kwargs["executable"] == "/bin/bash", "trap 9: /bin/bash is required; sh has no process substitution"


def test_the_child_output_is_captured_and_merged(tmp_path):
    """stderr must be folded into stdout so the log file holds the whole transcript."""
    log_file = tmp_path / "cmd.log"

    with patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process()) as popen:
        run_command("echo hi", str(log_file))

    _, call_kwargs = popen.call_args
    assert call_kwargs["stdout"] is subprocess.PIPE
    assert call_kwargs["stderr"] is subprocess.STDOUT


def test_the_cwd_is_passed_through_unchanged(tmp_path):
    """
    Trap 7: ``cwd`` must reach ``Popen`` exactly as the caller gave it.

    ``pipeline.py`` pins ``project_root = os.getcwd()`` and threads it into the Java
    and samtools calls. Every tool and reference path in ``config.json`` is relative to
    that directory, so normalising, resolving or dropping this value breaks reference
    lookup in a way no other test would catch.
    """
    project_root = tmp_path / "project_root"
    project_root.mkdir()
    log_file = tmp_path / "cmd.log"

    with patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process()) as popen:
        run_command("java -jar kestrel.jar", str(log_file), cwd=str(project_root))

    _, call_kwargs = popen.call_args
    assert call_kwargs["cwd"] == str(project_root)


def test_cwd_is_none_when_the_caller_omits_it(tmp_path):
    """Omitting ``cwd`` must let the child inherit the parent's working directory."""
    log_file = tmp_path / "cmd.log"

    with patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process()) as popen:
        run_command("samtools --version", str(log_file))

    _, call_kwargs = popen.call_args
    assert call_kwargs["cwd"] is None, "an omitted cwd must stay None, not become os.getcwd()"


def test_the_success_path_returns_true_and_logs_the_child_output(tmp_path):
    """A zero exit status returns True and the child's output lands in the log file."""
    log_file = tmp_path / "cmd.log"
    output = (b"[bwa] processing\n", b"[bwa] done\n")

    with patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process(0, output)):
        assert run_command("bwa mem ref.fa reads.fq", str(log_file)) is True

    written = log_file.read_text()
    assert "[bwa] processing" in written
    assert "[bwa] done" in written


def test_a_non_critical_failure_returns_false_and_still_writes_the_log(tmp_path):
    """A non-zero exit status with ``critical=False`` returns False rather than raising."""
    log_file = tmp_path / "cmd.log"

    with patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process(1, (b"E::main not found\n",))):
        assert run_command("samtools view missing.bam", str(log_file), critical=False) is False

    assert "E::main not found" in log_file.read_text(), "the diagnostic output must survive a failure"


def test_a_critical_failure_raises_runtime_error(tmp_path):
    """A non-zero exit status with ``critical=True`` aborts the pipeline."""
    log_file = tmp_path / "cmd.log"
    command = "samtools quickcheck -v broken.bam"

    with (
        patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process(1)),
        pytest.raises(RuntimeError) as excinfo,
    ):
        run_command(command, str(log_file), critical=True)

    assert command in str(excinfo.value), "the failing command must be named in the exception"


def test_a_non_critical_failure_is_logged_at_error_level(tmp_path, caplog):
    """
    A failed command must be reported at ERROR, so it is visible at the default level.

    This is the amplifier the rest of this effort keeps running into: a non-critical
    failure returns ``False`` and, when the failure was only logged at DEBUG, produced
    no output whatsoever at the pipeline's default INFO level. The stage carried on
    with missing input and the first sign of trouble appeared several stages later.

    The ``at_level(INFO)`` here is the point of the test, not incidental setup: it
    reproduces the level an ordinary run uses. ``pytest.ini`` sets ``log_level = DEBUG``
    globally, so without it caplog would capture a DEBUG record too and the test would
    pass while real runs stayed silent.
    """
    log_file = tmp_path / "cmd.log"
    command = "samtools view missing.bam"

    with (
        caplog.at_level(logging.INFO, logger=UTILS_LOGGER),
        patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process(1)),
    ):
        assert run_command(command, str(log_file), critical=False) is False

    errors = [record for record in caplog.records if record.levelno == logging.ERROR]
    assert errors, "a failed command must be logged at ERROR; at INFO a DEBUG log is silent"
    assert any(command in record.getMessage() for record in errors), "the ERROR log must name the failing command"


def test_a_critical_failure_is_logged_at_error_level_before_raising(tmp_path, caplog):
    """The critical path must log the failure at ERROR as well as raising."""
    log_file = tmp_path / "cmd.log"
    command = "java -jar kestrel.jar"

    with (
        caplog.at_level(logging.INFO, logger=UTILS_LOGGER),
        patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process(1)),
        pytest.raises(RuntimeError),
    ):
        run_command(command, str(log_file), critical=True)

    errors = [record for record in caplog.records if record.levelno == logging.ERROR]
    assert errors, "a critical failure must be logged at ERROR before the exception propagates"


def test_a_successful_command_logs_nothing_at_error_level(tmp_path, caplog):
    """The ERROR level stays reserved for failures; a clean run must not use it."""
    log_file = tmp_path / "cmd.log"

    with (
        caplog.at_level(logging.INFO, logger=UTILS_LOGGER),
        patch("vntyper.scripts.utils.subprocess.Popen", return_value=_fake_process(0)),
    ):
        assert run_command("samtools --version", str(log_file)) is True

    assert [record for record in caplog.records if record.levelno >= logging.ERROR] == []
