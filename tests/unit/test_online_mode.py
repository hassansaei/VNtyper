"""``vntyper online``: a failed remote job must not look like a successful run.

Two defects, both of which turn a failure into a green exit code:

* ``run_online_mode`` logged ``"Job failed or status unknown."`` and **returned**.
  ``handle_online`` returned, ``main`` returned, the process exited 0, and a wrapping
  ``subprocess.run(..., check=True)`` saw success. Nothing downstream could tell a
  failed job from a completed one except by looking for the zip.
* ``poll_job_status`` looped forever on any status other than ``completed`` or
  ``failed``. The web service maps Celery's ``PENDING``/``STARTED`` to
  ``pending``/``started`` and everything else through ``status.lower()``, so a
  ``REVOKED`` job - a terminal state - polled until the operator killed the client.

The tests also pin what the client reads out of ``/job-status/``. That endpoint now
returns a **generic** ``error`` message for a failed job rather than the exception
detail, and validates job ids as UUIDs. The client must not depend on the old shape.
"""

import logging
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts import cli_handlers, online_mode

pytestmark = pytest.mark.unit


class FakeResponse:
    """The subset of ``requests.Response`` the online client uses.

    Attributes:
        status_code: HTTP status.
        _payload: The decoded JSON body.
        text: The raw body, used only in error messages.
    """

    def __init__(self, payload: dict, status_code: int = 200) -> None:
        self.status_code = status_code
        self._payload = payload
        self.text = str(payload)
        self.content = b"zip-bytes"

    def json(self) -> dict:
        """Return the decoded body.

        Returns:
            dict: The payload.
        """
        return self._payload


def _statuses(*sequence: str):
    """Build a ``requests.get`` stub that walks ``sequence`` then repeats the last.

    Args:
        *sequence: The ``status`` values to return, in order.

    Returns:
        Callable: A stub suitable for ``mock.patch``.
    """
    remaining = list(sequence)

    def _get(url, timeout=None):
        value = remaining.pop(0) if len(remaining) > 1 else remaining[0]
        return FakeResponse({"job_id": "j", "status": value})

    return _get


# --------------------------------------------------------------------------------------
# poll_job_status
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("terminal", ["completed", "failed"])
def test_a_terminal_status_is_returned_immediately(terminal: str, monkeypatch) -> None:
    """The happy paths, unchanged.

    Args:
        terminal: The status the server reports.
        monkeypatch: Pytest monkeypatch fixture.
    """
    monkeypatch.setattr(online_mode.requests, "get", _statuses(terminal))
    assert online_mode.poll_job_status("http://api", "j") == terminal


def test_polling_continues_through_the_in_progress_statuses(monkeypatch) -> None:
    """``pending`` and ``started`` are not answers, so keep asking.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
    """
    slept: list[float] = []
    monkeypatch.setattr(online_mode.time, "sleep", slept.append)
    monkeypatch.setattr(online_mode.requests, "get", _statuses("pending", "started", "retry", "completed"))

    assert online_mode.poll_job_status("http://api", "j", poll_interval=3) == "completed"
    assert slept == [3, 3, 3]


def test_an_unrecognised_status_is_terminal_rather_than_an_endless_loop(monkeypatch, caplog) -> None:
    """A ``REVOKED`` job used to poll forever; it is terminal and must be reported.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
        caplog: Pytest log capture.
    """
    monkeypatch.setattr(online_mode.time, "sleep", mock.Mock())
    monkeypatch.setattr(online_mode.requests, "get", _statuses("revoked"))

    with caplog.at_level(logging.WARNING):
        assert online_mode.poll_job_status("http://api", "j") == "revoked"

    assert any(record.levelno >= logging.WARNING and "revoked" in record.message for record in caplog.records)


def test_polling_gives_up_instead_of_running_forever(monkeypatch, caplog) -> None:
    """The cap: a job that never leaves ``pending`` must end the client, not hang.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
        caplog: Pytest log capture.
    """
    slept: list[float] = []
    monkeypatch.setattr(online_mode.time, "sleep", slept.append)
    monkeypatch.setattr(online_mode.requests, "get", _statuses("pending"))

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="did not reach a terminal status"):
        online_mode.poll_job_status("http://api", "j", poll_interval=5, timeout=20)

    assert len(slept) <= 4, f"polled {len(slept)} times for a 20s timeout at 5s intervals"
    assert any(record.levelno >= logging.ERROR for record in caplog.records)


def test_a_non_200_status_response_still_raises(monkeypatch) -> None:
    """Unchanged: an HTTP error is a hard failure, not a status.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
    """
    monkeypatch.setattr(online_mode.requests, "get", lambda url, timeout=None: FakeResponse({}, status_code=404))
    with pytest.raises(RuntimeError, match="Failed to get job status"):
        online_mode.poll_job_status("http://api", "j")


def test_the_generic_failure_message_is_surfaced_not_parsed(monkeypatch, caplog) -> None:
    """``/job-status/`` returns a generic ``error`` for a failed job; log it verbatim.

    The endpoint deliberately no longer returns the task's exception detail. The
    client must read ``status`` for its decision and treat ``error`` as opaque text.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
        caplog: Pytest log capture.
    """
    generic = "Job failed. Please contact the administrator."

    def _get(url, timeout=None):
        return FakeResponse({"job_id": "j", "status": "failed", "error": generic})

    monkeypatch.setattr(online_mode.requests, "get", _get)
    with caplog.at_level(logging.ERROR):
        assert online_mode.poll_job_status("http://api", "j") == "failed"

    assert any(generic in record.message for record in caplog.records)


def test_the_client_does_not_read_any_field_the_endpoint_stopped_returning() -> None:
    """Pin the response contract: ``status`` decides, ``error`` is opaque, nothing else.

    A client that parsed the old exception detail would be a live regression against
    the hardened endpoint. This asserts the source names no other field.
    """
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(online_mode.poll_job_status))
    read = {
        node.args[0].value
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "get"
        and node.args
        and isinstance(node.args[0], ast.Constant)
        and isinstance(node.args[0].value, str)
    }
    assert read, "found no dict lookups in poll_job_status; this test would be vacuous"
    assert read <= {"status", "error"}, f"poll_job_status reads unexpected response fields: {sorted(read)}"


# --------------------------------------------------------------------------------------
# run_online_mode
# --------------------------------------------------------------------------------------


def _online_kwargs(tmp_path: Path) -> dict:
    """Return the arguments ``run_online_mode`` needs for a fresh submission.

    Args:
        tmp_path: Pytest temporary directory.

    Returns:
        dict: Keyword arguments.
    """
    return {
        "config": {"api": {"base_url": "http://api"}, "bam_processing": {"bam_region_hg19": "chr1:1-2"}},
        "bam": str(tmp_path / "in.bam"),
        "output_dir": str(tmp_path / "out"),
        "reference_assembly": "hg19",
        "threads": 2,
    }


@pytest.fixture
def _stubbed_online(monkeypatch):
    """Stub every network and samtools call ``run_online_mode`` makes.

    Args:
        monkeypatch: Pytest monkeypatch fixture.

    Yields:
        dict: The installed stubs, by name.
    """
    stubs = {
        "subset_bam": mock.Mock(),
        "submit_job": mock.Mock(return_value={"job_id": "job-1"}),
        "poll_job_status": mock.Mock(return_value="completed"),
        "download_results": mock.Mock(),
        "get_region_string_with_fallback": mock.Mock(return_value="chr1:1-2"),
    }
    for name, stub in stubs.items():
        monkeypatch.setattr(online_mode, name, stub)
    yield stubs


def test_a_completed_job_downloads_its_results(tmp_path: Path, _stubbed_online) -> None:
    """The happy path, unchanged.

    Args:
        tmp_path: Pytest temporary directory.
        _stubbed_online: The installed stubs.
    """
    online_mode.run_online_mode(**_online_kwargs(tmp_path))
    assert _stubbed_online["download_results"].call_count == 1


@pytest.mark.parametrize("status", ["failed", "revoked", ""])
def test_a_job_that_did_not_complete_raises(status: str, tmp_path: Path, _stubbed_online, caplog) -> None:
    """The defect: this used to log and return, so the process exited 0.

    Args:
        status: The terminal status the poller reports.
        tmp_path: Pytest temporary directory.
        _stubbed_online: The installed stubs.
        caplog: Pytest log capture.
    """
    _stubbed_online["poll_job_status"].return_value = status

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="job-1"):
        online_mode.run_online_mode(**_online_kwargs(tmp_path))

    assert _stubbed_online["download_results"].call_count == 0
    assert any(record.levelno >= logging.ERROR for record in caplog.records)


def test_a_submission_with_no_job_id_raises(tmp_path: Path, _stubbed_online) -> None:
    """A response the API should never send used to end the run quietly.

    Args:
        tmp_path: Pytest temporary directory.
        _stubbed_online: The installed stubs.
    """
    _stubbed_online["submit_job"].return_value = {"message": "queued"}
    with pytest.raises(RuntimeError, match="job_id"):
        online_mode.run_online_mode(**_online_kwargs(tmp_path))


def test_a_resumed_job_that_failed_also_raises(tmp_path: Path, _stubbed_online) -> None:
    """``--resume`` has its own copy of the completed/failed branch.

    Args:
        tmp_path: Pytest temporary directory.
        _stubbed_online: The installed stubs.
    """
    kwargs = _online_kwargs(tmp_path)
    output_path = Path(kwargs["output_dir"])
    output_path.mkdir(parents=True)
    (output_path / "job_id.txt").write_text("old-job", encoding="utf-8")
    _stubbed_online["poll_job_status"].return_value = "failed"

    with pytest.raises(RuntimeError, match="old-job"):
        online_mode.run_online_mode(**kwargs, resume=True)

    assert _stubbed_online["submit_job"].call_count == 0


def test_a_resumed_job_that_completed_downloads(tmp_path: Path, _stubbed_online) -> None:
    """The other half of the resume branch.

    Args:
        tmp_path: Pytest temporary directory.
        _stubbed_online: The installed stubs.
    """
    kwargs = _online_kwargs(tmp_path)
    output_path = Path(kwargs["output_dir"])
    output_path.mkdir(parents=True)
    (output_path / "job_id.txt").write_text("old-job", encoding="utf-8")

    online_mode.run_online_mode(**kwargs, resume=True)

    assert _stubbed_online["download_results"].call_count == 1
    assert _stubbed_online["submit_job"].call_count == 0


def test_the_poll_interval_and_timeout_come_from_the_config(tmp_path: Path, _stubbed_online) -> None:
    """An operator must be able to raise the cap without editing the source.

    Args:
        tmp_path: Pytest temporary directory.
        _stubbed_online: The installed stubs.
    """
    kwargs = _online_kwargs(tmp_path)
    kwargs["config"]["api"] = {"base_url": "http://api", "poll_interval_seconds": 2, "poll_timeout_seconds": 90}

    online_mode.run_online_mode(**kwargs)

    call = _stubbed_online["poll_job_status"].call_args
    assert call.kwargs["poll_interval"] == 2
    assert call.kwargs["timeout"] == 90


# --------------------------------------------------------------------------------------
# handle_online: the exit code a wrapping script actually sees
# --------------------------------------------------------------------------------------


def test_the_handler_exits_non_zero_when_the_job_failed(tmp_path: Path, monkeypatch) -> None:
    """The whole point: ``subprocess.run(..., check=True)`` must see the failure.

    Args:
        tmp_path: Pytest temporary directory.
        monkeypatch: Pytest monkeypatch fixture.
    """
    monkeypatch.setattr(
        cli_handlers,
        "run_online_mode",
        mock.Mock(side_effect=RuntimeError("Online job job-1 did not complete: status 'failed'.")),
    )
    args = mock.Mock(output_dir=str(tmp_path), reference_assembly="hg19", threads=2, bam="in.bam")
    args.email = args.cohort_id = args.passphrase = None
    args.resume = False

    with pytest.raises(SystemExit) as excinfo:
        cli_handlers.handle_online(args, config={}, parser=mock.Mock(), log_level_value=logging.INFO, log_file_str=None)

    assert excinfo.value.code == 1


def test_the_handler_returns_normally_when_the_job_completed(tmp_path: Path, monkeypatch) -> None:
    """A completed job must still exit 0.

    Args:
        tmp_path: Pytest temporary directory.
        monkeypatch: Pytest monkeypatch fixture.
    """
    monkeypatch.setattr(cli_handlers, "run_online_mode", mock.Mock())
    args = mock.Mock(output_dir=str(tmp_path), reference_assembly="hg19", threads=2, bam="in.bam")
    args.email = args.cohort_id = args.passphrase = None
    args.resume = False

    cli_handlers.handle_online(args, config={}, parser=mock.Mock(), log_level_value=logging.INFO, log_file_str=None)
