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
import stat
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
        self.closed = False

    def json(self) -> dict:
        """Return the decoded body.

        Returns:
            dict: The payload.
        """
        return self._payload

    def close(self) -> None:
        """Record that the client released this response."""
        self.closed = True


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


def test_status_get_retries_a_transient_connection_failure_and_closes_the_response(monkeypatch) -> None:
    """A dropped idempotent GET is retried, while its successful response is released."""
    response = FakeResponse({"job_id": "j", "status": "completed"})
    attempts = 0

    def _get(url, timeout=None):
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise online_mode.requests.ConnectionError("reset")
        return response

    monkeypatch.setattr(online_mode.requests, "get", _get)
    monkeypatch.setattr(online_mode.time, "sleep", lambda seconds: None)

    assert online_mode.poll_job_status("http://api", "j") == "completed"
    assert attempts == 2
    assert response.closed is True


def test_status_get_retries_are_bounded(monkeypatch) -> None:
    """A persistently unavailable status endpoint cannot loop beyond the retry bound."""
    attempts = 0

    def _get(url, timeout=None):
        nonlocal attempts
        attempts += 1
        raise online_mode.requests.Timeout("no answer")

    monkeypatch.setattr(online_mode.requests, "get", _get)
    waits: list[float] = []
    monkeypatch.setattr(online_mode.time, "sleep", waits.append)

    with pytest.raises(RuntimeError, match="failed after 3 attempts"):
        online_mode.poll_job_status("http://api", "j")

    assert attempts == 3
    assert waits == [2.0, 4.0]


def test_status_get_closes_a_discarded_5xx_before_retrying(monkeypatch) -> None:
    """Retrying a gateway response must not retain its connection in the pool."""
    unavailable = FakeResponse({}, status_code=503)
    completed = FakeResponse({"job_id": "j", "status": "completed"})
    responses = iter((unavailable, completed))
    monkeypatch.setattr(online_mode.requests, "get", lambda url, timeout=None: next(responses))
    monkeypatch.setattr(online_mode.time, "sleep", lambda seconds: None)

    assert online_mode.poll_job_status("http://api", "j") == "completed"
    assert unavailable.closed is True
    assert completed.closed is True


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


def test_a_non_200_status_response_still_raises_and_closes(monkeypatch, caplog) -> None:
    """Unchanged: an HTTP error is a hard failure, not a status.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
    """
    response = FakeResponse({}, status_code=404)
    monkeypatch.setattr(online_mode.requests, "get", lambda url, timeout=None: response)
    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="Failed to get job status") as excinfo:
        online_mode.poll_job_status("http://api", "j")

    assert response.closed is True
    assert any(record.levelno == logging.ERROR and record.message == str(excinfo.value) for record in caplog.records)


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
# submit_job and download_results
# --------------------------------------------------------------------------------------


def test_submission_posts_once_with_fresh_file_handles_at_byte_zero(tmp_path: Path, monkeypatch) -> None:
    """A response-loss exception cannot trigger a duplicate job submission."""
    bam_path = tmp_path / "subset.bam"
    bai_path = tmp_path / "subset.bam.bai"
    bam_path.write_bytes(b"bam-bytes")
    bai_path.write_bytes(b"bai-bytes")
    positions: list[tuple[int, int]] = []
    handles = []
    response = FakeResponse({"job_id": "job-1"})

    def _post(url, files, data, timeout):
        bam_handle = files["bam_file"][1]
        bai_handle = files["bai_file"][1]
        positions.append((bam_handle.tell(), bai_handle.tell()))
        handles.extend((bam_handle, bai_handle))
        bam_handle.read()
        bai_handle.read()
        return response

    post = mock.Mock(side_effect=_post)
    monkeypatch.setattr(online_mode.requests, "post", post)

    result = online_mode.submit_job("http://api", str(bam_path), str(bai_path), "hg19", 2)

    assert result == {"job_id": "job-1"}
    assert post.call_count == 1
    assert positions == [(0, 0)]
    assert all(handle.closed for handle in handles)
    assert response.closed is True


def test_submission_connection_failure_is_not_retried(tmp_path: Path, monkeypatch) -> None:
    """POST is unsafe to replay without a server idempotency key."""
    bam_path = tmp_path / "subset.bam"
    bam_path.write_bytes(b"bam-bytes")
    post = mock.Mock(side_effect=online_mode.requests.ConnectionError("response lost"))
    monkeypatch.setattr(online_mode.requests, "post", post)

    with pytest.raises(online_mode.requests.ConnectionError, match="response lost"):
        online_mode.submit_job("http://api", str(bam_path), None, "hg19", 2)

    assert post.call_count == 1


def test_submission_http_failure_is_logged_and_closed(tmp_path: Path, monkeypatch, caplog) -> None:
    """A rejected upload releases its response and records the same error it raises."""
    bam_path = tmp_path / "subset.bam"
    bam_path.write_bytes(b"bam-bytes")
    response = FakeResponse({"detail": "rejected"}, status_code=400)
    monkeypatch.setattr(online_mode.requests, "post", lambda url, files, data, timeout: response)

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="Failed to submit job") as excinfo:
        online_mode.submit_job("http://api", str(bam_path), None, "hg19", 2)

    assert response.closed is True
    assert any(record.levelno == logging.ERROR and record.message == str(excinfo.value) for record in caplog.records)


class _StreamingResponse(FakeResponse):
    """A response double whose body can fail between bounded chunks."""

    def __init__(self, chunks, *, status_code: int = 200, close_error: OSError | None = None) -> None:
        super().__init__({}, status_code=status_code)
        self._chunks = chunks
        self._close_error = close_error
        self.chunk_sizes: list[int] = []

    def iter_content(self, chunk_size: int):
        """Yield the configured chunks and record the requested upper bound."""
        self.chunk_sizes.append(chunk_size)
        yield from self._chunks

    def close(self) -> None:
        """Record the close and optionally model cleanup itself failing."""
        super().close()
        if self._close_error is not None:
            raise self._close_error


def test_download_streams_bounded_chunks_atomically_and_closes_response(tmp_path: Path, monkeypatch) -> None:
    """The result body is never buffered wholesale or exposed under its final name early."""
    response = _StreamingResponse([b"a" * 10, b"", b"b" * 5])
    seen: dict[str, object] = {}

    def _get(url, timeout=None, stream=False):
        seen["stream"] = stream
        return response

    monkeypatch.setattr(online_mode.requests, "get", _get)

    online_mode.download_results("http://api", "job1", tmp_path)

    final_path = tmp_path / "job1.zip"
    assert seen["stream"] is True
    assert response.chunk_sizes == [1 << 20]
    assert final_path.read_bytes() == b"a" * 10 + b"b" * 5
    assert stat.S_IMODE(final_path.stat().st_mode) & stat.S_IRUSR
    assert list(tmp_path.glob(".job1.*.tmp")) == []
    assert response.closed is True


def test_download_failure_removes_temp_and_preserves_existing_result(tmp_path: Path, monkeypatch) -> None:
    """A broken stream cannot replace a complete archive with partial bytes."""
    final_path = tmp_path / "job1.zip"
    final_path.write_bytes(b"previous-complete-result")

    def _chunks():
        yield b"partial"
        raise online_mode.requests.ConnectionError("stream reset")

    response = _StreamingResponse(_chunks())
    monkeypatch.setattr(online_mode.requests, "get", lambda url, timeout=None, stream=False: response)

    with pytest.raises(online_mode.requests.ConnectionError, match="stream reset"):
        online_mode.download_results("http://api", "job1", tmp_path)

    assert final_path.read_bytes() == b"previous-complete-result"
    assert list(tmp_path.glob(".job1.*.tmp")) == []
    assert response.closed is True


def test_download_close_failure_cannot_mask_stream_failure_or_skip_cleanup(tmp_path: Path, monkeypatch, caplog) -> None:
    """Cleanup stays best-effort when both the response stream and close fail."""
    final_path = tmp_path / "job1.zip"
    final_path.write_bytes(b"previous-complete-result")

    def _chunks():
        yield b"partial"
        raise online_mode.requests.ConnectionError("stream reset")

    response = _StreamingResponse(_chunks(), close_error=OSError("close refused"))
    monkeypatch.setattr(online_mode.requests, "get", lambda url, timeout=None, stream=False: response)

    with caplog.at_level(logging.WARNING), pytest.raises(online_mode.requests.ConnectionError, match="stream reset"):
        online_mode.download_results("http://api", "job1", tmp_path)

    assert final_path.read_bytes() == b"previous-complete-result"
    assert list(tmp_path.glob(".job1.*.tmp")) == []
    assert response.closed is True
    assert any(
        record.levelno == logging.WARNING
        and record.message == "Failed to close response from http://api/download/job1/: close refused"
        for record in caplog.records
    )


def test_download_replace_failure_removes_temp_and_preserves_existing_result(tmp_path: Path, monkeypatch) -> None:
    """Atomic publication failure leaves neither partial scratch data nor a damaged result."""
    final_path = tmp_path / "job1.zip"
    final_path.write_bytes(b"previous-complete-result")
    response = _StreamingResponse([b"new-complete-result"])
    monkeypatch.setattr(online_mode.requests, "get", lambda url, timeout=None, stream=False: response)
    monkeypatch.setattr(online_mode.os, "replace", mock.Mock(side_effect=OSError("replace refused")))

    with pytest.raises(OSError, match="replace refused"):
        online_mode.download_results("http://api", "job1", tmp_path)

    assert final_path.read_bytes() == b"previous-complete-result"
    assert list(tmp_path.glob(".job1.*.tmp")) == []
    assert response.closed is True


def test_download_retries_5xx_get_and_closes_every_response(tmp_path: Path, monkeypatch) -> None:
    """A transient gateway result is safely retried without leaking either response."""
    unavailable = _StreamingResponse([], status_code=502)
    completed = _StreamingResponse([b"zip"])
    responses = iter((unavailable, completed))
    get = mock.Mock(side_effect=lambda url, timeout=None, stream=False: next(responses))
    monkeypatch.setattr(online_mode.requests, "get", get)
    monkeypatch.setattr(online_mode.time, "sleep", lambda seconds: None)

    online_mode.download_results("http://api", "job1", tmp_path)

    assert get.call_count == 2
    assert all(call.kwargs["stream"] is True for call in get.call_args_list)
    assert unavailable.closed is True
    assert completed.closed is True


def test_download_retries_5xx_even_when_discarded_response_close_fails(tmp_path: Path, monkeypatch) -> None:
    """A 5xx close error cannot turn a safe retry into a download failure."""
    unavailable = _StreamingResponse([], status_code=502, close_error=OSError("close refused"))
    completed = _StreamingResponse([b"zip"])
    responses = iter((unavailable, completed))
    get = mock.Mock(side_effect=lambda url, timeout=None, stream=False: next(responses))
    monkeypatch.setattr(online_mode.requests, "get", get)
    monkeypatch.setattr(online_mode.time, "sleep", lambda seconds: None)

    online_mode.download_results("http://api", "job1", tmp_path)

    assert get.call_count == 2
    assert (tmp_path / "job1.zip").read_bytes() == b"zip"
    assert unavailable.closed is True
    assert completed.closed is True


def test_download_http_failure_is_logged_and_closed(tmp_path: Path, monkeypatch, caplog) -> None:
    """A factual 4xx is not retried, but its response is still released."""
    response = _StreamingResponse([], status_code=404)
    get = mock.Mock(return_value=response)
    monkeypatch.setattr(online_mode.requests, "get", get)

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="Failed to download results") as excinfo:
        online_mode.download_results("http://api", "job1", tmp_path)

    assert get.call_count == 1
    assert response.closed is True
    assert any(record.levelno == logging.ERROR and record.message == str(excinfo.value) for record in caplog.records)


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


@pytest.mark.parametrize("job_file_contents", [None, "", " \n\t"])
def test_jobless_resume_is_rejected_before_alignment_work(
    job_file_contents: str | None, tmp_path: Path, _stubbed_online
) -> None:
    """Missing and blank resume state cannot silently become a fresh local submission."""
    kwargs = _online_kwargs(tmp_path)
    output_path = Path(kwargs["output_dir"])
    output_path.mkdir(parents=True)
    if job_file_contents is not None:
        (output_path / "job_id.txt").write_text(job_file_contents, encoding="utf-8")

    with pytest.raises(ValueError, match="no job to resume"):
        online_mode.run_online_mode(**kwargs, resume=True)

    assert _stubbed_online["get_region_string_with_fallback"].call_count == 0
    assert _stubbed_online["subset_bam"].call_count == 0
    assert _stubbed_online["submit_job"].call_count == 0


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
