"""Response hygiene for the web service: cleanup after, and honest status codes.

Two properties are pinned here, both about what an endpoint leaves behind once it
has answered.

1. ``/cohort-download/`` builds its archive in a temporary file. That file is the
   service's, not the caller's, so nothing of it may survive the response — and
   the response itself must still carry the complete archive, so the cleanup has
   to happen *after* the body is sent rather than before it. Both halves are
   asserted together: an implementation that deletes too early fails the
   integrity check, and one that never deletes fails the filesystem check. The
   ways a response ends are covered one by one — delivered, failed while being
   built, and abandoned mid-body — because the archive holds a whole cohort's
   results and "usually removed" is not the guarantee being claimed.

2. ``/job-queue/`` must answer an unknown ``job_id`` with 404. A handler that
   catches its own ``HTTPException`` in a broad ``except Exception`` turns that
   404 into a 500, which tells the caller "the service is broken" when the truth
   is "no such job". The reverse must stay true as well: a genuinely unexpected
   failure inside the handler still has to surface as 500, so re-raising must not
   become a way of hiding real breakage.

The temporary directory is redirected at ``tempfile.tempdir`` so the leak check
is an assertion about an empty directory this test owns, rather than a diff of
the machine's shared ``/tmp``.

``docker`` is put on ``sys.path`` by ``tests/unit/web/conftest.py``, which pytest
imports before this module.
"""

import asyncio
import io
import os
import tempfile
import zipfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

pytestmark = pytest.mark.unit

from app.cohorts import COHORT_KEY_PREFIX  # noqa: E402

GOOD_PASSPHRASE = "correct-horse-battery-staple"

# A job identifier in the form the service issues, which is the only form the
# job routes accept -- see `test_job_identifiers.py`.
QUEUED_JOB_ID = "b7c41d90-2e58-4a63-9f07-15c8de3b6a24"

# Recognisable bytes so the delivered archive can be compared to its source.
RESULT_PAYLOAD = b"vntyper-result-bytes-for-job-a"
PATIENT_PAYLOAD = b"patient-bytes-that-must-not-be-served"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _spool_dir(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Point `tempfile` at an empty directory owned by the calling test.

    `tempfile.gettempdir()` caches into the module global patched here, so every
    `NamedTemporaryFile` the request opens lands in this directory and nothing
    else does.

    Args:
        monkeypatch: Standard pytest fixture; restores the global at teardown.
        tmp_path: Per-test scratch directory.

    Returns:
        Path: The now-empty directory temporary files will be created in.
    """
    spool = tmp_path / "spool"
    spool.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(spool))
    return spool


def _cohort_with_one_result(client, fake_redis, tmp_path: Path) -> str:
    """Create a protected cohort holding one job that has a result archive.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the `web_app` fixture configured as the job tree.

    Returns:
        str: The identifier of the created cohort.
    """
    response = client.post("/create-cohort/", data={"alias": "study", "passphrase": GOOD_PASSPHRASE})
    assert response.status_code == 200, response.text
    cohort_id = response.json()["cohort_id"]

    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")
    (tmp_path / "output" / "job-a.zip").write_bytes(RESULT_PAYLOAD)
    return cohort_id


def _send_response(response) -> bytes:
    """Run one response as ASGI and return its body bytes."""
    messages = []

    async def receive():
        return {"type": "http.request", "body": b"", "more_body": False}

    async def send(message):
        messages.append(message)

    scope = {
        "type": "http",
        "asgi": {"version": "3.0", "spec_version": "2.1"},
        "http_version": "1.1",
        "method": "GET",
        "scheme": "http",
        "path": "/",
        "raw_path": b"/",
        "query_string": b"",
        "root_path": "",
        "headers": [],
        "client": ("127.0.0.1", 50000),
        "server": ("testserver", 80),
    }
    asyncio.run(response(scope, receive, send))
    return b"".join(message.get("body", b"") for message in messages)


def test_individual_download_is_bound_before_its_path_is_replaced(web_app, tmp_path: Path) -> None:
    """A response streams the inode opened by the endpoint, not a later symlink target."""
    job_id = "b7c41d90-2e58-4a63-9f07-15c8de3b6a24"
    archive = tmp_path / "output" / f"{job_id}.zip"
    archive.write_bytes(RESULT_PAYLOAD)
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_PAYLOAD)

    response = web_app.download_result(job_id)
    archive.unlink()
    archive.symlink_to(patient)

    assert _send_response(response) == RESULT_PAYLOAD
    assert archive.is_symlink()
    assert patient.read_bytes() == PATIENT_PAYLOAD


def test_individual_download_rejects_a_hard_link_to_external_bytes(client, tmp_path: Path) -> None:
    """Existence alone cannot make a multiply-linked external inode downloadable."""
    job_id = "b7c41d90-2e58-4a63-9f07-15c8de3b6a24"
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_PAYLOAD)
    os.link(patient, tmp_path / "output" / f"{job_id}.zip")

    response = client.get(f"/download/{job_id}/")

    assert response.status_code == 404
    assert PATIENT_PAYLOAD not in response.content
    assert patient.read_bytes() == PATIENT_PAYLOAD


def test_individual_download_opens_a_fifo_nonblocking_before_rejecting_it(
    client, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A final-component FIFO swap cannot block the API worker before ``fstat``."""
    from app import archive_delivery

    job_id = "b7c41d90-2e58-4a63-9f07-15c8de3b6a24"
    fifo = tmp_path / "output" / f"{job_id}.zip"
    os.mkfifo(fifo)
    original_open = archive_delivery.os.open

    def require_nonblocking(path, flags, *args, **kwargs):
        if str(path) == fifo.name:
            assert flags & os.O_NONBLOCK
        return original_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(archive_delivery.os, "open", require_nonblocking)

    response = client.get(f"/download/{job_id}/")

    assert response.status_code == 404
    assert fifo.exists()


# ---------------------------------------------------------------------------
# /cohort-download/ cleans up after itself
# ---------------------------------------------------------------------------


def test_cohort_download_leaves_no_temporary_file_behind(
    client, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The archive the endpoint builds does not outlive the response.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        monkeypatch: Standard pytest fixture, used to redirect `tempfile`.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    spool = _spool_dir(monkeypatch, tmp_path)
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    response = client.get("/cohort-download/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200, response.text
    assert response.content  # the body was consumed, so any cleanup has run
    assert sorted(p.name for p in spool.iterdir()) == []


def test_cohort_download_still_delivers_the_complete_archive(
    client, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Cleanup runs after the body, so the caller still gets an intact zip.

    Deleting the temporary file before the response is sent would also leave the
    directory empty, which is why this test reads the archive back rather than
    checking the status code.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        monkeypatch: Standard pytest fixture, used to redirect `tempfile`.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    spool = _spool_dir(monkeypatch, tmp_path)
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    response = client.get("/cohort-download/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200, response.text
    assert response.headers["content-type"] == "application/zip"

    with zipfile.ZipFile(io.BytesIO(response.content)) as delivered:
        assert delivered.testzip() is None
        assert delivered.namelist() == ["job-a.zip"]
        assert delivered.read("job-a.zip") == RESULT_PAYLOAD

    # ... and the archive that produced those bytes is gone by now.
    assert sorted(p.name for p in spool.iterdir()) == []


def test_cohort_download_is_bound_before_its_scratch_path_is_replaced(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The cohort response streams its built archive, not a substituted patient alias."""
    spool = _spool_dir(monkeypatch, tmp_path)
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)
    response = web_app.cohort_download(
        cohort_id=cohort_id,
        alias=None,
        passphrase=GOOD_PASSPHRASE,
        header_passphrase=None,
    )
    scratch = next(spool.iterdir())
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_PAYLOAD)
    scratch.unlink()
    scratch.symlink_to(patient)

    delivered = _send_response(response)

    with zipfile.ZipFile(io.BytesIO(delivered)) as archive:
        assert archive.read("job-a.zip") == RESULT_PAYLOAD
    assert scratch.is_symlink()
    assert patient.read_bytes() == PATIENT_PAYLOAD


def test_cohort_download_does_not_package_a_hard_linked_external_member(client, fake_redis, tmp_path: Path) -> None:
    """A multiply-linked member is excluded instead of packaging external bytes."""
    response = client.post("/create-cohort/", data={"alias": "study", "passphrase": GOOD_PASSPHRASE})
    cohort_id = response.json()["cohort_id"]
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_PAYLOAD)
    os.link(patient, tmp_path / "output" / "job-a.zip")

    response = client.get("/cohort-download/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200
    assert PATIENT_PAYLOAD not in response.content
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        assert archive.namelist() == []
    assert patient.read_bytes() == PATIENT_PAYLOAD


def test_cohort_download_discards_a_member_that_gains_a_hard_link_during_copy(
    client, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Post-copy ownership validation happens before a ZIP member is committed."""
    from app import archive_delivery

    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)
    member = tmp_path / "output" / "job-a.zip"
    external_alias = tmp_path / "external-alias.zip"
    original_copy = archive_delivery.shutil.copyfileobj
    linked = False

    def copy_then_link(source, target, *args, **kwargs):
        nonlocal linked
        result = original_copy(source, target, *args, **kwargs)
        if not linked:
            linked = True
            os.link(member, external_alias)
        return result

    monkeypatch.setattr(archive_delivery.shutil, "copyfileobj", copy_then_link)

    response = client.get("/cohort-download/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        assert archive.namelist() == []
    assert external_alias.read_bytes() == RESULT_PAYLOAD


def test_cohort_download_opens_a_fifo_nonblocking_before_skipping_it(
    client, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A FIFO member is rejected without waiting for a writer."""
    from app import archive_delivery

    response = client.post("/create-cohort/", data={"alias": "study", "passphrase": GOOD_PASSPHRASE})
    cohort_id = response.json()["cohort_id"]
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")
    fifo = tmp_path / "output" / "job-a.zip"
    os.mkfifo(fifo)
    original_open = archive_delivery.os.open

    def require_nonblocking(path, flags, *args, **kwargs):
        if str(path) == fifo.name:
            assert flags & os.O_NONBLOCK
        return original_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(archive_delivery.os, "open", require_nonblocking)

    response = client.get("/cohort-download/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        assert archive.namelist() == []
    assert fifo.exists()


def test_cohort_download_leaves_no_temporary_file_when_the_archive_cannot_be_built(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A failure part-way through the build does not strand the scratch file.

    The archive is assembled member by member, and a member's result can age out
    between being seen and being read. The cleanup is attached to the response,
    so a build that never produces a response has to reclaim the file itself --
    otherwise the one path that leaves a temporary file behind is the one nobody
    is watching.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched API module whose archive writer is made to fail.
        fake_redis: The store backing the app's cohort client.
        monkeypatch: Standard pytest fixture, used to redirect `tempfile`.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    spool = _spool_dir(monkeypatch, tmp_path)
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    def _archive_write_failed(*args, **kwargs):
        """Stand in for destination archive I/O failing mid-build.

        Args:
            *args: Ignored.
            **kwargs: Ignored.

        Raises:
            RuntimeError: Always.
        """
        raise RuntimeError("archive write failed")

    monkeypatch.setattr(web_app, "write_owned_zip_member", _archive_write_failed)

    with pytest.raises(RuntimeError, match="archive write failed"):
        client.get("/cohort-download/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert sorted(p.name for p in spool.iterdir()) == []


def test_cohort_download_removes_the_archive_when_the_send_fails(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A connection dropped mid-body still takes the archive with it.

    Starlette awaits every body send *before* it runs a response's background
    task, and it does so without a `finally`. A client that disconnects part-way
    through the download therefore skipped the cleanup entirely, and what was
    left on disk was a ZIP of the whole cohort's results.

    Driven as raw ASGI rather than through TestClient, because a failing `send`
    is the condition under test and TestClient will not produce one.

    Args:
        client: TestClient fixture from conftest, used to create the cohort.
        web_app: The patched `app.main` module, holding the ASGI app.
        fake_redis: The store backing the app's cohort client.
        monkeypatch: Standard pytest fixture, used to redirect `tempfile`.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    import asyncio

    spool = _spool_dir(monkeypatch, tmp_path)
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    scope = {
        "type": "http",
        "asgi": {"version": "3.0", "spec_version": "2.1"},
        "http_version": "1.1",
        "method": "GET",
        "scheme": "http",
        "path": "/cohort-download/",
        "raw_path": b"/cohort-download/",
        "query_string": f"cohort_id={cohort_id}&passphrase={GOOD_PASSPHRASE}".encode(),
        "root_path": "",
        "headers": [(b"host", b"testserver")],
        "client": ("127.0.0.1", 50000),
        "server": ("testserver", 80),
    }

    async def receive():
        """Return the (empty) request body.

        Returns:
            dict: A single complete `http.request` message.
        """
        return {"type": "http.request", "body": b"", "more_body": False}

    async def send(message):
        """Accept the response start, then drop the connection.

        Args:
            message: The ASGI message being sent.

        Raises:
            ConnectionResetError: On the first body chunk.
        """
        if message["type"] == "http.response.body":
            msg = "the client went away mid-download"
            raise ConnectionResetError(msg)

    with pytest.raises(ConnectionResetError, match="went away"):
        asyncio.run(web_app.app(scope, receive, send))

    assert sorted(p.name for p in spool.iterdir()) == [], "a ZIP of the cohort's results was left on disk"


# ---------------------------------------------------------------------------
# /job-queue/ reports what actually happened
# ---------------------------------------------------------------------------


def test_job_queue_reports_an_unknown_job_id_as_not_found(client) -> None:
    """An identifier with no job behind it is a 404, not a server error.

    The identifier is well-formed and simply unknown, so the 404 has to come from
    the store lookup rather than from the format check that precedes it -- this
    is the handler re-raising its own `HTTPException` instead of catching it in
    its broad handler, which is what the test is for.

    Args:
        client: TestClient fixture from conftest.
    """
    response = client.get("/job-queue/", params={"job_id": QUEUED_JOB_ID})

    assert response.status_code == 404
    assert response.json()["detail"] == "Job ID not found"


def test_job_queue_still_reports_a_real_failure_as_a_server_error(client, web_app) -> None:
    """Re-raising the handler's own 404 must not hide genuine breakage.

    Args:
        client: TestClient fixture from conftest.
        web_app: Patched `app.main` module, whose Redis client is replaced here.
    """
    broken = MagicMock()
    broken.lrange.return_value = []
    broken.get.side_effect = RuntimeError("redis went away mid-request")
    web_app.redis_client = broken

    response = client.get("/job-queue/", params={"job_id": QUEUED_JOB_ID})

    assert response.status_code == 500
    assert response.json()["detail"] == "Error retrieving job position"


def test_job_queue_reports_the_position_of_a_queued_job(client, fake_redis) -> None:
    """The successful lookup path is unchanged.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's job client.
    """
    fake_redis.set(QUEUED_JOB_ID, "task-a")
    fake_redis.rpush("vntyper_job_queue", "task-z", "task-a")

    response = client.get("/job-queue/", params={"job_id": QUEUED_JOB_ID})

    assert response.status_code == 200
    body = response.json()
    assert body["job_id"] == QUEUED_JOB_ID
    assert body["position_in_queue"] == 2
    assert body["total_jobs_in_queue"] == 2


def test_job_queue_reports_a_known_job_that_has_left_the_queue(client, fake_redis) -> None:
    """A job that is running or finished is reported, not treated as missing.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's job client.
    """
    fake_redis.set(QUEUED_JOB_ID, "task-a")

    response = client.get("/job-queue/", params={"job_id": QUEUED_JOB_ID})

    assert response.status_code == 200
    body = response.json()
    assert body["position_in_queue"] is None
    assert body["total_jobs_in_queue"] == 0
    assert body["status"]


def test_job_queue_without_a_job_id_reports_the_queue_length(client, fake_redis) -> None:
    """The no-argument form is unaffected by the change.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's job client.
    """
    fake_redis.rpush("vntyper_job_queue", "task-a", "task-b")

    response = client.get("/job-queue/")

    assert response.status_code == 200
    assert response.json()["total_jobs_in_queue"] == 2
