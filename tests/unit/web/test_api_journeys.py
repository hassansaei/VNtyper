"""End-to-end journeys through the web service's own endpoints.

Five separate changes altered behaviour at this API's boundary: a mandatory cohort
passphrase, unique cohort aliases, an upload-filename allowlist, a ceiling on upload
and request size, a required Redis credential, a generic job-failure message and
cleanup of the archive `/cohort-download/` builds. Each of those has a test that
pins the property it introduced. None of them asks whether the service still does
the thing it exists to do.

That is what this module is for. It drives the real routes through the in-process
`client` fixture and follows the two journeys a user actually makes:

* **A job.** Submit an alignment, optionally with its index, and get an identifier
  back; watch it through `/job-status/` and `/job-queue/`; collect the result from
  `/download/` once it exists, and be told so honestly while it does not.
* **A cohort.** Create one, join a job to it, read it back, download its members'
  results as one archive, and run a joint analysis over them.

The assertions are about outcomes rather than status codes: the bytes that were
uploaded, under the name that was sent, in the directory the worker will look in;
the keyword arguments the Celery task was handed; the contents of the archive that
came back. A journey test that only checked for 200 would have passed against every
one of the defects this effort fixed.

The negative paths are asserted here too. They are the journeys' boundaries - the
conditions under which a journey must *not* start - and stating them alongside the
happy path is what stops a future regression from being "fixed" by loosening the
guard rather than by repairing the journey. The rules themselves are owned by
`test_cohort_authz.py`, `test_upload_safety.py` and `test_upload_limits.py`.

Nothing here runs: the Celery task objects are the `MagicMock`s the `web_app`
fixture installs, and `AsyncResult` is replaced where a journey needs a job to have
reached a particular state. `docker` is put on `sys.path` by
`tests/unit/web/conftest.py`, which pytest imports before this module.
"""

import io
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

pytestmark = pytest.mark.unit

from app.cohorts import COHORT_KEY_PREFIX  # noqa: E402
from app.version import API_VERSION  # noqa: E402

PASSPHRASE = "correct-horse-battery-staple"

# The endpoint never parses an upload, it only stores it, so recognisable bytes are
# more useful here than a valid BGZF stream: they make "the file that arrived is the
# file that was sent" an assertion rather than an assumption.
BAM_NAME = "sample.bam"
BAM_BYTES = b"alignment-bytes-for-sample.bam"

# `tasks.run_vntyper_job` looks for its index at exactly `<bam_path>.bai`, so an
# index uploaded as `sample.bam.bai` is the case where the worker finds one.
BAI_NAME = "sample.bam.bai"
BAI_BYTES = b"index-bytes-for-sample.bam.bai"

# What a finished job leaves in the output directory for `/download/` to serve.
RESULT_BYTES = b"zipped-vntyper-results"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _upload(name: str = BAM_NAME, payload: bytes = BAM_BYTES) -> dict[str, tuple[str, bytes, str]]:
    """Build the multipart file part for an alignment upload.

    Args:
        name: The filename to send, exactly as a client would.
        payload: The bytes to send under it.

    Returns:
        dict: A `files=` mapping for the test client.
    """
    return {"bam_file": (name, payload, "application/octet-stream")}


def _submit_job(client, **form: str) -> tuple[str, dict]:
    """Submit a valid job and return its identifier and the response body.

    Args:
        client: TestClient fixture from conftest.
        **form: Extra form fields, e.g. `cohort_id` and `passphrase`.

    Returns:
        tuple: The new job's identifier and the decoded response body.
    """
    response = client.post("/run-job/", files=_upload(), data=form)
    assert response.status_code == 200, response.text
    body = response.json()
    return body["job_id"], body


def _report_celery_status(monkeypatch: pytest.MonkeyPatch, web_app, status: str) -> None:
    """Make every Celery status lookup in the app report one state.

    Args:
        monkeypatch: Standard pytest fixture; restores `AsyncResult` at teardown.
        web_app: The patched `app.main` module.
        status: The Celery state to report, e.g. `"PENDING"` or `"SUCCESS"`.
    """
    monkeypatch.setattr(web_app, "AsyncResult", lambda task_id: SimpleNamespace(status=status, info=None))


def _create_cohort(client, alias: str = "family-a", passphrase: str = PASSPHRASE) -> str:
    """Create a cohort and return its identifier.

    Args:
        client: TestClient fixture from conftest.
        alias: The label to attach to the cohort.
        passphrase: The credential that will guard it.

    Returns:
        str: The new cohort's identifier.
    """
    response = client.post("/create-cohort/", data={"alias": alias, "passphrase": passphrase})
    assert response.status_code == 200, response.text
    return response.json()["cohort_id"]


# ---------------------------------------------------------------------------
# Journey 1: a job, from submission to download
# ---------------------------------------------------------------------------


def test_job_journey_from_submission_to_download(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A submitted alignment is stored, queued, tracked and finally downloadable.

    Walks the whole path a single job takes through the API and checks what each
    step actually produced, not merely that it answered.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        fake_redis: The store backing the app's Redis clients.
        monkeypatch: Standard pytest fixture, used to report Celery states.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    job_id, body = _submit_job(client)
    assert body["message"] == "Job submitted"
    assert body["cohort_id"] is None

    # The upload landed in the job's own input directory, under the name that was
    # sent and with the bytes that were sent.
    stored_bam = tmp_path / "handoff" / job_id / BAM_NAME
    assert stored_bam.read_bytes() == BAM_BYTES

    # ... and the worker was handed that path, plus the output directory it is
    # expected to write into. Nothing else can tell it where the upload went.
    web_app.run_vntyper_job.delay.assert_called_once()
    handed_over = web_app.run_vntyper_job.delay.call_args.kwargs
    assert handed_over["bam_path"] == str(stored_bam)
    assert handed_over["output_dir"] == str(tmp_path / "output" / job_id)
    assert handed_over["reference_assembly"] == "hg38"
    assert handed_over["advntr_mode"] is False
    assert handed_over["cohort_key"] is None
    assert handed_over["email"] is None

    # The job is findable afterwards: the identifier the caller was given resolves
    # to the task, and the task is in the queue.
    assert fake_redis.get(job_id) == "task-1"
    assert fake_redis.lrange("vntyper_job_queue", 0, -1) == ["task-1"]

    # While it waits, both tracking endpoints agree that it is waiting.
    _report_celery_status(monkeypatch, web_app, "PENDING")
    status = client.get(f"/job-status/{job_id}/")
    assert status.status_code == 200
    assert status.json() == {"job_id": job_id, "status": "pending", "error": None}

    queued = client.get("/job-queue/", params={"job_id": job_id})
    assert queued.status_code == 200
    assert queued.json()["position_in_queue"] == 1
    assert queued.json()["total_jobs_in_queue"] == 1

    # There is nothing to download yet, and the service says so rather than
    # serving an empty or partial file.
    not_ready = client.get(f"/download/{job_id}/")
    assert not_ready.status_code == 404
    assert not_ready.json() == {"detail": "File not found"}

    # Once the worker has finished and written the result archive, the same two
    # endpoints report completion and hand the archive over intact.
    (tmp_path / "output" / f"{job_id}.zip").write_bytes(RESULT_BYTES)
    _report_celery_status(monkeypatch, web_app, "SUCCESS")

    finished = client.get(f"/job-status/{job_id}/")
    assert finished.status_code == 200
    assert finished.json()["status"] == "completed"

    download = client.get(f"/download/{job_id}/")
    assert download.status_code == 200
    assert download.headers["content-type"] == "application/zip"
    assert f"{job_id}.zip" in download.headers["content-disposition"]
    assert download.content == RESULT_BYTES


def test_job_journey_stores_an_uploaded_index_beside_its_alignment(client, web_app, tmp_path: Path) -> None:
    """A submission carrying a `.bai` keeps both files, under the names sent.

    The worker looks for the index at `<bam_path>.bai` and builds one with samtools
    when it is absent, so an index stored anywhere else is silently discarded work.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    response = client.post(
        "/run-job/",
        files={
            "bam_file": (BAM_NAME, BAM_BYTES, "application/octet-stream"),
            "bai_file": (BAI_NAME, BAI_BYTES, "application/octet-stream"),
        },
    )

    assert response.status_code == 200, response.text
    job_id = response.json()["job_id"]

    job_input_dir = tmp_path / "handoff" / job_id
    assert sorted(path.name for path in job_input_dir.iterdir()) == [BAM_NAME, BAI_NAME]
    assert (job_input_dir / BAM_NAME).read_bytes() == BAM_BYTES
    assert (job_input_dir / BAI_NAME).read_bytes() == BAI_BYTES

    bam_path = web_app.run_vntyper_job.delay.call_args.kwargs["bam_path"]
    assert Path(f"{bam_path}.bai").read_bytes() == BAI_BYTES


def test_job_journey_routes_an_advntr_submission_to_the_long_queue(client, web_app, tmp_path: Path) -> None:
    """The optional adVNTR path still reaches its own queue with the same paths.

    adVNTR jobs are enqueued through `apply_async` rather than `delay` so they land
    on a separate worker queue; that is a second, otherwise untested submission
    path through the same upload handling.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    job_id, _ = _submit_job(client, advntr_mode="true")

    web_app.run_vntyper_job.delay.assert_not_called()
    web_app.run_vntyper_job.apply_async.assert_called_once()
    call = web_app.run_vntyper_job.apply_async.call_args
    assert call.kwargs["queue"] == "vntyper_long_queue"
    assert call.kwargs["kwargs"]["advntr_mode"] is True
    assert call.kwargs["kwargs"]["bam_path"] == str(tmp_path / "handoff" / job_id / BAM_NAME)


# ---------------------------------------------------------------------------
# Journey 2: a cohort, from creation to joint analysis
# ---------------------------------------------------------------------------


def test_cohort_journey_from_creation_to_joint_analysis(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A cohort can be created, joined, read, downloaded and analysed with one credential.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery tasks.
        fake_redis: The store backing the app's Redis clients.
        monkeypatch: Standard pytest fixture, used to report Celery states.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    credentials = {"cohort_id": _create_cohort(client), "passphrase": PASSPHRASE}
    cohort_id = credentials["cohort_id"]

    # Joining a job to the cohort is a write, so it needs the same credential.
    job_id, body = _submit_job(client, cohort_id=cohort_id, passphrase=PASSPHRASE)
    assert body["cohort_id"] == cohort_id

    # The membership is recorded, and the worker is told which cohort it is
    # contributing to, so its own bookkeeping lands in the right place.
    assert fake_redis.smembers(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs") == {job_id}
    assert web_app.run_vntyper_job.delay.call_args.kwargs["cohort_key"] == f"{COHORT_KEY_PREFIX}{cohort_id}"

    # Reading the cohort back reports its label and every member's state.
    _report_celery_status(monkeypatch, web_app, "PENDING")
    status = client.get("/cohort-status/", params=credentials)
    assert status.status_code == 200, status.text
    assert status.json() == {
        "cohort_id": cohort_id,
        "alias": "family-a",
        "jobs": [{"job_id": job_id, "status": "pending"}],
    }

    # Once the member has a result, the cohort download is one archive holding it,
    # readable and byte-identical to what the job produced.
    member_result = tmp_path / "output" / f"{job_id}.zip"
    member_result.write_bytes(RESULT_BYTES)

    download = client.get("/cohort-download/", params=credentials)
    assert download.status_code == 200, download.text
    assert download.headers["content-type"] == "application/zip"
    with zipfile.ZipFile(io.BytesIO(download.content)) as archive:
        assert archive.testzip() is None
        assert archive.namelist() == [f"{job_id}.zip"]
        assert archive.read(f"{job_id}.zip") == RESULT_BYTES

    # The joint analysis is enqueued over exactly the members' archives.
    analysis = client.post("/cohort-analysis/", data=credentials)
    assert analysis.status_code == 200, analysis.text
    analysis_job_id = analysis.json()["analysis_job_id"]
    assert analysis.json()["message"] == "Cohort analysis started"

    web_app.run_cohort_analysis_job.delay.assert_called_once()
    handed_over = web_app.run_cohort_analysis_job.delay.call_args.kwargs
    assert handed_over["cohort_id"] == cohort_id
    assert handed_over["zip_paths"] == [str(member_result)]
    assert handed_over["output_dir"] == str(tmp_path / "output" / analysis_job_id)

    # The analysis is trackable like any other job ...
    assert fake_redis.get(analysis_job_id) == "task-1"

    # ... and requesting it did not consume the members' results. What the task
    # itself may delete is pinned in `test_cohort_analysis_task.py`; the point
    # here is that a member is still downloadable after the cohort was analysed.
    assert member_result.read_bytes() == RESULT_BYTES
    assert client.get(f"/download/{job_id}/").content == RESULT_BYTES


def test_cohort_journey_gathers_several_members_into_one_archive(
    client, web_app, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A cohort with more than one finished job downloads as one archive of all of them.

    The single-member case cannot distinguish "collects the cohort" from "collects
    the first job it finds".

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        monkeypatch: Standard pytest fixture, used to report Celery states.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    credentials = {"cohort_id": _create_cohort(client), "passphrase": PASSPHRASE}

    results = {}
    for _ in range(2):
        job_id, _ = _submit_job(client, cohort_id=credentials["cohort_id"], passphrase=PASSPHRASE)
        payload = f"results-for-{job_id}".encode()
        (tmp_path / "output" / f"{job_id}.zip").write_bytes(payload)
        results[f"{job_id}.zip"] = payload

    _report_celery_status(monkeypatch, web_app, "SUCCESS")

    status = client.get("/cohort-status/", params=credentials)
    assert status.status_code == 200, status.text
    assert {job["status"] for job in status.json()["jobs"]} == {"completed"}
    assert len(status.json()["jobs"]) == 2

    download = client.get("/cohort-download/", params=credentials)
    assert download.status_code == 200, download.text
    with zipfile.ZipFile(io.BytesIO(download.content)) as archive:
        assert archive.testzip() is None
        assert sorted(archive.namelist()) == sorted(results)
        for name, payload in results.items():
            assert archive.read(name) == payload


# ---------------------------------------------------------------------------
# Journey 3: the boundaries - what must keep failing, and how
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("route", ["/cohort-status/", "/cohort-download/"])
def test_a_cohort_is_not_readable_with_the_wrong_passphrase(client, route: str) -> None:
    """The wrong credential is refused rather than answered.

    Args:
        client: TestClient fixture from conftest.
        route: The cohort read endpoint under test.
    """
    cohort_id = _create_cohort(client)

    response = client.get(route, params={"cohort_id": cohort_id, "passphrase": "not-the-passphrase"})

    assert response.status_code == 401
    assert response.json()["detail"] == "Incorrect passphrase"


def test_a_cohort_that_does_not_exist_is_reported_as_missing(client) -> None:
    """An unknown cohort is a 404, whatever passphrase accompanies it.

    Args:
        client: TestClient fixture from conftest.
    """
    response = client.get("/cohort-status/", params={"cohort_id": "no-such-cohort", "passphrase": PASSPHRASE})

    assert response.status_code == 404
    assert response.json()["detail"] == "Cohort ID not found"


def test_a_cohort_analysis_is_refused_without_the_passphrase(client) -> None:
    """The analysis route is a cohort read too, and is authorized on the same terms.

    Args:
        client: TestClient fixture from conftest.
    """
    cohort_id = _create_cohort(client)

    response = client.post("/cohort-analysis/", data={"cohort_id": cohort_id})

    assert response.status_code == 401
    assert response.json()["detail"] == "Passphrase required for this cohort"


def test_an_upload_the_service_does_not_accept_is_refused_before_anything_is_stored(
    client, web_app, tmp_path: Path
) -> None:
    """A disallowed filename is a 400, and leaves no directory behind.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
    """
    response = client.post("/run-job/", files=_upload("sample.txt", b"not an alignment"))

    assert response.status_code == 400
    assert "sample.txt" not in response.json()["detail"]
    assert list((tmp_path / "input").iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()


# Well-formed, and belonging to no job. Well-formed matters: the routes refuse an
# identifier they could not have issued before they look anything up, so an
# obviously malformed value would answer 404 without reaching the store at all.
# That path is `test_job_identifiers.py`'s subject; this one is about a real
# identifier for a job that is not there.
UNKNOWN_JOB_ID = "6f1e0c74-9a3b-4d52-8e17-2b95a4c60d38"


@pytest.mark.parametrize(
    ("path", "params", "detail"),
    [
        (f"/job-status/{UNKNOWN_JOB_ID}/", None, "Job ID not found"),
        (f"/download/{UNKNOWN_JOB_ID}/", None, "File not found"),
        ("/job-queue/", {"job_id": UNKNOWN_JOB_ID}, "Job ID not found"),
    ],
)
def test_an_unknown_job_id_is_reported_as_missing_not_as_a_server_error(
    client, path: str, params: dict | None, detail: str
) -> None:
    """Every job route answers an unrecognised identifier with 404.

    A 500 here would tell a caller polling a mistyped identifier that the service
    is broken, and would hide a real fault behind an ordinary one.

    Args:
        client: TestClient fixture from conftest.
        path: The route under test.
        params: Query parameters, if the route takes the identifier that way.
        detail: The message the route is expected to return.
    """
    response = client.get(path, params=params)

    assert response.status_code == 404
    assert response.json()["detail"] == detail


# ---------------------------------------------------------------------------
# Journey 4: the service still describes itself
# ---------------------------------------------------------------------------


def test_health_and_version_still_answer_with_their_real_payloads(client) -> None:
    """The two unauthenticated metadata routes respond, and say something true.

    Args:
        client: TestClient fixture from conftest.
    """
    health = client.get("/health/")
    assert health.status_code == 200
    assert health.json() == {"status": "ok"}

    version = client.get("/version/")
    assert version.status_code == 200
    assert version.json()["api_version"] == API_VERSION
    assert "tool_version" in version.json()


def test_usage_statistics_aggregate_the_recorded_jobs(client, fake_redis) -> None:
    """Usage statistics count jobs, distinct users and states from the usage store.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's usage client.
    """
    fake_redis.hset("usage:job-a", mapping={"user_hash": "user-1", "status": "completed"})
    fake_redis.hset("usage:job-b", mapping={"user_hash": "user-1", "status": "failed"})
    fake_redis.hset("usage:job-c", mapping={"user_hash": "user-2", "status": "completed"})

    response = client.get("/usage-statistics/")

    assert response.status_code == 200, response.text
    assert response.json() == {
        "total_jobs": 3,
        "unique_users": 2,
        "job_statuses": {"completed": 2, "failed": 1},
    }
