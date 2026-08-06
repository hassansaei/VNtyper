"""The shape a job identifier has to have before the service will use it.

Three routes take a job identifier from the caller -- `/job-status/{job_id}/`,
`/download/{job_id}/` and `/job-queue/?job_id=` -- and each of them uses it
directly: as a Redis key, and as part of a path in the output directory. The
service issues those identifiers itself, as version-4 UUIDs, so anything that is
not one cannot name a job. Accepting it anyway means addressing arbitrary keys
in the job database, where a name that happens to hold something other than a
string answers with a server error rather than "no such job".

So the format is checked before the value is used at all, in one place shared by
the three routes, and an identifier that is not one is answered exactly as an
identifier that names nothing: 404. Not 400 -- the distinction would tell a
caller which identifiers exist -- and not 500, which would report the service as
broken when the truth is that there is no such job.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module.
"""

from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

import pytest

pytestmark = pytest.mark.unit

from app.identifiers import is_job_id  # noqa: E402


def _report_celery_status(monkeypatch: pytest.MonkeyPatch, web_app, status: str) -> None:
    """Make every Celery status lookup in the app report one state.

    `/job-status/` builds an `AsyncResult`, which reaches for a result backend
    the unit tier does not have.

    Args:
        monkeypatch: Standard pytest fixture; restores `AsyncResult` at teardown.
        web_app: The patched `app.main` module.
        status: The Celery state to report.
    """
    monkeypatch.setattr(web_app, "AsyncResult", lambda task_id: SimpleNamespace(status=status, info=None))

# Values a caller could send that cannot have come from this service. The queue
# key is the interesting one: it lives in the same Redis database as the job
# mappings and holds a list, so reading it as a string is a type error.
NOT_JOB_IDS = [
    ("vntyper_job_queue", "the queue's own key, which holds a list rather than a string"),
    ("", "an empty identifier names nothing"),
    ("no-such-job", "plainly malformed"),
    ("../../etc/passwd", "path traversal: the identifier also names a file"),
    ("a/b", "a path separator has no place in a single key or path segment"),
    ("usage:job", "a key belonging to another part of the service"),
    ("*", "a glob, which a key scan would expand"),
    ("00000000-0000-0000-0000-00000000000", "35 hex digits: one short of a UUID"),
    ("00000000-0000-0000-0000-0000000000000", "37 hex digits: one over"),
    ("00000000000000000000000000000000", "unhyphenated: not the form this service issues"),
    ("{00000000-0000-0000-0000-000000000000}", "braced form"),
    ("urn:uuid:00000000-0000-0000-0000-000000000000", "URN form"),
    ("0000000g-0000-0000-0000-000000000000", "'g' is not a hex digit"),
    ("0000000 -0000-0000-0000-000000000000", "embedded space"),
    ("00000000-0000-0000-0000-000000000000\n", "trailing newline"),
    ("00000000-0000-0000-0000-000000000000\x00", "trailing NUL"),
]


@pytest.mark.parametrize("value,reason", NOT_JOB_IDS, ids=[repr(value) for value, _ in NOT_JOB_IDS])
def test_a_value_the_service_never_issued_is_not_a_job_id(value: str, reason: str) -> None:
    """Nothing outside the identifier format the service issues is accepted.

    Args:
        value: The candidate identifier.
        reason: Why it must be refused; reported on failure.
    """
    assert is_job_id(value) is False, reason


def test_an_identifier_this_service_issues_is_accepted() -> None:
    """The check must accept what the service itself hands out."""
    for _ in range(16):
        assert is_job_id(str(uuid4())) is True


def test_the_check_is_indifferent_to_hex_case() -> None:
    """An identifier echoed back in upper case is still the same identifier."""
    issued = str(uuid4())

    assert is_job_id(issued.upper()) is True


# ---------------------------------------------------------------------------
# Endpoint level: the check has to be wired into every route that takes one.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", ["vntyper_job_queue", "no-such-job", "*"])
def test_job_status_answers_an_unusable_identifier_with_404(client, fake_redis, value: str) -> None:
    """A value that cannot name a job is reported as missing, not as a fault.

    The queue key is seeded so it holds what it holds in a running service: a
    list. Reading it as a string is what turns a mistyped identifier into a
    server error.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's job client.
        value: The identifier under test.
    """
    fake_redis.rpush("vntyper_job_queue", "task-1")

    response = client.get(f"/job-status/{value}/")

    assert response.status_code == 404
    assert response.json()["detail"] == "Job ID not found"


@pytest.mark.parametrize("value", ["vntyper_job_queue", "no-such-job", "*"])
def test_job_queue_answers_an_unusable_identifier_with_404(client, fake_redis, value: str) -> None:
    """The queue route applies the same check to the identifier it is given.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's job client.
        value: The identifier under test.
    """
    fake_redis.rpush("vntyper_job_queue", "task-1")

    response = client.get("/job-queue/", params={"job_id": value})

    assert response.status_code == 404
    assert response.json()["detail"] == "Job ID not found"


@pytest.mark.parametrize("value", ["vntyper_job_queue", "no-such-job", "*"])
def test_download_answers_an_unusable_identifier_with_404(client, value: str) -> None:
    """The download route builds a path from the identifier and checks it first.

    Args:
        client: TestClient fixture from conftest.
        value: The identifier under test.
    """
    response = client.get(f"/download/{value}/")

    assert response.status_code == 404
    assert response.json()["detail"] == "File not found"


def test_download_does_not_serve_a_file_named_by_an_unusable_identifier(client, tmp_path: Path) -> None:
    """A 404 has to mean nothing was served, not merely that nothing was found.

    Args:
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the output tree.
    """
    (tmp_path / "output" / "vntyper_job_queue.zip").write_bytes(b"not-a-result")

    response = client.get("/download/vntyper_job_queue/")

    assert response.status_code == 404
    assert b"not-a-result" not in response.content


def test_the_identifiers_the_service_issues_are_still_accepted_everywhere(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A real job identifier reaches all three routes unchanged.

    Without this, every test above would also pass against routes that refused
    every identifier they were given.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module.
        fake_redis: The store backing the app's job client.
        monkeypatch: Standard pytest fixture, used to report a Celery state.
        tmp_path: The directory the fixture configured as the job tree.
    """
    _report_celery_status(monkeypatch, web_app, "PENDING")
    job_id = str(uuid4())
    fake_redis.set(job_id, "task-1")
    fake_redis.rpush("vntyper_job_queue", "task-1")
    (tmp_path / "output" / f"{job_id}.zip").write_bytes(b"result-bytes")

    assert client.get(f"/job-status/{job_id}/").status_code == 200
    assert client.get("/job-queue/", params={"job_id": job_id}).json()["position_in_queue"] == 1
    assert client.get(f"/download/{job_id}/").content == b"result-bytes"


def test_an_analysis_identifier_is_usable_on_the_routes_that_take_one(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A cohort analysis hands back an identifier the job routes accept.

    The analysis endpoint mints its own identifier, so it has to be minted in the
    form the check demands or the caller is handed something it cannot poll.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module.
        fake_redis: The store backing the app's clients.
        monkeypatch: Standard pytest fixture, used to report a Celery state.
        tmp_path: The directory the fixture configured as the job tree.
    """
    from app.cohorts import COHORT_KEY_PREFIX

    _report_celery_status(monkeypatch, web_app, "PENDING")
    passphrase = "correct-horse-battery-staple"
    created = client.post("/create-cohort/", data={"alias": "study", "passphrase": passphrase})
    cohort_id = created.json()["cohort_id"]
    member = str(uuid4())
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", member)
    (tmp_path / "output" / f"{member}.zip").write_bytes(b"member-result")

    analysis = client.post("/cohort-analysis/", data={"cohort_id": cohort_id, "passphrase": passphrase})
    assert analysis.status_code == 200, analysis.text
    analysis_job_id = analysis.json()["analysis_job_id"]

    assert is_job_id(analysis_job_id)
    assert client.get(f"/job-status/{analysis_job_id}/").status_code == 200
