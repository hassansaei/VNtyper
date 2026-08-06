"""What a job submission leaves on the shared volume when it is not accepted.

`/run-job/` writes two client-supplied files into a per-job directory on the
volume every job shares. A submission can still be refused after that point --
for a malformed email address, or for a cohort the caller is not authorized to
join -- and a refusal that has already written to the volume has to take what it
wrote with it. Otherwise a caller who is never given a job identifier still
leaves an alignment behind, addressable by nobody and removed by nothing: the
worker's own cleanup only runs for jobs that were actually enqueued.

Two properties are pinned here:

1. A submission refused for any reason leaves neither the upload nor the
   directories that were created to hold it.
2. An accepted submission is unaffected -- the same directories are created and
   kept, and the job is enqueued.

The assertions look at the filesystem rather than at status codes, because a
status code cannot distinguish "refused" from "refused and cleaned up".

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.job_workspace` is importable here.
"""

from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

from app.job_workspace import job_workspace  # noqa: E402

BAM_NAME = "sample.bam"
BAM_BYTES = b"alignment-bytes"
PASSPHRASE = "correct-horse-battery-staple"


def _submit(client, **form: str):
    """Post a valid alignment with the given extra form fields.

    Args:
        client: TestClient fixture from conftest.
        **form: Extra form fields to send alongside the upload.

    Returns:
        httpx.Response: The service's response.
    """
    return client.post(
        "/run-job/",
        files={"bam_file": (BAM_NAME, BAM_BYTES, "application/octet-stream")},
        data=form,
    )


def _tree(root: Path) -> list[Path]:
    """List everything under a directory, files and directories alike.

    Args:
        root: The directory to walk.

    Returns:
        list[Path]: Every entry found beneath it.
    """
    return sorted(root.rglob("*"))


# ---------------------------------------------------------------------------
# The workspace helper on its own, with no app in the way.
# ---------------------------------------------------------------------------


def test_workspace_creates_both_directories(tmp_path: Path) -> None:
    """A job gets an input and an output directory named after it.

    Args:
        tmp_path: Scratch directory standing in for the volume.
    """
    with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1") as (job_input, job_output):
        assert Path(job_input) == tmp_path / "input" / "job-1"
        assert Path(job_output) == tmp_path / "output" / "job-1"
        assert Path(job_input).is_dir()
        assert Path(job_output).is_dir()


def test_workspace_keeps_its_directories_when_the_submission_succeeds(tmp_path: Path) -> None:
    """Nothing is removed on the way out of a successful submission.

    Args:
        tmp_path: Scratch directory standing in for the volume.
    """
    with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1") as (job_input, _):
        (Path(job_input) / BAM_NAME).write_bytes(BAM_BYTES)

    assert (tmp_path / "input" / "job-1" / BAM_NAME).read_bytes() == BAM_BYTES
    assert (tmp_path / "output" / "job-1").is_dir()


def test_workspace_removes_everything_it_created_when_the_body_raises(tmp_path: Path) -> None:
    """A failure part-way through takes the upload and both directories with it.

    Args:
        tmp_path: Scratch directory standing in for the volume.
    """
    # noqa: SIM117 - kept nested on purpose. The point of the test is that the
    # workspace's teardown runs *while* the exception propagates, so pytest.raises has
    # to visibly wrap the whole `with job_workspace(...)` block, teardown included.
    with pytest.raises(RuntimeError, match="refused"):  # noqa: SIM117
        with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1") as (job_input, _):
            (Path(job_input) / BAM_NAME).write_bytes(BAM_BYTES)
            msg = "the submission was refused"
            raise RuntimeError(msg)

    assert _tree(tmp_path / "input") == []
    assert _tree(tmp_path / "output") == []


# ---------------------------------------------------------------------------
# Endpoint level: every refusal path has to reclaim what it wrote.
# ---------------------------------------------------------------------------


def test_a_submission_refused_for_its_email_leaves_nothing_behind(client, web_app, tmp_path: Path) -> None:
    """A malformed email address is a 400 that costs the volume nothing.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the fixture configured as the job tree.
    """
    response = _submit(client, email="not-an-email-address")

    assert response.status_code == 400
    assert "job_id" not in response.json()
    assert _tree(tmp_path / "input") == []
    assert _tree(tmp_path / "output") == []
    web_app.run_vntyper_job.delay.assert_not_called()
    web_app.run_vntyper_job.apply_async.assert_not_called()


def test_a_submission_refused_for_its_cohort_leaves_nothing_behind(client, web_app, tmp_path: Path) -> None:
    """An alias with no identifier beside it is a 400 that costs the volume nothing.

    An alias is a label rather than an identifier, so this submission cannot name
    a cohort and is refused -- after the point at which the upload used to have
    been written.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the fixture configured as the job tree.
    """
    response = _submit(client, alias="family-a", passphrase=PASSPHRASE)

    assert response.status_code == 400
    assert _tree(tmp_path / "input") == []
    assert _tree(tmp_path / "output") == []
    web_app.run_vntyper_job.delay.assert_not_called()


def test_a_submission_refused_for_its_passphrase_leaves_nothing_behind(client, tmp_path: Path) -> None:
    """Joining a cohort with the wrong credential is refused and cleaned up.

    Args:
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the job tree.
    """
    created = client.post("/create-cohort/", data={"alias": "family-b", "passphrase": PASSPHRASE})
    assert created.status_code == 200, created.text
    cohort_id = created.json()["cohort_id"]

    response = _submit(client, cohort_id=cohort_id, passphrase="not-the-passphrase")

    assert response.status_code == 401
    assert _tree(tmp_path / "input") == []
    assert _tree(tmp_path / "output") == []


def test_a_submission_refused_after_the_upload_leaves_nothing_behind(
    client, web_app, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A failure while the job is being enqueued is cleaned up too.

    The email and cohort refusals above are the two the endpoint decides for
    itself. This one stands in for everything else that can go wrong once the
    upload is on disk, so the guarantee does not depend on enumerating causes.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        monkeypatch: Standard pytest fixture; restores the task at teardown.
        tmp_path: The directory the fixture configured as the job tree.
    """
    monkeypatch.setattr(web_app.run_vntyper_job, "delay", _raise_enqueue_failure)

    with pytest.raises(RuntimeError, match="the broker is unreachable"):
        _submit(client)

    assert _tree(tmp_path / "input") == []
    assert _tree(tmp_path / "output") == []


def _raise_enqueue_failure(**_kwargs) -> None:
    """Stand in for a broker that refuses the job after the upload is stored.

    Raises:
        RuntimeError: Always.
    """
    msg = "the broker is unreachable"
    raise RuntimeError(msg)


def test_an_accepted_submission_still_keeps_its_directories(client, web_app, tmp_path: Path) -> None:
    """The happy path is unchanged: the upload stays and the job is enqueued.

    Without this, every test above would also pass against an endpoint that
    deleted the job directory unconditionally.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the fixture configured as the job tree.
    """
    response = _submit(client)

    assert response.status_code == 200, response.text
    job_id = response.json()["job_id"]
    assert (tmp_path / "input" / job_id / BAM_NAME).read_bytes() == BAM_BYTES
    assert (tmp_path / "output" / job_id).is_dir()
    web_app.run_vntyper_job.delay.assert_called_once()
