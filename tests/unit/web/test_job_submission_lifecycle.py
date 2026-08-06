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
   kept, and the job is enqueued. That holds from the moment the task is
   accepted, not only at the end of the request: a worker owns the input
   directory from then on, so a later failure must keep it.

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


def test_workspace_keeps_its_directories_when_the_body_raises_after_the_commit(tmp_path: Path) -> None:
    """Once the task is accepted the directories belong to the worker, not to the request.

    The reclaim above is scoped to submissions that were *not* accepted. After the
    commit point, a failure in the request is no longer evidence that nothing is
    running: the worker holds the input path and is reading it.

    Args:
        tmp_path: Scratch directory standing in for the volume.
    """
    with pytest.raises(RuntimeError, match="after acceptance"):  # noqa: SIM117
        with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1") as workspace:
            (Path(workspace.input_dir) / BAM_NAME).write_bytes(BAM_BYTES)
            workspace.commit()
            msg = "a bookkeeping write failed after acceptance"
            raise RuntimeError(msg)

    assert (tmp_path / "input" / "job-1" / BAM_NAME).read_bytes() == BAM_BYTES
    assert (tmp_path / "output" / "job-1").is_dir()


def test_workspace_still_unpacks_as_the_pair_of_directories(tmp_path: Path) -> None:
    """The commit point is additive: the two directories are still what it yields."""
    with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1") as workspace:
        job_input, job_output = workspace
        assert job_input == workspace.input_dir == str(tmp_path / "input" / "job-1")
        assert job_output == workspace.output_dir == str(tmp_path / "output" / "job-1")


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


@pytest.mark.parametrize("write", ["set", "rpush"])
def test_a_failure_after_acceptance_does_not_delete_the_running_jobs_input(
    write: str, client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Every bookkeeping write after the enqueue is injected with a failure in turn.

    Celery has accepted the task by then, so a worker is already reading the BAM at
    the path it was handed. Reclaiming the directory here deletes a running job's
    input underneath it -- a worse outcome than the leftover the reclaim exists to
    prevent, and one that produces a wrong result rather than a wasted inode.

    Args:
        write: Name of the Redis method to fail.
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        fake_redis: The in-process Redis every client in `web_app` shares.
        monkeypatch: Standard pytest fixture; restores the method at teardown.
        tmp_path: The directory the fixture configured as the job tree.
    """

    def _fail(*_args, **_kwargs):
        msg = "the connection to Redis dropped"
        raise ConnectionError(msg)

    monkeypatch.setattr(fake_redis, write, _fail)

    with pytest.raises(ConnectionError, match="dropped"):
        _submit(client)

    web_app.run_vntyper_job.delay.assert_called_once()
    accepted_input = Path(web_app.run_vntyper_job.delay.call_args.kwargs["bam_path"])
    assert accepted_input.read_bytes() == BAM_BYTES, "the worker's input was deleted underneath it"
    assert Path(web_app.run_vntyper_job.delay.call_args.kwargs["output_dir"]).is_dir()


def test_an_enqueue_failure_leaves_no_ghost_member_in_the_cohort(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Cohort membership is recorded before the enqueue, so it has to be undone with it.

    The member is added first because the task needs the cohort key. If the broker
    then refuses the task, the identifier stays in the cohort's Set naming a job
    that does not exist and whose directories have just been reclaimed: every
    cohort report afterwards looks for results that will never arrive.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        fake_redis: The in-process Redis every client in `web_app` shares.
        monkeypatch: Standard pytest fixture; restores the task at teardown.
        tmp_path: The directory the fixture configured as the job tree.
    """
    created = client.post("/create-cohort/", data={"alias": "family-c", "passphrase": PASSPHRASE})
    assert created.status_code == 200, created.text
    cohort_id = created.json()["cohort_id"]
    members_key = f"cohort:{cohort_id}:jobs"

    monkeypatch.setattr(web_app.run_vntyper_job, "delay", _raise_enqueue_failure)

    with pytest.raises(RuntimeError, match="the broker is unreachable"):
        _submit(client, cohort_id=cohort_id, passphrase=PASSPHRASE)

    assert fake_redis.smembers(members_key) == set()
    assert _tree(tmp_path / "input") == []


def test_an_accepted_cohort_submission_is_recorded_as_a_member(client, fake_redis) -> None:
    """The compensation above must not fire on the happy path.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The in-process Redis every client in `web_app` shares.
    """
    created = client.post("/create-cohort/", data={"alias": "family-d", "passphrase": PASSPHRASE})
    cohort_id = created.json()["cohort_id"]

    response = _submit(client, cohort_id=cohort_id, passphrase=PASSPHRASE)

    assert response.status_code == 200, response.text
    assert fake_redis.smembers(f"cohort:{cohort_id}:jobs") == {response.json()["job_id"]}


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
