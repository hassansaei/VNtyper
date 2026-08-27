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

import hashlib
from contextlib import suppress
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


def test_workspace_reclaims_nested_contents_through_its_bound_descriptors(tmp_path: Path) -> None:
    """Refusal cleanup recursively clears captured directories without path traversal."""
    with pytest.raises(RuntimeError, match="refused"):  # noqa: SIM117
        with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1") as workspace:
            nested_input = Path(workspace.input_dir) / "nested"
            nested_output = Path(workspace.output_dir) / "stage" / "deeper"
            nested_input.mkdir()
            nested_output.mkdir(parents=True)
            (nested_input / BAM_NAME).write_bytes(BAM_BYTES)
            (nested_output / "partial.txt").write_text("partial", encoding="utf-8")
            raise RuntimeError("the submission was refused")

    assert _tree(tmp_path / "input") == []
    assert _tree(tmp_path / "output") == []


def test_job_workspace_recursive_cleanup_preserves_primary_while_closing_ancestors(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Nested descriptor closes cannot replace the recursive cleanup error."""
    from app import job_workspace as workspace_module

    root = tmp_path / "root"
    trigger = root / "nested" / "deeper" / "trigger"
    trigger.parent.mkdir(parents=True)
    trigger.write_bytes(b"data")
    real_open = workspace_module.os.open
    real_stat = workspace_module.os.stat
    real_close = workspace_module.os.close
    root_descriptor = real_open(root, workspace_module.os.O_RDONLY | workspace_module.os.O_DIRECTORY)
    opened: list[int] = []
    closed: list[int] = []

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def failing_stat(path, *args, **kwargs):
        if path == "trigger":
            raise OSError("primary recursive cleanup failure")
        return real_stat(path, *args, **kwargs)

    def failing_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)
        raise OSError(f"secondary close failure {descriptor}")

    monkeypatch.setattr(workspace_module.os, "open", recording_open)
    monkeypatch.setattr(workspace_module.os, "stat", failing_stat)
    monkeypatch.setattr(workspace_module.os, "close", failing_close)
    try:
        with pytest.raises(OSError, match="primary recursive cleanup failure"):
            workspace_module._clear_directory(root_descriptor)
    finally:
        real_close(root_descriptor)
        for descriptor in opened:
            with suppress(OSError):
                real_close(descriptor)

    assert closed == list(reversed(opened))


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


def test_handoff_refuses_replaced_names_and_reclaims_only_the_bound_inodes(tmp_path: Path) -> None:
    """A pre-handoff swap leaves replacements alone and no uploaded bytes behind."""
    input_root = tmp_path / "input"
    output_root = tmp_path / "output"
    displaced_input = input_root / "job-1-displaced"
    displaced_output = output_root / "job-1-displaced"

    with pytest.raises(RuntimeError, match="input directory identity changed before handoff"):  # noqa: SIM117
        with job_workspace(str(input_root), str(output_root), "job-1") as workspace:
            alignment = Path(workspace.input_dir) / BAM_NAME
            alignment.write_bytes(BAM_BYTES)
            alignment_metadata = alignment.stat()
            Path(workspace.input_dir).rename(displaced_input)
            Path(workspace.output_dir).rename(displaced_output)
            Path(workspace.input_dir).mkdir()
            Path(workspace.output_dir).mkdir()
            (Path(workspace.input_dir) / "replacement.txt").write_text("input replacement", encoding="utf-8")
            (Path(workspace.output_dir) / "replacement.txt").write_text("output replacement", encoding="utf-8")
            workspace.handoff(
                (alignment_metadata.st_dev, alignment_metadata.st_ino),
                None,
                hashlib.sha256(BAM_BYTES).hexdigest(),
                None,
            )

    assert (Path(input_root) / "job-1" / "replacement.txt").read_text(encoding="utf-8") == "input replacement"
    assert (Path(output_root) / "job-1" / "replacement.txt").read_text(encoding="utf-8") == "output replacement"
    assert list(displaced_input.iterdir()) == []
    assert list(displaced_output.iterdir()) == []


def test_workspace_close_attempts_both_descriptors_without_masking_the_body_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The body error stays primary when one workspace descriptor close fails."""
    from app import job_workspace as workspace_module

    real_close = workspace_module.os.close
    closed: list[int] = []
    output_descriptor: int | None = None
    input_descriptor: int | None = None

    def close_with_one_failure(descriptor: int) -> None:
        closed.append(descriptor)
        if descriptor == output_descriptor:
            real_close(descriptor)
            raise OSError("output close failed")
        real_close(descriptor)

    with pytest.raises(RuntimeError, match="body failed"):  # noqa: SIM117
        with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1") as workspace:
            output_descriptor = workspace.output_descriptor
            input_descriptor = workspace.input_descriptor
            monkeypatch.setattr(workspace_module.os, "close", close_with_one_failure)
            raise RuntimeError("body failed")

    assert closed == [output_descriptor, input_descriptor]


def test_workspace_binding_failure_attempts_all_closes_and_preserves_the_primary_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Setup cleanup neither leaks a bound directory nor masks its bind failure."""
    from app import job_workspace as workspace_module

    real_open = workspace_module.os.open
    real_fstat = workspace_module.os.fstat
    real_close = workspace_module.os.close
    owned_descriptors: list[int] = []
    owned_close_attempts: list[int] = []
    input_close_failed = False

    def recording_open(path, *args, **kwargs):
        descriptor = real_open(path, *args, **kwargs)
        if len(owned_descriptors) < 2:
            owned_descriptors.append(descriptor)
        return descriptor

    def fail_output_binding(descriptor: int):
        if len(owned_descriptors) == 2 and descriptor == owned_descriptors[1]:
            raise OSError("primary bind failure")
        return real_fstat(descriptor)

    def close_with_secondary_error(descriptor: int) -> None:
        nonlocal input_close_failed
        if descriptor in owned_descriptors and descriptor not in owned_close_attempts:
            owned_close_attempts.append(descriptor)
        real_close(descriptor)
        if owned_descriptors and descriptor == owned_descriptors[0] and not input_close_failed:
            input_close_failed = True
            raise OSError("secondary close failure")

    monkeypatch.setattr(workspace_module.os, "open", recording_open)
    monkeypatch.setattr(workspace_module.os, "fstat", fail_output_binding)
    monkeypatch.setattr(workspace_module.os, "close", close_with_secondary_error)
    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    try:
        with pytest.raises(RuntimeError, match="primary bind failure"):  # noqa: SIM117
            with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1"):
                raise AssertionError("binding failure must prevent the context from yielding")

        assert owned_close_attempts == owned_descriptors
        assert not input_dir.exists()
        assert not output_dir.exists()
    finally:
        for path in (input_dir, output_dir):
            if path.exists():
                path.rmdir()


def test_workspace_construction_fstat_failure_closes_and_reclaims_every_bound_directory(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The constructor's second identity read remains inside setup cleanup."""
    from app import job_workspace as workspace_module

    real_fstat = workspace_module.os.fstat
    real_close = workspace_module.os.close
    fstat_calls = 0
    opened: list[int] = []
    closed: list[int] = []
    real_open = workspace_module.os.open

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def fail_constructor_fstat(descriptor: int):
        nonlocal fstat_calls
        fstat_calls += 1
        if fstat_calls == 3:
            raise OSError("constructor fstat failure")
        return real_fstat(descriptor)

    def recording_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)

    monkeypatch.setattr(workspace_module.os, "open", recording_open)
    monkeypatch.setattr(workspace_module.os, "fstat", fail_constructor_fstat)
    monkeypatch.setattr(workspace_module.os, "close", recording_close)

    with pytest.raises(RuntimeError, match="constructor fstat failure"):  # noqa: SIM117
        with job_workspace(str(tmp_path / "input"), str(tmp_path / "output"), "job-1"):
            raise AssertionError("construction failure must prevent the context from yielding")

    assert closed == opened
    assert _tree(tmp_path / "input") == []
    assert _tree(tmp_path / "output") == []


def test_workspace_documentation_distinguishes_the_spool_from_the_output_volume() -> None:
    """The helper documentation names the two different trust surfaces."""
    from app import job_workspace as workspace_module

    assert workspace_module.__doc__ is not None
    assert "protected handoff spool" in workspace_module.__doc__
    assert "service-private result store" in workspace_module.__doc__


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
    assert (tmp_path / "handoff" / job_id / BAM_NAME).read_bytes() == BAM_BYTES
    assert (tmp_path / "output" / job_id).is_dir()
    web_app.run_vntyper_job.delay.assert_called_once()
