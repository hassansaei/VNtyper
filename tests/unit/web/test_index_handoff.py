"""What happens to an alignment index between the endpoint and the worker.

`/run-job/` accepts an optional index part and stores it under the name the
client sent -- the documented input formats include `sample.bam.bai`,
`sample.bai` and, for CRAM, `sample.cram.crai`. The worker then has to find that
index, and afterwards has to remove it: the input directory holds
patient-derived data on a volume every job shares, and nothing else ever cleans
it up.

The two halves have to agree on one path. When they do not, an uploaded index is
ignored -- the worker rebuilds one it was already given -- and the file the
client sent stays on the volume for good. So the endpoint tells the worker where
it put the index, rather than the worker guessing, and the worker's lookup and
its cleanup both use what it was told.

The assertions are about the filesystem and the arguments handed over, not about
status codes: an index that is accepted, stored and then never mentioned again
produces exactly the same 200 as one that is used.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module.
"""

from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

BAM_BYTES = b"alignment-bytes"
INDEX_BYTES = b"index-bytes"


# ---------------------------------------------------------------------------
# The endpoint tells the worker where the index went
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("alignment_name", "index_name"),
    [
        ("sample.bam", "sample.bam.bai"),
        ("sample.bam", "sample.bai"),
        ("sample.cram", "sample.cram.crai"),
        ("sample.cram", "sample.crai"),
    ],
)
def test_the_worker_is_told_where_the_uploaded_index_was_stored(
    client, web_app, tmp_path: Path, alignment_name: str, index_name: str
) -> None:
    """Every accepted index name is handed over, not only the conventional one.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the fixture configured as the job tree.
        alignment_name: The alignment filename the client sends.
        index_name: The index filename the client sends alongside it.
    """
    response = client.post(
        "/run-job/",
        files={
            "bam_file": (alignment_name, BAM_BYTES, "application/octet-stream"),
            "bai_file": (index_name, INDEX_BYTES, "application/octet-stream"),
        },
    )

    assert response.status_code == 200, response.text
    job_input_dir = tmp_path / "input" / response.json()["job_id"]
    handed_over = web_app.run_vntyper_job.delay.call_args.kwargs
    assert handed_over["index_path"] == str(job_input_dir / index_name)
    assert (job_input_dir / index_name).read_bytes() == INDEX_BYTES


def test_a_submission_without_an_index_hands_over_none(client, web_app) -> None:
    """No index means no index path; the worker builds one for itself.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
    """
    response = client.post("/run-job/", files={"bam_file": ("sample.bam", BAM_BYTES, "application/octet-stream")})

    assert response.status_code == 200, response.text
    assert web_app.run_vntyper_job.delay.call_args.kwargs["index_path"] is None


def test_the_long_queue_submission_hands_over_the_index_too(client, web_app, tmp_path: Path) -> None:
    """The adVNTR path enqueues through `apply_async` and must carry the same paths.

    Args:
        client: TestClient fixture from conftest.
        web_app: The patched `app.main` module, holding the mocked Celery task.
        tmp_path: The directory the fixture configured as the job tree.
    """
    response = client.post(
        "/run-job/",
        files={
            "bam_file": ("sample.bam", BAM_BYTES, "application/octet-stream"),
            "bai_file": ("sample.bai", INDEX_BYTES, "application/octet-stream"),
        },
        data={"advntr_mode": "true"},
    )

    assert response.status_code == 200, response.text
    job_input_dir = tmp_path / "input" / response.json()["job_id"]
    handed_over = web_app.run_vntyper_job.apply_async.call_args.kwargs["kwargs"]
    assert handed_over["index_path"] == str(job_input_dir / "sample.bai")


# ---------------------------------------------------------------------------
# The worker uses what it was told, and clears the input directory
# ---------------------------------------------------------------------------


@pytest.fixture
def vntyper_task(monkeypatch: pytest.MonkeyPatch, fake_redis):
    """Run `run_vntyper_job`'s body with Redis and every subprocess neutralised.

    `subprocess.run` is replaced with a recorder that also does what `samtools
    index` would do, so "was an index built?" is an assertion about the commands
    the task issued rather than about a file that appeared by magic.

    Args:
        monkeypatch: Standard pytest fixture; restores every patch at teardown.
        fake_redis: In-process Redis stand-in from conftest.

    Returns:
        tuple: A callable invoking the task body, and the list of commands run.
    """
    from app import tasks

    for attr in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, attr, fake_redis)

    commands: list[list[str]] = []

    def _record(command, *args, **kwargs):
        """Record a command, standing in for `samtools index` where relevant.

        Args:
            command: The argument vector the task asked for.
            *args: Ignored.
            **kwargs: Ignored.
        """
        commands.append(list(command))
        if command[:2] == ["samtools", "index"]:
            Path(f"{command[2]}.bai").write_bytes(b"generated-index")

    monkeypatch.setattr(tasks.subprocess, "run", _record)

    def _invoke(**kwargs) -> None:
        """Call the task body as a worker would.

        Args:
            **kwargs: Arguments forwarded to `run_vntyper_job`.
        """
        tasks.run_vntyper_job.push_request(id="task-1")
        try:
            tasks.run_vntyper_job.run(**kwargs)
        finally:
            tasks.run_vntyper_job.pop_request()

    return _invoke, commands


def _job_input(tmp_path: Path, alignment_name: str, index_name: str | None) -> tuple[Path, Path | None]:
    """Lay out a job input directory the way the endpoint would have.

    Args:
        tmp_path: Scratch directory standing in for the input tree.
        alignment_name: The stored alignment's filename.
        index_name: The stored index's filename, or None for no index.

    Returns:
        tuple: The alignment path, and the index path or None.
    """
    job_input_dir = tmp_path / "input" / "job-1"
    job_input_dir.mkdir(parents=True)
    alignment = job_input_dir / alignment_name
    alignment.write_bytes(BAM_BYTES)
    if index_name is None:
        return alignment, None
    index = job_input_dir / index_name
    index.write_bytes(INDEX_BYTES)
    return alignment, index


def _run(invoke, alignment: Path, tmp_path: Path, index: Path | None) -> None:
    """Invoke the task for one alignment.

    Args:
        invoke: The callable from the `vntyper_task` fixture.
        alignment: The stored alignment path.
        tmp_path: Scratch directory standing in for the job tree.
        index: The stored index path, or None.
    """
    invoke(
        bam_path=str(alignment),
        output_dir=str(tmp_path / "output" / "job-1"),
        thread=1,
        reference_assembly="hg38",
        fast_mode=False,
        keep_intermediates=False,
        archive_results=False,
        index_path=None if index is None else str(index),
    )


@pytest.mark.parametrize(
    ("alignment_name", "index_name"),
    [
        ("sample.bam", "sample.bam.bai"),
        ("sample.bam", "sample.bai"),
        ("sample.cram", "sample.cram.crai"),
    ],
)
def test_an_index_that_was_supplied_is_not_rebuilt(
    vntyper_task, tmp_path: Path, alignment_name: str, index_name: str
) -> None:
    """An index the submission carried is used as it is.

    Rebuilding it is not merely wasted work: the rebuilt index lands under a
    different name, so the uploaded one is left on the volume afterwards.

    Args:
        vntyper_task: The task fixture above.
        tmp_path: Scratch directory standing in for the job tree.
        alignment_name: The stored alignment's filename.
        index_name: The stored index's filename.
    """
    invoke, commands = vntyper_task
    alignment, index = _job_input(tmp_path, alignment_name, index_name)

    _run(invoke, alignment, tmp_path, index)

    assert not any(command[:2] == ["samtools", "index"] for command in commands)


@pytest.mark.parametrize(
    ("alignment_name", "index_name"),
    [
        ("sample.bam", "sample.bam.bai"),
        ("sample.bam", "sample.bai"),
        ("sample.cram", "sample.cram.crai"),
        ("sample.cram", None),
        ("sample.bam", None),
    ],
)
def test_the_job_input_directory_is_empty_afterwards(
    vntyper_task, tmp_path: Path, alignment_name: str, index_name: str | None
) -> None:
    """No patient-derived file survives the job that was given it.

    Args:
        vntyper_task: The task fixture above.
        tmp_path: Scratch directory standing in for the job tree.
        alignment_name: The stored alignment's filename.
        index_name: The stored index's filename, or None for no index.
    """
    invoke, _ = vntyper_task
    alignment, index = _job_input(tmp_path, alignment_name, index_name)

    _run(invoke, alignment, tmp_path, index)

    assert not alignment.exists()
    assert not (tmp_path / "input" / "job-1").exists()


def test_an_index_is_still_built_when_the_submission_carried_none(vntyper_task, tmp_path: Path) -> None:
    """The fallback is unchanged: with no index supplied, the worker makes one.

    Without this, the tests above would also pass against a worker that had
    stopped indexing altogether.

    Args:
        vntyper_task: The task fixture above.
        tmp_path: Scratch directory standing in for the job tree.
    """
    invoke, commands = vntyper_task
    alignment, _ = _job_input(tmp_path, "sample.bam", None)

    _run(invoke, alignment, tmp_path, None)

    assert ["samtools", "index", str(alignment)] in commands
