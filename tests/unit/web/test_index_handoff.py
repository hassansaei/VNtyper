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

Two further ways a patient-derived file survives that cleanup are pinned at the
bottom of this module. The index the worker builds for itself carries a name the
client never sent, so cleanup has to derive that name rather than only remove
what it was handed; and the cleanup block shares a `finally` with Redis
bookkeeping that can fail before it is reached.

The pipeline used to be a third way: it reconstructed `f"{in_bam}.bai"`, missed
an uploaded `sample.bai`, and built a second index beside the alignment. Since
#210 it resolves both names htslib resolves and builds any index it needs in the
run's output directory, so it no longer writes into the input tree at all --
pinned below as an absence.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module.
"""

import logging
import subprocess
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


def _index_the_way_the_pipeline_does(command: list[str]) -> None:
    """Reproduce the index `vntyper pipeline` builds for a non-fast BAM run.

    Read off `vntyper/scripts/fastq_bam_processing.py`, not assumed, and it calls
    the pipeline's own `resolve_bam_index` -- now in
    `vntyper/scripts/alignment_index.py` -- rather than restating its rules: on the
    non-fast BAM path the pipeline looks for `<file>.bai` and then `<stem>.bai`,
    the two names htslib itself resolves, and indexes only when neither exists.

    Until #210 it reconstructed `f"{in_bam}.bai"` alone and wrote the index it
    built beside the alignment. Both halves of that were wrong for a job: it
    never saw an index the client had uploaded as `sample.bai`, so it built a
    second one, and it built it inside the job's own input directory on the
    volume every job shares. The index now goes into the run's output directory,
    which is what this stand-in reproduces -- a recorder that still wrote beside
    the alignment would keep the tests below passing while describing behaviour
    that no longer exists.

    The CRAM branch of the same function extracts unmapped reads with
    `samtools view` instead and builds no index at all, so nothing is created
    here for `--cram`.

    Args:
        command: The argument vector the task asked for.
    """
    from vntyper.scripts.alignment_index import resolve_bam_index

    if "pipeline" not in command or "--bam" not in command or "--fast-mode" in command:
        return
    alignment = command[command.index("--bam") + 1]
    if resolve_bam_index(alignment) is not None:
        return
    output_dir = Path(command[command.index("-o") + 1])
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "output_input.bam.bai").write_bytes(b"pipeline-index")


@pytest.fixture
def vntyper_task(monkeypatch: pytest.MonkeyPatch, fake_redis):
    """Run `run_vntyper_job`'s body with Redis and every subprocess neutralised.

    `subprocess.run` is replaced with a recorder that also does what `samtools
    index` would do, so "was an index built?" is an assertion about the commands
    the task issued rather than about a file that appeared by magic.

    The stand-in for `vntyper pipeline` has a side effect too, for the same
    reason: a recorder that only records makes the input directory look tidier
    than the real thing ever leaves it, and the second index this module tests
    for would never appear.

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
            # The suffix mirrors real `samtools index`, which writes `.crai`
            # beside a CRAM and `.bai` beside a BAM -- verified by conversion
            # and indexing under #188. Do not collapse this back to one suffix:
            # a stand-in that always writes `.bai` encodes the very defect #188
            # fixed, and made the CRAM-without-an-index case below pass against
            # a worker that was leaving the real index on the shared volume.
            alignment = command[2]
            suffix = ".crai" if alignment.lower().endswith(".cram") else ".bai"
            Path(f"{alignment}{suffix}").write_bytes(b"generated-index")
            return
        _index_the_way_the_pipeline_does(command)

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


def _run(invoke, alignment: Path, tmp_path: Path, index: Path | None, *, cohort_key: str | None = None) -> None:
    """Invoke the task for one alignment.

    Args:
        invoke: The callable from the `vntyper_task` fixture.
        alignment: The stored alignment path.
        tmp_path: Scratch directory standing in for the job tree.
        index: The stored index path, or None.
        cohort_key: The cohort the job belongs to, or None for a lone job.
    """
    invoke(
        bam_path=str(alignment),
        output_dir=str(tmp_path / "output" / "job-1"),
        thread=1,
        reference_assembly="hg38",
        fast_mode=False,
        keep_intermediates=False,
        archive_results=False,
        cohort_key=cohort_key,
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

    "No file" includes the index the worker builds for itself when the submission
    carried none, which lands under a name the client never sent -- see the
    dedicated test below. The pipeline stand-in
    (`_index_the_way_the_pipeline_does`) writes into the output directory, as the
    pipeline now does, so nothing it produces is in scope for this assertion.

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


# ---------------------------------------------------------------------------
# The index the pipeline builds for itself is cleaned up as well
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("alignment", "expected"),
    [
        (
            "/data/input/job-1/sample.bam",
            (
                "/data/input/job-1/sample.bam.bai",
                "/data/input/job-1/sample.bam.crai",
                "/data/input/job-1/sample.bai",
                "/data/input/job-1/sample.crai",
            ),
        ),
        (
            "/data/input/job-1/sample.cram",
            (
                "/data/input/job-1/sample.cram.bai",
                "/data/input/job-1/sample.cram.crai",
                "/data/input/job-1/sample.bai",
                "/data/input/job-1/sample.crai",
            ),
        ),
        # An extensionless name collapses `<alignment>` and `<stem>` onto the
        # same two paths; cleanup must not then try to delete either twice.
        (
            "/data/input/job-1/sample",
            ("/data/input/job-1/sample.bai", "/data/input/job-1/sample.crai"),
        ),
    ],
)
def test_the_derived_index_names_are_exactly_the_four_forms(alignment: str, expected: tuple[str, ...]) -> None:
    """Both the `<alignment>.*` and the `<stem>.*` forms are named, and nothing else.

    Args:
        alignment: The stored alignment path.
        expected: The paths cleanup is entitled to remove for it.
    """
    from app.tasks import derived_index_paths

    assert derived_index_paths(alignment) == expected


@pytest.mark.parametrize("alignment", ["/data/input/job-1/sample.bam", "/data/input/job-1/SAMPLE.CRAM"])
def test_every_derived_index_name_stays_inside_the_jobs_own_directory(alignment: str) -> None:
    """The containment property, asserted rather than assumed.

    This is what makes widening the removal set safe: the names are built from
    the alignment and joined back onto its own directory, so cleanup cannot
    reach a sibling job's input directory however the alignment is named.

    Args:
        alignment: The stored alignment path.
    """
    from app.tasks import derived_index_paths

    job_input_dir = "/data/input/job-1"
    assert [str(Path(path).parent) for path in derived_index_paths(alignment)] == [job_input_dir] * 4


def test_the_pipeline_no_longer_builds_an_index_beside_the_submission(vntyper_task, tmp_path: Path) -> None:
    """An uploaded `sample.bai` is now honoured by the pipeline as well (#210).

    The endpoint and the worker both accept `sample.bai`, so the worker skips its
    preflight `samtools index`. The pipeline used not to accept it: it
    reconstructed the name as `f"{in_bam}.bai"`, did not find it, and indexed the
    alignment itself -- leaving a second index beside the submission under a name
    neither the client nor the worker ever mentioned, on the volume every job
    shares. It now resolves both names htslib resolves, and when it does have to
    build one it builds it in the run's output directory.

    So there is nothing here for cleanup to find, which is the point: this test
    pins the *absence*, and the removal of an index that does exist under a name
    the client never sent is pinned by the test below.

    Args:
        vntyper_task: The task fixture above.
        tmp_path: Scratch directory standing in for the job tree.
    """
    invoke, _ = vntyper_task
    alignment, index = _job_input(tmp_path, "sample.bam", "sample.bai")

    _run(invoke, alignment, tmp_path, index)

    assert not Path(f"{alignment}.bai").exists(), "the pipeline built an index beside the patient's alignment"
    assert not (tmp_path / "output" / "job-1" / "output_input.bam.bai").exists(), (
        "the uploaded sample.bai was ignored and a redundant index was built anyway"
    )
    assert not (tmp_path / "input" / "job-1").exists()


def test_an_index_under_a_name_the_client_never_sent_is_still_removed(vntyper_task, tmp_path: Path) -> None:
    """Cleanup's reach over derived index names has to keep working.

    With no index uploaded, the worker builds `sample.bam.bai` itself before
    starting the pipeline -- a name the client never sent. Cleanup used to remove
    only the alignment and the exact path the client sent, then find the
    directory non-empty and deliberately leave it, so that file stayed on the
    shared volume for good. It is the derived-name removal (`derived_index_paths`)
    that closes it, and this is that removal observed end to end.

    Args:
        vntyper_task: The task fixture above.
        tmp_path: Scratch directory standing in for the job tree.
    """
    invoke, commands = vntyper_task
    alignment, _ = _job_input(tmp_path, "sample.bam", None)

    _run(invoke, alignment, tmp_path, None)

    assert ["samtools", "index", str(alignment)] in commands, "the worker never built the index this test is about"
    assert not Path(f"{alignment}.bai").exists(), "the index the worker built was left on the shared volume"
    assert not (tmp_path / "input" / "job-1").exists()


def test_cleanup_leaves_a_file_that_is_not_a_derived_index_name(
    vntyper_task, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Widening the removal set must not turn cleanup into a directory wipe.

    Only the alignment, the index the submission carried, and the index names the
    job's own tooling can deterministically produce for that alignment are
    removed. Anything else in the directory is still reported and left where it
    is -- the existing guarantee, which a glob over the job directory would have
    quietly replaced.

    Args:
        vntyper_task: The task fixture above.
        tmp_path: Scratch directory standing in for the job tree.
        caplog: Captures the task's log record.
    """
    invoke, _ = vntyper_task
    alignment, index = _job_input(tmp_path, "sample.bam", "sample.bai")
    bystander = alignment.parent / "unrelated.bai.txt"
    bystander.write_bytes(b"not this job's to delete")
    caplog.set_level(logging.WARNING, logger="app.tasks")

    _run(invoke, alignment, tmp_path, index)

    assert bystander.read_bytes() == b"not this job's to delete"
    assert f"Input directory {alignment.parent} still holds files and was left in place" in caplog.text


# ---------------------------------------------------------------------------
# Nothing in the bookkeeping around cleanup can stop it happening
# ---------------------------------------------------------------------------


def _redis_is_down(*_args, **_kwargs):
    """Stand in for a Redis that became unreachable as the task exits.

    Raises:
        ConnectionError: Always.
    """
    msg = "Error 111 connecting to redis:6379. Connection refused."
    raise ConnectionError(msg)


def test_a_redis_failure_on_the_way_out_still_clears_the_input_directory(
    vntyper_task, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The queue-list write and the file removal share a `finally`; one must not veto the other.

    `lrem` removing the task from `vntyper_job_queue` is display bookkeeping. The
    removals after it are the only thing that ever takes patient-derived data off
    the shared volume. Ordering them the other way round, unguarded, means a Redis
    outage at exactly the wrong moment keeps an alignment on disk indefinitely.

    Args:
        vntyper_task: The task fixture above.
        fake_redis: The in-process Redis the task's three clients share.
        monkeypatch: Standard pytest fixture; restores the method at teardown.
        tmp_path: Scratch directory standing in for the job tree.
    """
    invoke, _ = vntyper_task
    alignment, index = _job_input(tmp_path, "sample.bam", "sample.bam.bai")
    monkeypatch.setattr(fake_redis, "lrem", _redis_is_down)

    _run(invoke, alignment, tmp_path, index)

    assert not alignment.exists(), "the alignment survived a Redis failure in the cleanup block"
    assert not (tmp_path / "input" / "job-1").exists()


def test_a_redis_failure_on_the_way_out_does_not_replace_the_pipeline_failure(
    vntyper_task, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A job that failed must still report why, not what went wrong while tidying up.

    An exception raised inside `finally` supersedes the one propagating through
    it, so an unguarded bookkeeping call turns "the pipeline exited non-zero" into
    "the connection to Redis was refused" in the recorded task result.

    Args:
        vntyper_task: The task fixture above.
        fake_redis: The in-process Redis the task's three clients share.
        monkeypatch: Standard pytest fixture; restores every patch at teardown.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    invoke, _ = vntyper_task
    alignment, index = _job_input(tmp_path, "sample.bam", "sample.bam.bai")

    def _pipeline_exits_nonzero(command, *_args, **_kwargs):
        """Fail the `vntyper pipeline` call the way a bad exit code does.

        Args:
            command: The argument vector the task asked for.
            *_args: Ignored.
            **_kwargs: Ignored.

        Raises:
            subprocess.CalledProcessError: Always.
        """
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(tasks.subprocess, "run", _pipeline_exits_nonzero)
    monkeypatch.setattr(fake_redis, "lrem", _redis_is_down)

    with pytest.raises(subprocess.CalledProcessError):
        _run(invoke, alignment, tmp_path, index)

    assert not alignment.exists()
    assert not (tmp_path / "input" / "job-1").exists()


def test_a_retention_failure_on_the_way_out_still_clears_the_input_directory(
    vntyper_task, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Extending a cohort's TTL is the other unguarded call ahead of the removals.

    It only runs for cohort jobs, so a cohort submission is the one that loses its
    cleanup to it.

    Args:
        vntyper_task: The task fixture above.
        monkeypatch: Standard pytest fixture; restores the patch at teardown.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    invoke, _ = vntyper_task
    alignment, index = _job_input(tmp_path, "sample.bam", "sample.bam.bai")
    monkeypatch.setattr(tasks, "extend_cohort_retention", _redis_is_down)

    _run(invoke, alignment, tmp_path, index, cohort_key="cohort:family-a")

    assert not alignment.exists(), "the alignment survived a retention failure in the cleanup block"
    assert not (tmp_path / "input" / "job-1").exists()


def test_the_queue_entry_is_still_removed_when_nothing_goes_wrong(vntyper_task, fake_redis, tmp_path: Path) -> None:
    """Guarding the bookkeeping must not amount to skipping it.

    Without this, every test above would also pass against a worker that had
    stopped writing to `vntyper_job_queue` at all.

    Args:
        vntyper_task: The task fixture above.
        fake_redis: The in-process Redis the task's three clients share.
        tmp_path: Scratch directory standing in for the job tree.
    """
    invoke, _ = vntyper_task
    alignment, index = _job_input(tmp_path, "sample.bam", "sample.bam.bai")
    fake_redis.rpush("vntyper_job_queue", "task-1")

    _run(invoke, alignment, tmp_path, index)

    assert fake_redis.lrange("vntyper_job_queue", 0, -1) == []
