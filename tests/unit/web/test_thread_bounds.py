"""Thread-count form values are validated before a job is admitted."""

from pathlib import Path

import pytest

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("thread", ["0", "-1", "65", "999999"])
def test_out_of_range_thread_counts_are_refused_before_side_effects(
    client,
    web_app,
    tmp_path: Path,
    thread: str,
) -> None:
    """An invalid resource claim creates no input and enqueues no task.

    Args:
        client: In-process FastAPI test client.
        web_app: Patched web application module.
        tmp_path: Scratch directory containing the configured input root.
        thread: Out-of-range form token under test.
    """
    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"x", "application/octet-stream")},
        data={"thread": thread, "reference_assembly": "hg19"},
    )

    assert response.status_code == 422
    assert list((tmp_path / "input").iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()
    web_app.run_vntyper_job.apply_async.assert_not_called()


@pytest.mark.parametrize("thread", [1, 64])
def test_boundary_thread_counts_are_forwarded_unchanged(client, web_app, thread: int) -> None:
    """Both accepted endpoints of the interval reach the worker unchanged.

    Args:
        client: In-process FastAPI test client.
        web_app: Patched web application module.
        thread: Boundary value under test.
    """
    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"x", "application/octet-stream")},
        data={"thread": str(thread), "reference_assembly": "hg19"},
    )

    assert response.status_code == 200
    assert web_app.run_vntyper_job.delay.call_args.kwargs["thread"] == thread


def test_omitting_thread_preserves_the_default(client, web_app) -> None:
    """An omitted form field still forwards the shipped default of four.

    Args:
        client: In-process FastAPI test client.
        web_app: Patched web application module.
    """
    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"x", "application/octet-stream")},
        data={"reference_assembly": "hg19"},
    )

    assert response.status_code == 200
    assert web_app.run_vntyper_job.delay.call_args.kwargs["thread"] == 4
