"""Service-level hardening pins for docker/app.

The service must ship with debug behavior disabled. Nothing reads this setting
today, which is exactly why a ``True`` default would be dangerous: a future
reader would inherit debug behavior nobody chose.
"""

import inspect

import pytest

pytestmark = pytest.mark.unit

from app.config import Settings  # noqa: E402


def test_debug_is_not_a_shipped_default() -> None:
    """The service must never default to debug behavior."""
    assert Settings.DEBUG is False
    assert Settings().DEBUG is False


def test_run_job_handler_is_not_a_coroutine_function(web_app) -> None:
    """A blocking handler must be a ``def``, so FastAPI runs it off the loop.

    Args:
        web_app: Patched ``app.main`` module from conftest.
    """
    assert not inspect.iscoroutinefunction(web_app.run_vntyper)


def test_run_job_still_submits_after_the_dispatch_change(client, web_app, tmp_path) -> None:
    """The endpoint's observable behavior is unchanged by the dispatch fix.

    Args:
        client: TestClient fixture from conftest.
        web_app: Patched ``app.main`` module, sharing this test's ``tmp_path``.
        tmp_path: The same directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"bamdata", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 200
    job_id = response.json()["job_id"]
    assert (tmp_path / "handoff" / job_id / "sample.bam").read_bytes() == b"bamdata"
    web_app.run_vntyper_job.delay.assert_called_once()
