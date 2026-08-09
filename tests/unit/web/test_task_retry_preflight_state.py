"""Focused worker-retry tests for stale curated preflight state."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from tests.unit.web.test_tasks import (
    _invoke_vntyper_job,
    _job_kwargs,
    _make_job_input,
    _subprocess_stub,
)

pytestmark = pytest.mark.unit


@pytest.fixture
def redis_mocks(monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Replace queue, cohort and usage Redis clients with independent mocks."""
    from app import tasks

    mocks = SimpleNamespace(
        queue=MagicMock(name="redis_client"),
        cohort=MagicMock(name="redis_cohort_client"),
        usage=MagicMock(name="redis_usage_client"),
    )
    monkeypatch.setattr(tasks, "redis_client", mocks.queue)
    monkeypatch.setattr(tasks, "redis_cohort_client", mocks.cohort)
    monkeypatch.setattr(tasks, "redis_usage_client", mocks.usage)
    return mocks


@pytest.fixture
def no_email_task(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Prevent worker failure handling from contacting the Celery broker."""
    from app import tasks

    mocked = MagicMock(name="send_email_task")
    monkeypatch.setattr(tasks, "send_email_task", mocked)
    return mocked


@pytest.mark.parametrize("artifact_kind", ["regular", "symlink"])
def test_reused_job_clears_stale_preflight_state_before_a_new_attempt(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    fake_redis,
    tmp_path: Path,
    artifact_kind: str,
) -> None:
    """A later generic failure cannot surface a prior attempt's diagnosis.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The mocked queue and cohort Redis clients.
        no_email_task: The mocked email task.
        fake_redis: Stateful in-memory usage Redis client.
        tmp_path: Scratch directory standing in for the reused job tree.
        artifact_kind: Whether stale state is a regular file or an untrusted symlink.
    """
    from app import tasks

    del redis_mocks, no_email_task
    bam_path, _ = _make_job_input(tmp_path)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    stale_payload = {
        "code": "reference_unresolved",
        "message": "Unable to resolve CRAM reference: contig=chr1, M5=old-digest.",
        "candidates": [["cli", "full.fa", "probe exited non-zero"]],
    }
    stale_artifact = output_dir / "preflight_error.json"
    protected_artifact = tmp_path / "protected-preflight-error.json"
    if artifact_kind == "regular":
        stale_artifact.write_text(json.dumps(stale_payload), encoding="utf-8")
    else:
        protected_artifact.write_text(json.dumps(stale_payload), encoding="utf-8")
        stale_artifact.symlink_to(protected_artifact)

    fake_redis.hset(
        "usage:job-1",
        mapping={"status": "failed", "code": stale_payload["code"], "message": stale_payload["message"]},
    )
    monkeypatch.setattr(tasks, "redis_usage_client", fake_redis)
    pipeline_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])
    _subprocess_stub(monkeypatch, tasks, pipeline_error=pipeline_error)

    with pytest.raises(subprocess.CalledProcessError):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    assert not stale_artifact.exists()
    assert fake_redis.hget("usage:job-1", "status") == "failed"
    assert fake_redis.hget("usage:job-1", "code") is None
    assert fake_redis.hget("usage:job-1", "message") is None
    if artifact_kind == "symlink":
        assert protected_artifact.read_text(encoding="utf-8") == json.dumps(stale_payload)
