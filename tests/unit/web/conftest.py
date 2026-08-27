"""Fixtures for the FastAPI + Celery web service (docker/app) tests.

docker/app is not an installed package, so it is put on sys.path here rather
than in pyproject. Everything below is in-process: TestClient uses an ASGI
transport (no sockets) and fakeredis is pure Python, so these tests satisfy the
unit-tier purity rule and belong in the `unit` marker, not a new one.
"""

import socket
import sys
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

pytestmark = pytest.mark.unit

_DOCKER_DIR = Path(__file__).resolve().parents[3] / "docker"
if str(_DOCKER_DIR) not in sys.path:
    sys.path.insert(0, str(_DOCKER_DIR))


@pytest.fixture(autouse=True, scope="module")
def _no_real_sockets() -> Iterator[None]:
    """Block real socket connections for every test under tests/unit/web.

    docker/app talks to Redis at request time, and `validate_email()` in
    `main.run_vntyper` performs a live DNS lookup unless deliverability
    checking is disabled. Neither should ever fire in the unit tier: Redis
    access must go through the `fake_redis`/`web_app` fixtures, and any test
    that exercises the email path must stub `validate_email` itself. If
    either slips through, fail loudly instead of silently hitting a real
    service.

    Yields:
        None: This fixture only wraps test execution to install and remove
            the guard.
    """
    mp = pytest.MonkeyPatch()

    def _blocked_connect(self, address, *_args, **_kwargs):
        raise RuntimeError(
            f"tests/unit/web attempted a real socket connection to {address!r}. "
            "Use the fake_redis/web_app/client fixtures instead of live Redis or DNS."
        )

    mp.setattr(socket.socket, "connect", _blocked_connect)
    yield
    mp.undo()


@pytest.fixture(autouse=True)
def _ordinary_worker_storage_roots(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Point worker storage contracts at each test's isolated directories."""
    from app import tasks

    monkeypatch.setattr(tasks.settings, "DEFAULT_HANDOFF_SPOOL_DIR", str(tmp_path / "input"))
    monkeypatch.setattr(tasks.settings, "DEFAULT_INPUT_DIR", str(tmp_path / "input"))
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(tmp_path / "output"))


@pytest.fixture
def fake_redis():
    """Provide an in-process Redis stand-in.

    Returns:
        fakeredis.FakeRedis: A pure-Python client with the same interface as
            `redis.Redis(decode_responses=True)`, backed by an in-memory store
            instead of a socket.
    """
    import fakeredis

    return fakeredis.FakeRedis(decode_responses=True)


@pytest.fixture
def web_app(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fake_redis) -> Iterator[object]:
    """Provide the FastAPI app with Redis, Celery and rate limiting neutralised.

    Patches the real module-level names in `docker/app/main.py` (verified by
    reading the module, not guessed): the three `redis.Redis` clients, the two
    job-submission directories, and the two imported Celery task objects. Rate
    limiting is neutralised via `app.dependency_overrides` rather than by
    running the `startup` event, because `startup_event()` calls
    `FastAPILimiter.init()` against a real `aioredis` connection.

    Args:
        monkeypatch: Standard pytest fixture; also auto-restores every patch
            below at teardown.
        tmp_path: Per-test scratch directory used for the job input/output
            trees.
        fake_redis: In-process Redis stand-in from the `fake_redis` fixture.

    Yields:
        object: The imported `app.main` module, with `.app` holding the
            FastAPI instance.
    """
    from app import main as app_main

    for attr in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(app_main, attr, fake_redis)

    (tmp_path / "input").mkdir(exist_ok=True)
    (tmp_path / "handoff").mkdir(exist_ok=True)
    (tmp_path / "output").mkdir(exist_ok=True)
    monkeypatch.setattr(app_main, "DEFAULT_INPUT_DIR", str(tmp_path / "input"))
    monkeypatch.setattr(app_main, "DEFAULT_HANDOFF_SPOOL_DIR", str(tmp_path / "handoff"), raising=False)
    monkeypatch.setattr(app_main, "DEFAULT_OUTPUT_DIR", str(tmp_path / "output"))

    for attr in ("run_vntyper_job", "run_cohort_analysis_job"):
        task = MagicMock()
        task.delay.return_value = SimpleNamespace(id="task-1")
        task.apply_async.return_value = SimpleNamespace(id="task-1")
        monkeypatch.setattr(app_main, attr, task)

    async def _noop_rate_limiter() -> None:
        """Stand in for `RateLimiter.__call__` without touching Redis."""
        return None

    app_main.app.dependency_overrides[app_main.simple_rate_limiter] = _noop_rate_limiter
    app_main.app.dependency_overrides[app_main.high_rate_limiter] = _noop_rate_limiter

    yield app_main

    app_main.app.dependency_overrides.clear()


@pytest.fixture
def client(web_app):
    """Provide a TestClient that does NOT run the startup event.

    Entering the context manager (`with TestClient(...) as client`) would fire
    the `startup` event, which calls `FastAPILimiter.init()` against a real
    `aioredis` connection. Constructing without `with` skips lifespan
    entirely; rate limiting is neutralised separately in `web_app` via
    `dependency_overrides`.

    Args:
        web_app: The patched `app.main` module from the `web_app` fixture.

    Returns:
        starlette.testclient.TestClient: An in-process ASGI client for
            `web_app.app`.
    """
    from fastapi.testclient import TestClient

    return TestClient(web_app.app)
