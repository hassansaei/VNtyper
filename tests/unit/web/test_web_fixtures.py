"""Smoke tests proving the docker/app web-service fixtures actually work.

These are not feature tests for docker/app itself - G1 through G5 own the
security-relevant endpoint behaviour. This file only proves the `web_app` and
`client` fixtures in `conftest.py` wire up a real, working FastAPI app: Redis
calls land in `fake_redis`, Celery tasks never run, rate limiting never blocks
a request, and no test can slip past the socket guard unnoticed.
"""

import socket

import pytest
from fastapi import FastAPI

pytestmark = pytest.mark.unit


def test_web_app_is_a_real_fastapi_app(web_app) -> None:
    """`web_app` exposes the actual FastAPI instance from docker/app/main.py.

    Args:
        web_app: Patched `app.main` module from conftest.
    """
    assert isinstance(web_app.app, FastAPI)
    assert web_app.app.title == "VNtyper Online API"


def test_health_endpoint_round_trips_through_the_asgi_app(client) -> None:
    """The health endpoint responds with real content, not just <500.

    Args:
        client: TestClient fixture from conftest.
    """
    response = client.get("/health/")

    assert response.status_code == 200
    assert response.json() == {"status": "ok"}


def test_rate_limiting_is_neutralised_for_both_limiter_tiers(client) -> None:
    """Requests reach real endpoint logic instead of raising for a missing limiter.

    `/version/` sits behind `simple_rate_limiter`, `/download/{job_id}/` sits
    behind `high_rate_limiter`. Without the dependency_overrides in `web_app`,
    both raise "You must call FastAPILimiter.init" before reaching the
    handler; here they must reach it and return their real status codes.

    Args:
        client: TestClient fixture from conftest.
    """
    version_response = client.get("/version/")
    assert version_response.status_code == 200
    assert "api_version" in version_response.json()

    download_response = client.get("/download/does-not-exist/")
    assert download_response.status_code == 404
    assert download_response.json() == {"detail": "File not found"}


def test_cohort_create_and_status_round_trip_through_fake_redis(client, fake_redis) -> None:
    """A cohort created via the API is readable back out of `fake_redis`.

    Exercises the full write/read path (`create_cohort` -> `redis_cohort_client.hset`
    -> `get_cohort_status` -> `redis_cohort_client.hgetall`) to prove the
    `fake_redis` fixture is the same store both the app and the test see, not
    a separate no-op stand-in.

    The passphrase is not incidental: a cohort cannot be created without one,
    and cannot be read back without it either. Authorization itself is covered
    in `test_cohort_authz.py`; here the credential is only what makes the
    round trip reachable.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The fakeredis instance backing the app's Redis clients.
    """
    passphrase = "g0-smoke-passphrase"
    create_response = client.post("/create-cohort/", data={"alias": "g0-smoke-cohort", "passphrase": passphrase})
    assert create_response.status_code == 200
    cohort_id = create_response.json()["cohort_id"]

    assert fake_redis.hget(f"cohort:{cohort_id}", "alias") == "g0-smoke-cohort"

    status_response = client.get("/cohort-status/", params={"cohort_id": cohort_id, "passphrase": passphrase})
    assert status_response.status_code == 200
    body = status_response.json()
    assert body["cohort_id"] == cohort_id
    assert body["alias"] == "g0-smoke-cohort"
    assert body["jobs"] == []


def test_fake_redis_is_a_working_in_memory_client(fake_redis) -> None:
    """`fake_redis` behaves like a real Redis client without touching a socket.

    Args:
        fake_redis: The fakeredis instance under test.
    """
    assert fake_redis.ping() is True

    fake_redis.set("g0-smoke-key", "g0-smoke-value")
    assert fake_redis.get("g0-smoke-key") == "g0-smoke-value"


def test_run_vntyper_job_is_mocked_not_the_real_celery_task(web_app) -> None:
    """The Celery task object is a MagicMock, so no subprocess or broker is touched.

    Args:
        web_app: Patched `app.main` module from conftest.
    """
    result = web_app.run_vntyper_job.delay(bam_path="/tmp/fake.bam")

    assert result.id == "task-1"
    web_app.run_vntyper_job.delay.assert_called_once()


def test_socket_guard_blocks_real_network_connections() -> None:
    """The module-scoped `_no_real_sockets` guard rejects a real connect() call.

    This proves the guard from conftest.py is actually installed and active
    for this test module, rather than a comment nobody enforces.
    """
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        with pytest.raises(RuntimeError, match="real socket connection"):
            sock.connect(("example.invalid", 80))
    finally:
        sock.close()
