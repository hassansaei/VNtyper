"""FastAPI startup and shutdown own the rate-limiter Redis connection."""

import asyncio
import subprocess
from types import SimpleNamespace

import pytest
from app.lifecycle import AsyncRedisClient, close_rate_limiter_client, initialize_rate_limiter, probe_tool_version

pytestmark = pytest.mark.unit


class _RedisClient:
    def __init__(self, close_failure: BaseException | None = None) -> None:
        self.close_count = 0
        self.close_failure = close_failure

    async def aclose(self) -> None:
        self.close_count += 1
        if self.close_failure is not None:
            raise self.close_failure


class _Limiter:
    def __init__(self, failure: BaseException | None = None) -> None:
        self.client: AsyncRedisClient | None = None
        self.failure = failure

    async def init(self, client: AsyncRedisClient) -> None:
        self.client = client
        if self.failure is not None:
            raise self.failure


def test_rate_limiter_initialization_returns_the_owned_client() -> None:
    client = _RedisClient()
    factory_calls = []

    def redis_factory(url: str, **options: object) -> _RedisClient:
        factory_calls.append((url, options))
        return client

    limiter = _Limiter()

    result = asyncio.run(initialize_rate_limiter("redis://rate-limit/2", redis_factory=redis_factory, limiter=limiter))

    assert result is client
    assert limiter.client is client
    assert factory_calls == [
        ("redis://rate-limit/2", {"encoding": "utf8", "decode_responses": True}),
    ]
    assert client.close_count == 0


def test_failed_rate_limiter_initialization_closes_the_client_before_propagating() -> None:
    client = _RedisClient()
    limiter = _Limiter(RuntimeError("script load failed"))

    with pytest.raises(RuntimeError, match="script load failed"):
        asyncio.run(
            initialize_rate_limiter(
                "redis://rate-limit/2",
                redis_factory=lambda *_args, **_kwargs: client,
                limiter=limiter,
            )
        )

    assert client.close_count == 1


def test_initialization_failure_remains_primary_when_client_close_also_fails(caplog: pytest.LogCaptureFixture) -> None:
    client = _RedisClient(RuntimeError("close failed"))
    limiter = _Limiter(RuntimeError("script load failed"))

    with pytest.raises(RuntimeError, match="script load failed"):
        asyncio.run(
            initialize_rate_limiter(
                "redis://rate-limit/2",
                redis_factory=lambda *_args, **_kwargs: client,
                limiter=limiter,
            )
        )

    assert client.close_count == 1
    assert "Failed to close rate-limiting Redis after startup failure: close failed" in caplog.text


def test_rate_limiter_initialization_cancellation_closes_the_owned_client() -> None:
    client = _RedisClient()
    limiter = _Limiter(asyncio.CancelledError())

    with pytest.raises(asyncio.CancelledError):
        asyncio.run(
            initialize_rate_limiter(
                "redis://rate-limit/2",
                redis_factory=lambda *_args, **_kwargs: client,
                limiter=limiter,
            )
        )

    assert client.close_count == 1


def test_rate_limiter_shutdown_uses_aclose_exactly_once() -> None:
    client = _RedisClient()

    asyncio.run(close_rate_limiter_client(client))

    assert client.close_count == 1


def test_rate_limiter_shutdown_propagates_a_close_failure() -> None:
    client = _RedisClient(RuntimeError("close failed"))

    with pytest.raises(RuntimeError, match="close failed"):
        asyncio.run(close_rate_limiter_client(client))

    assert client.close_count == 1


@pytest.mark.parametrize(
    ("effect", "expected"),
    [
        ("VNtyper 2.0\n", "VNtyper 2.0"),
        (subprocess.CalledProcessError(1, ["vntyper", "-v"], output="bad version\n"), "error retrieving tool version"),
        (FileNotFoundError(), "VNtyper tool not installed"),
        (subprocess.TimeoutExpired(["vntyper", "-v"], 5), "timeout retrieving tool version"),
    ],
)
def test_tool_version_probe_preserves_every_public_outcome(effect: object, expected: str) -> None:
    def runner(*_args: object, **_kwargs: object) -> str:
        if isinstance(effect, BaseException):
            raise effect
        return str(effect)

    assert probe_tool_version(runner) == expected


def test_application_lifespan_closes_the_client_after_an_exception(web_app, monkeypatch: pytest.MonkeyPatch) -> None:
    client = _RedisClient()

    async def startup() -> _RedisClient:
        return client

    monkeypatch.setattr(web_app, "startup_event", startup)

    async def exercise() -> None:
        with pytest.raises(RuntimeError, match="request loop failed"):
            async with web_app.lifespan(SimpleNamespace()):
                raise RuntimeError("request loop failed")

    asyncio.run(exercise())

    assert client.close_count == 1


def test_startup_returns_the_client_and_caches_the_tool_version(web_app, monkeypatch: pytest.MonkeyPatch) -> None:
    client = _RedisClient()

    async def initialize(_url: str) -> _RedisClient:
        return client

    monkeypatch.setattr(web_app, "require_redis_password", lambda: "secret")
    monkeypatch.setattr(web_app, "initialize_rate_limiter", initialize)
    monkeypatch.setattr(web_app, "probe_tool_version", lambda: "vntyper 2.0.26")

    result = asyncio.run(web_app.startup_event())

    assert result is client
    assert web_app.TOOL_VERSION == "vntyper 2.0.26"
    assert client.close_count == 0


def test_startup_closes_the_client_if_the_version_probe_fails(web_app, monkeypatch: pytest.MonkeyPatch) -> None:
    client = _RedisClient()

    async def initialize(_url: str) -> _RedisClient:
        return client

    def fail_probe() -> str:
        raise RuntimeError("unexpected version failure")

    monkeypatch.setattr(web_app, "require_redis_password", lambda: "secret")
    monkeypatch.setattr(web_app, "initialize_rate_limiter", initialize)
    monkeypatch.setattr(web_app, "probe_tool_version", fail_probe)

    with pytest.raises(RuntimeError, match="unexpected version failure"):
        asyncio.run(web_app.startup_event())

    assert client.close_count == 1


def test_startup_cancellation_during_version_probe_closes_the_client(web_app, monkeypatch: pytest.MonkeyPatch) -> None:
    client = _RedisClient()

    async def initialize(_url: str) -> _RedisClient:
        return client

    def cancel_probe() -> str:
        raise asyncio.CancelledError

    monkeypatch.setattr(web_app, "require_redis_password", lambda: "secret")
    monkeypatch.setattr(web_app, "initialize_rate_limiter", initialize)
    monkeypatch.setattr(web_app, "probe_tool_version", cancel_probe)

    with pytest.raises(asyncio.CancelledError):
        asyncio.run(web_app.startup_event())

    assert client.close_count == 1


def test_version_probe_failure_remains_primary_when_startup_cleanup_fails(
    web_app, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    client = _RedisClient(RuntimeError("close failed"))

    async def initialize(_url: str) -> _RedisClient:
        return client

    def fail_probe() -> str:
        raise RuntimeError("unexpected version failure")

    monkeypatch.setattr(web_app, "require_redis_password", lambda: "secret")
    monkeypatch.setattr(web_app, "initialize_rate_limiter", initialize)
    monkeypatch.setattr(web_app, "probe_tool_version", fail_probe)

    with pytest.raises(RuntimeError, match="unexpected version failure"):
        asyncio.run(web_app.startup_event())

    assert client.close_count == 1
    assert "Failed to close rate-limiting Redis after version-probe failure: close failed" in caplog.text


def test_lifespan_body_failure_remains_primary_when_shutdown_cleanup_fails(
    web_app, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    client = _RedisClient(RuntimeError("close failed"))

    async def startup() -> _RedisClient:
        return client

    monkeypatch.setattr(web_app, "startup_event", startup)

    async def exercise() -> None:
        with pytest.raises(RuntimeError, match="request loop failed"):
            async with web_app.lifespan(SimpleNamespace()):
                raise RuntimeError("request loop failed")

    asyncio.run(exercise())

    assert client.close_count == 1
    assert "Failed to close rate-limiting Redis while preserving lifespan failure: close failed" in caplog.text


def test_application_registers_lifespan_without_deprecated_event_handlers(
    web_app, monkeypatch: pytest.MonkeyPatch
) -> None:
    client = _RedisClient()

    async def startup() -> _RedisClient:
        return client

    monkeypatch.setattr(web_app, "startup_event", startup)

    async def exercise_registered_lifespan() -> None:
        async with web_app.app.router.lifespan_context(web_app.app):
            assert client.close_count == 0

    asyncio.run(exercise_registered_lifespan())

    assert web_app.app.router.on_startup == []
    assert web_app.app.router.on_shutdown == []
    assert client.close_count == 1
