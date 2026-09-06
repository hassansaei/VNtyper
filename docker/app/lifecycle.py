"""Resource ownership for the VNtyper API startup and shutdown lifecycle."""

import logging
import subprocess
from collections.abc import Callable
from typing import Protocol

import redis.asyncio as aioredis
from fastapi_limiter import FastAPILimiter

logger = logging.getLogger(__name__)


class AsyncRedisClient(Protocol):
    """Redis connection behavior owned by the application lifecycle."""

    async def aclose(self) -> None:
        """Close the client and its connection pool."""


class RedisFactory(Protocol):
    """Constructor shape used by ``redis.asyncio.from_url``."""

    def __call__(self, url: str, *, encoding: str, decode_responses: bool) -> AsyncRedisClient:
        """Create one asynchronous Redis client."""


class LimiterInitializer(Protocol):
    """Subset of ``FastAPILimiter`` needed during startup."""

    async def init(self, client: AsyncRedisClient) -> None:
        """Load the limiter script into the supplied Redis client."""


def patch_fastapi_limiter() -> None:
    """Patch fastapi-limiter for compatibility with FastAPI >= 0.137.

    FastAPI 0.137+ wraps included routers into `_IncludedRouter` objects which
    do not have a `.path` attribute, causing `RateLimiter.__call__` to crash with
    AttributeError when scanning `request.app.routes`.
    """
    from starlette.requests import Request
    from starlette.responses import Response
    from fastapi_limiter.depends import RateLimiter
    import redis.exceptions

    async def _patched_call(self: RateLimiter, request: Request, response: Response):
        if not FastAPILimiter.redis:
            raise Exception("You must call FastAPILimiter.init in startup event of fastapi!")

        route = request.scope.get("route")
        route_index = getattr(route, "path", request.scope.get("path", "0"))
        dep_index = 0
        if route and hasattr(route, "dependencies"):
            for j, dependency in enumerate(route.dependencies):
                if self is dependency.dependency:
                    dep_index = j
                    break

        identifier = self.identifier or FastAPILimiter.identifier
        callback = self.callback or FastAPILimiter.http_callback
        rate_key = await identifier(request)
        key = f"{FastAPILimiter.prefix}:{rate_key}:{route_index}:{dep_index}"
        try:
            pexpire = await self._check(key)
        except redis.exceptions.NoScriptError:
            FastAPILimiter.lua_sha = await FastAPILimiter.redis.script_load(
                FastAPILimiter.lua_script
            )
            pexpire = await self._check(key)
        if pexpire != 0:
            return await callback(request, response, pexpire)

    RateLimiter.__call__ = _patched_call


patch_fastapi_limiter()


async def initialize_rate_limiter(
    redis_url: str,
    *,
    redis_factory: RedisFactory = aioredis.from_url,
    limiter: LimiterInitializer = FastAPILimiter,
) -> AsyncRedisClient:
    """Initialize rate limiting and return the Redis client the app now owns.

    Args:
        redis_url: Authenticated Redis URL for the rate-limiter database.
        redis_factory: Constructor for the asynchronous Redis client.
        limiter: Limiter implementation whose script is initialized in Redis.

    Returns:
        The initialized Redis client, to be closed during application shutdown.

    Raises:
        Exception: Propagates Redis or limiter initialization failures after
            closing the client. A close failure is logged without replacing the
            initialization failure.
    """
    patch_fastapi_limiter()
    client = redis_factory(redis_url, encoding="utf8", decode_responses=True)
    try:
        await limiter.init(client)
    except BaseException:
        logger.error("Failed to initialize rate limiting.")
        try:
            await client.aclose()
        except BaseException as close_error:
            logger.error(f"Failed to close rate-limiting Redis after startup failure: {close_error}")
        raise

    logger.info("Rate limiting initialized successfully.")
    return client


async def close_rate_limiter_client(client: AsyncRedisClient) -> None:
    """Close the rate-limiter Redis client during application shutdown.

    Args:
        client: Client returned by :func:`initialize_rate_limiter`.

    Raises:
        Exception: Propagates a Redis close failure after logging it.
    """
    try:
        await client.aclose()
    except BaseException as error:
        logger.error(f"Failed to close rate-limiting Redis: {error}")
        raise


def probe_tool_version(command_runner: Callable[..., str] = subprocess.check_output) -> str:
    """Return the VNtyper CLI version or its established fallback text.

    Args:
        command_runner: Callable compatible with ``subprocess.check_output``.

    Returns:
        The stripped CLI version, or the existing public fallback string for a
        failed, missing, or timed-out executable.
    """
    try:
        output = command_runner(
            ["vntyper", "-v"],
            stderr=subprocess.STDOUT,
            text=True,
            timeout=5,
        )
        version = output.strip()
        logger.info(f"VNtyper tool version: {version}")
        return version
    except subprocess.CalledProcessError as error:
        logger.error(f"Error retrieving tool version: {error.output.strip()}")
        return "error retrieving tool version"
    except FileNotFoundError:
        logger.error("VNtyper tool not found.")
        return "VNtyper tool not installed"
    except subprocess.TimeoutExpired:
        logger.error("Timeout expired while retrieving tool version.")
        return "timeout retrieving tool version"
