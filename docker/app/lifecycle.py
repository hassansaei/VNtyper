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
