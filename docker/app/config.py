import logging
import os
from urllib.parse import quote


class Settings:
    PROJECT_NAME: str = "VNtyper API"
    DEBUG: bool = True

    # Logging configuration
    LOG_LEVEL: str = os.getenv("LOG_LEVEL", "INFO")
    LOG_FORMAT: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Directory configurations
    DEFAULT_INPUT_DIR: str = os.getenv("DEFAULT_INPUT_DIR", "/opt/vntyper/input")
    DEFAULT_OUTPUT_DIR: str = os.getenv("DEFAULT_OUTPUT_DIR", "/opt/vntyper/output")

    # Upload size ceiling, in bytes. Caps the total an individual job submission
    # may write into the input directory above; the rate limits below bound how
    # many requests arrive, not how large each one is. The documented workflow
    # uploads a MUC1-region subset measured in megabytes, so 1 GiB leaves ample
    # headroom for an unsubsetted regional or exome-scale BAM while still
    # bounding what a single request can consume. Deployments with a different
    # volume budget override it with MAX_UPLOAD_BYTES.
    MAX_UPLOAD_BYTES: int = int(os.getenv("MAX_UPLOAD_BYTES", 1024 * 1024 * 1024))

    # Rate limiting configurations
    RATE_LIMIT_SIMPLE_TIMES: int = int(os.getenv("RATE_LIMIT_SIMPLE_TIMES", 100))
    RATE_LIMIT_SIMPLE_SECONDS: int = int(os.getenv("RATE_LIMIT_SIMPLE_SECONDS", 60))
    RATE_LIMIT_HIGH_TIMES: int = int(os.getenv("RATE_LIMIT_HIGH_TIMES", 10))
    RATE_LIMIT_HIGH_SECONDS: int = int(os.getenv("RATE_LIMIT_HIGH_SECONDS", 60))
    RATE_LIMITING_REDIS_DB: int = int(os.getenv("RATE_LIMITING_REDIS_DB", 2))  # Use DB 2 for rate limiting

    # SMTP Configuration
    SMTP_HOST: str = os.getenv("SMTP_HOST", "smtp.hoster.com")
    SMTP_PORT: int = int(os.getenv("SMTP_PORT", 587))
    SMTP_USERNAME: str = os.getenv("SMTP_USERNAME", "your_smtp_username")
    SMTP_PASSWORD: str = os.getenv("SMTP_PASSWORD", "your_smtp_password")
    EMAIL_FROM: str = os.getenv("EMAIL_FROM", "noreply@hoster.com")

    # API Base URL
    API_BASE_URL: str = os.getenv("API_BASE_URL", "http://localhost:8000")

    # How long a completed job's result archive stays retrievable, in days.
    #
    # `/download/{job_id}/` requires no credential (#189), so this value *is* the exposure
    # window for a completed job: the id is the capability. Shortening it does not fix that
    # design, and is not claimed to -- it reduces exposure duration, not exposure.
    #
    # 3 days rather than 7, chosen from what the code assumes about collection rather than by
    # feel: `vntyper online` polls for at most 4 hours
    # (`online_mode.DEFAULT_POLL_TIMEOUT_SECONDS`) and then gives up, and `config.json` has no
    # `api_settings` block overriding it, so the CLI path never needs days. The binding case is
    # a person opening the completion email later, and 3 days covers submit-Friday,
    # download-Monday.
    MAX_RESULT_AGE_DAYS: int = int(os.getenv("MAX_RESULT_AGE_DAYS", 3))

    # Cohort configurations
    #
    # Read `COHORT_RETENTION_DAYS` through `cohort_retention_days()` rather than directly: a
    # cohort must not outlive the archives it lists. See that method for why.
    COHORT_RETENTION_DAYS: int = int(os.getenv("COHORT_RETENTION_DAYS", 14))  # Default to 14 days
    PASSWORD_HASH_SCHEME: str = "bcrypt"

    # Usage statistics configurations
    USAGE_REDIS_DB: int = int(os.getenv("USAGE_REDIS_DB", 4))
    USAGE_DATA_RETENTION_DAYS: int = int(os.getenv("USAGE_DATA_RETENTION_DAYS", 30))
    USAGE_DATA_RETENTION_SECONDS: int = USAGE_DATA_RETENTION_DAYS * 86400

    def cohort_retention_days(self) -> int:
        """
        Effective cohort lifetime, clamped so a cohort never outlives its members' archives.

        A cohort is openable for ``COHORT_RETENTION_DAYS`` while each member's ``.zip`` is
        reclaimed at ``MAX_RESULT_AGE_DAYS`` (``tasks.delete_old_results``). With the shipped
        14 and 3 that leaves an 11-day window in which ``/cohort-status/`` lists member job
        ids whose ``/download/{job_id}/`` returns 404 -- a cohort page advertising downloads
        that do not work.

        It is worse than a static comparison suggests: ``extend_cohort_retention`` fires on
        every job completion and every cohort analysis, so an actively-used cohort rolls its
        TTL forward from its *most recent* job while each archive ages out on a fixed clock
        from *its own* creation. An active cohort can stay open with almost none of its
        members intact.

        Clamping is the honest resolution -- the alternative, documenting the gap, leaves the
        page listing broken downloads. A WARNING names both configured values when the clamp
        bites, so an operator who set 14 learns that 14 is not what they got.

        Returns:
            int: ``min(COHORT_RETENTION_DAYS, MAX_RESULT_AGE_DAYS)``.
        """
        if self.COHORT_RETENTION_DAYS <= self.MAX_RESULT_AGE_DAYS:
            return self.COHORT_RETENTION_DAYS

        # `logging.getLogger` rather than the module-level `logger`, which is defined below
        # this class: that name resolves at call time and so would work, but depending on
        # definition order from inside a class body is a trap for whoever moves either one.
        logging.getLogger(__name__).warning(
            "COHORT_RETENTION_DAYS=%s outlives MAX_RESULT_AGE_DAYS=%s, so a cohort would stay "
            "openable after its members' result archives were deleted. Clamping cohort retention "
            "to %s days.",
            self.COHORT_RETENTION_DAYS,
            self.MAX_RESULT_AGE_DAYS,
            self.MAX_RESULT_AGE_DAYS,
        )
        return self.MAX_RESULT_AGE_DAYS


settings = Settings()

# Configure logging
logging.basicConfig(level=settings.LOG_LEVEL, format=settings.LOG_FORMAT)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Redis credential
#
# The API, the Celery worker and the beat scheduler all authenticate against the
# same Redis instance, so they must resolve the same value from the same place.
# Read here, once, with no fallback: a default in the source tree is a shared
# secret, and a default that differs between modules lets the processes come up
# authenticating differently against one instance.
# ---------------------------------------------------------------------------

REDIS_PASSWORD_ENV_VAR = "REDIS_PASSWORD"


def get_redis_password() -> str | None:
    """Return the configured Redis password, or None when it is unset.

    Read from the environment on every call rather than captured at import time,
    so importing an application module never depends on the deployment
    environment. Callers that must not proceed without a credential should use
    `require_redis_password` instead.

    Returns:
        str | None: The configured password, or None if REDIS_PASSWORD is unset
            or empty.
    """
    return os.environ.get(REDIS_PASSWORD_ENV_VAR) or None


def require_redis_password() -> str:
    """Return the configured Redis password, failing fast when it is missing.

    Call this from process startup - the FastAPI startup event, the Celery
    worker and beat init signals - so a deployment that never set the variable
    stops there with an explicit message.

    Returns:
        str: The configured password.

    Raises:
        RuntimeError: If REDIS_PASSWORD is unset or empty.
    """
    password = get_redis_password()
    if not password:
        msg = (
            f"{REDIS_PASSWORD_ENV_VAR} is not set. Redis, the API and the Celery worker must all be "
            f"given the same password and there is no default. Set {REDIS_PASSWORD_ENV_VAR} to a "
            "freshly generated secret before starting the service."
        )
        logger.error(msg)
        raise RuntimeError(msg)
    return password


def build_redis_url(host: str, port: int, db: int, password: str | None = None) -> str:
    """Build a redis:// URL, percent-encoding the credential.

    A Redis password is an arbitrary string and routinely contains characters
    that are delimiters inside a URL authority (`@`, `:`, `/`, `#`).
    Interpolating one straight into an f-string yields a URL the client parses
    as a different host, port or path, so the encoding is done here rather than
    repeated at each call site.

    Args:
        host: Redis hostname.
        port: Redis port.
        db: Redis logical database number.
        password: Password to embed. When None or empty, the URL is built
            without a credential section.

    Returns:
        str: A URL safe to hand to redis-py, Celery or fastapi-limiter.
    """
    if not password:
        return f"redis://{host}:{port}/{db}"
    return f"redis://:{quote(password, safe='')}@{host}:{port}/{db}"
