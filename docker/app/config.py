import os
import logging


class Settings:
    PROJECT_NAME: str = "VNtyper API"
    DEBUG: bool = True

    # Logging configuration
    LOG_LEVEL: str = os.getenv("LOG_LEVEL", "INFO")
    LOG_FORMAT: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Directory configurations
    DEFAULT_INPUT_DIR: str = os.getenv("DEFAULT_INPUT_DIR", "/opt/vntyper/input")
    DEFAULT_OUTPUT_DIR: str = os.getenv("DEFAULT_OUTPUT_DIR", "/opt/vntyper/output")

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

    # Cohort configurations
    COHORT_RETENTION_DAYS: int = int(os.getenv("COHORT_RETENTION_DAYS", 14))  # Default to 14 days
    PASSWORD_HASH_SCHEME: str = "bcrypt"

    # Max result age for cleanup
    MAX_RESULT_AGE_DAYS: int = int(os.getenv("MAX_RESULT_AGE_DAYS", 7))

    # Usage statistics configurations
    USAGE_REDIS_DB: int = int(os.getenv("USAGE_REDIS_DB", 4))
    USAGE_DATA_RETENTION_DAYS: int = int(os.getenv("USAGE_DATA_RETENTION_DAYS", 30))
    USAGE_DATA_RETENTION_SECONDS: int = USAGE_DATA_RETENTION_DAYS * 86400


settings = Settings()

# Configure logging
logging.basicConfig(level=settings.LOG_LEVEL, format=settings.LOG_FORMAT)
logger = logging.getLogger(__name__)
