# docker/app/config.py

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
    RATE_LIMIT_TIMES: int = int(os.getenv("RATE_LIMIT_TIMES", 10))  # Default to 10 requests
    RATE_LIMIT_SECONDS: int = int(os.getenv("RATE_LIMIT_SECONDS", 60))  # Per 60 seconds
    RATE_LIMITING_REDIS_DB: int = int(os.getenv("RATE_LIMITING_REDIS_DB", 2))  # Use DB 2 for rate limiting

settings = Settings()

# Configure logging
logging.basicConfig(level=settings.LOG_LEVEL, format=settings.LOG_FORMAT)
logger = logging.getLogger(__name__)
