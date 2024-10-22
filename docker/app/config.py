# docker/app/config.py

import os
import logging

class Settings:
    PROJECT_NAME: str = "VNtyper API"
    PROJECT_VERSION: str = "2.0.0-alpha.7"
    DEBUG: bool = True

    # Logging configuration
    LOG_LEVEL: str = os.getenv("LOG_LEVEL", "INFO")
    LOG_FORMAT: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Directory configurations
    DEFAULT_INPUT_DIR: str = os.getenv("DEFAULT_INPUT_DIR", "/opt/vntyper/input")
    DEFAULT_OUTPUT_DIR: str = os.getenv("DEFAULT_OUTPUT_DIR", "/opt/vntyper/output")

settings = Settings()

# Configure logging
logging.basicConfig(level=settings.LOG_LEVEL, format=settings.LOG_FORMAT)
logger = logging.getLogger(__name__)
