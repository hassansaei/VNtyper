# docker/app/celery_app.py

from celery import Celery
from celery.schedules import crontab
from kombu import Queue
import os
import logging

# Configure logging for Celery
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Retrieve Redis password from environment variables
# Set a default password if REDIS_PASSWORD is not provided
DEFAULT_REDIS_PASSWORD = "qE3!#zjraRG*`X2g4%<x&J"
REDIS_PASSWORD = os.getenv("REDIS_PASSWORD", DEFAULT_REDIS_PASSWORD)

# Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB = int(os.getenv("REDIS_DB", 0))  # Default Celery DB

# Construct Redis URL with password
REDIS_URL = f"redis://:{REDIS_PASSWORD}@{REDIS_HOST}:{REDIS_PORT}/{REDIS_DB}"

# Initialize Celery
celery_app = Celery(
    "worker",
    broker=REDIS_URL,
    backend=REDIS_URL
)

# Autodiscover tasks from the 'app.tasks' module
celery_app.autodiscover_tasks(["app.tasks"])

# Configure Celery Beat schedule for periodic tasks
celery_app.conf.beat_schedule = {
    'delete-old-results-daily': {
        'task': 'app.tasks.delete_old_results',
        'schedule': crontab(hour=0, minute=0),  # Runs daily at midnight UTC
    },
}
celery_app.conf.timezone = 'UTC'

# Define task routes to direct 'run_vntyper_job' to 'vntyper_queue'
celery_app.conf.task_routes = {
    'app.tasks.run_vntyper_job': {'queue': 'vntyper_queue'},
}

# Define task queues
celery_app.conf.task_queues = (
    Queue('vntyper_queue'),
    # ------------------------------------------------
    # ADDED: Queue for long-running advntr jobs
    # ------------------------------------------------
    Queue('vntyper_long_queue'),
    Queue('celery'),  # Default queue for other tasks
)

# Remove any rate limits for Celery tasks (rate limiting is handled in FastAPI)
celery_app.conf.task_annotations = {
    'app.tasks.run_vntyper_job': {'rate_limit': None},
}

# Optional: Log Redis connection details (excluding password for security)
logger.info(f"Celery Broker URL: redis://{REDIS_HOST}:{REDIS_PORT}/{REDIS_DB}")
