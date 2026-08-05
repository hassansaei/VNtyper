# docker/app/celery_app.py

from celery import Celery
from celery.schedules import crontab
from celery.signals import beat_init, celeryd_init
from kombu import Queue
import os
import logging

from .config import build_redis_url, get_redis_password, require_redis_password

# Configure logging for Celery
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB = int(os.getenv("REDIS_DB", 0))  # Default Celery DB

# Broker and result-backend URL. The credential comes from the single accessor in
# config.py - the same one app/main.py and app/tasks.py use - and is
# percent-encoded by build_redis_url, because a Redis password may legitimately
# contain URL delimiters. check_redis_password below refuses to start a worker
# when the variable is unset.
REDIS_URL = build_redis_url(REDIS_HOST, REDIS_PORT, REDIS_DB, get_redis_password())

# Initialize Celery
celery_app = Celery("worker", broker=REDIS_URL, backend=REDIS_URL)

# Autodiscover tasks from the 'app.tasks' module
celery_app.autodiscover_tasks(["app.tasks"])

# Configure Celery Beat schedule for periodic tasks
celery_app.conf.beat_schedule = {
    "delete-old-results-daily": {
        "task": "app.tasks.delete_old_results",
        "schedule": crontab(hour=0, minute=0),  # Runs daily at midnight UTC
    },
}
celery_app.conf.timezone = "UTC"

# Define task routes to direct 'run_vntyper_job' to 'vntyper_queue'
celery_app.conf.task_routes = {
    "app.tasks.run_vntyper_job": {"queue": "vntyper_queue"},
}

# Define task queues
celery_app.conf.task_queues = (
    Queue("vntyper_queue"),
    # ------------------------------------------------
    # ADDED: Queue for long-running advntr jobs
    # ------------------------------------------------
    Queue("vntyper_long_queue"),
    Queue("celery"),  # Default queue for other tasks
)

# Remove any rate limits for Celery tasks (rate limiting is handled in FastAPI)
celery_app.conf.task_annotations = {
    "app.tasks.run_vntyper_job": {"rate_limit": None},
}

# Optional: Log Redis connection details (excluding password for security)
logger.info(f"Celery Broker URL: redis://{REDIS_HOST}:{REDIS_PORT}/{REDIS_DB}")


@celeryd_init.connect
@beat_init.connect
def check_redis_password(**_kwargs) -> None:
    """Refuse to start a worker or beat process without a Redis credential.

    Wired to the Celery startup signals rather than run at import time: the API
    and the test suite both import this module, and neither should depend on the
    deployment environment merely to be imported. `celeryd_init` fires from
    `Worker.on_before_init`, before the worker touches the broker.

    The failure is raised as `SystemExit`, not `RuntimeError`: Celery's
    `Signal.send` catches `Exception` from its receivers and only logs it, so a
    plain raise here would be swallowed and the process would start anyway.

    Args:
        **_kwargs: Celery signal payload, unused.

    Raises:
        SystemExit: If REDIS_PASSWORD is unset or empty. Exit code 1.
    """
    try:
        require_redis_password()
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
