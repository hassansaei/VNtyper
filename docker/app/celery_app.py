# docker/app/celery_app.py

from celery import Celery
from celery.schedules import crontab
from kombu import Queue
import os

# Initialize Celery
celery_app = Celery(
    "worker",
    broker=os.getenv("CELERY_BROKER_URL", "redis://redis:6379/0"),
    backend=os.getenv("CELERY_RESULT_BACKEND", "redis://redis:6379/0")
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
    Queue('celery'),  # Default queue for other tasks
)

# Remove any rate limits for Celery tasks (rate limiting is handled in FastAPI)
celery_app.conf.task_annotations = {
    'app.tasks.run_vntyper_job': {'rate_limit': None},
}
