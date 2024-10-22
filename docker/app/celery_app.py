# docker/app/celery_app.py

from celery import Celery
from celery.schedules import crontab
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
