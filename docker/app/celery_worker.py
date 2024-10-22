# docker/app/celery_worker.py

from celery import Celery
import os

# Initialize Celery
celery_app = Celery(
    "worker",
    broker=os.getenv("CELERY_BROKER_URL", "redis://redis:6379/0"),
    backend=os.getenv("CELERY_RESULT_BACKEND", "redis://redis:6379/0")
)

# Import tasks to ensure they're registered with Celery
import app.tasks
