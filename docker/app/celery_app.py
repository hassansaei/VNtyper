# docker/app/celery_app.py

from celery import Celery
import os

# Initialize Celery
celery_app = Celery(
    "worker",
    broker=os.getenv("CELERY_BROKER_URL", "redis://redis:6379/0"),
    backend=os.getenv("CELERY_RESULT_BACKEND", "redis://redis:6379/0")
)

# Autodiscover tasks from the 'app.tasks' module
celery_app.autodiscover_tasks(["app.tasks"])
