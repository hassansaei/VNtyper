"""Best-effort bookkeeping for a failed pipeline subprocess."""

from __future__ import annotations

import logging
import subprocess
from typing import Any

from .failure_email import build_failure_email_content
from .job_failures import read_preflight_failure


def report_pipeline_failure(
    *,
    redis_usage_client: Any,
    job_id: str,
    output_dir: str,
    email: str | None,
    cohort_key: str | None,
    error: subprocess.CalledProcessError | subprocess.TimeoutExpired,
    send_email_task: Any,
    logger: logging.Logger,
) -> None:
    """Record and notify about ``error`` without raising a secondary failure.

    Args:
        redis_usage_client: Redis client that owns per-job usage hashes.
        job_id: Public identifier of the failed job.
        output_dir: Job output directory that may hold curated preflight metadata.
        email: Notification recipient, or ``None`` when no email was requested.
        cohort_key: Internal cohort Redis key, or ``None`` outside a cohort.
        error: The pipeline subprocess exception that must remain primary.
        send_email_task: Celery task used to enqueue notification delivery.
        logger: Worker logger for secondary-failure diagnostics.
    """
    job_key = f"usage:{job_id}"
    try:
        redis_usage_client.hset(job_key, "status", "failed")
    except Exception as status_error:
        logger.error(f"Error recording failed status for {job_id}: {status_error}")

    try:
        preflight_failure = read_preflight_failure(output_dir)
        if preflight_failure is not None:
            redis_usage_client.hset(job_key, mapping=preflight_failure)
    except Exception as preflight_error:
        logger.error(f"Error recording preflight failure for {job_id}: {preflight_error}")

    if email:
        try:
            cohort_id = cohort_key.split(":", 1)[1] if cohort_key else None
            content = build_failure_email_content(job_id, cohort_id, error)
            send_email_task.delay(to_email=email, subject="VNtyper Job Failed", content=content)
        except Exception as email_error:
            logger.error(f"Error queueing failure email for {job_id}: {email_error}")
