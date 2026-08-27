"""Path-free HTML content for pipeline failure notifications."""

from __future__ import annotations

import html
import subprocess


def build_failure_email_content(
    job_id: str,
    cohort_id: str | None,
    error: subprocess.CalledProcessError | subprocess.TimeoutExpired,
) -> str:
    """Render escaped, recipient-actionable pipeline failure details.

    Command arguments and worker paths remain private diagnostics in the server
    log. A recipient only needs the identifiers and the failure class.

    Args:
        job_id: Public identifier of the failed job.
        cohort_id: Public cohort identifier, or ``None`` outside a cohort.
        error: Process exit or deadline failure raised by the pipeline runner.

    Returns:
        HTML containing escaped identifiers and a reconstructed failure reason.
    """
    if isinstance(error, subprocess.TimeoutExpired):
        failure_reason = f"The analysis did not finish within {error.timeout:g} seconds and was stopped."
    else:
        failure_reason = f"The analysis exited with status {error.returncode}."

    lines = [
        f"<p>Your VNtyper job with Job ID <strong>{html.escape(job_id)}</strong> has failed.</p>",
    ]
    if cohort_id is not None:
        lines.append(f"<p>Cohort ID: <strong>{html.escape(cohort_id)}</strong></p>")
    lines.append(f"<p>{html.escape(failure_reason)}</p>")
    lines.append("<p>Quote the job ID when reporting this so the failure can be looked up in the server logs.</p>")
    return "\n".join(lines)
