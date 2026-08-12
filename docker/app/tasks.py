# docker/app/tasks.py

import hashlib
import os
import shutil
import subprocess
from datetime import datetime, timedelta, timezone

import redis
from celery.utils.log import get_task_logger

from vntyper.scripts.archive_safety import (
    clear_stale_archive,
    create_safe_archive,
    quarantine_archive,
    revoke_public_archive,
)

from .archive_delivery import snapshot_owned_archives
from .celery_app import celery_app
from .cohorts import cohort_key, extend_cohort_retention
from .config import get_redis_password, settings
from .job_failures import clear_preflight_failure, read_preflight_failure
from .utils import send_email

logger = get_task_logger(__name__)

# Environment variables for Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

# Redis credential, from the single accessor in config.py so the worker resolves
# the same value as the API. There is no fallback; app/celery_app.py refuses to
# start a worker when the variable is unset.
REDIS_PASSWORD = get_redis_password()

# Redis DBs
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Job mappings
COHORT_REDIS_DB = int(os.getenv("COHORT_REDIS_DB", 3))  # Cohort data
USAGE_REDIS_DB = settings.USAGE_REDIS_DB  # Usage statistics

# Initialize Redis clients with authentication
redis_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True,
)
redis_cohort_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=COHORT_REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True,
)
redis_usage_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=USAGE_REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True,
)


def is_cram(alignment_path: str) -> bool:
    """Report whether a stored alignment is a CRAM.

    Matched without regard to case: the upload allowlist in `uploads.py`
    compiles its pattern with `re.IGNORECASE`, so `SAMPLE.CRAM` is a name the
    endpoint stores and enqueues, and a case-sensitive test here would send
    exactly those submissions back down the BAM path.

    Args:
        alignment_path: The stored alignment.

    Returns:
        bool: True for a CRAM, False for a BAM.
    """
    return str(alignment_path).lower().endswith(".cram")


def resolve_index_path(alignment_path: str, index_path: str | None) -> str:
    """Name the submitted or conventional index this job will clean up.

    The pipeline now builds missing indexes only in its run-local output view.
    The fallback still has to be format-aware so cleanup covers conventional
    index names left by older jobs without inventing a BAM suffix for CRAM.

    Args:
        alignment_path: The stored alignment.
        index_path: The index the submission carried, or None.

    Returns:
        str: The supplied index unchanged, or the conventional name beside the
            alignment when the submission carried none.
    """
    if index_path:
        return index_path
    return f"{alignment_path}.crai" if is_cram(alignment_path) else f"{alignment_path}.bai"


def derived_index_paths(alignment_path: str) -> tuple[str, ...]:
    """Name every conventional index entry cleanup may find beside an alignment.

    Current pipeline preflight never creates these in the input tree, but jobs
    may carry either accepted spelling and older workers could have generated
    one. Cleanup remains deliberately bounded to these exact names.

    The names are derived from the alignment rather than discovered, and every
    one of them is joined back onto the alignment's own directory, so this can
    only ever name files inside the per-job input directory. A glob would be
    shorter and would also match files this job did not create.

    Args:
        alignment_path: The stored alignment.

    Returns:
        tuple[str, ...]: The `<alignment>.bai`, `<alignment>.crai`,
            `<stem>.bai` and `<stem>.crai` paths beside the alignment, in that
            order and without repeats -- an alignment whose name carries no
            extension collapses the two forms onto the same pair.
    """
    directory, filename = os.path.split(alignment_path)
    stem = os.path.splitext(filename)[0]

    names: list[str] = []
    for name in (f"{filename}.bai", f"{filename}.crai", f"{stem}.bai", f"{stem}.crai"):
        if name not in names:
            names.append(name)
    return tuple(os.path.join(directory, name) for name in names)


def remove_job_input_files(alignment_path: str, index_path: str) -> None:
    """Delete the alignment and every index belonging to one job.

    Both are patient-derived and sit on the volume every job shares, so this is
    the only thing that ever removes them. Each removal is isolated: one that
    fails is logged and the rest still run, and nothing raises out of here,
    because this is called from a `finally` where an exception would replace
    whatever the task was already failing with.

    Args:
        alignment_path: The stored alignment.
        index_path: The index this job used, as resolved by
            `resolve_index_path`.
    """
    paths: list[str] = []
    for path in (alignment_path, index_path, *derived_index_paths(alignment_path)):
        if path not in paths:
            paths.append(path)

    for path in paths:
        try:
            if os.path.exists(path):
                os.remove(path)
                logger.info(f"Deleted input file: {path}")
        except Exception as e:
            logger.error(f"Error deleting input file {path}: {e}")


def build_vntyper_command(
    alignment_path: str,
    output_dir: str,
    thread: int,
    reference_assembly: str,
    fast_mode: bool = False,
    keep_intermediates: bool = False,
    archive_results: bool = False,
    advntr_mode: bool = False,
) -> list[str]:
    """Assemble the vntyper CLI invocation for one job.

    Extracted from run_vntyper_job so the flag selection is testable without a
    Celery worker. Behaviour is unchanged apart from the alignment flag (#188):
    the endpoint has always accepted `.cram`, but the command hardcoded `--bam`,
    so every accepted CRAM was handed to the CLI as a BAM and took the BAM code
    path.

    Args:
        alignment_path: The stored alignment.
        output_dir: The per-job output directory.
        thread: Threads to give the pipeline.
        reference_assembly: The assembly the alignment is against.
        fast_mode: Whether to append ``--fast-mode``.
        keep_intermediates: Whether to append ``--keep-intermediates``.
        archive_results: Whether to append ``--archive-results``.
        advntr_mode: Whether to include the adVNTR stage.

    Returns:
        list[str]: The argument vector to run.
    """
    command = [
        "conda",
        "run",
        # Without this, conda buffers the pipeline's stdout and stderr until it exits, so
        # a web job that runs for an hour logs nothing until the hour is up and a job that
        # hangs logs nothing at all (#213). `docker/entrypoint.sh` carries the same flag;
        # this is the second invocation, and it is the one every *web* job goes through.
        # It precedes `-n` because `conda run` parses its own options before the
        # environment selector.
        "--no-capture-output",
        "-n",
        "vntyper",
        "vntyper",
        "pipeline",
        "--cram" if is_cram(alignment_path) else "--bam",
        alignment_path,
        "-o",
        output_dir,
        "--threads",
        str(thread),
        "--reference-assembly",
        reference_assembly,
    ]

    if fast_mode:
        command.append("--fast-mode")
    if keep_intermediates:
        command.append("--keep-intermediates")
    if archive_results:
        command.append("--archive-results")
    if advntr_mode:
        command.extend(["--extra-modules", "advntr", "--advntr-max-coverage", "300"])

    return command


@celery_app.task(bind=True, max_retries=3, default_retry_delay=60)
def send_email_task(self, to_email: str, subject: str, content: str):
    """
    Celery task to send an email via SMTP.
    Retries up to 3 times in case of failure.
    """
    try:
        send_email(to_email=to_email, subject=subject, content=content)
        logger.info(f"Email sent to {to_email} with subject '{subject}'")
    except Exception as e:
        logger.error(f"Failed to send email to {to_email}: {e}")
        raise self.retry(exc=e) from e


@celery_app.task(bind=True)
def run_vntyper_job(
    self,
    bam_path: str,
    output_dir: str,
    thread: int,
    reference_assembly: str,
    fast_mode: bool,
    keep_intermediates: bool,
    archive_results: bool,
    email: str | None = None,
    cohort_key: str | None = None,
    client_ip: str | None = None,
    user_agent: str | None = None,
    advntr_mode: bool = False,
    index_path: str | None = None,
):
    """
    Celery task to run VNtyper pipeline with parameters.
    Sends an email upon job completion or failure if email is provided.

    `index_path` is where the submission's own index was stored, when it carried
    one. The endpoint accepts several index names, so cleanup uses that exact
    path rather than guessing. Missing-index construction belongs exclusively
    to pipeline preflight, which builds a run-local index under `output_dir`.
    """
    job_id = os.path.basename(output_dir)
    archive_published = False
    # Bound before the try block: the cleanup below runs whether or not the task
    # got as far as its first Redis call, and must remove exactly the index this
    # job used -- uploaded or generated -- rather than raise a NameError of its
    # own on the way.
    index_path = resolve_index_path(bam_path, index_path)
    try:
        job_key = f"usage:{job_id}"
        redis_usage_client.hdel(job_key, "code", "message")
        clear_preflight_failure(output_dir)
        logger.info(f"Starting VNtyper job for BAM file: {bam_path}")

        # Generate a unique hash for the user
        user_data = f"{client_ip}-{user_agent}"
        user_hash = hashlib.sha256(user_data.encode("utf-8")).hexdigest()

        # Store initial usage data
        usage_data = {
            "user_hash": user_hash,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "job_id": job_id,
            "status": "started",
        }
        redis_usage_client.hset(f"usage:{job_id}", mapping=usage_data)
        redis_usage_client.expire(f"usage:{job_id}", settings.USAGE_DATA_RETENTION_SECONDS)

        if archive_results:
            clear_stale_archive(
                output_dir, "zip", protected_paths=(bam_path, index_path) if index_path else (bam_path,)
            )

        # Build the base command for VNtyper
        command = build_vntyper_command(
            alignment_path=bam_path,
            output_dir=output_dir,
            thread=thread,
            reference_assembly=reference_assembly,
            fast_mode=fast_mode,
            keep_intermediates=keep_intermediates,
            archive_results=False,
            advntr_mode=advntr_mode,
        )

        # Run the VNtyper pipeline
        try:
            subprocess.run(command, check=True)
            logger.info(f"VNtyper job completed for {bam_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running VNtyper: {e}")
            # Update usage data on failure
            redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
            preflight_failure = read_preflight_failure(output_dir)
            if preflight_failure is not None:
                redis_usage_client.hset(f"usage:{job_id}", mapping=preflight_failure)
            # Send failure email if provided
            if email:
                subject = "VNtyper Job Failed"
                if cohort_key:
                    cohort_id = cohort_key.split(":", 1)[1]
                    content = f"""
                    <p>Your VNtyper job with Job ID <strong>{job_id}</strong> has failed.</p>
                    <p>Cohort ID: <strong>{cohort_id}</strong></p>
                    <p>Error Details:</p>
                    <pre>{str(e)}</pre>
                    """
                else:
                    content = f"""
                    <p>Your VNtyper job with Job ID <strong>{job_id}</strong> has failed.</p>
                    <p>Error Details:</p>
                    <pre>{str(e)}</pre>
                    """
                send_email_task.delay(to_email=email, subject=subject, content=content)
            raise

        if archive_results:
            create_safe_archive(
                output_dir,
                "zip",
                output_dir,
                protected_paths=(bam_path, index_path) if index_path else (bam_path,),
            )
            archive_published = True
            try:
                shutil.rmtree(output_dir)
                logger.info(f"Archived results to {output_dir}.zip and removed original directory")
            except Exception as e:
                archive_published = False
                try:
                    quarantine_path = quarantine_archive(
                        output_dir,
                        "zip",
                        protected_paths=(bam_path, index_path) if index_path else (bam_path,),
                    )
                    logger.error(f"Result cleanup failed; preserved complete archive at {quarantine_path}: {e}")
                except Exception as rollback_error:
                    logger.error(f"Error quarantining failed job's public archive: {rollback_error}")
                    archive_published = True
                raise

        # Update usage data on success
        redis_usage_client.hset(f"usage:{job_id}", "status", "completed")

        if cohort_key:
            redis_cohort_client.sadd(f"{cohort_key}:jobs", job_id)
            logger.info(f"Assigned job {job_id} to cohort {cohort_key}")

        # Construct the download URL
        download_url = f"{settings.API_BASE_URL}/api/download/{job_id}/"

        # Send success email if provided
        if email:
            subject = "VNtyper Job Completed Successfully"
            # The email used to state no deadline at all. That was already an omission; with a
            # shorter window it becomes a misleading one, because the recipient has no way to
            # know the results may be gone before they next look.
            #
            # It says "about N days", not "stops working after N days", because the latter
            # would be false. `MAX_RESULT_AGE_DAYS` is a *cleanup eligibility* threshold, not
            # an enforced deadline: `delete_old_results` runs once daily at 00:00 UTC
            # (celery_app.py) and removes an archive only once its timestamp is past the
            # cutoff, and `/download/{job_id}/` performs no age check at all -- it serves
            # whatever still exists. A job finishing at 00:01 is therefore retrievable for
            # nearly N+1 days, and a deletion that errors leaves the file in place. Promising
            # a hard deadline the service does not enforce would be a new false claim rather
            # than a fix for an old one.
            #
            # It also names the property of the link that matters: it is a capability, not an
            # authenticated route (#189).
            retention_days = settings.MAX_RESULT_AGE_DAYS
            if cohort_key:
                cohort_id = cohort_key.split(":", 1)[1]
                content = f"""
                <p>Your VNtyper job has been completed successfully.</p>
                <p>Job ID: <strong>{job_id}</strong></p>
                <p>Cohort ID: <strong>{cohort_id}</strong></p>
                <p>You can download your results <a href="{download_url}">here</a>.</p>
                <p>Results are deleted from the server about {retention_days} days after
                completion, so download them before then. Until they are deleted, anyone with
                this link can retrieve them -- it carries no password, so treat it as
                confidential.</p>
                """
            else:
                content = f"""
                <p>Your VNtyper job has been completed successfully.</p>
                <p>Job ID: <strong>{job_id}</strong></p>
                <p>You can download your results <a href="{download_url}">here</a>.</p>
                <p>Results are deleted from the server about {retention_days} days after
                completion, so download them before then. Until they are deleted, anyone with
                this link can retrieve them -- it carries no password, so treat it as
                confidential.</p>
                """
            send_email_task.delay(to_email=email, subject=subject, content=content)
            logger.info(f"Triggered email sending to {email} with download link")

    except Exception as e:
        logger.error(f"Error in VNtyper job: {e}")
        if archive_published:
            try:
                quarantine_path = quarantine_archive(
                    output_dir,
                    "zip",
                    protected_paths=(bam_path, index_path) if index_path else (bam_path,),
                )
                archive_published = False
                if quarantine_path is not None:
                    logger.error(f"Preserved failed job's complete archive at {quarantine_path}")
            except Exception as quarantine_error:
                logger.error(f"Error quarantining failed job's public archive: {quarantine_error}")
                try:
                    revoke_public_archive(
                        output_dir,
                        "zip",
                        protected_paths=(bam_path, index_path) if index_path else (bam_path,),
                    )
                    archive_published = False
                except Exception as rollback_error:
                    logger.error(f"Error removing failed job's public archive: {rollback_error}")
        try:
            redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
        except Exception as status_error:
            logger.error(f"Error recording failed status for {job_id}: {status_error}")
        raise
    finally:
        # Patient-derived inputs come off the shared volume before bookkeeping.
        remove_job_input_files(bam_path, index_path)

        # The per-job input directory holds nothing but this submission's own
        # files, so it goes with them. Removing it only when it is empty means
        # anything unexpected left in there is reported rather than deleted.
        job_input_dir = os.path.dirname(bam_path)
        try:
            if os.path.isdir(job_input_dir):
                if os.listdir(job_input_dir):
                    logger.warning(f"Input directory {job_input_dir} still holds files and was left in place")
                else:
                    os.rmdir(job_input_dir)
                    logger.info(f"Deleted empty input directory: {job_input_dir}")
        except Exception as e:
            logger.error(f"Error deleting input directory {job_input_dir}: {e}")

        # Bookkeeping last, and best-effort. Neither call owns anything on
        # disk: `lrem` maintains the queue-position display, and the retention
        # extension only pushes out a TTL that `delete_old_results()` would
        # otherwise act on. Raising from either would mask the task's own
        # exception, and failing a job whose pipeline completed and whose
        # results are already written reports a success as a failure.
        try:
            # Remove the task ID from the Redis list
            redis_client.lrem("vntyper_job_queue", 0, self.request.id)
            logger.info(f"Removed task ID {self.request.id} from vntyper_job_queue")
        except Exception as e:
            logger.error(f"Error removing task ID {self.request.id} from vntyper_job_queue: {e}")

        # Extend cohort TTL if necessary
        if cohort_key:
            try:
                extend_cohort_retention(redis_cohort_client, cohort_key, settings.cohort_retention_days() * 86400)
            except Exception as e:
                logger.error(f"Error extending retention for cohort {cohort_key}: {e}")


@celery_app.task
def delete_old_results():
    """
    Celery task to delete result ZIP files older than a specified age.
    Also deletes cohorts that have expired.
    """
    max_age_days = settings.MAX_RESULT_AGE_DAYS
    output_dir = settings.DEFAULT_OUTPUT_DIR
    # Both sides of the age comparison below are UTC-aware. They must stay that way
    # together: mixing an aware datetime with a naive one raises TypeError, and the
    # instants themselves are unchanged by the timezone, so the comparison result is
    # the same as it was with the two naive local-time values.
    cutoff_time = datetime.now(timezone.utc) - timedelta(days=max_age_days)

    logger.info(f"Running delete_old_results task. Deleting files older than {max_age_days} days.")

    for filename in os.listdir(output_dir):
        if filename.endswith(".zip"):
            file_path = os.path.join(output_dir, filename)
            if os.path.isfile(file_path):
                file_creation_time = datetime.fromtimestamp(os.path.getctime(file_path), tz=timezone.utc)
                if file_creation_time < cutoff_time:
                    try:
                        os.remove(file_path)
                        logger.info(f"Deleted old result file: {file_path}")
                    except Exception as e:
                        logger.error(f"Error deleting file {file_path}: {e}")

    # Delete expired cohorts
    for key in redis_cohort_client.scan_iter("cohort:*"):
        if not redis_cohort_client.exists(key):
            # Cohort has expired
            cohort_jobs_key = f"{key}:jobs"
            job_ids = redis_cohort_client.smembers(cohort_jobs_key)
            for job_id in job_ids:
                # Delete associated job results
                zip_path = os.path.join(output_dir, f"{job_id}.zip")
                if os.path.exists(zip_path):
                    try:
                        os.remove(zip_path)
                        logger.info(f"Deleted old result file: {zip_path}")
                    except Exception as e:
                        logger.error(f"Error deleting file {zip_path}: {e}")
                # Remove job ID from Redis
                redis_client.delete(job_id)
                redis_client.delete(f"celery-task-meta-{job_id}")  # Remove Celery task meta
            # Remove cohort jobs set
            redis_cohort_client.delete(cohort_jobs_key)
            logger.info(f"Deleted expired cohort data: {key}")

    logger.info("Completed delete_old_results task.")


# ----------------------------------------------------------------------
# Feature #82: Joint Cohort Analysis Task
# ----------------------------------------------------------------------
@celery_app.task(bind=True)
def run_cohort_analysis_job(
    self,
    cohort_id: str,
    zip_paths: list[str],
    output_dir: str,
    user_ip: str | None = None,
    user_agent: str | None = None,
):
    """
    Celery task to run a joint cohort analysis using 'vntyper cohort'.
    zip_paths: list of full paths to individual .zip result files for cohort samples
    output_dir: directory to store the joint analysis output and final zip.

    This task:
      1) Creates a file listing all .zip paths.
      2) Runs 'vntyper cohort --input-file <list> -o <output_dir>'.
      3) Zips the combined results and updates usage data.
    """
    task_id = self.request.id
    job_id = os.path.basename(output_dir)
    # Bound before the try block: the cleanup below runs whether or not the task
    # got as far as naming its scratch file, and must not mask the original
    # failure with a NameError of its own.
    input_file = None
    archive_published = False
    logger.info(f"Starting joint cohort analysis for Cohort ID: {cohort_id}")

    # Generate a unique hash for the user
    user_data = f"{user_ip}-{user_agent}"
    user_hash = hashlib.sha256(user_data.encode("utf-8")).hexdigest()

    # Store initial usage data
    usage_data = {
        "user_hash": user_hash,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "job_id": job_id,
        "status": "started",
        "analysis_type": "cohort_analysis",
        "cohort_id": cohort_id,
    }
    redis_usage_client.hset(f"usage:{job_id}", mapping=usage_data)
    redis_usage_client.expire(f"usage:{job_id}", settings.USAGE_DATA_RETENTION_SECONDS)

    try:
        clear_stale_archive(output_dir, "zip", protected_paths=zip_paths)
        # 1) Create directory, input file listing all .zip files
        os.makedirs(output_dir, exist_ok=True)
        snapshot_dir = os.path.join(output_dir, ".cohort-members")
        with snapshot_owned_archives(zip_paths, snapshot_dir) as snapshots:
            cohort_input_file = os.path.join(output_dir, "cohort_input.txt")
            input_file = cohort_input_file
            with open(cohort_input_file, "w") as f:
                f.writelines(f"{zpath}\n" for zpath in snapshots.paths)

            # 2) Run the "vntyper cohort" command
            command = [
                "conda",
                "run",
                # Same reason as build_vntyper_command: without this, conda buffers the child
                # until it exits, so a cohort analysis logs nothing while it runs (#213).
                "--no-capture-output",
                "-n",
                "vntyper",
                "vntyper",
                "cohort",
                "--input-file",
                cohort_input_file,
                "-o",
                output_dir,
            ]
            logger.info(f"Running command: {' '.join(command)}")
            subprocess.run(command, check=True)
            logger.info("Joint cohort analysis completed.")

        os.remove(cohort_input_file)
        input_file = None

        # 3) Zip the results
        create_safe_archive(output_dir, "zip", output_dir, protected_paths=zip_paths)
        archive_published = True
        logger.info(f"Zipped results to {output_dir}.zip")

        # Update usage data on success
        redis_usage_client.hset(f"usage:{job_id}", "status", "completed")

    except Exception as e:
        logger.error(f"Error in joint cohort analysis for {cohort_id}: {e}")
        if archive_published:
            try:
                quarantine_path = quarantine_archive(output_dir, "zip", protected_paths=zip_paths)
                archive_published = False
                if quarantine_path is not None:
                    logger.error(f"Preserved failed cohort's complete archive at {quarantine_path}")
            except Exception as quarantine_error:
                logger.error(f"Error quarantining failed cohort archive: {quarantine_error}")
                try:
                    revoke_public_archive(output_dir, "zip", protected_paths=zip_paths)
                    archive_published = False
                except Exception as rollback_error:
                    logger.error(f"Error removing failed cohort archive: {rollback_error}")
        try:
            redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
        except Exception as status_error:
            logger.error(f"Error recording failed status for cohort analysis {job_id}: {status_error}")
        raise
    finally:
        # Remove the task ID from the Redis list
        try:
            redis_client.lrem("vntyper_job_queue", 0, task_id)
            logger.info(f"Removed cohort analysis task ID {task_id} from queue")
        except Exception as e:
            logger.error(f"Error removing cohort analysis task ID {task_id} from queue: {e}")

        # Extend cohort TTL if necessary
        if cohort_id:
            try:
                extend_cohort_retention(
                    redis_cohort_client,
                    cohort_key(cohort_id),
                    settings.cohort_retention_days() * 86400,
                )
            except Exception as e:
                logger.error(f"Error extending retention for cohort {cohort_id}: {e}")

        # Delete the listing file this task wrote for itself. That file is the
        # only thing here the task owns; the .zip paths it names are the
        # members' own results, which /download/{job_id}/ serves and which a
        # later analysis reads again, so they are left alone. Result archives
        # are aged out centrally by delete_old_results().
        try:
            if input_file and os.path.exists(input_file):
                os.remove(input_file)
                logger.info(f"Deleted cohort input file: {input_file}")
        except Exception as e:
            logger.error(f"Error deleting cohort input file {input_file}: {e}")

        # Optionally, delete the output directory if it's empty
        try:
            if os.path.exists(output_dir) and not os.listdir(output_dir):
                os.rmdir(output_dir)
                logger.info(f"Deleted empty output directory: {output_dir}")
        except Exception as e:
            logger.error(f"Error deleting directory {output_dir}: {e}")
