# docker/app/tasks.py

import hashlib
import os
import shutil
import subprocess
from datetime import datetime, timedelta, timezone

import redis
from celery.utils.log import get_task_logger

from .celery_app import celery_app
from .cohorts import cohort_key, extend_cohort_retention
from .config import get_redis_password, settings
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
    """Name the index this job will use.

    `samtools index` writes `.crai` beside a CRAM and `.bai` beside a BAM, so
    the fallback has to be chosen by format (#188). Naming `.bai` for a CRAM
    named a file that is never created: the existence check never found the
    index the worker had itself just built, and cleanup then removed nothing
    while the real `.crai` stayed on the volume every job shares.

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
    """Name every index this job's own tooling can put beside its alignment.

    The worker is not the only thing that indexes the submission. Non-fast BAM
    processing in `vntyper/scripts/fastq_bam_processing.py` reconstructs the
    index name as `f"{in_bam}.bai"` and builds it whenever that exact path is
    missing -- it never sees an index the client uploaded as `sample.bai`. So a
    submission that carried `sample.bai` still ends up with `sample.bam.bai`
    beside it, under a name neither the client nor the worker ever mentioned,
    and cleanup that removes only the submitted paths leaves it on the volume
    every job shares.

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
    one. The endpoint accepts several index names, so the worker is told which
    one it got rather than guessing: guessing means rebuilding an index the job
    was already given, under a different name, and leaving the supplied file on
    the shared volume afterwards. With no index supplied it falls back to the
    conventional name beside the alignment, which is also where `samtools index`
    puts the one it builds.
    """
    job_id = os.path.basename(output_dir)
    # Bound before the try block: the cleanup below runs whether or not the task
    # got as far as its first Redis call, and must remove exactly the index this
    # job used -- uploaded or generated -- rather than raise a NameError of its
    # own on the way.
    index_path = resolve_index_path(bam_path, index_path)
    try:
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

        # Ensure the alignment has an index the pipeline can find
        if not os.path.exists(index_path):
            logger.info(f"Index not found for {bam_path}. Generating index.")
            try:
                subprocess.run(["samtools", "index", bam_path], check=True)
                logger.info(f"Successfully generated index at {index_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error generating index: {e}")
                # Update usage data on failure
                redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
                raise

        # Build the base command for VNtyper
        command = build_vntyper_command(
            alignment_path=bam_path,
            output_dir=output_dir,
            thread=thread,
            reference_assembly=reference_assembly,
            fast_mode=fast_mode,
            keep_intermediates=keep_intermediates,
            archive_results=archive_results,
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

        # Optionally, archive results
        if archive_results:
            try:
                shutil.make_archive(output_dir, "zip", output_dir)
                shutil.rmtree(output_dir)
                logger.info(f"Archived results to {output_dir}.zip and removed original directory")
            except Exception as e:
                logger.error(f"Error archiving results: {e}")
                # Update usage data on failure
                redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
                raise

        # Update usage data on success
        redis_usage_client.hset(f"usage:{job_id}", "status", "completed")

        # **Removed**: Cohort assignment from the task
        # if cohort_key:
        #     redis_cohort_client.sadd(f"{cohort_key}:jobs", job_id)

        # **Cohort assignment is now handled only upon successful completion**
        if cohort_key:
            redis_cohort_client.sadd(f"{cohort_key}:jobs", job_id)
            logger.info(f"Assigned job {job_id} to cohort {cohort_key}")

        # Construct the download URL
        download_url = f"{settings.API_BASE_URL}/api/download/{job_id}/"

        # Send success email if provided
        if email:
            subject = "VNtyper Job Completed Successfully"
            if cohort_key:
                cohort_id = cohort_key.split(":", 1)[1]
                content = f"""
                <p>Your VNtyper job has been completed successfully.</p>
                <p>Job ID: <strong>{job_id}</strong></p>
                <p>Cohort ID: <strong>{cohort_id}</strong></p>
                <p>You can download your results <a href="{download_url}">here</a>.</p>
                """
            else:
                content = f"""
                <p>Your VNtyper job has been completed successfully.</p>
                <p>Job ID: <strong>{job_id}</strong></p>
                <p>You can download your results <a href="{download_url}">here</a>.</p>
                """
            send_email_task.delay(to_email=email, subject=subject, content=content)
            logger.info(f"Triggered email sending to {email} with download link")

    except Exception as e:
        logger.error(f"Error in VNtyper job: {e}")
        # Update usage data on failure
        redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
        raise
    finally:
        # Patient-derived files come off the shared volume before anything else
        # in this block, and every bookkeeping call after them is isolated.
        # Both halves matter: this block used to open with two unguarded Redis
        # calls, so a broker that became unreachable as the task exited left the
        # alignment, its index and the whole input directory behind -- and
        # replaced the pipeline's own exception with a connection error, hiding
        # why the job failed.
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
                extend_cohort_retention(redis_cohort_client, cohort_key, settings.COHORT_RETENTION_DAYS * 86400)
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
        # 1) Create directory, input file listing all .zip files
        os.makedirs(output_dir, exist_ok=True)
        input_file = os.path.join(output_dir, "cohort_input.txt")
        with open(input_file, "w") as f:
            f.writelines(f"{zpath}\n" for zpath in zip_paths)

        # 2) Run the "vntyper cohort" command
        command = [
            "conda",
            "run",
            "-n",
            "vntyper",
            "vntyper",
            "cohort",
            "--input-file",
            input_file,
            "-o",
            output_dir,
        ]
        logger.info(f"Running command: {' '.join(command)}")
        subprocess.run(command, check=True)
        logger.info("Joint cohort analysis completed.")

        # 3) Zip the results
        try:
            shutil.make_archive(output_dir, "zip", output_dir)
            logger.info(f"Zipped results to {output_dir}.zip")
        except Exception as e:
            logger.error(f"Error zipping results for cohort analysis: {e}")
            redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
            raise

        # Update usage data on success
        redis_usage_client.hset(f"usage:{job_id}", "status", "completed")

    except Exception as e:
        logger.error(f"Error in joint cohort analysis for {cohort_id}: {e}")
        redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
        raise
    finally:
        # Remove the task ID from the Redis list
        redis_client.lrem("vntyper_job_queue", 0, task_id)
        logger.info(f"Removed cohort analysis task ID {task_id} from queue")

        # Extend cohort TTL if necessary
        if cohort_id:
            extend_cohort_retention(
                redis_cohort_client,
                cohort_key(cohort_id),
                settings.COHORT_RETENTION_DAYS * 86400,
            )

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
