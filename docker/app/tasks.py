# docker/app/tasks.py

import hashlib
import os
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

from . import admission
from .archive_delivery import snapshot_owned_archives
from .celery_app import celery_app
from .cohort_workspace import cohort_workspace
from .cohorts import cohort_key, extend_cohort_retention
from .config import get_redis_password, settings
from .failure_reporting import report_pipeline_failure
from .job_failures import clear_preflight_failure
from .pipeline_job_workspace import (
    PipelineJobWorkspace,
    open_pipeline_job_workspace,
    reclaim_unopened_spool_inputs,
)
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
    alignment_is_cram: bool | None = None,
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
        "--cram" if (is_cram(alignment_path) if alignment_is_cram is None else alignment_is_cram) else "--bam",
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
    capacity_reserved: bool = False,
    workspace_identity: dict[str, object] | None = None,
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
    heartbeat: admission.ReservationHeartbeat | None = None
    workspace: PipelineJobWorkspace | None = None
    input_is_cram = is_cram(bam_path)
    try:
        workspace = open_pipeline_job_workspace(
            bam_path,
            index_path,
            output_dir,
            workspace_identity,
            spool_root=settings.DEFAULT_HANDOFF_SPOOL_DIR,
            shared_roots=(
                settings.DEFAULT_INPUT_DIR,
                settings.DEFAULT_OUTPUT_DIR,
                settings.DEFAULT_HANDOFF_SPOOL_DIR,
            ),
        )
        if capacity_reserved:
            heartbeat = admission.ReservationHeartbeat(
                redis_client,
                job_id,
                active_lease_seconds=settings.ADMISSION_ACTIVE_LEASE_SECONDS,
                heartbeat_seconds=settings.ADMISSION_HEARTBEAT_SECONDS,
            )
            heartbeat.start()
        index_path = resolve_index_path(bam_path, index_path)
        job_key = f"usage:{job_id}"
        redis_usage_client.hdel(job_key, "code", "message")
        clear_preflight_failure(workspace.bound_output_path)
        logger.info(f"Starting VNtyper job for BAM file: {workspace.bound_alignment_path}")

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
            alignment_path=workspace.bound_alignment_path,
            output_dir=workspace.bound_output_path,
            thread=thread,
            reference_assembly=reference_assembly,
            fast_mode=fast_mode,
            keep_intermediates=keep_intermediates,
            archive_results=False,
            advntr_mode=advntr_mode,
            alignment_is_cram=input_is_cram,
        )

        # Run the VNtyper pipeline
        try:
            subprocess.run(command, check=True)
            workspace.require_current_output()
            workspace.remove_views()
            logger.info(f"VNtyper job completed for {workspace.bound_alignment_path}")
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            logger.error(f"Error running VNtyper: {e}")
            report_pipeline_failure(
                redis_usage_client=redis_usage_client,
                job_id=job_id,
                output_dir=workspace.bound_output_path,
                email=email,
                cohort_key=cohort_key,
                error=e,
                send_email_task=send_email_task,
                logger=logger,
            )
            raise

        if archive_results:
            create_safe_archive(
                output_dir,
                "zip",
                workspace.bound_output_path,
                protected_paths=(workspace.bound_alignment_path,)
                if workspace.bound_index_path is None
                else (workspace.bound_alignment_path, workspace.bound_index_path),
                root_descriptor=workspace.output_descriptor,
            )
            archive_published = True
            try:
                workspace.remove_output()
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
        if not archive_published:
            workspace.require_current_output("completion bookkeeping")

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

        if not archive_published:
            workspace.require_current_output("completion bookkeeping")

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
        if heartbeat is not None:
            heartbeat.stop()
        # Research inputs leave the protected handoff spool before bookkeeping.
        if workspace is not None:
            for cleanup_name, cleanup in (
                ("private input views", workspace.remove_views),
                ("uploaded inputs", workspace.remove_inputs),
                ("detached output", workspace.discard_detached_output),
                ("workspace descriptors", workspace.close),
            ):
                try:
                    cleanup()
                except Exception as cleanup_error:
                    logger.error(f"Error cleaning {cleanup_name} for {job_id}: {cleanup_error}")
        else:
            try:
                reclaim_unopened_spool_inputs(
                    bam_path,
                    index_path,
                    output_dir,
                    settings.DEFAULT_HANDOFF_SPOOL_DIR,
                )
            except Exception as cleanup_error:
                logger.error(f"Error cleaning unopened uploaded inputs for {job_id}: {cleanup_error}")

        # Bookkeeping last, and best-effort. Neither call owns anything on
        # disk: `lrem` maintains the queue-position display, and the retention
        # extension only pushes out a TTL that `delete_old_results()` would
        # otherwise act on. Raising from either would mask the task's own
        # exception, and failing a job whose pipeline completed and whose
        # results are already written reports a success as a failure.
        try:
            admission.release_worker_bookkeeping(
                redis_client, self.request.id, job_id, f"task ID {self.request.id}", logger
            )
        except Exception as e:
            logger.error(f"Error releasing worker bookkeeping for {job_id}: {e}")
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
    capacity_reserved: bool = False,
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
    heartbeat: admission.ReservationHeartbeat | None = None
    try:
        if capacity_reserved:
            heartbeat = admission.ReservationHeartbeat(
                redis_client,
                job_id,
                active_lease_seconds=settings.ADMISSION_ACTIVE_LEASE_SECONDS,
                heartbeat_seconds=settings.ADMISSION_HEARTBEAT_SECONDS,
            )
            heartbeat.start()
        logger.info(f"Starting joint cohort analysis for Cohort ID: {cohort_id}")

        user_data = f"{user_ip}-{user_agent}"
        usage_data = {
            "user_hash": hashlib.sha256(user_data.encode("utf-8")).hexdigest(),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "job_id": job_id,
            "status": "started",
            "analysis_type": "cohort_analysis",
            "cohort_id": cohort_id,
        }
        redis_usage_client.hset(f"usage:{job_id}", mapping=usage_data)
        redis_usage_client.expire(f"usage:{job_id}", settings.USAGE_DATA_RETENTION_SECONDS)

        # 1) Exclusively reserve this task's public output name before removing
        # an earlier archive. Every task write then stays in private staging,
        # anchored through its open descriptor.
        with cohort_workspace(output_dir) as workspace:
            clear_stale_archive(output_dir, "zip", protected_paths=zip_paths)
            try:
                snapshot_dir = os.path.join(workspace.staging_path, ".cohort-members")
                with snapshot_owned_archives(
                    zip_paths,
                    snapshot_dir,
                    parent_descriptor=workspace.descriptor,
                ) as snapshots:
                    cohort_input_file = workspace.child_path("cohort_input.txt")
                    input_file = cohort_input_file
                    input_descriptor = os.open(
                        "cohort_input.txt",
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
                        0o600,
                        dir_fd=workspace.descriptor,
                    )
                    with os.fdopen(input_descriptor, "w") as f:
                        f.writelines(f"{zpath}\n" for zpath in snapshots.paths)

                    # 2) Run the "vntyper cohort" command while the descriptors
                    # anchoring every member snapshot remain open.
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
                        workspace.bound_path,
                    ]
                    logger.info(f"Running command: {' '.join(command)}")
                    subprocess.run(command, check=True)
                    logger.info("Joint cohort analysis completed.")

                workspace.unlink("cohort_input.txt")
                input_file = None

                # 3) Zip the results from the descriptor-bound directory, not from
                # whatever may now occupy its original public name.
                create_safe_archive(
                    output_dir,
                    "zip",
                    workspace.staging_path,
                    protected_paths=zip_paths,
                    root_descriptor=workspace.descriptor,
                )
                archive_published = True
                logger.info(f"Zipped results to {output_dir}.zip")
            finally:
                if input_file is not None:
                    try:
                        workspace.unlink("cohort_input.txt")
                        logger.info(f"Deleted cohort input file: {input_file}")
                    except Exception as cleanup_error:
                        logger.error(f"Error deleting cohort input file {input_file}: {cleanup_error}")

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
        if heartbeat is not None:
            heartbeat.stop()
        try:
            admission.release_worker_bookkeeping(
                redis_client, task_id, job_id, f"cohort analysis task ID {task_id}", logger
            )
        except Exception as e:
            logger.error(f"Error releasing worker bookkeeping for {job_id}: {e}")

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
