# docker/app/main.py

from fastapi import FastAPI, UploadFile, File, Form
from fastapi.responses import FileResponse
from uuid import uuid4
import os
from .tasks import run_vntyper_job
import shutil
import logging

from .config import settings

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=settings.PROJECT_VERSION
)

logger = logging.getLogger(__name__)

@app.post("/run-job/")
async def run_vntyper(
    file: UploadFile = File(...),
    thread: int = Form(4),
    reference_assembly: str = Form("hg38"),
    fast_mode: bool = Form(False),
    keep_intermediates: bool = Form(False),
    archive_results: bool = Form(False)
):
    """
    Endpoint to run VNtyper job with additional parameters.
    """
    logger.info("Received job submission")
    # Ensure input and output directories exist
    os.makedirs(DEFAULT_INPUT_DIR, exist_ok=True)
    os.makedirs(DEFAULT_OUTPUT_DIR, exist_ok=True)

    # Generate a unique job ID
    job_id = str(uuid4())
    job_input_dir = os.path.join(DEFAULT_INPUT_DIR, job_id)
    job_output_dir = os.path.join(DEFAULT_OUTPUT_DIR, job_id)

    os.makedirs(job_input_dir, exist_ok=True)
    os.makedirs(job_output_dir, exist_ok=True)

    # Save the uploaded file
    file_path = os.path.join(job_input_dir, file.filename)
    with open(file_path, "wb") as f:
        shutil.copyfileobj(file.file, f)
    logger.info(f"Saved uploaded file to {file_path}")

    # Enqueue the Celery task
    run_vntyper_job.delay(
        bam_path=file_path,
        output_dir=job_output_dir,
        thread=thread,
        reference_assembly=reference_assembly,
        fast_mode=fast_mode,
        keep_intermediates=keep_intermediates,
        archive_results=archive_results
    )
    logger.info(f"Enqueued job {job_id}")

    return {"message": "Job submitted", "job_id": job_id}

@app.get("/download/{job_id}")
def download_result(job_id: str):
    """
    Endpoint to download the result (zipped file).
    """
    logger.info(f"Received download request for job {job_id}")
    zip_path = os.path.join(DEFAULT_OUTPUT_DIR, f"{job_id}.zip")
    if os.path.exists(zip_path):
        logger.info(f"Serving zip file {zip_path}")
        return FileResponse(
            zip_path,
            media_type="application/zip",
            filename=f"{job_id}.zip"
        )
    logger.warning(f"File not found: {zip_path}")
    return {"error": "File not found"}
