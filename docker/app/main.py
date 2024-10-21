# docker/app/main.py

from fastapi import FastAPI, UploadFile, File, Form, BackgroundTasks
from fastapi.responses import FileResponse
from .tasks import run_vntyper_job
import os

app = FastAPI()

# Environment variables for default directories
DEFAULT_INPUT_DIR = os.getenv("DEFAULT_INPUT_DIR", "/opt/vntyper/input")
DEFAULT_OUTPUT_DIR = os.getenv("DEFAULT_OUTPUT_DIR", "/opt/vntyper/output")

@app.post("/run-job/")
async def run_vntyper(
    file: UploadFile = File(...),
    thread: int = Form(4),
    reference_assembly: str = Form("hg38"),
    fast_mode: bool = Form(False),
    keep_intermediates: bool = Form(False),
    archive_results: bool = Form(False),
    background_tasks: BackgroundTasks = None
):
    """
    Endpoint to run VNtyper job with additional parameters.
    """
    # Ensure input and output directories exist
    os.makedirs(DEFAULT_INPUT_DIR, exist_ok=True)
    os.makedirs(DEFAULT_OUTPUT_DIR, exist_ok=True)
    
    # Save the uploaded file
    file_path = os.path.join(DEFAULT_INPUT_DIR, file.filename)
    with open(file_path, "wb") as f:
        f.write(await file.read())
    
    # Define output path
    output_dir = os.path.join(DEFAULT_OUTPUT_DIR, os.path.splitext(file.filename)[0])
    os.makedirs(output_dir, exist_ok=True)
    
    # Add background task
    background_tasks.add_task(
        run_vntyper_job,
        bam_path=file_path,
        output_dir=output_dir,
        thread=thread,
        reference_assembly=reference_assembly,
        fast_mode=fast_mode,
        keep_intermediates=keep_intermediates,
        archive_results=archive_results
    )
    
    return {"message": "Job started", "output_dir": f"/download/{os.path.basename(output_dir)}"}
    
@app.get("/download/{output_dir}")
async def download_result(output_dir: str):
    """
    Endpoint to download the result (zipped file).
    """
    zip_path = os.path.join(DEFAULT_OUTPUT_DIR, f"{output_dir}.zip")
    if os.path.exists(zip_path):
        return FileResponse(zip_path, media_type="application/zip", filename=f"{output_dir}.zip")
    return {"error": "File not found"}
