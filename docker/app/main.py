# docker/app/main.py
import os
from fastapi import FastAPI, UploadFile, File, BackgroundTasks
from fastapi.responses import FileResponse
from .tasks import run_vntyper_job

app = FastAPI()

# Use environment variables for default input/output paths
DEFAULT_INPUT_DIR = os.getenv("DEFAULT_INPUT_DIR", "/opt/vntyper/input")
DEFAULT_OUTPUT_DIR = os.getenv("DEFAULT_OUTPUT_DIR", "/opt/vntyper/output")

@app.post("/run-job/")
async def run_vntyper(file: UploadFile = File(...), background_tasks: BackgroundTasks = None):
    """
    Endpoint to run VNtyper job.
    Takes an input file, runs VNtyper in the background, and returns a zipped results file.
    """
    file_path = f"{DEFAULT_INPUT_DIR}/{file.filename}"
    
    # Save the uploaded file
    with open(file_path, "wb") as f:
        f.write(file.file.read())

    # Run the VNtyper job in the background
    output_zip = f"{DEFAULT_OUTPUT_DIR}/{file.filename}.zip"
    background_tasks.add_task(run_vntyper_job, file_path, output_zip)

    return {"message": "Job started", "output_zip": output_zip}

@app.get("/download/{file_name}")
async def download_result(file_name: str):
    """
    Endpoint to download the result (zipped file).
    """
    file_path = f"{DEFAULT_OUTPUT_DIR}/{file_name}"
    
    return FileResponse(file_path, media_type="application/zip", filename=file_name)
