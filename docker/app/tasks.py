# docker/app/tasks.py
import subprocess
import os

def run_vntyper_job(input_file: str, output_zip: str):
    """
    Function that runs the VNtyper command and zips the results.
    """
    output_dir = os.path.dirname(output_zip)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Run the VNtyper command
    vntyper_command = [
        "conda", "run", "-n", "vntyper", "vntyper",
        "--input", input_file,
        "--output", output_dir
    ]
    subprocess.run(vntyper_command, check=True)

    # Zip the output directory
    subprocess.run(["zip", "-r", output_zip, output_dir], check=True)
