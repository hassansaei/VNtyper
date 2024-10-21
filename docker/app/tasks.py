# docker/app/tasks.py

import subprocess
import os
import shutil

def run_vntyper_job(bam_path: str, output_dir: str, thread: int, reference_assembly: str,
                   fast_mode: bool, keep_intermediates: bool, archive_results: bool):
    """
    Function to run VNtyper pipeline with parameters.
    """
    command = [
        "conda", "run", "-n", "vntyper", "vntyper", "pipeline",
        "--bam", bam_path,
        "-o", output_dir,
        "--thread", str(thread),
        "--reference-assembly", reference_assembly
    ]
    
    if fast_mode:
        command.append("--fast-mode")
    if keep_intermediates:
        command.append("--keep-intermediates")
    if archive_results:
        command.append("--archive-results")
    
    # Run the command
    subprocess.run(command, check=True)
    
    # Optionally, archive results
    if archive_results:
        shutil.make_archive(output_dir, 'zip', output_dir)
        shutil.rmtree(output_dir)
