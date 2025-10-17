#!/usr/bin/env python3
# vntyper/modules/shark/shark_filtering.py

import logging
import os
from pathlib import Path

from vntyper.scripts.utils import load_config, run_command


def load_shark_config(config_path=None):
    """
    Load the SHARK configuration file.
    Similar to load_advntr_config, this function will load shark_config.json
    by default.
    """
    if config_path is None:
        # Default path to shark_config.json
        config_path = os.path.join(os.path.dirname(__file__), "shark_config.json")
    return load_config(config_path)


# Load the shark settings (if needed globally)
shark_config = load_shark_config()
shark_settings = shark_config.get("shark_settings", {})


def run_shark_filter(
    fastq_1,
    fastq_2,
    output_dir,
    config,  # shark_config dict
    main_config,  # main config dict (for tool paths)
    sample_name,
    reference_assembly="hg19",
    threads=4,
):
    """
    Run SHARK filtering on FASTQ files.

    Args:
        fastq_1 (str): Path to input FASTQ R1.
        fastq_2 (str): Path to input FASTQ R2.
        output_dir (str): Directory to store SHARK output.
        config (dict): The shark_config dictionary (contains muc1_region_fasta).
        main_config (dict): The main configuration dictionary (contains tool paths).
        sample_name (str): Sample name.
        reference_assembly (str): Reference assembly (hg19 or hg38).
        threads (int): Number of threads to use for SHARK.

    Returns:
        tuple: (filtered_fastq_1, filtered_fastq_2)

    Raises:
        ValueError: If muc1_region_fasta not defined in shark_config.json.
        RuntimeError: If SHARK filtering fails.
    """
    # shark_path should come from the main config since it contains the tool paths.
    shark_path = main_config.get("tools", {}).get("shark", "shark")

    # Get muc1_region_fasta from the shark_config
    muc1_region_fasta = config.get("shark_settings", {}).get("muc1_region_fasta")
    if not muc1_region_fasta:
        raise ValueError("muc1_region_fasta not defined in shark_config.json")

    filtered_fastq_1 = os.path.join(output_dir, f"{sample_name}_shark_R1.fastq")
    filtered_fastq_2 = os.path.join(output_dir, f"{sample_name}_shark_R2.fastq")

    # Add the threads parameter to the command using the -t option
    command = (
        f"{shark_path} -r {muc1_region_fasta} -1 {fastq_1} -2 {fastq_2} "
        f"-o {filtered_fastq_1} -p {filtered_fastq_2} -t {threads}"
    )

    log_file = Path(output_dir) / f"{sample_name}_shark.log"
    logging.info(f"Running SHARK filtering with command: {command}")

    success = run_command(command, str(log_file), critical=True)
    if not success:
        logging.error("SHARK filtering failed.")
        raise RuntimeError("SHARK filtering failed.")

    return filtered_fastq_1, filtered_fastq_2
