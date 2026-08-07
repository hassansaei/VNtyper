# vntyper/modules/shark/shark_filtering.py

import logging
import os
from pathlib import Path

from vntyper.scripts.command_builders import quote_path
from vntyper.scripts.utils import load_config, run_command

logger = logging.getLogger(__name__)


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
        reference_assembly (str): Accepted for API compatibility and **ignored**.
            SHARK matches k-mers against a single MUC1 region FASTA and does not
            select a region by coordinate, so the assembly does not change what
            it retains. Passing anything other than hg19/GRCh37 logs a warning.
            See issue #187.
        threads (int): Number of threads to use for SHARK.

    Returns:
        tuple: (filtered_fastq_1, filtered_fastq_2)

    Raises:
        ValueError: If muc1_region_fasta not defined in shark_config.json.
        RuntimeError: If SHARK filtering fails.
    """
    if reference_assembly and reference_assembly.lower() not in ("hg19", "grch37"):
        logger.warning(
            f"reference_assembly={reference_assembly!r} does not select a region for SHARK. "
            "SHARK is sequence-based, not coordinate-based: it filters against the single "
            "MUC1 region FASTA regardless of assembly. The parameter is retained for API "
            "compatibility only. See issue #187."
        )

    # shark_path should come from the main config since it contains the tool paths.
    shark_path = main_config.get("tools", {}).get("shark", "shark")

    # Get muc1_region_fasta from the shark_config
    muc1_region_fasta = config.get("shark_settings", {}).get("muc1_region_fasta")
    if not muc1_region_fasta:
        raise ValueError("muc1_region_fasta not defined in shark_config.json")

    filtered_fastq_1 = os.path.join(output_dir, f"{sample_name}_shark_R1.fastq")
    filtered_fastq_2 = os.path.join(output_dir, f"{sample_name}_shark_R2.fastq")

    # `run_command` runs this as one string under bash (trap 9), so quoting can only
    # happen here. Every operand is a path or a thread count and is quoted;
    # `shark_path` is not, because config.json holds a command *prefix* there --
    # "mamba run -n shark_env shark" (trap 6) -- which quoting would collapse into a
    # single token bash looks for as one binary.
    command = (
        f"{shark_path} -r {quote_path(muc1_region_fasta)} "
        f"-1 {quote_path(fastq_1)} -2 {quote_path(fastq_2)} "
        f"-o {quote_path(filtered_fastq_1)} -p {quote_path(filtered_fastq_2)} "
        f"-t {quote_path(threads)}"
    )

    log_file = Path(output_dir) / f"{sample_name}_shark.log"
    logger.info(f"Running SHARK filtering with command: {command}")

    success = run_command(command, str(log_file), critical=True)
    if not success:
        logger.error("SHARK filtering failed.")
        raise RuntimeError("SHARK filtering failed.")

    return filtered_fastq_1, filtered_fastq_2
