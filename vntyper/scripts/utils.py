#!/usr/bin/env python3
# vntyper/scripts/utils.py

import gzip
import importlib.resources as pkg_resources
import json
import logging
import os
import shlex
import subprocess
import sys


def run_command(command, log_file, critical=False):
    """
    Helper function to run a shell command and log its output.

    Args:
        command (str): The command to run.
        log_file (str): The path to the log file where stdout and stderr will be logged.
        critical (bool): If True, the pipeline will stop if the command fails.

    Returns:
        bool: True if the command succeeded, False otherwise.
    """
    logging.debug(f"Running command: {command}")
    with open(log_file, "w") as lf:
        process = subprocess.Popen(
            command,
            shell=True,
            executable="/bin/bash",  # Ensure Bash is used for process substitution
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        for line in process.stdout:
            decoded_line = line.decode()
            lf.write(decoded_line)
            logging.debug(decoded_line.strip())
        process.wait()

        if process.returncode != 0:
            logging.debug(f"Command failed: {command}")
            if critical:
                raise RuntimeError(f"Critical command failed: {command}")
            return False
    return True


def setup_logging(log_level=logging.INFO, log_file=None):
    """
    Sets up logging for the application.

    Args:
        log_level (int): Logging level (e.g., logging.INFO).
        log_file (str, optional): Path to a log file. If None, logs are printed to console.
    """
    logger = logging.getLogger()  # Get the root logger
    logger.setLevel(log_level)  # Set the overall logging level

    # Clear existing handlers so we don't duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    # If writing logs to a file, create a FileHandler, set its level, and attach
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Always attach a console handler at the same requested log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


def create_output_directories(base_output_dir):
    """
    Creates necessary output directories for the pipeline.

    Args:
        base_output_dir (str): Base directory for output files.

    Returns:
        dict: A dictionary mapping directory names to their paths.
    """
    dirs = {
        "base": base_output_dir,
        "kestrel": os.path.join(base_output_dir, "kestrel"),
        "advntr": os.path.join(base_output_dir, "advntr"),
        "fastq_bam_processing": os.path.join(base_output_dir, "fastq_bam_processing"),
        "alignment_processing": os.path.join(base_output_dir, "alignment_processing"),
        "coverage": os.path.join(base_output_dir, "coverage"),
    }

    for dir_path in dirs.values():
        try:
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
                logging.info(f"Created directory: {dir_path}")
            else:
                logging.info(f"Directory already exists: {dir_path}")
        except Exception as e:
            logging.error(f"Failed to create directory {dir_path}: {e}")
            raise

    return dirs


def get_tool_version(command, version_flag):
    """
    Runs a command to get the version of a tool and returns the parsed version string.

    Args:
        command (str): The command to run (e.g., "fastp").
        version_flag (str): The flag or argument to pass to the command to get its version
            (e.g., "-v").

    Returns:
        str: The parsed version string or 'unknown' if parsing fails.
    """
    try:
        # Split the command properly in case it's a compound command like "mamba run ..."
        full_command = shlex.split(command) + shlex.split(version_flag)
        result = subprocess.run(full_command, capture_output=True, text=True)
        output = result.stdout.strip() or result.stderr.strip()

        # Parse version from the output based on the command
        if command.startswith("fastp"):
            if "fastp" in output:
                return output.split("\n")[1].split(" ")[1]
            return "unknown"
        if command.startswith("samtools"):
            if "samtools" in output:
                return output.split("\n")[1].split(" ")[1]
            return "unknown"
        if command.startswith("bwa"):
            # Capture the second line which contains the version info
            lines = output.split("\n")
            if len(lines) >= 2 and "Version" in lines[1]:
                return lines[1].split(": ")[1]
            return "unknown"
        if "advntr" in command:
            lines = output.split("\n")
            if len(lines) >= 3 and "adVNTR" in lines[2]:
                return lines[2].split(": ")[0].split(" ")[1]
            return "unknown"
        if "java" in command and "kestrel" in command:
            # Handle Kestrel version parsing (Java + JAR execution)
            if "kestrel" in output:
                return output.split("\n")[-1].split(": ")[1]
            return "unknown"
        if command.startswith("java"):  # Handling java_path case
            return output.split("\n")[0]  # Return the first line of the Java version output
        return "unknown"

    except FileNotFoundError:
        logging.error(f"Command not found: {command}")
        return "unknown"
    except PermissionError:
        logging.error(f"Permission denied: {command}")
        return "unknown"
    except IndexError as e:
        logging.error(f"Failed to parse version for {command}: {e}")
        return "unknown"
    except Exception as e:
        logging.error(f"Failed to get version for {command}: {e}")
        return "unknown"


def get_tool_versions(config):
    """
    Retrieves the versions of the tools specified in the config and returns them as a
    dictionary.

    Args:
        config (dict): The configuration dictionary.

    Returns:
        dict: A dictionary with tool names as keys and their version strings as values.
    """
    tools = config.get("tools", {})
    versions = {}

    # Define version commands for each tool
    version_commands = {
        "fastp": "",
        "samtools": "",
        "bwa": "",
        "advntr": "",
        "java_path": "--version",
        "kestrel": '-jar "{kestrel_path}" -h'.format(
            kestrel_path=tools.get("kestrel", ""),
        ),
    }

    for tool, command in tools.items():
        version_flag = version_commands.get(tool, "")
        # Special handling for kestrel as it needs the java_path in front
        if tool == "kestrel":
            command = f"{tools.get('java_path', 'java')} {version_flag}"
        versions[tool] = get_tool_version(command, version_flag)

    return versions


def search(regex: str, df, case=False):
    """
    Searches for a regex pattern in a DataFrame.

    Args:
        regex (str): The regex pattern to search for.
        df (DataFrame): The pandas DataFrame to search.
        case (bool): Whether the search should be case-sensitive.

    Returns:
        DataFrame: A DataFrame containing the rows where the pattern was found.
    """
    logging.debug("Starting regex search in DataFrame.")
    try:
        textlikes = df.select_dtypes(include=[object, "object"])
        result_df = df[
            textlikes.apply(lambda column: column.str.contains(regex, regex=True, case=case, na=False)).any(axis=1)
        ]
        logging.debug("Regex search completed.")
        return result_df
    except Exception as e:
        logging.error(f"Error during regex search: {e}")
        raise


def load_config(config_path=None):
    """
    Load the configuration file with fallback to the default package config.

    Args:
        config_path (str or None): Path to the user-provided config file.

    Returns:
        dict: The loaded configuration dictionary.
    """
    if config_path is not None and os.path.exists(config_path):
        # User provided a config path
        try:
            with open(config_path) as config_file:
                config = json.load(config_file)
                logging.info(f"Configuration loaded from {config_path}")
                return config
        except json.JSONDecodeError as e:
            logging.error(f"Error decoding JSON from the config file: {e}")
            raise
        except Exception as e:
            logging.error(f"Unexpected error loading config file {config_path}: {e}")
            raise
    else:
        # No config path provided or file does not exist; use default config from package data
        try:
            with pkg_resources.open_text("vntyper", "config.json") as f:
                config = json.load(f)
                logging.info("Loaded default config from package data.")
                return config
        except Exception as e:
            logging.error("Error: Default config file not found in package data.")
            logging.error(e)
            sys.exit(1)


def validate_bam_file(file_path):
    """
    Validates the alignment file (BAM or CRAM) for existence, correct extension, and
    integrity using samtools quickcheck.

    This function was originally intended for BAM files, but we now extend it
    to handle CRAM as well. The logic remains the same; we just allow .cram
    extension in addition to .bam and run samtools quickcheck regardless.

    Args:
        file_path (str): Path to the BAM or CRAM file.

    Raises:
        ValueError: If any validation check fails.
    """
    if not file_path:
        logging.error("No alignment file provided.")
        raise ValueError("No alignment file provided.")

    if not os.path.isfile(file_path):
        logging.error(f"Alignment file does not exist: {file_path}")
        raise ValueError(f"Alignment file does not exist: {file_path}")

    # Modified to allow both .bam and .cram extensions
    if not (file_path.endswith(".bam") or file_path.endswith(".cram")):
        logging.error(f"Invalid alignment file extension for file: {file_path}. Must be .bam or .cram")
        raise ValueError(f"Invalid alignment file extension for file: {file_path}")

    # Perform samtools quickcheck
    command = f"samtools quickcheck -v {file_path}"
    log_file = f"{file_path}.quickcheck.log"
    success = run_command(command, log_file, critical=True)
    if not success:
        logging.error(f"Alignment file failed quickcheck: {file_path}")
        raise ValueError(f"Alignment file failed quickcheck: {file_path}")

    logging.info(f"Alignment file validated successfully: {file_path}")


def validate_fastq_file(file_path):
    """
    Validates the FASTQ file for existence, correct extension, and basic formatting.

    Args:
        file_path (str): Path to the FASTQ file.

    Raises:
        ValueError: If any validation check fails.
    """
    if not file_path:
        logging.error("No FASTQ file provided.")
        raise ValueError("No FASTQ file provided.")

    if not os.path.isfile(file_path):
        logging.error(f"FASTQ file does not exist: {file_path}")
        raise ValueError(f"FASTQ file does not exist: {file_path}")

    valid_extensions = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    if not file_path.endswith(valid_extensions):
        logging.error(f"Invalid FASTQ file extension for file: {file_path}")
        raise ValueError(f"Invalid FASTQ file extension for file: {file_path}")

    # Check basic FASTQ formatting by reading the first few lines
    try:
        open_func = gzip.open if file_path.endswith(".gz") else open
        with open_func(file_path, "rt") as f:
            for _ in range(4):  # Read first 4 lines of the first read
                line = f.readline()
                if not line:
                    logging.error(f"FASTQ file is incomplete or empty: {file_path}")
                    raise ValueError(f"FASTQ file is incomplete or empty: {file_path}")
        logging.info(f"FASTQ file validated successfully: {file_path}")
    except Exception as e:
        logging.error(f"Error validating FASTQ file {file_path}: {e}")
        raise
