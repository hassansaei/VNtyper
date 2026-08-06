# vntyper/scripts/utils.py

import gzip
import importlib.resources as pkg_resources
import json
import logging
import os
import shlex
import subprocess
import sys

from vntyper.scripts.command_builders import quote_path

logger = logging.getLogger(__name__)


def run_command(command, log_file, critical=False, cwd=None):
    """
    Helper function to run a shell command and log its output.

    The command is executed by bash via ``shell=True`` with
    ``executable="/bin/bash"``. That is deliberate and load-bearing: the CRAM
    unmapped-read path builds a command using bash process substitution, which no
    ``shell=False`` argv list can express. Callers are responsible for quoting the
    values they interpolate into ``command``.

    stdout and stderr are merged and streamed to ``log_file`` as the child runs. A
    non-zero exit status is logged at ERROR either way, so a failure is visible at the
    default log level; ``critical`` only decides whether it also aborts the pipeline.

    Args:
        command (str): The command to run, as a single shell string.
        log_file (str): The path to the log file where stdout and stderr will be logged.
        critical (bool): If True, the pipeline will stop if the command fails.
        cwd (str, optional): The working directory to execute the command in.
            If None, the subprocess will inherit the current working directory.
            Setting this explicitly is important for tools like Java that need to
            determine the working directory during initialization, and because every
            tool and reference path in ``config.json`` is relative to it.

    Returns:
        bool: True if the command succeeded, False if it failed and ``critical`` is False.

    Raises:
        RuntimeError: If the command exits non-zero and ``critical`` is True.
    """
    logger.debug(f"Running command: {command}")
    if cwd is not None:
        logger.debug(f"Working directory: {cwd}")

    with open(log_file, "w") as lf:
        process = subprocess.Popen(
            command,
            shell=True,
            executable="/bin/bash",  # Ensure Bash is used for process substitution
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=cwd,  # Explicitly set the working directory (or None to use inherited)
        )
        for line in process.stdout:
            decoded_line = line.decode()
            lf.write(decoded_line)
            logger.debug(decoded_line.strip())
        process.wait()

        if process.returncode != 0:
            # ERROR, not DEBUG: a non-critical failure returns False and the pipeline
            # carries on, so this line is the only signal the stage did not work. At
            # the default INFO level a DEBUG log here produced no output at all.
            logger.error(f"Command failed: {command}")
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
    root_logger = logging.getLogger()  # Get the root logger
    root_logger.setLevel(log_level)  # Set the overall logging level

    # Clear existing handlers so we don't duplicate logs
    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    # If writing logs to a file, create a FileHandler, set its level, and attach
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    # Always attach a console handler at the same requested log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)


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
                logger.info(f"Created directory: {dir_path}")
            else:
                logger.info(f"Directory already exists: {dir_path}")
        except Exception as e:
            logger.error(f"Failed to create directory {dir_path}: {e}")
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
        result = subprocess.run(full_command, capture_output=True, text=True, check=False)
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
        logger.error(f"Command not found: {command}")
        return "unknown"
    except PermissionError:
        logger.error(f"Permission denied: {command}")
        return "unknown"
    except IndexError as e:
        logger.error(f"Failed to parse version for {command}: {e}")
        return "unknown"
    except Exception as e:
        logger.error(f"Failed to get version for {command}: {e}")
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
        # Special handling for kestrel as it needs the java_path in front.
        # The flag is folded into ``command`` here, so it must not ALSO be passed to
        # get_tool_version: that builds its argv as
        # ``shlex.split(command) + shlex.split(version_flag)``, which used to hand
        # subprocess.run ``java -jar <jar> -h -jar <jar> -h``. Clearing version_flag,
        # rather than dropping the pre-built command, is deliberate: get_tool_version
        # selects its kestrel parse on ``"java" in command and "kestrel" in command``,
        # so ``command`` has to keep carrying both substrings.
        if tool == "kestrel":
            command = f"{tools.get('java_path', 'java')} {version_flag}"
            version_flag = ""
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
    logger.debug("Starting regex search in DataFrame.")
    try:
        textlikes = df.select_dtypes(include=[object, "object"])
        result_df = df[
            textlikes.apply(lambda column: column.str.contains(regex, regex=True, case=case, na=False)).any(axis=1)
        ]
        logger.debug("Regex search completed.")
        return result_df
    except Exception as e:
        logger.error(f"Error during regex search: {e}")
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
                logger.info(f"Configuration loaded from {config_path}")
                return config
        except json.JSONDecodeError as e:
            logger.error(f"Error decoding JSON from the config file: {e}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error loading config file {config_path}: {e}")
            raise
    else:
        # No config path provided or file does not exist; use default config from package data
        try:
            with pkg_resources.open_text("vntyper", "config.json") as f:
                config = json.load(f)
                logger.info("Loaded default config from package data.")
                return config
        except Exception as e:
            logger.error("Error: Default config file not found in package data.")
            logger.error(e)
            sys.exit(1)


def validate_bam_file(file_path, cwd=None):
    """
    Validates the alignment file (BAM or CRAM) for existence, correct extension, and
    integrity using samtools quickcheck.

    This function was originally intended for BAM files, but we now extend it
    to handle CRAM as well. The logic remains the same; we just allow .cram
    extension in addition to .bam and run samtools quickcheck regardless.

    Args:
        file_path (str): Path to the BAM or CRAM file.
        cwd (str, optional): Working directory to use when running samtools.
            Important for environments where the working directory may become invalid.

    Raises:
        ValueError: If any validation check fails.
    """
    if not file_path:
        logger.error("No alignment file provided.")
        raise ValueError("No alignment file provided.")

    if not os.path.isfile(file_path):
        logger.error(f"Alignment file does not exist: {file_path}")
        raise ValueError(f"Alignment file does not exist: {file_path}")

    # Modified to allow both .bam and .cram extensions
    if not (file_path.endswith((".bam", ".cram"))):
        logger.error(f"Invalid alignment file extension for file: {file_path}. Must be .bam or .cram")
        raise ValueError(f"Invalid alignment file extension for file: {file_path}")

    # Perform samtools quickcheck
    # Quoted because `run_command` runs with shell=True: an unquoted path containing a
    # space becomes two operands, and one containing a metacharacter is executed.
    #
    # critical=False, not True: `run_command`'s critical path raises RuntimeError and
    # never returns False, which made the `if not success:` branch below unreachable and
    # the documented `Raises: ValueError` untrue on a real quickcheck failure. This
    # function aborts the run itself, by raising, so `run_command` does not also need
    # to -- and the exception callers see is the ValueError every other check here (and
    # `validate_fastq_file`) raises for unusable input.
    command = f"samtools quickcheck -v {quote_path(file_path)}"
    log_file = f"{file_path}.quickcheck.log"
    success = run_command(command, log_file, critical=False, cwd=cwd)
    if not success:
        logger.error(f"Alignment file failed quickcheck: {file_path}")
        raise ValueError(f"Alignment file failed quickcheck: {file_path}")

    logger.info(f"Alignment file validated successfully: {file_path}")


def validate_fastq_file(file_path):
    """
    Validates the FASTQ file for existence, correct extension, and basic formatting.

    Args:
        file_path (str): Path to the FASTQ file.

    Raises:
        ValueError: If any validation check fails.
    """
    if not file_path:
        logger.error("No FASTQ file provided.")
        raise ValueError("No FASTQ file provided.")

    if not os.path.isfile(file_path):
        logger.error(f"FASTQ file does not exist: {file_path}")
        raise ValueError(f"FASTQ file does not exist: {file_path}")

    valid_extensions = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    if not file_path.endswith(valid_extensions):
        logger.error(f"Invalid FASTQ file extension for file: {file_path}")
        raise ValueError(f"Invalid FASTQ file extension for file: {file_path}")

    # Check basic FASTQ formatting by reading the first few lines
    try:
        open_func = gzip.open if file_path.endswith(".gz") else open
        with open_func(file_path, "rt") as f:
            for _ in range(4):  # Read first 4 lines of the first read
                line = f.readline()
                if not line:
                    logger.error(f"FASTQ file is incomplete or empty: {file_path}")
                    raise ValueError(f"FASTQ file is incomplete or empty: {file_path}")
        logger.info(f"FASTQ file validated successfully: {file_path}")
    except Exception as e:
        logger.error(f"Error validating FASTQ file {file_path}: {e}")
        raise
