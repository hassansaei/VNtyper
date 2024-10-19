import os
import json
import logging
import subprocess
import shlex
import subprocess as sp

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
    logging.info(f"Running command: {command}")
    with open(log_file, "w") as lf:
        process = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)
        for line in process.stdout:
            lf.write(line.decode())
            logging.info(line.decode().strip())
        process.wait()

        if process.returncode != 0:
            logging.error(f"Command failed: {command}")
            if critical:
                raise RuntimeError(f"Critical command failed: {command}")
            return False
    return True

def setup_logging(log_level=logging.INFO, log_file=None):
    logger = logging.getLogger()  # Get the root logger
    logger.setLevel(log_level)  # Set the overall logging level
    
    if logger.hasHandlers():
        logger.handlers.clear()
    
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

def create_output_directories(base_output_dir):
    dirs = {
        "base": base_output_dir,
        "kestrel": os.path.join(base_output_dir, "kestrel"),
        "advntr": os.path.join(base_output_dir, "advntr"),
        "fastq_bam_processing": os.path.join(base_output_dir, "fastq_bam_processing"),
        "alignment_processing": os.path.join(base_output_dir, "alignment_processing"),
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
        version_flag (str): The flag or argument to pass to the command to get its version (e.g., "-v").
    
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
            return output.split("\n")[1].split(" ")[1] if "fastp" in output else "unknown"
        elif command.startswith("samtools"):
            return output.split("\n")[1].split(" ")[1] if "samtools" in output else "unknown"
        elif command.startswith("bwa"):
            # Capture the third line which contains the version info
            lines = output.split("\n")
            return lines[1].split(": ")[1] if len(lines) >= 3 and "Version" in lines[1] else "unknown"
        elif "advntr" in command:
            # Capture the third line which contains the version info
            lines = output.split("\n")
            return lines[2].split(": ")[0].split(" ")[1] if len(lines) >= 3 and "adVNTR" in lines[2] else "unknown"
        elif "java" in command and "kestrel" in command:
            # Handle Kestrel version parsing (Java + JAR execution)
            # get the last line of the output
            return output.split("\n")[-1].split(": ")[1] if "kestrel" in output else "unknown"
        elif command.startswith("/usr/bin/java"):  # Handling java_path case
            return output.split("\n")[0]  # Return the first line of the Java version output
        else:
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
    Retrieves the versions of the tools specified in the config and returns them as a dictionary.

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
        "kestrel": "-jar \"{kestrel_path}\" -h".format(
            kestrel_path=config["tools"]["kestrel"],
        )
    }

    for tool, command in tools.items():
        version_flag = version_commands.get(tool, "")
        # Special handling for kestrel as it needs the java_path in front
        if tool == "kestrel":
            command = f"{tools['java_path']} {version_flag}"
        versions[tool] = get_tool_version(command, version_flag)

    return versions

def search(regex: str, df, case=False):
    logging.debug("Starting regex search in DataFrame.")
    try:
        textlikes = df.select_dtypes(include=[object, "object"])
        result_df = df[
            textlikes.apply(
                lambda column: column.str.contains(regex, regex=True, case=case, na=False)
            ).any(axis=1)
        ]
        logging.debug("Regex search completed.")
        return result_df
    except Exception as e:
        logging.error(f"Error during regex search: {e}")
        raise

def load_config(config_path=None):
    if config_path is None:
        logging.info("No config file path provided, proceeding without configuration.")
        return None

    if not os.path.exists(config_path):
        logging.error(f"Config file not found: {config_path}")
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        with open(config_path, 'r') as config_file:
            config = json.load(config_file)
            logging.info(f"Configuration loaded from {config_path}")
            return config
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from the config file: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error loading config file {config_path}: {e}")
        raise
