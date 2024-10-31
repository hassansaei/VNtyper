# vntyper/scripts/install_references.py

import json
import logging
import sys
from pathlib import Path
from typing import Optional, Dict, Any, List
from urllib.request import urlretrieve
import tarfile
import zipfile
import subprocess
import hashlib
from datetime import datetime


def load_install_config(config_path: Path) -> Dict[str, Any]:
    """
    Load the installation configuration from a JSON file.

    Args:
        config_path (Path): Path to the installation configuration JSON file.

    Returns:
        dict: Parsed configuration dictionary.

    Raises:
        SystemExit: If the configuration file is missing or malformed.
    """
    if not config_path.exists():
        logging.error(f"Installation config file not found at {config_path}")
        sys.exit(1)

    try:
        with config_path.open('r') as f:
            config = json.load(f)
        return config
    except json.JSONDecodeError as e:
        logging.error(f"Error parsing JSON config: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unexpected error reading config: {e}")
        sys.exit(1)


def download_file(url: str, dest_path: Path):
    """
    Download a file from a URL to the specified destination path.

    Args:
        url (str): URL to download the file from.
        dest_path (Path): Local path to save the downloaded file.

    Raises:
        SystemExit: If the download fails.
    """
    if dest_path.exists():
        logging.info(f"File already exists at {dest_path}. Skipping download.")
        return

    logging.info(f"Downloading from {url} to {dest_path}...")
    try:
        # Ensure parent directories exist
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(url, dest_path)
        logging.info(f"Successfully downloaded {dest_path.name}")
    except Exception as e:
        logging.error(f"Failed to download {url}: {e}")
        sys.exit(1)


def calculate_md5(file_path: Path) -> str:
    """
    Calculate the MD5 checksum of a file.

    Args:
        file_path (Path): Path to the file.

    Returns:
        str: MD5 checksum in hexadecimal format.

    Raises:
        SystemExit: If reading the file fails.
    """
    logging.debug(f"Calculating MD5 for {file_path}...")
    hash_md5 = hashlib.md5()
    try:
        with file_path.open("rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        md5_checksum = hash_md5.hexdigest()
        logging.debug(f"MD5 for {file_path} is {md5_checksum}")
        return md5_checksum
    except Exception as e:
        logging.error(f"Failed to calculate MD5 for {file_path}: {e}")
        sys.exit(1)


def execute_index_command(index_command: str, fasta_path: Path):
    """
    Execute the indexing command for a FASTA file.

    Args:
        index_command (str): The indexing command with a placeholder for the file path.
        fasta_path (Path): Path to the FASTA file to index.

    Raises:
        SystemExit: If the indexing fails.
    """
    command = index_command.format(path=str(fasta_path))
    logging.info(f"Executing indexing command: {command}")
    try:
        # Split the command into arguments
        args = command.split()
        result = subprocess.run(
            args,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        logging.info(f"Successfully executed: {command}")
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Indexing command failed for {fasta_path}: {e.stderr.decode().strip()}"
        )
        sys.exit(1)


def update_config(config_path: Path, references: Dict[str, Path]):
    """
    Update the main config.json with paths to the downloaded references.

    Args:
        config_path (Path): Path to the main config.json file.
        references (dict): Dictionary mapping reference keys to their paths.

    Raises:
        SystemExit: If updating the config fails.
    """
    if not config_path.exists():
        logging.error(
            f"Main config file {config_path} does not exist. Cannot update references."
        )
        sys.exit(1)

    try:
        with config_path.open('r') as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        logging.error(f"Error parsing main config.json: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unexpected error reading main config.json: {e}")
        sys.exit(1)

    # Ensure 'reference_data' section exists
    if "reference_data" not in config:
        config["reference_data"] = {}

    # Update the reference_data section
    for ref_key, ref_path in references.items():
        config["reference_data"][ref_key] = str(ref_path)

    try:
        with config_path.open('w') as f:
            json.dump(config, f, indent=2)
        logging.info(f"Successfully updated {config_path} with new reference paths.")
    except Exception as e:
        logging.error(f"Failed to write updated config.json: {e}")
        sys.exit(1)


def process_ucsc_references(ucsc_refs: Dict[str, Dict[str, str]],
                            output_dir: Path, bwa_path: str, skip_indexing: bool,
                            md5_dict: Dict[str, str]):
    """
    Process UCSC references by downloading and indexing.

    Args:
        ucsc_refs (dict): Dictionary of UCSC references.
        output_dir (Path): Base output directory.
        bwa_path (str): Path to the bwa executable.
        skip_indexing (bool): Whether to skip the indexing step.
        md5_dict (dict): Dictionary to store MD5 checksums.
    """
    for ref_name, ref_info in ucsc_refs.items():
        url = ref_info.get("url")
        target_path = output_dir / ref_info.get("target_path")
        index_command = ref_info.get("index_command", None)

        if not url or not ref_info.get("target_path"):
            logging.warning(
                f"Missing URL or target_path for UCSC reference {ref_name}. Skipping."
            )
            continue

        # Download the UCSC reference
        download_file(url, target_path)

        # Calculate MD5 checksum
        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logging.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        # Index the reference if required and not skipped
        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logging.info(f"Skipping indexing for {target_path}")


def process_vntyper_references(vntyper_refs: Dict[str, Dict[str, str]],
                               output_dir: Path, bwa_path: str, skip_indexing: bool,
                               md5_dict: Dict[str, str]):
    """
    Process VNtyper references by downloading and extracting.

    Args:
        vntyper_refs (dict): Dictionary of VNtyper references.
        output_dir (Path): Base output directory.
        bwa_path (str): Path to the bwa executable.
        skip_indexing (bool): Whether to skip the indexing step.
        md5_dict (dict): Dictionary to store MD5 checksums.
    """
    for ref_name, ref_info in vntyper_refs.items():
        url = ref_info.get("url")
        target_path = output_dir / ref_info.get("target_path")
        extract_to = ref_info.get("extract_to", None)
        index_command = ref_info.get("index_command", None)

        if not url or not ref_info.get("target_path"):
            logging.warning(
                f"Missing URL or target_path for VNtyper reference {ref_name}. Skipping."
            )
            continue

        # Download the VNtyper reference
        download_file(url, target_path)

        # Calculate MD5 checksum
        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logging.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        # Extract if required
        if extract_to:
            extract_dir = output_dir / extract_to
            extract_dir.mkdir(parents=True, exist_ok=True)
            if target_path.suffix == '.zip':
                try:
                    with zipfile.ZipFile(target_path, 'r') as zip_ref:
                        zip_ref.extractall(path=extract_dir)
                    logging.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    logging.error(f"Failed to extract {target_path}: {e}")
                    sys.exit(1)
            elif target_path.suffixes[-2:] == ['.tar', '.gz'] or target_path.suffix == '.tgz':
                try:
                    with tarfile.open(target_path, "r:gz") as tar:
                        tar.extractall(path=extract_dir)
                    logging.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    logging.error(f"Failed to extract {target_path}: {e}")
                    sys.exit(1)
            else:
                logging.warning(
                    f"Unsupported archive format for {target_path}. Skipping extraction."
                )

        # Index the reference if required and not skipped
        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logging.info(f"Skipping indexing for {target_path}")


def process_own_repository_references(own_repo_refs: Dict[str, Any],
                                      output_dir: Path, skip_indexing: bool,
                                      md5_dict: Dict[str, str]):
    """
    Process own repository references by downloading specific FASTA files.

    Args:
        own_repo_refs (dict): Dictionary of own repository references.
        output_dir (Path): Base output directory.
        skip_indexing (bool): Whether to skip the indexing step.
        md5_dict (dict): Dictionary to store MD5 checksums.
    """
    raw_files: List[Dict[str, str]] = own_repo_refs.get("raw_files", [])
    for file_info in raw_files:
        url = file_info.get("url")
        target_path = output_dir / file_info.get("target_path")
        index_command = file_info.get("index_command", None)

        if not url or not file_info.get("target_path"):
            logging.warning(
                "Missing URL or target_path for own repository raw file. Skipping."
            )
            continue

        # Download the raw file
        download_file(url, target_path)

        # Calculate MD5 checksum
        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logging.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        # Index the file if required and not skipped
        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logging.info(f"Skipping indexing for {target_path}")


def write_md5_checksums(md5_dict: Dict[str, str], output_dir: Path):
    """
    Write the MD5 checksums to a file in the output directory.

    Args:
        md5_dict (dict): Dictionary mapping file paths to their MD5 checksums.
        output_dir (Path): Base output directory.
    """
    checksum_file = output_dir / "md5_checksums.txt"
    try:
        with checksum_file.open('w') as f:
            for file_path, md5 in md5_dict.items():
                relative_path = file_path.replace(str(output_dir) + "/", "")
                f.write(f"{md5}  {relative_path}\n")
        logging.info(f"MD5 checksums written to {checksum_file}")
    except Exception as e:
        logging.error(f"Failed to write MD5 checksums to {checksum_file}: {e}")
        sys.exit(1)


def setup_logging(output_dir: Path):
    """
    Setup logging to output to both stdout and a log file in the output directory.

    Args:
        output_dir (Path): Directory where logs will be stored.
    """
    log_file = output_dir / "install_references.log"

    # Create a custom logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Remove all handlers associated with the root logger object
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Create handlers
    c_handler = logging.StreamHandler(sys.stdout)
    f_handler = logging.FileHandler(log_file)

    c_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.INFO)

    # Create formatters and add them to handlers
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    c_handler.setFormatter(formatter)
    f_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

    logging.info(f"Logging initialized. Logs will be saved to {log_file}")


def main(output_dir: Path, config_path: Optional[Path] = None, skip_indexing: bool = False):
    """
    Main function to execute the install_references process.

    Args:
        output_dir (Path): Directory where references will be installed.
        config_path (Optional[Path]): Path to the main config.json file to update.
        skip_indexing (bool): Whether to skip the indexing step.
    """
    # Define the path to the install_references_config.json
    script_dir = Path(__file__).parent
    install_config_path = script_dir / "install_references_config.json"

    # Load the installation configuration
    install_config = load_install_config(install_config_path)

    ucsc_refs = install_config.get("ucsc_references", {})
    vntyper_refs = install_config.get("vntyper_references", {})
    own_repo_refs = install_config.get("own_repository_references", {})
    bwa_path = install_config.get("bwa_path", "bwa")  # Default to 'bwa'

    # Set the base output directory
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logging.error(f"Failed to create output directory {output_dir}: {e}")
        sys.exit(1)

    # Setup logging after output_dir is created
    setup_logging(output_dir)

    # Initialize a dictionary to store MD5 checksums
    md5_dict = {}

    # Process UCSC references
    if ucsc_refs:
        logging.info("Processing UCSC references...")
        process_ucsc_references(ucsc_refs, output_dir, bwa_path, skip_indexing, md5_dict)

    # Process VNtyper references
    if vntyper_refs:
        logging.info("Processing VNtyper references...")
        process_vntyper_references(vntyper_refs, output_dir, bwa_path, skip_indexing, md5_dict)

    # Process own repository references
    if own_repo_refs:
        logging.info("Processing own repository references...")
        process_own_repository_references(own_repo_refs, output_dir, skip_indexing, md5_dict)

    # Write MD5 checksums to file
    if md5_dict:
        write_md5_checksums(md5_dict, output_dir)

    # Update the main config.json with new reference paths if config_path is provided
    if config_path and config_path.exists():
        # Prepare a dictionary with updated reference paths
        updated_references = {}

        # Collect all references from UCSC
        for ref_key, ref_info in ucsc_refs.items():
            ref_path = output_dir / ref_info.get("target_path")
            updated_references[f"ucsc_{ref_key}"] = ref_path.resolve()

        # Collect all references from VNtyper
        for ref_key, ref_info in vntyper_refs.items():
            ref_path = output_dir / ref_info.get("target_path")
            updated_references[f"vntyper_{ref_key}"] = ref_path.resolve()

        # Collect references from own repository raw_files
        raw_files: List[Dict[str, str]] = own_repo_refs.get("raw_files", [])
        for file_info in raw_files:
            ref_name = Path(file_info.get("target_path")).stem
            ref_path = output_dir / file_info.get("target_path")
            updated_references[f"own_repo_{ref_name}"] = ref_path.resolve()

        # Update config.json
        update_config(config_path, updated_references)
    else:
        if config_path:
            logging.warning(
                f"Config file {config_path} does not exist. Skipping config update."
            )
        else:
            logging.info("No config_path provided. Skipping config update.")

    logging.info("All references have been installed and configured successfully.")


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments for standalone execution
    parser = argparse.ArgumentParser(
        description="Install necessary reference files for vntyper."
    )
    parser.add_argument(
        '-d', '--output-dir',
        type=Path,
        required=True,
        help="Directory where references will be installed."
    )
    parser.add_argument(
        '-c', '--config-path',
        type=Path,
        default=None,
        help="Path to the main config.json file to update. If not provided, config update is skipped."
    )
    parser.add_argument(
        '--skip-indexing',
        action='store_true',
        help="Skip the indexing step during reference installation."
    )
    args = parser.parse_args()
    main(args.output_dir, args.config_path, args.skip_indexing)
