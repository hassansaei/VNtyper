#!/usr/bin/env python3
# vntyper/scripts/install_references.py

import gzip
import hashlib
import json
import logging
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path
from typing import Any, Optional
from urllib.request import urlretrieve


def load_install_config(config_path: Path) -> dict[str, Any]:
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
        with config_path.open("r") as f:
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
        args = command.split()
        subprocess.run(args, check=True, capture_output=True)
        logging.info(f"Successfully executed: {command}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Indexing command failed for {fasta_path}: {e.stderr.decode().strip()}")
        sys.exit(1)


###############################################################################
# Multi-Aligner Support Functions
###############################################################################


def check_executable_available(executable: str) -> bool:
    """
    Check if an executable is available in the system PATH.

    Args:
        executable (str): Name or path of the executable to check.

    Returns:
        bool: True if executable is available, False otherwise.
    """
    try:
        result = subprocess.run(["which", executable], capture_output=True, text=True, check=False)
        if result.returncode == 0:
            logging.debug(f"Found executable: {executable} at {result.stdout.strip()}")
            return True
        else:
            logging.debug(f"Executable not found: {executable}")
            return False
    except Exception as e:
        logging.debug(f"Error checking executable {executable}: {e}")
        return False


def get_enabled_aligners(aligner_config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """
    Get dictionary of enabled aligners from configuration.

    Args:
        aligner_config (dict): Dictionary of aligner configurations.

    Returns:
        dict: Dictionary of enabled aligners with their configurations.
    """
    enabled = {}
    for aligner_name, aligner_info in aligner_config.items():
        if aligner_info.get("enabled", False):
            executable = aligner_info.get("executable", aligner_name)
            if check_executable_available(executable):
                enabled[aligner_name] = aligner_info
                logging.info(f"  ✓ {aligner_name}: {aligner_info.get('description', 'No description')}")
            else:
                logging.warning(
                    f"  ✗ {aligner_name} is enabled in config but executable '{executable}' not found. "
                    f"Skipping this aligner."
                )
    return enabled


def detect_index_conflicts(aligners: dict[str, dict[str, Any]]) -> list[str]:
    """
    Detect potential conflicts between aligner index file extensions.

    Args:
        aligners (dict): Dictionary of aligner configurations.

    Returns:
        list: List of warning messages about conflicts (empty if none).
    """
    warnings = []
    extension_map = {}

    for aligner_name, aligner_info in aligners.items():
        index_files = aligner_info.get("index_files", [])
        for ext in index_files:
            if ext in extension_map:
                warnings.append(
                    f"Index file extension '{ext}' is used by both "
                    f"'{extension_map[ext]}' and '{aligner_name}'. "
                    f"This may cause conflicts."
                )
            else:
                extension_map[ext] = aligner_name

    return warnings


def check_index_exists(ref_path: Path, aligner_name: str, aligner_info: dict[str, Any]) -> bool:
    """
    Check if all required index files exist for an aligner.

    Args:
        ref_path (Path): Path to the reference FASTA file.
        aligner_name (str): Name of the aligner.
        aligner_info (dict): Aligner configuration dictionary.

    Returns:
        bool: True if all index files exist, False otherwise.
    """
    index_files = aligner_info.get("index_files", [])

    if aligner_info.get("index_dir_required", False):
        # DRAGMAP-style: check for directory with index files
        index_dir = ref_path.parent / f"{ref_path.stem}_{aligner_name}_index"
        if not index_dir.exists():
            return False
        return all((index_dir / index_file).exists() for index_file in index_files)
    elif aligner_info.get("requires_index_base", False):
        # Bowtie2-style: check for index base + extensions
        index_base = ref_path.parent / f"{ref_path.stem}_{aligner_name}"
        return all(Path(str(index_base) + ext).exists() for ext in index_files)
    else:
        # Standard: check for ref_path + extension
        for ext in index_files:
            index_file_path = Path(str(ref_path) + ext)
            if not index_file_path.exists():
                logging.debug(f"Missing index file: {index_file_path}")
                return False
        return True


def execute_aligner_index(ref_path: Path, aligner_name: str, aligner_info: dict[str, Any], threads: int = 4) -> bool:
    """
    Execute indexing for a specific aligner.

    Args:
        ref_path (Path): Path to the reference FASTA file.
        aligner_name (str): Name of the aligner.
        aligner_info (dict): Aligner configuration dictionary.
        threads (int): Number of threads to use (if supported).

    Returns:
        bool: True if indexing succeeded, False otherwise.
    """
    index_command_template = aligner_info.get("index_command", "")
    if not index_command_template:
        logging.error(f"No index_command specified for {aligner_name}")
        return False

    # Prepare command parameters
    params = {
        "ref_path": str(ref_path),
        "threads": threads if aligner_info.get("supports_threading", False) else aligner_info.get("threads_default", 4),
    }

    # Handle special cases
    if aligner_info.get("index_dir_required", False):
        # DRAGMAP: needs separate index directory
        index_dir = ref_path.parent / f"{ref_path.stem}_{aligner_name}_index"
        index_dir.mkdir(parents=True, exist_ok=True)
        params["index_dir"] = str(index_dir)

    if aligner_info.get("requires_index_base", False):
        # Bowtie2: needs separate index base name
        index_base = ref_path.parent / f"{ref_path.stem}_{aligner_name}"
        params["index_base"] = str(index_base)

    # Format command
    try:
        command = index_command_template.format(**params)
    except KeyError as e:
        logging.error(f"Missing parameter in index command for {aligner_name}: {e}")
        return False

    logging.info(f"  Indexing with {aligner_name}...")
    logging.debug(f"  Command: {command}")

    try:
        subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logging.info(f"  ✓ {aligner_name} indexing complete")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"  ✗ {aligner_name} indexing failed: {e.stderr.strip()}")
        return False


def index_reference_with_aligners(
    ref_path: Path, aligners: dict[str, dict[str, Any]], threads: int = 4, force_reindex: bool = False
) -> dict[str, bool]:
    """
    Index a reference file with multiple aligners.

    Args:
        ref_path (Path): Path to the reference FASTA file.
        aligners (dict): Dictionary of aligner configurations.
        threads (int): Number of threads to use for indexing.
        force_reindex (bool): Force re-indexing even if indices exist.

    Returns:
        dict: Dictionary mapping aligner names to success status (True/False).
    """
    results = {}

    logging.info(f"Indexing reference: {ref_path.name}")
    logging.info(f"  Using {len(aligners)} aligner(s) with {threads} threads")

    for aligner_name, aligner_info in aligners.items():
        # Check if index already exists
        if not force_reindex and check_index_exists(ref_path, aligner_name, aligner_info):
            logging.info(f"  ✓ {aligner_name} index already exists, skipping")
            results[aligner_name] = True
            continue

        # Execute indexing
        success = execute_aligner_index(ref_path, aligner_name, aligner_info, threads)
        results[aligner_name] = success

    return results


def update_config(config_path: Path, references: dict[str, Path]):
    """
    Update the main config.json with paths to the downloaded references.

    Args:
        config_path (Path): Path to the main config.json file.
        references (dict): Dictionary mapping reference keys to their paths.

    Raises:
        SystemExit: If updating the config fails.
    """
    if not config_path.exists():
        logging.error(f"Main config file {config_path} does not exist. Cannot update references.")
        sys.exit(1)

    try:
        with config_path.open("r") as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        logging.error(f"Error parsing main config.json: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unexpected error reading main config.json: {e}")
        sys.exit(1)

    if "reference_data" not in config:
        config["reference_data"] = {}

    for ref_key, ref_path in references.items():
        config["reference_data"][ref_key] = str(ref_path)

    try:
        with config_path.open("w") as f:
            json.dump(config, f, indent=2)
        logging.info(f"Successfully updated {config_path} with new reference paths.")
    except Exception as e:
        logging.error(f"Failed to write updated config.json: {e}")
        sys.exit(1)


def process_ucsc_references(
    ucsc_refs: dict[str, dict[str, str]],
    output_dir: Path,
    bwa_path: str,
    skip_indexing: bool,
    md5_dict: dict[str, str],
    aligners: Optional[dict[str, dict[str, Any]]] = None,
    index_threads: int = 4,
):
    """
    Process UCSC references by downloading and indexing.

    Args:
        ucsc_refs (dict): Dictionary of UCSC references.
        output_dir (Path): Base output directory.
        bwa_path (str): Path to the bwa executable (legacy, kept for compatibility).
        skip_indexing (bool): Whether to skip the indexing step.
        md5_dict (dict): Dictionary to store MD5 checksums.
        aligners (dict, optional): Dictionary of aligner configurations for multi-aligner indexing.
        index_threads (int): Number of threads to use for indexing.
    """
    for ref_name, ref_info in ucsc_refs.items():
        url = ref_info.get("url")
        target_path = output_dir / ref_info.get("target_path")
        index_command = ref_info.get("index_command", None)

        if not url or not ref_info.get("target_path"):
            logging.warning(f"Missing URL or target_path for UCSC reference {ref_name}. Skipping.")
            continue

        download_file(url, target_path)

        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logging.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        if target_path.suffix == ".zip":
            try:
                with zipfile.ZipFile(target_path, "r") as zip_ref:
                    zip_ref.extractall(path=target_path)
                logging.info(f"Successfully extracted {target_path.name}")
            except Exception as e:
                logging.error(f"Failed to extract {target_path}: {e}")
                sys.exit(1)
        elif target_path.suffix == ".gz":
            try:
                output_path = target_path.with_suffix("")
                with gzip.open(target_path, "rb") as f_in, open(output_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                logging.info(f"Successfully extracted {target_path.name} to {output_path.name}")
            except Exception as e:
                logging.error(f"Failed to extract {target_path}: {e}")
                sys.exit(1)
        elif target_path.suffixes[-2:] == [".tar", ".gz"] or target_path.suffix == ".tgz":
            try:
                with tarfile.open(target_path, "r:gz") as tar:
                    tar.extractall(path=target_path)
                logging.info(f"Successfully extracted {target_path.name}")
            except Exception as e:
                logging.error(f"Failed to extract {target_path}: {e}")
                sys.exit(1)
        else:
            logging.warning(f"Unsupported archive format for {target_path}. Skipping extraction.")

        # Multi-aligner indexing
        if not skip_indexing:
            output_path = target_path.with_suffix("")

            if aligners:
                # Use multi-aligner indexing
                index_reference_with_aligners(output_path, aligners, threads=index_threads, force_reindex=False)
            elif index_command:
                # Fall back to legacy single indexing command
                logging.warning(f"No aligners configured, using legacy index_command for {output_path.name}")
                execute_index_command(index_command, output_path)
        elif skip_indexing:
            logging.info(f"Skipping indexing for {target_path.with_suffix('')}")


def process_vntyper_references(
    vntyper_refs: dict[str, dict[str, str]],
    output_dir: Path,
    bwa_path: str,
    skip_indexing: bool,
    md5_dict: dict[str, str],
):
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
            logging.warning(f"Missing URL or target_path for VNtyper reference {ref_name}. Skipping.")
            continue

        download_file(url, target_path)

        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logging.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        if extract_to:
            extract_dir = output_dir / extract_to
            extract_dir.mkdir(parents=True, exist_ok=True)
            if target_path.suffix == ".zip":
                try:
                    with zipfile.ZipFile(target_path, "r") as zip_ref:
                        zip_ref.extractall(path=extract_dir)
                    logging.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    logging.error(f"Failed to extract {target_path}: {e}")
                    sys.exit(1)
            elif target_path.suffixes[-2:] == [".tar", ".gz"] or target_path.suffix == ".tgz":
                try:
                    with tarfile.open(target_path, "r:gz") as tar:
                        tar.extractall(path=extract_dir)
                    logging.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    logging.error(f"Failed to extract {target_path}: {e}")
                    sys.exit(1)
            else:
                logging.warning(f"Unsupported archive format for {target_path}. Skipping extraction.")

        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logging.info(f"Skipping indexing for {target_path}")


def process_own_repository_references(
    own_repo_refs: dict[str, Any],
    output_dir: Path,
    skip_indexing: bool,
    md5_dict: dict[str, str],
):
    """
    Process own repository references by downloading specific FASTA files.

    Args:
        own_repo_refs (dict): Dictionary of own repository references.
        output_dir (Path): Base output directory.
        skip_indexing (bool): Whether to skip the indexing step.
        md5_dict (dict): Dictionary to store MD5 checksums.
    """
    raw_files: list[dict[str, str]] = own_repo_refs.get("raw_files", [])
    for file_info in raw_files:
        url = file_info.get("url")
        target_path = output_dir / file_info.get("target_path")
        index_command = file_info.get("index_command", None)

        if not url or not file_info.get("target_path"):
            logging.warning("Missing URL or target_path for own repository raw file. Skipping.")
            continue

        download_file(url, target_path)

        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logging.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logging.info(f"Skipping indexing for {target_path}")


def write_md5_checksums(md5_dict: dict[str, str], output_dir: Path):
    """
    Write the MD5 checksums to a file in the output directory.

    Args:
        md5_dict (dict): Dictionary mapping file paths to their MD5 checksums.
        output_dir (Path): Base output directory.
    """
    checksum_file = output_dir / "md5_checksums.txt"
    try:
        with checksum_file.open("w") as f:
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

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    c_handler = logging.StreamHandler(sys.stdout)
    f_handler = logging.FileHandler(log_file)

    c_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    c_handler.setFormatter(formatter)
    f_handler.setFormatter(formatter)

    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

    logging.info(f"Logging initialized. Logs will be saved to {log_file}")


def main(
    output_dir: Path,
    config_path: Optional[Path] = None,
    skip_indexing: bool = False,
    index_threads: int = 4,
    aligners_to_use: Optional[list[str]] = None,
    references_to_process: Optional[list[str]] = None,
):
    """
    Main function to execute the install_references process.

    Args:
        output_dir (Path): Directory where references will be installed.
        config_path (Optional[Path]): Path to the main config.json file to update.
        skip_indexing (bool): Whether to skip the indexing step.
        index_threads (int): Number of threads to use for indexing.
        aligners_to_use (list, optional): List of specific aligners to use (overrides config).
        references_to_process (list, optional): List of specific references to process (e.g., ['hg19', 'hg38']).
    """
    script_dir = Path(__file__).parent
    install_config_path = script_dir / "install_references_config.json"

    install_config = load_install_config(install_config_path)

    ucsc_refs = install_config.get("ucsc_references", {})
    ncbi_refs = install_config.get("ncbi_references", {})
    vntyper_refs = install_config.get("vntyper_references", {})
    own_repo_refs = install_config.get("own_repository_references", {})
    bwa_path = install_config.get("bwa_path", "bwa")  # Default to 'bwa'

    # Load aligner configurations
    aligner_config = install_config.get("aligners", {})

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logging.error(f"Failed to create output directory {output_dir}: {e}")
        sys.exit(1)

    setup_logging(output_dir)

    # Filter references (after logging is set up)
    # Default to hg19 and hg38 (UCSC) for backward compatibility
    if references_to_process is None:
        references_to_process = ["hg19", "hg38"]
        logging.info("No references specified, using default: hg19, hg38")

    all_available_refs = set(ucsc_refs.keys()) | set(ncbi_refs.keys()) | set(vntyper_refs.keys())
    requested_refs = set(references_to_process)
    found_refs = requested_refs & all_available_refs
    missing_refs = requested_refs - all_available_refs

    if missing_refs:
        logging.warning(f"Requested references not found in config: {', '.join(sorted(missing_refs))}")
        logging.warning(f"Available references: {', '.join(sorted(all_available_refs))}")

    if not found_refs:
        logging.error("None of the requested references were found in the configuration.")
        sys.exit(1)

    logging.info(f"Processing references: {', '.join(sorted(found_refs))}")

    ucsc_refs = {k: v for k, v in ucsc_refs.items() if k in references_to_process}
    ncbi_refs = {k: v for k, v in ncbi_refs.items() if k in references_to_process}
    vntyper_refs = {k: v for k, v in vntyper_refs.items() if k in references_to_process}

    # Initialize aligners
    enabled_aligners = {}
    if not skip_indexing and aligner_config:
        logging.info("=" * 80)
        logging.info("ALIGNER CONFIGURATION")
        logging.info("=" * 80)
        logging.info("Checking available aligners:")

        # Get enabled aligners
        all_enabled = get_enabled_aligners(aligner_config)

        # Filter by user-specified aligners if provided, otherwise default to BWA only
        if aligners_to_use:
            for aligner_name in aligners_to_use:
                if aligner_name in all_enabled:
                    enabled_aligners[aligner_name] = all_enabled[aligner_name]
                elif aligner_name in aligner_config:
                    logging.warning(f"  ✗ {aligner_name} was specified but is not available or not enabled")
                else:
                    logging.error(f"  ✗ Unknown aligner: {aligner_name}")
        else:
            # Default to BWA only
            if "bwa" in all_enabled:
                enabled_aligners["bwa"] = all_enabled["bwa"]
                logging.info("  Using default aligner: bwa")
            else:
                logging.warning("  Default aligner 'bwa' not available")

        if not enabled_aligners:
            logging.warning("No aligners available. Indexing will be skipped.")
            logging.warning(
                "To enable aligners, install them and ensure they are in your PATH, "
                "or set enabled:true in install_references_config.json"
            )
        else:
            logging.info(f"\nWill use {len(enabled_aligners)} aligner(s) for indexing:")
            for name in enabled_aligners:
                logging.info(f"  • {name}")

            # Detect index file conflicts
            conflicts = detect_index_conflicts(enabled_aligners)
            if conflicts:
                logging.warning("\n⚠ Index file conflicts detected:")
                for warning in conflicts:
                    logging.warning(f"  • {warning}")
                logging.warning("  These conflicts may cause issues if aligners overwrite each other's index files.")
            else:
                logging.info("  ✓ No index file conflicts detected")

        logging.info("=" * 80)
        logging.info("")

    md5_dict = {}

    # Process UCSC references
    if ucsc_refs:
        logging.info("Processing UCSC references...")
        process_ucsc_references(
            ucsc_refs,
            output_dir,
            bwa_path,
            skip_indexing,
            md5_dict,
            aligners=enabled_aligners,
            index_threads=index_threads,
        )

    # Process NCBI references
    if ncbi_refs:
        logging.info("Processing NCBI references...")
        process_ucsc_references(
            ncbi_refs,
            output_dir,
            bwa_path,
            skip_indexing,
            md5_dict,
            aligners=enabled_aligners,
            index_threads=index_threads,
        )

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
        updated_references = {}

        # Collect all references from UCSC
        for ref_key, ref_info in ucsc_refs.items():
            ref_path = output_dir / ref_info.get("target_path")
            updated_references[f"ucsc_{ref_key}"] = ref_path.resolve()

        # Collect all references from NCBI
        for ref_key, ref_info in ncbi_refs.items():
            ref_path = output_dir / ref_info.get("target_path")
            updated_references[f"ncbi_{ref_key}"] = ref_path.resolve()

        # Collect all references from VNtyper
        for ref_key, ref_info in vntyper_refs.items():
            ref_path = output_dir / ref_info.get("target_path")
            updated_references[f"vntyper_{ref_key}"] = ref_path.resolve()

        # Collect references from own repository
        raw_files: list[dict[str, str]] = own_repo_refs.get("raw_files", [])
        for file_info in raw_files:
            ref_name = Path(file_info.get("target_path")).stem
            ref_path = output_dir / file_info.get("target_path")
            updated_references[f"own_repo_{ref_name}"] = ref_path.resolve()

        update_config(config_path, updated_references)
    else:
        if config_path:
            logging.warning(f"Config file {config_path} does not exist. Skipping config update.")
        else:
            logging.info("No config_path provided. Skipping config update.")

    logging.info("All references have been installed and configured successfully.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Install necessary reference files for vntyper with multi-aligner support.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Install all references with all enabled aligners
  python install_references.py -d reference/

  # Install specific reference only (e.g., hg19)
  python install_references.py -d reference/ --references hg19

  # Install with specific aligners only
  python install_references.py -d reference/ --aligners bwa minimap2

  # Install hg19 with BWA aligner and 8 threads
  python install_references.py -d reference/ --references hg19 --aligners bwa --threads 8

  # Install with 8 threads
  python install_references.py -d reference/ --threads 8

  # Skip indexing
  python install_references.py -d reference/ --skip-indexing
        """,
    )
    parser.add_argument(
        "-d",
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where references will be installed.",
    )
    parser.add_argument(
        "-c",
        "--config-path",
        type=Path,
        default=None,
        help="Path to the main config.json file to update. If not provided, config update is skipped.",
    )
    parser.add_argument(
        "--skip-indexing",
        action="store_true",
        help="Skip the indexing step during reference installation.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use for indexing (default: 4).",
    )
    parser.add_argument(
        "--aligners",
        nargs="+",
        default=None,
        metavar="ALIGNER",
        help="Specific aligners to use (e.g., bwa bwa-mem2 minimap2). "
        "If not specified, only BWA will be used (default).",
    )
    parser.add_argument(
        "--references",
        nargs="+",
        default=None,
        metavar="REFERENCE",
        help="Specific references to process (e.g., hg19 hg38 GRCh37 GRCh38). "
        "Default: hg19 hg38 (UCSC references only).",
    )

    args = parser.parse_args()
    main(args.output_dir, args.config_path, args.skip_indexing, args.threads, args.aligners, args.references)
