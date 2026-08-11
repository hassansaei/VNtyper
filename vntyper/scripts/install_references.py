# vntyper/scripts/install_references.py

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import logging
import shlex
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path
from typing import Any
from urllib.request import urlretrieve

# `docker/Dockerfile.base` copies this module and `reference_bundle.py` alone into a
# build stage and runs them without installing the package, and every module imported
# here joins the base-image content hash. `reference_bundle` is therefore the only
# `vntyper` import allowed at module scope; anything else is inlined below or imported
# inside the function that needs it.
from vntyper.scripts.reference_bundle import verify_sha256

logger = logging.getLogger(__name__)

#: Config sections holding one physical chromosome FASTA per entry, in the order
#: `--from-source` walks them.
GENOME_SECTIONS = ("ucsc_references", "ncbi_references", "ensembl_references")


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
        logger.error(f"Installation config file not found at {config_path}")
        sys.exit(1)

    try:
        with config_path.open("r") as f:
            config = json.load(f)
        return config
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing JSON config: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error reading config: {e}")
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
        logger.info(f"File already exists at {dest_path}. Skipping download.")
        return

    logger.info(f"Downloading from {url} to {dest_path}...")
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(url, dest_path)
        logger.info(f"Successfully downloaded {dest_path.name}")
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
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
    logger.debug(f"Calculating MD5 for {file_path}...")
    hash_md5 = hashlib.md5()
    try:
        with file_path.open("rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        md5_checksum = hash_md5.hexdigest()
        logger.debug(f"MD5 for {file_path} is {md5_checksum}")
        return md5_checksum
    except Exception as e:
        logger.error(f"Failed to calculate MD5 for {file_path}: {e}")
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
    logger.info(f"Executing indexing command: {command}")
    try:
        args = command.split()
        subprocess.run(args, check=True, capture_output=True)
        logger.info(f"Successfully executed: {command}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Indexing command failed for {fasta_path}: {e.stderr.decode().strip()}")
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
            logger.debug(f"Found executable: {executable} at {result.stdout.strip()}")
            return True
        else:
            logger.debug(f"Executable not found: {executable}")
            return False
    except Exception as e:
        logger.debug(f"Error checking executable {executable}: {e}")
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
                logger.info(f"  ✓ {aligner_name}: {aligner_info.get('description', 'No description')}")
            else:
                logger.warning(
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
    extension_map: dict[str, str] = {}

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
                logger.debug(f"Missing index file: {index_file_path}")
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
        logger.error(f"No index_command specified for {aligner_name}")
        return False

    # Prepare command parameters. Paths are quoted with `_quote` before they are
    # interpolated into the shell=True command below; `threads` is numeric and is never
    # a candidate for quoting.
    params = {
        "ref_path": _quote(ref_path),
        "threads": threads if aligner_info.get("supports_threading", False) else aligner_info.get("threads_default", 4),
    }

    # Handle special cases
    if aligner_info.get("index_dir_required", False):
        # DRAGMAP: needs separate index directory
        index_dir = ref_path.parent / f"{ref_path.stem}_{aligner_name}_index"
        index_dir.mkdir(parents=True, exist_ok=True)
        params["index_dir"] = _quote(index_dir)

    if aligner_info.get("requires_index_base", False):
        # Bowtie2: needs separate index base name
        index_base = ref_path.parent / f"{ref_path.stem}_{aligner_name}"
        params["index_base"] = _quote(index_base)

    # Format command
    try:
        command = index_command_template.format(**params)
    except KeyError as e:
        logger.error(f"Missing parameter in index command for {aligner_name}: {e}")
        return False

    logger.info(f"  Indexing with {aligner_name}...")
    logger.debug(f"  Command: {command}")

    try:
        subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logger.info(f"  ✓ {aligner_name} indexing complete")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ {aligner_name} indexing failed: {e.stderr.strip()}")
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

    logger.info(f"Indexing reference: {ref_path.name}")
    logger.info(f"  Using {len(aligners)} aligner(s) with {threads} threads")

    for aligner_name, aligner_info in aligners.items():
        # Check if index already exists
        if not force_reindex and check_index_exists(ref_path, aligner_name, aligner_info):
            logger.info(f"  ✓ {aligner_name} index already exists, skipping")
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
        logger.error(f"Main config file {config_path} does not exist. Cannot update references.")
        sys.exit(1)

    try:
        with config_path.open("r") as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing main config.json: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error reading main config.json: {e}")
        sys.exit(1)

    if "reference_data" not in config:
        config["reference_data"] = {}

    for ref_key, ref_path in references.items():
        config["reference_data"][ref_key] = str(ref_path)

    try:
        with config_path.open("w") as f:
            json.dump(config, f, indent=2)
        logger.info(f"Successfully updated {config_path} with new reference paths.")
    except Exception as e:
        logger.error(f"Failed to write updated config.json: {e}")
        sys.exit(1)


def process_ucsc_references(
    ucsc_refs: dict[str, dict[str, str]],
    output_dir: Path,
    bwa_path: str,
    skip_indexing: bool,
    md5_dict: dict[str, str],
    aligners: dict[str, dict[str, Any]] | None = None,
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
        target_path_str = ref_info.get("target_path")
        index_command = ref_info.get("index_command", None)

        if not url or not target_path_str:
            logger.warning(f"Missing URL or target_path for UCSC reference {ref_name}. Skipping.")
            continue

        target_path = output_dir / target_path_str

        download_file(url, target_path)

        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logger.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        if target_path.suffix == ".zip":
            try:
                with zipfile.ZipFile(target_path, "r") as zip_ref:
                    zip_ref.extractall(path=target_path)
                logger.info(f"Successfully extracted {target_path.name}")
            except Exception as e:
                logger.error(f"Failed to extract {target_path}: {e}")
                sys.exit(1)
        elif target_path.suffix == ".gz":
            try:
                output_path = target_path.with_suffix("")
                with gzip.open(target_path, "rb") as f_in, open(output_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                logger.info(f"Successfully extracted {target_path.name} to {output_path.name}")
            except Exception as e:
                logger.error(f"Failed to extract {target_path}: {e}")
                sys.exit(1)
        elif target_path.suffixes[-2:] == [".tar", ".gz"] or target_path.suffix == ".tgz":
            try:
                with tarfile.open(target_path, "r:gz") as tar:
                    tar.extractall(path=target_path)
                logger.info(f"Successfully extracted {target_path.name}")
            except Exception as e:
                logger.error(f"Failed to extract {target_path}: {e}")
                sys.exit(1)
        else:
            logger.warning(f"Unsupported archive format for {target_path}. Skipping extraction.")

        # Multi-aligner indexing
        if not skip_indexing:
            output_path = target_path.with_suffix("")

            if aligners:
                # Use multi-aligner indexing
                index_reference_with_aligners(output_path, aligners, threads=index_threads, force_reindex=False)
            elif index_command:
                # Fall back to legacy single indexing command
                logger.warning(f"No aligners configured, using legacy index_command for {output_path.name}")
                execute_index_command(index_command, output_path)
        elif skip_indexing:
            logger.info(f"Skipping indexing for {target_path.with_suffix('')}")


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
        target_path_str = ref_info.get("target_path")
        extract_to = ref_info.get("extract_to", None)
        index_command = ref_info.get("index_command", None)

        if not url or not target_path_str:
            logger.warning(f"Missing URL or target_path for VNtyper reference {ref_name}. Skipping.")
            continue

        target_path = output_dir / target_path_str

        download_file(url, target_path)

        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logger.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        if extract_to:
            extract_dir = output_dir / extract_to
            extract_dir.mkdir(parents=True, exist_ok=True)
            if target_path.suffix == ".zip":
                try:
                    with zipfile.ZipFile(target_path, "r") as zip_ref:
                        zip_ref.extractall(path=extract_dir)
                    logger.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    logger.error(f"Failed to extract {target_path}: {e}")
                    sys.exit(1)
            elif target_path.suffixes[-2:] == [".tar", ".gz"] or target_path.suffix == ".tgz":
                try:
                    with tarfile.open(target_path, "r:gz") as tar:
                        tar.extractall(path=extract_dir)
                    logger.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    logger.error(f"Failed to extract {target_path}: {e}")
                    sys.exit(1)
            else:
                logger.warning(f"Unsupported archive format for {target_path}. Skipping extraction.")

        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logger.info(f"Skipping indexing for {target_path}")


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
        target_path_str = file_info.get("target_path")
        index_command = file_info.get("index_command", None)

        if not url or not target_path_str:
            logger.warning("Missing URL or target_path for own repository raw file. Skipping.")
            continue

        target_path = output_dir / target_path_str

        download_file(url, target_path)

        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logger.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logger.info(f"Skipping indexing for {target_path}")


###############################################################################
# Building references from their upstream sources
###############################################################################


def _quote(value: object) -> str:
    """Quote one operand for a shell command.

    A three-line copy of :func:`vntyper.scripts.samtools_command_fragments.quote_path`
    (which `command_builders` only re-exports), deliberately duplicated rather than
    imported: either module would join the Docker base-image content-hash set, where
    every future edit costs a 25-120 minute rebuild. See the import comment at the top
    of this file. ``TestQuote`` in `tests/unit/test_install_references_derivations.py`
    pins that the two helpers behave identically.

    Args:
        value: The operand to quote. Anything with a ``str()``.

    Returns:
        str: The operand as exactly one shell word.
    """
    return shlex.quote(str(value))


def _run_shell(command: str) -> subprocess.CompletedProcess[str]:
    """Run one shell command the way `utils.run_command` would, without importing it.

    ``shell=True`` with an explicit ``/bin/bash`` is the repository-wide convention:
    the commands built here contain a redirect, and every interpolated *path* is
    quoted with :func:`_quote` while command *prefixes* from the config are not, so
    that ``mamba run -n vntyper samtools`` stays several words.

    Args:
        command: The complete shell command.

    Returns:
        subprocess.CompletedProcess[str]: The finished process, stdout/stderr captured.
    """
    return subprocess.run(command, shell=True, executable="/bin/bash", capture_output=True, text=True, check=False)


def _discard_failed_output(destination: Path) -> None:
    """Remove a derivation output that must not survive its own failure.

    Installation merges rather than replaces, so anything left on disk is copied into the
    next run's staging directory and treated as real. A zero-byte FASTA from a failed
    ``samtools faidx`` redirect, or a FASTA that just failed its digest, would be carried
    forward as though it had been produced correctly. Its ``.fai`` goes too, so no index
    is orphaned beside a file that no longer exists.

    Args:
        destination: The derivation output to discard.
    """
    destination.unlink(missing_ok=True)
    Path(f"{destination}.fai").unlink(missing_ok=True)


def derive_region_fasta(source_fasta: Path, region: str, destination: Path, samtools: str) -> Path:
    """Cut a MUC1 region out of a chromosome FASTA.

    This is how both SHARK references are produced. The hg19 one reproduces the FASTA
    this repository used to track byte-for-byte, which is what makes deriving rather
    than shipping them safe.

    Args:
        source_fasta: Indexed chromosome FASTA to cut from.
        region: ``chr1:start-end``, in the source's own contig naming.
        destination: File to write.
        samtools: samtools command prefix from the install config.

    Returns:
        Path: ``destination``.

    Raises:
        RuntimeError: If samtools fails. The truncated output is discarded first - see
            :func:`_discard_failed_output`.
    """
    command = f"{samtools} faidx {_quote(source_fasta)} {_quote(region)} > {_quote(destination)}"
    completed = _run_shell(command)
    if completed.returncode != 0:
        # Bash creates and truncates the redirect target before samtools runs, so a failed
        # derivation has already left a zero-byte FASTA behind.
        _discard_failed_output(destination)
        message = (
            f"faidx derivation failed for {destination.name} "
            f"({region} of {source_fasta.name}): {completed.stderr.strip()}"
        )
        logger.error(message)
        raise RuntimeError(message)
    return destination


def index_fasta_with_samtools(fasta: Path, samtools: str) -> None:
    """Write the ``.fai`` beside a FASTA.

    Both the chromosome FASTAs (so :func:`derive_region_fasta` can cut from them) and
    the derived outputs (so the pipeline can read them) need one.

    Args:
        fasta: FASTA to index in place.
        samtools: samtools command prefix from the install config.

    Raises:
        RuntimeError: If samtools fails.
    """
    completed = _run_shell(f"{samtools} faidx {_quote(fasta)}")
    if completed.returncode != 0:
        message = f"samtools faidx failed for {fasta.name}: {completed.stderr.strip()}"
        logger.error(message)
        raise RuntimeError(message)


def merge_pairwise_motifs(seed_fasta: Path, filter_config: Path, destination: Path) -> Path:
    """Build the pairwise-and-self merged MUC1 motif reference from its two seeds.

    Ported from `reference/generate_vntr_reference.py`, which hardcodes three
    CWD-relative filenames and moves out of this repository with the rest of the seed
    data. The bundle build stages the seeds and then runs ``--from-source``, so this
    has to produce the file with no script beside it to call.

    The logic is reproduced exactly because the output digest is asserted: every
    ordered pair from ``itertools.product`` over an insertion-ordered dict, self-pairs
    included, written as ``>h1-h2`` carrying ``seq(h2) + seq(h1)``. That inversion is
    deliberate in the original and load-bearing for the digest, as is the absence of
    line wrapping.

    Args:
        seed_fasta: Motif seed with one line of sequence per record.
        filter_config: JSON mapping a contig name to the partners disallowed after it.
        destination: File to write.

    Returns:
        Path: ``destination``.

    Raises:
        ValueError: If the seed holds no records, or its headers and sequence lines do
            not pair up one-to-one - which would silently emit a shorter reference
            under the expected name.
    """
    lines = seed_fasta.read_text(encoding="utf-8").splitlines()
    headers = [line.strip() for line in lines if line.startswith(">")]
    sequences = [line.strip() for line in lines if line.strip() and not line.startswith(">")]

    if not headers:
        message = f"{seed_fasta.name} contains no records; cannot build {destination.name}"
        logger.error(message)
        raise ValueError(message)
    if len(headers) != len(sequences):
        message = (
            f"{seed_fasta.name} is not one sequence line per record: "
            f"{len(headers)} headers against {len(sequences)} sequence lines"
        )
        logger.error(message)
        raise ValueError(message)

    contigs = {header[1:]: sequence for header, sequence in zip(headers, sequences, strict=True)}
    filters = json.loads(filter_config.read_text(encoding="utf-8"))
    disallowed = {name: set(partners) for name, partners in filters.items()}

    with destination.open("w", encoding="utf-8") as handle:
        for first, second in itertools.product(contigs, repeat=2):
            if second in disallowed.get(first, ()):
                continue
            handle.write(f">{first}-{second}\n{contigs[second]}{contigs[first]}\n")

    logger.info(f"Merged {len(contigs)} motifs from {seed_fasta.name} into {destination.name}")
    return destination


def _missing_seed_message(output: str, missing: list[str], output_dir: Path) -> str:
    """Compose the one message both the preflight and the derivation loop raise.

    Args:
        output: The derivation output that cannot be produced.
        missing: Seed file names that are not in ``output_dir``.
        output_dir: Reference tree being populated.

    Returns:
        str: The operator-facing message.
    """
    return (
        f"cannot derive {output}: seed(s) {', '.join(missing)} are not present in {output_dir}. "
        "The bundle build stages the seeds into the reference tree before running --from-source."
    )


def _preflight_literal_seeds(install_config: dict[str, Any], output_dir: Path) -> None:
    """Fail on a missing literal-derivation seed before a single genome is downloaded.

    :func:`run_derivations` already names a missing ``filter_config.json`` well, but only
    after every genome has been downloaded, decompressed and BWA-indexed - three hours
    into a build that a one-line edit to the workflow's staging step could have broken.
    This file's philosophy is preflight first; this is the same check, moved to where it
    costs a second, and it uses the same message so the two cannot drift apart.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree being populated.

    Raises:
        RuntimeError: If a ``literal`` derivation names a seed that is not present.
    """
    for spec in install_config.get("derivations", []):
        if spec.get("kind") != "literal":
            continue
        missing = [name for name in spec.get("from_seeds", []) if not (output_dir / name).exists()]
        if missing:
            message = _missing_seed_message(spec.get("output", "<unnamed derivation>"), missing, output_dir)
            logger.error(message)
            raise RuntimeError(message)


def run_derivations(
    install_config: dict[str, Any],
    output_dir: Path,
    references: dict[str, Path],
    selected: set[str] | None = None,
) -> None:
    """Run every in-scope derivation and verify it against its committed digest.

    Two kinds are configured. ``shark`` cuts a region out of an installed chromosome
    FASTA; ``literal`` merges the two MUC1 motif seeds. Every output is checked against
    ``expected_sha256`` before it is indexed, so a silent change upstream turns the
    bundle build red instead of publishing different sequence under an unchanged name.

    A ``shark`` derivation whose ``from`` is outside ``selected`` is skipped, not failed:
    ``--references hg19`` legitimately builds no hg38 reference, and because installation
    merges rather than replaces, a later ``--references hg38`` run fills it in beside the
    hg19 tree. A source that *is* selected but did not arrive remains a hard error. The
    ``literal`` derivation depends on seeds rather than on a genome, so it always runs.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree being populated.
        references: Physical id to installed chromosome FASTA, for ``from`` lookups.
        selected: The reference ids this run was asked for. None means every configured
            source is in scope, which is what a caller with no selection wants.

    Raises:
        RuntimeError: If a selected source reference, or a literal seed, is missing.
        ValueError: If a derivation declares no digest or an unknown kind, if a literal
            derivation does not name exactly two seeds, or if a derived file does not
            match ``expected_sha256``.
    """
    samtools = install_config.get("samtools_path", "samtools")

    for spec in install_config.get("derivations", []):
        kind = spec.get("kind")
        output = spec["output"]
        destination = output_dir / output
        expected = spec.get("expected_sha256")
        if not expected:
            # Checked before the scope test on purpose: a derivation missing its digest is
            # a malformed config, and that should surface on any run, not only the one that
            # happens to select its source.
            message = f"derivation {output} declares no expected_sha256; refusing to produce it unverified"
            logger.error(message)
            raise ValueError(message)

        if kind == "shark":
            source_id = spec["from"]
            if selected is not None and source_id not in selected:
                logger.info(
                    f"Skipping {output}: its source '{source_id}' is not in this run's reference "
                    f"selection. Install it later with --references {source_id}."
                )
                continue
            source = references.get(source_id)
            if source is None or not source.exists():
                message = f"cannot derive {output}: source reference '{source_id}' is not installed"
                logger.error(message)
                raise RuntimeError(message)
            derive_region_fasta(source, spec["region"], destination, samtools)
            provenance = f"{source.name} at {spec['region']}"
        elif kind == "literal":
            seed_names = spec.get("from_seeds", [])
            if len(seed_names) != 2:
                message = (
                    f"derivation {output} must name exactly two seeds (motif FASTA, filter config); "
                    f"got {len(seed_names)}"
                )
                logger.error(message)
                raise ValueError(message)
            seeds = [output_dir / name for name in seed_names]
            missing = [seed.name for seed in seeds if not seed.exists()]
            if missing:
                message = _missing_seed_message(output, missing, output_dir)
                logger.error(message)
                raise RuntimeError(message)
            merge_pairwise_motifs(seeds[0], seeds[1], destination)
            provenance = " and ".join(seed.name for seed in seeds)
        else:
            message = f"derivation {output} declares unknown kind '{kind}'; expected 'shark' or 'literal'"
            logger.error(message)
            raise ValueError(message)

        try:
            verify_sha256(destination, expected)
        except ValueError as mismatch:
            _discard_failed_output(destination)
            message = f"{mismatch}; discarded {output} rather than leave a wrong reference in the tree"
            logger.error(message)
            raise ValueError(message) from mismatch
        index_fasta_with_samtools(destination, samtools)
        logger.info(f"Derived {output} from {provenance}")


def resolve_source_location(ref_id: str, entry: dict[str, Any], release_spec: dict[str, Any] | None) -> tuple[str, str]:
    """Decide where one genome is fetched from and what it must hash to.

    A release spec's ``sources`` block may **fill in** a field the shipped config leaves
    blank, but it may not contradict one. The trust anchor is
    install_references_config.json, which lives in this repository and is a base-image
    content-hash input; a spec lives beside the assets it describes, so letting it win
    would let a release publish bytes that VNtyper's own committed provenance does not
    describe. Moving a future release to a newer upstream therefore takes a reviewed
    VNtyper commit. ``scripts/bundle_release.py`` runs the same comparison in its
    ``--check-spec-only`` preflight, so the disagreement surfaces in minute one rather
    than three hours in.

    Args:
        ref_id: Physical reference id, e.g. ``hg19``.
        entry: That reference's section of the install config.
        release_spec: Parsed ``--release-spec`` contents, or None.

    Returns:
        tuple[str, str]: The URL to download and the expected SHA-256.

    Raises:
        ValueError: If the spec and the config both declare a field and disagree, or if
            either the URL or the digest is missing. ``--from-source`` never installs
            bytes it cannot verify.
    """
    pinned: dict[str, Any] = (release_spec or {}).get("sources", {}).get(ref_id, {})

    for field in ("url", "source_sha256"):
        declared = pinned.get(field)
        committed = entry.get(field)
        if declared and committed and declared != committed:
            message = (
                f"reference '{ref_id}': the release spec pins {field} {declared!r} but "
                f"install_references_config.json pins {committed!r}. The committed config is the "
                "trust anchor; change it in a reviewed commit and match the spec to it rather than "
                "publishing a release whose provenance does not describe its bytes."
            )
            logger.error(message)
            raise ValueError(message)

    url = pinned.get("url") or entry.get("url")
    digest = pinned.get("source_sha256") or entry.get("source_sha256")

    if not url:
        message = f"no source URL configured for reference '{ref_id}'; --from-source cannot build it"
        logger.error(message)
        raise ValueError(message)
    if not digest:
        message = (
            f"no source_sha256 configured for reference '{ref_id}' ({url}); "
            "--from-source refuses to install unverified bytes"
        )
        logger.error(message)
        raise ValueError(message)
    return url, digest


def decompress_source(archive: Path) -> Path:
    """Expand a downloaded genome, or hand back a file that needs no expanding.

    Decompression is unconditional rather than skipped when the target exists: the
    archive has just been digest-verified, so rewriting from it is the only way to be
    sure the FASTA beside it is the file that digest describes.

    Args:
        archive: The verified download.

    Returns:
        Path: The FASTA to index and derive from.
    """
    if archive.suffix != ".gz":
        logger.info(f"{archive.name} is not gzip-compressed; using it where it landed")
        return archive

    destination = archive.with_suffix("")
    with gzip.open(archive, "rb") as compressed, destination.open("wb") as expanded:
        shutil.copyfileobj(compressed, expanded)
    logger.info(f"Expanded {archive.name} to {destination.name}")
    return destination


def install_from_source(
    install_config: dict[str, Any],
    output_dir: Path,
    references: list[str],
    aligners: dict[str, dict[str, Any]],
    index_threads: int,
    release_spec: dict[str, Any] | None = None,
    skip_indexing: bool = False,
) -> dict[str, Path]:
    """Build every requested reference from its upstream source.

    This is what the bundle build workflow runs, so it is the definition of what a
    bundle contains - there is no second derivation implementation to drift from it.

    When ``release_spec`` is given (release builds), the URLs and digests it names take
    precedence and the downloaded bytes are verified **before** decompression or
    indexing. Without it, the URLs and ``source_sha256`` values in
    install_references_config.json are used; those pin Ensembl to an explicit release
    rather than the mutable ``current`` path.

    The two downloadable MUC1 seeds are installed too, minus anything a derivation
    produces, as is the adVNTR database. Those are common assets rather than selectable
    ones, so ``--references`` does not gate them - see the comment at the call site.
    ``filter_config.json`` is not downloadable and must already be in the tree - the
    bundle build stages the seeds before calling this, and a checkout has them tracked.
    :func:`_preflight_literal_seeds` says so by name **before** the first download, so a
    staging-step mistake costs a second rather than three hours.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree to populate.
        references: Physical reference ids to build.
        aligners: Enabled aligner configurations; empty to skip aligner indexing.
        index_threads: Threads handed to the aligners.
        release_spec: Parsed ``--release-spec`` contents, or None.
        skip_indexing: Skip the optional indexing of the seeds and the adVNTR database.
            It does **not** skip ``samtools faidx`` on the chromosome FASTAs: the region
            derivations cut from those indexes, so that one is a prerequisite rather than
            an optional extra.

    Returns:
        dict[str, Path]: physical id -> installed chromosome FASTA, for
        :func:`run_derivations`.

    Raises:
        RuntimeError: If a download fails, if any aligner index fails, or if a selected
            derivation source or literal seed is missing.
        ValueError: If an entry has no ``target_path``, no URL or no digest, if the
            release spec contradicts a committed pin, or if a download or derivation does
            not match its expected digest.
    """
    samtools = install_config.get("samtools_path", "samtools")
    selected = set(references)
    installed: dict[str, Path] = {}
    _preflight_literal_seeds(install_config, output_dir)

    for section in GENOME_SECTIONS:
        for ref_id, entry in install_config.get(section, {}).items():
            if ref_id not in selected:
                continue

            target_path_str = entry.get("target_path")
            if not target_path_str:
                message = f"reference '{ref_id}' declares no target_path; cannot install it from source"
                logger.error(message)
                raise ValueError(message)

            url, expected = resolve_source_location(ref_id, entry, release_spec)
            archive = output_dir / target_path_str

            logger.info(f"Building {ref_id} from source")
            try:
                download_file(url, archive)
            except SystemExit as exit_signal:
                # `download_file` ends the process, which is right for the legacy CLI path
                # but wrong here: this function is called by the bundle build and by
                # `staged_install`, both of which have cleanup to do first.
                message = f"failed to download {ref_id} from {url}"
                logger.error(message)
                raise RuntimeError(message) from exit_signal
            # Before decompression and before any indexing: unverified bytes are never
            # expanded, never indexed and never derived from. A mismatch also takes the
            # archive with it - `download_file` skips a destination that already exists, so
            # leaving truncated or tampered bytes on disk would make every later run skip
            # the download, re-hash the same bad file and fail identically forever.
            try:
                verify_sha256(archive, expected)
            except ValueError as mismatch:
                archive.unlink(missing_ok=True)
                message = f"{mismatch}; removed {archive.name} so a retry downloads it again"
                logger.error(message)
                raise ValueError(message) from mismatch

            fasta = decompress_source(archive)
            index_fasta_with_samtools(fasta, samtools)
            if aligners:
                results = index_reference_with_aligners(fasta, aligners, threads=index_threads, force_reindex=False)
                # The legacy path logs a failure and carries on, which is right for an
                # operator installing locally. It is wrong here. `execute_aligner_index`
                # returns False on a CalledProcessError, and `bundle_release.verify_tree`
                # only checks that each sidecar `is_file()` - so `bwa index` running out
                # of disk 7 GB into a 14 GB runner would leave some sidecar names present,
                # let the build continue, and archive a truncated index into an immutable
                # release that cannot be patched.
                failed = sorted(name for name, ok in results.items() if not ok)
                if failed:
                    message = (
                        f"indexing {ref_id} ({fasta.name}) failed for aligner(s) {', '.join(failed)}. "
                        "--from-source refuses to package an incomplete index: the sidecar files a failed "
                        "run leaves behind still exist, so nothing downstream would notice."
                    )
                    logger.error(message)
                    raise RuntimeError(message)

            installed[ref_id] = fasta

    _install_source_seeds(install_config, output_dir, skip_indexing)

    # Unconditional, unlike the legacy path, which filters this by `--references` and so
    # installs no adVNTR database at all unless someone asks for `vntr_db_advntr` by name -
    # and nobody does, because it is not an assembly. The bundle build runs this function
    # with a list of six genomes, so keeping that filter here would publish an immutable
    # refs release with no adVNTR databases in it. The adVNTR database is a common asset
    # every install needs, exactly like the MUC1 seeds above; Task 10 makes that explicit
    # for both paths.
    vntyper_refs = install_config.get("vntyper_references", {})
    if vntyper_refs:
        logger.info("Installing common VNtyper references (not selectable: every install needs them)...")
        process_vntyper_references(vntyper_refs, output_dir, install_config.get("bwa_path", "bwa"), skip_indexing, {})

    run_derivations(install_config, output_dir, installed, selected)
    return installed


def _install_source_seeds(install_config: dict[str, Any], output_dir: Path, skip_indexing: bool = False) -> None:
    """Fetch the common seed files, skipping anything a derivation produces.

    ``All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`` is still listed as a
    downloadable raw file, but it is derived now, so fetching it first would be a
    wasted round trip over a file that is about to be overwritten.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree to populate.
        skip_indexing: Skip the seeds' ``samtools faidx`` step.
    """
    own_repo_refs = install_config.get("own_repository_references", {})
    raw_files = own_repo_refs.get("raw_files", [])
    if not raw_files:
        return

    derived = {spec.get("output") for spec in install_config.get("derivations", [])}
    seeds = [entry for entry in raw_files if entry.get("target_path") not in derived]
    if not seeds:
        return

    logger.info("Installing common seed files...")
    process_own_repository_references({"raw_files": seeds}, output_dir, skip_indexing, {})


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
        logger.info(f"MD5 checksums written to {checksum_file}")
    except Exception as e:
        logger.error(f"Failed to write MD5 checksums to {checksum_file}: {e}")
        sys.exit(1)


def setup_logging(output_dir: Path):
    """
    Setup logging to output to both stdout and a log file in the output directory.

    Args:
        output_dir (Path): Directory where logs will be stored.
    """
    log_file = output_dir / "install_references.log"

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    c_handler = logging.StreamHandler(sys.stdout)
    f_handler = logging.FileHandler(log_file)

    c_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    c_handler.setFormatter(formatter)
    f_handler.setFormatter(formatter)

    root_logger.addHandler(c_handler)
    root_logger.addHandler(f_handler)

    logger.info(f"Logging initialized. Logs will be saved to {log_file}")


def main(
    output_dir: Path,
    config_path: Path | None = None,
    skip_indexing: bool = False,
    index_threads: int = 4,
    aligners_to_use: list[str] | None = None,
    references_to_process: list[str] | None = None,
    from_source: bool = False,
    release_spec_path: Path | None = None,
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
        from_source (bool): Build every reference from its upstream source and run the
            configured derivations, instead of the legacy per-section download. This is
            the path the bundle build workflow runs. Defaults to False, which leaves the
            existing behaviour untouched.
        release_spec_path (Optional[Path]): Release builds only - take every source URL
            and digest from this file rather than from the shipped install config.
    """
    script_dir = Path(__file__).parent
    install_config_path = script_dir / "install_references_config.json"

    install_config = load_install_config(install_config_path)

    ucsc_refs = install_config.get("ucsc_references", {})
    ncbi_refs = install_config.get("ncbi_references", {})
    ensembl_refs = install_config.get("ensembl_references", {})
    vntyper_refs = install_config.get("vntyper_references", {})
    own_repo_refs = install_config.get("own_repository_references", {})
    bwa_path = install_config.get("bwa_path", "bwa")  # Default to 'bwa'

    # Load aligner configurations
    aligner_config = install_config.get("aligners", {})

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create output directory {output_dir}: {e}")
        sys.exit(1)

    setup_logging(output_dir)

    # Filter references (after logging is set up)
    # Default to hg19 and hg38 (UCSC) for backward compatibility
    if references_to_process is None:
        references_to_process = ["hg19", "hg38"]
        logger.info("No references specified, using default: hg19, hg38")

    all_available_refs = (
        set(ucsc_refs.keys()) | set(ncbi_refs.keys()) | set(ensembl_refs.keys()) | set(vntyper_refs.keys())
    )
    requested_refs = set(references_to_process)
    found_refs = requested_refs & all_available_refs
    missing_refs = requested_refs - all_available_refs

    if missing_refs:
        logger.warning(f"Requested references not found in config: {', '.join(sorted(missing_refs))}")
        logger.warning(f"Available references: {', '.join(sorted(all_available_refs))}")

    if not found_refs:
        logger.error("None of the requested references were found in the configuration.")
        sys.exit(1)

    logger.info(f"Processing references: {', '.join(sorted(found_refs))}")

    ucsc_refs = {k: v for k, v in ucsc_refs.items() if k in references_to_process}
    ncbi_refs = {k: v for k, v in ncbi_refs.items() if k in references_to_process}
    ensembl_refs = {k: v for k, v in ensembl_refs.items() if k in references_to_process}
    vntyper_refs = {k: v for k, v in vntyper_refs.items() if k in references_to_process}

    # Initialize aligners
    enabled_aligners = {}
    if not skip_indexing and aligner_config:
        logger.info("=" * 80)
        logger.info("ALIGNER CONFIGURATION")
        logger.info("=" * 80)
        logger.info("Checking available aligners:")

        # Get enabled aligners
        all_enabled = get_enabled_aligners(aligner_config)

        # Filter by user-specified aligners if provided, otherwise default to BWA only
        if aligners_to_use:
            for aligner_name in aligners_to_use:
                if aligner_name in all_enabled:
                    enabled_aligners[aligner_name] = all_enabled[aligner_name]
                elif aligner_name in aligner_config:
                    logger.warning(f"  ✗ {aligner_name} was specified but is not available or not enabled")
                else:
                    logger.error(f"  ✗ Unknown aligner: {aligner_name}")
        else:
            # Default to BWA only
            if "bwa" in all_enabled:
                enabled_aligners["bwa"] = all_enabled["bwa"]
                logger.info("  Using default aligner: bwa")
            else:
                logger.warning("  Default aligner 'bwa' not available")

        if not enabled_aligners:
            logger.warning("No aligners available. Indexing will be skipped.")
            logger.warning(
                "To enable aligners, install them and ensure they are in your PATH, "
                "or set enabled:true in install_references_config.json"
            )
        else:
            logger.info(f"\nWill use {len(enabled_aligners)} aligner(s) for indexing:")
            for name in enabled_aligners:
                logger.info(f"  • {name}")

            # Detect index file conflicts
            conflicts = detect_index_conflicts(enabled_aligners)
            if conflicts:
                logger.warning("\n⚠ Index file conflicts detected:")
                for warning in conflicts:
                    logger.warning(f"  • {warning}")
                logger.warning("  These conflicts may cause issues if aligners overwrite each other's index files.")
            else:
                logger.info("  ✓ No index file conflicts detected")

        logger.info("=" * 80)
        logger.info("")

    md5_dict: dict[str, str] = {}

    if from_source:
        # The source path owns every section at once: it must verify each download
        # before expanding it, and it must know which chromosome FASTAs landed where so
        # the derivations can cut from them. It writes no md5_checksums.txt - the
        # digests it enforces are SHA-256 values committed in the install config.
        release_spec = load_install_config(release_spec_path) if release_spec_path else None
        if release_spec is not None:
            logger.info(f"Taking source URLs and digests from {release_spec_path}")
        install_from_source(
            install_config,
            output_dir,
            sorted(found_refs),
            enabled_aligners,
            index_threads,
            release_spec,
            skip_indexing,
        )
    else:
        # Process UCSC references
        if ucsc_refs:
            logger.info("Processing UCSC references...")
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
            logger.info("Processing NCBI references...")
            process_ucsc_references(
                ncbi_refs,
                output_dir,
                bwa_path,
                skip_indexing,
                md5_dict,
                aligners=enabled_aligners,
                index_threads=index_threads,
            )

        # Process ENSEMBL references
        if ensembl_refs:
            logger.info("Processing ENSEMBL references...")
            process_ucsc_references(
                ensembl_refs,
                output_dir,
                bwa_path,
                skip_indexing,
                md5_dict,
                aligners=enabled_aligners,
                index_threads=index_threads,
            )

        # Process VNtyper references
        if vntyper_refs:
            logger.info("Processing VNtyper references...")
            process_vntyper_references(vntyper_refs, output_dir, bwa_path, skip_indexing, md5_dict)

        # Process own repository references
        if own_repo_refs:
            logger.info("Processing own repository references...")
            process_own_repository_references(own_repo_refs, output_dir, skip_indexing, md5_dict)

    # Write MD5 checksums to file
    if md5_dict:
        write_md5_checksums(md5_dict, output_dir)

    # Update the main config.json with new reference paths if config_path is provided
    if config_path and config_path.exists():
        updated_references = {}

        # Collect all references from UCSC
        for ref_key, ref_info in ucsc_refs.items():
            ucsc_target = ref_info.get("target_path")
            if ucsc_target:
                ref_path = output_dir / ucsc_target
                updated_references[f"ucsc_{ref_key}"] = ref_path.resolve()

        # Collect all references from NCBI
        for ref_key, ref_info in ncbi_refs.items():
            ncbi_target = ref_info.get("target_path")
            if ncbi_target:
                ref_path = output_dir / ncbi_target
                updated_references[f"ncbi_{ref_key}"] = ref_path.resolve()

        # Collect all references from ENSEMBL
        for ref_key, ref_info in ensembl_refs.items():
            ensembl_target = ref_info.get("target_path")
            if ensembl_target:
                ref_path = output_dir / ensembl_target
                updated_references[f"ensembl_{ref_key}"] = ref_path.resolve()

        # Collect all references from VNtyper
        for ref_key, ref_info in vntyper_refs.items():
            vntyper_target = ref_info.get("target_path")
            if vntyper_target:
                ref_path = output_dir / vntyper_target
                updated_references[f"vntyper_{ref_key}"] = ref_path.resolve()

        # Collect references from own repository
        raw_files: list[dict[str, str]] = own_repo_refs.get("raw_files", [])
        for file_info in raw_files:
            own_target = file_info.get("target_path")
            if own_target:
                ref_name = Path(own_target).stem
                ref_path = output_dir / own_target
                updated_references[f"own_repo_{ref_name}"] = ref_path.resolve()

        update_config(config_path, updated_references)
    else:
        if config_path:
            logger.warning(f"Config file {config_path} does not exist. Skipping config update.")
        else:
            logger.info("No config_path provided. Skipping config update.")

    logger.info("All references have been installed and configured successfully.")


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
        help="Specific references to process (e.g., hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl). "
        "Default: hg19 hg38 (UCSC references only).",
    )
    parser.add_argument(
        "--from-source",
        action="store_true",
        help="Build references from their upstream sources and run the configured derivations. "
        "Needs the MUC1 seed files already present in the output directory.",
    )
    parser.add_argument(
        "--release-spec",
        type=Path,
        default=None,
        help="Release builds only: take every source URL and digest from this file. Requires --from-source.",
    )

    args = parser.parse_args()
    if args.release_spec and not args.from_source:
        parser.error("--release-spec requires --from-source; it only affects the from-source build")
    main(
        args.output_dir,
        args.config_path,
        args.skip_indexing,
        args.threads,
        args.aligners,
        args.references,
        args.from_source,
        args.release_spec,
    )
