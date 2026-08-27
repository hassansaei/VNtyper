# vntyper/scripts/install_references.py

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# `docker/Dockerfile.base` copies this module and its module-scope siblings into a
# build stage and runs them without installing the package, and every module imported
# here joins the base-image content hash. Anything not provisioned there is imported
# inside the function that needs it.
from vntyper.scripts.install_references_logging import InstallLogHandler, attach_install_log, finish_install_log
from vntyper.scripts.reference_bundle import safe_extract, safe_extract_zip, sha256_of, staged_install, verify_sha256
from vntyper.scripts.reference_download import download_file
from vntyper.scripts.reference_integrity import fetch_verified_asset, verify_existing_asset

logger = logging.getLogger(__name__)

#: Config sections holding one physical chromosome FASTA per entry, in the order
#: `--from-source` and `install_from_bundle`'s `--references` selection both walk them.
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

    No shell is involved - the command is parsed into an argv list, not run through
    `/bin/bash` - so this is a correctness bug, not an injection one. The path is
    substituted quoted (:func:`_quote`) and the result parsed with `shlex.split`
    rather than `str.split`, so a path containing whitespace (e.g.
    ``--output-dir "/data/my refs"``) still reaches the aligner as one argv element
    instead of being shredded across two.

    Args:
        index_command (str): The indexing command with a placeholder for the file path.
        fasta_path (Path): Path to the FASTA file to index.

    Raises:
        SystemExit: If the indexing fails.
    """
    command = index_command.format(path=_quote(fasta_path))
    logger.info(f"Executing indexing command: {command}")
    try:
        args = shlex.split(command)
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


def canonical_reference_keys(install_config: dict[str, Any], output_dir: Path) -> dict[str, Path]:
    """Map every installed, *verified* reference onto the config key the pipeline reads.

    Reads exactly the schema fixed for the installer and nothing else. Genome keys are
    derived from the registry rather than written by hand, so the writer and the readers
    cannot drift: neither of them owns the name. `reference_registry` is imported here,
    function-locally, rather than at module scope: `docker/Dockerfile.base` copies this
    module alone into a build stage and runs it by path without installing the package,
    and `tests/unit/test_docker_installer_imports.py` fails the build guard if a sibling
    the Dockerfile does not COPY is imported at module scope.

    Parked finding (was #244): a path being *present* under `output_dir` used to be
    enough. `staged_install` seeds a new run's staging directory from whatever tree is
    already there (`reference_bundle.staged_install`, `seed_from_existing=True`) so that
    a partial install (`--references hg19`) does not erase a different assembly a
    previous run installed - that is deliberate accumulation, not a bug. But the same
    seeding step carried forward a file this or any earlier `install-references` run
    never verified at all: an old, hand-copied or tampered FASTA sitting in the output
    directory before the very first real install. This now requires a
    :mod:`reference_provenance` record for the path instead - written at install time by
    `_install_bundle_asset`/`_record_bundle_provenance` (bundle path) and
    `_record_source_provenance` (`--from-source` path) once that specific file's bytes
    have actually been checked against a committed or upstream digest. Presence with no
    record is reported by name, together with the `--references <label>` (or plain
    `install-references`, for a common asset installed on every run) that installs and
    verifies it, and left out of the mapping - neither silently dropped nor silently
    trusted.

    This never re-hashes a file to decide: it only checks whether a ledger key exists,
    so a 700 MB genome FASTA costs nothing extra on every `config.json` write, no matter
    how many times this runs.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Directory the references were installed into.

    Returns:
        dict[str, Path]: Absolute paths, keyed by `reference_data` key.

    Raises:
        KeyError: If an entry is missing a required schema field. Failing here is
            correct - a silently skipped entry is a silently missing reference.
    """
    from vntyper.scripts.reference_provenance import has_record, load_provenance, relative_posix
    from vntyper.scripts.reference_registry import REFERENCE_KINDS, reference_keys

    prefix = REFERENCE_KINDS["bwa"]["prefix"]
    written: dict[str, Path] = {}
    remedies: dict[str, str] = {}
    common_remedy = "vntyper install-references (installs common assets on every run)"

    for section in ("ucsc_references", "ncbi_references", "ensembl_references"):
        for physical_id, entry in install_config.get(section, {}).items():
            key = f"{prefix}_{physical_id}"
            written[key] = (output_dir / entry["installed_path"]).resolve()
            remedies[key] = f"vntyper install-references --references {physical_id}"

    for entry in install_config.get("common_references", []):
        key = entry["config_key"]
        written[key] = (output_dir / entry["installed_path"]).resolve()
        remedies[key] = common_remedy

    for spec in install_config.get("derivations", []):
        if spec["kind"] == "shark":
            (key,) = reference_keys("shark", spec["assembly"])
            remedies[key] = f"vntyper install-references --references {spec['from']}"
        else:
            key = spec["config_key"]
            remedies[key] = common_remedy
        written[key] = (output_dir / spec["output"]).resolve()

    # Only consider files that are actually there. A partial install (`--references
    # hg38`) must not write a key pointing at an hg19 FASTA nobody downloaded - that is
    # the same class of defect as #163 itself, just from the other direction.
    provenance = load_provenance(output_dir)
    result: dict[str, Path] = {}
    for key, path in written.items():
        if not path.exists():
            continue
        try:
            relative = relative_posix(path, output_dir)
        except ValueError:
            relative = None
        if relative is not None and has_record(provenance, relative):
            result[key] = path
        else:
            logger.warning(
                f"{path} exists but has no install provenance record; not writing "
                f"config.json[{key!r}] to point at an unverified file. Run "
                f"`{remedies.get(key, 'vntyper install-references')}` to install and verify it."
            )
    return result


def update_config(config_path: Path, references: dict[str, Path]):
    """
    Update the main config.json with paths to the downloaded references.

    Writes atomically: the new document is written to a sibling `.json.tmp` file and
    then renamed over `config_path` with `os.replace`, which POSIX and Windows both
    guarantee is atomic. A write failure (disk full, a `json.dump` error on a bad value)
    therefore never leaves `config_path` truncated or half-written - the previous config
    survives untouched, and only the `.tmp` file is corrupted.

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

    tmp_path = config_path.with_suffix(".json.tmp")
    try:
        with tmp_path.open("w") as f:
            json.dump(config, f, indent=2)
        os.replace(tmp_path, config_path)
        logger.info(f"Successfully updated {config_path} with new reference paths.")
    except Exception as e:
        logger.error(f"Failed to write updated config.json: {e}")
        tmp_path.unlink(missing_ok=True)
        sys.exit(1)


def process_vntyper_references(
    vntyper_refs: dict[str, dict[str, str]],
    output_dir: Path,
    bwa_path: str,
    skip_indexing: bool,
    md5_dict: dict[str, str],
    release_spec: dict[str, Any] | None = None,
):
    """
    Process VNtyper references by downloading and extracting.

    Args:
        vntyper_refs (dict): Dictionary of VNtyper references.
        output_dir (Path): Base output directory.
        bwa_path (str): Path to the bwa executable.
        skip_indexing (bool): Whether to skip the indexing step.
        md5_dict (dict): Dictionary to store MD5 checksums.
        release_spec (dict, optional): Parsed ``--release-spec`` contents, or None. A
            digest it names must agree with the one committed in
            install_references_config.json (see :func:`resolve_seed_digest`); it is
            never allowed to override it.

    Raises:
        ValueError: If a seed - ``vntr_db_advntr.zip`` is the one this function
            downloads - has no committed or spec-pinned digest, or if the downloaded or
            already-present bytes do not match it (MAJOR 3, milestone-5 PR-2 review).
        RuntimeError: If a verified archive cannot be extracted safely.
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

        expected = resolve_seed_digest(target_path_str, ref_info, release_spec)
        expected = fetch_verified_asset(target_path_str, url, target_path, expected, download_file)
        _record_source_provenance(output_dir, target_path, expected, url)

        md5_checksum = calculate_md5(target_path)
        md5_dict[str(target_path)] = md5_checksum
        logger.info(f"MD5 checksum for {target_path.name}: {md5_checksum}")

        if extract_to:
            extract_dir = output_dir / extract_to
            extract_dir.mkdir(parents=True, exist_ok=True)
            # `canonical_reference_keys` looks common assets like the adVNTR VNTR
            # databases up by their *extracted* path, not the archive's - the archive's
            # own digest is what was actually verified above, so every member this
            # archive is known (from its own listing) to have produced inherits it as
            # its provenance record rather than going unrecorded.
            extracted_members: list[str] = []
            if target_path.suffix == ".zip":
                try:
                    extracted_members = safe_extract_zip(target_path, extract_dir)
                    logger.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    message = f"Failed to extract {target_path}: {e}"
                    logger.error(message)
                    raise RuntimeError(message) from e
            elif target_path.suffixes[-2:] == [".tar", ".gz"] or target_path.suffix == ".tgz":
                try:
                    extracted_members = safe_extract(target_path, extract_dir)
                    logger.info(f"Successfully extracted {target_path.name}")
                except Exception as e:
                    message = f"Failed to extract {target_path}: {e}"
                    logger.error(message)
                    raise RuntimeError(message) from e
            else:
                logger.warning(f"Unsupported archive format for {target_path}. Skipping extraction.")

            for member in extracted_members:
                member_path = extract_dir / member
                if member_path.is_file():
                    _record_source_provenance(output_dir, member_path, expected, url)

        if index_command and not skip_indexing:
            execute_index_command(index_command, target_path)
        elif index_command and skip_indexing:
            logger.info(f"Skipping indexing for {target_path}")


def process_own_repository_references(
    own_repo_refs: dict[str, Any],
    output_dir: Path,
    skip_indexing: bool,
    md5_dict: dict[str, str],
    release_spec: dict[str, Any] | None = None,
):
    """
    Process own repository references by downloading specific FASTA files.

    Args:
        own_repo_refs (dict): Dictionary of own repository references.
        output_dir (Path): Base output directory.
        skip_indexing (bool): Whether to skip the indexing step.
        md5_dict (dict): Dictionary to store MD5 checksums.
        release_spec (dict, optional): Parsed ``--release-spec`` contents, or None. A
            digest it names must agree with the one committed in
            install_references_config.json (see :func:`resolve_seed_digest`); it is
            never allowed to override it.

    Raises:
        ValueError: If a seed this function downloads (``MUC1_motifs_Rev_com.fa``,
            ``code-adVNTR_RUs.fa`` or ``filter_config.json``) has no committed or
            spec-pinned digest, or if the downloaded or already-present bytes do not
            match it (MAJOR 3, milestone-5 PR-2 review).
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

        expected = resolve_seed_digest(target_path_str, file_info, release_spec)
        expected = fetch_verified_asset(target_path_str, url, target_path, expected, download_file)
        _record_source_provenance(output_dir, target_path, expected, url)

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

    try:
        with destination.open("w", encoding="utf-8") as handle:
            for first, second in itertools.product(contigs, repeat=2):
                if second in disallowed.get(first, ()):
                    continue
                handle.write(f">{first}-{second}\n{contigs[second]}{contigs[first]}\n")
    except OSError:
        # ``open("w")`` truncates before the first write, so a failure part-way through --
        # a full disk is the realistic one -- leaves a short FASTA under the final name.
        # Nothing downstream catches that: ``run_derivations`` verifies against
        # ``expected_sha256`` only after this function *returns*, so an exception routes
        # around the check entirely. Discard it here, for the same reason and in the same
        # way :func:`derive_region_fasta` discards its own truncated output.
        _discard_failed_output(destination)
        raise

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


def _downloadable_seed_names(install_config: dict[str, Any]) -> set[str]:
    """Return every seed name `_install_source_seeds` can fetch on its own.

    `own_repository_references.raw_files` is downloaded later in the same
    :func:`install_from_source` run, after the genome loop and before
    :func:`run_derivations`. A seed named there does not need to exist yet when the
    preflight below runs - only a seed with no download source at all does.

    Args:
        install_config: The parsed install_references_config.json.

    Returns:
        set[str]: `target_path` of every `own_repository_references.raw_files` entry.
    """
    own_repo_refs = install_config.get("own_repository_references", {})
    return {entry.get("target_path") for entry in own_repo_refs.get("raw_files", []) if entry.get("target_path")}


def _preflight_literal_seeds(install_config: dict[str, Any], output_dir: Path) -> None:
    """Fail on a missing, non-downloadable literal-derivation seed before any genome runs.

    :func:`run_derivations` already names a missing ``filter_config.json`` well, but only
    after every genome has been downloaded, decompressed and BWA-indexed - three hours
    into a build that a one-line edit to the workflow's staging step could have broken.
    This file's philosophy is preflight first; this is the same check, moved to where it
    costs a second, and it uses the same message so the two cannot drift apart.

    A seed such as ``MUC1_motifs_Rev_com.fa`` is a downloadable ``raw_files`` entry that
    :func:`_install_source_seeds` fetches later in the *same* run, after the genome loop.
    Requiring it here too would reject a fresh, unstaged ``output_dir`` that is about to
    succeed. In the shipped config every literal-derivation seed, including
    ``filter_config.json``, now has a download source (a pinned commit under
    ``berntpopp/vntyper-data``'s ``seeds/``), so this preflight passes for a fresh,
    unstaged ``output_dir`` too - it only still raises for a seed that names no download
    source at all in whatever ``install_config`` a caller supplies.

    ``--derive-only`` does not call this. It downloads and stages nothing, so an absent seed
    there is the tree's shape rather than a fault, and :func:`derive_only` skips such a
    derivation and verifies whatever is already at its path instead.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree being populated.

    Raises:
        RuntimeError: If a ``literal`` derivation names a non-downloadable seed that is
            not present.
    """
    downloadable = _downloadable_seed_names(install_config)
    for spec in install_config.get("derivations", []):
        if spec.get("kind") != "literal":
            continue
        missing = [
            name for name in spec.get("from_seeds", []) if name not in downloadable and not (output_dir / name).exists()
        ]
        if missing:
            message = _missing_seed_message(spec.get("output", "<unnamed derivation>"), missing, output_dir)
            logger.error(message)
            raise RuntimeError(message)


def run_derivations(
    install_config: dict[str, Any],
    output_dir: Path,
    references: dict[str, Path],
    selected: set[str] | None = None,
) -> list[str]:
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

    Returns:
        list[str]: The outputs actually derived and verified, in config order. A caller
        that reports on the run needs this: the skipped ones are the difference between it
        and the configured list, and claiming a file is verified when its derivation was
        skipped is the failure this return value exists to make impossible.

    Raises:
        RuntimeError: If a selected source reference, or a literal seed, is missing.
        ValueError: If a derivation declares no digest or an unknown kind, if a literal
            derivation does not name exactly two seeds, or if a derived file does not
            match ``expected_sha256``.
    """
    samtools = install_config.get("samtools_path", "samtools")
    derived: list[str] = []

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
        _record_source_provenance(output_dir, destination, expected)
        index_fasta_with_samtools(destination, samtools)
        logger.info(f"Derived {output} from {provenance}")
        # Appended only past verify_sha256, so membership means "verified", not "attempted".
        derived.append(output)

    return derived


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


def resolve_seed_digest(name: str, entry: dict[str, Any], release_spec: dict[str, Any] | None) -> str:
    """Decide the SHA-256 a seed file must match, refusing a spec that disagrees.

    Mirrors :func:`resolve_source_location`'s trust model for the four seeds
    ``--from-source`` cannot derive (``MUC1_motifs_Rev_com.fa``, ``code-adVNTR_RUs.fa``,
    ``vntr_db_advntr.zip``, ``filter_config.json``, i.e. ``bundle_release.REQUIRED_SEEDS``):
    the digest committed in install_references_config.json is the trust anchor, and a
    ``--release-spec``'s ``seeds`` block - read directly here rather than via
    ``bundle_release.spec_seed_digests``, since this module may import only siblings
    provisioned by the Docker refs stage - may only corroborate it, never override it.

    Args:
        name: Seed file name, keyed the same way in the install config's
            ``target_path`` and in a release spec's ``seeds`` block (e.g.
            ``"code-adVNTR_RUs.fa"``).
        entry: This seed's section of the install config.
        release_spec: Parsed ``--release-spec`` contents, or None.

    Returns:
        str: Lowercase hex SHA-256 the seed must match, whether freshly downloaded or
        already present in the output tree.

    Raises:
        ValueError: If the spec and the config both declare a digest and disagree, or if
            neither declares one at all - ``--from-source`` never installs a seed it
            cannot verify.
    """
    seeds_block = (release_spec or {}).get("seeds", {})
    spec_entry = seeds_block.get(name) if isinstance(seeds_block, dict) else None
    declared = spec_entry.get("sha256") if isinstance(spec_entry, dict) else spec_entry
    committed = entry.get("source_sha256")

    if declared and committed and declared != committed:
        message = (
            f"seed '{name}': the release spec pins sha256 {declared!r} but "
            f"install_references_config.json pins {committed!r}. The committed config is the "
            "trust anchor; change it in a reviewed commit and match the spec to it rather than "
            "publishing a release whose provenance does not describe its bytes."
        )
        logger.error(message)
        raise ValueError(message)

    digest = committed or declared
    if not digest:
        message = f"no source_sha256 configured for seed '{name}'; --from-source refuses to install unverified bytes"
        logger.error(message)
        raise ValueError(message)
    return digest


def verify_seed(name: str, target_path: Path, entry: dict[str, Any], release_spec: dict[str, Any] | None) -> str:
    """Verify a seed against its committed digest, whether just downloaded or reused.

    ``download_file`` skips a destination that already exists, so without this a stale
    or corrupted seed left over from an earlier partial run - or planted there - would be
    extracted, indexed and activated without complaint (MAJOR 3, milestone-5 PR-2
    review). Calling this unconditionally after every ``download_file`` closes that gap
    for both the freshly-downloaded and the already-present case alike.

    Args:
        name: Seed file name, as used by :func:`resolve_seed_digest`.
        target_path: Where the seed landed in the output tree.
        entry: This seed's section of the install config.
        release_spec: Parsed ``--release-spec`` contents, or None.

    Returns:
        str: The digest `target_path` was verified against - both callers record it as
        this file's install provenance (see :mod:`reference_provenance`) rather than
        re-deriving it.

    Raises:
        ValueError: If no digest is configured, or if the file's bytes do not match it.
            On a mismatch the file is removed. The download call sites use
            :func:`reference_integrity.fetch_verified_asset` when a stale pre-existing
            file should be replaced within the same run.
    """
    expected = resolve_seed_digest(name, entry, release_spec)
    return verify_existing_asset(name, target_path, expected)


def _record_source_provenance(output_dir: Path, path: Path, sha256: str, source_url: str | None = None) -> None:
    """Record one ``--from-source`` install provenance entry.

    A thin wrapper around :mod:`reference_provenance`, imported function-locally like
    every other sibling `install_references.py` uses only inside a function body (see
    that module's docstring for why). Called once per verified file rather than
    threaded through a shared mutable accumulator: each call is an independent,
    idempotent read-modify-write of the small ledger, which is cheap next to the
    network and hashing cost already paid to produce `sha256`.

    Args:
        output_dir: Reference-tree root the file was installed under (the
            `staged_install` staging directory, mid-install).
        path: The installed file.
        sha256: The digest already verified for this file (or, for a decompressed
            genome FASTA, for the compressed source archive it was expanded from).
        source_url: The upstream URL the digest was verified against, when there is
            one. None for a locally derived output (a SHARK region cut, or the merged
            motif FASTA), which has no single upstream URL of its own.
    """
    from vntyper.scripts.reference_provenance import build_record, merge, relative_posix

    relative = relative_posix(path, output_dir)
    merge(output_dir, {relative: build_record(sha256=sha256, source="from-source", source_url=source_url)})


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
                fetch_verified_asset(archive.name, url, archive, expected, download_file)
            except RuntimeError as failure:
                # `download_file` raises rather than exiting; this re-wrap only adds the
                # reference id for the operator.
                message = f"failed to download {ref_id} from {url}"
                logger.error(message)
                raise RuntimeError(message) from failure
            fasta = decompress_source(archive)
            index_fasta_with_samtools(fasta, samtools)
            if aligners:
                # Parked finding (was #245): `force_reindex=True`, not the module
                # default. `decompress_source` above unconditionally rewrites `fasta`
                # from the just-verified archive every time this loop reaches it (see
                # its own docstring), but `staged_install` seeds a fresh staging
                # directory from whatever aligner sidecars (`.amb`/`.ann`/`.bwt`/`.pac`/
                # `.sa`) an earlier install left beside the *old* FASTA at this same
                # path. Without this, `check_index_exists` sees all five present and
                # skips indexing entirely - so after a maintainer repins this genome's
                # `source_sha256`, the tree ends up with a new FASTA and a `bwa mem`
                # index still built from the old one, both named identically by
                # `config.json` and by this run's own install provenance record. This
                # call site is exactly "--from-source actually replaced the FASTA in
                # this run" (unconditionally, right above), so forcing a reindex here
                # is defensible even though it costs a redundant index on the common
                # case where the bytes did not change.
                results = index_reference_with_aligners(fasta, aligners, threads=index_threads, force_reindex=True)
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

            # `expected` verified the compressed *archive*; `decompress_source` is a
            # deterministic gunzip, but the provenance record should describe the actual
            # installed artifact, not assume the expansion produced it byte-for-byte -
            # hence one fresh digest of `fasta` itself rather than reusing `expected`.
            _record_source_provenance(output_dir, fasta, sha256_of(fasta), url)
            installed[ref_id] = fasta

    _install_source_seeds(install_config, output_dir, skip_indexing, release_spec)

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
        process_vntyper_references(
            vntyper_refs,
            output_dir,
            install_config.get("bwa_path", "bwa"),
            skip_indexing,
            {},
            release_spec=release_spec,
        )

    run_derivations(install_config, output_dir, installed, selected)
    return installed


def installed_reference_map(install_config: dict[str, Any], output_dir: Path) -> dict[str, Path]:
    """Map every configured reference id to its installed FASTA, keeping only what exists.

    The install paths are declared per reference as ``installed_path``, so this reads the
    same field the installers write rather than reconstructing a filename.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree to inspect.

    Returns:
        dict[str, Path]: Reference id to installed FASTA, for the ids actually present.
    """
    found: dict[str, Path] = {}
    for group in ("ucsc_references", "ncbi_references", "ensembl_references"):
        for ref_id, spec in install_config.get(group, {}).items():
            installed = spec.get("installed_path")
            if not installed:
                continue
            fasta = output_dir / installed
            if fasta.is_file():
                found[ref_id] = fasta
    return found


def _verify_files_that_could_not_be_rebuilt(skipped: list[dict[str, Any]], output_dir: Path) -> None:
    """Check what is already at the path of a derivation this tree cannot rebuild.

    A skipped derivation writes nothing, so nothing unverified is *produced*. But a file
    may still be at that path from an earlier install, and leaving it unread would make
    ``--derive-only`` answer "are my derived files right?" with silence for exactly the
    files it could not rebuild -- while exiting 0. The digest is committed and the file is
    small, so there is no reason not to look.

    A mismatch is discarded rather than left in place, for the reason
    :func:`run_derivations` gives when it discards its own: a wrong reference produces a
    plausible result rather than an obvious failure, and it cannot be rebuilt here.

    Args:
        skipped: Derivation specs that were not rebuilt on this run.
        output_dir: Reference tree being inspected.

    Raises:
        ValueError: If a file present at a skipped derivation's path does not match its
            committed digest.
    """
    matching: list[str] = []
    absent: list[str] = []

    for spec in skipped:
        output = spec["output"]
        destination = output_dir / output
        if not destination.is_file():
            absent.append(output)
            continue
        try:
            verify_sha256(destination, spec["expected_sha256"])
        except ValueError as mismatch:
            _discard_failed_output(destination)
            message = (
                f"{mismatch}; discarded {output} rather than leave a wrong reference in the tree. "
                "It could not be rebuilt on this run because its source genome is not installed."
            )
            logger.error(message)
            raise ValueError(message) from mismatch
        matching.append(output)

    if matching:
        logger.info(
            "Of those, already present and matching their committed digests: %s",
            ", ".join(matching),
        )
    if absent:
        logger.warning(
            "Of those, missing from the tree: %s. Install what it is built from -- the source "
            "genome with `vntyper install-references --references <id>`, or the seeds -- and "
            "run this again, or install the published bundle, which ships them pre-built.",
            ", ".join(absent),
        )


def derive_only(install_config: dict[str, Any], output_dir: Path) -> None:
    """Rebuild the derived reference files from what is already on disk. No network.

    Three files in the reference tree are *derived* rather than downloaded: the two MUC1
    region FASTAs, cut out of an installed chromosome with ``samtools faidx``, and the
    merged MUC1 motif FASTA, built from two seeds. The published bundle ships them
    pre-built and ``--from-source`` builds them at the end of its run, so on the two
    ordinary paths they are already correct.

    This exists for the case in between, which had no command at all: a tree whose genomes
    and seeds are present but whose derived files are missing or suspect. Recovering that
    previously meant a full ``--from-source`` run, which re-downloads and BWA-indexes six
    chromosome FASTAs to regenerate three small files -- so in practice it was done by hand,
    with `samtools faidx` typed from the config, and by hand is exactly where an unverified
    reference file comes from.

    Every output is still checked against its committed ``expected_sha256`` by
    :func:`run_derivations`, so this is not a way to produce unverified bytes; it is the
    same verification on a cheaper path.

    **A derivation this tree cannot rebuild is skipped, not failed** -- for either reason,
    a missing source genome or missing seeds. A tree holding only hg19 legitimately derives
    only the hg19 region; and the published bundle ships the merged motif FASTA pre-built
    without staging ``filter_config.json`` beside it, so on the tree the Docker image
    actually carries, the literal derivation has no seeds and never will. Failing there
    would make this command unusable on the commonest tree in existence -- which is what it
    did until a run inside the image showed it.

    Skipping is safe because it is not silence. The closing summary names what was not
    rebuilt, and whatever is already at those paths is verified against the same committed
    digest by :func:`_verify_files_that_could_not_be_rebuilt`, so a stale file left by an
    earlier install is caught here rather than by a genotyping run months later.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: An existing reference tree.

    Raises:
        ValueError: If a derivation's output, or a file already present at the path of a
            derivation this tree cannot rebuild, does not match its committed digest.
    """
    references = installed_reference_map(install_config, output_dir)
    logger.info(
        "Deriving reference files from %d installed genome(s) in %s: %s",
        len(references),
        output_dir,
        ", ".join(sorted(references)) or "none",
    )

    # `run_derivations` skips a *shark* derivation whose source is outside `selected` and
    # raises for a *literal* one whose seeds are absent -- correct for `--from-source`, which
    # stages the seeds itself and so knows a missing one is a real fault. Nothing here
    # downloads or stages anything, so an absent seed is this tree's shape rather than a
    # fault. Filter those out first; the scope rule for genomes stays where it is, because
    # duplicating it here is how the two would drift.
    specs = install_config.get("derivations", [])
    runnable = []
    for spec in specs:
        missing_seeds = [name for name in spec.get("from_seeds", []) if not (output_dir / name).is_file()]
        if spec.get("kind") == "literal" and missing_seeds:
            logger.info(
                "Skipping %s: seed(s) %s are not in this tree. The published bundle ships this "
                "file pre-built and does not stage its seeds beside it, so this is the normal "
                "shape of a bundle-installed tree rather than a fault.",
                spec["output"],
                ", ".join(missing_seeds),
            )
            continue
        runnable.append(spec)

    # `selected` is the set of genomes actually present, which routes a derivation whose
    # source is absent down run_derivations' skip path instead of its hard-error path.
    derived = run_derivations({**install_config, "derivations": runnable}, output_dir, references, set(references))

    skipped = [spec for spec in specs if spec["output"] not in derived]
    if not skipped:
        logger.info(
            "Derived and verified all %d reference file(s) against their committed digests.",
            len(derived),
        )
        return

    logger.warning(
        "Derived and verified %d of %d reference file(s) in %s. Not rebuilt in this tree: %s.",
        len(derived),
        len(specs),
        output_dir,
        ", ".join(spec["output"] for spec in skipped),
    )
    _verify_files_that_could_not_be_rebuilt(skipped, output_dir)


def _install_source_seeds(
    install_config: dict[str, Any],
    output_dir: Path,
    skip_indexing: bool = False,
    release_spec: dict[str, Any] | None = None,
) -> None:
    """Fetch the common seed files, skipping anything a derivation produces.

    ``All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`` is derived, not downloaded, so
    fetching it would be a wasted round trip over a file that is about to be overwritten.
    Its ``raw_files`` entry was removed along with the dead URL it carried (#253); the filter
    below stays, because it enforces the rule rather than the one entry.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Reference tree to populate.
        skip_indexing: Skip the seeds' ``samtools faidx`` step.
        release_spec: Parsed ``--release-spec`` contents, or None. Forwarded to
            :func:`process_own_repository_references` so a seed's digest can be
            corroborated, never overridden, by the spec.
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
    process_own_repository_references({"raw_files": seeds}, output_dir, skip_indexing, {}, release_spec=release_spec)


###############################################################################
# Installing the published, checksummed reference bundle
###############################################################################

#: Matches bwa's own banner line, e.g. "Version: 0.7.18-r1243-dirty". The identical
#: pattern lives in `scripts/bundle_release.py::_BWA_VERSION`, which is not importable
#: here (see the module-scope import comment at the top of this file) - a three-line
#: literal is cheaper to keep in sync than a cross-file dependency neither side needs.
_BWA_VERSION_RE = re.compile(r"^Version:\s*(\S+)", re.MULTILINE)

#: The two files `scripts/bundle_release.py` writes inside every asset (never fetched
#: separately - see `install_from_bundle`'s docstring). They describe exactly one asset
#: each, so they are deleted from the shared staging tree right after being read: left
#: in place, the *last* asset's copy would silently overwrite every earlier one's, and
#: both would otherwise land inside the installed reference tree they were never part of.
_BUNDLE_METADATA_MEMBERS = ("release-manifest.json", "BUILD_INFO.json")


@dataclass(frozen=True)
class _BundleAsset:
    """One published `.tar.gz` this run needs, and what it becomes once extracted.

    Attributes:
        name: Asset file name, e.g. ``vntyper-references-refs-v1-ucsc-hg19.tar.gz``.
        sha256: The trust anchor committed in install_references_config.json.
        ref_id: Physical reference id (``hg19``, ...) for a genome asset, or None for
            the common asset, which covers several `common_references` entries at once
            and carries no single `kind`/`index_command` of its own.
        entry: That reference's section of the install config, or None for the common
            asset. Carries `kind`, `installed_path` and `index_command` - the fields the
            bwa-version check below needs.
    """

    name: str
    sha256: str
    ref_id: str | None
    entry: dict[str, Any] | None


def _local_bwa_version() -> str | None:
    """Read the locally installed bwa's version banner.

    Returns:
        str | None: The version string (e.g. ``0.7.18-r1243-dirty``), or None if `bwa`
        is not on PATH or its banner does not match the expected format. Either way, the
        caller treats "unknown" the same as "matches" - it re-indexes only on a positive,
        confirmed mismatch, never on the absence of information.
    """
    try:
        completed = subprocess.run(["bwa"], capture_output=True, text=True, check=False)
    except OSError:
        return None
    match = _BWA_VERSION_RE.search(f"{completed.stdout}\n{completed.stderr}")
    return match.group(1) if match else None


def _plan_bundle_assets(install_config: dict[str, Any], references: list[str]) -> list[_BundleAsset]:
    """Decide which assets one `install_from_bundle` run needs.

    Every selected genome's asset, plus the common asset unconditionally: the MUC1
    seeds and the adVNTR databases it carries are not selectable by `--references`, the
    same rule `install_from_source` already applies to `vntyper_references`.

    Args:
        install_config: The parsed install_references_config.json.
        references: Physical reference ids this run was asked for.

    Returns:
        list[_BundleAsset]: Genome assets first (in `GENOME_SECTIONS` order), then the
        common asset last.

    Raises:
        ValueError: If the config's `bundle` section, or a selected genome entry, is
            missing the asset name or digest this function needs to plan a fetch.
    """
    bundle = install_config.get("bundle", {})
    for field in ("repository", "release_tag", "common_asset", "common_asset_sha256"):
        if not bundle.get(field):
            message = (
                f"install_references_config.json's 'bundle' section declares no {field}; "
                "cannot fetch the published reference release"
            )
            logger.error(message)
            raise ValueError(message)

    selected = set(references)
    assets: list[_BundleAsset] = []
    for section in GENOME_SECTIONS:
        for ref_id, entry in install_config.get(section, {}).items():
            if ref_id not in selected:
                continue
            asset = entry.get("asset")
            digest = entry.get("asset_sha256")
            if not asset or not digest:
                message = (
                    f"reference '{ref_id}' declares no bundle asset/asset_sha256; "
                    "cannot install it from the published release"
                )
                logger.error(message)
                raise ValueError(message)
            assets.append(_BundleAsset(asset, digest, ref_id, entry))

    assets.append(_BundleAsset(bundle["common_asset"], bundle["common_asset_sha256"], None, None))
    return assets


def _read_bundle_json(path: Path, asset_name: str) -> dict[str, Any]:
    """Parse one of an extracted asset's two metadata members.

    Args:
        path: Path the member was extracted to.
        asset_name: The asset it came from, for the error message.

    Returns:
        dict[str, Any]: The parsed document.

    Raises:
        ValueError: If the file is absent (extraction did not produce it) or is not
            valid JSON.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        message = f"{asset_name}: {path.name} was not found after extraction: {error}"
        logger.error(message)
        raise ValueError(message) from error
    try:
        return json.loads(text)
    except json.JSONDecodeError as error:
        message = f"{asset_name}: {path.name} is not valid JSON: {error}"
        logger.error(message)
        raise ValueError(message) from error


def _verify_manifest_files(staging: Path, manifest: dict[str, Any], asset_name: str) -> None:
    """Re-check every file the asset's own manifest says it extracted.

    The asset's whole-archive digest (`verify_sha256` before `safe_extract`) is the
    security-critical check - a tampered archive never reaches extraction at all. This
    is defence in depth on top of it: independent per-file digests, taken from metadata
    that is itself inside the already-verified archive, catch an extraction that silently
    produced the wrong bytes for one member without the archive as a whole disagreeing.

    Args:
        staging: Directory the asset was just extracted into.
        manifest: Parsed `release-manifest.json` from that asset.
        asset_name: The asset's file name, for the error message.

    Raises:
        ValueError: If a file the manifest names was not extracted, or does not match
            the digest the manifest records for it.
    """
    for file_entry in manifest.get("files", []):
        relative = file_entry["path"]
        extracted = staging / relative
        if not extracted.is_file():
            message = f"{asset_name}: release-manifest.json names '{relative}', but extraction did not produce it"
            logger.error(message)
            raise ValueError(message)
        verify_sha256(extracted, file_entry["sha256"])


def _reindex_if_bwa_version_differs(
    staging: Path, entry: dict[str, Any], build_info: dict[str, Any], ref_id: str | None
) -> None:
    """Re-run bwa locally when the bundled index was built with a different bwa.

    A BWA index is not guaranteed compatible across bwa versions. Comparing
    `BUILD_INFO.json`'s `bwa_version` against what is actually on this machine's PATH
    catches the mismatch before it becomes a silent alignment problem; re-indexing
    locally is cheap next to a 250 MB chromosome download, and `bwa` is already a
    required dependency of every VNtyper install.

    Only called for an entry whose `kind` is `"bwa"` - the common asset carries no BWA
    index and has no `entry` at all. This is the one place today that reads `kind`: it
    tells this function "this asset's manifest describes a bwa-indexed genome, so a
    version mismatch here is actionable", as opposed to a hypothetical future asset kind
    this check would not know how to re-index.

    Args:
        staging: Directory the genome asset was extracted into.
        entry: That reference's section of the install config - `installed_path` and
            `index_command`.
        build_info: Parsed `BUILD_INFO.json` from the same asset.
        ref_id: Physical reference id, for the log message.

    Raises:
        ValueError: If a mismatch is detected but the entry names no `index_command` to
            recover with.
        RuntimeError: If the re-index command itself fails.
    """
    bundled_version = build_info.get("bwa_version")
    local_version = _local_bwa_version()
    if not bundled_version or not local_version or bundled_version == local_version:
        return

    fasta = staging / entry["installed_path"]
    logger.warning(
        f"{ref_id}: the bundled BWA index was built with bwa {bundled_version}, but the local bwa is "
        f"{local_version}; re-indexing {fasta.name} locally"
    )

    index_command_template = entry.get("index_command")
    if not index_command_template:
        message = f"{ref_id}: no index_command configured; cannot re-index {fasta.name} after a bwa version mismatch"
        logger.error(message)
        raise ValueError(message)

    # `_quote` + `shlex.split` rather than `str()` + `str.split`: an unquoted path
    # containing whitespace (e.g. `--output-dir "/data/my refs"`) would otherwise be
    # shredded across two argv elements and abort an otherwise valid install. No shell
    # is involved here, so this is a correctness fix, not an injection one.
    command = index_command_template.format(path=_quote(fasta))
    completed = subprocess.run(shlex.split(command), capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        message = f"{ref_id}: re-indexing {fasta.name} failed: {completed.stderr.strip()}"
        logger.error(message)
        raise RuntimeError(message)
    logger.info(f"  ✓ re-indexed {fasta.name} with the local bwa {local_version}")


def _record_bundle_provenance(staging: Path, manifest: dict[str, Any], asset: _BundleAsset, release_tag: str) -> None:
    """Record install provenance for the file(s) `canonical_reference_keys` looks up.

    Only the FASTA itself for a `kind == "bwa"` genome asset - never a BWA index
    sidecar. `_reindex_if_bwa_version_differs` may rewrite the sidecars with a
    locally-built index moments after this runs, which would make a sidecar's recorded
    manifest digest describe bytes the tree no longer has; `canonical_reference_keys`
    never looks a sidecar up by its own key anyway, only the FASTA's `installed_path`
    is one. The common asset carries no BWA index and is never reindexed, so every file
    its manifest lists is recorded.

    Args:
        staging: Directory the asset was extracted into.
        manifest: Parsed `release-manifest.json`, already verified file-by-file by
            `_verify_manifest_files`.
        asset: The asset just installed.
        release_tag: The release tag the asset was fetched from.
    """
    from vntyper.scripts.reference_provenance import build_record, merge

    files_by_path = {entry["path"]: entry["sha256"] for entry in manifest.get("files", [])}
    records: dict[str, dict[str, Any]] = {}

    if asset.entry is not None and asset.entry.get("kind") == "bwa":
        relative = asset.entry["installed_path"]
        digest = files_by_path.get(relative)
        if digest is not None:
            records[relative] = build_record(sha256=digest, source="bundle", asset=asset.name, release_tag=release_tag)
    else:
        for relative, digest in files_by_path.items():
            records[relative] = build_record(sha256=digest, source="bundle", asset=asset.name, release_tag=release_tag)

    merge(staging, records)


def _install_bundle_asset(
    repository: str, release_tag: str, asset: _BundleAsset, download_dir: Path, staging: Path
) -> None:
    """Fetch, verify and extract one asset, then react to its metadata.

    Args:
        repository: `owner/repo` the release lives in.
        release_tag: The immutable release tag, e.g. ``refs-v1``.
        asset: The asset to install.
        download_dir: Scratch directory for the downloaded archive; not part of the
            installed tree.
        staging: The `staged_install` staging directory to extract into.

    Raises:
        RuntimeError: If the download fails.
        ValueError: If the digest does not match, if `safe_extract` refuses the archive,
            or if the asset's own manifest disagrees with what was extracted.
    """
    url = f"https://github.com/{repository}/releases/download/{release_tag}/{asset.name}"
    archive = download_dir / asset.name
    logger.info(f"Fetching {asset.name} from {repository}@{release_tag}")
    try:
        download_file(url, archive)
    except RuntimeError as failure:
        message = f"failed to download {asset.name} from repository {repository}, release {release_tag} ({url})"
        logger.error(message)
        raise RuntimeError(message) from failure

    # Carry-forward from the Task 3 review: the digest is checked before a single byte
    # of the archive is written into the tree.
    verify_sha256(archive, asset.sha256)
    safe_extract(archive, staging)
    archive.unlink(missing_ok=True)

    manifest = _read_bundle_json(staging / "release-manifest.json", asset.name)
    build_info = _read_bundle_json(staging / "BUILD_INFO.json", asset.name)
    try:
        _verify_manifest_files(staging, manifest, asset.name)
        # Provenance is recorded from the now-verified manifest before any reindex, not
        # after: a bwa version mismatch below rewrites the sidecar index files with a
        # locally-built index, and this only ever records the FASTA itself for a `kind
        # == "bwa"` asset (see `_record_bundle_provenance`) - never a sidecar - so the
        # order does not matter for correctness, but recording first keeps "verified"
        # and "recorded" atomic in the log even if reindexing raises.
        _record_bundle_provenance(staging, manifest, asset, release_tag)
        if asset.entry is not None and asset.entry.get("kind") == "bwa":
            _reindex_if_bwa_version_differs(staging, asset.entry, build_info, asset.ref_id)
    finally:
        # Removed unconditionally, success or failure: metadata for a *different* asset
        # must never be mistaken for this one's, and neither file is part of the
        # reference tree the pipeline reads. See `_BUNDLE_METADATA_MEMBERS`.
        for name in _BUNDLE_METADATA_MEMBERS:
            (staging / name).unlink(missing_ok=True)


def install_from_bundle(install_config: dict[str, Any], output_dir: Path, references: list[str]) -> None:
    """Install references by fetching the published, checksummed release bundle.

    This is the default `install-references` path (see `main`); `--from-source` is the
    exception, not this. Each selected genome's asset, plus the common asset that
    carries the MUC1 seeds and the adVNTR databases, is downloaded from
    ``https://github.com/{repository}/releases/download/{release_tag}/{asset}``,
    verified against the SHA-256 committed in install_references_config.json, and
    extracted - all into one `staged_install` staging directory, so a failure partway
    through (a bad download, a digest mismatch, a manifest disagreement) leaves any
    previously installed tree untouched rather than half-replaced.

    That guarantee covers every failure `staged_install` can catch, which is every
    ordinary Python exception - it does not cover an uncatchable one. `staged_install`
    activates with two renames (the live tree aside to `.output_dir.previous.*`, then
    the finished staging tree into `output_dir`) that cannot be made atomic together
    without an OS-level directory-swap primitive this codebase does not use. A SIGKILL
    landing in the narrow window between those two renames skips its `except
    BaseException` handler entirely - there is no Python code left running to catch
    anything - and leaves **`output_dir` absent** while **both** `.output_dir.previous.*`
    (the old tree, moved aside but never restored) **and** `.output_dir.staging.*` (the
    new tree, fully built but never activated) survive on disk as orphans, with nothing
    logged. This is deliberately left as an operator-diagnosable state rather than
    auto-recovered on the next run's startup: the next full `install-references` call
    self-heals by rebuilding `output_dir` from scratch (`main`'s `output_dir.mkdir(...,
    exist_ok=True)` tolerates its absence), but the two siblings are not swept or
    reported by this codebase, and nothing else on the machine can tell which of them
    (if either) is safe to delete without an operator looking at both.

    ``release-manifest.json`` and ``BUILD_INFO.json`` are read from *inside* each
    already-verified archive, never fetched separately: a loose metadata file sitting
    beside the asset has no committed digest to check it against, so trusting it would
    reopen exactly the trust hole `asset_sha256` closes. Only after an asset's own
    digest has passed does this function read its manifest, to re-check every file it
    says it extracted, and its `BUILD_INFO.json`, to detect a bwa version mismatch and
    re-index locally rather than install a BWA index this machine's bwa may not agree
    with byte-for-byte.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Directory to install references into.
        references: Physical reference ids to install (e.g. ``["hg19", "hg38"]``).

    Raises:
        RuntimeError: If a download fails, or a bwa re-index fails.
        ValueError: If the config's `bundle` section or a selected genome entry has no
            asset/digest configured, if an asset's digest does not match, if
            `safe_extract` refuses an archive, or if an asset's own manifest disagrees
            with what was extracted from it.
    """
    assets = _plan_bundle_assets(install_config, references)
    bundle = install_config["bundle"]

    with (
        staged_install(output_dir) as staging,
        tempfile.TemporaryDirectory(prefix="vntyper-bundle-download-") as download_dir,
    ):
        for asset in assets:
            _install_bundle_asset(bundle["repository"], bundle["release_tag"], asset, Path(download_dir), staging)


def setup_logging(output_dir: Path, log_file: Path | None = None) -> InstallLogHandler:
    """Attach the install run's log file without rebuilding root logging.

    The CLI is the sole owner of root logging configuration.  If it has already run,
    its handlers and level remain untouched.  Direct module execution has no handlers,
    so :func:`attach_install_log` supplies the standalone INFO console in that case.

    Args:
        output_dir (Path): Directory where logs will be stored, when `log_file` is not
            given explicitly.
        log_file (Path, optional): Exact path to write the log to. Defaults to
            ``output_dir / "install_references.log"``. `main` passes a path *outside*
            `output_dir` (see `_staging_safe_log_path`) and moves the finished file into
            place itself once `output_dir` is no longer subject to being replaced by a
            `staged_install` activation.
    Returns:
        InstallLogHandler: Installer-owned handler for `main` to detach and close.
    """
    log_file = log_file or (output_dir / "install_references.log")
    file_handler = attach_install_log(output_dir, log_file)
    logger.info(f"Logging initialized. Logs will be saved to {log_file}")
    return file_handler


def _staging_safe_log_path(output_dir: Path) -> Path:
    """Choose a log location that a `staged_install` activation cannot lose.

    Both install paths `main` can take stage their work in a temporary directory beside
    `output_dir` and activate by renaming it into place. A log file written directly
    into `output_dir` would be caught by that rename dance in a way that loses most of
    a run's own log, regardless of which side of the swap it ends up on:

    - If `staged_install`'s `seed_from_existing` copies `output_dir` into staging (the
      default, so an install of one assembly does not erase another), the log file gets
      seeded too - but only as of whatever it contained at that exact moment, typically
      just "Logging initialized..." and "Processing references: ...". Every line logged
      afterwards keeps going to the *original* file descriptor, which is what gets
      renamed to `.previous.*` and then deleted once activation succeeds - taking the
      complete log with it. The tree left at `output_dir` ends up with only the seeded,
      truncated copy.
    - If the log were excluded from seeding instead, `output_dir` would end up with no
      log file at all after activation, for the same reason: nothing ever writes into
      the *new* directory at that path, and the complete one is still the file that gets
      deleted.

    Writing beside `output_dir` instead sidesteps both failure modes: nothing this
    function's log depends on ever needs to survive a `staged_install` rename. `main`
    moves the finished file into `output_dir / "install_references.log"` itself, in
    `_finalize_install_log`, once `output_dir` has settled.

    Args:
        output_dir: The install's final destination.

    Returns:
        Path: A sibling of `output_dir`, named after it so an operator can tell which
        run a stray file belongs to if `_finalize_install_log` is ever unable to move it
        into place.
    """
    return output_dir.parent / f".{output_dir.name}.install-references.log"


def _finalize_install_log(log_file: Path, output_dir: Path) -> None:
    """Move a run's log into its final home now that `output_dir` has settled.

    Called from a `finally` block in `main`, so this runs whether the install
    succeeded, raised, or hit an early `sys.exit` - the log is worth keeping either way,
    and a failed run is exactly when an operator most wants to read it.

    The installer-owned handler is detached and closed before this function runs.  A
    rename is therefore both cheaper than copy-then-delete and cannot leave that handler
    writing through a path which has moved.

    Args:
        log_file: Where `setup_logging` was pointed - `_staging_safe_log_path`'s return
            value.
        output_dir: The install's final destination.
    """
    if not log_file.exists():
        return
    if not output_dir.exists():
        # The run failed before anything was ever activated, or activation itself
        # failed and `staged_install` could not even restore `output_dir`. Either way
        # there is nowhere to move the log into; leave it at its sibling path and say
        # so, rather than silently dropping a log that is still complete up to the
        # point of failure.
        logger.error(f"install run did not produce {output_dir}; its log remains at {log_file}")
        return
    destination = output_dir / "install_references.log"
    try:
        log_file.rename(destination)
    except OSError as error:
        logger.error(f"could not move install log {log_file} into {destination}: {error}")


def main(
    output_dir: Path,
    config_path: Path | None = None,
    skip_indexing: bool = False,
    index_threads: int = 4,
    aligners_to_use: list[str] | None = None,
    references_to_process: list[str] | None = None,
    from_source: bool = False,
    release_spec_path: Path | None = None,
    derive_only_mode: bool = False,
):
    """
    Main function to execute the install_references process.

    Args:
        output_dir (Path): Directory where references will be installed.
        config_path (Optional[Path]): Path to the main config.json file to update.
        skip_indexing (bool): Whether to skip the indexing step. Only affects
            ``--from-source``; the default bundle path re-indexes only in reaction to a
            detected bwa version mismatch (see `install_from_bundle`), and that safety
            check is not something a speed flag should be able to turn off.
        index_threads (int): Number of threads to use for indexing. Only affects
            ``--from-source``.
        aligners_to_use (list, optional): List of specific aligners to use (overrides
            config). Only affects ``--from-source``; the bundle path always installs the
            BWA index the release already carries.
        references_to_process (list, optional): List of specific references to process (e.g., ['hg19', 'hg38']).
        from_source (bool): Build every reference from its upstream source and run the
            configured derivations, instead of fetching the published, checksummed
            release bundle. This is the path the bundle build workflow itself runs.
            Defaults to False, which now takes `install_from_bundle`'s path rather than
            the legacy per-section download this replaced: references come from a
            versioned release in ``berntpopp/vntyper-data``, not from six third-party
            hosts at install time.
        release_spec_path (Optional[Path]): Release builds only - take every source URL
            and digest from this file rather than from the shipped install config.
            Requires ``from_source``.
    """
    script_dir = Path(__file__).parent
    install_config_path = script_dir / "install_references_config.json"

    install_config = load_install_config(install_config_path)

    if derive_only_mode:
        # Nothing is downloaded and nothing is indexed: this rebuilds the three derived
        # files from genomes and seeds already on disk, verifying each against its
        # committed digest. Returns before any installer runs.
        derive_only(install_config, output_dir)
        return

    ucsc_refs = install_config.get("ucsc_references", {})
    ncbi_refs = install_config.get("ncbi_references", {})
    ensembl_refs = install_config.get("ensembl_references", {})
    vntyper_refs = install_config.get("vntyper_references", {})

    aligner_config = install_config.get("aligners", {})

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create output directory {output_dir}: {e}")
        sys.exit(1)

    # Written outside `output_dir` for the whole run and moved into place at the very
    # end, in `finally` below - see `_staging_safe_log_path`'s docstring for why a log
    # written directly into `output_dir` would lose most of its own content to the
    # `staged_install` rename dance both branches below go through.
    log_file = _staging_safe_log_path(output_dir)
    install_log_handler = setup_logging(output_dir, log_file=log_file)

    try:
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

        _install_references(
            install_config,
            output_dir,
            config_path,
            skip_indexing,
            index_threads,
            aligners_to_use,
            found_refs,
            from_source,
            release_spec_path,
            aligner_config,
        )
    finally:
        primary_failure = sys.exc_info()[0] is not None
        try:
            finish_install_log(
                install_log_handler,
                lambda: _finalize_install_log(log_file, output_dir),
            )
        except Exception as error:
            logger.error(f"Failed to finish install logging: {error}")
            if not primary_failure:
                raise


def _install_references(
    install_config: dict[str, Any],
    output_dir: Path,
    config_path: Path | None,
    skip_indexing: bool,
    index_threads: int,
    aligners_to_use: list[str] | None,
    found_refs: set[str],
    from_source: bool,
    release_spec_path: Path | None,
    aligner_config: dict[str, Any],
) -> None:
    """Run the selected install path and update config.json - the part of `main` that
    needs no `finally`, split out so the log-finalisation wrapper in `main` stays small.

    Args:
        install_config: The parsed install_references_config.json.
        output_dir: Directory where references are installed.
        config_path: Path to the main config.json to update, or None to skip.
        skip_indexing: Whether to skip indexing. Only affects ``--from-source``.
        index_threads: Threads for indexing. Only affects ``--from-source``.
        aligners_to_use: Specific aligners to use. Only affects ``--from-source``.
        found_refs: Physical reference ids this run will install.
        from_source: Build from upstream sources instead of fetching the release bundle.
        release_spec_path: Release builds only - source URLs/digests to use instead of
            the shipped config. Requires `from_source`.
        aligner_config: The install config's ``aligners`` section.
    """
    if from_source:
        # Initialize aligners. Only `--from-source` reads this: the bundle path in the
        # `else` branch below always installs the BWA index the release already carries,
        # and re-indexes on its own only when `install_from_bundle` detects a bwa version
        # mismatch - computing an "enabled aligners" set for it would be dead work, and
        # logging a whole "ALIGNER CONFIGURATION" section for a run it cannot affect
        # would tell the operator a flag mattered when it did not.
        enabled_aligners: dict[str, dict[str, Any]] = {}
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

        # The source path owns every section at once: it must verify each download
        # before expanding it, and it must know which chromosome FASTAs landed where so
        # the derivations can cut from them. It writes no md5_checksums.txt - the
        # digests it enforces are SHA-256 values committed in the install config.
        #
        # Wrapped in `staged_install` here rather than inside `install_from_source`
        # itself (carry-forward from PR-1 Codex finding M5): `install_from_source`'s own
        # direct-unit-test suite in test_install_references_derivations.py asserts its
        # returned paths equal `output_dir`-relative paths for the literal directory it
        # was called with, so redirecting its writes has to happen at this call site, by
        # substituting `staging` for `output_dir`, rather than by changing what the
        # function does with the `output_dir` it is given. A late failure - a bad
        # download, a failed derivation, a `bwa index` that runs out of disk mid-genome -
        # now leaves any previously installed tree untouched instead of a partially
        # rebuilt one, the same guarantee the bundle path gets from `install_from_bundle`
        # below - including the same inherent exception to it: a SIGKILL landing between
        # `staged_install`'s two activation renames leaves `output_dir` absent and both
        # `.previous.*`/`.staging.*` siblings orphaned with nothing logged, deliberately
        # left for manual inspection rather than auto-recovered. See
        # `install_from_bundle`'s docstring for the full explanation.
        release_spec = load_install_config(release_spec_path) if release_spec_path else None
        if release_spec is not None:
            logger.info(f"Taking source URLs and digests from {release_spec_path}")
        with staged_install(output_dir) as staging:
            install_from_source(
                install_config,
                staging,
                sorted(found_refs),
                enabled_aligners,
                index_threads,
                release_spec,
                skip_indexing,
            )
    else:
        # The default path: fetch the published, checksummed release bundle rather than
        # rebuilding from six third-party hosts. `install_from_bundle` stages its own
        # install (see its docstring), so it is handed `output_dir` directly rather than
        # a staging directory computed here - unlike the `--from-source` branch above,
        # its own unit tests call it exactly this way.
        install_from_bundle(install_config, output_dir, sorted(found_refs))

    # Update the main config.json with new reference paths if config_path is provided.
    # `canonical_reference_keys` re-derives every key from `install_config` (the full,
    # unfiltered config loaded above) and `output_dir` directly, rather than from the
    # per-section dicts this run happened to process - so the keys it writes are exactly
    # the ones `vntyper/` reads, and a partial run never writes a key pointing at a file
    # that was never installed.
    if config_path and config_path.exists():
        update_config(config_path, canonical_reference_keys(install_config, output_dir))
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
