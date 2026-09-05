# vntyper/modules/shark/shark_filtering.py

import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.command_builders import quote_path
from vntyper.scripts.pipeline_read_routing import count_fastq_records
from vntyper.scripts.run_configuration import resolve_compatibility_component
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


def select_muc1_region_fasta(config: dict, main_config: dict, reference_assembly: str) -> str:
    """Pick the MUC1 region FASTA for an assembly, most authoritative source first.

    1. ``main_config["reference_data"]`` - what `install-references` writes, so a
       custom ``--output-dir`` is honoured. `shark_config.json` is a separate file
       that `--config-path` never touches, so without this layer an installed tree
       would be ignored.
    2. ``config["shark_settings"]`` keyed by assembly - the shipped default.
    3. the legacy flat ``muc1_region_fasta`` key, **only** when the config carries no
       keyed entry at all. A partly-populated keyed config is incomplete, not legacy.

    Resolution is by key *membership*: a key present with value ``None`` is a
    deliberate "disabled" and does not fall through.

    Args:
        config: The shark_config dictionary.
        main_config: The main configuration dictionary.
        reference_assembly: Supported assembly label.

    Returns:
        str: Path to the region FASTA.

    Raises:
        ValueError: If nothing resolves, naming the assembly and every key tried.
    """
    from vntyper.scripts.reference_registry import reference_keys
    from vntyper.scripts.reference_resolution import resolve_from_mapping

    (key,) = reference_keys("shark", reference_assembly)
    settings = config.get("shark_settings", {})

    resolved = resolve_from_mapping("shark", reference_assembly, main_config.get("reference_data", {}))
    if resolved is not None:
        if resolved.value:
            return resolved.value
        raise ValueError(f"reference_data[{resolved.key!r}] is null; SHARK is disabled for {reference_assembly!r}")

    if key in settings:
        if settings[key]:
            return settings[key]
        raise ValueError(f"shark_settings[{key!r}] is null; SHARK is disabled for {reference_assembly!r}")

    keyed = [name for name in settings if name.startswith("muc1_region_fasta_")]
    if not keyed:
        legacy_key = "muc1_region_fasta"
        if legacy_key in settings:
            if settings[legacy_key]:
                return settings[legacy_key]
            raise ValueError(f"shark_settings[{legacy_key!r}] is null; SHARK is disabled for {reference_assembly!r}")

    raise ValueError(
        f"No SHARK MUC1 region FASTA for reference_assembly {reference_assembly!r}. "
        f"Tried reference_data[{key!r}] and shark_settings[{key!r}]."
    )


@dataclass(frozen=True)
class SharkSearchParameters:
    """Runtime search parameters for SHARK filtering (#312)."""

    kmer_size: int | str = 17
    confidence: float | str = 0.6

    @property
    def k_str(self) -> str:
        """String representation of k-mer size."""
        return str(self.kmer_size)

    @property
    def c_str(self) -> str:
        """String representation of confidence threshold."""
        return str(self.confidence)


def shark_search_parameters(config: dict | None) -> SharkSearchParameters:
    """Resolve the effective SHARK search parameters from the runtime sidecar.

    Args:
        config: The shark_config or runtime sidecar dictionary.

    Returns:
        Resolved SharkSearchParameters with defaults k=17, c=0.6.
    """
    settings = (config or {}).get("shark_settings", {})
    raw_k = settings.get("kmer_size", 17)
    raw_c = settings.get("confidence", 0.6)
    return SharkSearchParameters(kmer_size=raw_k, confidence=raw_c)


def run_shark_filter(
    fastq_1,
    fastq_2,
    output_dir,
    config,  # shark_config dict
    main_config,  # main config dict (for tool paths)
    sample_name,
    reference_assembly="hg19",
    threads=4,
    *,
    resolved_component=None,
    custom_context_active=False,
):
    """
    Run SHARK filtering on FASTQ files.

    Args:
        fastq_1 (str): Path to input FASTQ R1.
        fastq_2 (str): Path to input FASTQ R2.
        output_dir (str): Directory to store SHARK output.
        config (dict): The shark_config dictionary (contains muc1_region_fasta_* keys).
        main_config (dict): The main configuration dictionary (contains tool paths and,
            when populated by `install-references`, `reference_data`).
        sample_name (str): Sample name.
        reference_assembly (str): Supported assembly label. Selects which MUC1 region
            FASTA (hg19- or hg38-coordinate) SHARK filters against -- see
            `select_muc1_region_fasta`. See issue #152.
        threads (int): Number of threads to use for SHARK.
        resolved_component: Explicit empty SHARK decision component for this run.
        custom_context_active: Whether an explicit custom profile owns this run.

    Returns:
        tuple: (filtered_fastq_1, filtered_fastq_2)

    Raises:
        ValueError: If no MUC1 region FASTA resolves for `reference_assembly` -- see
            `select_muc1_region_fasta`.
        RuntimeError: If SHARK filtering fails.
    """
    resolve_compatibility_component(
        "shark",
        resolved_component,
        custom_context_active=custom_context_active,
    )

    # shark_path should come from the main config since it contains the tool paths.
    shark_path = main_config.get("tools", {}).get("shark", "shark")

    # Get the assembly-selected muc1_region_fasta.
    muc1_region_fasta = select_muc1_region_fasta(config, main_config, reference_assembly)

    filtered_fastq_1 = os.path.join(output_dir, f"{sample_name}_shark_R1.fastq")
    filtered_fastq_2 = os.path.join(output_dir, f"{sample_name}_shark_R2.fastq")

    search_params = shark_search_parameters(config)

    # `run_command` runs this as one string under bash (trap 9), so quoting can only
    # happen here. Every operand is a path or a thread count and is quoted;
    # `shark_path` is not, because config.json holds a command *prefix* there --
    # "mamba run -n shark_env shark" (trap 6) -- which quoting would collapse into a
    # single token bash looks for as one binary.
    command = (
        f"{shark_path} -r {quote_path(muc1_region_fasta)} "
        f"-1 {quote_path(fastq_1)} -2 {quote_path(fastq_2)} "
        f"-o {quote_path(filtered_fastq_1)} -p {quote_path(filtered_fastq_2)} "
        f"-t {quote_path(threads)} "
        f"-k {quote_path(search_params.k_str)} "
        f"-c {quote_path(search_params.c_str)}"
    )

    log_file = Path(output_dir) / f"{sample_name}_shark.log"
    logger.info(f"Running SHARK filtering with command: {command}")

    success = run_command(command, str(log_file), critical=True)
    if not success:
        logger.error("SHARK filtering failed.")
        raise RuntimeError("SHARK filtering failed.")

    return filtered_fastq_1, filtered_fastq_2


def write_shark_step_summary(
    filtered_fastq_1: str | Path,
    filtered_fastq_2: str | Path,
    output_path: str | Path,
) -> dict[str, str]:
    """Write readable metadata for the FASTQs produced by SHARK.

    The pipeline previously recorded a nonexistent ``filtered_R1.fastq.gz`` with an
    unsupported ``fastq`` summary type. This JSON sidecar instead names the two real,
    uncompressed SHARK outputs and records how many complete reads each contains.

    Args:
        filtered_fastq_1: SHARK's filtered R1 FASTQ.
        filtered_fastq_2: SHARK's filtered R2 FASTQ.
        output_path: Destination for the JSON sidecar.

    Returns:
        The string-valued payload written to ``output_path``.

    Raises:
        ValueError: If either FASTQ cannot be counted or has an incomplete record.
        OSError: If the summary cannot be written.
    """
    # The configured FASTQ shape is validated elsewhere as exactly four lines per
    # record, so this literal cannot represent a second configurable policy.
    kept_reads_r1 = count_fastq_records(filtered_fastq_1, lines_per_record=4)
    kept_reads_r2 = count_fastq_records(filtered_fastq_2, lines_per_record=4)
    payload = {
        "filtered_fastq_1": str(filtered_fastq_1),
        "filtered_fastq_2": str(filtered_fastq_2),
        "kept_reads_r1": str(kept_reads_r1),
        "kept_reads_r2": str(kept_reads_r2),
    }
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=4)
    logger.info(
        "SHARK kept %s R1 / %s R2 reads; recorded in %s",
        kept_reads_r1,
        kept_reads_r2,
        output_path,
    )
    return payload
