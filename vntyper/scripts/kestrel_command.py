"""Build the command line for the vendored Kestrel genotyper."""

from __future__ import annotations

import logging
import shlex
from collections.abc import Sequence
from pathlib import Path

from vntyper.scripts.command_builders import quote_path

logger = logging.getLogger(__name__)


def construct_kestrel_command(
    *,
    kmer_size: int,
    kestrel_path: str | Path,
    reference_vntr: str | Path,
    output_dir: str | Path,
    fastq_files: Sequence[str | Path],
    vcf_out: str | Path,
    java_path: str,
    java_memory: str,
    max_align_states: int,
    max_hap_states: int,
    log_level: str,
    sample_name: str,
    additional_settings: str = "",
) -> str:
    """
    Constructs the command for running Kestrel's mapping-free genotyping.

    This primarily sets Kestrel parameters:
      - k-mer size
      - maximum alignment/haplotype states
      - input FASTQs + reference MUC1-VNTR
      - output VCF + intermediate SAM
      - logging level, memory allocation

    Args:
        kmer_size (int): Size of the kmer to use in Kestrel alignment.
        kestrel_path (str): Path to the Kestrel jar file.
        reference_vntr (str): Path to the reference VNTR file (FASTA).
        output_dir (str): Directory for the intermediate & final output files.
        fastq_files: Non-empty ordered FASTQ paths for one Kestrel sample.
        vcf_out (str): Path to the VCF output file from Kestrel.
        java_path (str): Path to the Java executable.
        java_memory (str): Memory allocated to the JVM (e.g., "12g").
        max_align_states (int): Kestrel param for alignment states.
        max_hap_states (int): Kestrel param for haplotype states.
        log_level (str): Logging level (DEBUG, INFO, etc.).
        sample_name (str): Sample name for labeling in the VCF.
        additional_settings (str, optional): Extra command-line options to append.
            Defaults to an empty string.

    Raises:
        ValueError: If FASTQ operands or protected additional settings are invalid.

    Returns:
        str: The constructed command line string to run Kestrel.
    """
    if isinstance(fastq_files, (str, Path)) or not isinstance(fastq_files, Sequence):
        raise ValueError("fastq_files must be a non-scalar sequence of FASTQ paths.")
    normalized_fastqs: list[str] = []
    for fastq in fastq_files:
        if not isinstance(fastq, (str, Path)) or not str(fastq):
            raise ValueError("FASTQ input files are missing or invalid.")
        fastq_text = str(fastq)
        if Path(fastq_text).parent == Path(".") and fastq_text.startswith("-"):
            raise ValueError(f"Relative FASTQ input looks like an option: {fastq_text}")
        normalized_fastqs.append(fastq_text)
    if not normalized_fastqs:
        raise ValueError("FASTQ input files are missing or invalid.")
    if len(set(normalized_fastqs)) != len(normalized_fastqs):
        raise ValueError("FASTQ input files contain duplicate paths.")

    if not isinstance(additional_settings, str):
        raise ValueError("additional_settings must be a string.")
    try:
        additional_tokens = shlex.split(additional_settings)
    except ValueError as exc:
        raise ValueError(f"additional_settings could not be parsed: {exc}") from exc
    protected = ("--sample", "--filespersample")
    if any(
        token == "-s"
        or (token.startswith("-s") and not token.startswith("--"))
        or token in protected
        or any(token.startswith(f"{option}=") for option in protected)
        for token in additional_tokens
    ):
        raise ValueError("additional_settings cannot override Kestrel sample or input grouping.")

    fastq_inputs = " ".join(quote_path(path) for path in normalized_fastqs)

    # `run_command` runs this as one string under bash (trap 9), so quoting can only
    # happen here. Paths and the sample name are quoted -- the web service accepts
    # filenames that reach this line, and a space alone splits one argument into two.
    # `java_path` is left raw: config.json holds a command *prefix* there, which
    # quoting would collapse into a single token bash looks for as one binary. The
    # numeric and memory settings are operator-controlled config, not user input, and
    # are interpolated as before.
    base_command = (
        f"{java_path} -Xmx{java_memory} -jar {quote_path(kestrel_path)} -k {kmer_size} "
        f"--maxalignstates {max_align_states} --maxhapstates {max_hap_states} "
        f"-r {quote_path(reference_vntr)} -o {quote_path(vcf_out)} "
        f"-s{quote_path(sample_name)} {fastq_inputs} "
        f"--hapfmt sam -p {quote_path(f'{output_dir}/output.sam')} --logstderr --logstdout "
        f"--loglevel {log_level.upper()} --temploc {quote_path(output_dir)}"
    )
    if additional_settings:
        base_command += " " + additional_settings
    return base_command
