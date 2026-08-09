"""Build the command line for the vendored Kestrel genotyper."""

import logging

from vntyper.scripts.command_builders import quote_path

logger = logging.getLogger(__name__)


def construct_kestrel_command(
    kmer_size,
    kestrel_path,
    reference_vntr,
    output_dir,
    fastq_1,
    fastq_2,
    vcf_out,
    java_path,
    java_memory,
    max_align_states,
    max_hap_states,
    log_level,
    sample_name,
    additional_settings="",
):
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
        fastq_1 (str): Path to the first FASTQ (R1).
        fastq_2 (str | None): Optional path to the second FASTQ (R2).
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
        ValueError: If fastq_1 is missing.

    Returns:
        str: The constructed command line string to run Kestrel.
    """
    if not fastq_1:
        raise ValueError("FASTQ input files are missing or invalid.")

    fastq_inputs = quote_path(fastq_1)
    if fastq_2:
        fastq_inputs += f" {quote_path(fastq_2)}"

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
