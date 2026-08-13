"""Build the command line for the vendored Kestrel genotyper."""

from __future__ import annotations

import logging
import shlex
from collections.abc import Sequence
from pathlib import Path

from vntyper.scripts.command_builders import quote_path

logger = logging.getLogger(__name__)

#: Kestrel options that provably cannot reach the k-mer counter and do not collide with
#: a value this builder already sets, mapped to how many values each takes.
#:
#: This is an **allowlist**, and it is enforced only when Kestrel is handed a pre-built
#: IKC. Once counting runs as its own KAnalyze step, an option that changes counting has
#: to be applied to both commands or to neither, and a deny-list is fail-open at exactly
#: the boundary where a mistake silently changes genotypes. Every name here was checked
#: against ``java -jar kestrel.jar -h`` for the pinned 1.0.1 build.
#:
#: Four classes are deliberately absent, and a rejected option is a loud configuration
#: error the operator can see, while an accepted count-affecting one is a silent wrong
#: answer:
#:
#: * **Counter-affecting.** ``--mincount``, ``--minsize``, ``--minmask``, ``--charset``,
#:   ``--seqfilter``, ``--quality``, ``--noseqfilter``, ``-k``/``--ksize``,
#:   ``--temploc``, ``--free``/``--nofree`` (which reaches ``setFreeSegment``), and
#:   ``--lib``/``--liburl`` (``getCountModule`` forwards library URLs into the counter).
#: * **Count-map selection.** ``--memcount``/``--nomemcount`` choose an in-memory count
#:   map instead of the IKC one, so a supplied IKC would be read as sequence.
#: * **IKC lifecycle.** ``--rmikc``/``--normikc``, which this module owns.
#: * **Values this builder sets.** ``-f``/``--format``, ``-o``/``--out``,
#:   ``-p``/``--hapout``, ``-m``/``--outfmt``, ``--hapfmt``, ``-r``/``--ref``,
#:   ``-s``/``--sample``, ``--filespersample``, and the logging options.
#:
#: ``--noautoflank`` and ``--maxrepeat`` are **not** Kestrel options at all, despite
#: looking like the negations and siblings of ones that are.
CALL_ONLY_OPTIONS: dict[str, int] = {
    "--alpha": 1,
    "--ambiregions": 0,
    "--ambivar": 0,
    "--anchorboth": 0,
    "--autoflank": 0,
    "--byreference": 0,
    "--byregion": 0,
    "--countrev": 0,
    "--decaymin": 1,
    "--diffq": 1,
    "--flank": 1,
    "--interval": 1,
    "--maxalignstates": 1,
    "--maxhapstates": 1,
    "--mindiff": 1,
    "--noambigregions": 0,
    "--noambivar": 0,
    "--noanchorboth": 0,
    "--nocountrev": 0,
    "--norevregion": 0,
    "--normrefdesc": 0,
    "--peakscan": 1,
    "--revregion": 0,
    "--rmrefdesc": 0,
    "--scanlimitfactor": 1,
    "--varfilter": 1,
    "--weight": 1,
}

#: Allowlisted options whose value may legitimately begin with ``-``.
#:
#: The value check exists so ``--ambivar --memcount`` cannot smuggle ``--memcount``
#: through as ``--ambivar``'s value. These options take signed or negative-leading
#: numbers, so for them an option-looking value is accepted: ``--weight`` takes a
#: comma-separated list whose first element is routinely negative, and the decay and
#: quantile settings are floats an operator may write as ``-0.5`` by mistake, which
#: Kestrel should reject rather than this builder.
NUMERIC_VALUE_OPTIONS: frozenset[str] = frozenset(
    {"--alpha", "--decaymin", "--diffq", "--mindiff", "--scanlimitfactor", "--weight"}
)


def reject_non_call_options(tokens: list[str]) -> None:
    """Reject anything in ``additional_settings`` that could reach the k-mer counter.

    Applied only when Kestrel is being handed a pre-built IKC. Without the split there
    is no second command to desynchronise from, and enforcing it anyway would reject
    configurations that are perfectly valid for stock Kestrel.

    Args:
        tokens: The ``shlex``-split ``additional_settings``.

    Raises:
        ValueError: If any token is not an allowlisted call-only option or the value of
            one.
    """
    expecting: str | None = None
    for token in tokens:
        if expecting is not None:
            # A value is required here. Reject an option-looking token unless this
            # option's value may legitimately start with '-': without this,
            # "--ambivar --memcount" swallows --memcount as a value and the guard is
            # bypassed entirely.
            if token.startswith("-") and expecting not in NUMERIC_VALUE_OPTIONS:
                msg = f"additional_settings: {expecting!r} requires a value but was followed by {token!r}."
                logger.error(msg)
                raise ValueError(msg)
            expecting = None
            continue
        if token == "--":
            msg = (
                "additional_settings cannot contain '--': everything after it is a bare input "
                "operand, and Kestrel requires a supplied IKC to be its sole input."
            )
            logger.error(msg)
            raise ValueError(msg)
        if not token.startswith("-"):
            msg = (
                f"additional_settings cannot contain the bare operand {token!r}: Kestrel is "
                "given a pre-built IKC, which must be its sole input."
            )
            logger.error(msg)
            raise ValueError(msg)
        name, separator, _attached = token.partition("=")
        arity = CALL_ONLY_OPTIONS.get(name)
        if arity is None:
            msg = (
                f"additional_settings may only carry Kestrel call options; {name!r} is not one. "
                "An option that reaches the k-mer counter would apply to Kestrel but not to the "
                "separate KAnalyze count step, silently changing genotypes."
            )
            logger.error(msg)
            raise ValueError(msg)
        if arity == 0 and separator:
            msg = f"additional_settings: {name!r} is a flag and takes no value."
            logger.error(msg)
            raise ValueError(msg)
        expecting = None if (arity == 0 or separator) else name
    if expecting is not None:
        msg = f"additional_settings: {expecting!r} requires a value but none was given."
        logger.error(msg)
        raise ValueError(msg)


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
    java_opts: str = "",
    ikc_path: str | Path | None = None,
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
            Defaults to an empty string. When ``ikc_path`` is given these are
            restricted to :data:`CALL_ONLY_OPTIONS`.
        java_opts (str, optional): Extra JVM options, rendered between ``-Xmx`` and
            ``-jar``. The call step is single-threaded, so ``-XX:+UseSerialGC`` here
            avoids G1 spawning a GC worker per core against one application thread.
        ikc_path (str | Path | None, optional): A pre-built indexed k-mer count file.
            When given, Kestrel is handed ``-f ikc <path>`` instead of the FASTQ
            operands and does not count anything itself:
            ``IkcCountMap.preModuleRun`` adopts the file and returns false. The FASTQ
            paths are still validated, because they are what the count step consumed.

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

    # The allowlist applies only to the split path. Without a separate count step there
    # is no second command for a count-affecting option to desynchronise from, so
    # enforcing it unconditionally would reject configurations that are valid for stock
    # Kestrel -- which is why the guard and the split are inseparable rather than merely
    # committed together.
    if ikc_path is not None:
        reject_non_call_options(additional_tokens)

    if ikc_path is None:
        sample_inputs = " ".join(quote_path(path) for path in normalized_fastqs)
    else:
        # `-f ikc` sets the format for the files that follow it. Kestrel requires a
        # supplied IKC to be the sample's sole source: `IkcCountMap.preModuleRun` adopts
        # it only when `sources.length == 1`.
        sample_inputs = f"-f ikc {quote_path(ikc_path)}"
    jvm_options = f" {java_opts}" if java_opts else ""

    # `run_command` runs this as one string under bash (trap 9), so quoting can only
    # happen here. Paths and the sample name are quoted -- the web service accepts
    # filenames that reach this line, and a space alone splits one argument into two.
    # `java_path` is left raw: config.json holds a command *prefix* there, which
    # quoting would collapse into a single token bash looks for as one binary. The
    # numeric and memory settings are operator-controlled config, not user input, and
    # are interpolated as before.
    #
    # `--logstderr` is gone, and it never did anything: `OptLogStderr` and `OptLogStdout`
    # call the same setter on a single `logFile` field, so last-one-wins and
    # `--logstdout` -- emitted second -- has always won. `run_command` merges stderr into
    # stdout regardless. Removing it is a clarity fix with no behavioural change.
    base_command = (
        f"{java_path} -Xmx{java_memory}{jvm_options} -jar {quote_path(kestrel_path)} -k {kmer_size} "
        f"--maxalignstates {max_align_states} --maxhapstates {max_hap_states} "
        f"-r {quote_path(reference_vntr)} -o {quote_path(vcf_out)} "
        f"-s{quote_path(sample_name)} {sample_inputs} "
        f"--hapfmt sam -p {quote_path(f'{output_dir}/output.sam')} --logstdout "
        f"--loglevel {log_level.upper()} --temploc {quote_path(output_dir)}"
    )
    if additional_settings:
        base_command += " " + additional_settings
    return base_command
