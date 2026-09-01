# vntyper/modules/advntr/advntr_genotyping.py

import logging
import os
import re
import shlex
from pathlib import Path

import numpy as np
import pandas as pd

from vntyper.modules.advntr.advntr_result_io import invalidate_advntr_artifact, publish_advntr_result
from vntyper.modules.advntr.advntr_variant_annotations import (
    DELETION_PATTERN,
    INSERTION_PATTERN,
    derive_ru_and_pos,
)
from vntyper.modules.advntr.artifact_evidence import EVIDENCE_DISPOSITION_COLUMN
from vntyper.scripts.command_builders import quote_path
from vntyper.scripts.flagging import ADVNTR_FLAG_COLUMNS, compile_flag_rules
from vntyper.scripts.nomenclature_annotate import NOMENCLATURE_COLUMNS, annotate_advntr_frame
from vntyper.scripts.utils import load_config, run_command

logger = logging.getLogger(__name__)

#: Every ``LEN<n>`` token in an adVNTR ``State`` string.
#:
#: #192 defines the inserted length of a compound state as the **sum** of these, not the
#: first and not zero. @hassansaei, 2026-08-06: "use the sum of inserted lengths [...]
#: ``I9_2_A_LEN9&I50_2_A_LEN3`` to Insertion_len = 9 + 3 = 12 (not first-LEN-only). [...]
#: Do not keep first-LEN-wins as the defined semantics."
#:
#: What this replaced, recorded because #192 and the #198 PR body both diagnose it wrongly
#: as "first-LEN-wins": the historic expression was ``(LEN.*)``, greedy to the end of the
#: string, whose remainder was then split on ``LEN`` and coerced with
#: ``pd.to_numeric(errors="coerce")``. A single terminal ``LEN`` left a bare number and
#: parsed (``I22_2_G_LEN1`` -> 1), but *any* material after the first ``LEN`` -- a second
#: ``LEN``, a further ``&`` part, even a trailing ``&`` -- left a non-numeric remainder
#: that became **zero**. ``I9_2_A_LEN9&I50_2_A_LEN3`` gave 0, not 9. No input ever
#: behaved as "first-LEN-wins".
#:
#: ``Insertion_len`` feeds ``Net_indel_length = Insertion_len - Deletion_length``, and a
#: net change of 0 is not the pathogenic frame, so a collapsed pure-insertion compound was
#: dropped in silence. Summation makes those calls visible; it also drops a net-in-frame
#: compound such as ``I9_2_A_LEN2&D50_2&D51_2`` that the zero had been keeping. Both
#: directions are pinned by
#: ``tests/unit/test_advntr_output_parsing.py::TestSummedInsertionLength``.
LEN_TOKEN_PATTERN = re.compile(r"LEN(\d+)")


def sum_insertion_lengths(variant: str) -> int:
    """Total inserted bases named by every ``LEN`` token in an adVNTR ``State`` string.

    Implements the #192 decision: a compound state's inserted length is the sum over all
    of its parts, so ``I9_2_A_LEN9&I50_2_A_LEN3`` is 12. A state naming no ``LEN`` token,
    or naming one that is malformed (``LENX``), contributes 0.

    Args:
        variant (str): An adVNTR ``State``/``Variant`` string, possibly compound
            (parts joined by ``&``).

    Returns:
        int: The sum of every ``LEN<n>``; 0 when there is none. A non-string input
            (a ``NaN`` from an empty ``State`` column) is 0, as the historic
            ``str.extract`` pipeline it replaces also was.
    """
    if not isinstance(variant, str):
        return 0
    return sum(int(match) for match in LEN_TOKEN_PATTERN.findall(variant))


# -------------------------------------------------------------------------
def load_advntr_config(config_path=None):
    """
    Loads the adVNTR configuration file.
    """
    if config_path is None:
        # Default path to advntr_config.json
        config_path = os.path.join(os.path.dirname(__file__), "advntr_config.json")
    return load_config(config_path)


# Load the adVNTR settings
advntr_config = load_advntr_config()
advntr_settings = advntr_config.get("advntr_settings", {})


def advntr_output_extension(settings: dict) -> str:
    """
    Return the file extension adVNTR writes its results to, for a settings mapping.

    Args:
        settings (dict): An ``advntr_settings`` mapping.

    Returns:
        str: ``".vcf"`` when ``output_format`` is ``"vcf"``, otherwise ``".tsv"``.

    Raises:
        KeyError: If ``output_format`` is absent. The shipped ``advntr_config.json`` sets it;
            a partial configuration is not supported input. See :func:`run_advntr`.

    Note:
        The mapping is a **parameter, not the module global**, deliberately.
        :func:`run_advntr` reads the import-time :data:`advntr_settings`, while
        ``pipeline.py`` calls :func:`load_advntr_config` a second time and derives its own
        local mapping. Those are two independently loaded states, and this extension is
        derived by both -- the producer of the output path and its consumer.

        Taking the mapping as a parameter does **not** merge those two states, and this
        function should not be read as claiming it does. What it removes is the *duplicated
        derivation*: the rule turning ``output_format`` into an extension existed twice, once
        here and once inlined in ``pipeline.py``, so the two could drift apart in a way no
        test would catch. Now there is one rule, applied to whichever mapping the caller
        actually used, and each caller's mapping is visible at its call site. Both loads read
        the same packaged file, so they agree in practice; unifying the loads themselves would
        change which configuration ``pipeline.py`` reads and is a separate change (#247).
    """
    return ".vcf" if settings["output_format"] == "vcf" else ".tsv"


def resolve_advntr_threads(settings, pipeline_threads):
    """Return the ``-t`` value for a settings mapping and the pipeline's ``--threads``.

    ``null`` means *inherit the pipeline value*; an integer overrides it. Both are
    validated here rather than being interpolated into a command string and left for
    adVNTR to reject, because adVNTR's own message (``threads cannot be less than 1``,
    ``advntr_commands.py:72-73``) names neither the key nor the file the bad value came
    from.

    Args:
        settings (dict): An ``advntr_settings`` mapping.
        pipeline_threads (int): The pipeline's resolved ``--threads``.

    Returns:
        int: A positive thread count.

    Raises:
        KeyError: If ``threads`` is absent. A partial mapping is not supported input --
            the same rule ``output_format`` follows -- and ``null`` is not the same as
            missing: ``null`` is an explicit request to inherit.
        ValueError: If either value is not a positive integer.
    """
    if "threads" not in settings:
        raise KeyError(
            "advntr_settings is missing 'threads'. A partial mapping is not supported "
            "input: set an integer to pin the thread count, or null to inherit the "
            "pipeline's --threads. See vntyper/modules/advntr/advntr_config.json."
        )

    configured = settings["threads"]
    source = "advntr_settings['threads']"
    if configured is None:
        configured = pipeline_threads
        source = "the pipeline's --threads (advntr_settings['threads'] is null)"

    # bool is an int subclass, and True would silently become one thread.
    if isinstance(configured, bool) or not isinstance(configured, int):
        raise ValueError(f"{source} must be an integer, got {configured!r} ({type(configured).__name__})")
    if configured < 1:
        raise ValueError(f"{source} must be at least 1, got {configured}")
    return configured


#: Every spelling adVNTR declares for an option ``run_advntr`` sets itself
#: (``advntr/__main__.py`` at the pinned revision: ``-a/--alignment_file``, ``-o/--outfile``,
#: ``-fs/--frameshift``, ``-m/--models``, ``-t/--threads``, ``-vid/--vntr_id`` and
#: ``--working_directory``).
MANAGED_ADVNTR_OPTIONS = frozenset(
    {
        "-a",
        "--alignment_file",
        "-o",
        "--outfile",
        "-fs",
        "--frameshift",
        "-m",
        "--models",
        "-t",
        "--threads",
        "-vid",
        "--vntr_id",
        "--working_directory",
    }
)

#: What owns each managed option, so a refusal can say what to set instead of only what to
#: remove.
MANAGED_ADVNTR_OPTION_OWNERS = {
    "-t": "advntr_settings['threads']",
    "--threads": "advntr_settings['threads']",
    "-vid": "advntr_settings['vid']",
    "--vntr_id": "advntr_settings['vid']",
    "-o": "the output directory and sample name run_advntr is called with",
    "--outfile": "the output directory and sample name run_advntr is called with",
    "-m": "run_advntr's db_file argument",
    "--models": "run_advntr's db_file argument",
    "-a": "run_advntr's sorted_bam argument",
    "--alignment_file": "run_advntr's sorted_bam argument",
    "--working_directory": "run_advntr's output argument",
    "-fs": "the fixed `genotype -fs` mode this module runs",
    "--frameshift": "the fixed `genotype -fs` mode this module runs",
}

#: Every other option adVNTR's ``genotype`` subparser declares -- the ones this module does
#: not set, and which ``additional_commands`` therefore exists to carry.
#:
#: ``-h``/``--help`` are deliberately absent: adVNTR would print help and exit zero without
#: genotyping anything.
ADVNTR_EXTRA_OPTIONS = frozenset(
    {
        "-r",
        "--reference_filename",
        "-f",
        "--fasta",
        "-p",
        "--pacbio",
        "-n",
        "--nanopore",
        "--outfmt",
        "--vid_file",
        "--append",
        "--noref_aln",
        "--min_read_length",
        "-e",
        "--expansion",
        "-c",
        "--coverage",
        "--haploid",
        "-naive",
        "--naive",
        "-u",
        "--update",
        "--fullru",
        "-aln",
        "--aln",
    }
)

#: What argparse treats as a negative number rather than an option string.
_NEGATIVE_NUMBER = re.compile(r"^-\d+$|^-\d*\.\d+$")


def advntr_option_token(token):
    """Return the option spelling ``token`` is, or ``None`` if it is a value.

    Only the exact spelling is returned -- ``--opt=value`` yields ``--opt``, and nothing
    else is normalised. That is the point: see :func:`resolve_additional_commands` for why
    this is an allow-list rather than a pattern match.

    Args:
        token (str): One ``shlex``-split word of ``additional_commands``.

    Returns:
        str | None: The option spelling, or ``None`` if the word is a value.
    """
    if not token.startswith("-") or _NEGATIVE_NUMBER.match(token):
        return None
    return token.split("=", 1)[0] if token.startswith("--") else token


def resolve_additional_commands(settings):
    """Return ``additional_commands``, refusing any option this module already sets.

    ``additional_commands`` is interpolated verbatim *after* every option
    :func:`run_advntr` builds, and adVNTR parses the result with argparse, which lets the
    last occurrence of an option win. A fragment repeating a managed option therefore
    replaces a value VNtyper chose, silently and after validation: measured against a
    parser built from adVNTR's own declarations, each of ``--threads 3``, ``-t3``,
    ``--threads=3``, ``--thr 3`` and ``--thread 3`` overrides an earlier ``-t 12``.

    That is what made :func:`resolve_advntr_threads` decorative -- it rejects 0, negatives,
    bools and non-integers, and then a fragment saying ``--threads 0`` hands adVNTR exactly
    the value it refused. ``-o`` and ``-m`` are worse: ``pipeline.py`` reconstructs the
    output path independently and reads the file back, so a redirected artefact is a
    missing result rather than a relocated one.

    The rule is an **allow-list**, and that is not fastidiousness. The first version of this
    guard pattern-matched against the managed options -- exact spellings, ``--opt=value``,
    long abbreviations, and a one-character short option carrying an attached value. It
    therefore had to reimplement argparse's option matching, and it missed two whole
    classes, both measured against a parser built from adVNTR's own declarations:

    * **single-dash flag groups.** ``consume_optional`` decomposes them, so a zero-argument
      short option consumes its first character and re-forms ``'-' + tail``. Every one of
      ``-p``, ``-n``, ``-e``, ``-u`` is a free prefix: ``-pt3`` set threads to 3, and
      ``-po /elsewhere.vcf`` redirected the artefact.
    * **single-dash abbreviation.** ``allow_abbrev`` applies to these too, and no ``-v`` is
      declared, so ``-v 99999`` and ``-vi 99999`` both reach ``-vid``.

    So a word starting with ``-`` must now be an option adVNTR declares *and* one this
    module does not set. Abbreviations, attached values and flag groups are all refused
    because none of them is an exact declared spelling -- fail-closed by construction
    rather than by enumeration. It over-rejects only spellings argparse itself would
    accept for a permitted option (``-c30``); those fail loudly and are re-spelled.

    Unlike ``threads`` and ``output_format``, a missing key is not an error here. The default
    is an empty extension. Operators can explicitly opt into adVNTR-owned flags such as
    ``-aln``, but VNtyper does not consume its alignment sidecar. adVNTR 2.0.4 can finish the
    genotype VCF and then crash in that optional post-processing path, so enabling it by
    default turns a completed genotype into a correctly propagated subprocess failure.

    Args:
        settings (dict): An ``advntr_settings`` mapping.

    Returns:
        str: The fragment, unchanged. Validated as ``shlex`` reads it -- ``run_command``
        hands the assembled string to bash, which additionally expands ``$``, backticks and
        ``$'...'``, so this is not a claim about shell metacharacters. The fragment is
        operator-controlled configuration, not user input.

    Raises:
        ValueError: If the fragment does not parse as shell words, if any word is an option
            :func:`run_advntr` sets itself, or if any word is not an option adVNTR declares.
    """
    additional = settings.get("additional_commands", "")

    try:
        tokens = shlex.split(additional)
    except ValueError as exc:
        raise ValueError(
            f"advntr_settings['additional_commands'] is not a parseable command fragment ({exc}): "
            f"{additional!r}. It is interpolated into a string run under `bash -c`, so an "
            "unbalanced quote would surface as a shell syntax error inside adVNTR's log."
        ) from exc

    for token in tokens:
        option = advntr_option_token(token)
        if option is None:
            continue
        if option in MANAGED_ADVNTR_OPTIONS:
            owner = MANAGED_ADVNTR_OPTION_OWNERS[option]
            raise ValueError(
                f"advntr_settings['additional_commands'] contains {token!r}, which is adVNTR's "
                f"{option} -- an option run_advntr already sets. adVNTR parses with argparse, where "
                f"the last occurrence wins, so this would silently override {owner}. Set that instead, "
                "and keep additional_commands for flags adVNTR alone owns (such as -aln)."
            )
        if option not in ADVNTR_EXTRA_OPTIONS:
            raise ValueError(
                f"advntr_settings['additional_commands'] contains {token!r}, which adVNTR's `genotype` "
                "does not declare as an option. Abbreviations, attached values and single-dash flag "
                "groups are refused even when argparse would accept them, because each is a way to "
                "reach an option run_advntr sets -- `-pt3` sets the thread count and `-v` reaches "
                "-vid. Spell the option out in full, one word per option."
            )

    return additional


def run_advntr(db_file, sorted_bam, output, output_name, config, cwd=None, pipeline_threads=1):
    """
    Run adVNTR genotyping using the specified database file and BAM file, fetching settings from advntr_config.

    Args:
        db_file (str): Path to the adVNTR VNTR database file.
        sorted_bam (str): Path to the sorted BAM file.
        output (str): Directory where the results will be saved.
        output_name (str): Base name for the output files.
        config (dict): Main configuration dictionary.
        pipeline_threads (int): The pipeline's ``--threads``. Used only when
            ``advntr_settings['threads']`` is ``null``, which means "inherit".

    Returns:
        int: 0 on success, or 1 for the pre-command input-validation failures.

    Raises:
        RuntimeError: When the critical adVNTR command exits non-zero or reports
            failure. This propagates so no partial output can reach result parsing.
        ValueError: From :func:`resolve_advntr_threads` or
            :func:`resolve_additional_commands`, for a configured value this module refuses.
            Deliberately not turned into a return code so every caller must handle the
            configuration failure explicitly; ``cli.py`` catches ``ValueError`` and exits
            1 with the message.
        KeyError: From :func:`resolve_advntr_threads` or :func:`advntr_output_extension`,
            for a *missing* configuration key. ``cli.py`` does not catch this, so a partial
            mapping surfaces as a traceback rather than a clean exit -- the same shape
            ``output_format`` already had.
    """
    advntr_path = config["tools"]["advntr"]

    # Configuration is authoritative for `threads` and `output_format`. Both used to be read
    # with `.get` fallbacks that contradicted the shipped file -- `threads` defaulted to 8
    # while advntr_config.json sets 1, and `output_format` defaulted to "tsv" while it sets
    # "vcf" -- so dropping either key changed the emitted command with no error at all. That
    # is the pattern already rejected for the calibration constants in
    # confidence_assignment.py:108-111. A partial mapping is not supported input, so failing
    # loudly beats emitting a command nobody chose. (Note this sidecar is NOT what `--config-path`
    # replaces -- that flag loads the main config.json, while `load_advntr_config` defaults to the
    # packaged advntr_config.json independently of it.) (#247)
    #
    # `-t` is real as of adVNTR 2.0.0. This comment replaces one that said the opposite,
    # and the opposite was true at the time: `-t` set only `settings.CORES`
    # (advntr_commands.py:72-74), whose readers were VNTR model construction
    # (models.py:66,69,72,74) and two PacBio-only functions (vntr_finder.py:870,889) --
    # none on the `genotype -fs` path VNtyper uses. Wiring `--threads` to it then would
    # have changed the command without changing the runtime, which is why #215 was closed
    # as resolved-by-investigation and why #254 pinned the commit.
    #
    # adVNTR 2.0.0 moved the Viterbi DP into a `nogil` block and threaded the read loop,
    # so `-t N` now genuinely parallelises decoding: measured 203.3 s -> 1.71 s at `-t 16`
    # on example_7a61_hg19_subset.bam, with an identical selection digest at every thread
    # count. Do not restore the old wording without also re-pinning GIT_COMMIT to a
    # pre-2.0.0 revision.
    #
    # `null` means "inherit the pipeline's --threads"; an explicit integer overrides it.
    # A missing key still raises, because a partial mapping is not supported input.
    #
    # Note `--jobs J` x `-t T` oversubscribes. The golden-cohort runner already models
    # that with a separate `advntr_case_threads` knob, which remains the control.
    threads = resolve_advntr_threads(advntr_settings, pipeline_threads)

    # Optional operator-supplied flags for adVNTR, refusing any option this function
    # already sets. The shipped and missing-key defaults are empty: VNtyper does not read
    # adVNTR's `-aln` sidecar, and adVNTR 2.0.4 can crash while generating it after the VCF
    # is complete. Explicit opt-in remains supported through this extension point. See
    # `resolve_additional_commands`: argparse lets the last occurrence win, so an unchecked
    # fragment overrides the thread count validated one line above.
    additional_commands = resolve_additional_commands(advntr_settings)

    # Determine the output extension. Shared with pipeline.py, which reconstructs this path.
    output_ext = advntr_output_extension(advntr_settings)

    # Set the VNTR ID from config file or default to 25561
    vid = advntr_settings.get("vid", 25561)

    # ---------------------------------------------------------------------
    # Validate input paths before proceeding
    # ---------------------------------------------------------------------
    if not os.path.isfile(db_file):
        logger.critical(f"VNTR database file not found: {db_file}")
        return 1
    if not os.path.isfile(sorted_bam):
        logger.critical(f"Sorted BAM file not found: {sorted_bam}")
        return 1
    if not os.path.isdir(output):
        logger.warning(f"Output directory does not exist, creating: {output}")
        try:
            os.makedirs(output, exist_ok=True)
        except Exception as e:
            logger.critical(f"Could not create output directory {output}: {e}")
            return 1

    # `run_command` runs this as one string under bash (trap 9), so quoting can only
    # happen here. Paths, the sample-derived output name and the thread count are
    # quoted; `advntr_path` and non-empty `additional_commands` are not, because both hold
    # *command fragments* from config.json - `advntr` is "mamba run -n envadvntr advntr"
    # (trap 6) and `additional_commands` can be an explicit flag list such as "-aln".
    # Quoting either would collapse it into a single token bash then looks for as one
    # binary or one argument. They are operator-controlled configuration, not user input.
    additional_fragment = f" {additional_commands}" if additional_commands else ""
    advntr_command = (
        f"{advntr_path} genotype -fs -vid {quote_path(vid)} "
        f"--alignment_file {quote_path(sorted_bam)} "
        f"-o {quote_path(f'{output}/{output_name}_adVNTR{output_ext}')} "
        f"-m {quote_path(db_file)} --working_directory {quote_path(output)} "
        f"-t {quote_path(threads)}{additional_fragment}"
    )

    # Define log file for adVNTR output
    log_file = os.path.join(output, f"{output_name}_advntr.log")

    logger.info("Launching adVNTR genotyping...")
    logger.debug(f"Command: {advntr_command}")

    try:
        succeeded = run_command(advntr_command, log_file, critical=True, cwd=cwd)
    except RuntimeError:
        logger.error(f"adVNTR genotyping failed; adVNTR output is in {log_file}.")
        raise
    if not succeeded:
        msg = "adVNTR genotyping failed. Check the log for details."
        logger.error(msg)
        raise RuntimeError(msg)

    logger.info("adVNTR genotyping of MUC1-VNTR completed successfully.")
    return 0


#: Offset of the insertion arm's accepted-magnitude series: a **net gain** of ``3n+1``
#: bases. See :func:`accepted_frame_magnitudes`.
INSERTION_FRAME_OFFSET = 1

#: Offset of the deletion arm's accepted-magnitude series: a **net loss** of ``3n+2``
#: bases. See :func:`accepted_frame_magnitudes`.
DELETION_FRAME_OFFSET = 2


def accepted_frame_magnitudes(offset: int) -> np.ndarray:
    """The ``|Net_indel_length|`` values one arm of the pathogenic-frame filter accepts.

    The series is built from the two settings ``advntr_config.json`` exposes, read from
    the module-level ``advntr_settings`` global at call time (trap 1 in AGENTS.md: the
    filter reads the *derived* global, so a test must patch that one):

    * ``frameshift_multiplier`` (default 3) -- the codon width;
    * ``max_frameshift`` (default 100) -- how many terms of the series are accepted, so
      the largest accepted magnitude is ``(max_frameshift - 1) * multiplier + offset``.

    Args:
        offset (int): :data:`INSERTION_FRAME_OFFSET` for the insertion arm,
            :data:`DELETION_FRAME_OFFSET` for the deletion arm.

    Returns:
        np.ndarray: The accepted magnitudes as strings, matching the ``str`` dtype of the
            ``frame`` column they are tested against.
    """
    max_frameshift = advntr_settings.get("max_frameshift", 100)
    frameshift_multiplier = advntr_settings.get("frameshift_multiplier", 3)
    return (np.arange(max_frameshift) * frameshift_multiplier + offset).astype(str)


def derive_indel_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename adVNTR's raw columns and derive the indel lengths the filter decides on.

    Shared by :func:`advntr_processing_del` and :func:`advntr_processing_ins` so the two
    arms cannot derive the decision quantity differently. Adds:

    * ``Deletion_length`` -- ``Variant.str.count("D")``, the number of ``D`` events;
    * ``Insertion_length`` -- ``Variant.str.count("I")``, the number of ``I`` events
      (derived for debuggability; nothing reads it);
    * ``Insertion_len`` -- inserted bases, the sum over every ``LEN`` token (#192);
    * ``Net_indel_length`` -- the **signed** net change,
      ``Insertion_len - Deletion_length``. This is the quantity the pathogenic-frame
      decision is made on;
    * ``frame`` -- ``|Net_indel_length|`` as a string, the magnitude tested against
      :func:`accepted_frame_magnitudes`. It deliberately does **not** carry the sign; the
      sign is tested separately, by each arm, on ``Net_indel_length``.

    Args:
        df (pd.DataFrame): A raw adVNTR frame carrying a ``State`` column.

    Returns:
        pd.DataFrame: A copy with the columns above added and ``State``/``Pvalue\\n``
            renamed. No rows are dropped.
    """
    df1 = df.copy()
    df1.rename(columns={"State": "Variant", "Pvalue\n": "Pvalue"}, inplace=True)
    logger.debug("Renamed columns: 'State' -> 'Variant', 'Pvalue\\n' -> 'Pvalue'.")
    df1["Deletion_length"] = df1["Variant"].str.count("D")
    df1["Insertion_length"] = df1["Variant"].str.count("I")
    logger.debug("Calculated 'Deletion_length' and 'Insertion_length'.")
    df1["Insertion_len"] = df1["Variant"].map(sum_insertion_lengths).astype(int)
    logger.debug("Summed every 'LEN' token in 'Variant' into 'Insertion_len' (#192).")
    df1["Deletion_length"] = df1["Deletion_length"].astype(int)
    logger.debug("Converted 'Insertion_len' and 'Deletion_length' to integers.")
    df1["Net_indel_length"] = df1["Insertion_len"] - df1["Deletion_length"]
    df1["frame"] = df1["Net_indel_length"].abs().astype(str)
    logger.debug(f"Computed net indel lengths; sample: {df1['Net_indel_length'].head().tolist()}")
    return df1


# ---------------------------------------------------------------------------
# The pathogenic-frame filter (#182, decided 2026-08-06)
# ---------------------------------------------------------------------------
#
# The decision quantity is the SIGNED net change in bases,
#
#     Delta = Net_indel_length = Insertion_len - Deletion_length,
#
# and the pathogenic ADTKD-MUC1 condition is Delta = +1 (mod 3) -- the frame that yields
# the toxic MUC1-fs neo-protein, classically exemplified by dupC. Written out:
#
#     a net INSERTION of 3n+1 bases   ->  Delta = +1 (mod 3)   PATHOGENIC FRAME
#     a net DELETION  of 3n+2 bases   ->  Delta = -(3n+2) = +1 (mod 3)   PATHOGENIC FRAME
#     a net INSERTION of 3n+2 bases   ->  Delta = +2 (mod 3)   frameshift, other frame
#     a net DELETION  of 3n+1 bases   ->  Delta = +2 (mod 3)   frameshift, other frame
#     a net change of 3n bases        ->  Delta =  0 (mod 3)   in frame, not a frameshift
#
# So VNtyper does NOT filter on "changes the reading frame at all" -- it filters on
# "enters the pathogenic +1 frame", which is the narrower of the two. Rows in the other
# frame are genuine frameshifts and are dropped on purpose; that frame has not been
# established as pathogenic in ADTKD-MUC1 patients. @hassansaei on #182: "keep the same
# 3n+1 / 3n+2 rule for adVNTR as for Kestrel (#181). This is intentional shared
# convention, not something to relax independently."
#
# THE SIGN IS LOAD-BEARING AND USED TO BE DISCARDED. Both arms previously tested
# abs(Insertion_len - Deletion_length) against their series and guarded only on
# "has at least one deleted base" / "has at least one inserted base". A *mixed* state
# satisfies both guards, so either arm could admit it on the magnitude alone and the
# opposite, non-pathogenic frame was accepted -- e.g. I9_2_A_LEN3&D50_2 (3 inserted,
# 1 deleted, Delta = +2) was reported via the DELETION arm because |+2| = 2 is in the
# 3n+2 series. Each arm therefore now tests the sign of Net_indel_length as well, which
# is what makes the pair equivalent to the single signed rule Delta % 3 == 1. The sign
# test also subsumes the old presence guards: Delta < 0 implies Deletion_length >= 1 and
# Delta > 0 implies Insertion_len >= 1.
#
# This is the SAME rule Kestrel applies in scoring.extract_frameshifts (#181), which
# derives its own signed delta from the REF/ALT lengths -- direction = sign(Delta),
# frameshift_amount = |Delta| % 3 -- and was already correct for mixed indels. The two
# are asserted to agree, on pure AND mixed states, by
# tests/unit/test_frameshift_convention_parity.py.
#
# WHY Deletion_length = Variant.str.count("D") IS A COUNT OF DELETED BASES (#202, settled
# 2026-08-12 against adVNTR at berntpopp/adVNTR@05fd98a, branch enhanced_hmm -- the exact
# revision vntyper/dependencies/advntr/install_advntr.cfg pins, verified byte-identical to
# the installed advntr-1.3.3 egg).
#
# The subtraction above mixes two quantities that LOOK like different units: Insertion_len
# is a summed bp count (every LEN token, #192) and Deletion_length is a count of D
# characters. They are the same unit, because one D token always removes exactly one base:
#
#   1. One delete state per reference column, by construction. profile_hmm.py:39-47 emits
#      one 'D' + index per single gap character; hmm_utils.py:507 builds one
#      State(None, name='D%s_%s') per match column. State(None, ...) is silent -- it emits
#      nothing and skips exactly one column.
#   2. Delete states have NO self-loop, insert states do. hmm_utils.py:542,674 add
#      insert->itself transitions, which is precisely why an I token needs _LEN{k}. The
#      only D transition, hmm_utils.py:551,683, is delete_states[i-1] -> delete_states[i],
#      i.e. it advances a column. Deleting k bases therefore visits k distinct states.
#   3. The State string encodes a k-base deletion as k '&'-joined D tokens.
#      vntr_finder.py:606-612 joins consecutive D events; every _LEN write site
#      (vntr_finder.py:587,594,625,633,656; hmm_alignment.py:322,330,372,385,412) is
#      guarded by startswith("I"). No branch anywhere appends a length to a D.
#   4. adVNTR itself counts it this way: hmm_alignment.py:63-64 renders exactly one '-'
#      per D state, and its own commented-out reference implementation at
#      hmm_alignment.py:128-132 is literally `detected_mutation.count("D")`.
#
# So the two arms are consistent and no bp/event conversion is missing. Do not "fix" this
# by parsing a length out of a D token -- there is none to parse.
#
# KNOWN LIMITATION, opposite direction, and NOT a VNtyper defect: adVNTR's
# mutation_count_temp (vntr_finder.py:487) is keyed by state name and reset per read, not
# per repeat unit. One read spanning two RUs of the same pattern that deletes the same
# column in both emits that D token once, so count("D") can UNDERSTATE such a case. That is
# a representational limit of the State string, upstream of anything here.


def advntr_processing_del(df):
    """
    Keep the adVNTR calls whose **net** change is a deletion in the pathogenic frame.

    A row survives when ``Net_indel_length < 0`` (more bases deleted than inserted) and
    ``|Net_indel_length|`` is one of the ``3n+2`` magnitudes
    :func:`accepted_frame_magnitudes` yields, i.e. ``Net_indel_length % 3 == 1``. See the
    block comment above for why the sign test is not optional.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only those net deletions that enter
            the pathogenic frame.
    """
    logger.debug("Starting deletion processing.")
    df1 = derive_indel_columns(df)
    del_frame = accepted_frame_magnitudes(DELETION_FRAME_OFFSET)
    logger.debug(f"Accepted net-deletion magnitudes (first 5): {del_frame[:5].tolist()}")
    df1 = df1[(df1["Net_indel_length"] < 0) & df1["frame"].isin(del_frame)]
    logger.debug(f"Filtered DataFrame shape after deletion processing: {df1.shape}")
    return df1


def advntr_processing_ins(df):
    """
    Keep the adVNTR calls whose **net** change is an insertion in the pathogenic frame.

    A row survives when ``Net_indel_length > 0`` (more bases inserted than deleted) and
    ``Net_indel_length`` is one of the ``3n+1`` magnitudes
    :func:`accepted_frame_magnitudes` yields, i.e. ``Net_indel_length % 3 == 1``. See the
    block comment above for why the sign test is not optional.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only those net insertions that enter
            the pathogenic frame.
    """
    logger.debug("Starting insertion processing.")
    df1 = derive_indel_columns(df)
    ins_frame = accepted_frame_magnitudes(INSERTION_FRAME_OFFSET)
    logger.debug(f"Accepted net-insertion magnitudes (first 5): {ins_frame[:5].tolist()}")
    df1 = df1[(df1["Net_indel_length"] > 0) & df1["frame"].isin(ins_frame)]
    logger.debug(f"Filtered DataFrame shape after insertion processing: {df1.shape}")
    return df1


def load_ru_sequences(ru_fasta_path):
    """
    Load repeat unit (RU) sequences from a FASTA file.

    Args:
        ru_fasta_path (str): Path to the RU FASTA file.

    Returns:
        dict: A dictionary mapping RU identifier (as string) to its sequence.
    """
    ru_dict = {}
    with open(ru_fasta_path) as f:
        current_ru = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_ru and seq_lines:
                    ru_dict[current_ru] = "".join(seq_lines)
                header = line[1:]
                current_ru = header.removeprefix("RU")
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_ru and seq_lines:
            ru_dict[current_ru] = "".join(seq_lines)
    return ru_dict


def annotate_advntr_variants(variant_series, ru_fasta_path):
    """Annotate adVNTR variants with RU, POS, REF, and ALT.

    RU and POS are derived from the variant state string. Only REF and ALT require the
    repeat-unit FASTA.

    Args:
        variant_series (pd.Series): Series of variant strings (possibly with multiple parts separated by '&').
        ru_fasta_path (str): Path to the RU FASTA file.

    Returns:
        tuple: Four lists corresponding to RU, POS, REF, and ALT annotations.
    """
    variant_values = list(variant_series)
    ru_dict = load_ru_sequences(ru_fasta_path)
    ru_annotations, pos_annotations = derive_ru_and_pos(variant_values)
    ref_annotations = []
    alt_annotations = []

    for variant in variant_values:
        ref_parts = []
        alt_parts = []
        parts = variant.split("&") if isinstance(variant, str) else [""]
        for part in parts:
            part = part.strip()
            ins_match = INSERTION_PATTERN.match(part)
            del_match = DELETION_PATTERN.match(part)
            if ins_match:
                pos_val = int(ins_match.group(1))
                ru_val = ins_match.group(2)
                inserted_base = ins_match.group(3)
                ins_len = int(ins_match.group(4))
                ru_seq = ru_dict.get(ru_val, "")
                ref_base = ru_seq[pos_val - 1] if ru_seq and pos_val - 1 < len(ru_seq) else "."
                alt_val = ref_base + inserted_base * ins_len
                ref_parts.append(ref_base)
                alt_parts.append(alt_val)
            elif del_match:
                pos_val = int(del_match.group(1))
                ru_val = del_match.group(2)
                ru_seq = ru_dict.get(ru_val, "")
                if ru_seq and pos_val > 1 and pos_val - 1 < len(ru_seq):
                    prev_base = ru_seq[pos_val - 2]
                    del_base = ru_seq[pos_val - 1]
                    ref_allele = prev_base + del_base
                    alt_allele = prev_base
                else:
                    ref_allele = "."
                    alt_allele = "."
                ref_parts.append(ref_allele)
                alt_parts.append(alt_allele)
            else:
                ref_parts.append(".")
                alt_parts.append(".")
        ref_annotations.append(",".join(ref_parts))
        alt_annotations.append(",".join(alt_parts))

    return ru_annotations, pos_annotations, ref_annotations, alt_annotations


def process_advntr_output(output_path, output, output_name, config=None):
    """
    Process the adVNTR output to extract relevant information and generate final results.

    Repeat unit (RU) and position (POS) are always derived from each variant's state
    string. When configuration includes a valid 'reference_data.code_adVNTR_RUs' FASTA,
    the function also annotates REF and ALT values.

    Every final output retains the legacy caller and nomenclature columns. Positive
    results additionally append the four molecular-identity publication columns;
    the negative placeholder retains its unchanged pre-identity schema.

    If the VCF data is empty, a negative result is generated immediately with
    'VID' set to "Negative" and all other columns set to "None", and further processing is skipped.

    Args:
        output_path (str): Path to the adVNTR output file.
        output (str): Directory where the final results will be saved.
        output_name (str): Base name for the output files.
        config (dict, optional): Main configuration dictionary.

    Raises:
        ValueError: If the configured flag rules are invalid for the adVNTR result schema.
        RuntimeError: If the source is absent or unreadable, processing fails, or the
            complete result cannot be published atomically.
    """
    output_result_path = Path(output) / f"{output_name}_adVNTR_result.tsv"
    invalidate_advntr_artifact(output_result_path)

    flag_columns = set(ADVNTR_FLAG_COLUMNS)
    if config:
        ru_fasta_path = config.get("reference_data", {}).get("code_adVNTR_RUs")
        if ru_fasta_path and os.path.exists(ru_fasta_path):
            flag_columns.update({"REF", "ALT"})
    compiled_flag_rules = compile_flag_rules(advntr_config.get("flagging_rules", {}), flag_columns)

    if not os.path.exists(output_path):
        message = f"adVNTR output file {output_path} not found!"
        logger.error(message)
        raise RuntimeError(message)

    logger.info("Processing adVNTR result...")

    try:
        with open(output_path) as file:
            content = file.readlines()

        # Replace header to ensure consistency
        content = [line.replace("#VID", "VID") if line.startswith("#VID") else line for line in content]
        with open(output_path, "w") as file:
            file.writelines(content)
    except OSError as e:
        message = f"Error reading adVNTR output: {e}"
        logger.error(message)
        raise RuntimeError(message) from e

    try:
        logger.info("Loading data into DataFrame...")
        df = pd.read_csv(output_path, sep="\t", comment="#")
        logger.info(f"Data loaded successfully with shape: {df.shape}")
        logger.debug(f"First few rows of the DataFrame:\n{df.head()}")
    except Exception as e:
        message = f"Error loading data into DataFrame: {e}"
        logger.error(message)
        raise RuntimeError(message) from e

    # Immediately check if the loaded DataFrame is empty
    final_columns = [
        "VID",
        "Variant",
        "NumberOfSupportingReads",
        "MeanCoverage",
        "Pvalue",
        "RU",
        "POS",
        "REF",
        "ALT",
        "Flag",
        *NOMENCLATURE_COLUMNS,
    ]
    if df.empty:
        logger.warning("VCF file is empty. Generating default negative result.")
        advntr_concat = pd.DataFrame(
            [
                {
                    "VID": "Negative",
                    "Variant": "Not applicable",
                    "NumberOfSupportingReads": "Not applicable",
                    "MeanCoverage": "Not applicable",
                    "Pvalue": "Not applicable",
                    "RU": "Not applicable",
                    "POS": "Not applicable",
                    "REF": "Not applicable",
                    "ALT": "Not applicable",
                    "Flag": "Not applicable",
                    # Empty, not "Not applicable": these five are nullable by
                    # contract and are never padded with a placeholder string.
                    **dict.fromkeys(NOMENCLATURE_COLUMNS, ""),
                }
            ]
        )
        advntr_concat = advntr_concat[final_columns]
        publish_advntr_result(advntr_concat, output_result_path)
        logger.info(f"Processed adVNTR results saved to {output_result_path}")
        cleanup_files(output, output_name)
        return

    try:
        logger.info("Processing deletions...")
        df_del = advntr_processing_del(df)

        logger.info("Processing insertions...")
        df_ins = advntr_processing_ins(df)

        logger.info("Concatenating deletions and insertions...")
        advntr_concat = pd.concat([df_del, df_ins], axis=0)

        positive_result = not advntr_concat.empty
        if advntr_concat.empty:
            logger.warning("No pathogenic variant found after filtering. Generating default negative result.")
            advntr_concat = pd.DataFrame(
                [
                    {
                        "VID": "Negative",
                        "Variant": "Not applicable",
                        "NumberOfSupportingReads": "Not applicable",
                        "MeanCoverage": "Not applicable",
                        "Pvalue": "Not applicable",
                        "RU": "Not applicable",
                        "POS": "Not applicable",
                        "REF": "Not applicable",
                        "ALT": "Not applicable",
                        "Flag": "Not applicable",
                        # Empty, not "Not applicable": these five are nullable by
                        # contract and are never padded with a placeholder string.
                        **dict.fromkeys(NOMENCLATURE_COLUMNS, ""),
                    }
                ]
            )
        else:
            from vntyper.scripts.identity_candidates import translation_component_from_config
            from vntyper.scripts.molecular_identity_presentation import IDENTITY_COLUMNS
            from vntyper.scripts.nomenclature import load_nomenclature_config

            base_columns = [
                "VID",
                "Variant",
                "NumberOfSupportingReads",
                "MeanCoverage",
                "Pvalue",
                "Net_indel_length",
            ]
            advntr_concat = advntr_concat[base_columns]
            logger.info("Removing duplicates...")
            advntr_concat.drop_duplicates(subset=["VID", "Variant", "NumberOfSupportingReads"], inplace=True)

            # RU and POS come from the state string and must exist before flagging. Only
            # REF and ALT require the repeat-unit FASTA.
            ru_ann, pos_ann = derive_ru_and_pos(advntr_concat["Variant"])
            advntr_concat["RU"] = ru_ann
            advntr_concat["POS"] = pos_ann
            if config:
                ru_fasta_path = config.get("reference_data", {}).get("code_adVNTR_RUs")
                if ru_fasta_path and os.path.exists(ru_fasta_path):
                    logger.info("Annotating variants with RU-level REF and ALT bases.")
                    _ru_ann, _pos_ann, ref_ann, alt_ann = annotate_advntr_variants(
                        advntr_concat["Variant"], ru_fasta_path
                    )
                    advntr_concat["REF"] = ref_ann
                    advntr_concat["ALT"] = alt_ann

            # Apply flagging rules if available
            if compiled_flag_rules.rules:
                logger.info("Applying flagging rules to adVNTR output.")
                from vntyper.scripts.flagging import add_flags

                advntr_concat = add_flags(advntr_concat, compiled_flag_rules)

            # Name the variants. Done before the "ensure all columns present" sweep
            # below so the five land as computed values, not as "Not applicable".
            identity_component = translation_component_from_config(load_nomenclature_config())
            advntr_concat = annotate_advntr_frame(advntr_concat, identity_component=identity_component)

            # Ensure all final columns are present
            for col in final_columns:
                if col not in advntr_concat.columns:
                    advntr_concat[col] = "" if col in NOMENCLATURE_COLUMNS else "Not applicable"

        publication_columns = (
            [*final_columns, *IDENTITY_COLUMNS, EVIDENCE_DISPOSITION_COLUMN] if positive_result else final_columns
        )
        advntr_concat = advntr_concat[publication_columns]
    except Exception as e:
        message = f"Error during processing of deletions and insertions: {e}"
        logger.error(message)
        raise RuntimeError(message) from e

    publish_advntr_result(advntr_concat, output_result_path)
    logger.info(f"Processed adVNTR results saved to {output_result_path}")

    cleanup_files(output, output_name)


def cleanup_files(output, output_name):
    """
    Clean up intermediate files.

    Args:
        output (str): The output directory.
        output_name (str): The base name for the output files.
    """
    logger.info("Intermediate files cleaned up.")
