"""
cli_handlers.py

Module Purpose:
---------------
One handler function per ``vntyper`` subcommand, so that ``vntyper/cli.py``'s
``main()`` is nothing but parse-then-dispatch.

Every subcommand body used to live inline in ``main()``, which meant any two
changes to two different subcommands were changes to the same function. Splitting
them apart lets each subcommand be edited, tested and reverted on its own.

The ``report`` handler deliberately lives in its own module
(:mod:`vntyper.scripts.cli_report`) rather than here: it is the one subcommand
known to be broken (AGENTS.md trap 11, fixed in #179) and is owned separately.

Logging is **not** configured here. ``cli.py`` remains the sole place that calls
``setup_logging()``; handlers receive the already-resolved log level and log file
so they can forward them to the stage they invoke.
"""

import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Any, Protocol

from vntyper.scripts.artifact_names import validate_output_name

# Imported names are module attributes and module globals alike, which is what keeps
# the suite's monkeypatch seams and the handlers' own lookups working (#262).
from vntyper.scripts.cli_lazy_imports import (
    aggregate_cohort,
    install_references_main,
    run_online_mode,
)
from vntyper.scripts.cohort_pseudonyms import _mapping_at
from vntyper.scripts.pipeline import run_pipeline
from vntyper.scripts.reference_registry import get_reference_source, physical_reference_id, reference_keys
from vntyper.scripts.reference_resolution import ResolvedReference, resolve_from_mapping
from vntyper.scripts.report_assets import DEFAULT_REPORT_IGV

logger = logging.getLogger(__name__)


#: The optional stages ``--extra-modules`` can name. ``pipeline.py`` tests membership
#: against exactly these two strings (``"advntr" in extra_modules``,
#: ``"shark" in extra_modules``), so a name outside this set is a silent no-op rather
#: than an error - which is what makes validating here worth doing.
#: ``tests/unit/test_cli_handlers.py`` asserts the two lists stay in step.
KNOWN_EXTRA_MODULES: tuple[str, ...] = ("advntr", "shark")


def normalise_extra_modules(values: list[Any] | None) -> list[str]:
    """Turn every accepted spelling of ``--extra-modules`` into the pipeline's own.

    ``--extra-modules`` is an ``append`` action, so repeating it yields
    ``["advntr", "shark"]`` - but the documented comma form yields the single string
    ``"advntr,shark"``, and ``pipeline.py`` tests membership by exact string. That
    made ``--extra-modules advntr,shark`` match *neither* module and produce a
    Kestrel-only run with exit code 0 and a report that reads as a plain negative.
    The same held for a typo: nothing validated the names.

    Splits on commas, trims, drops empties, lower-cases and de-duplicates while
    preserving first-seen order.

    Args:
        values: ``args.extra_modules`` - a list of strings, or of lists of strings
            if a caller built it programmatically. None is treated as empty.

    Returns:
        list[str]: The normalised module names, in the order first named.

    Raises:
        ValueError: If any name is not in :data:`KNOWN_EXTRA_MODULES`.
    """
    modules: list[str] = []
    as_written: dict[str, str] = {}
    for value in values or []:
        for item in value if isinstance(value, list) else [value]:
            for part in str(item).split(","):
                written = part.strip()
                if not written:
                    continue
                name = written.lower()
                as_written.setdefault(name, written)
                if name not in modules:
                    modules.append(name)

    unknown = [as_written[name] for name in modules if name not in KNOWN_EXTRA_MODULES]
    if unknown:
        msg = (
            f"Unknown --extra-modules value(s): {', '.join(sorted(unknown))}. "
            f"Available modules: {', '.join(KNOWN_EXTRA_MODULES)}. "
            "An unrecognised name used to be accepted and then ignored, producing a Kestrel-only run "
            "that reports as a negative genotype."
        )
        logger.error(msg)
        raise ValueError(msg)

    logger.debug(f"extra_modules normalised to {modules}")
    return modules


def get_conf(config: dict[str, Any], key: str, fallback: Any) -> Any:
    """Read a CLI default out of ``config["default_values"]``.

    Args:
        config: The loaded configuration mapping.
        key: The ``default_values`` key to read.
        fallback: Value to return when the key is absent.

    Returns:
        Any: The configured value, or ``fallback``.
    """
    return config.get("default_values", {}).get(key, fallback)


class CommandHandler(Protocol):
    """The uniform signature every subcommand handler is dispatched through.

    The signature is uniform rather than minimal so that ``main()`` can dispatch
    through a table and never needs editing when a handler's needs change. A
    handler that does not need a parameter simply ignores it.
    """

    def __call__(
        self,
        args: argparse.Namespace,
        config: dict[str, Any],
        parser: argparse.ArgumentParser,
        log_level_value: int,
        log_file_str: str | None,
    ) -> None:
        """Run the subcommand.

        Args:
            args: The parsed arguments.
            config: The loaded configuration mapping.
            parser: The parser, for ``parser.error()`` on invalid combinations.
            log_level_value: The already-resolved ``logging`` level.
            log_file_str: The already-resolved log file path, or None.
        """
        ...  # pragma: no cover - structural typing only


def handle_install_references(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``install-references`` subcommand.

    Args:
        args: The parsed arguments.
        config: Unused; present for the uniform handler signature.
        parser: Used to reject ``--release-spec`` without ``--from-source``.
        log_level_value: Already applied by `cli.py`; the installer preserves it.
        log_file_str: Already applied by `cli.py`; the installer preserves its handler.
    """
    if args.derive_only and args.from_source:
        # Both build the derived files, but --derive-only deliberately downloads nothing.
        # Accepting both would make one of them silently meaningless.
        parser.error("--derive-only cannot be combined with --from-source; --from-source already derives")

    if args.release_spec and not args.from_source:
        # Silently ignoring a flag the user typed is worse than refusing it: only the
        # source path reads a release spec, so without --from-source the pinned URLs and
        # digests would have no effect at all. Usage errors exit 2, as argparse does.
        parser.error("--release-spec requires --from-source; it only affects the from-source build")

    install_references_main(
        output_dir=args.output_dir,
        config_path=args.config_path,
        skip_indexing=args.skip_indexing,
        index_threads=args.threads,
        aligners_to_use=args.aligners,
        references_to_process=args.references,
        from_source=args.from_source,
        release_spec_path=args.release_spec,
        derive_only_mode=args.derive_only,
    )
    sys.exit(0)


def _resolve_bwa_reference(
    config: dict[str, Any], reference_assembly: str, *, required: bool = True
) -> ResolvedReference | None:
    """Resolve the BWA reference config entry for an assembly, failing closed.

    Shared by :func:`select_bwa_reference`, which reduces this to a path, and by
    :func:`handle_pipeline`, which also needs the matched key and fallback status to
    record in the run summary (``reference_key_used``, ``reference_source_effective``) -
    re-deriving them independently at the call site would risk disagreeing with what was
    actually resolved.

    Args:
        config: Pipeline configuration.
        reference_assembly: Supported assembly label.
        required: Whether an unresolved reference raises (True, the run path) or returns
            None (False, the logging-safety guard).

    Returns:
        ResolvedReference | None: The matched key, its value and whether it was a
        UCSC-family fallback, or None when nothing resolves and `required` is False.

    Raises:
        ValueError: If no configured key resolves, if the key that resolved is present
            with an explicit ``null`` (deliberately disabled), or if a resolved path
            names no file on disk, and `required` is True in each case.
        ValueError: If `reference_assembly` is not a supported label. This must NOT be
            swallowed by callers - an unknown assembly is a configuration error, not a
            missing file.
    """
    resolved = resolve_from_mapping("bwa", reference_assembly, config.get("reference_data", {}))
    keys = ", ".join(reference_keys("bwa", reference_assembly))

    # A present-but-non-string value (an int, a dict, ...) is malformed config, not a
    # usable path. Treating it as unresolved - rather than handing it to `Path()` further
    # down the call chain - is what lets `required=False` return None instead of the
    # caller crashing on a type it never asked to handle; `required=True` still fails
    # closed, just with a clear message instead of a `Path`/`TypeError` two frames later.
    if resolved is not None and isinstance(resolved.value, str) and resolved.value:
        if not Path(resolved.value).exists():
            # Presence beat truthiness to get here, but a configured path that names no
            # file on disk is not a usable reference either. A default `install-references`
            # run only installs hg19/hg38, so a shipped config.json that also declares
            # bwa_reference_GRCh38 as a real (but not-yet-installed) relative path used to
            # resolve here, is_fallback=False, no warning - and the run died several stages
            # later with a message naming neither the assembly nor the remedy. Falling
            # through to the next tier (the UCSC-family key, when one exists) instead of
            # raising here would silently align a GRCh38-labelled run against UCSC sequence
            # - precisely the defect this milestone exists to kill - so this fails closed at
            # whichever tier `resolve_from_mapping` already matched, rather than re-walking
            # `reference_keys` to find one whose file happens to exist.
            if not required:
                return None
            message = (
                f"BWA reference for --reference-assembly {reference_assembly!r} is configured at "
                f"reference_data[{resolved.key!r}] = {resolved.value!r}, but that file does not exist. "
                f"Run `vntyper install-references --references {physical_reference_id(reference_assembly)}` "
                "to install it."
            )
            logger.error(message)
            raise ValueError(message)
        if resolved.is_fallback:
            # Name the effective source plainly. Deriving it from the key suffix would be
            # fragile; the fallback key is always the UCSC family key, so the source is ucsc.
            logger.warning(
                f"--reference-assembly {reference_assembly!r} has no {keys.split(', ')[0]!r} entry; "
                f"falling back to {resolved.key!r}. This run therefore uses 'ucsc' sequence, "
                f"not {get_reference_source(reference_assembly)!r}. "
                f"Run `vntyper install-references` to install the requested reference."
            )
        return resolved

    if not required:
        return None

    if resolved is not None and resolved.value is None:
        # Present with an explicit null is a deliberate "disabled", not an absence - say
        # so, distinctly from the "nothing configured at all" message below, the way
        # `shark_filtering.select_muc1_region_fasta` already distinguishes the two.
        message = (
            f"reference_data[{resolved.key!r}] is null; BWA reference for --reference-assembly "
            f"{reference_assembly!r} is disabled. Tried: {keys}. "
            "Run `vntyper install-references` or set one of those keys."
        )
    else:
        message = (
            f"No BWA reference configured for --reference-assembly {reference_assembly!r}. "
            f"Tried: {keys}. Run `vntyper install-references` or set one of those keys."
        )
    logger.error(message)
    raise ValueError(message)


def select_bwa_reference(config: dict[str, Any], reference_assembly: str, *, required: bool = True) -> str | None:
    """Resolve the BWA reference for an assembly, failing closed.

    Both this and `cli_logging_safety` must agree, or the guard that refuses to let
    `--log-file` name an operator input inspects a different file from the one the run
    opens for writing.

    Args:
        config: Pipeline configuration.
        reference_assembly: Supported assembly label.
        required: Whether an unresolved reference raises (True) or returns None
            (False).

    Returns:
        str | None: Path to the reference FASTA, or None when `required` is False and
        nothing resolves.

    Raises:
        ValueError: If no configured key resolves and `required` is True.
        ValueError: If `reference_assembly` is not a supported label. This must NOT be
            swallowed by callers - an unknown assembly is a configuration error, not a
            missing file.
    """
    resolved = _resolve_bwa_reference(config, reference_assembly, required=required)
    return resolved.value if resolved is not None else None


def handle_pipeline(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``pipeline`` subcommand.

    Args:
        args: The parsed arguments.
        config: The loaded configuration mapping.
        parser: The parser, used to reject invalid input combinations.
        log_level_value: The already-resolved ``logging`` level.
        log_file_str: The already-resolved log file path, or None.
    """
    # If log_file was not explicitly provided and output_dir is set, ensure log_file is correctly set
    if not args.log_file and args.output_dir:
        log_file = Path(args.output_dir) / "pipeline.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Setting log file to {log_file}")

    if args.output_dir is None:
        args.output_dir = get_conf(config, "output_dir", "out")
        logger.debug(f"output_dir set to {args.output_dir}")
    if args.threads is None:
        args.threads = get_conf(config, "threads", 4)
        logger.debug(f"threads set to {args.threads}")
    if args.reference_assembly is None:
        args.reference_assembly = get_conf(config, "reference_assembly", "hg19")
        logger.debug(f"reference_assembly set to {args.reference_assembly}")
    # A value supplied on the command line has already been refused in `cli.main`,
    # before the output directory and the log file were created. This call still
    # runs, because the value can also come from the configuration file, which the
    # CLI check deliberately does not reach into.
    # --output-name is validated, not forwarded: run_pipeline takes no such parameter
    # and three of its consumers name their files from a literal it cannot reach.
    # Refusing a value we cannot honour is the point - silently dropping it produced a
    # run whose artefacts were not where the caller was told to look.
    # See vntyper/scripts/artifact_names.py.
    args.output_name = validate_output_name(args.output_name or get_conf(config, "output_name", None))
    logger.debug(f"output_name is fixed at {args.output_name}")
    if args.archive_format is None:
        args.archive_format = get_conf(config, "archive_format", "zip")
        logger.debug(f"archive_format set to {args.archive_format}")

    flattened_modules = normalise_extra_modules(args.extra_modules)

    # Validate single input type
    input_types = sum(
        [
            1 if args.bam else 0,
            1 if args.cram else 0,
            1 if (args.fastq1 or args.fastq2) else 0,
        ]
    )
    if input_types > 1:
        parser.error("Provide either BAM, CRAM, or FASTQ files (not multiples).")
        logger.debug("Multiple input types detected.")
        sys.exit(1)

    if not args.bam and not args.cram and args.fastq1 is None:
        parser.error("When not providing BAM/CRAM, --fastq1 must be specified; --fastq2 is optional.")
        logger.debug("Missing required FASTQ1 input.")
        sys.exit(1)

    if args.fastq1 is not None and args.fastq2 is None and "shark" in flattened_modules:
        parser.error("SHARK requires paired-end FASTQ input; provide --fastq2 or remove the shark module.")

    # Construct module_args_dict for advntr, etc.
    module_args_dict: dict[str, dict[str, Any]] = {}
    if "advntr" in flattened_modules:
        module_args_dict["advntr"] = {}

        # If the user set advntr_reference somewhere
        if hasattr(args, "advntr_reference"):
            module_args_dict["advntr"]["advntr_reference"] = args.advntr_reference
            delattr(args, "advntr_reference")
            logger.debug(f"advntr_reference set to {args.advntr_reference}")

        # The new coverage parameter:
        if args.advntr_max_coverage:
            module_args_dict["advntr"]["max_coverage"] = args.advntr_max_coverage
            logger.debug(f"advntr_max_coverage set to {args.advntr_max_coverage}")

    else:
        module_args_dict["advntr"] = {}
        logger.debug("advntr module not included.")

    # (#62) If user tries to use 'shark' in BAM/CRAM mode, exit with a warning
    if (args.bam or args.cram) and ("shark" in flattened_modules):
        logger.warning("Shark is not supported in BAM mode; please use FASTQ mode or remove the shark flag.")
        logger.debug("Shark module detected with BAM/CRAM input; exiting.")
        sys.exit(1)

    # Determine which BWA reference to use, resolving membership-first through the
    # shared registry (#163): `cli_logging_safety`'s pre-open guard must resolve to the
    # same file this does, or `--log-file` can end up appending into whichever reference
    # the guard did not check. Only a FASTQ run actually aligns with it -
    # `pipeline_inputs.py` raises for a missing BWA reference only when
    # `input_type == "FASTQ"`, and the BAM/CRAM branches never read it - so `required` is
    # tied to that, not to a blanket True. A BAM/CRAM run resolves best-effort (its result
    # still feeds input-ownership protection and archiving) so a config with no BWA keys
    # at all does not abort a run that never needed one; a FASTQ run still fails closed,
    # with `pipeline_inputs`'s own message one layer down if this resolves to None anyway.
    is_fastq_input = not args.bam and not args.cram
    resolved_bwa_reference = _resolve_bwa_reference(config, args.reference_assembly, required=is_fastq_input)
    bwa_reference = resolved_bwa_reference.value if resolved_bwa_reference is not None else None
    reference_key_used = resolved_bwa_reference.key if resolved_bwa_reference is not None else None
    reference_source_effective = None
    if resolved_bwa_reference is not None:
        reference_source_effective = (
            "ucsc" if resolved_bwa_reference.is_fallback else get_reference_source(args.reference_assembly)
        )
    logger.debug(f"Using BWA reference {reference_key_used}: {bwa_reference}")

    # Where the sample name came from, recorded alongside it (#242). This is the only
    # place that knows: `Path(...).stem` of `S1_R1.fastq.gz` is `S1_R1.fastq`, and the
    # report's documented rule finishes that into `S1` - but only if it can tell a stem
    # from a name the operator typed, which the string cannot say. Deriving the name
    # *here* instead would also close the hole and would rename every FASTQ run's
    # Kestrel output files, because `run_kestrel` builds them from this same value.
    sample_name_is_explicit = args.sample_name is not None
    sample_name_val = args.sample_name
    if sample_name_val is None:
        if args.bam:
            sample_name_val = Path(args.bam).stem
            logger.debug(f"sample_name set from BAM file: {sample_name_val}")
        elif args.cram:
            # CRAM is derived from here for the same reason BAM is, and the arm is
            # not optional (#188). The web worker used to hand every accepted
            # alignment to --bam, so a CRAM reached this branch as a BAM and got
            # its stem; now that it arrives as --cram, without this arm every CRAM
            # run with no explicit --sample-name would fall through to the literal
            # "sample" below -- in the report and in the output filenames.
            sample_name_val = Path(args.cram).stem
            logger.debug(f"sample_name set from CRAM file: {sample_name_val}")
        elif args.fastq1:
            sample_name_val = Path(args.fastq1).stem
            logger.debug(f"sample_name set from FASTQ1 file: {sample_name_val}")
        else:
            sample_name_val = "sample"
            logger.debug(f"sample_name defaulted to: {sample_name_val}")

    # Process the new --summary-formats argument (comma-separated)
    summary_formats = []
    if args.summary_formats:
        summary_formats = [fmt.strip().lower() for fmt in args.summary_formats.split(",") if fmt.strip()]

    run_pipeline(
        bwa_reference=bwa_reference,
        output_dir=Path(args.output_dir),
        extra_modules=flattened_modules,
        module_args=module_args_dict,
        config=config,
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        bam=args.bam,
        cram=args.cram,
        reference_fasta=args.reference_fasta,
        threads=args.threads,
        reference_assembly=args.reference_assembly,
        reference_key_used=reference_key_used,
        reference_source_effective=reference_source_effective,
        fast_mode=args.fast_mode,
        delete_intermediates=args.delete_intermediates,
        archive_results=args.archive_results,
        archive_format=args.archive_format,
        custom_regions=args.custom_regions,
        bed_file=args.bed_file,
        log_level=log_level_value,  # Pass log_level to run_pipeline
        sample_name=sample_name_val,
        sample_name_is_explicit=sample_name_is_explicit,
        log_file=log_file_str,  # Pass the correctly determined log_file
        summary_formats=summary_formats,  # New parameter passed
        # An ordinary run never goes through `vntyper report`, so this is the only
        # route `--report-igv` has to the generator for the reports most people get.
        # `getattr` guards a direct `handle_pipeline()` call whose namespace predates
        # the option (#242).
        report_igv=getattr(args, "report_igv", DEFAULT_REPORT_IGV),
        run_configuration=getattr(args, "run_configuration", None),
        resume=getattr(args, "resume", False),
    )


def handle_cohort(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``cohort`` subcommand.

    Args:
        args: The parsed arguments.
        config: The loaded configuration mapping.
        parser: Unused; present for the uniform handler signature.
        log_level_value: Unused; present for the uniform handler signature.
        log_file_str: Unused; present for the uniform handler signature.
    """
    if args.summary_file is None:
        args.summary_file = get_conf(config, "summary_file", "cohort_summary.html")
        logger.debug(f"summary_file set to {args.summary_file}")

    # Prepare the list of input paths
    input_paths = []
    if args.input_dirs:
        input_paths.extend(args.input_dirs)
        logger.debug(f"Added input_dirs: {args.input_dirs}")
    if args.input_file:
        if not args.input_file.exists():
            logger.error(f"The input file {args.input_file} does not exist.")
            sys.exit(1)
        with open(args.input_file) as f:
            file_lines = [line.strip() for line in f if line.strip()]
            input_paths.extend(file_lines)
        logger.debug(f"Added input_file entries: {file_lines}")

    rare_max_freq = args.rare_allele_max_frequency
    if rare_max_freq is None:
        cohort_cfg = _mapping_at(config, "cohort", "cohort")
        raw_val = cohort_cfg.get("rare_allele_max_frequency", 0.05)
        try:
            val = float(raw_val)
            if not (0.0 <= val <= 1.0) or not math.isfinite(val):
                raise ValueError(f"Value {raw_val!r} is not between 0 and 1.")
            rare_max_freq = val
        except (ValueError, TypeError) as e:
            msg = f"Invalid cohort.rare_allele_max_frequency in config: {raw_val!r} ({e})"
            logger.error(msg)
            raise ValueError(msg) from e

    aggregate_cohort(
        input_paths=input_paths,
        output_dir=Path(args.output_dir),
        summary_file=args.summary_file,
        config=config,
        additional_formats=args.summary_formats,  # Pass the new parameter for extra formats
        pseudonymize_samples=args.pseudonymize_samples,  # Pass the new pseudonymize flag (value or None)
        rare_allele_max_frequency=rare_max_freq,
    )


def handle_online(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``online`` subcommand.

    Args:
        args: The parsed arguments.
        config: The loaded configuration mapping.
        parser: Unused; present for the uniform handler signature.
        log_level_value: Unused; present for the uniform handler signature.
        log_file_str: Unused; present for the uniform handler signature.

    Raises:
        SystemExit: With code 1 when the remote job did not complete.
    """
    # No need to set up logging here; it's already set up in cli.py

    if args.output_dir is None:
        args.output_dir = get_conf(config, "output_dir", "out")
        logger.debug(f"output_dir set to {args.output_dir}")
    if args.reference_assembly is None:
        args.reference_assembly = get_conf(config, "reference_assembly", "hg19")
        logger.debug(f"reference_assembly set to {args.reference_assembly}")
    if args.threads is None:
        args.threads = get_conf(config, "threads", 4)
        logger.debug(f"threads set to {args.threads}")

    try:
        run_online_mode(
            config=config,
            bam=args.bam,
            output_dir=args.output_dir,
            reference_assembly=args.reference_assembly,
            threads=args.threads,
            email=args.email,
            cohort_id=args.cohort_id,
            passphrase=args.passphrase,
            resume=args.resume,
        )
    except RuntimeError:
        # run_online_mode has already logged the reason at ERROR. Exit non-zero so a
        # wrapping `subprocess.run(..., check=True)` sees the failure: this used to
        # return normally and exit 0 on a failed remote job.
        sys.exit(1)
