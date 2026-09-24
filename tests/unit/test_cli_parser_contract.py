"""The full ``vntyper`` argument-parser contract, one table, one comparison.

``tests/unit/test_cli_parser.py`` pins the parser's *shape* - that five subcommands
exist and that a handful of options parse. It cannot catch a changed default, a
changed type or a lost ``required=True``, and those are the changes that break a
caller silently: a ``--threads`` that stops being ``int`` reaches ``run_pipeline`` as
a string, and a ``required`` that is dropped turns a usage error into a ``None``
travelling into the pipeline.

This module states the contract declaratively and compares the whole table at once,
so a change to any option is a diff on one line rather than an absent assertion.
Adding an option is *meant* to fail this test: the new row is the review.

Three things it also pins, which the shape tests do not:

* ``--reference-assembly`` offers exactly what :func:`reference_registry.list_assemblies`
  knows, on **both** subcommands that take it. A registry entry the CLI cannot name is
  unreachable; a CLI choice the registry does not know raises deep inside the run.
* Every subcommand rejects an unknown argument with argparse's exit code 2, rather
  than absorbing it into a positional.
* The global options live on the **top-level** parser only. ``cli.main``'s docstring
  used to claim they worked after the subcommand as well; they never have, because
  ``parent_parser`` is not passed to ``add_parser``.
"""

import argparse

import pytest

from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.reference_registry import list_assemblies

pytestmark = pytest.mark.unit


def test_keep_intermediates_remains_an_accepted_compatibility_flag(capsys: pytest.CaptureFixture[str]) -> None:
    """The append-only CLI contract survives after its inert Python plumbing is removed."""
    parser = build_parser()
    args = parser.parse_args(["pipeline", "--keep-intermediates", "--delete-intermediates"])

    assert args.keep_intermediates is True
    assert args.delete_intermediates is True
    with pytest.raises(SystemExit) as help_exit:
        parser.parse_args(["pipeline", "--help"])
    assert help_exit.value.code == 0
    help_text = " ".join(capsys.readouterr().out.split())
    assert "Compatibility flag: intermediate files" in help_text
    assert "wins when --keep-intermediates is also given" in help_text


#: One row per option: ``dest -> (option strings, action class, type name, default,
#: required, choices, nargs)``. ``type name`` is ``None`` for options argparse stores
#: verbatim. This is the contract; the parser is compared against it wholesale.
ParserRow = tuple

TOP_LEVEL_CONTRACT: dict[str, ParserRow] = {
    "config_path": (("--config-path",), "_StoreAction", "Path", None, False, None, None),
    "log_file": (("-f", "--log-file"), "_StoreAction", None, None, False, None, None),
    "log_level": (
        ("-l", "--log-level"),
        "_StoreAction",
        None,
        None,
        False,
        ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        None,
    ),
    "version": (("-v", "--version"), "_VersionAction", None, "==SUPPRESS==", False, None, 0),
}

ASSEMBLY_CHOICES = ("GRCh37", "GRCh38", "hg19", "hg19_ensembl", "hg19_ncbi", "hg38", "hg38_ensembl", "hg38_ncbi")

#: ``--report-igv``'s choices, stated here as literals rather than imported from
#: ``report_assets``. This is the contract table: a row that asks the implementation
#: what it implements agrees with whatever it implements, including with a mode
#: silently added or removed.
REPORT_IGV_CHOICES = ("embedded", "sidecar", "off")

#: ``--target`` on the operations that still default to the dominance workflow.
DOMINANCE_TARGETS = ("dominance", "callers", "length")

#: ``calibrate optimize``'s axis and objective choices, stated as literals for the same
#: reason as ``REPORT_IGV_CHOICES``. Separate tests below compare them against the axis
#: constants and the objective inventory they are meant to mirror, so a value added in
#: one place and not the other is a failure rather than an agreement with itself.
CUTOFF_AXIS_CHOICES = (
    "depth_floor_linked",
    "gg_gate_independent",
    "depth_score_high",
    "alt_depth_band",
    "var_active_region",
    "advntr_cutoff",
    "advntr_min_support",
)
CUTOFF_OBJECTIVE_CHOICES = (
    "max-sensitivity-at-specificity",
    "youden-j",
    "max-f1",
    "balanced-accuracy",
    "sensitivity",
    "specificity",
)

SUBCOMMAND_CONTRACT: dict[str, dict[str, ParserRow]] = {
    "calibrate": {},
    "pipeline": {
        "advntr_additional_commands": (
            ("--advntr-additional-commands",),
            "_StoreAction",
            "str",
            None,
            False,
            None,
            None,
        ),
        "advntr_max_coverage": (("--advntr-max-coverage",), "_StoreAction", "int", None, False, None, None),
        "archive_format": (("--archive-format",), "_StoreAction", "str", None, False, ("zip", "tar.gz"), None),
        "archive_results": (("--archive-results",), "_StoreTrueAction", None, False, False, None, 0),
        "bam": (("--bam",), "_StoreAction", "str", None, False, None, None),
        "bed_file": (("--bed-file",), "_StoreAction", "Path", None, False, None, None),
        "cram": (("--cram",), "_StoreAction", "str", None, False, None, None),
        "custom_regions": (("--custom-regions",), "_StoreAction", "str", None, False, None, None),
        "decision_profile": (("--decision-profile",), "_StoreAction", "Path", None, False, None, None),
        "research_decision_profile": (
            ("--research-decision-profile",),
            "_StoreAction",
            "Path",
            None,
            False,
            None,
            None,
        ),
        "calibration_bundle": (("--calibration-bundle",), "_StoreAction", "Path", None, False, None, None),
        "calibration_context": (("--calibration-context",), "_StoreAction", "Path", None, False, None, None),
        "delete_intermediates": (("--delete-intermediates",), "_StoreTrueAction", None, False, False, None, 0),
        "extra_modules": (("--extra-modules",), "_AppendAction", None, [], False, None, None),
        "fast_mode": (("--fast-mode",), "_StoreTrueAction", None, False, False, None, 0),
        "fastq1": (("--fastq1",), "_StoreAction", "str", None, False, None, None),
        "fastq2": (("--fastq2",), "_StoreAction", "str", None, False, None, None),
        "keep_intermediates": (("--keep-intermediates",), "_StoreTrueAction", None, False, False, None, 0),
        "length_annotation": (("--length-annotation",), "_StoreAction", "Path", None, False, None, None),
        "length_context": (("--length-context",), "_StoreAction", "Path", None, False, None, None),
        "length_model": (("--length-model",), "_StoreAction", "Path", None, False, None, None),
        "standard_length_model": (("--standard-length-model",), "_StoreAction", "Path", None, False, None, None),
        "estimate_vntr_length": (
            ("--estimate-vntr-length", "--no-estimate-vntr-length"),
            "BooleanOptionalAction",
            None,
            None,
            False,
            None,
            0,
        ),
        "measure_vntr_length_features": (
            ("--measure-vntr-length-features",),
            "_StoreTrueAction",
            None,
            False,
            False,
            None,
            0,
        ),
        "output_dir": (("-o", "--output-dir"), "_StoreAction", "str", None, False, None, None),
        "output_name": (("-n", "--output-name"), "_StoreAction", "str", None, False, None, None),
        "reference_assembly": (
            ("--reference-assembly",),
            "_StoreAction",
            "str",
            None,
            False,
            ASSEMBLY_CHOICES,
            None,
        ),
        "reference_fasta": (("--reference-fasta",), "_StoreAction", "Path", None, False, None, None),
        # Declared on `pipeline` as well as on `report`, and with the same default,
        # because an ordinary run reaches the report generator through `run_pipeline`
        # and never through the `report` subcommand. Two defaults would mean the same
        # run produced different artifacts depending on which subcommand rendered it.
        "report_igv": (("--report-igv",), "_StoreAction", "str", "embedded", False, REPORT_IGV_CHOICES, None),
        "resume": (("--resume",), "_StoreTrueAction", None, False, False, None, 0),
        "sample_name": (("-s", "--sample-name"), "_StoreAction", "str", None, False, None, None),
        "summary_formats": (("--summary-formats",), "_StoreAction", "str", "", False, None, None),
        "threads": (("--threads",), "_StoreAction", "positive_int", None, False, None, None),
    },
    "report": {
        "bam_file": (("--bam-file",), "_StoreAction", "Path", None, False, None, None),
        "bed_file": (("--bed-file",), "_StoreAction", "Path", None, False, None, None),
        "flanking": (("--flanking",), "_StoreAction", "int", None, False, None, None),
        "input_dir": (("--input-dir",), "_StoreAction", "Path", None, False, None, None),
        "output_dir": (("-o", "--output-dir"), "_StoreAction", "str", None, True, None, None),
        "reference_fasta": (("--reference-fasta",), "_StoreAction", "Path", None, False, None, None),
        "report_file": (("--report-file",), "_StoreAction", "str", None, False, None, None),
        # Spelled and defaulted exactly as `pipeline`'s; see the note there.
        "report_igv": (("--report-igv",), "_StoreAction", "str", "embedded", False, REPORT_IGV_CHOICES, None),
        # Spelled exactly as `pipeline`'s, deliberately: it is the same fact about
        # a run, and two spellings would be two options to remember.
        "sample_name": (("-s", "--sample-name"), "_StoreAction", "str", None, False, None, None),
        "vcf_file": (("--vcf-file",), "_StoreAction", "Path", None, False, None, None),
    },
    "cohort": {
        "input_dirs": (("-i", "--input-dirs"), "_StoreAction", None, None, False, None, "+"),
        "input_file": (("--input-file",), "_StoreAction", "Path", None, False, None, None),
        "output_dir": (("-o", "--output-dir"), "_StoreAction", "str", None, True, None, None),
        "pseudonymize_samples": (("--pseudonymize-samples",), "_StoreAction", None, None, False, None, "?"),
        "rare_allele_max_frequency": (
            ("--rare-allele-max-frequency",),
            "_StoreAction",
            "unit_fraction",
            None,
            False,
            None,
            None,
        ),
        "summary_file": (("--summary-file",), "_StoreAction", "str", None, False, None, None),
        "summary_formats": (("--summary-formats",), "_StoreAction", "str", "", False, None, None),
    },
    "install-references": {
        "aligners": (("--aligners",), "_StoreAction", None, None, False, None, "+"),
        "derive_only": (("--derive-only",), "_StoreTrueAction", None, False, False, None, 0),
        "from_source": (("--from-source",), "_StoreTrueAction", None, False, False, None, 0),
        # Note the odd one out: every other subcommand spells this -o.
        "output_dir": (("-d", "--output-dir"), "_StoreAction", "Path", None, True, None, None),
        "references": (("--references",), "_StoreAction", None, None, False, None, "+"),
        "release_spec": (("--release-spec",), "_StoreAction", "Path", None, False, None, None),
        "skip_indexing": (("--skip-indexing",), "_StoreTrueAction", None, False, False, None, 0),
        "threads": (("-t", "--threads"), "_StoreAction", "positive_int", 4, False, None, None),
    },
    "online": {
        "bam": (("--bam",), "_StoreAction", "str", None, True, None, None),
        "cohort_id": (("--cohort-id",), "_StoreAction", "str", None, False, None, None),
        "email": (("--email",), "_StoreAction", "str", None, False, None, None),
        "output_dir": (("-o", "--output-dir"), "_StoreAction", "str", None, False, None, None),
        "passphrase": (("--passphrase",), "_StoreAction", "str", None, False, None, None),
        "reference_assembly": (
            ("--reference-assembly",),
            "_StoreAction",
            "str",
            None,
            False,
            ASSEMBLY_CHOICES,
            None,
        ),
        "resume": (("--resume",), "_StoreTrueAction", None, False, False, None, 0),
        "threads": (("--threads",), "_StoreAction", "positive_int", None, False, None, None),
    },
}

#: Smallest argv that parses, per subcommand - the prefix an "unknown argument" test
#: has to satisfy before argparse gets as far as complaining about the unknown one.
MINIMAL_ARGV: dict[str, list[str]] = {
    "calibrate": [
        "calibrate",
        "fit",
        "--evidence",
        "evidence",
        "--objective",
        "lexicographic-safety-v1",
        "--output",
        "candidate",
    ],
    "pipeline": ["pipeline"],
    "report": ["report", "-o", "results"],
    "cohort": ["cohort", "-i", "a", "-o", "results"],
    "install-references": ["install-references", "-d", "refs"],
    "online": ["online", "--bam", "in.bam"],
}

#: ``calibrate``'s own options are empty on purpose - everything it accepts lives on a
#: nested subparser, one per operation. ``SUBCOMMAND_CONTRACT["calibrate"]`` therefore
#: cannot see a single one of them, and adding an operation or an option to an existing
#: operation used to be invisible to this file. This second table closes that hole with
#: the same contract in the same shape: one row per option, compared wholesale.
CALIBRATE_OPERATION_CONTRACT: dict[str, dict[str, ParserRow]] = {
    "assess": {
        "exposure_ledger": (("--exposure-ledger",), "_StoreAction", "Path", None, True, None, None),
        "intake": (("--intake",), "_StoreAction", "Path", None, True, None, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "profile": (("--profile",), "_StoreAction", "Path", None, True, None, None),
        "runs": (("--runs",), "_StoreAction", "Path", None, True, None, None),
        "target": (("--target",), "_StoreAction", None, None, True, ("callers", "length"), None),
    },
    "cohort": {
        "caller_policies": (("--caller-policies",), "_StoreAction", "Path", None, False, None, None),
        "caller_runs": (("--caller-runs",), "_StoreAction", "Path", None, False, None, None),
        "count_convention": (
            ("--count-convention",),
            "_StoreAction",
            None,
            "source-reported",
            False,
            ("source-reported", "complete", "canonical-only"),
            None,
        ),
        "folds": (("--folds",), "_StoreAction", "int", 5, False, None, None),
        "manifest": (("--manifest",), "_StoreAction", "Path", None, True, None, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "reference": (("--reference",), "_StoreAction", "Path", None, False, None, None),
        "seed": (("--seed",), "_StoreAction", "int", 20260915, False, None, None),
        "target": (("--target",), "_StoreAction", None, "auto", False, ("auto", "length", "callers", "both"), None),
    },
    "evaluate": {
        "authority": (("--authority",), "_StoreAction", "Path", None, False, None, None),
        "custody": (("--custody",), "_StoreAction", "Path", None, False, None, None),
        "evidence": (("--evidence",), "_StoreAction", "Path", None, True, None, None),
        "exposure_ledger": (("--exposure-ledger",), "_StoreAction", "Path", None, False, None, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "profile": (("--profile",), "_StoreAction", "Path", None, True, None, None),
        "target": (("--target",), "_StoreAction", None, "dominance", False, DOMINANCE_TARGETS, None),
        "validation": (("--validation",), "_StoreAction", "Path", None, False, None, None),
    },
    "export": {
        "authority": (("--authority",), "_StoreAction", "Path", None, True, None, None),
        "completion": (("--completion",), "_StoreAction", "Path", None, True, None, None),
        "evaluation": (("--evaluation",), "_StoreAction", "Path", None, True, None, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "profile": (("--profile",), "_StoreAction", "Path", None, True, None, None),
        "target": (("--target",), "_StoreAction", None, None, True, ("callers", "length"), None),
        "validation": (("--validation",), "_StoreAction", "Path", None, True, None, None),
    },
    "extract": {
        "length_annotation": (("--length-annotation",), "_StoreAction", "Path", None, False, None, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "partitions": (("--partitions",), "_StoreAction", "Path", None, False, None, None),
        "runs": (("--runs",), "_StoreAction", "Path", None, True, None, None),
        "sources": (("--sources",), "_StoreAction", "Path", None, False, None, None),
        "study": (("--study",), "_StoreAction", "Path", None, False, None, None),
        "target": (("--target",), "_StoreAction", None, "dominance", False, DOMINANCE_TARGETS, None),
        "truth": (("--truth",), "_StoreAction", "Path", None, False, None, None),
    },
    "fit": {
        "advntr_executable": (("--advntr-executable",), "_StoreAction", "Path", None, False, None, None),
        "evidence": (("--evidence",), "_StoreAction", "Path", None, True, None, None),
        "exposure_ledger": (("--exposure-ledger",), "_StoreAction", "Path", None, False, None, None),
        "objective": (
            ("--objective",),
            "_StoreAction",
            None,
            None,
            True,
            ("lexicographic-safety-v1", "caller-safety-v1", "length-total-v1"),
            None,
        ),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "target": (("--target",), "_StoreAction", None, "dominance", False, DOMINANCE_TARGETS, None),
    },
    "intake": {
        "cram_references": (("--cram-references",), "_StoreAction", "Path", None, False, None, None),
        "manifest": (("--manifest",), "_StoreAction", "Path", None, True, None, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "preprocessing_priority": (("--preprocessing-priority",), "_StoreAction", None, None, True, None, "+"),
        "temporary_directory": (("--temporary-directory",), "_StoreAction", "Path", None, False, None, None),
    },
    "optimize": {
        "advntr_executable": (("--advntr-executable",), "_StoreAction", "Path", None, False, None, None),
        "axes": (("--axis",), "_AppendAction", None, None, False, CUTOFF_AXIS_CHOICES, None),
        "caller": (("--caller",), "_StoreAction", None, "kestrel", False, ("kestrel", "advntr", "both"), None),
        "captures": (("--captures",), "_StoreAction", "Path", None, True, None, None),
        "folds": (("--folds",), "_StoreAction", "int", 5, False, None, None),
        "manifest": (("--manifest",), "_StoreAction", "Path", None, True, None, None),
        "max_breakpoints": (("--max-breakpoints",), "_StoreAction", "int", None, False, None, None),
        "min_sensitivity": (("--min-sensitivity",), "_StoreAction", "float", None, False, None, None),
        "min_specificity": (("--min-specificity",), "_StoreAction", "float", None, False, None, None),
        "objective": (("--objective",), "_StoreAction", None, None, True, CUTOFF_OBJECTIVE_CHOICES, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "seed": (("--seed",), "_StoreAction", "int", 20260915, False, None, None),
        "workers": (("--workers",), "_StoreAction", "int", 1, False, None, None),
    },
    "validate": {
        "authority": (("--authority",), "_StoreAction", "Path", None, False, None, None),
        "custody": (("--custody",), "_StoreAction", "Path", None, False, None, None),
        "evidence": (("--evidence",), "_StoreAction", "Path", None, True, None, None),
        "exposure_ledger": (("--exposure-ledger",), "_StoreAction", "Path", None, False, None, None),
        "output": (("--output",), "_StoreAction", "Path", None, True, None, None),
        "profile": (("--profile",), "_StoreAction", "Path", None, True, None, None),
        "target": (("--target",), "_StoreAction", None, "dominance", False, DOMINANCE_TARGETS, None),
        "validation": (("--validation",), "_StoreAction", "Path", None, False, None, None),
    },
}

#: Smallest argv that parses, per calibration operation, mirroring ``MINIMAL_ARGV``.
MINIMAL_CALIBRATE_ARGV: dict[str, list[str]] = {
    "assess": ["--target", "callers", "--profile", "p", "--intake", "i", "--runs", "r", "--exposure-ledger", "e"],
    "cohort": ["--manifest", "m"],
    "evaluate": ["--profile", "p", "--evidence", "e"],
    "export": [
        "--target",
        "callers",
        "--profile",
        "p",
        "--validation",
        "v",
        "--evaluation",
        "x",
        "--authority",
        "a",
        "--completion",
        "c",
    ],
    "extract": ["--runs", "r"],
    "fit": ["--evidence", "e", "--objective", "lexicographic-safety-v1"],
    "intake": ["--manifest", "m", "--preprocessing-priority", "id"],
    "optimize": ["--manifest", "m", "--captures", "c", "--objective", "youden-j"],
    "validate": ["--profile", "p", "--evidence", "e"],
}


def _describe(action: argparse.Action) -> ParserRow:
    """Reduce an argparse action to the tuple the contract table states.

    Args:
        action: The action to describe.

    Returns:
        ParserRow: Option strings, action class name, type name, default,
        required flag, choices and nargs.
    """
    return (
        tuple(action.option_strings),
        type(action).__name__,
        getattr(action.type, "__name__", None),
        action.default,
        action.required,
        tuple(action.choices) if action.choices else None,
        action.nargs,
    )


def _options(parser: argparse.ArgumentParser) -> dict[str, ParserRow]:
    """Describe every option on ``parser``, minus ``-h`` and the subparsers action.

    Args:
        parser: The parser to inspect.

    Returns:
        dict[str, ParserRow]: dest -> contract row.
    """
    return {
        action.dest: _describe(action)
        for action in parser._actions
        if not isinstance(action, argparse._SubParsersAction | argparse._HelpAction)
    }


def _subparsers(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    """Return the registered subparsers by name.

    Args:
        parser: The top-level parser.

    Returns:
        dict[str, argparse.ArgumentParser]: Subcommand name -> its parser.
    """
    actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    assert len(actions) == 1, f"expected exactly one subparsers action, found {len(actions)}"
    return dict(actions[0].choices)


def test_the_contract_table_covers_every_subcommand() -> None:
    """A new subcommand must arrive with its own contract row, not silently."""
    assert set(_subparsers(build_parser())) == set(SUBCOMMAND_CONTRACT)


def _calibrate_operations() -> dict[str, argparse.ArgumentParser]:
    """Return the ``calibrate`` operation parsers by name.

    Returns:
        dict[str, argparse.ArgumentParser]: Operation name -> its parser.
    """
    calibrate = _subparsers(build_parser())["calibrate"]
    actions = [action for action in calibrate._actions if isinstance(action, argparse._SubParsersAction)]
    assert len(actions) == 1, f"expected exactly one calibrate operation action, found {len(actions)}"
    return dict(actions[0].choices)


def test_the_calibrate_operation_table_covers_every_operation() -> None:
    """A new calibration operation must arrive with its own contract rows."""
    assert set(_calibrate_operations()) == set(CALIBRATE_OPERATION_CONTRACT)
    assert set(MINIMAL_CALIBRATE_ARGV) == set(CALIBRATE_OPERATION_CONTRACT)


@pytest.mark.parametrize("operation", sorted(CALIBRATE_OPERATION_CONTRACT))
def test_every_calibrate_operation_option_matches_the_contract(operation: str) -> None:
    """Flags, types, defaults and required-ness of one operation, as one table.

    Args:
        operation: The calibration operation under test.
    """
    actual = _options(_calibrate_operations()[operation])
    expected = CALIBRATE_OPERATION_CONTRACT[operation]
    assert actual == expected, (
        f"calibrate {operation}'s options drifted from the contract; "
        f"only in parser: {sorted(set(actual) - set(expected))}; "
        f"only in contract: {sorted(set(expected) - set(actual))}"
    )


@pytest.mark.parametrize("operation", sorted(CALIBRATE_OPERATION_CONTRACT))
def test_calibrate_operation_defaults_survive_a_minimal_parse(operation: str) -> None:
    """The declared defaults are what a minimal invocation actually produces.

    Args:
        operation: The calibration operation under test.
    """
    argv = ["calibrate", operation, *MINIMAL_CALIBRATE_ARGV[operation], "--output", "out"]
    args = build_parser().parse_args(argv)
    assert args.calibration_operation == operation
    for dest, row in CALIBRATE_OPERATION_CONTRACT[operation].items():
        if any(flag in argv for flag in row[0]):
            continue  # supplied on the command line, so its default is not observable
        assert getattr(args, dest) == row[3], f"calibrate {operation} --{dest} defaulted to {getattr(args, dest)!r}"


@pytest.mark.parametrize("operation", sorted(CALIBRATE_OPERATION_CONTRACT))
def test_every_calibrate_operation_rejects_an_unknown_argument(operation: str) -> None:
    """An unknown flag on an operation is a usage error, never an absorbed positional.

    Args:
        operation: The calibration operation under test.
    """
    argv = ["calibrate", operation, *MINIMAL_CALIBRATE_ARGV[operation], "--output", "out", "--not-a-real-option"]
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(argv)
    assert excinfo.value.code == 2


@pytest.mark.parametrize("operation", sorted(CALIBRATE_OPERATION_CONTRACT))
def test_every_calibrate_operation_requires_its_output(operation: str) -> None:
    """Dropping ``required=True`` on ``--output`` would send ``None`` to ``atomic_output``.

    Args:
        operation: The calibration operation under test.
    """
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(["calibrate", operation, *MINIMAL_CALIBRATE_ARGV[operation]])
    assert excinfo.value.code == 2


def test_the_cutoff_axis_choices_are_exactly_the_declared_axes() -> None:
    """An axis the CLI cannot name is unreachable; a name the axes reject fails deep."""
    from vntyper.scripts import calibration_cutoff_axes as axes
    from vntyper.scripts.cli_calibration_parser import CUTOFF_AXIS_CHOICES as declared

    constants = {
        axes.DEPTH_FLOOR_LINKED,
        axes.GG_GATE_INDEPENDENT,
        axes.DEPTH_SCORE_HIGH,
        axes.ALT_DEPTH_BAND,
        axes.ACTIVE_REGION,
        axes.ADVNTR_CUTOFF,
        axes.ADVNTR_MIN_SUPPORT,
    }
    assert set(CUTOFF_AXIS_CHOICES) == constants
    assert declared == CUTOFF_AXIS_CHOICES


def test_the_cutoff_objective_choices_are_exactly_the_supported_objectives() -> None:
    """The CLI offers every selectable objective, and nothing the selector refuses."""
    from vntyper.scripts.calibration_cutoff_selection import OBJECTIVES

    assert CUTOFF_OBJECTIVE_CHOICES == OBJECTIVES


def test_the_top_level_options_match_the_contract() -> None:
    """``-l``, ``-f``, ``-v`` and ``--config-path``, exactly as ``cli.main`` reads them."""
    assert _options(build_parser()) == TOP_LEVEL_CONTRACT


@pytest.mark.parametrize("command", sorted(SUBCOMMAND_CONTRACT))
def test_every_subcommand_option_matches_the_contract(command: str) -> None:
    """Flags, types, defaults and required-ness, compared as one table.

    Args:
        command: The subcommand under test.
    """
    actual = _options(_subparsers(build_parser())[command])
    expected = SUBCOMMAND_CONTRACT[command]
    assert actual == expected, (
        f"{command}'s options drifted from the contract; "
        f"only in parser: {sorted(set(actual) - set(expected))}; "
        f"only in contract: {sorted(set(expected) - set(actual))}"
    )


@pytest.mark.parametrize("command", sorted(SUBCOMMAND_CONTRACT))
def test_defaults_survive_a_minimal_parse(command: str) -> None:
    """The declared defaults must be what a minimal invocation actually produces.

    The contract table reads ``action.default``; this asserts argparse then puts
    that value on the namespace, which is what every handler goes on to read.

    Args:
        command: The subcommand under test.
    """
    args = build_parser().parse_args(MINIMAL_ARGV[command])
    for dest, row in SUBCOMMAND_CONTRACT[command].items():
        default = row[3]
        if any(flag in MINIMAL_ARGV[command] for flag in row[0]):
            continue  # supplied on the command line, so its default is not observable
        assert getattr(args, dest) == default, f"{command} --{dest} defaulted to {getattr(args, dest)!r}"


@pytest.mark.parametrize("command", sorted(SUBCOMMAND_CONTRACT))
def test_every_subcommand_rejects_an_unknown_argument(command: str) -> None:
    """An unknown flag is a usage error (exit 2), never an absorbed positional.

    Args:
        command: The subcommand under test.
    """
    argv = [*MINIMAL_ARGV[command], "--definitely-not-a-real-option"]
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(argv)
    assert excinfo.value.code == 2


@pytest.mark.parametrize("command", ["report", "cohort", "install-references", "online"])
def test_the_required_options_are_enforced(command: str) -> None:
    """Dropping a ``required=True`` would turn a usage error into a ``None``.

    Args:
        command: The subcommand under test.
    """
    required = [dest for dest, row in SUBCOMMAND_CONTRACT[command].items() if row[4]]
    assert required, f"{command} declares no required option; this test would be vacuous"
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args([command])
    assert excinfo.value.code == 2


@pytest.mark.parametrize("command", ["pipeline", "online"])
def test_reference_assembly_offers_exactly_what_the_registry_knows(command: str) -> None:
    """The CLI's choices and the registry's assemblies are the same set.

    A registry entry the CLI cannot name is unreachable; a CLI choice the registry
    does not know raises deep inside the run instead of at parse time.

    Args:
        command: The subcommand under test.
    """
    assemblies = list_assemblies()
    assert assemblies, "the registry listed no assemblies; this test would be vacuous"
    action = next(
        action for action in _subparsers(build_parser())[command]._actions if action.dest == "reference_assembly"
    )
    assert action.choices is not None, f"{command} --reference-assembly stopped constraining its values"
    assert sorted(action.choices) == sorted(assemblies)


@pytest.mark.parametrize("command", ["pipeline", "online"])
def test_an_unknown_reference_assembly_is_a_usage_error(command: str) -> None:
    """``choices`` must actually bite, not merely decorate the help text.

    Args:
        command: The subcommand under test.
    """
    argv = [*MINIMAL_ARGV[command], "--reference-assembly", "hg17"]
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(argv)
    assert excinfo.value.code == 2


@pytest.mark.parametrize("command", sorted(SUBCOMMAND_CONTRACT))
@pytest.mark.parametrize("option", ["--log-level", "--log-file", "--config-path"])
def test_the_global_options_are_rejected_after_the_subcommand(command: str, option: str) -> None:
    """The global options are top-level only, and the docstring now says so.

    ``build_parser`` never passes ``parent_parser`` to ``add_parser``, so
    ``vntyper pipeline --log-level DEBUG`` has always been a usage error. This
    pins the real behaviour so that adding the parent later is a visible change,
    not an accident.

    Args:
        command: The subcommand under test.
        option: The global option being tried in the wrong position.
    """
    value = "DEBUG" if option == "--log-level" else "x"
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args([*MINIMAL_ARGV[command], option, value])
    assert excinfo.value.code == 2


@pytest.mark.parametrize("command", sorted(SUBCOMMAND_CONTRACT))
@pytest.mark.parametrize("option", ["--log-level", "--log-file", "--config-path"])
def test_the_global_options_are_accepted_before_the_subcommand(command: str, option: str) -> None:
    """The position that does work, pinned alongside the one that does not.

    Args:
        command: The subcommand under test.
        option: The global option.
    """
    value = "DEBUG" if option == "--log-level" else "x"
    args = build_parser().parse_args([option, value, *MINIMAL_ARGV[command]])
    assert args.command == command
    assert getattr(args, option.lstrip("-").replace("-", "_")) is not None
