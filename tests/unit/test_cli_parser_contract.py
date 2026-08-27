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

SUBCOMMAND_CONTRACT: dict[str, dict[str, ParserRow]] = {
    "pipeline": {
        "advntr_max_coverage": (("--advntr-max-coverage",), "_StoreAction", "int", None, False, None, None),
        "archive_format": (("--archive-format",), "_StoreAction", "str", None, False, ("zip", "tar.gz"), None),
        "archive_results": (("--archive-results",), "_StoreTrueAction", None, False, False, None, 0),
        "bam": (("--bam",), "_StoreAction", "str", None, False, None, None),
        "bed_file": (("--bed-file",), "_StoreAction", "Path", None, False, None, None),
        "cram": (("--cram",), "_StoreAction", "str", None, False, None, None),
        "custom_regions": (("--custom-regions",), "_StoreAction", "str", None, False, None, None),
        "delete_intermediates": (("--delete-intermediates",), "_StoreTrueAction", None, False, False, None, 0),
        "extra_modules": (("--extra-modules",), "_AppendAction", None, [], False, None, None),
        "fast_mode": (("--fast-mode",), "_StoreTrueAction", None, False, False, None, 0),
        "fastq1": (("--fastq1",), "_StoreAction", "str", None, False, None, None),
        "fastq2": (("--fastq2",), "_StoreAction", "str", None, False, None, None),
        "keep_intermediates": (("--keep-intermediates",), "_StoreTrueAction", None, False, False, None, 0),
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
    "pipeline": ["pipeline"],
    "report": ["report", "-o", "results"],
    "cohort": ["cohort", "-i", "a", "-o", "results"],
    "install-references": ["install-references", "-d", "refs"],
    "online": ["online", "--bam", "in.bam"],
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
