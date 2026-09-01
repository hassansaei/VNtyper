"""The CLI argument-parsing seam: ``build_parser()`` and ``main(argv)``.

``vntyper/cli.py`` used to build its ``argparse`` tree inline inside ``main()`` and call
``parser.parse_args()`` with no argument, which made the CLI untestable in two separate
ways: there was no way to obtain a parser without running the program, and no way to
drive ``main`` with anything other than the real ``sys.argv``. ``cli.py`` sat at 0%
coverage as a result.

``build_parser()`` now lives in ``vntyper.scripts.cli_parser`` and ``main`` takes an
explicit ``argv``. These tests pin that seam: the parser's shape (every subcommand, the
global options, the mutually exclusive groups, the defaults) and ``main``'s argv
handling, including that ``argv=None`` still reads ``sys.argv`` so the console entry
point is unaffected.

This file deliberately stops at the parser. It does not exercise what the subcommands
go on to do - ``vntyper report`` in particular was known to be broken (AGENTS.md
trap 11, fixed in #179) and is owned elsewhere. The point here is only that those
paths are *reachable* from a test.
"""

import argparse
import sys

import pytest

from vntyper.cli import main
from vntyper.scripts.cli_parser import build_parser
from vntyper.version import __version__ as VERSION

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit

EXPECTED_SUBCOMMANDS = {"pipeline", "report", "cohort", "install-references", "online", "calibrate"}


def _subcommand_choices(parser: argparse.ArgumentParser) -> set[str]:
    """
    Return the names of every subcommand registered on ``parser``.

    Args:
        parser (argparse.ArgumentParser): The parser to inspect.

    Returns:
        set[str]: The registered subcommand names.
    """
    # argparse exposes no public accessor for registered subparsers, so this reaches
    # into `_actions`. The length assertion below is the guard against that private
    # detail changing shape and quietly turning the caller's assertion vacuous.
    subparser_actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    assert len(subparser_actions) == 1, f"expected exactly one subparsers action, found {len(subparser_actions)}"
    return set(subparser_actions[0].choices)


def test_build_parser_returns_an_argument_parser():
    """The seam every CLI test needs: a parser without running the program."""
    parser = build_parser()
    assert isinstance(parser, argparse.ArgumentParser)


def test_build_parser_returns_a_fresh_parser_on_every_call():
    """
    Each call must build a new parser rather than hand out a module-level singleton.

    argparse parsers carry mutable state, and a shared one would let a test that
    provokes an error leak that state into the next test.
    """
    assert build_parser() is not build_parser()


def test_every_subcommand_is_registered():
    """All six documented subcommands must survive the extraction."""
    choices = _subcommand_choices(build_parser())
    assert choices == EXPECTED_SUBCOMMANDS


def test_the_global_options_are_available_before_the_subcommand():
    """``-l``, ``-f`` and ``--config-path`` come from the shared parent parser."""
    args = build_parser().parse_args(["--log-level", "DEBUG", "--log-file", "run.log", "pipeline"])
    assert args.command == "pipeline"
    assert args.log_level == "DEBUG"
    assert args.log_file == "run.log"
    assert args.config_path is None


def test_the_pipeline_subcommand_keeps_its_defaults():
    """A bare ``pipeline`` invocation must parse with the defaults ``main`` relies on."""
    args = build_parser().parse_args(["pipeline"])
    assert args.command == "pipeline"
    assert args.output_dir is None
    assert args.threads is None
    assert args.extra_modules == []
    assert args.summary_formats == ""
    assert args.fast_mode is False


def test_the_report_subcommand_parses_an_output_directory():
    """
    The seam for the ``report`` tests: ``report -o DIR`` must reach ``args``.

    What ``report`` then does with those arguments is covered in
    ``tests/unit/test_cli_report.py``; this only pins that the arguments parse.
    """
    args = build_parser().parse_args(["report", "-o", "results"])
    assert args.command == "report"
    assert args.output_dir == "results"


def test_the_pipeline_region_options_stay_mutually_exclusive():
    """``--custom-regions`` and ``--bed-file`` must not be accepted together."""
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(["pipeline", "--custom-regions", "chr1:1-2", "--bed-file", "regions.bed"])
    assert excinfo.value.code == 2


def test_the_cohort_subcommand_still_requires_an_input_source():
    """``cohort`` needs one of ``--input-dirs`` / ``--input-file``; the group is required."""
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(["cohort", "-o", "results"])
    assert excinfo.value.code == 2


THREAD_SUBCOMMAND_ARGV = [
    ["pipeline"],
    ["install-references", "--output-dir", "references"],
    ["online", "--bam", "sample.bam"],
]


@pytest.mark.parametrize("prefix", THREAD_SUBCOMMAND_ARGV)
@pytest.mark.parametrize("value", ["0", "-1", "not-a-number"])
def test_every_threads_option_refuses_invalid_values(prefix: list[str], value: str) -> None:
    """Each thread-bearing subcommand reports malformed values as usage errors.

    Args:
        prefix: Minimal valid invocation for the subcommand under test.
        value: Non-positive or non-integer token under test.
    """
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args([*prefix, "--threads", value])

    assert excinfo.value.code == 2


@pytest.mark.parametrize("prefix", THREAD_SUBCOMMAND_ARGV)
@pytest.mark.parametrize("value", [1, 64, 999999])
def test_every_threads_option_preserves_positive_values(prefix: list[str], value: int) -> None:
    """Every previously accepted positive thread count remains unchanged.

    Args:
        prefix: Minimal valid invocation for the subcommand under test.
        value: Positive value under test, including one above the web bound.
    """
    args = build_parser().parse_args([*prefix, "--threads", str(value)])

    assert args.threads == value


def test_main_with_an_empty_argv_prints_help_and_exits_zero(capsys):
    """No subcommand is not an error: print the help and exit 0."""
    with pytest.raises(SystemExit) as excinfo:
        main([])
    assert excinfo.value.code == 0
    assert "usage:" in capsys.readouterr().out


def test_main_help_exits_zero(capsys):
    """``--help`` is handled by argparse and lists the subcommands."""
    with pytest.raises(SystemExit) as excinfo:
        main(["--help"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert "usage:" in out
    for subcommand in EXPECTED_SUBCOMMANDS:
        assert subcommand in out


def test_main_version_exits_zero_and_prints_the_package_version(capsys):
    """``--version`` reports the single authoritative version."""
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert VERSION in capsys.readouterr().out


def test_main_rejects_an_unknown_option_with_exit_code_two():
    """An unparsable argv exits 2, argparse's usage-error code."""
    with pytest.raises(SystemExit) as excinfo:
        main(["--definitely-not-a-real-option"])
    assert excinfo.value.code == 2


def test_main_reads_sys_argv_when_argv_is_none(monkeypatch, capsys):
    """
    ``argv=None`` must keep reading ``sys.argv[1:]``, so the console script still works.

    Two invocations with different ``sys.argv`` must produce different outcomes;
    otherwise a ``main`` that ignored ``sys.argv`` entirely would pass this test.
    """
    monkeypatch.setattr(sys, "argv", ["vntyper"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    assert "usage:" in capsys.readouterr().out

    monkeypatch.setattr(sys, "argv", ["vntyper", "--definitely-not-a-real-option"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 2
