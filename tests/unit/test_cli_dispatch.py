"""``main()`` is parse-then-dispatch: one handler per subcommand, nothing inline.

Every subcommand body used to live inside ``vntyper.cli.main``. Two people changing
two different subcommands were therefore changing the same function, and no
subcommand could be exercised without exercising all the argument plumbing around
it. The bodies now live in :mod:`vntyper.scripts.cli_handlers` (and, for ``report``
alone, :mod:`vntyper.scripts.cli_report`).

These tests pin the dispatch seam that split creates:

* ``HANDLERS`` maps each registered subcommand to the extracted function -- so a
  handler cannot be silently orphaned by a rename;
* ``main([...])`` reaches that handler for every subcommand, with the parsed
  namespace and the resolved log level/log file;
* logging is configured exactly once, in ``cli.py``, *before* any handler runs.

They deliberately stop at the handler boundary: what each handler goes on to do is
covered where that handler lives. ``vntyper report`` in particular was known broken
(AGENTS.md trap 11, fixed in #179) and is owned elsewhere -- the point here is
that it is reachable, and separately editable, from a test.
"""

import argparse
import logging

import pytest

from vntyper import cli
from vntyper.scripts import cli_calibrate, cli_handlers, cli_report
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


@pytest.fixture(autouse=True)
def _restore_root_logger():
    """Undo ``setup_logging``'s reconfiguration of the root logger.

    ``setup_logging`` clears every root handler, which would take pytest's own
    logging handlers with it and leave later tests unable to capture anything.
    """
    root = logging.getLogger()
    saved_handlers = list(root.handlers)
    saved_level = root.level
    try:
        yield
    finally:
        root.handlers[:] = saved_handlers
        root.setLevel(saved_level)


def _subcommand_choices() -> set[str]:
    """Return the subcommand names the parser registers.

    Returns:
        set[str]: The registered subcommand names.
    """
    parser = build_parser()
    actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    assert len(actions) == 1, f"expected exactly one subparsers action, found {len(actions)}"
    return set(actions[0].choices)


# Smallest argv that parses for each subcommand. Directory-valued options are
# formatted in with tmp_path so that nothing is written outside the test.
MINIMAL_ARGV: dict[str, list[str]] = {
    "pipeline": ["pipeline", "-o", "{tmp}/outputs", "--bam", "{tmp}/inputs/in.bam"],
    "report": ["report", "-o", "{tmp}"],
    "cohort": ["cohort", "-i", "{tmp}", "-o", "{tmp}"],
    "install-references": ["install-references", "-d", "{tmp}"],
    "online": ["online", "--bam", "{tmp}/in.bam", "-o", "{tmp}"],
    "calibrate": [
        "calibrate",
        "fit",
        "--evidence",
        "{tmp}/evidence",
        "--objective",
        "lexicographic-safety-v1",
        "--output",
        "{tmp}/candidate",
    ],
}

EXPECTED_HANDLERS = {
    "pipeline": cli_handlers.handle_pipeline,
    "report": cli_report.handle_report,
    "cohort": cli_handlers.handle_cohort,
    "install-references": cli_handlers.handle_install_references,
    "online": cli_handlers.handle_online,
    "calibrate": cli_calibrate.handle_calibrate,
}


def test_every_subcommand_has_a_handler() -> None:
    """A subcommand with no entry in HANDLERS falls through to the error path."""
    choices = _subcommand_choices()
    assert choices, "the parser registered no subcommands; this test would be vacuous"
    assert choices == set(cli.HANDLERS), (
        f"subcommands without a handler: {sorted(choices - set(cli.HANDLERS))}; "
        f"handlers without a subcommand: {sorted(set(cli.HANDLERS) - choices)}"
    )


def test_handlers_point_at_the_extracted_functions() -> None:
    """Pins where each handler lives, which is what makes them separately ownable."""
    assert cli.HANDLERS == EXPECTED_HANDLERS


def test_the_report_handler_lives_in_its_own_module() -> None:
    """``report`` is owned separately from every other subcommand, so it is alone."""
    assert cli_report.handle_report.__module__ == "vntyper.scripts.cli_report"
    assert not hasattr(cli_handlers, "handle_report")


def test_main_body_no_longer_contains_a_subcommand_body() -> None:
    """``main`` is parse-then-dispatch; the stage calls belong to the handlers."""
    import inspect

    source = inspect.getsource(cli.main)
    assert "HANDLERS" in source, "the dispatch table vanished; the rest of this assertion is vacuous"
    for leaked in ("run_pipeline(", "aggregate_cohort(", "run_online_mode(", "generate_summary_report("):
        assert leaked not in source, f"{leaked} is still inline in main()"


@pytest.mark.parametrize("command", sorted(MINIMAL_ARGV))
def test_main_dispatches_each_subcommand_to_its_handler(command, tmp_path, monkeypatch) -> None:
    """Drive the real ``main([...])`` and assert the right handler is reached."""
    calls = []

    def _record(args, config, parser, log_level_value, log_file_str):
        calls.append(
            {
                "args": args,
                "config": config,
                "parser": parser,
                "log_level_value": log_level_value,
                "log_file_str": log_file_str,
            }
        )

    monkeypatch.setitem(cli.HANDLERS, command, _record)
    argv = [part.format(tmp=tmp_path) for part in MINIMAL_ARGV[command]]

    cli.main(argv)

    assert len(calls) == 1, f"{command} reached its handler {len(calls)} times, expected once"
    call = calls[0]
    assert call["args"].command == command
    assert isinstance(call["config"], dict) and call["config"], "the handler got no configuration"
    assert isinstance(call["parser"], argparse.ArgumentParser)
    assert isinstance(call["log_level_value"], int)


def test_only_the_dispatched_handler_runs(tmp_path, monkeypatch) -> None:
    """A subcommand must not fall through into another subcommand's handler."""
    reached = []
    for name in cli.HANDLERS:
        monkeypatch.setitem(
            cli.HANDLERS,
            name,
            lambda args, config, parser, log_level_value, log_file_str, _name=name: reached.append(_name),
        )

    cli.main(["cohort", "-i", str(tmp_path), "-o", str(tmp_path)])

    assert reached == ["cohort"]


def test_the_handler_receives_the_resolved_log_file(tmp_path, monkeypatch) -> None:
    """``cli.py`` resolves the log file; the handler only forwards it onward."""
    seen = {}

    def _record(args, config, parser, log_level_value, log_file_str):
        seen["log_file_str"] = log_file_str
        seen["log_level_value"] = log_level_value

    monkeypatch.setitem(cli.HANDLERS, "pipeline", _record)
    output_dir = tmp_path / "outputs"
    cli.main(
        [
            "--log-level",
            "WARNING",
            "pipeline",
            "-o",
            str(output_dir),
            "--bam",
            str(tmp_path / "inputs" / "in.bam"),
        ]
    )

    assert seen["log_file_str"] == str(output_dir / "pipeline.log")
    assert seen["log_level_value"] == logging.WARNING


def test_logging_is_configured_before_the_handler_runs(tmp_path, monkeypatch) -> None:
    """``cli.py`` is the sole place logging is set up, and it happens first."""
    order = []

    def _fake_setup_logging(log_level=logging.INFO, log_file=None):
        order.append("setup_logging")

    monkeypatch.setattr(cli, "setup_logging", _fake_setup_logging)
    monkeypatch.setitem(
        cli.HANDLERS,
        "cohort",
        lambda args, config, parser, log_level_value, log_file_str: order.append("handler"),
    )

    cli.main(["cohort", "-i", str(tmp_path), "-o", str(tmp_path)])

    assert order == ["setup_logging", "handler"]


def test_a_handler_value_error_exits_one_instead_of_escaping_as_a_traceback(tmp_path, monkeypatch, caplog) -> None:
    """Item 10: a handler-level `ValueError` (AGENTS.md's `logger.error(msg)` then
    `raise ValueError(msg)` convention, no custom exception classes) must get the same
    clean `sys.exit(1)` every other CLI-level validation error already gets here, rather
    than escape `main()` uncaught. `_resolve_bwa_reference`'s new fail-closed error on a
    present-but-uninstalled BWA reference is what surfaced this, but the fix and this
    test exercise the real end-to-end path: a genuine FASTQ run through the real
    `handle_pipeline`, with nothing configured for a required BWA reference.
    """
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)
    monkeypatch.setattr(cli, "load_config", lambda config_path=None: {"reference_data": {}})

    fastq1 = tmp_path / "r1.fastq.gz"
    fastq1.write_bytes(b"fastq")

    with (
        caplog.at_level(logging.ERROR),
        pytest.raises(SystemExit) as excinfo,
    ):
        cli.main(["pipeline", "--fastq1", str(fastq1), "-o", str(tmp_path / "out")])

    assert excinfo.value.code == 1
    assert any(
        record.levelno >= logging.ERROR and "No BWA reference configured" in record.message for record in caplog.records
    ), "the handler's own diagnostic must not be swallowed by the clean exit"


def _called_function_names(module) -> set[str]:
    """Return every function name called anywhere in ``module``'s source.

    Uses the AST rather than a substring scan so that a docstring naming a
    function is not mistaken for a call to it.

    Args:
        module: The imported module to parse.

    Returns:
        set[str]: The called names, attribute calls reduced to the attribute.
    """
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(module))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name):
                names.add(func.id)
            elif isinstance(func, ast.Attribute):
                names.add(func.attr)
    return names


def test_no_handler_module_configures_logging() -> None:
    """Handlers must not re-run ``setup_logging``; ``cli.py`` owns it (AGENTS.md)."""
    assert "setup_logging" in _called_function_names(cli), "cli.py stopped configuring logging; nothing else may"

    for module in (cli_handlers, cli_report, cli_calibrate):
        called = _called_function_names(module)
        assert called, f"parsed no calls out of {module.__name__}; this assertion would be vacuous"
        assert "setup_logging" not in called, f"{module.__name__} configures logging; only cli.py may"
        assert "basicConfig" not in called, f"{module.__name__} calls basicConfig; only cli.py configures logging"


def test_an_unroutable_command_exits_one(tmp_path, monkeypatch, caplog) -> None:
    """The fall-through path is still reachable once a handler goes missing."""
    monkeypatch.setattr(cli, "HANDLERS", {})
    # setup_logging() clears every root handler, caplog's included, so a caplog
    # assertion made after main() has run sees nothing at all. Stub it out: this
    # test is about dispatch fall-through, not about logging configuration.
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)

    with caplog.at_level(logging.ERROR, logger="vntyper.cli"), pytest.raises(SystemExit) as excinfo:
        cli.main(["cohort", "-i", str(tmp_path), "-o", str(tmp_path)])

    assert excinfo.value.code == 1
    assert any(record.levelno >= logging.ERROR and "Unknown command" in record.message for record in caplog.records)


def test_an_unhonourable_output_name_exits_one_before_anything_is_created(tmp_path, monkeypatch, caplog) -> None:
    """A flag the pipeline cannot honour is refused before the run leaves a trace.

    ``--output-name`` was validated inside the handler, which runs after logging
    has been configured -- and the log file for a ``pipeline`` run is
    ``<output_dir>/pipeline.log``, so the output directory was created, written to,
    and only then did the run die on input that was never acceptable. It died with
    an unhandled ``ValueError`` as well: a traceback for a usage error.
    """
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)
    called = []
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *a, **k: called.append(True))

    output_dir = tmp_path / "new-dir"
    with caplog.at_level(logging.CRITICAL, logger="vntyper.cli"), pytest.raises(SystemExit) as excinfo:
        cli.main(["pipeline", "--bam", str(tmp_path / "x.bam"), "-o", str(output_dir), "--output-name", "sample"])

    assert excinfo.value.code == 1
    assert not called, "the handler ran on input the CLI had already decided to refuse"
    assert not output_dir.exists(), "the refused run created its output directory"
    assert any("--output-name" in record.message for record in caplog.records)


def test_the_default_output_name_still_dispatches(tmp_path, monkeypatch) -> None:
    """The check is a guard on one flag, not a new way for ordinary runs to fail."""
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)
    called = []
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *a, **k: called.append(True))

    cli.main(
        [
            "pipeline",
            "--bam",
            str(tmp_path / "inputs" / "x.bam"),
            "-o",
            str(tmp_path / "outputs"),
            "--output-name",
            "output",
        ]
    )

    assert called == [True]


def test_no_command_prints_help_and_exits_zero(capsys) -> None:
    """Unchanged behaviour, pinned here because dispatch now sits right after it."""
    with pytest.raises(SystemExit) as excinfo:
        cli.main([])

    assert excinfo.value.code == 0
    assert "usage:" in capsys.readouterr().out
