import builtins
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from golden_cohort import artifacts, launcher  # noqa: E402

pytestmark = pytest.mark.unit


def _resolution(tree: Path, marker_present: bool, *, in_tree: bool = True) -> dict[str, object]:
    return {
        "tree": str(tree),
        "vntyper_file": str(tree / "vntyper/__init__.py"),
        "in_tree": in_tree,
        "marker": "vntyper.scripts.marker",
        "marker_present": marker_present,
        "error": None,
    }


def test_modules_use_the_single_golden_cohort_import_identity() -> None:
    assert artifacts.__name__ == "golden_cohort.artifacts"
    assert launcher.__name__ == "golden_cohort.launcher"


def test_launch_line_reports_actual_and_expected_marker_states(tmp_path: Path) -> None:
    line = launcher._launch_line("after", _resolution(tmp_path, True), expect_marker=False)
    assert line.startswith("GATE-LAUNCH")
    assert "marker_state=present" in line
    assert "expected_marker=absent" in line


@pytest.mark.parametrize(
    ("info", "expect_marker"),
    [
        ({"in_tree": False, "marker_present": True}, True),
        ({"in_tree": True, "marker_present": False}, True),
    ],
)
def test_launch_returns_abort_when_tree_or_marker_expectation_disagrees(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    info: dict[str, bool],
    expect_marker: bool,
) -> None:
    resolved = _resolution(tmp_path, info["marker_present"], in_tree=info["in_tree"])
    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: resolved)

    assert (
        launcher.launch(
            tree=tmp_path,
            side="after",
            marker="vntyper.scripts.marker",
            expect_marker=expect_marker,
            commands_log=None,
            argv=["pipeline"],
        )
        == 97
    )


def test_launch_returns_unhandled_code_when_cli_raises(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def raise_unhandled() -> None:
        raise RuntimeError("boom")

    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: _resolution(tmp_path, True))
    monkeypatch.setitem(sys.modules, "vntyper.cli", SimpleNamespace(main=raise_unhandled))

    assert (
        launcher.launch(
            tree=tmp_path,
            side="after",
            marker="vntyper.scripts.marker",
            expect_marker=True,
            commands_log=None,
            argv=["pipeline"],
        )
        == 98
    )


@pytest.mark.parametrize("exit_code", [0, 1])
def test_launch_returns_cli_exit_codes_unchanged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, exit_code: int
) -> None:
    def exit_cli() -> None:
        raise SystemExit(exit_code)

    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: _resolution(tmp_path, True))
    monkeypatch.setitem(sys.modules, "vntyper.cli", SimpleNamespace(main=exit_cli))

    assert (
        launcher.launch(
            tree=tmp_path,
            side="after",
            marker="vntyper.scripts.marker",
            expect_marker=True,
            commands_log=None,
            argv=["pipeline"],
        )
        == exit_code
    )


def test_resolve_reports_tree_and_marker_state() -> None:
    info = launcher.resolve(REPO_ROOT, "vntyper.scripts.pipeline_guards")
    assert info["in_tree"] is True
    assert info["marker_present"] is True
    assert info["error"] is None
    assert Path(info["vntyper_file"]).is_relative_to(REPO_ROOT)


def test_resolve_reports_an_unimportable_marker(tmp_path: Path) -> None:
    info = launcher.resolve(tmp_path, "missing_golden_cohort_parent.marker")
    assert info["in_tree"] is False
    assert info["marker_present"] is False
    assert str(info["error"]).startswith("find_spec(missing_golden_cohort_parent.marker): ModuleNotFoundError:")


def test_resolve_reports_an_import_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    original_import = builtins.__import__

    def fail_vntyper_import(name: str, *args: Any, **kwargs: Any) -> Any:
        if name == "vntyper":
            raise RuntimeError("import failed")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_vntyper_import)
    info = launcher.resolve(tmp_path, "vntyper.scripts.marker")
    assert info == {
        "tree": str(tmp_path),
        "vntyper_file": None,
        "in_tree": False,
        "marker": "vntyper.scripts.marker",
        "marker_present": None,
        "error": "RuntimeError: import failed",
    }


def test_resolve_returns_structured_import_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An unavailable candidate is returned as an unstartable resolution."""
    original_import = builtins.__import__

    def fail_vntyper_import(name: str, *args: Any, **kwargs: Any) -> Any:
        if name == "vntyper":
            raise ImportError("candidate import failed")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_vntyper_import)

    assert launcher.resolve(tmp_path, "vntyper.scripts.marker") == {
        "tree": str(tmp_path),
        "vntyper_file": None,
        "in_tree": False,
        "marker": "vntyper.scripts.marker",
        "marker_present": None,
        "error": "ImportError: candidate import failed",
    }


def test_launch_returns_abort_when_import_failed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    resolved = _resolution(tmp_path, False)
    resolved.update({"vntyper_file": None, "in_tree": False, "marker_present": None, "error": "ImportError: missing"})
    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: resolved)

    assert (
        launcher.launch(
            tree=tmp_path,
            side="before",
            marker="vntyper.scripts.marker",
            expect_marker=False,
            commands_log=None,
            argv=[],
        )
        == 97
    )


@pytest.mark.parametrize(("system_exit_code", "expected"), [(None, 0), ("usage", 1)])
def test_launch_normalizes_non_integer_system_exit_codes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    system_exit_code: object,
    expected: int,
) -> None:
    def exit_cli() -> None:
        raise SystemExit(system_exit_code)

    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: _resolution(tmp_path, True))
    monkeypatch.setitem(sys.modules, "vntyper.cli", SimpleNamespace(main=exit_cli))

    assert (
        launcher.launch(
            tree=tmp_path,
            side="after",
            marker="vntyper.scripts.marker",
            expect_marker=True,
            commands_log=None,
            argv=[],
        )
        == expected
    )


def test_launch_records_commands_and_returns_zero_when_cli_returns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    commands_log = tmp_path / "commands.jsonl"
    recorded: list[Path] = []
    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: _resolution(tmp_path, True))
    monkeypatch.setattr(launcher, "_record_commands", recorded.append)
    monkeypatch.setitem(sys.modules, "vntyper.cli", SimpleNamespace(main=lambda: None))

    assert (
        launcher.launch(
            tree=tmp_path,
            side="after",
            marker="vntyper.scripts.marker",
            expect_marker=True,
            commands_log=commands_log,
            argv=["report"],
        )
        == 0
    )
    assert recorded == [commands_log]


@pytest.mark.parametrize(
    "argv",
    [
        ["pipeline", "--bam", "/inputs/sample.bam", "--output-dir", "/outputs/sample"],
        ["cohort", "--input-dir", "/outputs", "--output-dir", "/outputs/cohort"],
    ],
)
def test_launch_uses_resolved_tree_and_exact_cli_argv_then_restores_process_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    argv: list[str],
) -> None:
    repository = tmp_path / "repository"
    repository.mkdir()
    alias = tmp_path / "repository-link"
    alias.symlink_to(repository, target_is_directory=True)
    original_cwd = Path.cwd()
    original_path_object = sys.path
    original_path = sys.path.copy()
    original_argv_object = sys.argv
    original_argv = sys.argv.copy()
    observed: dict[str, object] = {}

    def observe_resolution(tree: Path, marker: str) -> dict[str, object]:
        observed["resolve_tree"] = tree
        observed["resolve_cwd"] = Path.cwd()
        observed["resolve_path_head"] = sys.path[0]
        return _resolution(tree, True)

    def observe_cli() -> None:
        observed["cli_cwd"] = Path.cwd()
        observed["cli_path_head"] = sys.path[0]
        observed["cli_argv"] = sys.argv.copy()

    monkeypatch.setattr(launcher, "resolve", observe_resolution)
    monkeypatch.setitem(sys.modules, "vntyper.cli", SimpleNamespace(main=observe_cli))

    try:
        assert (
            launcher.launch(
                tree=alias,
                side="after",
                marker="vntyper.scripts.marker",
                expect_marker=True,
                commands_log=None,
                argv=argv,
            )
            == 0
        )
        assert observed == {
            "resolve_tree": repository.resolve(),
            "resolve_cwd": repository.resolve(),
            "resolve_path_head": str(repository.resolve()),
            "cli_cwd": repository.resolve(),
            "cli_path_head": str(repository.resolve()),
            "cli_argv": ["vntyper", *argv],
        }
        assert Path.cwd() == original_cwd
        assert sys.path is original_path_object
        assert sys.path == original_path
        assert sys.argv is original_argv_object
        assert sys.argv == original_argv
    finally:
        os.chdir(original_cwd)
        sys.path[:] = original_path
        sys.argv[:] = original_argv


def test_launch_restores_command_recorder_process_state(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    original_init = subprocess.Popen.__init__
    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: _resolution(tmp_path, True))
    monkeypatch.setitem(sys.modules, "vntyper.cli", SimpleNamespace(main=lambda: None))

    try:
        assert (
            launcher.launch(
                tree=tmp_path,
                side="after",
                marker="vntyper.scripts.marker",
                expect_marker=True,
                commands_log=tmp_path / "commands.jsonl",
                argv=["report"],
            )
            == 0
        )
        assert subprocess.Popen.__init__ is original_init
        assert (tmp_path / "commands.jsonl").parent.is_dir()
    finally:
        subprocess.Popen.__init__ = original_init  # type: ignore[method-assign]


def test_command_recording_is_confined_to_a_child_process(tmp_path: Path) -> None:
    commands_log = tmp_path / "commands.jsonl"
    child = f"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, {str(REPO_ROOT / "scripts")!r})
from golden_cohort import launcher

def fake_init(self, args, *rest, **kwargs):
    self.args = args
    self.returncode = 0

subprocess.Popen.__init__ = fake_init
launcher._record_commands(Path({str(commands_log)!r}))
subprocess.Popen(["tool", "a b"], shell=True)

def fail_open(*args, **kwargs):
    raise OSError("log unavailable")

launcher.open = fail_open
second = subprocess.Popen(["second"])
records = [json.loads(line) for line in Path({str(commands_log)!r}).read_text().splitlines()]
print(json.dumps({{"records": records, "survived": second.returncode == 0}}))
"""
    completed = subprocess.run([sys.executable, "-c", child], check=True, capture_output=True, text=True)
    result = json.loads(completed.stdout)
    assert result == {
        "records": [{"command": "tool 'a b'", "shell": True}],
        "survived": True,
    }


# --- Attribute markers -------------------------------------------------------------------
#
# The marker has to be present on one side and absent on the other, and `admissibility`
# refuses two sides that expected the same state. A branch that only *modifies* files
# therefore has nothing to name: #259 adds no module at all -- even
# `tests/unit/test_advntr_command.py`, which it grows by 116 lines, already exists on main.
# Both candidate module markers resolve on both sides, so neither can distinguish them.
#
# `module:attribute` gives such a branch a marker without inventing a throwaway module for
# the PR to carry, and it is a *stronger* witness than a sibling module would be: it is
# reached through the same import machinery as the code under test.

#: Present on the candidate side of #259 and absent on `origin/main`.
ATTRIBUTE_MARKER = "vntyper.modules.advntr.advntr_genotyping:resolve_advntr_threads"


def test_parse_marker_splits_an_attribute_from_its_module() -> None:
    assert launcher.parse_marker(ATTRIBUTE_MARKER) == (
        "vntyper.modules.advntr.advntr_genotyping",
        "resolve_advntr_threads",
    )


def test_parse_marker_leaves_a_bare_module_alone() -> None:
    assert launcher.parse_marker("vntyper.scripts.pipeline_guards") == ("vntyper.scripts.pipeline_guards", None)


@pytest.mark.parametrize("marker", ["vntyper.scripts.pipeline_guards:", ":resolve_advntr_threads", ":"])
def test_parse_marker_refuses_a_half_written_marker(marker: str) -> None:
    """A marker with an empty half must not quietly degrade into a module check.

    That degradation is the dangerous one: `--marker mod:` would check only `mod`, which
    exists on both sides, and the gate would then compare two sides it had "verified".
    """
    with pytest.raises(ValueError, match="marker"):
        launcher.parse_marker(marker)


def test_resolve_reports_a_present_attribute_marker() -> None:
    info = launcher.resolve(REPO_ROOT, ATTRIBUTE_MARKER)

    assert info["marker_present"] is True
    assert info["error"] is None


def test_resolve_reports_an_absent_attribute_as_absent_and_not_as_an_error() -> None:
    """This is the baseline side's answer, and it must be a clean False.

    An absent attribute reported as an error would make every `before` run abort with
    EXIT_ABORT rather than run, which reads as a broken harness rather than a working one.
    """
    module, _ = launcher.parse_marker(ATTRIBUTE_MARKER)

    info = launcher.resolve(REPO_ROOT, f"{module}:no_such_function")

    assert info["marker_present"] is False
    assert info["error"] is None
    # The attribute is what made it absent, not a mistyped module: the same module answers
    # present on its own, and would answer present on both sides of the comparison.
    assert launcher.resolve(REPO_ROOT, module)["marker_present"] is True


def test_resolve_reports_an_attribute_whose_module_cannot_be_imported() -> None:
    info = launcher.resolve(REPO_ROOT, "missing_golden_cohort_parent.marker:thing")

    assert info["marker_present"] is False
    assert str(info["error"]).startswith("import_module(missing_golden_cohort_parent.marker): ModuleNotFoundError:")


def test_resolve_reports_a_malformed_marker_rather_than_treating_it_as_a_module() -> None:
    info = launcher.resolve(REPO_ROOT, "vntyper.scripts.pipeline_guards:")

    assert info["marker_present"] is False
    assert "ValueError" in str(info["error"])


def test_a_bare_module_marker_is_still_answered_without_importing_the_module() -> None:
    """`find_spec` executes nothing, and bare markers must keep that property.

    Runs 1-3 used `vntyper.scripts.pipeline_guards`; importing a marker module would change
    what those runs measured, because the import happens before `vntyper.cli.main`
    configures logging.
    """
    module = "vntyper.scripts.pipeline_guards"
    sys.modules.pop(module, None)

    assert launcher.marker_is_present(module) is True
    assert module not in sys.modules, "a bare module marker must not be imported"


def test_launch_aborts_when_the_marker_itself_cannot_be_resolved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """`resolve` recorded marker-resolution failures and nothing consumed the record.

    A marker whose module raises on import reported `marker_state=absent`, which is
    indistinguishable from a genuinely absent attribute. On the side that expects the
    marker absent -- the baseline -- the run then proceeded with its witness void while the
    launch line asserted it held. The attribute form makes this newly reachable: a bare
    module's `find_spec` cannot fail this way for a module that exists.
    """
    resolution = _resolution(tmp_path, False)
    resolution["error"] = "import_module(boom): RuntimeError: deps missing"
    monkeypatch.setattr(launcher, "resolve", lambda *_args, **_kwargs: resolution)

    code = launcher.launch(
        tree=tmp_path, side="before", marker="boom:anything", expect_marker=False, commands_log=None, argv=[]
    )

    assert code == launcher.EXIT_ABORT
    assert "marker-unresolvable" in capsys.readouterr().out


def test_launch_still_runs_when_the_marker_resolved_cleanly(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The abort above must key on the error, not on the marker being absent."""
    monkeypatch.setattr(launcher, "resolve", lambda *_args, **_kwargs: _resolution(tmp_path, False))
    monkeypatch.setattr(launcher, "_run_cli", lambda _argv: 0, raising=False)

    code = launcher.launch(
        tree=tmp_path, side="before", marker="mod:attr", expect_marker=False, commands_log=None, argv=["--help"]
    )

    assert code != launcher.EXIT_ABORT


def test_parse_marker_refuses_more_than_one_colon() -> None:
    """`partition` would give ('a.b', 'c:d'), and `hasattr(mod, 'c:d')` is always False --
    so a mistyped marker reads as absent rather than as the error it is."""
    with pytest.raises(ValueError, match="marker"):
        launcher.parse_marker("a.b:c:d")
