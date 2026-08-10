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

    def replacement_init(self: Any, args: Any, *rest: Any, **kwargs: Any) -> None:
        del self, args, rest, kwargs

    def install_recorder(log_path: Path) -> Any:
        assert log_path == tmp_path / "commands.jsonl"
        subprocess.Popen.__init__ = replacement_init  # type: ignore[method-assign]
        return original_init

    monkeypatch.setattr(sys, "path", sys.path.copy())
    monkeypatch.setattr(launcher.os, "chdir", lambda tree: None)
    monkeypatch.setattr(launcher, "resolve", lambda tree, marker: _resolution(tmp_path, True))
    monkeypatch.setattr(launcher, "_record_commands", install_recorder)
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
