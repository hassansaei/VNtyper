"""Adversarial source-write and pytest-capability tests for the mutation harness."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_source
import mutation_test

pytestmark = pytest.mark.unit

_SOURCE = "def f(value):\n    return value >= 1 and value\n"
_CANARY_SOURCE = "\n" * 74 + "VALUE = 6 / 3\n"
_CANARY_NODE = "tests/unit/test_scoring.py::test_split_depth_and_calculate_frame_score_no_frameshift"
_CANARY_ASSERTION = "E       AssertionError: Frame_Score should be 1.0"
_CANARY_FAILURE = f"FAILED {_CANARY_NODE}\n{_CANARY_ASSERTION}\n1 failed in 0.01s\n"


def _write_source(root: Path, relative: str = "sample.py", content: str = _SOURCE) -> Path:
    target = root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(content, encoding="utf-8")
    return target


def test_sweep_replaces_a_hardlink_entry_without_mutating_its_external_alias(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sweep_root = tmp_path / "sweep"
    sweep_root.mkdir()
    protected = tmp_path / "real-source.py"
    protected.write_text(_SOURCE, encoding="utf-8")
    target = sweep_root / "sample.py"
    target.hardlink_to(protected)
    observed_mutants: list[str] = []

    def kill_mutant(_paths: tuple[str, ...], _timeout: int, _parallel: bool = False, **_kwargs: object) -> bool:
        mutant = target.read_text(encoding="utf-8")
        assert mutant != _SOURCE, "pytest must observe the installed mutant"
        assert protected.read_text(encoding="utf-8") == _SOURCE, "a hardlink alias must retain the original inode bytes"
        observed_mutants.append(mutant)
        return False

    monkeypatch.setattr(mutation_test, "run_tests", kill_mutant)

    result = mutation_test.sweep_module("sample.py", ("tests",), 5, False, repo_root=sweep_root)

    assert result.killed == result.total > 0
    assert len(observed_mutants) == result.total
    assert target.read_text(encoding="utf-8") == _SOURCE
    assert protected.read_text(encoding="utf-8") == _SOURCE
    assert target.stat().st_ino != protected.stat().st_ino


def test_sweep_rejects_a_symlink_target_without_reading_or_mutating_its_external_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sweep_root = tmp_path / "sweep"
    sweep_root.mkdir()
    protected = tmp_path / "real-source.py"
    protected.write_text(_SOURCE, encoding="utf-8")
    target = sweep_root / "sample.py"
    target.symlink_to(protected)
    runs: list[object] = []

    def record_run(*_args: object, **_kwargs: object) -> bool:
        runs.append(object())
        return False

    monkeypatch.setattr(mutation_test, "run_tests", record_run)

    with pytest.raises(RuntimeError, match="mutation target is not a regular file"):
        mutation_test.sweep_module("sample.py", ("tests",), 5, False, repo_root=sweep_root)

    assert runs == []
    assert target.is_symlink()
    assert protected.read_text(encoding="utf-8") == _SOURCE


def test_partial_mutant_temp_write_preserves_the_installed_source_and_cleans_the_temp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = _write_source(tmp_path)
    real_write = os.write
    source_write_calls = 0

    def partial_then_fail(descriptor: int, content: bytes) -> int:
        nonlocal source_write_calls
        source_write_calls += 1
        if source_write_calls == 1:
            return real_write(descriptor, content[: max(1, len(content) // 2)])
        if source_write_calls == 2:
            raise OSError("partial source write")
        return real_write(descriptor, content)

    monkeypatch.setattr(os, "write", partial_then_fail)
    monkeypatch.setattr(mutation_test, "run_tests", lambda *_a, **_k: False)

    with pytest.raises(OSError, match="partial source write"):
        mutation_test.sweep_module("sample.py", ("tests",), 5, False, repo_root=tmp_path)

    assert source_write_calls >= 2, "source installation must use checked descriptor writes"
    assert target.read_text(encoding="utf-8") == _SOURCE
    assert list(tmp_path.iterdir()) == [target]


def test_source_replacement_retries_a_temp_collision_without_removing_the_existing_entry(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = _write_source(tmp_path)
    collision = tmp_path / ".sample.py.mutation-collision"
    collision.write_text("not ours\n", encoding="utf-8")
    identities = iter((SimpleNamespace(hex="collision"), SimpleNamespace(hex="owned")))
    monkeypatch.setattr(mutation_source.uuid, "uuid4", lambda: next(identities))

    mutation_source.atomic_replace_source(tmp_path, "sample.py", b"replacement\n")

    assert target.read_bytes() == b"replacement\n"
    assert collision.read_text(encoding="utf-8") == "not ours\n"


def test_source_replacement_rejects_non_regular_targets_and_zero_progress_writes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    directory = tmp_path / "directory.py"
    directory.mkdir()
    with pytest.raises(RuntimeError, match="not a regular file"):
        mutation_source.read_source_bytes(tmp_path, "directory.py")
    with pytest.raises(RuntimeError, match="not a regular file"):
        mutation_source.atomic_replace_source(tmp_path, "directory.py", b"replacement\n")

    target = _write_source(tmp_path)
    monkeypatch.setattr(mutation_source.os, "write", lambda *_a, **_k: 0)
    with pytest.raises(OSError, match="made no progress"):
        mutation_source.atomic_replace_source(tmp_path, "sample.py", b"replacement\n")

    assert target.read_text(encoding="utf-8") == _SOURCE
    assert set(tmp_path.iterdir()) == {directory, target}


@pytest.mark.parametrize("error", [RuntimeError("pytest failed"), KeyboardInterrupt("SIGINT")])
def test_sweep_exposes_the_mutant_during_pytest_and_restores_after_base_exceptions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    error: BaseException,
) -> None:
    target = _write_source(tmp_path)
    observed: list[str] = []

    def interrupt(_paths: tuple[str, ...], _timeout: int, _parallel: bool = False, **_kwargs: object) -> bool:
        observed.append(target.read_text(encoding="utf-8"))
        raise error

    monkeypatch.setattr(mutation_test, "run_tests", interrupt)

    with pytest.raises(type(error), match=str(error)):
        mutation_test.sweep_module("sample.py", ("tests",), 5, False, repo_root=tmp_path)

    assert observed and observed[0] != _SOURCE
    assert target.read_text(encoding="utf-8") == _SOURCE


def test_pytest_uses_the_pinned_root_descriptor_and_explicit_required_plugins(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    real_root.mkdir()
    sweep_root.mkdir()
    captured: dict[str, object] = {}

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        captured.update(command=command, **kwargs)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(mutation_test, "REAL_REPO_ROOT", real_root)
    monkeypatch.setattr(mutation_test.subprocess, "run", fake_run)
    monkeypatch.setenv("PYTEST_ADDOPTS", "-p ambient_plugin")
    monkeypatch.setenv("PYTEST_PLUGINS", "ambient_plugin")

    run = mutation_test.run_pytest(("tests/unit/test_scoring.py",), 7, parallel=True, repo_root=sweep_root)

    assert run.passed is True
    pass_fds = captured["pass_fds"]
    assert isinstance(pass_fds, tuple) and len(pass_fds) == 1
    assert captured["cwd"] == Path("/proc/self/fd") / str(pass_fds[0])
    environment = captured["env"]
    assert isinstance(environment, dict)
    assert environment["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] == "1"
    assert environment["PYTHONDONTWRITEBYTECODE"] == "1"
    assert "PYTEST_ADDOPTS" not in environment
    assert "PYTEST_PLUGINS" not in environment
    assert environment["PYTHONPATH"].split(os.pathsep)[0] == str(captured["cwd"])
    command = captured["command"]
    assert isinstance(command, list)
    assert command[command.index("pytest_timeout") - 1 : command.index("pytest_timeout") + 1] == [
        "-p",
        "pytest_timeout",
    ]
    assert command[command.index("xdist.plugin") - 1 : command.index("xdist.plugin") + 1] == ["-p", "xdist.plugin"]
    assert command[-1] == "tests/unit/test_scoring.py"


def test_pytest_executes_the_opened_root_when_the_public_directory_is_replaced(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sweep_root = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    sweep_root.mkdir()
    (sweep_root / "sentinel.txt").write_text("opened-root\n", encoding="utf-8")
    observed: list[str] = []

    def replace_public_root(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        cwd_argument = kwargs["cwd"]
        assert isinstance(cwd_argument, (str, os.PathLike))
        cwd = Path(cwd_argument)
        sweep_root.rename(displaced)
        sweep_root.mkdir()
        (sweep_root / "sentinel.txt").write_text("replacement-root\n", encoding="utf-8")
        try:
            observed.append((cwd / "sentinel.txt").read_text(encoding="utf-8"))
        finally:
            shutil.rmtree(sweep_root)
            displaced.rename(sweep_root)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(mutation_test.subprocess, "run", replace_public_root)

    run = mutation_test.run_pytest(("tests/unit",), 7, repo_root=sweep_root)

    assert run.passed is True
    assert observed == ["opened-root\n"]


def test_pytest_rejects_an_absolute_node_before_opening_a_child(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[object] = []

    def record_child(*_args: object, **_kwargs: object) -> None:
        calls.append(object())

    monkeypatch.setattr(mutation_test.subprocess, "run", record_child)

    with pytest.raises(ValueError, match="workspace-relative"):
        mutation_test.run_pytest((str(tmp_path / "test_replacement.py"),), 7, repo_root=tmp_path)

    assert calls == []


@pytest.mark.parametrize(
    ("returncode", "output", "expected"),
    [
        (1, _CANARY_FAILURE, None),
        (1, "ERROR: plugin failed during pytest configuration\n", "canary infrastructure failure"),
        (1, f"{_CANARY_FAILURE}INTERNALERROR> plugin teardown failed\n", "canary infrastructure failure"),
        (1, _CANARY_FAILURE.replace("1 failed", "11 failed"), "canary infrastructure failure"),
        (
            1,
            f"FAILED tests/unit/test_other.py::test_other - AssertionError\n{_CANARY_FAILURE}",
            "canary infrastructure failure",
        ),
        (1, f"FAILED {_CANARY_NODE}\n", "canary infrastructure failure"),
        (1, f"{_CANARY_ASSERTION}\n", "canary infrastructure failure"),
        (0, "1 passed\n", "canary survived"),
        (2, _CANARY_FAILURE, "canary infrastructure failure"),
    ],
)
def test_canary_accepts_only_the_expected_scoring_node_and_assertion(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    returncode: int,
    output: str,
    expected: str | None,
) -> None:
    target = _write_source(tmp_path, mutation_test.CANARY_KEY[0], _CANARY_SOURCE)
    original = target.read_bytes()
    observed: list[bytes] = []

    def fake_run_pytest(*_args: object, **_kwargs: object) -> mutation_test.TestRun:
        observed.append(target.read_bytes())
        return mutation_test.TestRun(returncode == 0, output, returncode, False)

    monkeypatch.setattr(mutation_test, "run_pytest", fake_run_pytest)

    refusal = mutation_test.canary_refusal(17, repo_root=tmp_path)

    assert observed == [original.replace(b" / ", b" * ")]
    if expected is None:
        assert refusal is None
    else:
        assert refusal is not None and expected in refusal
    assert target.read_bytes() == original


def test_canary_breaks_a_hardlink_before_installing_and_restores_its_own_name(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protected = tmp_path / "real-scoring.py"
    protected.write_text(_CANARY_SOURCE, encoding="utf-8")
    target = tmp_path / mutation_test.CANARY_KEY[0]
    target.parent.mkdir(parents=True)
    target.hardlink_to(protected)

    def fake_run_pytest(*_args: object, **_kwargs: object) -> mutation_test.TestRun:
        assert b"VALUE = 6 * 3" in target.read_bytes()
        assert protected.read_text(encoding="utf-8") == _CANARY_SOURCE
        return mutation_test.TestRun(False, _CANARY_FAILURE, 1, False)

    monkeypatch.setattr(mutation_test, "run_pytest", fake_run_pytest)

    assert mutation_test.canary_refusal(17, repo_root=tmp_path) is None
    assert target.read_text(encoding="utf-8") == _CANARY_SOURCE
    assert protected.read_text(encoding="utf-8") == _CANARY_SOURCE
    assert target.stat().st_ino != protected.stat().st_ino
