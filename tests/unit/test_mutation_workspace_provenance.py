from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_provenance
import mutation_workspace

pytestmark = pytest.mark.unit


def _workspace(real_root: Path, sweep_root: Path) -> mutation_workspace.MutationWorkspace:
    return mutation_workspace.MutationWorkspace(
        real_root=real_root,
        sweep_root=sweep_root,
        head="a" * 40,
        overlay_changes=(),
        baseline_manifest=(),
        baseline_status=b"",
        baseline_digests={},
    )


def test_provenance_probe_returns_paths_below_the_disposable_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    package = sweep_root / "vntyper/__init__.py"
    target = sweep_root / "vntyper/scripts/scoring.py"
    package.parent.mkdir(parents=True)
    target.parent.mkdir(parents=True)
    package.write_text("", encoding="utf-8")
    target.write_text("", encoding="utf-8")
    real_root.mkdir()
    seen: dict[str, object] = {}

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        seen.update(command=command, cwd=kwargs["cwd"], env=kwargs["env"])
        return subprocess.CompletedProcess(command, 0, f"{package}\n{target}\n", "")

    real_subdirectory = real_root / "scripts"
    monkeypatch.setenv(
        "PYTHONPATH",
        os.pathsep.join((str(real_root), str(real_subdirectory), str(tmp_path / "dependency"))),
    )
    monkeypatch.setattr(mutation_provenance.subprocess, "run", fake_run)

    paths = mutation_workspace.verify_import_provenance(
        _workspace(real_root, sweep_root),
        ("vntyper/scripts/scoring.py",),
    )

    assert paths == (package.resolve(), target.resolve())
    assert seen["cwd"] == sweep_root
    environment = seen["env"]
    assert isinstance(environment, dict)
    assert environment["PYTHONDONTWRITEBYTECODE"] == "1"
    assert environment["PYTHONPATH"].split(os.pathsep)[0] == str(sweep_root)
    assert str(real_root) not in environment["PYTHONPATH"].split(os.pathsep)
    assert str(real_subdirectory) not in environment["PYTHONPATH"].split(os.pathsep)
    command = seen["command"]
    assert isinstance(command, list)
    assert command[:2] == [sys.executable, "-B"]
    assert command[-1] == "vntyper.scripts.scoring"


def test_provenance_probe_rejects_an_editable_real_root_import(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    real_package = real_root / "vntyper/__init__.py"
    sweep_target = sweep_root / "vntyper/scripts/scoring.py"
    real_package.parent.mkdir(parents=True)
    sweep_target.parent.mkdir(parents=True)
    real_package.write_text("", encoding="utf-8")
    sweep_target.write_text("", encoding="utf-8")

    monkeypatch.setattr(
        mutation_provenance.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(
            command,
            0,
            f"{real_package}\n{sweep_target}\n",
            "",
        ),
    )

    with pytest.raises(RuntimeError) as exc_info:
        mutation_workspace.verify_import_provenance(
            _workspace(real_root, sweep_root),
            ("vntyper/scripts/scoring.py",),
        )

    message = str(exc_info.value)
    assert "import escaped disposable worktree" in message
    assert str(real_root) in message
    assert str(sweep_root) in message


def test_provenance_probe_rejects_a_non_python_target_before_spawning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    real_root.mkdir()
    sweep_root.mkdir()
    calls: list[list[str]] = []
    monkeypatch.setattr(
        mutation_provenance.subprocess,
        "run",
        lambda command, **_kwargs: calls.append(command),
    )

    with pytest.raises(RuntimeError, match="provenance target must be a Python module"):
        mutation_workspace.verify_import_provenance(_workspace(real_root, sweep_root), ("config.json",))

    assert calls == []


def test_provenance_probe_normalizes_spawn_and_child_failures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    real_root.mkdir()
    sweep_root.mkdir()
    workspace = _workspace(real_root, sweep_root)

    monkeypatch.setattr(
        mutation_provenance.subprocess,
        "run",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("spawn blocked")),
    )
    with pytest.raises(RuntimeError, match="import provenance probe failed: spawn blocked"):
        mutation_workspace.verify_import_provenance(workspace, ("vntyper/scripts/scoring.py",))

    monkeypatch.setattr(
        mutation_provenance.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(command, 2, "partial\n", "import failed"),
    )
    with pytest.raises(RuntimeError) as exc_info:
        mutation_workspace.verify_import_provenance(workspace, ("vntyper/scripts/scoring.py",))
    assert str(exc_info.value) == "import provenance probe failed: partial\nimport failed"


def test_provenance_probe_rejects_an_incomplete_path_list(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    package = sweep_root / "vntyper/__init__.py"
    package.parent.mkdir(parents=True)
    package.write_text("", encoding="utf-8")
    real_root.mkdir()
    monkeypatch.setattr(
        mutation_provenance.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(command, 0, f"{package}\n", ""),
    )

    with pytest.raises(RuntimeError, match="unexpected path count"):
        mutation_workspace.verify_import_provenance(
            _workspace(real_root, sweep_root),
            ("vntyper/scripts/scoring.py",),
        )


def test_provenance_probe_rejects_a_missing_reported_import(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    package = sweep_root / "vntyper/__init__.py"
    missing_target = sweep_root / "vntyper/scripts/missing.py"
    package.parent.mkdir(parents=True)
    package.write_text("", encoding="utf-8")
    real_root.mkdir()
    monkeypatch.setattr(
        mutation_provenance.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(
            command,
            0,
            f"{package}\n{missing_target}\n",
            "",
        ),
    )

    with pytest.raises(RuntimeError, match="import provenance probe returned an invalid path"):
        mutation_workspace.verify_import_provenance(
            _workspace(real_root, sweep_root),
            ("vntyper/scripts/scoring.py",),
        )
