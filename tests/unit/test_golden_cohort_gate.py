"""Unit tests for the golden-cohort command dispatcher."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))


def _load_gate() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "golden_cohort_gate_under_test", REPO_ROOT / "scripts" / "golden_cohort_gate.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = _load_gate()


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (["launch", "--", "--bam", "x"], (["launch"], ["--bam", "x"])),
        (["matrix", "--data-dir", "tests/data"], (["matrix", "--data-dir", "tests/data"], [])),
    ],
)
def test_split_launch_argv_preserves_each_side(argv: list[str], expected: tuple[list[str], list[str]]) -> None:
    assert gate._split_launch_argv(argv) == expected


def test_parser_requires_a_subcommand() -> None:
    with pytest.raises(SystemExit) as exc_info:
        gate.build_parser().parse_args([])
    assert exc_info.value.code == 2


def test_cmd_matrix_returns_one_when_matrix_construction_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_matrix(_args: object) -> dict[str, object]:
        raise ValueError("matrix drift")

    monkeypatch.setattr(gate, "_matrix_from_args", fail_matrix)
    assert gate.cmd_matrix(SimpleNamespace(out=None)) == 1


def test_cmd_matrix_rejects_an_empty_matrix(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(gate, "_matrix_from_args", lambda _args: {"cases": [], "probes": []})
    assert gate.cmd_matrix(SimpleNamespace(out=None)) == 1


def test_cmd_matrix_writes_the_exact_built_matrix(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    built = {"cases": [{"case_id": "case-a"}], "probes": [], "cohort_cases": [], "check": {"ok": True}}
    output = tmp_path / "nested" / "matrix.json"
    monkeypatch.setattr(gate, "_matrix_from_args", lambda _args: built)

    assert gate.cmd_matrix(SimpleNamespace(out=output)) == 0
    assert output.read_text(encoding="utf-8") == json.dumps(built, indent=2)


def test_cmd_probe_reports_a_pinned_success_and_unpinned_leak(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    tree = tmp_path / "tree"
    tree.mkdir()
    results = iter(
        [
            {"file": str(tree / "vntyper/__init__.py"), "sys_path_0": "", "marker": True},
            {"file": "/editable/vntyper/__init__.py", "sys_path_0": "/tmp", "marker": True},
            {"file": str(tree / "vntyper/__init__.py"), "sys_path_0": "/tmp", "marker": True},
        ]
    )

    def fake_run(*args: Any, **kwargs: Any) -> subprocess.CompletedProcess[str]:
        del args, kwargs
        return subprocess.CompletedProcess(["python"], 0, stdout=json.dumps(next(results)) + "\n", stderr="")

    monkeypatch.setattr(gate.subprocess, "run", fake_run)
    assert gate.cmd_probe(SimpleNamespace(tree=tree, marker="vntyper.scripts.marker")) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["tree"] == str(tree)
    assert payload["results"]["pythonpath_pinned"]["file"] == str(tree / "vntyper/__init__.py")
    assert payload["unpinned_script_leaked"] is True


def test_cmd_probe_returns_one_when_the_pinned_probe_has_no_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    tree = tmp_path / "tree"
    tree.mkdir()
    responses = iter(
        [
            subprocess.CompletedProcess(["python"], 0, stdout="{}\n", stderr=""),
            subprocess.CompletedProcess(["python"], 0, stdout="{}\n", stderr=""),
            subprocess.CompletedProcess(["python"], 1, stdout="", stderr="pinned import failed"),
        ]
    )
    monkeypatch.setattr(gate.subprocess, "run", lambda *args, **kwargs: next(responses))

    assert gate.cmd_probe(SimpleNamespace(tree=tree, marker="vntyper.scripts.marker")) == 1
    payload = json.loads(capsys.readouterr().out)
    assert payload["results"]["pythonpath_pinned"] == {"error": "pinned import failed"}


def test_cmd_launch_forwards_the_wrapper_contract(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    observed: dict[str, object] = {}

    def fake_launch(**kwargs: object) -> int:
        observed.update(kwargs)
        return 7

    monkeypatch.setattr(gate.launcher, "launch", fake_launch)
    args = SimpleNamespace(
        tree=tmp_path,
        side="after",
        marker="vntyper.scripts.marker",
        expect_marker="present",
        commands_log=tmp_path / "commands.jsonl",
    )

    assert gate.cmd_launch(args, ["pipeline", "--bam", "sample.bam"]) == 7
    assert observed == {
        "tree": tmp_path,
        "side": "after",
        "marker": "vntyper.scripts.marker",
        "expect_marker": True,
        "commands_log": tmp_path / "commands.jsonl",
        "argv": ["pipeline", "--bam", "sample.bam"],
    }


@pytest.mark.parametrize(
    ("record", "expected"),
    [
        ({"launch_verified": True, "expectations_met": True, "cohort_results": {}}, 0),
        ({"launch_verified": False, "expectations_met": True, "cohort_results": {}}, 1),
        ({"launch_verified": True, "expectations_met": False, "cohort_results": {}}, 1),
        (
            {
                "launch_verified": True,
                "expectations_met": True,
                "cohort_results": {"cohort": {"blocked": True}},
            },
            1,
        ),
    ],
)
def test_cmd_run_converts_the_side_record_to_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, record: dict[str, object], expected: int
) -> None:
    monkeypatch.setattr(gate, "_matrix_from_args", lambda _args: {"cases": [{"case_id": "a"}]})
    monkeypatch.setattr(gate.runner, "run_side", lambda **kwargs: record)
    args = SimpleNamespace(
        tree=tmp_path / "tree",
        run_root=tmp_path / "run",
        side="after",
        marker="vntyper.scripts.marker",
        expect_marker="present",
        threads=4,
        advntr_threads=8,
        jobs=2,
        timeout=60,
        no_cohort=False,
    )
    assert gate.cmd_run(args) == expected


def test_cmd_run_returns_one_when_matrix_construction_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_matrix(_args: object) -> dict[str, object]:
        raise ValueError("matrix drift")

    monkeypatch.setattr(gate, "_matrix_from_args", fail_matrix)
    assert gate.cmd_run(SimpleNamespace()) == 1


@pytest.mark.parametrize(
    ("argv", "selected", "forwarded"),
    [
        (["matrix", "--data-dir", "tests/data"], "matrix", None),
        (["probe", "--tree", "/tree", "--marker", "marker"], "probe", None),
        (
            [
                "launch",
                "--tree",
                "/tree",
                "--side",
                "after",
                "--marker",
                "marker",
                "--expect-marker",
                "present",
                "--",
                "--bam",
                "x",
            ],
            "launch",
            ["--bam", "x"],
        ),
        (
            [
                "run",
                "--side",
                "after",
                "--tree",
                "/tree",
                "--run-root",
                "/run",
                "--marker",
                "marker",
                "--expect-marker",
                "present",
            ],
            "run",
            None,
        ),
        (["compare", "--before-root", "/before", "--after-root", "/after"], "compare", None),
    ],
)
def test_main_dispatches_only_the_selected_handler(
    monkeypatch: pytest.MonkeyPatch, argv: list[str], selected: str, forwarded: list[str] | None
) -> None:
    calls: list[tuple[str, object]] = []

    def handler(name: str, status: int):
        def invoke(*args: object) -> int:
            calls.append((name, args[1] if len(args) == 2 else None))
            return status

        return invoke

    statuses = {"matrix": 3, "probe": 4, "launch": 5, "run": 6, "compare": 7}
    for name, status in statuses.items():
        monkeypatch.setattr(gate, f"cmd_{name}", handler(name, status))

    assert gate.main(argv) == statuses[selected]
    assert calls == [(selected, forwarded)]
