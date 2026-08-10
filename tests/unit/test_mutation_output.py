"""Unit tests for atomic mutation-report installation."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest import mock

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_output
import mutation_test

pytestmark = pytest.mark.unit


def test_atomic_write_replaces_complete_text(tmp_path: Path) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")

    mutation_output.atomic_write_text(destination, "new\n")

    assert destination.read_text(encoding="utf-8") == "new\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_failed_replace_preserves_old_text_and_removes_temp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    monkeypatch.setattr(mutation_output.os, "replace", mock.Mock(side_effect=OSError("blocked")), raising=False)

    with pytest.raises(OSError, match="blocked"):
        mutation_output.atomic_write_text(destination, "new\n")

    assert destination.read_text(encoding="utf-8") == "old\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_failed_stream_open_closes_descriptor_and_removes_temp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    real_close = mutation_output.os.close
    closed: list[int] = []

    def record_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)

    monkeypatch.setattr(mutation_output.os, "close", record_close)
    monkeypatch.setattr(mutation_output.os, "fdopen", mock.Mock(side_effect=OSError("cannot open stream")))

    with pytest.raises(OSError, match="cannot open stream"):
        mutation_output.atomic_write_text(destination, "new\n")

    assert len(closed) == 1
    assert destination.read_text(encoding="utf-8") == "old\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_second_output_failure_leaves_complete_json_and_old_complete_markdown(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    report = tmp_path / "report.md"
    results_json = tmp_path / "results.json"
    report.write_text("old report\n", encoding="utf-8")
    results_json.write_text('{"old": true}\n', encoding="utf-8")
    real_atomic_write = mutation_output.atomic_write_text
    destinations: list[Path] = []

    def fail_second(destination: Path, content: str) -> None:
        destinations.append(destination)
        if len(destinations) == 2:
            raise OSError("second replacement blocked")
        real_atomic_write(destination, content)

    monkeypatch.setattr(mutation_test, "atomic_write_text", fail_second, raising=False)
    result = mutation_test.ModuleResult(path="sample.py", killed=1)

    with pytest.raises(OSError, match="second replacement blocked"):
        mutation_test.write_outputs([result], 2.0, report, results_json)

    assert destinations == [results_json, report]
    assert json.loads(results_json.read_text(encoding="utf-8")) == mutation_test.results_to_dict([result], 2.0)
    assert report.read_text(encoding="utf-8") == "old report\n"


def test_render_failure_happens_before_either_output_is_replaced(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    report = tmp_path / "report.md"
    results_json = tmp_path / "results.json"
    report.write_text("old report\n", encoding="utf-8")
    results_json.write_text('{"old": true}\n', encoding="utf-8")
    writes = mock.Mock()
    monkeypatch.setattr(mutation_test, "atomic_write_text", writes)
    monkeypatch.setattr(
        mutation_test,
        "format_markdown",
        mock.Mock(side_effect=ValueError("render failed")),
    )

    with pytest.raises(ValueError, match="render failed"):
        mutation_test.write_outputs(
            [mutation_test.ModuleResult(path="sample.py", killed=1)],
            2.0,
            report,
            results_json,
        )

    writes.assert_not_called()
    assert results_json.read_text(encoding="utf-8") == '{"old": true}\n'
    assert report.read_text(encoding="utf-8") == "old report\n"


def test_render_only_resolves_outputs_at_the_real_root_without_a_workspace(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_root = tmp_path / "real"
    outside = tmp_path / "outside"
    real_root.mkdir()
    outside.mkdir()
    saved = real_root / "results.json"
    saved.write_text(json.dumps({"elapsed": 1.0, "modules": []}) + "\n", encoding="utf-8")
    relative_report = Path("docs/report.md")
    monkeypatch.chdir(outside)
    monkeypatch.setattr(mutation_test, "REAL_REPO_ROOT", real_root)
    monkeypatch.setattr(mutation_test, "_refuse_if_dirty", lambda *_a, **_k: None)
    monkeypatch.setattr(
        mutation_test,
        "detached_head_workspace",
        lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("render-only requested a workspace")),
        raising=False,
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["mutation_test.py", "--render-only", str(saved), "--output", str(relative_report)],
    )

    assert mutation_test.main() == 0
    report = real_root / relative_report
    assert report.is_file()
    assert report.read_text(encoding="utf-8").startswith("# Mutation Testing")
    assert not (outside / relative_report).exists()
