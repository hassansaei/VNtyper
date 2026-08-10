"""Unit tests for atomic mutation-report installation."""

from __future__ import annotations

import json
import os
import stat
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest import mock

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_output
import mutation_test

pytestmark = pytest.mark.unit


class _FailingStream:
    """Proxy a real text stream while failing one requested durability phase."""

    def __init__(self, stream: Any, phase: str) -> None:
        self._stream = stream
        self._phase = phase

    def __enter__(self) -> _FailingStream:
        return self

    def __exit__(self, *_args: object) -> None:
        self._stream.close()

    def write(self, content: str) -> int:
        if self._phase == "write":
            raise OSError("write blocked")
        return self._stream.write(content)

    def flush(self) -> None:
        if self._phase == "flush":
            raise OSError("flush blocked")
        self._stream.flush()

    def fileno(self) -> int:
        return self._stream.fileno()


def test_atomic_write_replaces_complete_text(tmp_path: Path) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")

    mutation_output.atomic_write_text(destination, "new\n")

    assert destination.read_text(encoding="utf-8") == "new\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_atomic_write_preserves_the_prior_destination_mode(tmp_path: Path) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    destination.chmod(0o640)

    mutation_output.atomic_write_text(destination, "new\n")

    assert stat.S_IMODE(destination.stat().st_mode) == 0o640


def test_atomic_write_fsyncs_the_complete_file_and_the_replaced_parent_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    real_fsync = os.fsync
    synced_kinds: list[str] = []

    def record_fsync(descriptor: int) -> None:
        mode = os.fstat(descriptor).st_mode
        synced_kinds.append("directory" if stat.S_ISDIR(mode) else "regular")
        real_fsync(descriptor)

    monkeypatch.setattr(mutation_output.os, "fsync", record_fsync)

    mutation_output.atomic_write_text(destination, "new\n")

    assert synced_kinds == ["regular", "directory"]


@pytest.mark.parametrize("phase", ["write", "flush"])
def test_stream_failure_preserves_the_old_report_and_removes_the_temp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    phase: str,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    real_fdopen = mutation_output.os.fdopen

    def failing_fdopen(descriptor: int, *args: Any, **kwargs: Any) -> _FailingStream:
        return _FailingStream(real_fdopen(descriptor, *args, **kwargs), phase)

    monkeypatch.setattr(mutation_output.os, "fdopen", failing_fdopen)

    with pytest.raises(OSError, match=f"{phase} blocked"):
        mutation_output.atomic_write_text(destination, "new\n")

    assert destination.read_text(encoding="utf-8") == "old\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_file_fsync_failure_preserves_the_old_report_and_removes_the_temp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    monkeypatch.setattr(mutation_output.os, "fsync", mock.Mock(side_effect=OSError("fsync blocked")))

    with pytest.raises(OSError, match="fsync blocked"):
        mutation_output.atomic_write_text(destination, "new\n")

    assert destination.read_text(encoding="utf-8") == "old\n"
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


def test_cleanup_failure_is_visible_and_uses_the_opened_parent_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    cleanup_dir_fds: list[int | None] = []

    def fail_cleanup(_path: object, *args: object, **kwargs: object) -> None:
        del args
        directory_fd = kwargs.get("dir_fd")
        assert directory_fd is None or isinstance(directory_fd, int)
        cleanup_dir_fds.append(directory_fd)
        raise OSError("cleanup blocked")

    monkeypatch.setattr(mutation_output.os, "replace", mock.Mock(side_effect=OSError("replace blocked")))
    monkeypatch.setattr(mutation_output.os, "unlink", fail_cleanup)

    with pytest.raises(OSError, match="cleanup blocked"):
        mutation_output.atomic_write_text(destination, "new\n")

    assert cleanup_dir_fds and all(isinstance(descriptor, int) for descriptor in cleanup_dir_fds)
    assert destination.read_text(encoding="utf-8") == "old\n"


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

    assert len(closed) == 3
    assert len(set(closed)) == 3
    assert destination.read_text(encoding="utf-8") == "old\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_temp_name_collision_is_retried_without_removing_the_existing_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    collision = tmp_path / ".report.md.mutation-collision"
    collision.write_text("not ours\n", encoding="utf-8")
    identities = iter((SimpleNamespace(hex="collision"), SimpleNamespace(hex="owned")))
    monkeypatch.setattr(mutation_output.uuid, "uuid4", lambda: next(identities))

    mutation_output.atomic_write_text(destination, "new\n")

    assert destination.read_text(encoding="utf-8") == "new\n"
    assert collision.read_text(encoding="utf-8") == "not ours\n"


@pytest.mark.parametrize("interrupt_transition", [2, 3])
def test_interrupt_during_stream_descriptor_handoff_closes_every_descriptor_and_removes_temp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    interrupt_transition: int,
) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    transitions = 0

    @contextmanager
    def interrupt_second_transition() -> Iterator[None]:
        nonlocal transitions
        transitions += 1
        current = transitions
        yield
        if current == interrupt_transition:
            raise KeyboardInterrupt("SIGINT during descriptor handoff")

    monkeypatch.setattr(mutation_output, "_defer_termination", interrupt_second_transition)

    with pytest.raises(KeyboardInterrupt, match="descriptor handoff"):
        mutation_output.atomic_write_text(destination, "new\n")

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


def test_write_outputs_refuses_one_path_for_both_artifacts_before_installing_anything(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "combined.txt"
    destination.write_text("old\n", encoding="utf-8")
    writes = mock.Mock()
    monkeypatch.setattr(mutation_test, "atomic_write_text", writes)

    with pytest.raises(ValueError, match="distinct"):
        mutation_test.write_outputs(
            [mutation_test.ModuleResult(path="sample.py", killed=1)], 2.0, destination, destination
        )

    writes.assert_not_called()
    assert destination.read_text(encoding="utf-8") == "old\n"


def _prepare_path_validation_main(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    argv: list[str],
) -> tuple[Path, list[str]]:
    real_root = tmp_path / "real"
    real_root.mkdir(exist_ok=True)
    target = real_root / "sample.py"
    target.write_text("VALUE = 1\n", encoding="utf-8")
    events: list[str] = []
    monkeypatch.setattr(mutation_test, "REAL_REPO_ROOT", real_root)
    monkeypatch.setattr(mutation_test, "TARGETS", {"sample.py": ("tests/unit/x.py",)})
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(mutation_test, "_refuse_if_dirty", lambda *_a, **_k: None)
    monkeypatch.setattr(mutation_test.signal, "signal", lambda *_a, **_k: None)
    monkeypatch.setattr(mutation_test, "clear_bytecode_caches", lambda *_a, **_k: 0)
    monkeypatch.setattr(mutation_test, "baseline_refusal", lambda *_a, **_k: None)

    def record_sweep(*_args: object, **_kwargs: object) -> mutation_test.ModuleResult:
        events.append("sweep")
        return mutation_test.ModuleResult(path="sample.py", killed=1)

    def record_write(*_args: object, **_kwargs: object) -> None:
        events.append("write")

    monkeypatch.setattr(mutation_test, "sweep_module", record_sweep)
    monkeypatch.setattr(mutation_test, "write_outputs", record_write)
    monkeypatch.setattr(mutation_test, "dirty_paths", lambda *_a, **_k: [])
    return target, events


def test_main_refuses_an_output_that_is_the_selected_source_before_sweeping(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    target, events = _prepare_path_validation_main(tmp_path, monkeypatch, ["mutation_test.py", "--output", "sample.py"])
    original = target.read_bytes()

    assert mutation_test.main() == 1

    assert events == []
    assert target.read_bytes() == original
    assert "overlap" in capsys.readouterr().err


def test_main_refuses_output_and_results_that_resolve_to_the_same_path_before_sweeping(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _target, events = _prepare_path_validation_main(
        tmp_path,
        monkeypatch,
        ["mutation_test.py", "--output", "result.txt", "--results-json", "result.txt"],
    )

    assert mutation_test.main() == 1

    assert events == []
    assert "distinct" in capsys.readouterr().err


def test_main_refuses_an_existing_hardlink_output_alias_of_the_selected_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    target, events = _prepare_path_validation_main(tmp_path, monkeypatch, ["mutation_test.py", "--output", "report.md"])
    report = target.parent / "report.md"
    report.hardlink_to(target)

    assert mutation_test.main() == 1

    assert events == []
    assert target.read_text(encoding="utf-8") == "VALUE = 1\n"
    assert report.read_text(encoding="utf-8") == "VALUE = 1\n"
    assert "overlap" in capsys.readouterr().err


def test_render_only_refuses_to_replace_its_own_input_before_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    real_root = tmp_path / "real"
    real_root.mkdir()
    saved = real_root / "results.json"
    payload = json.dumps({"elapsed": 1.0, "modules": []}) + "\n"
    saved.write_text(payload, encoding="utf-8")
    events: list[str] = []
    monkeypatch.chdir(real_root)
    monkeypatch.setattr(mutation_test, "REAL_REPO_ROOT", real_root)
    monkeypatch.setattr(mutation_test, "_refuse_if_dirty", lambda *_a, **_k: None)
    monkeypatch.setattr(mutation_test, "write_outputs", lambda *_a, **_k: events.append("write"))
    monkeypatch.setattr(
        sys,
        "argv",
        ["mutation_test.py", "--render-only", "results.json", "--output", "results.json"],
    )

    assert mutation_test.main() == 1

    assert events == []
    assert saved.read_text(encoding="utf-8") == payload
    assert "distinct" in capsys.readouterr().err
