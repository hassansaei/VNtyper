"""End-to-end orchestration tests for the isolated mutation sweep."""

from __future__ import annotations

import contextlib
import hashlib
import os
import signal
import sys
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import mutation_test  # noqa: E402
import mutation_workspace_fs  # noqa: E402


@dataclass
class OrchestrationHarness:
    """Observable state owned by one fully isolated ``main`` invocation."""

    real_root: Path
    sweep_root: Path
    real_target: Path
    output: Path
    results_json: Path
    capability: mutation_workspace_fs.RootCapability
    events: list[str] = field(default_factory=list)
    execution_roots: list[object] = field(default_factory=list)
    cache_roots: list[Path] = field(default_factory=list)


def _install_orchestration_harness(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    failure: str | None = None,
    stub_signals: bool = True,
) -> OrchestrationHarness:
    """Install deterministic boundaries around one real ``main`` orchestration."""
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    real_target = real_root / "sample.py"
    sweep_target = sweep_root / "sample.py"
    real_target.parent.mkdir(parents=True)
    sweep_target.parent.mkdir(parents=True)
    real_target.write_text("VALUE = 1\n", encoding="utf-8")
    sweep_target.write_text("VALUE = 1\n", encoding="utf-8")
    capability = mutation_workspace_fs.open_root_capability(sweep_root)
    harness = OrchestrationHarness(
        real_root=real_root,
        sweep_root=sweep_root,
        real_target=real_target,
        output=real_root / "report.md",
        results_json=real_root / "results.json",
        capability=capability,
    )
    verify_labels = ["verify-after-canary", "verify-after-sweep"]

    @dataclass
    class FakeWorkspace:
        real_root: Path
        sweep_root: Path
        root_capability: mutation_workspace_fs.RootCapability
        head: str = "a" * 40
        overlay_changes: tuple[object, ...] = ()
        baseline_manifest: tuple[object, ...] = ()

        def verify_baseline(self) -> None:
            label = verify_labels.pop(0)
            harness.events.append(label)
            if failure == label:
                message = "canary restore" if label == "verify-after-canary" else "ordinary restore"
                raise RuntimeError(message)

    workspace = FakeWorkspace(real_root, sweep_root, capability)

    @contextlib.contextmanager
    def fake_detached(
        requested_root: Path,
        targets: tuple[str, ...],
        excluded_outputs: tuple[Path, ...],
    ) -> Iterator[FakeWorkspace]:
        assert requested_root == real_root
        assert targets == ("sample.py",)
        assert excluded_outputs == (harness.output, harness.results_json)
        if failure == "workspace":
            os.close(capability.descriptor)
            raise RuntimeError("workspace add failed")
        harness.events.append("workspace-enter")
        try:
            yield workspace
        finally:
            harness.events.append("workspace-exit")
            os.close(capability.descriptor)
            if failure == "cleanup":
                raise RuntimeError(f"orphaned worktree: {sweep_root}")

    def fake_dirty_guard(targets: object, outputs: object) -> str | None:
        assert tuple(targets) == ("sample.py",)
        assert tuple(outputs) == (harness.output, harness.results_json)
        harness.events.append("dirty-guard")
        if failure == "dirty":
            os.close(capability.descriptor)
        return "dirty selected target" if failure == "dirty" else None

    def fake_capture(root: Path, targets: object) -> dict[str, str]:
        assert root == real_root
        assert tuple(targets) == ("sample.py",)
        harness.events.append("real-digest-captured")
        if failure == "digest-capture":
            os.close(capability.descriptor)
            raise OSError("digest read failed")
        return {"sample.py": "startup-digest"}

    def fake_verify_real(root: Path, expected: object) -> None:
        assert root == real_root
        assert expected == {"sample.py": "startup-digest"}
        assert real_target.read_text(encoding="utf-8") == "VALUE = 1\n"
        harness.events.append("real-digest-verified")
        if failure in {"digest-verify", "sweep-and-digest-verify"}:
            raise RuntimeError("real target digest changed")

    def fake_provenance(actual_workspace: object, modules: tuple[str, ...]) -> tuple[Path, ...]:
        assert actual_workspace is workspace
        assert modules == ("sample.py",)
        harness.events.append("provenance")
        if failure == "provenance":
            raise RuntimeError("escaped import")
        return (sweep_root / "vntyper/__init__.py", sweep_target)

    def fake_cache_clear(root: Path) -> int:
        expected = Path("/proc/self/fd") / str(capability.descriptor) / "vntyper"
        assert root == expected
        harness.cache_roots.append(root)
        harness.events.append("cache-clear")
        if failure == "cache-clear":
            raise OSError("cache clear failed")
        return 0

    def fake_baseline(targets: object, timeout: int, *, repo_root: object) -> str | None:
        assert tuple(targets) == ("sample.py",)
        assert timeout == 600
        assert repo_root is capability
        harness.execution_roots.append(repo_root)
        harness.events.append("baseline")
        return "baseline red" if failure == "baseline" else None

    def fake_canary(timeout: int, *, repo_root: object) -> str | None:
        assert timeout == 600
        assert repo_root is capability
        harness.execution_roots.append(repo_root)
        harness.events.append("canary")
        return "canary survived" if failure == "canary" else None

    def fake_sweep(
        path: str,
        tests: tuple[str, ...],
        timeout: int,
        verbose: bool,
        *,
        repo_root: object,
    ) -> mutation_test.ModuleResult:
        assert (path, tests, timeout, verbose) == ("sample.py", ("tests/unit/x.py",), 600, False)
        assert repo_root is capability
        harness.execution_roots.append(repo_root)
        harness.events.append("sweep")
        if failure == "interrupt":
            raise KeyboardInterrupt("operator interrupt")
        if failure in {"sweep", "sweep-and-digest-verify"}:
            raise RuntimeError("sweep failed")
        return mutation_test.ModuleResult(path=path, killed=1)

    def fake_outputs(
        results: list[mutation_test.ModuleResult],
        elapsed: float,
        output: Path | None,
        results_json: Path | None,
    ) -> None:
        assert len(results) == 1 and results[0].path == "sample.py"
        assert elapsed == 2.0
        assert output == harness.output
        assert results_json == harness.results_json
        harness.events.append("outputs")
        if failure == "outputs":
            raise OSError("json replace")
        harness.results_json.write_bytes(b"new-json")
        harness.output.write_bytes(b"new-markdown")

    monotonic = iter((10.0, 12.0))
    monkeypatch.setattr(mutation_test, "REAL_REPO_ROOT", real_root)
    monkeypatch.setattr(
        mutation_test,
        "TARGETS",
        {"sample.py": ("tests/unit/x.py",), "other.py": ("tests/unit/y.py",)},
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "mutation_test.py",
            "--module",
            "sample.py",
            "--output",
            str(harness.output),
            "--results-json",
            str(harness.results_json),
        ],
    )
    if stub_signals:
        monkeypatch.setattr(
            mutation_test,
            "_install_signal_handlers",
            lambda: harness.events.append("signal-install") or (),
            raising=False,
        )
    monkeypatch.setattr(mutation_test, "_refuse_if_dirty", fake_dirty_guard)
    monkeypatch.setattr(mutation_test, "_real_target_digests", fake_capture, raising=False)
    monkeypatch.setattr(mutation_test, "_verify_real_target_digests", fake_verify_real, raising=False)
    monkeypatch.setattr(mutation_test, "detached_head_workspace", fake_detached, raising=False)
    monkeypatch.setattr(mutation_test, "verify_import_provenance", fake_provenance, raising=False)
    monkeypatch.setattr(mutation_test, "clear_bytecode_caches", fake_cache_clear)
    monkeypatch.setattr(mutation_test, "baseline_refusal", fake_baseline)
    monkeypatch.setattr(mutation_test, "canary_refusal", fake_canary)
    monkeypatch.setattr(mutation_test, "sweep_module", fake_sweep)
    monkeypatch.setattr(mutation_test, "write_outputs", fake_outputs)
    monkeypatch.setattr(mutation_test.time, "monotonic", lambda: next(monotonic))
    return harness


def test_main_orders_capabilities_measurement_cleanup_and_digest_verification(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Removing or reordering any safety stage must fail before publication."""
    harness = _install_orchestration_harness(tmp_path, monkeypatch)

    assert mutation_test.main() == 0

    assert harness.events == [
        "dirty-guard",
        "real-digest-captured",
        "signal-install",
        "workspace-enter",
        "provenance",
        "cache-clear",
        "baseline",
        "canary",
        "verify-after-canary",
        "sweep",
        "verify-after-sweep",
        "outputs",
        "workspace-exit",
        "real-digest-verified",
    ]
    assert harness.execution_roots == [harness.capability] * 3
    assert harness.cache_roots
    assert harness.output.read_bytes() == b"new-markdown"
    assert harness.results_json.read_bytes() == b"new-json"
    assert harness.real_target.read_bytes() == b"VALUE = 1\n"


def test_real_target_digests_detect_a_changed_selected_source(tmp_path: Path) -> None:
    """A real-source byte change across the sweep must force a non-zero outcome."""
    first = tmp_path / "first.py"
    second = tmp_path / "second.py"
    first.write_bytes(b"FIRST\n")
    second.write_bytes(b"SECOND\n")
    expected = mutation_test._real_target_digests(tmp_path, {"first.py": object(), "second.py": object()})

    assert expected == {
        "first.py": hashlib.sha256(b"FIRST\n").hexdigest(),
        "second.py": hashlib.sha256(b"SECOND\n").hexdigest(),
    }

    second.write_bytes(b"CHANGED\n")
    with pytest.raises(RuntimeError, match="second.py"):
        mutation_test._verify_real_target_digests(tmp_path, expected)


@pytest.mark.parametrize(
    ("failure", "required_text", "forbidden_event"),
    [
        ("workspace", "workspace add failed", "provenance"),
        ("provenance", "escaped import", "baseline"),
        ("cache-clear", "cache clear failed", "baseline"),
        ("baseline", "baseline red", "canary"),
        ("canary", "canary survived", "sweep"),
        ("verify-after-canary", "canary restore", "sweep"),
        ("sweep", "sweep failed", "outputs"),
        ("interrupt", "operator interrupt", "outputs"),
        ("verify-after-sweep", "ordinary restore", "outputs"),
    ],
)
def test_orchestration_unwinds_without_publishing_on_stage_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    failure: str,
    required_text: str,
    forbidden_event: str,
) -> None:
    """Every pre-publication failure must clean up and preserve prior artifacts."""
    harness = _install_orchestration_harness(tmp_path, monkeypatch, failure=failure)
    harness.output.write_bytes(b"old-markdown")
    harness.results_json.write_bytes(b"old-json")

    assert mutation_test.main() == 1

    assert required_text in capsys.readouterr().err
    if failure != "workspace":
        assert "workspace-exit" in harness.events
    assert "real-digest-verified" in harness.events
    assert forbidden_event not in harness.events
    assert harness.output.read_bytes() == b"old-markdown"
    assert harness.results_json.read_bytes() == b"old-json"
    assert harness.real_target.read_bytes() == b"VALUE = 1\n"


def test_main_reports_output_failure_and_still_unwinds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """An output exception is non-zero after context and digest cleanup."""
    harness = _install_orchestration_harness(tmp_path, monkeypatch, failure="outputs")
    harness.output.write_bytes(b"old-markdown")
    harness.results_json.write_bytes(b"old-json")

    assert mutation_test.main() == 1

    assert "json replace" in capsys.readouterr().err
    assert harness.events[-2:] == ["workspace-exit", "real-digest-verified"]
    assert harness.output.read_bytes() == b"old-markdown"
    assert harness.results_json.read_bytes() == b"old-json"
    assert harness.real_target.read_bytes() == b"VALUE = 1\n"


def test_cleanup_failure_is_nonzero_and_a_fresh_retry_succeeds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A retained orphan never becomes the workspace for the next invocation."""
    failed = _install_orchestration_harness(tmp_path / "first", monkeypatch, failure="cleanup")

    assert mutation_test.main() == 1

    assert f"orphaned worktree: {failed.sweep_root}" in capsys.readouterr().err
    assert failed.output.read_bytes() == b"new-markdown"
    assert failed.results_json.read_bytes() == b"new-json"
    assert failed.real_target.read_bytes() == b"VALUE = 1\n"
    assert "real-digest-verified" in failed.events

    retry = _install_orchestration_harness(tmp_path / "retry", monkeypatch)
    assert retry.sweep_root != failed.sweep_root
    assert mutation_test.main() == 0
    assert "outputs" in retry.events
    assert retry.real_target.read_bytes() == b"VALUE = 1\n"


def test_digest_capture_failure_is_nonzero_without_a_false_verification(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Verification runs only after a complete startup digest set exists."""
    harness = _install_orchestration_harness(tmp_path, monkeypatch, failure="digest-capture")

    assert mutation_test.main() == 1

    assert "digest read failed" in capsys.readouterr().err
    assert harness.events == ["dirty-guard", "real-digest-captured"]


def test_digest_verification_failure_forces_nonzero_after_complete_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A real-target mismatch invalidates an otherwise completed measurement."""
    harness = _install_orchestration_harness(tmp_path, monkeypatch, failure="digest-verify")

    assert mutation_test.main() == 1

    assert "real target digest changed" in capsys.readouterr().err
    assert harness.output.read_bytes() == b"new-markdown"
    assert harness.results_json.read_bytes() == b"new-json"


def test_digest_verification_does_not_hide_the_original_sweep_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Both failures remain diagnostic when the sweep and final proof fail."""
    _install_orchestration_harness(tmp_path, monkeypatch, failure="sweep-and-digest-verify")

    assert mutation_test.main() == 1

    error = capsys.readouterr().err
    assert "sweep failed" in error
    assert "real target digest changed" in error


@pytest.mark.parametrize("signal_name", ["SIGTERM", "SIGHUP", "SIGQUIT"])
def test_handled_signal_unwinds_workspace(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    signal_name: str,
) -> None:
    """Every catchable termination signal must enter the common cleanup path."""
    signum = getattr(signal, signal_name, None)
    if signum is None:
        pytest.skip(f"{signal_name} is unavailable on this platform")
    harness = _install_orchestration_harness(tmp_path, monkeypatch, stub_signals=False)
    registered: dict[int, Callable[[int, object], None]] = {}
    monkeypatch.setattr(mutation_test.signal, "signal", lambda number, handler: registered.__setitem__(number, handler))

    def interrupting_sweep(*_args: object, **_kwargs: object) -> mutation_test.ModuleResult:
        registered[signum](signum, None)
        raise AssertionError("signal handler must raise")

    monkeypatch.setattr(mutation_test, "sweep_module", interrupting_sweep)

    assert mutation_test.main() == 1
    assert "workspace-exit" in harness.events
    assert "real-digest-verified" in harness.events
    assert "terminated by signal" in capsys.readouterr().err


def test_sigint_is_not_replaced(monkeypatch: pytest.MonkeyPatch) -> None:
    """SIGINT retains Python's native ``KeyboardInterrupt`` semantics."""
    registered: dict[int, object] = {}
    monkeypatch.setattr(mutation_test.signal, "signal", lambda number, handler: registered.__setitem__(number, handler))

    mutation_test._install_signal_handlers()

    assert signal.SIGINT not in registered
