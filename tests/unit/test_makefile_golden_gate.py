"""The optional golden gate runs real selection only when both roots are supplied."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_gate(
    tmp_path: Path, target: str, roots: dict[str, str], status: int = 0
) -> tuple[subprocess.CompletedProcess[str], list[str] | None]:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    capture = tmp_path / "argv.json"
    python = fake_bin / "python"
    python.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "Path(os.environ['GOLDEN_CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n"
        "sys.exit(int(os.environ['GOLDEN_STATUS']))\n",
        encoding="utf-8",
    )
    python.chmod(0o755)
    env = os.environ.copy()
    env.pop("VNTYPER_SIM_ROOT", None)
    env.pop("VNTYPER_ADVNTR_ROOT", None)
    env.update(roots)
    env.update(PATH=f"{fake_bin}{os.pathsep}{env['PATH']}", GOLDEN_CAPTURE=str(capture), GOLDEN_STATUS=str(status))
    result = subprocess.run(
        ["make", "--no-print-directory", target], cwd=REPO_ROOT, env=env, capture_output=True, text=True, check=False
    )
    return result, json.loads(capture.read_text(encoding="utf-8")) if capture.exists() else None


@pytest.mark.parametrize(
    "roots",
    [
        {},
        {"VNTYPER_SIM_ROOT": "/simulation"},
        {"VNTYPER_ADVNTR_ROOT": "/advntr"},
        {"VNTYPER_SIM_ROOT": "", "VNTYPER_ADVNTR_ROOT": "/advntr"},
    ],
)
def test_optional_gate_reports_absent_roots_without_running_pytest(tmp_path: Path, roots: dict[str, str]) -> None:
    result, argv = _run_gate(tmp_path, "test-golden-if-configured", roots)

    assert result.returncode == 0
    assert "Golden tier not run" in result.stdout
    assert argv is None


@pytest.mark.parametrize("status", [0, 17])
def test_optional_gate_selects_full_golden_and_propagates_failure(tmp_path: Path, status: int) -> None:
    # Supplied but nonexistent paths must reach pytest, whose root validation fails.
    result, argv = _run_gate(
        tmp_path,
        "test-golden-if-configured",
        {"VNTYPER_SIM_ROOT": str(tmp_path / "missing sim"), "VNTYPER_ADVNTR_ROOT": str(tmp_path / "missing advntr")},
        status,
    )

    assert (result.returncode == 0) == (status == 0)
    assert argv == ["-m", "pytest", "-m", "golden", "tests/golden", "-q", "-rs"]
    assert "Golden tier not run" not in result.stdout


@pytest.mark.parametrize("roots", [{}, {"VNTYPER_SIM_ROOT": "/simulation"}, {"VNTYPER_ADVNTR_ROOT": "/advntr"}])
def test_explicit_gate_requires_both_roots(tmp_path: Path, roots: dict[str, str]) -> None:
    result, argv = _run_gate(tmp_path, "test-golden", roots)

    assert result.returncode != 0
    assert "set" in result.stderr
    assert argv is None


def test_check_all_reaches_optional_golden_gate() -> None:
    declaration = next(
        line
        for line in (REPO_ROOT / "Makefile").read_text(encoding="utf-8").splitlines()
        if line.startswith("check-all:")
    )

    assert declaration.partition(":")[2].split().count("test-golden-if-configured") == 1


def test_agent_instructions_describe_conditional_golden_gate() -> None:
    instructions = (REPO_ROOT / "AGENTS.md").read_text(encoding="utf-8")

    assert "also runs golden when both corpus roots are set" in instructions
    assert "`make check-all` calls `test-golden-if-configured`" in instructions
    assert "The tier deliberately remains outside\n  `make check-all`" not in instructions
