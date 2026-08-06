#!/usr/bin/env python3
"""
tests/unit/test_marker_hygiene.py

Guards the test suite against silently losing coverage in CI.

CI runs `pytest -m unit`. A test file under tests/unit/ that forgets
`pytestmark = pytest.mark.unit` is collected, deselected, and never run - with no
warning and a green build. That is not hypothetical: test_haplo_count_and_selection.py
sat unmarked with 30 passing tests covering the Issue #136/#145 genotype tie-breaking
logic, invisible to CI.

`--strict-markers` (pytest.ini) catches a *misspelled* marker but cannot catch a
*missing* one. This test closes that gap.

It works off the live session's collected items, so it costs no subprocess and no
measurable time: any file that contributes zero items to a `-m unit` run is either
unmarked or has no tests, and we distinguish those by source inspection.
"""

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

UNIT_DIR = Path(__file__).resolve().parent

# A module defines tests if it has a top-level `def test_*` or a `class Test*`.
_DEFINES_TESTS = re.compile(r"^\s*(?:async\s+def\s+test_|def\s+test_|class\s+Test)", re.MULTILINE)


def test_every_unit_file_is_selected_by_the_unit_marker(request: pytest.FixtureRequest) -> None:
    """Every test file under tests/unit/ must be selected by `pytest -m unit`.

    Args:
        request: Pytest request, used to reach the current session's collected items.

    Raises:
        AssertionError: If a file defining tests contributed no items to this run.
    """
    selected_files = {
        Path(str(item.fspath)).resolve() for item in request.session.items if str(item.fspath).endswith(".py")
    }

    unmarked = []
    for path in sorted(UNIT_DIR.rglob("test_*.py")):
        if path.resolve() in selected_files:
            continue
        if _DEFINES_TESTS.search(path.read_text(encoding="utf-8")):
            unmarked.append(str(path.relative_to(UNIT_DIR)))

    assert not unmarked, (
        "These files under tests/unit/ define tests but contributed nothing to a "
        f"`pytest -m unit` run, so CI never executes them: {unmarked}. "
        "Add `pytestmark = pytest.mark.unit` immediately after the imports."
    )


def test_no_test_modules_outside_the_known_tiers() -> None:
    """Test modules must live under tests/unit, tests/integration or tests/docker.

    A ``test_*.py`` anywhere else is either a mis-placed test that no tier runs,
    or a helper module whose name makes a real collection failure invisible.

    Raises:
        AssertionError: If a stray test module is found.
    """
    tests_root = UNIT_DIR.parent
    allowed = {"unit", "integration", "docker"}
    stray = [
        str(path.relative_to(tests_root))
        for path in tests_root.rglob("test_*.py")
        if path.relative_to(tests_root).parts[0] not in allowed
    ]
    assert not stray, (
        f"These test modules live outside tests/{{unit,integration,docker}}: {stray}. "
        "Move real tests into a tier; rename helpers so they do not start with 'test_'."
    )
