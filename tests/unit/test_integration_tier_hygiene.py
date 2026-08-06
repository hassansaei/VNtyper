"""
tests/unit/test_integration_tier_hygiene.py

Guards the properties that made the integration tier untrustworthy as regression
evidence. Each assertion here corresponds to a defect that was actually present:

* ``tests/integration/test_pipeline_integration.py`` carried its own copies of
  ``test_config`` and ``ensure_test_data`` plus ``compute_md5``/``download_file``.
  The copies shadowed the canonical versions in ``tests/conftest.py`` and
  ``tests/support/data_utils.py``, and the copy of ``ensure_test_data`` ignored
  ``VNTYPER_TEST_DATA_SKIP_DOWNLOAD`` - so CI's fail-fast path was silently
  bypassed for the one tier that downloads 1.2 GB.
* It re-implemented the Kestrel assertions inline, and that inline copy skipped the
  ``Confidence`` check whenever a case expected ``Negative`` - exactly the class most
  likely to move under a filtering change.
* It called ``logging.basicConfig`` at import, which AGENTS.md forbids anywhere in
  this repository.
* No workflow ran the tier on any trigger, so none of the above was ever observed.

These are cheap source-level checks: no subprocess, no Docker, no test data.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

TESTS_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TESTS_ROOT.parent
INTEGRATION_MODULE = TESTS_ROOT / "integration" / "test_pipeline_integration.py"
ROOT_CONFTEST = TESTS_ROOT / "conftest.py"
DATA_UTILS = TESTS_ROOT / "support" / "data_utils.py"
WORKFLOWS = REPO_ROOT / ".github" / "workflows"


def _parse(path: Path) -> ast.Module:
    """Parse a Python source file into an AST.

    Args:
        path: File to parse.

    Returns:
        ast.Module: The parsed module.
    """
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _is_fixture_decorator(node: ast.expr) -> bool:
    """Is this decorator ``@pytest.fixture`` or ``@pytest.fixture(...)``?

    Args:
        node: A decorator expression.

    Returns:
        bool: True if it marks a pytest fixture.
    """
    if isinstance(node, ast.Call):
        node = node.func
    if isinstance(node, ast.Attribute):
        return node.attr == "fixture"
    return isinstance(node, ast.Name) and node.id == "fixture"


def _fixture_names(path: Path) -> set[str]:
    """Names of the pytest fixtures a module defines at top level.

    Args:
        path: Module to inspect.

    Returns:
        set[str]: Fixture names.
    """
    return {
        node.name
        for node in _parse(path).body
        if isinstance(node, ast.FunctionDef) and any(_is_fixture_decorator(dec) for dec in node.decorator_list)
    }


def _public_function_names(path: Path) -> set[str]:
    """Names of the top-level, non-underscore functions a module defines.

    Args:
        path: Module to inspect.

    Returns:
        set[str]: Function names.
    """
    return {
        node.name for node in _parse(path).body if isinstance(node, ast.FunctionDef) and not node.name.startswith("_")
    }


def test_the_integration_module_redefines_no_canonical_fixture() -> None:
    """The integration tier must use the fixtures in tests/conftest.py, not copies.

    A locally redefined fixture wins over the conftest one silently, so the two can
    drift indefinitely without any test failing.

    Raises:
        AssertionError: If a canonical fixture name is redefined.
    """
    duplicated = _fixture_names(INTEGRATION_MODULE) & _fixture_names(ROOT_CONFTEST)
    assert not duplicated, (
        f"{INTEGRATION_MODULE.name} redefines fixture(s) {sorted(duplicated)} that "
        f"already exist in {ROOT_CONFTEST.name}. The local copy shadows the canonical "
        "one, so the two tiers stop agreeing without anything going red. Delete the copy."
    )


def test_the_integration_module_redefines_no_data_utility() -> None:
    """Test-data fetching lives in tests/support/data_utils.py and nowhere else.

    Raises:
        AssertionError: If a data utility is re-implemented in the integration module.
    """
    duplicated = _public_function_names(INTEGRATION_MODULE) & _public_function_names(DATA_UTILS)
    assert not duplicated, (
        f"{INTEGRATION_MODULE.name} re-implements {sorted(duplicated)} from "
        f"tests/support/{DATA_UTILS.name}. Import them instead - the copies drifted on "
        "download timeouts, MD5 handling and archive prefix stripping."
    )


def test_the_integration_module_does_not_assert_kestrel_fields_itself() -> None:
    """Kestrel validation must be delegated to the shared orchestration.

    The inline copy this replaced skipped the ``Confidence`` assertion for every case
    expecting ``Negative``, so the two negative-control samples asserted nothing about
    the call itself. Delegating keeps the local and Docker tiers literally identical.

    Raises:
        AssertionError: If the module names Kestrel result columns directly, or does not
            call ``run_bam_test_case``.
    """
    source = INTEGRATION_MODULE.read_text(encoding="utf-8")

    assert "run_bam_test_case" in source, (
        f"{INTEGRATION_MODULE.name} must validate BAM cases through "
        "tests.support.orchestration.run_bam_test_case, which the Docker tier also uses."
    )

    inline = [column for column in ("Confidence", "Depth_Score", "Estimated_Depth_") if column in source]
    assert not inline, (
        f"{INTEGRATION_MODULE.name} names Kestrel result column(s) {inline} directly, "
        "which means it is asserting them inline again. Validation belongs in "
        "tests/helpers.py::validate_kestrel_output so both tiers share one implementation."
    )


def _calls_basic_config(path: Path) -> bool:
    """Does this module call ``logging.basicConfig(...)``?

    Matched on the AST, not on the text, so prose and assertions that merely name the
    function - including this module's own - do not count.

    Args:
        path: Module to inspect.

    Returns:
        bool: True if the module contains a call to ``basicConfig``.
    """
    return any(
        isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "basicConfig"
        for node in ast.walk(_parse(path))
    )


def test_no_collected_test_module_configures_the_root_logger() -> None:
    """`logging.basicConfig` is forbidden in every module pytest imports.

    It reconfigures the root logger for whatever process imports it - which, for a test
    module or a conftest, is every pytest run that collects it. The integration module
    called it at import, so it fired during collection of tiers that never run it.

    ``tests/benchmark/`` is excluded on purpose: those are standalone scripts with their
    own ``__main__``, never imported by pytest, and configuring logging is what a script
    entry point is supposed to do.

    Raises:
        AssertionError: If a collected module calls basicConfig.
    """
    candidates = [
        path for path in sorted(TESTS_ROOT.rglob("*.py")) if "benchmark" not in path.relative_to(TESTS_ROOT).parts
    ]
    offenders = [str(path.relative_to(REPO_ROOT)) for path in candidates if _calls_basic_config(path)]
    assert not offenders, (
        f"These pytest-imported modules call logging.basicConfig: {offenders}. AGENTS.md "
        "forbids it - declare `logger = logging.getLogger(__name__)` and let it propagate."
    )


def test_the_integration_tier_is_not_in_the_unit_tier() -> None:
    """Integration tests do real I/O and must never be selected by `-m unit`.

    Raises:
        AssertionError: If the integration module claims the unit marker.
    """
    source = INTEGRATION_MODULE.read_text(encoding="utf-8")
    assert "pytest.mark.unit" not in source, (
        f"{INTEGRATION_MODULE.name} carries the unit marker, but it runs the real "
        "pipeline against 1.2 GB of sample data. `make test-unit` must stay pure."
    )
    assert "pytest.mark.integration" in source, (
        f"{INTEGRATION_MODULE.name} defines no integration marker, so `-m integration` "
        "selects nothing and the tier silently does not exist."
    )


def test_a_scheduled_workflow_runs_the_integration_tier() -> None:
    """Some workflow must run the integration tier on a timer.

    The tier's dependencies - the Zenodo archive, the UCSC/NCBI/ENSEMBL reference URLs,
    the vendored Kestrel JAR - rot on their own schedule, not on the repository's. Until
    a cron looked for it, the first person to notice was whoever opened the next PR.

    Raises:
        AssertionError: If no scheduled workflow invokes the integration tier.
    """
    if not WORKFLOWS.is_dir():
        pytest.skip("no GitHub Actions workflows present in this tree")

    scheduled_with_integration = [
        path.name
        for path in sorted(WORKFLOWS.glob("*.yml"))
        for text in [path.read_text(encoding="utf-8")]
        if re.search(r"^\s*schedule:\s*$", text, re.MULTILINE)
        and re.search(r"^\s*-\s*cron:", text, re.MULTILINE)
        and "make test-integration" in text
    ]
    assert scheduled_with_integration, (
        "No workflow runs `make test-integration` on a `schedule:` trigger. "
        "ci-tests.yml gates on `-m unit` only, so without this the pipeline is never "
        "exercised end to end by CI at all."
    )
