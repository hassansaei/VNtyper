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
import shlex
from pathlib import Path, PurePosixPath

import pytest
import yaml

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


def _call_nodes(path: Path, name: str) -> list[ast.Call]:
    """Every invocation of ``name`` in a module, as AST call nodes.

    Matched on ``ast.Call``, so a bare name reference (``_ = run_bam_test_case``), an
    import, or a mention in a docstring does not count - only a real invocation. Both
    ``f(...)`` and ``mod.f(...)`` forms match, keyed on the final attribute name.

    Args:
        path: Module to inspect.
        name: Function name to look for.

    Returns:
        list[ast.Call]: The matching call nodes, in traversal order.
    """
    return [
        node
        for node in ast.walk(_parse(path))
        if isinstance(node, ast.Call)
        and (
            (isinstance(node.func, ast.Name) and node.func.id == name)
            or (isinstance(node.func, ast.Attribute) and node.func.attr == name)
        )
    ]


def test_the_integration_module_does_not_assert_kestrel_fields_itself() -> None:
    """Kestrel validation must be delegated to the shared orchestration.

    The inline copy this replaced skipped the ``Confidence`` assertion for every case
    expecting ``Negative``, so the two negative-control samples asserted nothing about
    the call itself. Delegating keeps the local and Docker tiers literally identical.

    The delegation is checked on the AST rather than with a substring search: an import,
    a docstring mention or a bare name reference (``_ = run_bam_test_case``) all contain
    the string while validating nothing, and a substring guard stays green through
    exactly that edit.

    Raises:
        AssertionError: If the module names Kestrel result columns directly, or does not
            actually call ``run_bam_test_case`` with arguments.
    """
    source = INTEGRATION_MODULE.read_text(encoding="utf-8")

    calls = _call_nodes(INTEGRATION_MODULE, "run_bam_test_case")
    assert calls, (
        f"{INTEGRATION_MODULE.name} must validate BAM cases through "
        "tests.support.orchestration.run_bam_test_case, which the Docker tier also uses. "
        "No call to it was found in the AST - importing the name, or referencing it "
        "without calling it, validates nothing."
    )
    assert any(call.args or call.keywords for call in calls), (
        f"{INTEGRATION_MODULE.name} calls run_bam_test_case() without passing the case, "
        "the runner or the output directory, so it cannot be validating anything."
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


# Shell tokens that may precede a command without changing which command runs.
_TRANSPARENT_PREFIXES = frozenset({"sudo", "time", "env", "nice", "exec", "command", "xvfb-run"})


def _is_cron_scheduled(workflow: dict) -> bool:
    """Does this parsed workflow declare a ``schedule:`` trigger with a cron entry?

    PyYAML resolves the bare key ``on`` to the boolean ``True`` under YAML 1.1, so the
    trigger block has to be looked up under both spellings.

    Args:
        workflow: A parsed GitHub Actions workflow document.

    Returns:
        bool: True if at least one cron schedule is declared.
    """
    triggers = workflow.get("on", workflow.get(True))
    if not isinstance(triggers, dict):
        return False
    schedule = triggers.get("schedule")
    if not isinstance(schedule, list):
        return False
    return any(isinstance(entry, dict) and entry.get("cron") for entry in schedule)


def _run_scripts(workflow: dict) -> list[str]:
    """Every ``run:`` script body in a parsed workflow.

    Args:
        workflow: A parsed GitHub Actions workflow document.

    Returns:
        list[str]: The shell scripts, one per step that declares ``run:``.
    """
    scripts: list[str] = []
    jobs = workflow.get("jobs")
    if not isinstance(jobs, dict):
        return scripts
    for job in jobs.values():
        steps = job.get("steps") if isinstance(job, dict) else None
        if not isinstance(steps, list):
            continue
        scripts.extend(step["run"] for step in steps if isinstance(step, dict) and isinstance(step.get("run"), str))
    return scripts


def _executes_make_target(script: str, target: str) -> bool:
    """Does this shell script actually *run* ``make <target>``?

    Splits the script into commands, drops comments and transparent prefixes, then
    requires a command whose executable is ``make`` and whose arguments include the
    target as a token in its own right. ``echo 'make test-integration'`` therefore does
    not count: the executable is ``echo`` and the make invocation is a quoted string it
    prints. Neither does a commented-out line.

    Args:
        script: The body of a workflow step's ``run:`` key.
        target: The make target that must be invoked.

    Returns:
        bool: True if the script invokes the target.
    """
    for line in script.splitlines():
        for segment in re.split(r"&&|\|\||;|\||\bthen\b|\belse\b|\bdo\b", line):
            try:
                tokens = shlex.split(segment, comments=True)
            except ValueError:
                continue
            # Drop leading `VAR=value` assignments and prefixes that just wrap a command.
            while tokens and (
                re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", tokens[0]) or tokens[0] in _TRANSPARENT_PREFIXES
            ):
                tokens = tokens[1:]
            if tokens and PurePosixPath(tokens[0]).name == "make" and target in tokens[1:]:
                return True
    return False


def test_a_scheduled_workflow_runs_the_integration_tier() -> None:
    """Some workflow must run the integration tier on a timer.

    The tier's dependencies - the Zenodo archive, the UCSC/NCBI/ENSEMBL reference URLs,
    the vendored Kestrel JAR - rot on their own schedule, not on the repository's. Until
    a cron looked for it, the first person to notice was whoever opened the next PR.

    The workflow is parsed rather than grepped, and the step's shell is parsed rather
    than substring-matched. A grep for ``make test-integration`` is satisfied by a
    comment, by a job name, or by ``echo 'make test-integration'`` - none of which run
    anything, which is precisely the failure this test exists to make impossible.

    Raises:
        AssertionError: If no scheduled workflow actually executes the tier.
    """
    if not WORKFLOWS.is_dir():
        pytest.skip("no GitHub Actions workflows present in this tree")

    paths = sorted(set(WORKFLOWS.glob("*.yml")) | set(WORKFLOWS.glob("*.yaml")))
    scheduled_with_integration = []
    for path in paths:
        try:
            workflow = yaml.safe_load(path.read_text(encoding="utf-8"))
        except yaml.YAMLError:
            continue
        if not isinstance(workflow, dict) or not _is_cron_scheduled(workflow):
            continue
        if any(_executes_make_target(script, "test-integration") for script in _run_scripts(workflow)):
            scheduled_with_integration.append(path.name)

    assert scheduled_with_integration, (
        "No workflow *executes* `make test-integration` in a step of a workflow with a "
        "`schedule:` cron trigger. ci-tests.yml gates on `-m unit` only, so without this "
        "the pipeline is never exercised end to end by CI at all. Note that mentioning "
        "the command - in a comment, a step name, or an `echo` - does not run it."
    )
