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
KESTREL_MULTIFILE_MODULE = TESTS_ROOT / "integration" / "test_kestrel_multifile.py"
ROOT_CONFTEST = TESTS_ROOT / "conftest.py"
DATA_UTILS = TESTS_ROOT / "support" / "data_utils.py"
WORKFLOWS = REPO_ROOT / ".github" / "workflows"
MAKEFILE = REPO_ROOT / "Makefile"


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


#: Statements after which the rest of their block cannot run.
_TERMINATING_STATEMENTS = (ast.Return, ast.Raise, ast.Break, ast.Continue)

#: Bare calls that end the test the same way a ``return`` does.
_TERMINATING_CALLS = frozenset({"exit", "skip", "xfail"})


def _terminates_the_block(statement: ast.stmt) -> bool:
    """Does this statement stop the statements after it in its block from running?

    Args:
        statement: The statement to classify.

    Returns:
        bool: True for ``return``/``raise``/``break``/``continue`` and for a bare
            ``sys.exit(...)``, ``pytest.skip(...)`` or ``pytest.xfail(...)``.
    """
    if isinstance(statement, _TERMINATING_STATEMENTS):
        return True
    if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Call):
        func = statement.value.func
        called = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
        return called in _TERMINATING_CALLS
    return False


def _is_constant(test: ast.expr, truth: bool) -> bool:
    """Is this condition a literal that always evaluates to ``truth``?"""
    return isinstance(test, ast.Constant) and bool(test.value) is truth


def _unreachable(tree: ast.AST) -> set[int]:
    """Identify every AST node that cannot run.

    Two shapes are recognised, and they are the two that make the delegation guard
    satisfiable without delegating anything: a statement placed after an unconditional
    ``return``/``raise``/``sys.exit()`` in the same block, and the dead branch of an
    ``if``/``while`` whose condition is a literal.

    Args:
        tree: The parsed module.

    Returns:
        set[int]: ``id()`` of every node inside a region that never executes.
    """
    dead: set[int] = set()

    def kill(node: ast.AST) -> None:
        for child in ast.walk(node):
            dead.add(id(child))

    for node in ast.walk(tree):
        for field in ("body", "orelse", "finalbody"):
            block = getattr(node, field, None)
            if not isinstance(block, list):
                continue
            terminated = False
            for statement in block:
                if terminated:
                    kill(statement)
                elif isinstance(statement, ast.stmt) and _terminates_the_block(statement):
                    terminated = True
        if isinstance(node, ast.If | ast.While):
            if _is_constant(node.test, False):
                for statement in node.body:
                    kill(statement)
            elif _is_constant(node.test, True):
                for statement in node.orelse:
                    kill(statement)

    return dead


def _call_nodes(path: Path, name: str) -> list[ast.Call]:
    """Every *reachable* invocation of ``name`` in a module, as AST call nodes.

    Matched on ``ast.Call``, so a bare name reference (``_ = run_bam_test_case``), an
    import, or a mention in a docstring does not count - only a real invocation. Both
    ``f(...)`` and ``mod.f(...)`` forms match, keyed on the final attribute name.

    Calls that cannot run are excluded, which is the difference between "the module
    mentions the shared validator" and "the module uses it". A call parked under
    ``if False:`` or after a ``return`` satisfies the plain AST test while validating
    exactly as much as a comment does.

    Args:
        path: Module to inspect.
        name: Function name to look for.

    Returns:
        list[ast.Call]: The matching call nodes, in traversal order.
    """
    tree = _parse(path)
    dead = _unreachable(tree)
    return [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and id(node) not in dead
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


def test_real_kestrel_multifile_proof_is_integration_only() -> None:
    """The real Java/samtools proof must never leak into the pure unit tier."""
    tree = _parse(KESTREL_MULTIFILE_MODULE)
    marker_attributes = {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Attribute)
        and isinstance(node.value.value, ast.Name)
        and node.value.value.id == "pytest"
        and node.value.attr == "mark"
    }

    assert "integration" in marker_attributes
    assert "unit" not in marker_attributes


def test_real_kestrel_multifile_proof_uses_registered_data() -> None:
    """The real samples must come through the canonical manifest and data validator."""
    tree = _parse(KESTREL_MULTIFILE_MODULE)
    tests = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name.startswith("test_")]
    assert tests, "The multifile module contains no executable test."
    assert all({"test_config", "ensure_test_data"} <= {arg.arg for arg in node.args.args} for node in tests)

    string_literals = {
        node.value for node in ast.walk(tree) if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    assert "file_resources" in string_literals
    assert "example_b178_hg19_subset.bam" in string_literals
    assert "example_40cf_hg38_subset.bam" in string_literals
    assert "tests/data/example_b178_hg19_subset.bam" not in string_literals
    assert "tests/data/example_40cf_hg38_subset.bam" not in string_literals


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


#: make options that read the makefile and report, without running a single recipe.
#: ``make -n test-integration`` prints the commands the tier would run and exits 0, so a
#: step carrying one satisfies "the workflow mentions the target" while executing nothing.
_MAKE_NON_EXECUTING_FLAGS = frozenset(
    {
        "-n",
        "--dry-run",
        "--just-print",
        "--recon",
        "-q",
        "--question",
        "-t",
        "--touch",
        "-p",
        "--print-data-base",
    }
)

#: A line that is nothing but ``exit`` or ``exit <status>``. Anything after it in the
#: same step is never reached; a conditional exit (``if ...; then exit 0; fi``) is not
#: this shape, so it does not disable the rest of the script.
_UNCONDITIONAL_EXIT = re.compile(r"^\s*exit(\s+\d+)?\s*$")


def _executes_make_target(script: str, target: str) -> bool:
    """Does this shell script actually *run* ``make <target>``?

    Splits the script into commands, drops comments and transparent prefixes, then
    requires a command whose executable is ``make`` and whose arguments include the
    target as a token in its own right. ``echo 'make test-integration'`` therefore does
    not count: the executable is ``echo`` and the make invocation is a quoted string it
    prints. Neither does a commented-out line.

    Two further ways to satisfy the shape without running anything are rejected here,
    because a guard that accepts either is a guard that can be turned off without going
    red:

    * ``make -n test-integration`` (and ``-q``, ``-t``, ``-p``) reads the makefile,
      reports, and runs no recipe. It exits 0, so the step is green and the tier never
      ran. Any of :data:`_MAKE_NON_EXECUTING_FLAGS` disqualifies the invocation.
    * a ``make`` call placed after a line that is nothing but ``exit`` is unreachable.
      Scanning stops there. A *conditional* exit is not this shape and is left alone.

    Args:
        script: The body of a workflow step's ``run:`` key.
        target: The make target that must be invoked.

    Returns:
        bool: True if the script invokes the target in a way that runs its recipe.
    """
    for line in script.splitlines():
        if _UNCONDITIONAL_EXIT.match(line):
            return False
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
            if not tokens or PurePosixPath(tokens[0]).name != "make" or target not in tokens[1:]:
                continue
            if any(argument in _MAKE_NON_EXECUTING_FLAGS for argument in tokens[1:]):
                continue
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


# ---------------------------------------------------------------------------
# The guards' own guards
#
# Both helpers above answer "is this really happening?" about a file in the repository,
# so a helper that is too generous does not fail - it silently stops checking. These
# pin the cases each one was accepting.
# ---------------------------------------------------------------------------


def _write_module(tmp_path: Path, source: str) -> Path:
    """Write a throwaway module and return its path."""
    path = tmp_path / "sample.py"
    path.write_text(source, encoding="utf-8")
    return path


def test_a_call_that_runs_is_found(tmp_path) -> None:
    """The ordinary case, so the reachability filter cannot pass by rejecting everything."""
    path = _write_module(tmp_path, "def test_x():\n    run_bam_test_case(case, runner)\n")

    assert len(_call_nodes(path, "run_bam_test_case")) == 1


@pytest.mark.parametrize(
    ("source", "reason"),
    [
        ("def test_x():\n    return\n    run_bam_test_case(case)\n", "after an unconditional return"),
        ("def test_x():\n    raise SystemExit\n    run_bam_test_case(case)\n", "after a raise"),
        ("def test_x():\n    pytest.skip('later')\n    run_bam_test_case(case)\n", "after pytest.skip"),
        ("def test_x():\n    if False:\n        run_bam_test_case(case)\n", "in a literally false branch"),
        ("def test_x():\n    if True:\n        pass\n    else:\n        run_bam_test_case(case)\n", "in a dead else"),
    ],
)
def test_a_call_that_cannot_run_is_not_counted(tmp_path, source, reason) -> None:
    """
    Delegation that never executes is not delegation.

    The AST guard was satisfied by the *presence* of a call, so each of these turned the
    shared Kestrel validation off while leaving the guard green.
    """
    path = _write_module(tmp_path, source)

    assert _call_nodes(path, "run_bam_test_case") == [], reason


def test_a_plain_make_invocation_executes_the_target() -> None:
    """The ordinary case: the guard must still recognise a real invocation."""
    assert _executes_make_target("make test-integration", "test-integration")
    assert _executes_make_target("cd repo && sudo make -j4 test-integration", "test-integration")


@pytest.mark.parametrize("target_name", ["test-integration", "test-integration-parallel"])
def test_the_integration_target_materializes_its_declared_derived_fixtures(target_name: str) -> None:
    """A clean scheduled run must not depend on ignored CRAM and single-end files."""
    target = next(
        line for line in MAKEFILE.read_text(encoding="utf-8").splitlines() if line.startswith(f"{target_name}:")
    )

    assert target.split()[1:] == ["cram-fixtures"]


@pytest.mark.parametrize(
    "flag",
    ["-n", "--dry-run", "--just-print", "--recon", "-q", "--question", "-t", "--touch"],
)
def test_a_make_invocation_that_runs_no_recipe_does_not_count(flag) -> None:
    """
    `make -n test-integration` prints what it would do, runs nothing, and exits 0.

    The guard token-matched the target, so any of these switched the scheduled tier off
    while the step stayed green and the guard stayed satisfied.
    """
    assert not _executes_make_target(f"make {flag} test-integration", "test-integration")


def test_a_make_invocation_after_an_unconditional_exit_does_not_count() -> None:
    """A step that has already exited 0 never reaches the line below it."""
    script = "echo skipping for now\nexit 0\nmake test-integration\n"

    assert not _executes_make_target(script, "test-integration")


def test_a_conditional_exit_does_not_disable_the_rest_of_the_script() -> None:
    """`if ...; then exit 0; fi` is an early return for one case, not a dead script."""
    script = 'if [ -z "$DATA" ]; then exit 0; fi\nmake test-integration\n'

    assert _executes_make_target(script, "test-integration")
