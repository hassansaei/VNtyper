"""
Unit tests for the advisory mutation harness, ``scripts/mutation_test.py``.

The harness rewrites production source in place, so two of its properties are worth
pinning hard: it must only ever produce mutants that compile, and it must put the
original file back. Both are asserted below against a throwaway module in ``tmp_path``;
nothing here touches a real ``vntyper`` module.

The third property is the one that made two earlier sweeps fictional and is the reason
this file exists at all: byte-length-preserving mutants plus CPython's ``(mtime, size)``
``.pyc`` validation means a stale cache silently runs the *unmutated* code. The defence
is ``PYTHONDONTWRITEBYTECODE=1`` in the child plus deleting every ``__pycache__``, and
both halves are pinned here so a future edit cannot quietly drop one.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import mutation_test  # noqa: E402

# ---------------------------------------------------------------------------
# Numeric literal mutation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("literal", "expected"),
    [("0", "1"), ("1", "2"), ("41", "42"), ("-1", "0")],
)
def test_an_integer_literal_is_incremented(literal, expected) -> None:
    """Off-by-one on a threshold is the defect this operator exists to model."""
    assert mutation_test.mutate_number(literal) == expected


def test_a_float_literal_is_incremented() -> None:
    """Floats round-trip through ``repr`` so the mutant is still a valid literal."""
    assert mutation_test.mutate_number("1.5") == "2.5"


@pytest.mark.parametrize("literal", ["0x1f", "0b1010", "0o17", "1_000", "3j"])
def test_literals_that_cannot_be_incremented_safely_are_skipped(literal) -> None:
    """
    Hex, binary, octal, underscored and imaginary literals are declined, not guessed.

    Rewriting ``0x1f`` as ``32`` would still compile, so nothing downstream would
    reject it - it would just silently stop being the literal the author wrote.
    """
    assert mutation_test.mutate_number(literal) is None


# ---------------------------------------------------------------------------
# Mutant generation
# ---------------------------------------------------------------------------


def _write_module(tmp_path: Path, source: str) -> Path:
    """Write a throwaway module and return its path."""
    path = tmp_path / "sample.py"
    path.write_text(source, encoding="utf-8")
    return path


def test_every_generated_mutant_compiles(tmp_path, monkeypatch) -> None:
    """
    The headline safety property: a syntax error is not a defect a test could catch.

    ``*args`` -> ``/args`` is the case that makes this non-trivial - the ``*`` is a
    mutable operator token, but the result does not parse, so it must be dropped rather
    than counted as an uncaught mutant.
    """
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    path = _write_module(
        tmp_path,
        "def f(*args, **kwargs):\n    return len(args) + len(kwargs)\n",
    )

    mutants = mutation_test.generate_mutants(path)

    assert mutants, "the `+` in the return should still yield a mutant"
    for mutant in mutants:
        compile(mutant.source, "sample.py", "exec")


def test_comparison_and_boolean_operators_are_mutated(tmp_path, monkeypatch) -> None:
    """The decision-flipping operators are the point of the exercise."""
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    path = _write_module(tmp_path, "def f(a, b):\n    return a >= 1 and b\n")

    swaps = {(m.original, m.replacement) for m in mutation_test.generate_mutants(path)}

    assert (">=", "<") in swaps
    assert ("and", "or") in swaps
    assert ("1", "2") in swaps


def test_not_is_deleted_rather_than_replaced(tmp_path, monkeypatch) -> None:
    """
    Removing ``not`` inverts a guard while leaving the statement's shape intact.

    That is the highest-signal operator for this codebase's "stages mark, they do not
    filter" logic, and the empty replacement has to still produce parseable source.
    """
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    path = _write_module(tmp_path, "def f(a):\n    if not a:\n        return 0\n    return 1\n")

    mutants = [m for m in mutation_test.generate_mutants(path) if m.original == "not"]

    assert len(mutants) == 1
    assert mutants[0].replacement == ""
    assert "if  a:" in mutants[0].source
    compile(mutants[0].source, "sample.py", "exec")


def test_strings_and_comments_are_never_mutated(tmp_path, monkeypatch) -> None:
    """
    A mutation inside a docstring or a comment is a guaranteed survivor and pure noise.

    Tokenizing gives this for free - both arrive as their own token types - and this
    test is what stops someone "simplifying" the generator into a regex over the source.
    """
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    path = _write_module(
        tmp_path,
        '"""Docstring with and, or, True and 1 < 2."""\n# comment with and or True 1 < 2\nX = 5\n',
    )

    mutants = mutation_test.generate_mutants(path)

    assert {(m.original, m.replacement) for m in mutants} == {("5", "6")}


def test_the_reported_line_number_is_the_mutated_token_s_line(tmp_path, monkeypatch) -> None:
    """The survivor list is only actionable if the line numbers are right."""
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    path = _write_module(tmp_path, "A = 1\n\n\nB = 2 and 3\n")

    lines = {m.line for m in mutation_test.generate_mutants(path) if m.original == "and"}

    assert lines == {4}


# ---------------------------------------------------------------------------
# The .pyc defence
# ---------------------------------------------------------------------------


def test_clear_bytecode_caches_removes_every_pycache(tmp_path) -> None:
    """
    Half the ``.pyc`` defence: a cache left by an earlier run shadows the mutant.

    ``PYTHONDONTWRITEBYTECODE`` only stops caches being *written*. A ``.pyc`` already on
    disk whose recorded (mtime, size) still matches is loaded regardless, which is
    exactly the stale-cache failure this deletes its way out of.
    """
    (tmp_path / "pkg" / "__pycache__").mkdir(parents=True)
    (tmp_path / "pkg" / "__pycache__" / "mod.pyc").write_bytes(b"stale")
    (tmp_path / "pkg" / "sub" / "__pycache__").mkdir(parents=True)
    (tmp_path / "pkg" / "keep.py").write_text("X = 1\n", encoding="utf-8")

    removed = mutation_test.clear_bytecode_caches(tmp_path)

    assert removed == 2
    assert not list(tmp_path.rglob("__pycache__"))
    assert (tmp_path / "pkg" / "keep.py").exists(), "only caches are removed"


def test_the_test_subprocess_disables_bytecode_writing(monkeypatch) -> None:
    """
    The other half of the defence, pinned so it cannot be dropped in a refactor.

    Without both ``-B`` and ``PYTHONDONTWRITEBYTECODE=1`` a byte-length-preserving
    mutant written in the same second as the file it replaces is invisible to the cache
    validator, the interpreter runs the unmutated bytecode, and every mutant "survives".
    """
    captured: dict = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["env"] = kwargs["env"]
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(mutation_test.subprocess, "run", fake_run)

    mutation_test.run_tests(("tests/unit/test_scoring.py",), timeout=60)

    assert captured["env"]["PYTHONDONTWRITEBYTECODE"] == "1"
    assert "-B" in captured["command"]


def test_a_timed_out_test_run_counts_as_a_kill(monkeypatch) -> None:
    """A mutant that hangs the suite has been detected, not missed."""

    def fake_run(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(mutation_test.subprocess, "run", fake_run)

    assert mutation_test.run_tests(("tests/unit",), timeout=1) is False


# ---------------------------------------------------------------------------
# Sweeping restores the file
# ---------------------------------------------------------------------------


def test_the_original_file_is_restored_after_a_sweep(tmp_path, monkeypatch) -> None:
    """
    The harness writes over real production source, so this is non-negotiable.

    Every mutant here is reported as killed, so the sweep runs its full loop and still
    has to leave the file byte-identical.
    """
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(mutation_test, "run_tests", lambda *_, **__: False)
    original = "def f(a):\n    return a >= 1 and a\n"
    _write_module(tmp_path, original)

    result = mutation_test.sweep_module("sample.py", ("tests",), timeout=5, verbose=False)

    assert result.killed == result.total > 0
    assert (tmp_path / "sample.py").read_text(encoding="utf-8") == original


def test_the_original_file_is_restored_when_a_sweep_raises(tmp_path, monkeypatch) -> None:
    """
    An interrupted sweep must not leave a mutated module on disk.

    This is the path that matters in practice - a maintainer pressing Ctrl-C - and
    without the ``finally`` it would leave production source silently wrong in a
    working tree someone is about to commit from.
    """
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)

    def explode(*_args, **_kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(mutation_test, "run_tests", explode)
    original = "def f(a):\n    return a >= 1\n"
    _write_module(tmp_path, original)

    with pytest.raises(KeyboardInterrupt):
        mutation_test.sweep_module("sample.py", ("tests",), timeout=5, verbose=False)

    assert (tmp_path / "sample.py").read_text(encoding="utf-8") == original


def test_a_surviving_mutant_is_confirmed_against_the_full_tier(tmp_path, monkeypatch) -> None:
    """
    Scoping the tests would bias the score down; the escalation is what removes the bias.

    A mutant that the module's own tests miss is re-run against the whole unit tier
    before it is recorded as a survivor, so "survived" means survived everything.
    """
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    seen: list[tuple[tuple[str, ...], bool]] = []

    def fake_run_tests(paths, _timeout, parallel=False):
        seen.append((paths, parallel))
        # Pass the scoped run, fail the full tier: the mutant is killed on escalation.
        return paths != ("tests/unit",)

    monkeypatch.setattr(mutation_test, "run_tests", fake_run_tests)
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")

    result = mutation_test.sweep_module("sample.py", ("scoped.py",), timeout=5, verbose=False)

    assert (("tests/unit",), True) in seen, "a scoped survivor must be escalated, in parallel"
    assert (("scoped.py",), False) in seen, "the scoped run comes first and is serial"
    assert result.killed > 0
    assert result.survived == 0, "the full tier killed it, so it is not a survivor"
    assert result.survivors == []


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def test_a_module_with_no_mutants_scores_zero_rather_than_dividing_by_zero() -> None:
    """``--module`` can select a module that yields nothing; that must not raise."""
    assert mutation_test.ModuleResult(path="x.py").score == 0.0
    assert mutation_test.ModuleResult(path="x.py").total == 0


def test_the_report_names_every_surviving_mutant() -> None:
    """
    The survivor list is the actionable half of the artefact.

    A score alone says the suite is weak; the named survivors say where.
    """
    survivor = mutation_test.Mutant(
        path="vntyper/scripts/scoring.py",
        line=42,
        original=">=",
        replacement="<",
        source="",
        killed=False,
    )
    result = mutation_test.ModuleResult(path="vntyper/scripts/scoring.py", killed=3, survived=1, survivors=[survivor])

    report = mutation_test.format_report([result], elapsed=60.0)

    assert "75.0%" in report
    assert "line   42" in report
    assert "'>=' -> '<'" in report


def _survivor(line: int = 42, path: str = "vntyper/scripts/scoring.py") -> "mutation_test.Mutant":
    """Build a surviving mutant for the classification tests."""
    return mutation_test.Mutant(path=path, line=line, original=">=", replacement="<", source="", killed=False)


def test_an_unclassified_survivor_is_reported_as_a_genuine_gap() -> None:
    """The default must be "this is a real gap", never "probably equivalent"."""
    survivor = _survivor()
    result = mutation_test.ModuleResult(path=survivor.path, killed=3, survived=1, survivors=[survivor])

    page = mutation_test.format_markdown([result], elapsed=1.0)

    assert survivor.equivalence_reason is None
    assert "### Genuine gaps" in page
    assert "None - every survivor is classified equivalent below." not in page


def test_a_classified_survivor_is_separated_out_and_its_reason_is_published(monkeypatch) -> None:
    """
    An equivalent mutant is unkillable, so counting it against the suite is misleading.

    It is only defensible if the reason travels with it, which is why the reason is
    rendered into the page rather than kept in someone's head.
    """
    survivor = _survivor()
    monkeypatch.setitem(mutation_test.EQUIVALENT_MUTANTS, survivor.key, "the config always supplies this default")
    result = mutation_test.ModuleResult(path=survivor.path, killed=3, survived=1, survivors=[survivor])

    page = mutation_test.format_markdown([result], elapsed=1.0)

    assert survivor.equivalence_reason == "the config always supplies this default"
    assert "the config always supplies this default" in page
    # 3 killed of 4 raw = 75%; excluding the one equivalent mutant, 3 of 3 = 100%.
    assert "75.0%" in page
    assert "100.0%" in page


def test_both_the_raw_and_the_adjusted_score_are_published(monkeypatch) -> None:
    """
    Publishing only one number would overstate the result in one direction or the other.

    The raw score counts unkillable mutants against the suite; the adjusted score rests
    on human judgement. The page has to show both and say which is which.
    """
    survivor = _survivor()
    monkeypatch.setitem(mutation_test.EQUIVALENT_MUTANTS, survivor.key, "unreachable bound")
    result = mutation_test.ModuleResult(path=survivor.path, killed=1, survived=1, survivors=[survivor])

    page = mutation_test.format_markdown([result], elapsed=1.0)

    assert "raw mutation score of 50.0%" in page
    assert "**100.0%**" in page


def test_a_stale_equivalence_entry_is_surfaced_rather_than_ignored(monkeypatch) -> None:
    """
    A classification that no longer matches anything is excusing a mutant that is gone.

    Line numbers shift when a module is edited, so entries rot silently. Left unreported,
    a stale entry becomes a permanent excuse attached to whatever now sits on that line.
    """
    stale_key = ("vntyper/scripts/scoring.py", 999, ">=", "<")
    monkeypatch.setitem(mutation_test.EQUIVALENT_MUTANTS, stale_key, "no longer present")
    result = mutation_test.ModuleResult(path="vntyper/scripts/scoring.py", killed=1, survived=0)

    assert mutation_test.stale_equivalence_keys([result]) == [stale_key]
    assert "STALE EQUIVALENCE ENTRIES" in mutation_test.format_report([result], elapsed=1.0)


def test_a_scoped_run_does_not_report_other_modules_entries_as_stale(monkeypatch) -> None:
    """`--module scoring.py` must not declare every other module's classification rotten."""
    other_key = ("vntyper/scripts/flagging.py", 10, "and", "or")
    monkeypatch.setitem(mutation_test.EQUIVALENT_MUTANTS, other_key, "not swept in this run")
    result = mutation_test.ModuleResult(path="vntyper/scripts/scoring.py", killed=1, survived=0)

    assert mutation_test.stale_equivalence_keys([result]) == []


def test_every_committed_equivalence_entry_targets_a_real_module() -> None:
    """A typo'd path in the classification table would silently never apply."""
    for path, _line, _original, _replacement in mutation_test.EQUIVALENT_MUTANTS:
        assert path in mutation_test.TARGETS, path


def test_the_page_records_that_branch_coverage_is_measured_but_not_enabled() -> None:
    """
    Branch coverage was measured at 63.80% against a floor of 66 and left disabled.

    That is a real finding with a real decision attached, and the next person to reach
    for `branch = true` needs to find it rather than re-derive it.
    """
    result = mutation_test.ModuleResult(path="vntyper/scripts/scoring.py", killed=4, survived=0)

    page = mutation_test.format_markdown([result], elapsed=1.0)

    assert "63.80%" in page
    assert "not" in page and "lowered" in page
    assert "144 more covered units" in page


def test_the_page_warns_against_building_from_the_tree_during_a_sweep() -> None:
    """
    The harness rewrites `vntyper/scripts/*.py` in place, so a build started mid-sweep
    bakes a live mutant into its artefact.

    This is not hypothetical: a Docker image built during a sweep crashed in the
    container with a pandas `KeyError` in `motif_processing.py`, which reads exactly
    like a production bug and cost a diagnosis cycle before it was traced back. The
    `finally` restore protects the repository, not anything already built from it, so
    the warning has to be on the page next to the `.pyc` one - both are "the result you
    got is not the result you think you got".
    """
    result = mutation_test.ModuleResult(path="vntyper/scripts/scoring.py", killed=4, survived=0)

    page = mutation_test.format_markdown([result], elapsed=1.0)

    assert "in place" in page
    assert "docker build" in page.lower()
    assert "git diff --quiet -- vntyper/" in page


def test_the_page_says_the_43_5_percent_baseline_is_not_directly_comparable() -> None:
    """
    Different mutant population, different modules - the totals do not line up.

    Presenting them as a like-for-like improvement would be exactly the overstatement
    this artefact exists to avoid.
    """
    result = mutation_test.ModuleResult(path="vntyper/scripts/confidence_assignment.py", killed=4, survived=0)

    page = mutation_test.format_markdown([result], elapsed=1.0)

    assert "43.5%" in page
    assert "not directly comparable" in page


def test_the_markdown_page_embeds_the_raw_output_and_the_pyc_warning() -> None:
    """
    The committed docs page has to be the evidence, not a summary of it.

    It also has to carry the ``.pyc`` warning: the page is where the next person looks
    before re-running the sweep, and that trap cost two fictional sweeps.
    """
    result = mutation_test.ModuleResult(path="vntyper/scripts/scoring.py", killed=4, survived=0)

    page = mutation_test.format_markdown([result], elapsed=120.0)

    assert page.startswith("# Mutation Testing")
    assert "PYTHONDONTWRITEBYTECODE=1" in page
    assert "```text" in page, "the raw console report is embedded verbatim"
    assert mutation_test.format_report([result], elapsed=120.0) in page


def test_every_declared_target_and_its_tests_exist() -> None:
    """
    A typo in TARGETS would silently mutate nothing, or score against the wrong tests.

    Both failures look like a healthy run, so the paths are checked against the repo.
    """
    for module, test_paths in mutation_test.TARGETS.items():
        assert (mutation_test.REPO_ROOT / module).is_file(), module
        for test_path in test_paths:
            assert (mutation_test.REPO_ROOT / test_path).is_file(), test_path


# ---------------------------------------------------------------------------
# The clean-tree preflight
# ---------------------------------------------------------------------------

#: Commands whose effect is to throw the working tree away. None of them may appear in
#: anything the harness prints: the one situation where a user reads these messages is
#: the one where their own uncommitted work is on the line.
DESTRUCTIVE_ADVICE = ("git checkout --", "git checkout -- ", "git restore", "git reset --hard", "git clean")


def _fake_subprocess_run(git_outputs: list[str], calls: list[list[str]]):
    """
    Build a ``subprocess.run`` double that scripts ``git`` and fails every pytest run.

    Args:
        git_outputs (list[str]): ``git status --porcelain`` stdout, one entry per call.
        calls (list[list[str]]): Mutated in place with every command issued.

    Returns:
        Callable: A drop-in replacement for ``subprocess.run``.
    """
    pending = list(git_outputs)

    def fake_run(command, *_args, **_kwargs):
        calls.append(list(command))
        if command[0] == "git":
            return subprocess.CompletedProcess(command, 0, stdout=pending.pop(0), stderr="")
        # A non-zero pytest exit means "mutant killed", so sweeps finish immediately.
        return subprocess.CompletedProcess(command, 1, stdout="", stderr="")

    return fake_run


def _prepare_main(tmp_path, monkeypatch, git_outputs: list[str], argv: list[str] | None = None) -> list[list[str]]:
    """Point ``main()`` at ``tmp_path``, script ``git``, and neutralise signal setup."""
    monkeypatch.setattr(mutation_test, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(mutation_test, "TARGETS", {"sample.py": ("tests/unit/x.py",)})
    monkeypatch.setattr(mutation_test.signal, "signal", lambda *_args: None)
    monkeypatch.setattr(sys, "argv", argv or ["mutation_test.py"])
    calls: list[list[str]] = []
    monkeypatch.setattr(subprocess, "run", _fake_subprocess_run(git_outputs, calls))
    return calls


def test_a_dirty_target_file_stops_the_sweep_before_anything_is_rewritten(tmp_path, monkeypatch, capsys) -> None:
    """
    The harness restores the text it read at the start, not the committed text.

    So a sweep begun on a modified target measures *uncommitted* code and publishes it
    as the committed baseline, and then reports the maintainer's own edit as a failed
    restoration. Refusing up front is the only honest outcome.
    """
    original = "def f(a):\n    return a >= 1  # uncommitted work\n"
    _write_module(tmp_path, original)
    _prepare_main(tmp_path, monkeypatch, [" M sample.py\n"])
    swept: list[str] = []
    monkeypatch.setattr(mutation_test, "sweep_module", lambda path_str, *_a, **_k: swept.append(path_str))

    exit_code = mutation_test.main()

    assert exit_code != 0, "a sweep that never ran must not report success"
    assert swept == [], "the sweep must not rewrite a file that holds uncommitted work"
    assert (tmp_path / "sample.py").read_text(encoding="utf-8") == original
    assert "sample.py" in capsys.readouterr().err, "the message must name the dirty file"


def test_the_preflight_refusal_never_tells_the_user_to_discard_their_work(tmp_path, monkeypatch, capsys) -> None:
    """``git checkout --`` on a file the user just edited destroys the edit, silently."""
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    _prepare_main(tmp_path, monkeypatch, [" M sample.py\n"])
    monkeypatch.setattr(
        mutation_test,
        "sweep_module",
        lambda path_str, *_a, **_k: mutation_test.ModuleResult(path=path_str, killed=1),
    )

    mutation_test.main()

    err = capsys.readouterr().err
    for advice in DESTRUCTIVE_ADVICE:
        assert advice not in err, f"the refusal must not recommend {advice!r}"


def test_the_unrestored_file_error_never_tells_the_user_to_discard_their_work(tmp_path, monkeypatch, capsys) -> None:
    """
    The end-of-run check fires on any diff, and cannot tell a mutant from a real edit.

    It runs after a clean preflight here - the genuine "restore failed" case - and even
    then it must not name a command that throws work away.
    """
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    _prepare_main(tmp_path, monkeypatch, ["", " M sample.py\n"])

    exit_code = mutation_test.main()

    assert exit_code != 0
    err = capsys.readouterr().err
    assert "sample.py" in err
    for advice in DESTRUCTIVE_ADVICE:
        assert advice not in err, f"the restore failure must not recommend {advice!r}"


def test_a_clean_target_tree_is_swept_normally(tmp_path, monkeypatch) -> None:
    """The preflight must gate on uncommitted work only - a clean tree runs as before."""
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    _prepare_main(tmp_path, monkeypatch, ["", ""])

    exit_code = mutation_test.main()

    assert exit_code == 0


def test_the_preflight_checks_only_the_targets_the_run_will_rewrite(tmp_path, monkeypatch) -> None:
    """
    ``--module`` narrows the sweep, so it must narrow the guard too.

    Blocking on a file this run will never touch would make the harness unusable in any
    working tree with an edit in flight.
    """
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    calls = _prepare_main(
        tmp_path,
        monkeypatch,
        ["", ""],
        argv=["mutation_test.py", "--module", "sample.py"],
    )
    monkeypatch.setattr(mutation_test, "TARGETS", {"sample.py": ("tests/unit/x.py",), "other.py": ("tests/unit/y.py",)})

    mutation_test.main()

    git_calls = [c for c in calls if c[0] == "git"]
    assert git_calls, "the preflight must consult git before mutating anything"
    assert "sample.py" in git_calls[0]
    assert "other.py" not in git_calls[0], "an unselected target must not gate the run"


def test_the_preflight_also_guards_the_files_it_will_overwrite(tmp_path, monkeypatch) -> None:
    """
    `write_outputs()` rewrites the report and the raw results wholesale, every run.

    Guarding only the mutated sources meant an uncommitted edit to
    ``docs/development/mutation-testing.md`` was destroyed by the next `make mutation`
    with no warning at all.
    """
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    report = tmp_path / "docs/development/mutation-testing.md"
    results = tmp_path / "docs/development/mutation-results.json"
    calls = _prepare_main(
        tmp_path,
        monkeypatch,
        ["", ""],
        argv=["mutation_test.py", "--output", str(report), "--results-json", str(results)],
    )

    mutation_test.main()

    preflight = [c for c in calls if c[0] == "git"][0]
    assert "docs/development/mutation-testing.md" in preflight
    assert "docs/development/mutation-results.json" in preflight


def test_a_dirty_report_page_stops_the_run_before_it_is_overwritten(tmp_path, monkeypatch, capsys) -> None:
    """The docs page is generated, but an edit in flight is still somebody's work."""
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    report = tmp_path / "docs/development/mutation-testing.md"
    _prepare_main(
        tmp_path,
        monkeypatch,
        [" M docs/development/mutation-testing.md\n"],
        argv=["mutation_test.py", "--output", str(report)],
    )
    swept: list[str] = []
    monkeypatch.setattr(mutation_test, "sweep_module", lambda path_str, *_a, **_k: swept.append(path_str))

    exit_code = mutation_test.main()

    assert exit_code != 0
    assert swept == []
    assert "docs/development/mutation-testing.md" in capsys.readouterr().err


def test_render_only_runs_the_preflight_over_the_page_it_rewrites(tmp_path, monkeypatch, capsys) -> None:
    """
    ``--render-only`` mutates nothing, which is exactly why it looked safe to exempt.

    It still overwrites the docs page, so `make mutation-render` over an uncommitted edit
    destroyed it just as surely as a full sweep would have.
    """
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    saved = tmp_path / "results.json"
    saved.write_text(json.dumps({"elapsed": 1.0, "modules": []}) + "\n", encoding="utf-8")
    report = tmp_path / "docs/development/mutation-testing.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("# uncommitted edit\n", encoding="utf-8")
    _prepare_main(
        tmp_path,
        monkeypatch,
        [" M docs/development/mutation-testing.md\n"],
        argv=["mutation_test.py", "--render-only", str(saved), "--output", str(report)],
    )

    exit_code = mutation_test.main()

    assert exit_code != 0
    assert report.read_text(encoding="utf-8") == "# uncommitted edit\n"
    assert "docs/development/mutation-testing.md" in capsys.readouterr().err


def test_render_only_over_a_clean_page_still_renders(tmp_path, monkeypatch) -> None:
    """The guard must gate on uncommitted work only - the ordinary re-render is unchanged."""
    saved = tmp_path / "results.json"
    saved.write_text(json.dumps({"elapsed": 1.0, "modules": []}) + "\n", encoding="utf-8")
    report = tmp_path / "docs/development/mutation-testing.md"
    _prepare_main(
        tmp_path,
        monkeypatch,
        [""],
        argv=["mutation_test.py", "--render-only", str(saved), "--output", str(report)],
    )

    assert mutation_test.main() == 0
    assert report.exists()


def test_an_indeterminate_git_answer_stops_the_run(tmp_path, monkeypatch, capsys) -> None:
    """
    Unknown is not clean. The guard's own failure must not become permission to proceed.

    This is the fail-open pattern the branch exists to remove, and it had shipped inside
    the fix for it: `dirty_paths` returned "nothing is dirty" whenever git could not be
    run at all.
    """
    _write_module(tmp_path, "def f(a):\n    return a >= 1\n")
    _prepare_main(tmp_path, monkeypatch, [""])
    monkeypatch.setattr(subprocess, "run", _raise_missing_git)
    swept: list[str] = []
    monkeypatch.setattr(mutation_test, "sweep_module", lambda path_str, *_a, **_k: swept.append(path_str))

    exit_code = mutation_test.main()

    assert exit_code != 0, "a guard that cannot check must not wave the sweep through"
    assert swept == []
    assert "cannot determine" in capsys.readouterr().err


def _raise_missing_git(command, *_args, **_kwargs):
    """Stand in for a checkout with no ``git`` on PATH."""
    if command[0] == "git":
        raise FileNotFoundError("git")
    return subprocess.CompletedProcess(command, 1, stdout="", stderr="")
