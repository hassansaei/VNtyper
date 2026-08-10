"""Unit tests for the mutation harness's selected-target and output guards."""

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import mutation_guard  # noqa: E402

#: Every command whose effect is to throw uncommitted changes away.
DESTRUCTIVE = ("git checkout", "git restore", "git reset", "git clean", "git stash drop")


# ---------------------------------------------------------------------------
# Reading git's answer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        (" M vntyper/scripts/scoring.py", "vntyper/scripts/scoring.py"),
        ("M  vntyper/scripts/scoring.py", "vntyper/scripts/scoring.py"),
        ("?? vntyper/scripts/new.py", "vntyper/scripts/new.py"),
        ("MM vntyper/scripts/flagging.py", "vntyper/scripts/flagging.py"),
        ("A  vntyper/scripts/added.py", "vntyper/scripts/added.py"),
    ],
)
def test_every_porcelain_status_code_yields_its_path(line, expected) -> None:
    """Staged, unstaged, both and untracked all mean "do not sweep this file"."""
    assert mutation_guard.parse_porcelain(f"{line}\n") == [expected]


def test_a_rename_reports_the_path_that_is_on_disk() -> None:
    """``R  old -> new``: the file the sweep would overwrite is the destination."""
    assert mutation_guard.parse_porcelain("R  old.py -> vntyper/scripts/new.py\n") == ["vntyper/scripts/new.py"]


def test_a_quoted_path_is_unquoted() -> None:
    """git quotes paths containing special characters; the quotes are not the name."""
    assert mutation_guard.parse_porcelain(' M "vntyper/scripts/odd name.py"\n') == ["vntyper/scripts/odd name.py"]


@pytest.mark.parametrize("output", ["", "\n", "   \n\n"])
def test_a_clean_tree_reports_nothing(output) -> None:
    """A clean tree is the normal case and must never produce a phantom path."""
    assert mutation_guard.parse_porcelain(output) == []


def test_every_dirty_file_is_reported_not_just_the_first() -> None:
    """The refusal names files, so missing one would send someone hunting."""
    output = " M vntyper/scripts/scoring.py\n?? vntyper/scripts/flagging.py\n"

    assert mutation_guard.parse_porcelain(output) == [
        "vntyper/scripts/scoring.py",
        "vntyper/scripts/flagging.py",
    ]


def test_git_is_asked_only_about_the_paths_it_was_given(tmp_path, monkeypatch) -> None:
    """
    An unscoped ``git status`` would block the sweep on any unrelated edit in the tree.

    The pathspec after ``--`` is what keeps the guard proportionate: it refuses only
    when a file the sweep is about to rewrite is dirty.
    """
    seen: dict[str, object] = {}

    def fake_run(command, **kwargs):
        seen["command"] = command
        seen["cwd"] = kwargs.get("cwd")
        return subprocess.CompletedProcess(command, 0, stdout=" M a.py\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert mutation_guard.dirty_paths(tmp_path, ["a.py", "b.py"]) == ["a.py"]
    assert seen["command"] == ["git", "status", "--porcelain", "--", "a.py", "b.py"]
    assert seen["cwd"] == tmp_path


def test_a_failing_git_call_blocks_the_sweep(tmp_path, monkeypatch) -> None:
    """
    Non-zero from ``git status`` means *unknown*, not *clean* - so it refuses.

    This reverses the choice made when the guard was first written, which returned ``[]``
    here on the reasoning that "outside a repository the guard has nothing to protect".
    That reasoning was wrong twice over. It confuses "there is nothing to protect" with
    "I could not find out", and the two are indistinguishable from inside this function:
    a corrupted index, a `safe.directory` refusal and a genuine non-repository all exit
    non-zero. And a guard that waves the sweep through whenever its own check breaks is
    exactly the fail-open pattern the rest of this branch exists to remove - it would
    have shipped inside the fix for it.

    The cost of being wrong is asymmetric. Refusing wrongly prints a message and wastes a
    maintainer's minute; proceeding wrongly rewrites production source, publishes
    uncommitted code as the committed baseline, and overwrites the generated report.
    """

    def fake_run(command, **_kwargs):
        return subprocess.CompletedProcess(command, 128, stdout="", stderr="not a git repository")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="not a git repository"):
        mutation_guard.dirty_paths(tmp_path, ["a.py"])


def test_a_missing_git_binary_blocks_the_sweep(tmp_path, monkeypatch) -> None:
    """
    No git, no answer, no sweep. Also a reversal - see the test above for why.

    The guard was described as "a safety net, not a dependency". It is a dependency:
    without it the harness cannot establish the one precondition that makes its
    read-mutate-restore contract safe.
    """

    def fake_run(*_args, **_kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="git"):
        mutation_guard.dirty_paths(tmp_path, ["a.py"])


def test_the_refusal_to_guess_says_what_went_wrong(tmp_path, monkeypatch) -> None:
    """A bare "refusing" would leave the reader unable to tell a broken repo from a dirty one."""

    def fake_run(command, **_kwargs):
        return subprocess.CompletedProcess(command, 128, stdout="", stderr="detected dubious ownership")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(RuntimeError) as excinfo:
        mutation_guard.dirty_paths(tmp_path, ["a.py"])

    message = str(excinfo.value)
    assert "detected dubious ownership" in message
    assert "Nothing has been modified." in message


# ---------------------------------------------------------------------------
# What the run would overwrite
# ---------------------------------------------------------------------------


def test_the_generated_outputs_are_guarded_alongside_the_sources(tmp_path) -> None:
    """
    `write_outputs()` overwrites the docs page unconditionally, so it is a target too.

    Checking only the mutated sources left an uncommitted edit to
    ``docs/development/mutation-testing.md`` to be destroyed by the next `make mutation`.
    """
    guarded = mutation_guard.writable_paths(
        tmp_path,
        ["vntyper/scripts/scoring.py"],
        [tmp_path / "docs/development/mutation-testing.md", tmp_path / "docs/development/mutation-results.json"],
    )

    assert guarded == [
        "docs/development/mutation-results.json",
        "docs/development/mutation-testing.md",
        "vntyper/scripts/scoring.py",
    ]


def test_an_output_outside_the_repository_is_not_asked_about(tmp_path) -> None:
    """`--output /tmp/x.md` is legal; git has nothing to say about it, so it is dropped."""
    guarded = mutation_guard.writable_paths(tmp_path, ["a.py"], [Path("/tmp/elsewhere.md")])

    assert guarded == ["a.py"]


def test_an_unrequested_output_contributes_nothing(tmp_path) -> None:
    """`--results-json` is optional, and `--render-only` writes no JSON at all."""
    assert mutation_guard.writable_paths(tmp_path, ["a.py"], [None, None]) == ["a.py"]


def test_a_path_is_guarded_once_even_if_named_twice(tmp_path) -> None:
    """A duplicate would name the same file twice in the refusal."""
    output = tmp_path / "a.py"

    assert mutation_guard.writable_paths(tmp_path, ["a.py"], [output, output]) == ["a.py"]


# ---------------------------------------------------------------------------
# Wording the two failures
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "formatter",
    [
        mutation_guard.format_dirty_tree_refusal,
        mutation_guard.format_unrestored_warning,
        mutation_guard.format_indeterminate_refusal,
    ],
)
def test_no_message_recommends_a_command_that_discards_work(formatter) -> None:
    """
    This is the finding that produced this module, and it is the regression to hold.

    The harness used to answer a dirty target with "Run `git checkout --` on: ...",
    which deletes the uncommitted work it had just mistaken for a failed restore.
    """
    message = formatter(["vntyper/scripts/scoring.py"])

    for command in DESTRUCTIVE:
        assert command not in message, f"{formatter.__name__} must not recommend {command!r}"


@pytest.mark.parametrize(
    "formatter",
    [mutation_guard.format_dirty_tree_refusal, mutation_guard.format_unrestored_warning],
)
def test_every_affected_file_is_named(formatter) -> None:
    """A message that says "some file" sends the reader to `git status` to find out."""
    message = formatter(["vntyper/scripts/scoring.py", "vntyper/scripts/flagging.py"])

    assert "vntyper/scripts/scoring.py" in message
    assert "vntyper/scripts/flagging.py" in message


def test_the_refusal_says_the_tree_was_left_alone() -> None:
    """
    Someone who sees an error from a tool that rewrites source needs to know, first,
    whether their files were touched. The refusal happens before any write, so it says so.
    """
    message = mutation_guard.format_dirty_tree_refusal(["vntyper/scripts/scoring.py"])

    assert "Nothing has been modified." in message
    assert message.startswith("ERROR:")


def test_indeterminate_refusal_names_the_two_guarded_boundaries() -> None:
    """Unknown Git state must not be confused with real-source mutation."""
    message = mutation_guard.format_indeterminate_refusal("git unavailable")

    assert "cannot verify selected targets and requested outputs" in message
    assert "rewrites production source in place" not in message


def test_dirty_refusal_explains_why_targets_and_outputs_must_be_clean() -> None:
    """The refusal distinguishes snapshot identity from real-output replacement."""
    message = mutation_guard.format_dirty_tree_refusal(["vntyper/scripts/scoring.py"])

    assert "selected targets define the committed measurement baseline" in message
    assert "requested outputs are replaced in the real checkout" in message
    assert "rewrites each target in place" not in message


def test_the_unrestored_warning_points_at_a_read_only_inspection() -> None:
    """`git diff` shows the damage without touching it - that is the whole advice."""
    message = mutation_guard.format_unrestored_warning(["vntyper/scripts/scoring.py"])

    assert "git diff" in message
    assert "do not commit" in message
