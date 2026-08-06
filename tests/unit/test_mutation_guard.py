"""
Unit tests for ``scripts/mutation_guard.py``, the mutation harness's working-tree guard.

The harness rewrites production source in place and restores whatever it read at the
start, so it is only correct on a clean tree. This module is what enforces that, and the
property worth pinning hardest is negative: **no message it produces may name a command
that discards uncommitted work**. The one moment a user reads these messages is the
moment their own edits are on the line, and the harness's previous advice - ``git
checkout --`` - would have deleted exactly the work it was misdiagnosing.
"""

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


def test_a_failing_git_call_does_not_block_the_sweep(tmp_path, monkeypatch) -> None:
    """Outside a repository the guard has nothing to protect, so it stays out of the way."""

    def fake_run(command, **_kwargs):
        return subprocess.CompletedProcess(command, 128, stdout="", stderr="not a git repository")

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert mutation_guard.dirty_paths(tmp_path, ["a.py"]) == []


def test_a_missing_git_binary_does_not_block_the_sweep(tmp_path, monkeypatch) -> None:
    """The guard is a safety net, not a dependency: no git, no refusal."""

    def fake_run(*_args, **_kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert mutation_guard.dirty_paths(tmp_path, ["a.py"]) == []


# ---------------------------------------------------------------------------
# Wording the two failures
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "formatter",
    [mutation_guard.format_dirty_tree_refusal, mutation_guard.format_unrestored_warning],
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


def test_the_unrestored_warning_points_at_a_read_only_inspection() -> None:
    """`git diff` shows the damage without touching it - that is the whole advice."""
    message = mutation_guard.format_unrestored_warning(["vntyper/scripts/scoring.py"])

    assert "git diff" in message
    assert "do not commit" in message
