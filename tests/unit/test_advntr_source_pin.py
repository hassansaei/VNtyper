"""The adVNTR source pin, and why it has to be a commit (#254).

`vntyper/modules/advntr/advntr_genotyping.py` records *why* adVNTR runs at `-t 1` by citing
exact upstream lines: `advntr_commands.py:72-74` (`-t` sets only `settings.CORES`),
`models.py:66,69,72,74` and `vntr_finder.py:870,889` (its only three readers, none on
VNtyper's `genotype -fs` path), and `genome_analyzer.py:211`. That comment exists to stop
#215 being re-litigated in two years.

Line numbers are only true of a specific revision. `install_advntr.cfg` pinned
`GIT_BRANCH="enhanced_hmm"` -- a **mutable branch** -- so a later install could fetch a
different tree and silently invalidate the evidence, while every other upstream artefact here
is pinned immutably: the seed URLs name a commit, the bundles carry a sha256, the Kestrel JAR
is vendored.

These tests pin the shape of the configuration. Whether the cited lines still say what they
are claimed to say is verified against the pinned revision itself, which is recorded in the
commit that introduced the pin.
"""

import re
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

ADVNTR_DIR = Path(__file__).resolve().parents[2] / "vntyper" / "dependencies" / "advntr"
CFG = ADVNTR_DIR / "install_advntr.cfg"
INSTALLER = ADVNTR_DIR / "install_advntr.sh"


def _cfg_value(name: str) -> str | None:
    """Read a shell-assignment value out of the cfg, ignoring comments."""
    for line in CFG.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("#") or "=" not in stripped:
            continue
        key, _, value = stripped.partition("=")
        if key.strip() == name:
            return value.strip().strip('"').strip("'")
    return None


def test_the_source_revision_is_pinned_to_a_full_commit_sha():
    """A branch is mutable; the evidence in advntr_genotyping.py is not.

    A short SHA is refused too: it is ambiguous in principle and unverifiable at a glance
    against the revision recorded alongside the evidence.
    """
    commit = _cfg_value("GIT_COMMIT")

    assert commit, "GIT_COMMIT must be set - a branch alone cannot pin evidence"
    assert re.fullmatch(r"[0-9a-f]{40}", commit), f"GIT_COMMIT must be a full 40-character SHA, got {commit!r}"


def test_the_installer_checks_out_the_pinned_revision():
    """Declaring a pin that nothing acts on is worse than no pin: it reads as a guarantee."""
    installer = INSTALLER.read_text(encoding="utf-8")

    assert "GIT_COMMIT" in installer, "the installer must consume the pin"
    assert "git checkout" in installer, "the installer must actually check the revision out"


def test_the_branch_is_the_forks_default():
    """The fork's default branch used to be `master` while VNtyper installed `enhanced_hmm`,
    so the default branch was not the code that runs. `main` now carries that work.

    The branch is only a clone hint -- GIT_COMMIT determines the tree -- but a default that
    disagrees with what is installed is a trap for anyone reading the fork.
    """
    assert _cfg_value("GIT_BRANCH") == "main"


def test_the_pinned_repository_is_the_fork_the_evidence_came_from():
    assert _cfg_value("GIT_REPO") == "https://github.com/berntpopp/adVNTR.git"


def test_the_config_is_resolved_beside_the_script_not_in_the_working_directory():
    """The image copies this directory to /tmp/advntr/ and runs the script from elsewhere.

    A bare relative `install_advntr.cfg` was therefore never sourced there: the build used
    the script's own fallbacks, so the GIT_COMMIT pin would have had no effect on the very
    image it exists to pin. A config read only when you happen to be standing in the right
    directory is worse than no config, because it looks like it applied.
    """
    installer = INSTALLER.read_text(encoding="utf-8")

    assert "BASH_SOURCE" in installer, "the config must be resolved relative to the script"
    assert 'CONFIG_FILE="install_advntr.cfg"' not in installer, "a bare relative name resolves against the CWD"


# --- Behaviour, by running the installer ------------------------------------------------
#
# Everything above reads the script as text, which cannot tell a live `git checkout` from one
# sitting in a comment. These run it. The installer prints its settings and then refuses an
# existing --install-dir without --overwrite, so it reports what it resolved and exits before
# it would clone: no network, no git, deterministic.

requires_bash = pytest.mark.skipif(shutil.which("bash") is None, reason="bash is not available")


def _resolved_settings(tmp_path: Path, *args: str) -> str:
    """Run the installer far enough to print its settings, and no further."""
    stop_here = tmp_path / "already-installed"
    stop_here.mkdir()

    result = subprocess.run(
        ["bash", str(INSTALLER), "-d", str(stop_here), *args],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=60,
        check=False,
    )

    assert result.returncode == 1, f"expected the existing-directory refusal, got {result.returncode}: {result.stderr}"
    assert "already exists" in result.stdout, "the run must stop before cloning"
    return result.stdout


@requires_bash
def test_the_shipped_pin_reaches_a_run_started_from_another_directory(tmp_path):
    """The whole point of resolving beside the script, checked by running it from elsewhere.

    ``docker/Dockerfile.base`` runs ``bash /tmp/advntr/install_advntr.sh`` with the build's
    working directory somewhere else entirely. Reading the script's text cannot show that the
    shipped config is actually picked up there; running it from an unrelated cwd can.
    """
    settings = _resolved_settings(tmp_path)

    assert _cfg_value("GIT_COMMIT") in settings, "the shipped pin must reach a run started elsewhere"
    assert _cfg_value("GIT_BRANCH") in settings
    assert _cfg_value("GIT_REPO") in settings


@requires_bash
def test_an_explicit_config_replaces_the_shipped_one_rather_than_layering_over_it(tmp_path):
    """``-c`` has to win completely, not partially.

    Sourcing the shipped config first and letting ``-c`` override afterwards leaks every
    value the caller's file does not set. GIT_COMMIT is the one that matters: a ``-c`` that
    selects a different branch would clone that branch and then check out *this* repository's
    pinned commit on top of it -- silently building a tree the caller did not choose -- or
    abort claiming their pinned revision is missing when they pinned nothing.

    This is only a hazard because the default now resolves beside the script and therefore
    always exists. Before that it was a bare relative name that usually matched nothing.
    """
    override = tmp_path / "local.cfg"
    override.write_text('GIT_REPO="https://example.invalid/fork.git"\nGIT_BRANCH="my-branch"\n', encoding="utf-8")

    settings = _resolved_settings(tmp_path, "-c", str(override))

    assert "my-branch" in settings, "the override's branch must be used"
    assert "example.invalid" in settings, "the override's repository must be used"
    assert _cfg_value("GIT_COMMIT") not in settings, (
        "the shipped commit pin leaked into a run configured by --config; "
        "it would be checked out on top of the caller's branch"
    )


@requires_bash
def test_a_config_path_that_does_not_exist_is_refused(tmp_path):
    """Silently falling back to the shipped config would build something other than asked."""
    result = subprocess.run(
        ["bash", str(INSTALLER), "-c", str(tmp_path / "nope.cfg")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=60,
        check=False,
    )

    assert result.returncode == 1
    assert "not found" in result.stdout + result.stderr


@requires_bash
def test_a_symlinked_invocation_still_finds_the_shipped_pin(tmp_path):
    """``dirname "${BASH_SOURCE[0]}"`` is the directory of the *name*, not of the file.

    Linking the script onto PATH -- ``/usr/local/bin/install_advntr.sh`` pointing here -- is
    an ordinary thing to do. Without resolving the link, the config is looked for beside the
    symlink, is not found, and every value falls through to the built-in defaults. GIT_COMMIT
    has no built-in default, so the install silently becomes unpinned and builds the branch
    tip: the exact failure the pin exists to prevent, reached by a route that looks like
    normal use.
    """
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    link = bin_dir / "install_advntr.sh"
    link.symlink_to(INSTALLER)
    stop_here = tmp_path / "already-installed"
    stop_here.mkdir()

    result = subprocess.run(
        ["bash", str(link), "-d", str(stop_here)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=60,
        check=False,
    )

    assert _cfg_value("GIT_COMMIT") in result.stdout, "the pin must survive a symlinked invocation"
    assert "no GIT_COMMIT pin is set" not in result.stderr


@requires_bash
def test_an_unpinned_install_warns_rather_than_proceeding_quietly(tmp_path):
    """An unpinned install stays possible -- the script must work for someone who copied it
    out on its own -- but never quietly. A build that takes the branch tip without saying so
    is indistinguishable from one that honoured the pin.
    """
    lone = tmp_path / "install_advntr.sh"
    lone.write_text(INSTALLER.read_text(encoding="utf-8"), encoding="utf-8")
    stop_here = tmp_path / "already-installed"
    stop_here.mkdir()

    result = subprocess.run(
        ["bash", str(lone), "-d", str(stop_here)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=60,
        check=False,
    )

    assert "no GIT_COMMIT pin is set" in result.stderr
    assert "No configuration file was found" in result.stderr


@requires_bash
@pytest.mark.skipif(shutil.which("git") is None, reason="git is not available")
def test_an_absent_pinned_revision_actually_aborts_the_install(tmp_path):
    """The abort path, executed rather than grepped for.

    This replaced a test that sliced the script from the first textual ``GIT_COMMIT`` to EOF
    and asserted ``exit 1`` appeared somewhere in it -- which any unrelated ``exit 1`` below
    that point satisfied. Here a real (local, offline) repository is cloned and a revision
    that is not in it is pinned, so the only way to reach a zero exit is for the
    fallback-to-branch-tip failure to be real.
    """
    origin = tmp_path / "origin"
    origin.mkdir()
    env = {
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@example.invalid",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@example.invalid",
        "PATH": "/usr/bin:/bin",
        "HOME": str(tmp_path),
    }
    (origin / "README").write_text("x\n", encoding="utf-8")
    for command in (["init", "-b", "main"], ["add", "-A"], ["commit", "-m", "c"]):
        subprocess.run(["git", *command], cwd=origin, env=env, check=True, capture_output=True)

    absent = "0" * 40
    override = tmp_path / "local.cfg"
    override.write_text(
        f'GIT_REPO="{origin}"\nGIT_BRANCH="main"\nGIT_COMMIT="{absent}"\n',
        encoding="utf-8",
    )

    result = subprocess.run(
        ["bash", str(INSTALLER), "-c", str(override), "-d", str(tmp_path / "build")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env=env,
        timeout=120,
        check=False,
    )

    assert result.returncode != 0, "an absent pinned revision must not build the branch tip"
    assert absent in result.stderr, "the abort must name the revision it could not find"
    assert "Refusing to build against the branch tip" in result.stderr
