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


#: A shell assignment: an identifier, ``=``, then the value. Anchored and identifier-shaped
#: so that prose containing ``=`` cannot be mistaken for one -- the cfg's comments cite
#: upstream lines like ``settings.CORES = args.threads``, and a looser split would read that
#: as an assignment to ``settings.CORES`` the moment a comment marker moved.
_ASSIGNMENT = re.compile(r"^(?P<key>[A-Za-z_][A-Za-z0-9_]*)=(?P<value>.*)$")


def _cfg_value(name: str) -> str | None:
    """Read a shell-assignment value out of the cfg, ignoring comments.

    Quotes are stripped as one matching surrounding pair rather than character by
    character, so a value that legitimately contains a quote survives intact.
    """
    for line in CFG.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        match = _ASSIGNMENT.match(stripped)
        if match is None or match.group("key") != name:
            continue
        value = match.group("value").strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            value = value[1:-1]
        return value
    return None


def test_the_cfg_parser_reads_assignments_and_not_prose():
    """The helper every assertion below rests on, checked directly.

    If it silently returned None the config tests would still pass their `== "main"` style
    comparisons only by luck, and `GIT_COMMIT in stdout` assertions would compare against
    nothing. Its failure mode is quiet, so it gets its own test.
    """
    assert _cfg_value("GIT_BRANCH") == "main"
    assert _cfg_value("NOT_A_REAL_KEY") is None
    # The cfg's comments cite `settings.CORES = args.threads`; that must not parse.
    assert _cfg_value("settings.CORES") is None
    assert "#" not in (_cfg_value("GIT_REPO") or ""), "a comment must not leak into a value"


def test_the_source_revision_is_pinned_to_a_full_commit_sha():
    """A branch is mutable; the evidence in advntr_genotyping.py is not.

    A short SHA is refused too: it is ambiguous in principle and unverifiable at a glance
    against the revision recorded alongside the evidence.
    """
    commit = _cfg_value("GIT_COMMIT")

    assert commit, "GIT_COMMIT must be set - a branch alone cannot pin evidence"
    assert re.fullmatch(r"[0-9a-f]{40}", commit), f"GIT_COMMIT must be a full 40-character SHA, got {commit!r}"


def test_the_branch_is_the_forks_default():
    """The fork's default branch used to be `master` while VNtyper installed `enhanced_hmm`,
    so the default branch was not the code that runs. `main` now carries that work.

    The branch is only a clone hint -- GIT_COMMIT determines the tree -- but a default that
    disagrees with what is installed is a trap for anyone reading the fork.
    """
    assert _cfg_value("GIT_BRANCH") == "main"


def test_the_pinned_repository_is_the_fork_the_evidence_came_from():
    assert _cfg_value("GIT_REPO") == "https://github.com/berntpopp/adVNTR.git"


# Everything above reads the cfg, which is data. Nothing here greps the *installer* text:
# reading a shell script as a string cannot tell a live `git checkout` from one sitting in a
# comment, and every such assertion this file used to carry has a behavioural equivalent
# below. Two are gone entirely -- one asserted `exit 1` appeared somewhere below the first
# mention of GIT_COMMIT, which any unrelated `exit 1` satisfied.


# --- Behaviour, by running the installer ------------------------------------------------
#
# Everything above reads the script as text, which cannot tell a live `git checkout` from one
# sitting in a comment. These run it. The installer prints its settings and then refuses an
# existing --install-dir without --overwrite, so it reports what it resolved and exits before
# it would clone: no network, no git, deterministic.

requires_bash = pytest.mark.skipif(shutil.which("bash") is None, reason="bash is not available")

#: A minimal environment for the tests that drive real git. Identity is supplied so the run
#: does not depend on the machine's git config, and PATH is pinned so it does not depend on
#: what else happens to be installed. Callers add HOME, which must be the test's tmp_path so
#: git cannot read the developer's ~/.gitconfig.
_GIT_ENV = {
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.invalid",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.invalid",
    "PATH": "/usr/bin:/bin",
}


def _resolved_settings(tmp_path: Path, *args: str) -> str:
    """Run the installer far enough to print its settings, and no further.

    ``-e ""`` clears CONDA_ENV. The shipped config sets it to ``envadvntr``, and the script
    activates that environment -- exiting non-zero if conda is absent *or* the environment
    does not exist -- before it reaches the existing-directory refusal this relies on to
    stop. Without clearing it the test asserts something different depending on which conda
    environments the machine happens to have, which is how it first went green locally and
    red on every CI runner.
    """
    stop_here = tmp_path / "already-installed"
    stop_here.mkdir()

    result = subprocess.run(
        ["bash", str(INSTALLER), "-e", "", "-d", str(stop_here), *args],
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
        ["bash", str(link), "-e", "", "-d", str(stop_here)],
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
        ["bash", str(lone), "-e", "", "-d", str(stop_here)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=60,
        check=False,
    )

    assert "no GIT_COMMIT pin is set" in result.stderr
    assert "No configuration file was found" in result.stderr


def _local_repo_with_two_commits(tmp_path: Path) -> tuple[Path, str, str]:
    """Build an offline git repository and return it with its two revisions, oldest first."""
    origin = tmp_path / "origin"
    origin.mkdir()
    env = _GIT_ENV | {"HOME": str(tmp_path)}
    revisions = []
    for index, command in enumerate((["init", "-b", "main"], ["add", "-A"], ["commit", "-m", "one"])):
        if index == 0:
            (origin / "README").write_text("one\n", encoding="utf-8")
        subprocess.run(["git", *command], cwd=origin, env=env, check=True, capture_output=True)
    revisions.append(
        subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=origin, env=env, check=True, capture_output=True, text=True
        ).stdout.strip()
    )
    (origin / "README").write_text("two\n", encoding="utf-8")
    for command in (["add", "-A"], ["commit", "-m", "two"]):
        subprocess.run(["git", *command], cwd=origin, env=env, check=True, capture_output=True)
    revisions.append(
        subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=origin, env=env, check=True, capture_output=True, text=True
        ).stdout.strip()
    )
    return origin, revisions[0], revisions[1]


@requires_bash
@pytest.mark.skipif(shutil.which("git") is None, reason="git is not available")
def test_the_pinned_revision_is_what_gets_built_not_the_branch_tip(tmp_path):
    """The positive half: a pin that *is* present must actually be checked out.

    ``test_an_absent_pinned_revision_actually_aborts_the_install`` proves the failure path,
    which a checkout that silently did nothing would also satisfy -- there would be nothing
    to fail. This pins the parent commit of a two-commit branch and reads back the revision
    the installer reports, so the branch tip and the pin are different values and only one
    of them is right. The run then dies at ``python setup.py install`` against a repository
    that has none; that is after the point being tested, so its exit status is ignored.
    """
    origin, first, second = _local_repo_with_two_commits(tmp_path)
    assert first != second
    override = tmp_path / "local.cfg"
    override.write_text(f'GIT_REPO="{origin}"\nGIT_BRANCH="main"\nGIT_COMMIT="{first}"\n', encoding="utf-8")

    result = subprocess.run(
        ["bash", str(INSTALLER), "-e", "", "-c", str(override), "-d", str(tmp_path / "build")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env=_GIT_ENV | {"HOME": str(tmp_path)},
        timeout=120,
        check=False,
    )

    assert f"adVNTR revision: {first}" in result.stdout, (
        f"the pinned revision was not what got built; stdout was:\n{result.stdout}"
    )
    assert second not in result.stdout, "the branch tip was built instead of the pin"


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
    env = _GIT_ENV | {"HOME": str(tmp_path)}
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


#: The module whose comment the cfg's own comment cross-references.
GENOTYPING = Path(__file__).resolve().parents[2] / "vntyper" / "modules" / "advntr" / "advntr_genotyping.py"

#: A leading comment marker, with the single space that conventionally follows it.
_COMMENT_MARKER = re.compile(r"(?m)^[ \t]*#[ \t]?")

#: A sentence terminator followed by whitespace. ``203.3 s`` and ``subset.bam,`` do not
#: match, because the character after the dot is not whitespace in either.
_SENTENCE_SPLIT = re.compile(r"(?<=[.;])\s")


def _claims_about(text: str, needle: str) -> list[str]:
    """Return the sentences of ``text`` that mention ``needle``.

    Comment markers are stripped and wrapped lines rejoined, because the claim this checks
    is written across several commented lines in both documents. A blank line ends a claim,
    so a sentence can never absorb the paragraph after it -- without that, a ``.py`` file
    would flatten into one string and a version number from unrelated code could satisfy
    the assertion.

    Args:
        text: The whole file.
        needle: The substring a sentence must contain to be returned.

    Returns:
        list[str]: The matching sentences, whitespace-normalised.
    """
    without_markers = _COMMENT_MARKER.sub("", text)
    claims: list[str] = []
    for block in re.split(r"\n\s*\n", without_markers):
        joined = " ".join(block.split())
        claims.extend(sentence for sentence in _SENTENCE_SPLIT.split(joined) if needle in sentence)
    return claims


class TestTheTwoDocumentsAgreeOnWhichReleaseMadeThreadingReal:
    """`install_advntr.cfg` and `advntr_genotyping.py` both state which adVNTR release moved
    the Viterbi DP into a ``nogil`` block, and the cfg explicitly tells a reader to keep the
    two in step. They disagreed: the cfg credited 2.0.1, the module credited 2.0.0.

    2.0.0 is right. `f76d326` ("perf(hmm): bit-exact nogil Viterbi rewrite") and `82b1c2b`
    ("perf: parallelise read decoding, making -t genuinely work") both precede `d66b839`
    ("release: 2.0.0"). Tag `v2.0.1` is `59479ae`, whose one substantive commit is `2e7a3d0`
    ("build: stop shipping the test harness and dev scripts in the egg") -- packaging, not
    decoding.

    Nothing else pins this. The pin tests above check the cfg's *shape*, and a wrong version
    in a comment survives every one of them -- which is exactly how it shipped.
    """

    @pytest.mark.parametrize("path", [CFG, GENOTYPING], ids=["cfg", "module"])
    def test_the_nogil_claim_names_the_release_that_carries_it(self, path):
        claims = _claims_about(path.read_text(encoding="utf-8"), "nogil")

        assert claims, f"{path.name} no longer states which release moved the DP into a nogil block"
        assert any("2.0.0" in claim for claim in claims), (
            f"{path.name} must credit 2.0.0 with the nogil rewrite; found: {claims}"
        )

    @pytest.mark.parametrize("path", [CFG, GENOTYPING], ids=["cfg", "module"])
    def test_no_document_credits_the_packaging_release_with_the_rewrite(self, path):
        claims = _claims_about(path.read_text(encoding="utf-8"), "nogil")

        offenders = [claim for claim in claims if "2.0.1" in claim]
        assert not offenders, (
            f"{path.name} credits 2.0.1 with the nogil rewrite, which is a packaging release: {offenders}"
        )

    def test_the_checker_reads_a_claim_split_across_commented_lines(self):
        """Non-vacuity: the wording that actually shipped must be what this rejects."""
        shipped = (
            "# adVNTR 2.0.1 (tag v2.0.1) moved the Viterbi DP into a `nogil` block and threaded the\n"
            "# read loop. Verified at this revision: 61.1 -> 3.12 ms per decode attempt serially.\n"
        )

        claims = _claims_about(shipped, "nogil")

        assert claims == [
            "adVNTR 2.0.1 (tag v2.0.1) moved the Viterbi DP into a `nogil` block and threaded the read loop."
        ]
