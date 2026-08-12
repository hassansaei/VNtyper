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


def test_a_missing_pinned_revision_aborts_rather_than_using_the_branch_tip():
    """Silently falling back to the branch tip is the exact failure the pin prevents.

    The build would report success while having compiled a tree nobody chose -- and the
    evidence comment would then be describing code that is not installed.
    """
    installer = INSTALLER.read_text(encoding="utf-8")
    checkout_block = installer[installer.index("GIT_COMMIT") :]

    assert "exit 1" in checkout_block, "a failed checkout of the pinned revision must abort the install"


def test_the_branch_is_the_forks_default():
    """The fork's default branch used to be `master` while VNtyper installed `enhanced_hmm`,
    so the default branch was not the code that runs. `main` now carries that work.

    The branch is only a clone hint -- GIT_COMMIT determines the tree -- but a default that
    disagrees with what is installed is a trap for anyone reading the fork.
    """
    assert _cfg_value("GIT_BRANCH") == "main"


def test_the_pinned_repository_is_the_fork_the_evidence_came_from():
    assert _cfg_value("GIT_REPO") == "https://github.com/berntpopp/adVNTR.git"
