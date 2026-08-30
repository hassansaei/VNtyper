"""CLI and Git-boundary tests for integration compatibility checking."""

import copy
import importlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any

import pytest

from tests.unit.test_integration_compatibility import _contract, _live_case, _manifest, _resources

pytestmark = pytest.mark.unit


def _module() -> Any:
    try:
        return importlib.import_module("scripts.check_integration_compatibility")
    except ModuleNotFoundError:
        return None


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True, text=True)
    return result.stdout.strip()


def _repository(tmp_path: Path) -> tuple[Path, str, str]:
    repo = tmp_path / "repo"
    repo.mkdir(parents=True)
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "tests@example.invalid")
    _git(repo, "config", "user.name", "Tests")
    contract = _contract()
    paths = {
        "tests/compatibility/real_success_baseline.json": _manifest(contract),
        "tests/test_data_config.json": {
            **_resources("tests/data/input.bam"),
            "integration_tests": {"bam_tests": [_live_case(contract)]},
        },
    }
    for relative, value in paths.items():
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
    version_path = repo / "vntyper/version.py"
    version_path.parent.mkdir(parents=True)
    version_path.write_text('__version__ = "2.0.24"\n')
    _git(repo, "add", ".")
    _git(repo, "commit", "-qm", "base")
    base = _git(repo, "rev-parse", "HEAD")
    second = _contract(name="case-b", path="tests/data/second.bam")
    manifest_path = repo / "tests/compatibility/real_success_baseline.json"
    manifest_path.write_text(json.dumps(_manifest(contract, second)))
    config_path = repo / "tests/test_data_config.json"
    config = json.loads(config_path.read_text())
    config["file_resources"].append({"local_path": "tests/data", "filename": "second.bam", "md5sum": "a" * 32})
    config["integration_tests"]["bam_tests"].append(_live_case(second))
    config_path.write_text(json.dumps(config))
    _git(repo, "add", ".")
    _git(repo, "commit", "-qm", "append")
    head = _git(repo, "rev-parse", "HEAD")
    return repo, base, head


def test_main_accepts_explicit_pr_or_direct_push_base(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Catch comparing against HEAD/origin-main instead of the supplied event base."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, _ = _repository(tmp_path)

    exit_code = module.main(["--base-revision", base, "--repo-root", str(repo)], environ={"CI": "true"})

    assert exit_code == 0
    assert "passed" in capsys.readouterr().out.lower()


def test_main_passes_authoritative_checkout_version_to_compatibility_engine(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Catch selecting observations from mutable test JSON instead of package metadata."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, _ = _repository(tmp_path)
    captured: dict[str, object] = {}

    def capture(*_args: object, **kwargs: object) -> None:
        captured.update(kwargs)

    monkeypatch.setattr(module, "check_compatibility", capture)

    assert module.main(["--base-revision", base, "--repo-root", str(repo)], environ={"CI": "true"}) == 0
    assert captured["observation_version"] == "2.0.24"


def test_checkout_version_reader_fails_closed_for_missing_or_ambiguous_assignment(tmp_path: Path) -> None:
    """Never fall back to whichever VNtyper package is installed in the environment."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    with pytest.raises(ValueError, match="cannot read current package version"):
        module._read_checkout_version(tmp_path)
    version_path = tmp_path / "vntyper/version.py"
    version_path.parent.mkdir()
    version_path.write_text('__version__ = "2.0.24"\n__version__ = "2.0.25"\n')
    with pytest.raises(ValueError, match="missing or ambiguous"):
        module._read_checkout_version(tmp_path)


def test_main_requires_explicit_base_in_ci(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Catch silently using a local merge-base fallback as CI evidence."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, _, _ = _repository(tmp_path)

    assert module.main(["--repo-root", str(repo)], environ={"CI": "true"}) == 1
    assert "--base-revision is required in CI" in capsys.readouterr().err


def test_local_fallback_prints_local_only_evidence(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Catch presenting an inferred local merge base as authoritative CI evidence."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, _ = _repository(tmp_path)
    _git(repo, "remote", "add", "origin", str(repo))
    _git(repo, "update-ref", "refs/remotes/origin/main", base)

    assert module.main(["--repo-root", str(repo)], environ={}) == 0
    assert "local-only evidence" in capsys.readouterr().err


@pytest.mark.parametrize("revision", ["", "0" * 40, "missing-revision"])
def test_main_fails_closed_for_invalid_base(tmp_path: Path, capsys: pytest.CaptureFixture[str], revision: str) -> None:
    """Catch treating missing, unreachable, or all-zero event bases as an empty bootstrap."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, _, _ = _repository(tmp_path)

    assert module.main(["--base-revision", revision, "--repo-root", str(repo)], environ={"CI": "true"}) == 1
    assert "invalid base revision" in capsys.readouterr().err


def test_read_json_at_revision_fails_on_malformed_json(tmp_path: Path) -> None:
    """Catch treating malformed Git content as an absent manifest."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, _, _ = _repository(tmp_path)
    path = repo / "broken.json"
    path.write_text("{")
    _git(repo, "add", "broken.json")
    _git(repo, "commit", "-qm", "broken")
    revision = _git(repo, "rev-parse", "HEAD")

    with pytest.raises(ValueError, match="malformed JSON"):
        module.read_json_at_revision(repo, revision, "broken.json")


def test_read_json_at_revision_distinguishes_absent_path(tmp_path: Path) -> None:
    """Catch conflating a legitimate pre-manifest revision with Git command failure."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, _, head = _repository(tmp_path)

    assert module.read_json_at_revision(repo, head, "absent.json", allow_absent=True) is None
    with pytest.raises(ValueError, match="path is absent"):
        module.read_json_at_revision(repo, head, "absent.json")


def test_main_rejects_absent_custom_manifest_when_default_exists_at_base(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Catch treating a missing custom manifest as permission to re-enter bootstrap."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, _ = _repository(tmp_path)
    custom_manifest = "tests/compatibility/custom.json"
    custom_path = repo / custom_manifest
    custom_path.write_text((repo / module.DEFAULT_MANIFEST).read_text())

    exit_code = module.main(
        ["--base-revision", base, "--repo-root", str(repo), "--manifest", custom_manifest],
        environ={"CI": "true"},
    )

    assert exit_code == 1
    error = capsys.readouterr().err
    assert custom_manifest in error
    assert "default manifest already exists at base" in error


def test_resolve_base_rejects_shallow_repository(tmp_path: Path) -> None:
    """Catch accepting incomplete history where append-only evidence cannot be proven."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    source, base, _ = _repository(tmp_path / "source")
    shallow = tmp_path / "shallow"
    subprocess.run(
        ["git", "clone", "-q", "--depth", "1", f"file://{source}", str(shallow)],
        check=True,
        env={**os.environ, "GIT_ALLOW_PROTOCOL": "file"},
    )

    with pytest.raises(ValueError, match="shallow"):
        module.resolve_base_revision(shallow, base, ci=True)


def test_main_fails_closed_when_branch_is_behind_supplied_base(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Catch replacing an ahead PR base with a merge-base that hides missing contracts."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, head = _repository(tmp_path)
    _git(repo, "checkout", "-q", base)

    assert module.main(["--base-revision", head, "--repo-root", str(repo)], environ={"CI": "true"}) == 1
    assert "base revision is not an ancestor of HEAD" in capsys.readouterr().err


def test_resolve_base_uses_merge_base_when_supplied_base_tip_advanced(tmp_path: Path) -> None:
    """Catch rejecting a feature branch only because its explicit base branch advanced."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, feature_head = _repository(tmp_path)
    _git(repo, "branch", "feature", feature_head)
    _git(repo, "checkout", "-q", "-b", "advanced-base", base)
    _git(repo, "commit", "--allow-empty", "-qm", "advance base")
    advanced_base = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "feature")

    assert module.resolve_base_revision(repo, advanced_base, ci=True, allow_merge_base=True) == base


def test_main_rejects_non_fast_forward_push_that_removes_before_contract(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Catch hiding a removed push-before contract by falling back to an older common ancestor."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, before = _repository(tmp_path)
    _git(repo, "checkout", "-q", "-b", "rewritten", base)
    _git(repo, "commit", "--allow-empty", "-qm", "rewrite without appended contract")

    exit_code = module.main(
        ["--base-revision", before, "--repo-root", str(repo)],
        environ={"CI": "true", "EVENT_NAME": "push"},
    )

    assert exit_code == 1
    assert "base revision is not an ancestor of HEAD" in capsys.readouterr().err


def test_main_rejects_ambiguous_criss_cross_pull_request_merge_bases(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Catch silently choosing one of multiple best pull-request merge bases."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, base, _ = _repository(tmp_path)
    _git(repo, "checkout", "-q", "-b", "left", base)
    _git(repo, "commit", "--allow-empty", "-qm", "left base")
    left = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "-b", "right", base)
    _git(repo, "commit", "--allow-empty", "-qm", "right base")
    right = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "-b", "pr-head", left)
    _git(repo, "merge", "-q", "--no-ff", "-m", "merge right into PR", right)
    _git(repo, "checkout", "-q", "-b", "advanced-base", right)
    _git(repo, "merge", "-q", "--no-ff", "-m", "merge left into base", left)
    advanced_base = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "pr-head")
    assert set(_git(repo, "merge-base", "--all", advanced_base, "HEAD").splitlines()) == {left, right}

    exit_code = module.main(
        ["--base-revision", advanced_base, "--repo-root", str(repo)],
        environ={"CI": "true", "EVENT_NAME": "pull_request"},
    )

    assert exit_code == 1
    assert "no unique exact merge base" in capsys.readouterr().err


def test_resolve_base_rejects_unrelated_history(tmp_path: Path) -> None:
    """Catch treating a resolved commit with no common history as usable evidence."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo, _, feature_head = _repository(tmp_path)
    _git(repo, "checkout", "-q", "--orphan", "unrelated")
    _git(repo, "commit", "-qm", "unrelated root")
    unrelated = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", feature_head)

    with pytest.raises(ValueError, match="no unique exact merge base"):
        module.resolve_base_revision(repo, unrelated, ci=True, allow_merge_base=True)


@pytest.mark.parametrize(
    ("failure", "message"),
    [
        ("initial-ancestry", "Git ancestry check failed"),
        ("reverse-ancestry", "Git ancestry check failed for HEAD"),
        ("merge-output", "no unique exact merge base"),
        ("merge-ancestry", "Git merge-base ancestry check failed"),
    ],
)
def test_resolve_base_fails_closed_for_inconsistent_git_results(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
    message: str,
) -> None:
    """Catch accepting subprocess errors or a merge base Git cannot prove is an ancestor."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    supplied = "a" * 40
    common = "b" * 40

    def fake_git(_repo_root: Path, arguments: list[str]) -> subprocess.CompletedProcess[str]:
        result = subprocess.CompletedProcess(["git", *arguments], 0, "", "")
        if arguments == ["rev-parse", "--is-shallow-repository"]:
            result.stdout = "false\n"
        elif arguments == ["rev-parse", "--verify", f"{supplied}^{{commit}}"]:
            result.stdout = f"{supplied}\n"
        elif arguments == ["merge-base", "--is-ancestor", supplied, "HEAD"]:
            result.returncode = 2 if failure == "initial-ancestry" else 1
        elif arguments == ["merge-base", "--is-ancestor", "HEAD", supplied]:
            result.returncode = 2 if failure == "reverse-ancestry" else 1
        elif arguments == ["merge-base", "--all", supplied, "HEAD"]:
            result.stdout = "invalid\n" if failure == "merge-output" else f"{common}\n"
        elif arguments == ["merge-base", "--is-ancestor", common, "HEAD"]:
            result.returncode = 1 if failure == "merge-ancestry" else 0
        else:  # pragma: no cover - the assertion exposes an unexpected production command
            raise AssertionError(f"unexpected Git command: {arguments}")
        return result

    monkeypatch.setattr(module, "_git", fake_git)

    with pytest.raises(ValueError, match=message):
        module.resolve_base_revision(tmp_path, supplied, ci=True, allow_merge_base=True)


def test_git_boundaries_fail_for_missing_worktree_and_inconsistent_tree_output(tmp_path: Path) -> None:
    """Catch subprocess setup failures and accepting a directory listing as one JSON file."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    with pytest.raises(ValueError, match="Git subprocess failed"):
        module.resolve_base_revision(tmp_path / "absent", "HEAD", ci=True)
    repo, _, head = _repository(tmp_path)
    with pytest.raises(ValueError, match="inconsistent Git output"):
        module.read_json_at_revision(repo, head, "tests")


def test_observation_provenance_resolves_to_matching_ancestral_version(tmp_path: Path) -> None:
    """Bind each versioned observation to reachable package metadata at its exact commit."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "tests@example.invalid")
    _git(repo, "config", "user.name", "Tests")
    version_path = repo / "vntyper/version.py"
    version_path.parent.mkdir()
    version_path.write_text('__version__ = "2.0.24"\n')
    _git(repo, "add", ".")
    _git(repo, "commit", "-qm", "release")
    provenance = _git(repo, "rev-parse", "HEAD")
    manifest: dict[str, Any] = {
        "schema_version": 2,
        "contracts": [],
        "observation_sets": [
            {
                "version": "2.0.24",
                "provenance_commit": provenance,
                "extends": None,
                "report_overrides": [],
            }
        ],
    }

    module.verify_observation_provenance(repo, manifest)

    wrong_version = copy.deepcopy(manifest)
    wrong_version["observation_sets"][0]["version"] = "2.0.25"
    with pytest.raises(ValueError, match="does not match package version"):
        module.verify_observation_provenance(repo, wrong_version)
    missing = copy.deepcopy(manifest)
    missing["observation_sets"][0]["provenance_commit"] = "a" * 40
    with pytest.raises(ValueError, match="does not resolve"):
        module.verify_observation_provenance(repo, missing)


def test_observation_provenance_rejects_non_ancestral_commit(tmp_path: Path) -> None:
    """A detached observation commit cannot attest to the current branch."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "tests@example.invalid")
    _git(repo, "config", "user.name", "Tests")
    version_path = repo / "vntyper/version.py"
    version_path.parent.mkdir()
    version_path.write_text('__version__ = "2.0.24"\n')
    _git(repo, "add", ".")
    _git(repo, "commit", "-qm", "main")
    main = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "--orphan", "detached")
    _git(repo, "commit", "--allow-empty", "-qm", "detached")
    detached = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", "main")
    assert _git(repo, "rev-parse", "HEAD") == main
    manifest = {
        "schema_version": 2,
        "contracts": [],
        "observation_sets": [
            {
                "version": "2.0.24",
                "provenance_commit": detached,
                "extends": None,
                "report_overrides": [],
            }
        ],
    }

    with pytest.raises(ValueError, match="not an ancestor"):
        module.verify_observation_provenance(repo, manifest)


@pytest.mark.parametrize(
    ("observations", "message"),
    [
        ({}, "observation_sets must be a list"),
        (["not-an-object"], r"observation_sets\[0\] must be an object"),
    ],
)
def test_observation_provenance_rejects_malformed_containers(
    tmp_path: Path, observations: object, message: str
) -> None:
    """The Git boundary also fails closed if called before core schema validation."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    with pytest.raises(ValueError, match=message):
        module.verify_observation_provenance(
            tmp_path,
            {"schema_version": 2, "contracts": [], "observation_sets": observations},
        )


@pytest.mark.parametrize(
    ("failure", "message"),
    [
        ("ancestry-error", "Git ancestry check failed"),
        ("version-unreadable", "cannot read package version"),
        ("version-ambiguous", "package version .* is ambiguous"),
    ],
)
def test_observation_provenance_fails_closed_for_git_and_version_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
    message: str,
) -> None:
    """Subprocess errors and ambiguous metadata cannot become valid provenance."""
    module = _module()
    assert module is not None, "integration compatibility CLI is not implemented"
    provenance = "a" * 40
    manifest = {
        "schema_version": 2,
        "contracts": [],
        "observation_sets": [
            {
                "version": "2.0.24",
                "provenance_commit": provenance,
                "extends": None,
                "report_overrides": [],
            }
        ],
    }

    def fake_git(_repo_root: Path, arguments: list[str]) -> subprocess.CompletedProcess[str]:
        if arguments[0] == "rev-parse":
            return subprocess.CompletedProcess(arguments, 0, f"{provenance}\n", "")
        if arguments[0] == "merge-base":
            return subprocess.CompletedProcess(arguments, 2 if failure == "ancestry-error" else 0, "", "")
        if failure == "version-unreadable":
            return subprocess.CompletedProcess(arguments, 1, "", "missing")
        output = '__version__ = "2.0.24"\n'
        if failure == "version-ambiguous":
            output += '__version__ = "2.0.24"\n'
        return subprocess.CompletedProcess(arguments, 0, output, "")

    monkeypatch.setattr(module, "_git", fake_git)

    with pytest.raises(ValueError, match=message):
        module.verify_observation_provenance(tmp_path, manifest)
