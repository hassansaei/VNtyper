"""Command-line boundary for the inactive integration compatibility engine."""

from __future__ import annotations

import argparse
import importlib
import json
import os
import subprocess
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

_core_name = "scripts.integration_compatibility" if __package__ else "integration_compatibility"
_core = importlib.import_module(_core_name)
BOOTSTRAP_REVISION: str = _core.BOOTSTRAP_REVISION
check_compatibility: Any = _core.check_compatibility

DEFAULT_MANIFEST = "tests/compatibility/real_success_baseline.json"
DEFAULT_TEST_CONFIG = "tests/test_data_config.json"


def _git(repo_root: Path, arguments: Sequence[str]) -> subprocess.CompletedProcess[str]:
    """Run one read-only Git command and capture its complete result."""
    try:
        return subprocess.run(
            ["git", *arguments],
            cwd=repo_root,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as error:
        raise ValueError(f"Git subprocess failed: {error}") from error


def resolve_base_revision(repo_root: Path, revision: str | None, *, ci: bool) -> str:
    """Resolve an explicit event base or a clearly labelled local fallback.

    Args:
        repo_root: Git worktree root.
        revision: Explicit PR base or direct-push ``before`` revision.
        ci: Whether the caller is running in CI.

    Returns:
        Full base commit SHA.

    Raises:
        ValueError: If history or the requested base cannot prove the comparison.
    """
    shallow = _git(repo_root, ["rev-parse", "--is-shallow-repository"])
    if shallow.returncode != 0:
        raise ValueError(f"Git shallow-history check failed: {shallow.stderr.strip()}")
    if shallow.stdout.strip() != "false":
        raise ValueError("shallow Git history cannot prove append-only compatibility")
    if revision is None:
        if ci:
            raise ValueError("--base-revision is required in CI")
        fallback = _git(repo_root, ["merge-base", "origin/main", "HEAD"])
        if fallback.returncode != 0 or not fallback.stdout.strip():
            raise ValueError(f"local merge-base lookup failed: {fallback.stderr.strip()}")
        revision = fallback.stdout.strip()
        print(
            f"Using local-only evidence from merge-base origin/main HEAD: {revision}",
            file=sys.stderr,
        )
    if revision == "0" * 40:
        raise ValueError("invalid base revision: all-zero SHA")
    resolved = _git(repo_root, ["rev-parse", "--verify", f"{revision}^{{commit}}"])
    sha = resolved.stdout.strip()
    if resolved.returncode != 0 or len(sha) != 40:
        raise ValueError(f"invalid base revision {revision!r}: {resolved.stderr.strip()}")
    ancestry = _git(repo_root, ["merge-base", "--is-ancestor", sha, "HEAD"])
    if ancestry.returncode == 1:
        behind = _git(repo_root, ["merge-base", "--is-ancestor", "HEAD", sha])
        if behind.returncode == 0:
            raise ValueError(f"base revision is not an ancestor of HEAD: {sha}")
        if behind.returncode != 1:
            raise ValueError(f"Git ancestry check failed for HEAD and {sha}: {behind.stderr.strip()}")
        merged = _git(repo_root, ["merge-base", sha, "HEAD"])
        merge_bases = [line for line in merged.stdout.splitlines() if line]
        if merged.returncode != 0 or len(merge_bases) != 1 or len(merge_bases[0]) != 40:
            raise ValueError(f"base revision has no exact merge base with HEAD: {sha}: {merged.stderr.strip()}")
        sha = merge_bases[0]
        merge_ancestry = _git(repo_root, ["merge-base", "--is-ancestor", sha, "HEAD"])
        if merge_ancestry.returncode != 0:
            raise ValueError(f"Git merge-base ancestry check failed for {sha}: {merge_ancestry.stderr.strip()}")
    elif ancestry.returncode != 0:
        raise ValueError(f"Git ancestry check failed for {sha}: {ancestry.stderr.strip()}")
    return sha


def read_json_at_revision(
    repo_root: Path,
    revision: str,
    relative_path: str,
    *,
    allow_absent: bool = False,
) -> dict[str, Any] | None:
    """Load a JSON object from one exact Git revision.

    Args:
        repo_root: Git worktree root.
        revision: Verified commit revision.
        relative_path: Repository-relative path.
        allow_absent: Return ``None`` only when the path genuinely is not in the tree.

    Returns:
        Parsed JSON object, or ``None`` for an allowed absent path.

    Raises:
        ValueError: If Git fails, output is inconsistent, or JSON is malformed.
    """
    listed = _git(repo_root, ["ls-tree", "-r", "--name-only", revision, "--", relative_path])
    if listed.returncode != 0:
        raise ValueError(f"Git tree lookup failed for {revision}:{relative_path}: {listed.stderr.strip()}")
    entries = [line for line in listed.stdout.splitlines() if line]
    if not entries:
        if allow_absent:
            return None
        raise ValueError(f"path is absent at {revision}: {relative_path}")
    if entries != [relative_path]:
        raise ValueError(f"inconsistent Git output for {revision}:{relative_path}: {entries}")
    shown = _git(repo_root, ["show", f"{revision}:{relative_path}"])
    if shown.returncode != 0:
        raise ValueError(f"Git show failed for {revision}:{relative_path}: {shown.stderr.strip()}")
    try:
        value = json.loads(shown.stdout)
    except json.JSONDecodeError as error:
        raise ValueError(f"malformed JSON at {revision}:{relative_path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"JSON at {revision}:{relative_path} must be an object")
    return value


def _read_current_json(repo_root: Path, relative_path: str) -> dict[str, Any]:
    path = repo_root / relative_path
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot load current JSON {relative_path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"current JSON {relative_path} must be an object")
    return value


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-revision")
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", default=DEFAULT_MANIFEST)
    parser.add_argument("--test-config", default=DEFAULT_TEST_CONFIG)
    parser.add_argument("--resource-config", default=DEFAULT_TEST_CONFIG)
    return parser


def main(argv: Sequence[str] | None = None, *, environ: Mapping[str, str] | None = None) -> int:
    """Run append-only and bidirectional compatibility validation.

    Args:
        argv: Optional command arguments excluding the executable name.
        environ: Environment used only to determine whether CI requires an explicit base.

    Returns:
        Zero on success and one on any fail-closed validation or Git error.
    """
    args = _parser().parse_args(argv)
    environment = os.environ if environ is None else environ
    repo_root = args.repo_root.resolve()
    try:
        base_revision = resolve_base_revision(repo_root, args.base_revision, ci=bool(environment.get("CI")))
        current_manifest = _read_current_json(repo_root, args.manifest)
        live_config = _read_current_json(repo_root, args.test_config)
        resource_config = _read_current_json(repo_root, args.resource_config)
        base_manifest = read_json_at_revision(repo_root, base_revision, args.manifest, allow_absent=True)
        historical = None
        if base_manifest is None:
            if args.manifest != DEFAULT_MANIFEST:
                default_base_manifest = read_json_at_revision(
                    repo_root,
                    base_revision,
                    DEFAULT_MANIFEST,
                    allow_absent=True,
                )
                if default_base_manifest is not None:
                    raise ValueError(
                        f"custom manifest {args.manifest!r} is absent at {base_revision}; "
                        "default manifest already exists at base, so bootstrap is forbidden"
                    )
            historical = read_json_at_revision(repo_root, BOOTSTRAP_REVISION, args.test_config)
        check_compatibility(
            base_manifest,
            current_manifest,
            live_config,
            resource_config,
            historical_test_config=historical,
        )
    except ValueError as error:
        print(f"Integration compatibility check failed: {error}", file=sys.stderr)
        return 1
    print(f"Integration compatibility check passed against base {base_revision}")
    return 0


if __name__ == "__main__":  # pragma: no cover - direct execution bootstrap
    raise SystemExit(main())
