"""Parsed and executable contracts for the PyPI environment preflight."""

import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest
import yaml

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_PATH = ROOT / ".github" / "workflows" / "publish-pypi.yml"
ISSUE_URL = "https://github.com/hassansaei/VNtyper/issues/236"
ENVIRONMENT_ENDPOINT = "repos/hassansaei/VNtyper/environments/pypi"
POLICIES_ENDPOINT = f"{ENVIRONMENT_ENDPOINT}/deployment-branch-policies"
CUSTOM_RULES_ENDPOINT = f"{ENVIRONMENT_ENDPOINT}/deployment_protection_rules"
MAINTAINER_GUIDANCE_PATHS = (
    ROOT / "AGENTS.md",
    ROOT / "docs" / "superpowers" / "specs" / "2026-08-11-automated-pypi-oidc-design.md",
    ROOT / "docs" / "superpowers" / "plans" / "2026-08-11-automated-pypi-oidc-plan.md",
)

VALID_ENVIRONMENT: dict[str, object] = {
    "name": "pypi",
    "protection_rules": [{"id": 23, "node_id": "BR_x", "type": "branch_policy"}],
    "deployment_branch_policy": {"protected_branches": False, "custom_branch_policies": True},
}
VALID_POLICIES: dict[str, object] = {
    "total_count": 1,
    "branch_policies": [{"id": 17, "node_id": "BP_x", "name": "main", "type": "branch"}],
}
VALID_CUSTOM_RULES: dict[str, object] = {
    "total_count": 0,
    "custom_deployment_protection_rules": [],
}


def _workflow() -> dict[str, Any]:
    """Load the release workflow without YAML 1.1 coercing the ``on`` key."""
    return yaml.load(WORKFLOW_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)


def _preflight() -> dict[str, Any]:
    return next(
        step for step in _workflow()["jobs"]["validate-release"]["steps"] if step.get("id") == "pypi-environment"
    )


def _write_executable(path: Path, source: str) -> None:
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)


def _run_preflight(
    tmp_path: Path,
    payloads: dict[str, object],
    *,
    failed_endpoint: str = "",
) -> tuple[subprocess.CompletedProcess[str], list[list[str]], Path]:
    controller_scripts = tmp_path / "controller" / "scripts"
    controller_scripts.mkdir(parents=True)
    shutil.copy2(ROOT / "scripts" / "pypi_environment_contract.py", controller_scripts)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    calls_path = tmp_path / "gh-calls.jsonl"
    payloads_path = tmp_path / "payloads.json"
    payloads_path.write_text(json.dumps(payloads), encoding="utf-8")
    _write_executable(
        fake_bin / "gh",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

args = sys.argv[1:]
with pathlib.Path(os.environ["FAKE_GH_CALLS"]).open("a", encoding="utf-8") as handle:
    handle.write(json.dumps(args) + "\\n")
if len(args) != 2 or args[0] != "api":
    raise SystemExit(97)
endpoint = args[1]
if endpoint == os.environ.get("FAKE_GH_FAILED_ENDPOINT"):
    print("simulated API failure", file=sys.stderr)
    raise SystemExit(1)
payloads = json.loads(pathlib.Path(os.environ["FAKE_GH_PAYLOADS"]).read_text(encoding="utf-8"))
if endpoint not in payloads:
    raise SystemExit(98)
print(json.dumps(payloads[endpoint]))
""",
    )
    sentinel = tmp_path / "resolve-ran"
    script = f'{_preflight()["run"]}\nprintf resolve > "{sentinel}"\n'
    completed = subprocess.run(
        ["bash", "-c", script],
        cwd=tmp_path / "controller",
        env=os.environ
        | {
            "FAKE_GH_CALLS": str(calls_path),
            "FAKE_GH_FAILED_ENDPOINT": failed_endpoint,
            "FAKE_GH_PAYLOADS": str(payloads_path),
            "GH_TOKEN": "fake-token",
            "GITHUB_REPOSITORY": "hassansaei/VNtyper",
            "PATH": f"{fake_bin}:{os.environ['PATH']}",
            "RUNNER_TEMP": str(tmp_path),
        },
        check=False,
        capture_output=True,
        text=True,
    )
    calls = [json.loads(line) for line in calls_path.read_text(encoding="utf-8").splitlines()]
    return completed, calls, sentinel


def _valid_payloads() -> dict[str, object]:
    return {
        ENVIRONMENT_ENDPOINT: VALID_ENVIRONMENT,
        POLICIES_ENDPOINT: VALID_POLICIES,
        CUSTOM_RULES_ENDPOINT: VALID_CUSTOM_RULES,
    }


def test_maintainer_guidance_records_the_pypi_environment_contract() -> None:
    """Release guidance must preserve the reviewer-free, main-only OIDC policy."""
    required = (
        "reviewer-free",
        "exact branch `main`",
        "fails before package or registry writes",
        "#236",
        "never reintroduce `PYPI_API_TOKEN`",
    )

    for path in MAINTAINER_GUIDANCE_PATHS:
        normalized_guidance = " ".join(path.read_text(encoding="utf-8").split())
        for phrase in required:
            assert phrase in normalized_guidance, f"{path} is missing: {phrase}"


def test_preflight_is_unprivileged_and_ordered_before_candidate_resolution() -> None:
    """No release identity may be resolved before the live environment is proven safe."""
    workflow = _workflow()
    job = workflow["jobs"]["validate-release"]
    steps = job["steps"]
    preflight = _preflight()
    preflight_index = steps.index(preflight)
    resolve_index = next(index for index, step in enumerate(steps) if step.get("id") == "resolve")
    checkout_indexes = [index for index, step in enumerate(steps) if step.get("uses") == "actions/checkout@v7"]

    assert job["permissions"] == {"actions": "read", "contents": "read"}
    assert preflight["env"] == {
        "GH_TOKEN": "${{ github.token }}",
        "GITHUB_REPOSITORY": "${{ github.repository }}",
    }
    assert "continue-on-error" not in preflight
    assert checkout_indexes == [0, resolve_index + 1]
    assert [0, preflight_index, resolve_index, checkout_indexes[1]] == sorted(
        [0, preflight_index, resolve_index, checkout_indexes[1]]
    )


def test_preflight_fetches_only_the_three_read_contract_endpoints() -> None:
    """The controller must collect exactly the three read-only policy observations."""
    run = _preflight()["run"]
    calls = re.findall(r'^fetch_contract "([^"]+)" "\$[A-Z_]+"$', run, re.MULTILINE)

    assert calls == [
        "repos/${GITHUB_REPOSITORY}/environments/pypi",
        "repos/${GITHUB_REPOSITORY}/environments/pypi/deployment-branch-policies",
        "repos/${GITHUB_REPOSITORY}/environments/pypi/deployment_protection_rules",
    ]
    assert "set -euo pipefail" in run
    assert (
        'python scripts/pypi_environment_contract.py "$ENVIRONMENT_PATH" "$POLICIES_PATH" "$CUSTOM_RULES_PATH"' in run
    )
    assert "|| true" not in run
    assert "/pending_deployments" not in run
    assert "review_pending_deployments" not in run
    assert "--method POST" not in run


def test_preflight_accepts_the_exact_live_environment_contract(tmp_path: Path) -> None:
    """The checked-in shell preflight must accept only the valid three-response contract."""
    completed, calls, sentinel = _run_preflight(tmp_path, _valid_payloads())

    assert completed.returncode == 0, completed.stderr
    assert calls == [
        ["api", ENVIRONMENT_ENDPOINT],
        ["api", POLICIES_ENDPOINT],
        ["api", CUSTOM_RULES_ENDPOINT],
    ]
    assert sentinel.read_text(encoding="utf-8") == "resolve"


@pytest.mark.parametrize(
    ("failed_endpoint", "environment", "policies", "custom_rules"),
    (
        (ENVIRONMENT_ENDPOINT, VALID_ENVIRONMENT, VALID_POLICIES, VALID_CUSTOM_RULES),
        (
            "",
            {**VALID_ENVIRONMENT, "protection_rules": [{"type": "required_reviewers"}]},
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
        ),
        (
            "",
            VALID_ENVIRONMENT,
            {"total_count": 1, "branch_policies": [{"name": "release", "type": "branch"}]},
            VALID_CUSTOM_RULES,
        ),
        (
            "",
            VALID_ENVIRONMENT,
            {"total_count": 1, "branch_policies": [{"name": "main", "type": "tag"}]},
            VALID_CUSTOM_RULES,
        ),
        (
            "",
            VALID_ENVIRONMENT,
            {
                "total_count": 2,
                "branch_policies": [
                    {"name": "main", "type": "branch"},
                    {"name": "release", "type": "branch"},
                ],
            },
            VALID_CUSTOM_RULES,
        ),
        (
            "",
            VALID_ENVIRONMENT,
            VALID_POLICIES,
            {"total_count": 1, "custom_deployment_protection_rules": [{"id": 19}]},
        ),
    ),
    ids=("api-failure", "reviewer-rule", "wrong-branch", "tag-policy", "extra-policy", "custom-rule"),
)
def test_preflight_fails_closed_before_resolve_for_every_api_or_policy_error(
    tmp_path: Path,
    failed_endpoint: str,
    environment: dict[str, object],
    policies: dict[str, object],
    custom_rules: dict[str, object],
) -> None:
    """API and hostile policy observations must stop release resolution with #236 guidance."""
    payloads: dict[str, object] = {
        ENVIRONMENT_ENDPOINT: environment,
        POLICIES_ENDPOINT: policies,
        CUSTOM_RULES_ENDPOINT: custom_rules,
    }

    completed, _calls, sentinel = _run_preflight(tmp_path, payloads, failed_endpoint=failed_endpoint)

    assert completed.returncode != 0
    assert not sentinel.exists()
    assert "#236" in completed.stdout + completed.stderr
    assert ISSUE_URL in completed.stdout + completed.stderr
