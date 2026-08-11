"""Contracts for the live GitHub PyPI environment policy validator."""

import importlib
import json
from collections.abc import Callable
from copy import deepcopy
from pathlib import Path
from typing import cast

import pytest

pytestmark = pytest.mark.unit

contract = importlib.import_module("scripts.pypi_environment_contract")
ISSUE_URL = cast(str, contract.ISSUE_URL)
main = cast(Callable[[list[str]], int], contract.main)
validate_pypi_environment = cast(
    Callable[[dict[str, object], dict[str, object], dict[str, object]], None], contract.validate_pypi_environment
)


VALID_ENVIRONMENT: dict[str, object] = {
    "name": "pypi",
    "protection_rules": [{"id": 23, "node_id": "BR_x", "type": "branch_policy"}],
    "deployment_branch_policy": {
        "protected_branches": False,
        "custom_branch_policies": True,
    },
}
VALID_POLICIES: dict[str, object] = {
    "total_count": 1,
    "branch_policies": [{"id": 17, "node_id": "BP_x", "name": "main", "type": "branch"}],
}
VALID_CUSTOM_RULES: dict[str, object] = {"total_count": 0, "custom_deployment_protection_rules": []}


def test_exact_reviewer_free_main_only_policy_is_valid() -> None:
    """Only the reviewer-free, main-only deployment policy is acceptable."""
    assert validate_pypi_environment(VALID_ENVIRONMENT, VALID_POLICIES, VALID_CUSTOM_RULES) is None


@pytest.mark.parametrize(
    ("environment", "policies", "custom_rules", "field"),
    (
        (
            {**VALID_ENVIRONMENT, "protection_rules": [{"type": "required_reviewers"}]},
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "protection_rules",
        ),
        (
            {**VALID_ENVIRONMENT, "protection_rules": [{"type": "wait_timer"}]},
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "protection_rules",
        ),
        (
            {key: value for key, value in VALID_ENVIRONMENT.items() if key != "protection_rules"},
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "protection_rules",
        ),
        (
            {
                **VALID_ENVIRONMENT,
                "protection_rules": [
                    {"id": 23, "node_id": "BR_x", "type": "branch_policy"},
                    {"id": 24, "node_id": "BR_y", "type": "branch_policy"},
                ],
            },
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "protection_rules",
        ),
        (
            {**VALID_ENVIRONMENT, "protection_rules": [{"type": "unknown"}]},
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "protection_rules",
        ),
        ({**VALID_ENVIRONMENT, "protection_rules": {}}, VALID_POLICIES, VALID_CUSTOM_RULES, "protection_rules"),
        (
            {**VALID_ENVIRONMENT, "deployment_branch_policy": None},
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "deployment_branch_policy",
        ),
        (
            {
                **VALID_ENVIRONMENT,
                "deployment_branch_policy": {"protected_branches": True, "custom_branch_policies": False},
            },
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "deployment_branch_policy",
        ),
        (
            {
                **VALID_ENVIRONMENT,
                "deployment_branch_policy": {"protected_branches": False, "custom_branch_policies": False},
            },
            VALID_POLICIES,
            VALID_CUSTOM_RULES,
            "deployment_branch_policy",
        ),
        (VALID_ENVIRONMENT, {"total_count": 0, "branch_policies": []}, VALID_CUSTOM_RULES, "total_count"),
        (
            VALID_ENVIRONMENT,
            {
                "total_count": 2,
                "branch_policies": [
                    {"id": 17, "node_id": "BP_x", "name": "main", "type": "branch"},
                    {"id": 18, "node_id": "BP_y", "name": "release", "type": "branch"},
                ],
            },
            VALID_CUSTOM_RULES,
            "total_count",
        ),
        (
            VALID_ENVIRONMENT,
            {"total_count": 1, "branch_policies": [{"id": 17, "node_id": "BP_x", "name": "master", "type": "branch"}]},
            VALID_CUSTOM_RULES,
            "sole policy",
        ),
        (
            VALID_ENVIRONMENT,
            {"total_count": 1, "branch_policies": [{"id": 17, "node_id": "BP_x", "name": "main*", "type": "branch"}]},
            VALID_CUSTOM_RULES,
            "sole policy",
        ),
        (
            VALID_ENVIRONMENT,
            {"total_count": 1, "branch_policies": [{"id": 17, "node_id": "BP_x", "name": "main", "type": "tag"}]},
            VALID_CUSTOM_RULES,
            "sole policy",
        ),
        (
            VALID_ENVIRONMENT,
            {"total_count": True, "branch_policies": [{"name": "main", "type": "branch"}]},
            VALID_CUSTOM_RULES,
            "total_count",
        ),
        (
            VALID_ENVIRONMENT,
            VALID_POLICIES,
            {"total_count": 1, "custom_deployment_protection_rules": []},
            "custom deployment protection rules",
        ),
        (
            VALID_ENVIRONMENT,
            VALID_POLICIES,
            {"total_count": 1, "custom_deployment_protection_rules": [{"id": 19, "type": "required_reviewers"}]},
            "custom deployment protection rules",
        ),
        (
            VALID_ENVIRONMENT,
            VALID_POLICIES,
            {"total_count": False, "custom_deployment_protection_rules": []},
            "custom deployment protection rules",
        ),
    ),
)
def test_hostile_environment_or_policy_schema_is_rejected(
    environment: dict[str, object], policies: dict[str, object], custom_rules: dict[str, object], field: str
) -> None:
    """Every relaxed, malformed, or non-main policy must fail closed."""
    with pytest.raises(ValueError) as error:
        validate_pypi_environment(deepcopy(environment), deepcopy(policies), deepcopy(custom_rules))

    assert field in str(error.value)
    assert ISSUE_URL in str(error.value)


def test_main_accepts_valid_environment_and_policy_files(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """The CLI accepts decoded files whose policy is exactly the release contract."""
    environment_path = tmp_path / "environment.json"
    policies_path = tmp_path / "policies.json"
    custom_rules_path = tmp_path / "custom-rules.json"
    environment_path.write_text(json.dumps(VALID_ENVIRONMENT), encoding="utf-8")
    policies_path.write_text(json.dumps(VALID_POLICIES), encoding="utf-8")
    custom_rules_path.write_text(json.dumps(VALID_CUSTOM_RULES), encoding="utf-8")

    assert main([str(environment_path), str(policies_path), str(custom_rules_path)]) == 0
    assert "reviewer-free and restricted to branch main" in capsys.readouterr().out


@pytest.mark.parametrize("target", ("environment", "policies", "custom_rules"))
def test_main_reports_a_path_and_issue_for_invalid_json(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], target: str
) -> None:
    """Unreadable JSON reports the exact input path without a traceback."""
    environment_path = tmp_path / "environment.json"
    policies_path = tmp_path / "policies.json"
    custom_rules_path = tmp_path / "custom-rules.json"
    environment_path.write_text(json.dumps(VALID_ENVIRONMENT), encoding="utf-8")
    policies_path.write_text(json.dumps(VALID_POLICIES), encoding="utf-8")
    custom_rules_path.write_text(json.dumps(VALID_CUSTOM_RULES), encoding="utf-8")
    paths = {
        "environment": environment_path,
        "policies": policies_path,
        "custom_rules": custom_rules_path,
    }
    failed_path = paths[target]
    failed_path.write_text("{not-json", encoding="utf-8")

    assert main([str(environment_path), str(policies_path), str(custom_rules_path)]) == 1

    diagnostic = capsys.readouterr().err
    assert str(failed_path) in diagnostic
    assert "#236" in diagnostic
    assert "Traceback" not in diagnostic


def test_main_reports_policy_failure_separately_from_json_parsing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A valid payload with reviewers fails policy validation rather than JSON parsing."""
    environment_path = tmp_path / "environment.json"
    policies_path = tmp_path / "policies.json"
    custom_rules_path = tmp_path / "custom-rules.json"
    environment = {**VALID_ENVIRONMENT, "protection_rules": [{"type": "required_reviewers"}]}
    environment_path.write_text(json.dumps(environment), encoding="utf-8")
    policies_path.write_text(json.dumps(VALID_POLICIES), encoding="utf-8")
    custom_rules_path.write_text(json.dumps(VALID_CUSTOM_RULES), encoding="utf-8")

    assert main([str(environment_path), str(policies_path), str(custom_rules_path)]) == 1

    diagnostic = capsys.readouterr().err
    assert "protection_rules" in diagnostic
    assert "#236" in diagnostic
    assert "JSON" not in diagnostic
    assert "Traceback" not in diagnostic
