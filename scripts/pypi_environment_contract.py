"""Validate the live GitHub PyPI environment's release policy."""

from __future__ import annotations

import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import NoReturn

ISSUE_URL = "https://github.com/hassansaei/VNtyper/issues/236"


def _reject(message: str) -> NoReturn:
    """Raise a policy error with the corrective issue reference.

    Args:
        message: The contract clause the live GitHub configuration violated.

    Raises:
        ValueError: Always, with the project issue that documents the required policy.
    """
    raise ValueError(f"{message}; see {ISSUE_URL}")


def _has_exact_count(payload: Mapping[str, object], expected: int) -> bool:
    """Return whether a payload has an integer ``total_count`` equal to ``expected``.

    Args:
        payload: Decoded GitHub API response containing a count.
        expected: Required count for the endpoint contract.

    Returns:
        ``True`` only for an actual integer count with the requested value.
    """
    return type(payload.get("total_count")) is int and payload["total_count"] == expected


def validate_pypi_environment(
    environment: Mapping[str, object], policies: Mapping[str, object], custom_rules: Mapping[str, object]
) -> None:
    """Fail unless the live PyPI environment is reviewer-free and main-only.

    Args:
        environment: Decoded GitHub environment endpoint payload.
        policies: Decoded deployment-branch-policy endpoint payload.
        custom_rules: Decoded custom deployment-protection-rule endpoint payload.

    Raises:
        ValueError: If a protection rule, branch policy, or custom rule violates #236.
    """
    if environment.get("name") != "pypi":
        _reject("environment name must be 'pypi'")

    protection_rules = environment.get("protection_rules")
    if not isinstance(protection_rules, list) or len(protection_rules) != 1:
        _reject("protection_rules must contain exactly the built-in branch_policy rule")
    rule = protection_rules[0]
    if not isinstance(rule, Mapping) or rule.get("type") != "branch_policy":
        _reject("protection_rules must contain exactly the built-in branch_policy rule")

    if environment.get("deployment_branch_policy") != {
        "protected_branches": False,
        "custom_branch_policies": True,
    }:
        _reject("deployment_branch_policy must be custom and main-only")

    if not _has_exact_count(policies, 1):
        _reject("total_count must be exactly 1")
    rows = policies.get("branch_policies")
    if not isinstance(rows, list) or len(rows) != 1:
        _reject("branch_policies must contain exactly one row")
    row = rows[0]
    if not isinstance(row, Mapping) or row.get("name") != "main" or row.get("type") != "branch":
        _reject("the sole policy must be branch main")

    if not _has_exact_count(custom_rules, 0) or custom_rules.get("custom_deployment_protection_rules") != []:
        _reject("custom deployment protection rules must be empty")


def _read_object(path: Path) -> Mapping[str, object]:
    """Read one JSON object from an endpoint-response file.

    Args:
        path: JSON file captured from a GitHub API endpoint.

    Returns:
        The decoded top-level JSON object.

    Raises:
        OSError: If the input file cannot be read.
        json.JSONDecodeError: If the input is not valid JSON.
        ValueError: If the decoded JSON is not an object.
    """
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError("top-level JSON value must be an object")
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    """Validate environment-policy API responses supplied as JSON files.

    Args:
        argv: Three JSON paths: environment, branch policies, and custom rules.

    Returns:
        Zero when the live policy meets #236, otherwise one.
    """
    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) != 3:
        print(
            f"PyPI environment policy validation failed (#236): expected three JSON file paths; see {ISSUE_URL}",
            file=sys.stderr,
        )
        return 1

    paths = [Path(argument) for argument in arguments]
    payloads: list[Mapping[str, object]] = []
    for path in paths:
        try:
            payloads.append(_read_object(path))
        except (OSError, json.JSONDecodeError, ValueError) as error:
            print(f"PyPI environment policy validation failed for {path} (#236): {error}", file=sys.stderr)
            return 1

    try:
        environment, policies, custom_rules = payloads
        validate_pypi_environment(environment, policies, custom_rules)
    except ValueError as error:
        print(f"PyPI environment policy validation failed (#236): {error}", file=sys.stderr)
        return 1

    print("PyPI environment is reviewer-free and restricted to branch main.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
