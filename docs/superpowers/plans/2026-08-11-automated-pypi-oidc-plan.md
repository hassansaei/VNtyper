# Automated PyPI OIDC Publishing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make VNtyper releases fail early when the live `pypi` environment is not reviewer-free and main-only, while retaining tokenless PyPI Trusted Publishing.

**Architecture:** A new pure Python adapter validates the three GitHub REST payloads that define the live environment, its deployment branch policies, and its custom deployment protection rules. The default-branch release controller downloads those payloads in `validate-release` and invokes the adapter before candidate checkout, package building, GHCR promotion, or PyPI publication. Unit tests own the schema and hostile mutations; workflow tests own permissions, ordering, propagation, and the continued absence of credentials or approval bypasses.

**Tech Stack:** Python 3.10, `argparse`, `json`, GitHub Actions YAML, GitHub REST API, pytest, PyYAML, Ruff, mypy, actionlint.

## Global Constraints

- Keep the GitHub environment named exactly `pypi` and keep the PyPI OIDC publisher tuple `hassansaei/VNtyper`, `publish-pypi.yml`, environment `pypi`.
- Permit the live environment only from the exact branch policy `{name: "main", type: "branch"}`.
- Require the built-in `branch_policy` protection rule and reject every required reviewer, wait timer, custom deployment protection rule, unrestricted policy, missing/extra branch policy, and tag policy.
- Do not add `PYPI_API_TOKEN`, Twine credentials, GitHub App approval, token fallback, environment deletion, tag creation, or tag movement.
- Keep the environment reviewer-free with no timer or environment secrets, custom policies for
  the exact branch `main`, and no custom deployment-protection rules. Its preflight fails
  before package or registry writes and cites #236 on a mismatch; OIDC remains the only
  publisher, so never reintroduce `PYPI_API_TOKEN`. After the first successful OIDC release,
  the owner separately deletes that obsolete repository secret.
- Keep `publish-pypi` as the only job with `id-token: write`; keep GHCR permissions, ten release checks, tag/ancestry validation, immutable evidence, concurrency, and retries unchanged.
- Preserve Python 3.10 compatibility and the repository's 120-character Ruff line length.
- Every behavioral edit follows RED, causal failure, minimal GREEN, refactor, verification, commit, and independent review.

---

### Task 1: Pure live-environment policy validator

**Files:**
- Create: `scripts/pypi_environment_contract.py`
- Create: `tests/unit/test_pypi_environment_contract.py`

**Interfaces:**
- Consumes: decoded GitHub environment, deployment-branch-policy, and custom deployment-protection-rule JSON objects as `Mapping[str, object]`.
- Produces: `validate_pypi_environment(environment: Mapping[str, object], policies: Mapping[str, object], custom_rules: Mapping[str, object]) -> None` and `main(argv: Sequence[str] | None = None) -> int`.

- [ ] **Step 1: Write the successful-policy and hostile-schema tests**

Create a literal valid fixture rather than importing production constants:

```python
VALID_ENVIRONMENT = {
    "name": "pypi",
    "protection_rules": [{"id": 23, "node_id": "BR_x", "type": "branch_policy"}],
    "deployment_branch_policy": {
        "protected_branches": False,
        "custom_branch_policies": True,
    },
}
VALID_POLICIES = {
    "total_count": 1,
    "branch_policies": [{"id": 17, "node_id": "BP_x", "name": "main", "type": "branch"}],
}
VALID_CUSTOM_RULES = {"total_count": 0, "custom_deployment_protection_rules": []}

def test_exact_reviewer_free_main_only_policy_is_valid() -> None:
    assert validate_pypi_environment(VALID_ENVIRONMENT, VALID_POLICIES, VALID_CUSTOM_RULES) is None
```

Parametrize independent mutations for `required_reviewers`, `wait_timer`, missing/duplicate/unknown built-in protection rules, malformed `protection_rules`, `deployment_branch_policy=None`, protected-branch mode, unrestricted mode, zero/two policies, `master`, glob patterns, `{name: "main", type: "tag"}`, nonzero custom-rule count, and any custom-rule row. Assert `ValueError` includes the rejected field and issue URL `https://github.com/hassansaei/VNtyper/issues/236`.

- [ ] **Step 2: Write CLI tests for file, JSON, and policy failures**

Use `tmp_path` JSON files and assert:

```python
assert main([str(environment_path), str(policies_path), str(custom_rules_path)]) == 0
assert "reviewer-free and restricted to branch main" in capsys.readouterr().out
```

Then corrupt each of the three files and assert exit `1`, the exact failed path, an actionable `#236` diagnostic, and no traceback. Test a valid JSON payload with a forbidden reviewer separately from invalid JSON so parsing and policy failures cannot collapse into one branch.

- [ ] **Step 3: Run the new tests to verify RED**

Run: `mamba run -n vntyper pytest -q tests/unit/test_pypi_environment_contract.py`

Expected: collection fails because `scripts.pypi_environment_contract` does not exist.

- [ ] **Step 4: Implement the minimal pure validator and CLI**

Implement strict type helpers that reject booleans where integers are required. Project server-owned identifiers out of the policy row but require the contract fields exactly:

```python
ISSUE_URL = "https://github.com/hassansaei/VNtyper/issues/236"

def validate_pypi_environment(
    environment: Mapping[str, object],
    policies: Mapping[str, object],
    custom_rules: Mapping[str, object],
) -> None:
    if environment.get("name") != "pypi":
        _reject("environment name must be 'pypi'", environment)
    protection_rules = environment.get("protection_rules")
    if not isinstance(protection_rules, list) or len(protection_rules) != 1:
        _reject("protection_rules must contain exactly the built-in branch_policy rule", environment)
    rule = protection_rules[0]
    if not isinstance(rule, Mapping) or rule.get("type") != "branch_policy":
        _reject("protection_rules must contain exactly the built-in branch_policy rule", environment)
    if environment.get("deployment_branch_policy") != {
        "protected_branches": False,
        "custom_branch_policies": True,
    }:
        _reject("deployment_branch_policy must be custom and main-only", environment)
    rows = policies.get("branch_policies")
    if type(policies.get("total_count")) is not int or policies["total_count"] != 1:
        _reject("total_count must be exactly 1", policies)
    if not isinstance(rows, list) or len(rows) != 1:
        _reject("branch_policies must contain exactly one row", policies)
    row = rows[0]
    if not isinstance(row, Mapping) or row.get("name") != "main" or row.get("type") != "branch":
        _reject("the sole policy must be branch main", policies)
    if custom_rules.get("total_count") != 0 or custom_rules.get("custom_deployment_protection_rules") != []:
        _reject("custom deployment protection rules must be empty", custom_rules)
```

Load each of the three inputs with `Path.read_text(encoding="utf-8")` and `json.loads`; require top-level objects; print one success line on `stdout`, and on `OSError`, `json.JSONDecodeError`, or `ValueError` print one `stderr` diagnostic and return `1`.

- [ ] **Step 5: Run focused GREEN and static checks**

Run:

```bash
mamba run -n vntyper pytest -q tests/unit/test_pypi_environment_contract.py
mamba run -n vntyper ruff format --check scripts/pypi_environment_contract.py tests/unit/test_pypi_environment_contract.py
mamba run -n vntyper ruff check scripts/pypi_environment_contract.py tests/unit/test_pypi_environment_contract.py
mamba run -n vntyper mypy scripts/pypi_environment_contract.py tests/unit/test_pypi_environment_contract.py
```

Expected: all commands exit `0`.

- [ ] **Step 6: Commit Task 1**

```bash
git add scripts/pypi_environment_contract.py tests/unit/test_pypi_environment_contract.py
git commit -m "feat(release): validate live PyPI environment policy" -m "Closes #236"
```

### Task 2: Default-branch workflow preflight

**Files:**
- Modify: `.github/workflows/publish-pypi.yml:35-225`
- Create: `tests/unit/test_pypi_environment_workflow.py`
- Modify: `tests/unit/test_release_workflow_contract.py:549-578`

**Interfaces:**
- Consumes: `scripts/pypi_environment_contract.py` from Task 1 and the Actions-provided `github.repository`/`github.token`.
- Produces: a `validate-release` step with id `pypi-environment` that writes no registry/package state and fails before candidate validation on every API or policy error.

- [ ] **Step 1: Write parsed workflow RED tests**

Load the workflow using `yaml.BaseLoader` and assert:

```python
job = workflow["jobs"]["validate-release"]
assert job["permissions"] == {"actions": "read", "contents": "read"}
step = next(item for item in job["steps"] if item.get("id") == "pypi-environment")
assert step["env"] == {
    "GH_TOKEN": "${{ github.token }}",
    "GITHUB_REPOSITORY": "${{ github.repository }}",
}
assert step["continue-on-error"] is absent
```

Pin ordering: controller checkout, `pypi-environment`, `resolve`, candidate checkout. Assert the run script fetches exactly `/repos/${GITHUB_REPOSITORY}/environments/pypi`, `/deployment-branch-policies`, and `/deployment_protection_rules`, passes all three files to `controller/scripts/pypi_environment_contract.py`, uses `set -euo pipefail`, and contains neither `|| true` nor a reviewer-approval API.

- [ ] **Step 2: Write executable workflow RED tests**

Extract the checked-in step and run it with a fake `gh` executable. The fake records argv and returns literal payloads. Assert a valid triple exits `0` and records exactly three GETs. Parametrize API failure, reviewer payload, wrong branch, tag type, extra policy, and custom protection rule; each must exit nonzero before a fake `resolve` sentinel can run and must mention #236.

- [ ] **Step 3: Strengthen privilege and credential regression tests**

Update the existing exact permission assertion for `validate-release` only. Add whole-workflow assertions that:

```python
assert [name for name, job in jobs.items() if job.get("permissions", {}).get("id-token") == "write"] == [
    "publish-pypi"
]
assert "PYPI_API_TOKEN" not in raw
assert "TWINE_PASSWORD" not in raw
assert "/pending_deployments" not in raw
assert "review_pending_deployments" not in raw
assert "--method POST" not in preflight["run"]
```

Also assert `publish-pypi.environment == {"name": "pypi"}` and its permission remains exactly `{"id-token": "write"}`.

- [ ] **Step 4: Run the workflow tests to verify RED**

Run:

```bash
mamba run -n vntyper pytest -q \
  tests/unit/test_pypi_environment_workflow.py \
  tests/unit/test_release_workflow_contract.py::test_release_trigger_uses_default_branch_repository_dispatch_for_production \
  tests/unit/test_release_workflow_contract.py::test_pypi_publish_is_protected_oidc_only_and_rerun_safe
```

Expected: failures identify the absent preflight step and missing `actions: read` permission.

- [ ] **Step 5: Implement the minimal preflight step**

Add `actions: read` to `validate-release`, then insert the step immediately after controller checkout:

```yaml
      - id: pypi-environment
        working-directory: controller
        env:
          GH_TOKEN: ${{ github.token }}
          GITHUB_REPOSITORY: ${{ github.repository }}
        run: |
          set -euo pipefail
          ENVIRONMENT_PATH="$RUNNER_TEMP/pypi-environment.json"
          POLICIES_PATH="$RUNNER_TEMP/pypi-deployment-branch-policies.json"
          CUSTOM_RULES_PATH="$RUNNER_TEMP/pypi-deployment-protection-rules.json"
          fetch_contract() {
            ENDPOINT=$1
            DESTINATION=$2
            TEMPORARY="${DESTINATION}.tmp"
            if ! gh api "$ENDPOINT" > "$TEMPORARY"; then
              rm -f -- "$TEMPORARY"
              echo "::error::Failed to read GitHub API endpoint ${ENDPOINT}; see https://github.com/hassansaei/VNtyper/issues/236"
              return 1
            fi
            mv -- "$TEMPORARY" "$DESTINATION"
          }
          fetch_contract "repos/${GITHUB_REPOSITORY}/environments/pypi" "$ENVIRONMENT_PATH"
          fetch_contract "repos/${GITHUB_REPOSITORY}/environments/pypi/deployment-branch-policies" "$POLICIES_PATH"
          fetch_contract "repos/${GITHUB_REPOSITORY}/environments/pypi/deployment_protection_rules" "$CUSTOM_RULES_PATH"
          python scripts/pypi_environment_contract.py "$ENVIRONMENT_PATH" "$POLICIES_PATH" "$CUSTOM_RULES_PATH"
```

Do not add an `if`, soft-failure path, approval API, environment reference, or write permission to this job.

- [ ] **Step 6: Run focused GREEN, mutations, and workflow checks**

Run:

```bash
mamba run -n vntyper pytest -q \
  tests/unit/test_pypi_environment_contract.py \
  tests/unit/test_pypi_environment_workflow.py \
  tests/unit/test_release_workflow_contract.py
mamba run -n vntyper actionlint .github/workflows/publish-pypi.yml
mamba run -n vntyper ruff format --check scripts/pypi_environment_contract.py tests/unit/test_pypi_environment_contract.py tests/unit/test_pypi_environment_workflow.py
mamba run -n vntyper ruff check scripts/pypi_environment_contract.py tests/unit/test_pypi_environment_contract.py tests/unit/test_pypi_environment_workflow.py
mamba run -n vntyper mypy scripts/pypi_environment_contract.py tests/unit/test_pypi_environment_contract.py tests/unit/test_pypi_environment_workflow.py
```

Then make these one-at-a-time temporary mutations and require the focused suite to fail before restoring each: remove `actions: read`; move the preflight after candidate checkout; add `continue-on-error`; change `main` to `v*`; change policy type to `tag`; add `PYPI_API_TOKEN`; remove `environment: pypi`.

- [ ] **Step 7: Commit Task 2**

```bash
git add .github/workflows/publish-pypi.yml tests/unit/test_pypi_environment_workflow.py tests/unit/test_release_workflow_contract.py
git commit -m "ci(release): fail fast on unsafe PyPI environment" -m "Closes #236"
```

### Task 3: Maintainer contract, full verification, review, and PR

**Files:**
- Modify: `AGENTS.md` release workflow section
- Modify: `docs/development/ci-followups.md`
- Modify: `docs/superpowers/specs/2026-08-11-automated-pypi-oidc-design.md`
- Modify: `docs/superpowers/plans/2026-08-11-automated-pypi-oidc-plan.md`

**Interfaces:**
- Consumes: Tasks 1-2 and live issue #236.
- Produces: executable maintainer guidance, a reviewed feature branch, and a PR against `main`.

- [ ] **Step 1: Write documentation-contract RED assertions**

In `tests/unit/test_pypi_environment_workflow.py`, assert normalized guidance contains all of:

```python
required = (
    "reviewer-free",
    "exact branch `main`",
    "fails before package or registry writes",
    "#236",
    "never reintroduce `PYPI_API_TOKEN`",
)
```

Require the same selected policy in AGENTS, the design, and the plan. Run the test and record the missing phrases as RED.

- [ ] **Step 2: Update the maintainer guidance minimally**

State that the `pypi` environment must have zero reviewers/timer, custom branch policies, exactly branch `main`, and no environment secrets; current controller preflights that state; mismatches point to #236 and fail before writes; OIDC remains the only publisher; after first successful OIDC release the owner separately deletes the obsolete repository secret.

- [ ] **Step 3: Run documentation and focused GREEN**

Run:

```bash
mamba run -n vntyper pytest -q tests/unit/test_pypi_environment_workflow.py
mamba run -n vntyper make docs-build
```

Expected: both exit `0`.

- [ ] **Step 4: Run all repository gates**

Run sequentially:

```bash
mamba run -n vntyper make format-check
mamba run -n vntyper make lint
mamba run -n vntyper make type-check-all
mamba run -n vntyper make test-unit
mamba run -n vntyper make test-unit-cov
mamba run -n vntyper make test-scripts-cov
mamba run -n vntyper make patch-coverage
mamba run -n vntyper make check-all
mamba run -n vntyper make ci-local
```

Because `.github/workflows/` changed, `ci-local` is mandatory. Do not change any threshold, marker, tier, tolerance, or expected value to make a gate pass; apply systematic debugging to every unexpected failure.

- [ ] **Step 5: Inspect the complete diff and commit documentation**

Run:

```bash
git diff --check
git diff --stat origin/main...HEAD
git diff origin/main...HEAD
git status --short
```

Then commit:

```bash
git add AGENTS.md docs/development/ci-followups.md \
  docs/superpowers/specs/2026-08-11-automated-pypi-oidc-design.md \
  docs/superpowers/plans/2026-08-11-automated-pypi-oidc-plan.md
git commit -m "docs(release): record automated OIDC environment" -m "Closes #236"
```

- [ ] **Step 6: Request and receive independent code review**

Review `origin/main...HEAD` against the approved design. Require explicit Critical/Important/Minor findings, security review of permissions and ordering, Python 3.10 review, mutation sensitivity, and a PASS/NEEDS_FIXES verdict. Apply every valid finding with `superpowers:receiving-code-review`, TDD, and a separate follow-up commit, then re-review until Critical=0 and Important=0.

- [ ] **Step 7: Re-run completion verification and publish the PR**

Run fresh `make check-all`, `make ci-local`, `git diff --check`, and clean-status checks. Then use `github:yeet` to push without force and open a PR against `main` with `Closes #236`, a security rationale, RED/GREEN evidence, mutation verdicts, and all gate results. Monitor required GitHub checks to terminal green; address actionable review feedback before merge.

### Task 4: Live migration and ordered release recovery

**Files:**
- External state only: GitHub environment `pypi`, Actions runs, PyPI, GHCR, GitHub Releases.

**Interfaces:**
- Consumes: merged Task 3 PR, `@hassansaei` Administration-write action from #236, signed tags `v2.0.11` and `v2.0.12`.
- Produces: PyPI 2.0.11 followed by 2.0.12, verified GHCR aliases, GitHub Releases, closed issue #236, and removal of only stale completed worktrees/branches.

- [ ] **Step 1: Verify the administrator migration before dispatch**

Read the environment, branch-policy, and custom-protection-rule endpoints. Require exactly the built-in `branch_policy` protection rule, custom policies enabled, exactly `{name: "main", type: "branch"}`, and zero custom deployment protection rules. Require issue #236 to record who changed it and when. Do not dispatch while the policy is wrong.

- [ ] **Step 2: Complete v2.0.11 first**

Inspect Actions run `31421072639`. If GitHub resumes it, monitor to terminal success. If it remains waiting after the policy change, have an administrator select **Start all waiting jobs**, or send one authenticated `vntyper_release` dispatch with `client_payload.tag=v2.0.11`. Verify PyPI JSON/file metadata exposes 2.0.11 and the release summary is successful. On `invalid-publisher`, ask `@hassansaei` to verify the PyPI tuple; do not add a token fallback.

- [ ] **Step 3: Re-verify the existing v2.0.12 tag and dry run**

Require tag object `f41b3aec14c37d19d01a64a2de136c68ec9bc88d`, good SSH signature, peeled commit `7f8583dd60565fec5d0297cbb26ccd2d7f439b22`, and successful no-write run `31463886477`. Never move or force-push it.

- [ ] **Step 4: Dispatch and monitor v2.0.12 production**

Send authenticated repository dispatch JSON `{"event_type":"vntyper_release","client_payload":{"tag":"v2.0.12"}}`. Monitor validation, gates, build, GHCR promotion, PyPI publication, and summary to terminal success. Do not start a second dispatch while the first is pending.

- [ ] **Step 5: Verify every public release surface**

Require PyPI 2.0.12 wheel and sdist, GHCR immutable aliases `v2.0.12`/`2.0.12`, floating aliases `2.0`/`2`/`latest`, a single evidence-bound digest, and OCI revision/version labels equal to the full release SHA and `2.0.12`. Require the GitHub Release to reference the existing signed tag and contain no substituted source archive claims.

- [ ] **Step 6: Close and clean up**

After the first successful OIDC publish, comment on #236 with run/PyPI/GHCR evidence and ask `@hassansaei` to delete repository secret `PYPI_API_TOKEN`; verify deletion when permission permits. Close #236 only after the live policy and both versions are verified. Delete only branches/worktrees already merged or completed; preserve unrelated work. Never delete or move release tags.
