"""Pure release-version policy for VNtyper release automation."""

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

REQUIRED_CHECK_NAMES: tuple[str, ...] = (
    "Lint (Ruff)",
    "Type Check (mypy)",
    "Unit Tests (Python 3.10)",
    "Unit Tests (Python 3.11)",
    "Unit Tests (Python 3.12)",
    "Unit Tests (Python 3.13)",
    "Docs build (strict)",
    "CI Success",
    "Build and test image",
    "Docker Success",
)


@dataclass(frozen=True)
class ReleaseVersion:
    """A strict release tag and its numeric components."""

    tag: str
    plain: str
    major: int
    minor: int
    patch: int


@dataclass(frozen=True)
class CheckVerdict:
    """Release verdict for one required workflow check."""

    name: str
    state: Literal["missing", "pending", "success", "failure"]
    conclusion: str | None
    check_run_id: int | None
    url: str | None
    reason: str


@dataclass(frozen=True)
class CheckPoll:
    """Aggregate result from one release-gate polling attempt."""

    action: Literal["wait", "success", "fail", "timeout"]
    verdicts: tuple[CheckVerdict, ...]
    attempt: int
    elapsed_seconds: int


@dataclass(frozen=True)
class AliasState:
    """Observed digest and package version for one registry alias."""

    digest: str
    version: str | None


@dataclass(frozen=True)
class AliasUpdate:
    """Planned decision for one release alias."""

    alias: str
    decision: Literal["create", "advance", "no-op", "skip-newer", "skip-unorderable", "fail-conflict"]
    execute: bool
    reason: str


def parse_release_tag(tag: str) -> ReleaseVersion:
    """Parse a strict VNtyper release tag.

    Args:
        tag: Candidate tag to parse.

    Returns:
        The parsed release version.

    Raises:
        ValueError: If the tag is not strict ``vMAJOR.MINOR.PATCH`` syntax.
    """
    match = re.fullmatch(r"v(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)", tag)
    if match is None:
        msg = f"Release tag {tag!r} must be strict vMAJOR.MINOR.PATCH with no prerelease or build suffix."
        raise ValueError(msg)

    major, minor, patch = (int(component) for component in match.groups())
    return ReleaseVersion(tag=tag, plain=tag[1:], major=major, minor=minor, patch=patch)


def required_aliases(version: ReleaseVersion) -> tuple[str, ...]:
    """Return the container aliases required for a release version.

    Args:
        version: Parsed release version.

    Returns:
        Aliases in promotion order, from immutable release tag to ``latest``.
    """
    return (
        version.tag,
        version.plain,
        f"{version.major}.{version.minor}",
        str(version.major),
        "latest",
    )


def classify_check_runs(
    candidate_sha: str,
    check_runs: Sequence[Mapping[str, object]],
    *,
    attempt: int,
    max_attempts: int = 120,
) -> CheckPoll:
    """Classify exact-SHA GitHub Actions checks for release gating.

    Args:
        candidate_sha: Full commit SHA whose checks are required.
        check_runs: Check Runs API records to classify.
        attempt: One-based polling attempt number.
        max_attempts: Last attempt before a pending poll times out.

    Returns:
        The aggregate action and ordered verdict for every required check.

    Raises:
        ValueError: If the polling bounds are not positive integers or the attempt is outside them.
    """
    if type(max_attempts) is not int or max_attempts < 1:
        msg = "max_attempts must be a positive integer."
        raise ValueError(msg)
    if type(attempt) is not int or attempt < 1 or attempt > max_attempts:
        msg = "attempt must be a one-based integer that does not exceed max_attempts."
        raise ValueError(msg)

    latest_runs: dict[str, Mapping[str, object]] = {}
    latest_ids: dict[str, int] = {}
    for check_run in check_runs:
        if check_run.get("head_sha") != candidate_sha:
            continue
        app = check_run.get("app")
        if not isinstance(app, Mapping) or app.get("slug") != "github-actions":
            continue
        name = check_run.get("name")
        check_run_id = check_run.get("id")
        if not isinstance(name, str) or name not in REQUIRED_CHECK_NAMES or type(check_run_id) is not int:
            continue
        if check_run_id > latest_ids.get(name, -1):
            latest_runs[name] = check_run
            latest_ids[name] = check_run_id

    verdicts: list[CheckVerdict] = []
    for name in REQUIRED_CHECK_NAMES:
        selected_run = latest_runs.get(name)
        if selected_run is None:
            verdicts.append(
                CheckVerdict(
                    name=name,
                    state="missing",
                    conclusion=None,
                    check_run_id=None,
                    url=None,
                    reason="No eligible exact-SHA GitHub Actions check run was found.",
                )
            )
            continue

        status_value = selected_run.get("status")
        status = status_value if isinstance(status_value, str) else None
        conclusion_value = selected_run.get("conclusion")
        conclusion = conclusion_value if isinstance(conclusion_value, str) else None
        url_value = selected_run.get("html_url")
        url = url_value if isinstance(url_value, str) else None
        check_run_id = latest_ids[name]

        if status != "completed":
            state: Literal["pending", "success", "failure"] = "pending"
            reason = f"Newest eligible check run has status {status!r}."
        elif conclusion == "success":
            state = "success"
            reason = "Newest eligible check run completed successfully."
        else:
            state = "failure"
            reason = f"Newest eligible check run completed with conclusion {conclusion!r}."
        verdicts.append(
            CheckVerdict(
                name=name,
                state=state,
                conclusion=conclusion,
                check_run_id=check_run_id,
                url=url,
                reason=reason,
            )
        )

    frozen_verdicts = tuple(verdicts)
    if any(verdict.state == "failure" for verdict in frozen_verdicts):
        action: Literal["wait", "success", "fail", "timeout"] = "fail"
    elif all(verdict.state == "success" for verdict in frozen_verdicts):
        action = "success"
    elif attempt >= max_attempts:
        action = "timeout"
    else:
        action = "wait"
    return CheckPoll(
        action=action,
        verdicts=frozen_verdicts,
        attempt=attempt,
        elapsed_seconds=(attempt - 1) * 30,
    )


def plan_alias_updates(
    version: ReleaseVersion,
    source_digest: str,
    current_aliases: Mapping[str, AliasState | None],
    *,
    dry_run: bool,
) -> tuple[AliasUpdate, ...]:
    """Plan immutable and monotonic release alias updates.

    Args:
        version: Parsed candidate release version.
        source_digest: Verified immutable source manifest digest.
        current_aliases: Current state for exactly the five required aliases.
        dry_run: Whether to suppress all executable writes.

    Returns:
        Alias decisions in exact-first promotion order.

    Raises:
        ValueError: If the digest is not canonical or the alias key set is incomplete.
    """
    if re.fullmatch(r"sha256:[0-9a-f]{64}", source_digest) is None:
        msg = f"Source digest {source_digest!r} must be canonical sha256:<64 lowercase hex characters>."
        raise ValueError(msg)

    aliases = required_aliases(version)
    if set(current_aliases) != set(aliases):
        msg = "Current aliases must contain exactly the five required release alias keys."
        raise ValueError(msg)

    candidate_parts = (version.major, version.minor, version.patch)
    updates: list[AliasUpdate] = []
    for index, alias in enumerate(aliases):
        current = current_aliases[alias]
        if current is None:
            decision: Literal["create", "advance", "no-op", "skip-newer", "skip-unorderable", "fail-conflict"] = (
                "create"
            )
            reason = "Alias is absent and will be created."
        elif index < 2:
            if current.digest == source_digest:
                decision = "no-op"
                reason = "Exact alias already resolves to the source digest."
            else:
                decision = "fail-conflict"
                reason = "Exact alias resolves to a different digest and is immutable."
        else:
            observed_version: ReleaseVersion | None = None
            if current.version is not None:
                try:
                    observed_version = parse_release_tag(f"v{current.version}")
                except ValueError:
                    observed_version = None

            if observed_version is None:
                decision = "skip-unorderable"
                reason = f"Existing alias version {current.version!r} is not strict MAJOR.MINOR.PATCH."
            else:
                observed_parts = (observed_version.major, observed_version.minor, observed_version.patch)
                if observed_parts < candidate_parts:
                    decision = "advance"
                    reason = f"Existing version {current.version!r} is older than candidate {version.plain!r}."
                elif observed_parts > candidate_parts:
                    decision = "skip-newer"
                    reason = f"Existing version {current.version!r} is newer than candidate {version.plain!r}."
                elif current.digest == source_digest:
                    decision = "no-op"
                    reason = "Alias already resolves to the candidate version and source digest."
                else:
                    decision = "fail-conflict"
                    reason = "Alias reports the candidate version but resolves to a different digest."

        execute = decision in {"create", "advance"} and not dry_run
        updates.append(AliasUpdate(alias=alias, decision=decision, execute=execute, reason=reason))
    return tuple(updates)
