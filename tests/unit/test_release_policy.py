import pytest

from scripts.release_policy import (
    REQUIRED_CHECK_NAMES,
    AliasState,
    classify_check_runs,
    parse_release_tag,
    plan_alias_updates,
    required_aliases,
)

pytestmark = pytest.mark.unit

SHA = "a" * 40


def _check(
    name: str,
    *,
    sha: str = SHA,
    run_id: int = 1,
    status: str = "completed",
    conclusion: str | None = "success",
    app: str = "github-actions",
) -> dict[str, object]:
    return {
        "name": name,
        "head_sha": sha,
        "id": run_id,
        "status": status,
        "conclusion": conclusion,
        "html_url": f"https://example.test/check/{run_id}",
        "app": {"slug": app},
    }


def test_release_tag_produces_all_five_aliases_in_promotion_order() -> None:
    version = parse_release_tag("v2.10.3")
    assert (version.plain, version.major, version.minor, version.patch) == ("2.10.3", 2, 10, 3)
    assert required_aliases(version) == ("v2.10.3", "2.10.3", "2.10", "2", "latest")


@pytest.mark.parametrize(
    "tag",
    ("2.0.10", "v2.0", "v2.0.10-rc1", "v2.0.10+local", "v02.0.10", "v2.00.10", "v2.0.010", "v2.0.10;echo pwned"),
)
def test_release_tag_rejects_every_non_strict_form(tag: str) -> None:
    with pytest.raises(ValueError, match="strict vMAJOR.MINOR.PATCH"):
        parse_release_tag(tag)


def test_all_ten_latest_exact_sha_github_actions_checks_are_required() -> None:
    runs = [_check(name) for name in REQUIRED_CHECK_NAMES]
    poll = classify_check_runs(SHA, runs, attempt=1)
    assert poll.action == "success"
    assert tuple(verdict.name for verdict in poll.verdicts) == REQUIRED_CHECK_NAMES
    assert (poll.attempt, poll.elapsed_seconds) == (1, 0)


def test_newest_same_name_check_run_controls_the_verdict() -> None:
    successful = [_check(name) for name in REQUIRED_CHECK_NAMES]
    poll = classify_check_runs(
        SHA,
        successful + [_check("Docker Success", run_id=99, conclusion="failure")],
        attempt=1,
    )
    docker = poll.verdicts[-1]
    assert poll.action == "fail"
    assert (docker.state, docker.conclusion, docker.check_run_id) == ("failure", "failure", 99)


def test_wrong_sha_and_external_checks_cannot_satisfy_a_missing_check() -> None:
    runs = [_check(name) for name in REQUIRED_CHECK_NAMES if name != "Docs build (strict)"]
    runs.extend(
        (
            _check("Docs build (strict)", sha="b" * 40, run_id=98),
            _check("Docs build (strict)", app="external-ci", run_id=99),
        )
    )
    poll = classify_check_runs(SHA, runs, attempt=119)
    docs = poll.verdicts[6]
    assert poll.action == "wait"
    assert (docs.state, docs.check_run_id, docs.url) == ("missing", None, None)


def test_missing_check_times_out_at_the_polling_bound() -> None:
    runs = [_check(name) for name in REQUIRED_CHECK_NAMES if name != "Docs build (strict)"]
    poll = classify_check_runs(SHA, runs, attempt=120)
    assert poll.action == "timeout"
    assert (poll.attempt, poll.elapsed_seconds) == (120, 3570)


@pytest.mark.parametrize("conclusion", ("skipped", "neutral"))
def test_completed_non_success_component_is_a_terminal_failure(conclusion: str) -> None:
    runs = [
        _check(name, conclusion=conclusion if name == "Build and test image" else "success")
        for name in REQUIRED_CHECK_NAMES
    ]
    poll = classify_check_runs(SHA, runs, attempt=1)
    build = poll.verdicts[8]
    assert poll.action == "fail"
    assert (build.state, build.conclusion) == ("failure", conclusion)


def test_incomplete_check_remains_pending() -> None:
    runs = [
        _check(name, status="in_progress", conclusion=None) if name == "CI Success" else _check(name)
        for name in REQUIRED_CHECK_NAMES
    ]
    poll = classify_check_runs(SHA, runs, attempt=1)
    aggregate = poll.verdicts[7]
    assert poll.action == "wait"
    assert (aggregate.state, aggregate.conclusion, aggregate.check_run_id) == ("pending", None, 1)


def test_exact_aliases_create_noop_or_fail_but_never_advance() -> None:
    version = parse_release_tag("v2.0.10")
    current = dict.fromkeys(required_aliases(version))
    source = "sha256:" + "a" * 64
    current["v2.0.10"] = AliasState(source, "2.0.10")
    current["2.0.10"] = AliasState("sha256:" + "b" * 64, "2.0.10")
    updates = plan_alias_updates(version, source, current, dry_run=False)
    assert [(update.alias, update.decision, update.execute) for update in updates[:2]] == [
        ("v2.0.10", "no-op", False),
        ("2.0.10", "fail-conflict", False),
    ]


def test_floating_aliases_create_advance_or_skip_newer_without_downgrade() -> None:
    version = parse_release_tag("v2.0.10")
    source = "sha256:" + "a" * 64
    current = dict.fromkeys(required_aliases(version))
    current["2"] = AliasState("sha256:" + "b" * 64, "2.0.9")
    current["latest"] = AliasState("sha256:" + "c" * 64, "2.1.0")
    updates = plan_alias_updates(version, source, current, dry_run=False)
    assert [(update.alias, update.decision, update.execute) for update in updates] == [
        ("v2.0.10", "create", True),
        ("2.0.10", "create", True),
        ("2.0", "create", True),
        ("2", "advance", True),
        ("latest", "skip-newer", False),
    ]


@pytest.mark.parametrize("observed_version", (None, "2.0", "v2.0.9"))
def test_floating_alias_with_unorderable_version_is_skipped(observed_version: str | None) -> None:
    version = parse_release_tag("v2.0.10")
    source = "sha256:" + "a" * 64
    current = dict.fromkeys(required_aliases(version))
    current["2.0"] = AliasState("sha256:" + "b" * 64, observed_version)
    update = plan_alias_updates(version, source, current, dry_run=False)[2]
    assert (update.alias, update.decision, update.execute) == ("2.0", "skip-unorderable", False)


def test_equal_floating_version_with_a_different_digest_is_a_hard_conflict() -> None:
    version = parse_release_tag("v2.0.10")
    source = "sha256:" + "a" * 64
    current = dict.fromkeys(required_aliases(version))
    current["latest"] = AliasState("sha256:" + "b" * 64, "2.0.10")
    latest = plan_alias_updates(version, source, current, dry_run=False)[-1]
    assert (latest.decision, latest.execute) == ("fail-conflict", False)


def test_each_completed_alias_prefix_converges_to_noop_on_rerun() -> None:
    version = parse_release_tag("v2.0.10")
    source = "sha256:" + "a" * 64
    aliases = ("v2.0.10", "2.0.10", "2.0", "2", "latest")
    for prefix_length in range(len(aliases) + 1):
        current = dict.fromkeys(aliases)
        for alias in aliases[:prefix_length]:
            current[alias] = AliasState(source, "2.0.10")
        updates = plan_alias_updates(version, source, current, dry_run=False)
        expected = [
            (alias, "no-op", False) if index < prefix_length else (alias, "create", True)
            for index, alias in enumerate(aliases)
        ]
        assert [(update.alias, update.decision, update.execute) for update in updates] == expected


def test_dry_run_preserves_decisions_but_disables_every_write() -> None:
    version = parse_release_tag("v2.0.10")
    source = "sha256:" + "a" * 64
    current = dict.fromkeys(required_aliases(version))
    current["2"] = AliasState("sha256:" + "b" * 64, "2.0.9")
    current["latest"] = AliasState("sha256:" + "c" * 64, "2.1.0")
    live = plan_alias_updates(version, source, current, dry_run=False)
    dry = plan_alias_updates(version, source, current, dry_run=True)
    assert [(update.alias, update.decision) for update in dry] == [
        ("v2.0.10", "create"),
        ("2.0.10", "create"),
        ("2.0", "create"),
        ("2", "advance"),
        ("latest", "skip-newer"),
    ]
    assert [(update.alias, update.decision) for update in dry] == [(update.alias, update.decision) for update in live]
    assert all(update.execute is False for update in dry)


@pytest.mark.parametrize("source", ("not-a-digest", "sha256:" + "a" * 63, "sha256:" + "A" * 64))
def test_alias_plan_rejects_noncanonical_source_digest(source: str) -> None:
    version = parse_release_tag("v2.0.10")
    current = dict.fromkeys(required_aliases(version))
    with pytest.raises(ValueError, match="sha256"):
        plan_alias_updates(version, source, current, dry_run=False)


def test_alias_plan_requires_exactly_all_five_alias_keys() -> None:
    version = parse_release_tag("v2.0.10")
    source = "sha256:" + "a" * 64
    current = dict.fromkeys(required_aliases(version))
    current.pop("latest")
    current["main"] = None
    with pytest.raises(ValueError, match="exactly"):
        plan_alias_updates(version, source, current, dry_run=False)
