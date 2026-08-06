"""The lifetime of a cohort alias: when it is claimed, and how long it holds.

An alias is a label a person chooses for a cohort, and the service keeps it
unique: the claim is one Redis key whose existence *is* the claim, taken with a
single atomic `SET NX EX` so two simultaneous creations cannot both believe they
won it. Uniqueness is only as good as the bookkeeping around that key, and two
properties are pinned here because neither follows from the claim itself.

1. **A refused creation must not hold an alias.** The claim is worthless to a
   cohort that was never created, and it lasts a full retention period, so
   claiming it before the request is known to be acceptable takes a name out of
   circulation on behalf of nothing.
2. **A live cohort must not outlive its own alias.** The cohort record and its
   member set are extended every time one of its jobs runs. If the claim is not
   extended with them, an actively used cohort eventually loses its name and a
   second cohort can take it -- so the label stops identifying what it labels.

Authorization itself is `test_cohort_authz.py`'s subject; this module is only
about the alias key. `docker` is put on `sys.path` by
`tests/unit/web/conftest.py`, which pytest imports before this module.
"""

from unittest.mock import MagicMock

import pytest

pytestmark = pytest.mark.unit

from app.cohorts import (  # noqa: E402
    ALIAS_KEY_PREFIX,
    COHORT_KEY_PREFIX,
    MAX_ALIAS_CHARS,
    alias_key,
    cohort_key,
    create_cohort_record,
    extend_cohort_retention,
    resolve_cohort,
)
from app.utils import MAX_PASSPHRASE_BYTES  # noqa: E402

PASSPHRASE = "correct-horse-battery-staple"
ALIAS = "family-a"
RETENTION = 3600


# ---------------------------------------------------------------------------
# A refused creation leaves the alias free
# ---------------------------------------------------------------------------


def test_an_over_long_passphrase_does_not_claim_the_alias(fake_redis) -> None:
    """A passphrase the service cannot store leaves the requested alias free.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    over_long = "p" * (MAX_PASSPHRASE_BYTES + 1)

    with pytest.raises(ValueError, match=str(MAX_PASSPHRASE_BYTES)):
        create_cohort_record(fake_redis, alias=ALIAS, passphrase=over_long, retention_seconds=RETENTION)

    assert fake_redis.get(alias_key(ALIAS)) is None
    assert fake_redis.keys(f"{COHORT_KEY_PREFIX}*") == []


def test_a_missing_passphrase_does_not_claim_the_alias(fake_redis) -> None:
    """A creation with no credential at all leaves the requested alias free.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    with pytest.raises(ValueError, match="passphrase is required"):
        create_cohort_record(fake_redis, alias=ALIAS, passphrase="   ", retention_seconds=RETENTION)

    assert fake_redis.get(alias_key(ALIAS)) is None


def test_a_failure_after_the_claim_releases_it_again(fake_redis) -> None:
    """The claim is given back when the cohort behind it cannot be written.

    The claim has to be taken before the record exists -- that is what makes it
    atomic -- so the only way it can be kept honest is to release it when what it
    was taken for does not happen.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    store = MagicMock(wraps=fake_redis)
    store.pipeline.side_effect = RuntimeError("the store went away mid-write")

    with pytest.raises(RuntimeError, match="went away"):
        create_cohort_record(store, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=RETENTION)

    assert fake_redis.get(alias_key(ALIAS)) is None


class _HalfAppliedPipeline:
    """A pipeline whose queued write lands but whose EXEC never returns.

    The transaction is what normally makes the record and its TTL one step, but a
    client that applied the HSET and then lost the connection before EXEC looks
    exactly like this from the caller's side. It is the state the rollback exists
    for, so it is built here rather than assumed impossible.
    """

    def __init__(self, store) -> None:
        """Record the store the queued write is applied to.

        Args:
            store: The Redis stand-in to write through to.
        """
        self._store = store

    def hset(self, *args, **kwargs):
        """Apply the write immediately, as a server that received it would.

        Args:
            *args: Passed to the store.
            **kwargs: Passed to the store.

        Returns:
            _HalfAppliedPipeline: Self, so calls chain as redis-py's do.
        """
        self._store.hset(*args, **kwargs)
        return self

    def expire(self, *_args, **_kwargs):
        """Queue an expiry that is never applied.

        Args:
            *_args: Ignored.
            **_kwargs: Ignored.

        Returns:
            _HalfAppliedPipeline: Self, so calls chain as redis-py's do.
        """
        return self

    def execute(self):
        """Fail the way a dropped connection does.

        Raises:
            RuntimeError: Always.
        """
        msg = "the store went away between writes"
        raise RuntimeError(msg)


def test_a_failure_after_the_record_leaves_no_unaddressable_cohort(fake_redis) -> None:
    """A cohort whose TTL was never set must not be left behind.

    The identifier is minted inside the creation and only reaches the caller in
    the response, so a creation that raises never reveals it. A hash written
    without its expiry is therefore permanent and unreachable at once: nothing
    can open it, nothing can delete it, and it still holds an alias and a
    passphrase hash.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    store = MagicMock(wraps=fake_redis)
    store.pipeline.return_value = _HalfAppliedPipeline(fake_redis)

    with pytest.raises(RuntimeError, match="went away"):
        create_cohort_record(store, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=RETENTION)

    assert fake_redis.keys(f"{COHORT_KEY_PREFIX}*") == [], "a cohort record survived with no way to reach it"
    assert fake_redis.get(alias_key(ALIAS)) is None


def test_a_created_cohort_carries_a_ttl(fake_redis) -> None:
    """The happy path still sets the retention the rollback above exists for.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    created = create_cohort_record(fake_redis, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=RETENTION)

    assert 0 < fake_redis.ttl(cohort_key(created["cohort_id"])) <= RETENTION


def test_the_alias_is_reusable_after_a_refused_creation(client) -> None:
    """End to end: a refused request does not take the name out of circulation.

    Args:
        client: TestClient fixture from conftest.
    """
    refused = client.post(
        "/create-cohort/",
        data={"alias": ALIAS, "passphrase": "p" * (MAX_PASSPHRASE_BYTES + 1)},
    )
    assert refused.status_code == 400
    assert str(MAX_PASSPHRASE_BYTES) in refused.json()["detail"]

    accepted = client.post("/create-cohort/", data={"alias": ALIAS, "passphrase": PASSPHRASE})

    assert accepted.status_code == 200, accepted.text
    assert accepted.json()["alias"] == ALIAS


def test_a_successful_creation_still_claims_the_alias(fake_redis) -> None:
    """The uniqueness guarantee is unchanged: a second cohort cannot take it.

    Without this, every test above would also pass against an implementation
    that had stopped claiming aliases altogether.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    created = create_cohort_record(fake_redis, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=RETENTION)

    assert fake_redis.get(alias_key(ALIAS)) == created["cohort_id"]
    with pytest.raises(ValueError, match="already in use"):
        create_cohort_record(fake_redis, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=RETENTION)


# ---------------------------------------------------------------------------
# A live cohort does not outlive its alias
# ---------------------------------------------------------------------------


def test_retention_extends_the_alias_alongside_the_cohort(fake_redis) -> None:
    """All three keys a cohort owns are pushed out together.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    created = create_cohort_record(fake_redis, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=60)
    key = cohort_key(created["cohort_id"])
    fake_redis.sadd(f"{key}:jobs", "a-member")
    fake_redis.expire(f"{key}:jobs", 60)

    extend_cohort_retention(fake_redis, key, RETENTION)

    assert fake_redis.ttl(key) > 60
    assert fake_redis.ttl(f"{key}:jobs") > 60
    assert fake_redis.ttl(alias_key(ALIAS)) > 60


def test_retention_of_a_cohort_without_an_alias_is_not_an_error(fake_redis) -> None:
    """An unlabelled cohort has no claim to extend and must not fail trying.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    created = create_cohort_record(fake_redis, alias=None, passphrase=PASSPHRASE, retention_seconds=60)
    key = cohort_key(created["cohort_id"])

    extend_cohort_retention(fake_redis, key, RETENTION)

    assert fake_redis.ttl(key) > 60
    assert fake_redis.keys(f"{ALIAS_KEY_PREFIX}*") == []


def test_a_finished_job_extends_its_cohorts_alias(fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    """The worker's own bookkeeping keeps the alias alive with the cohort.

    Extending the cohort but not its alias is what lets a cohort outlive its own
    name, so the worker is exercised here rather than only the helper it calls.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
        monkeypatch: Standard pytest fixture; restores every patch at teardown.
        tmp_path: Scratch directory standing in for the job input directory.
    """
    from app import tasks

    for attr in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, attr, fake_redis)
    monkeypatch.setattr(tasks.subprocess, "run", lambda *args, **kwargs: None)

    created = create_cohort_record(fake_redis, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=60)
    key = cohort_key(created["cohort_id"])
    fake_redis.sadd(f"{key}:jobs", "a-member")
    fake_redis.expire(f"{key}:jobs", 60)

    bam_path = tmp_path / "sample.bam"
    bam_path.write_bytes(b"alignment")
    (tmp_path / "sample.bam.bai").write_bytes(b"index")

    tasks.run_vntyper_job.push_request(id="task-1")
    try:
        tasks.run_vntyper_job.run(
            bam_path=str(bam_path),
            output_dir=str(tmp_path / "out"),
            thread=1,
            reference_assembly="hg38",
            fast_mode=False,
            keep_intermediates=False,
            archive_results=False,
            cohort_key=key,
        )
    finally:
        tasks.run_vntyper_job.pop_request()

    assert fake_redis.ttl(alias_key(ALIAS)) > 60
    assert fake_redis.ttl(key) > 60
    assert fake_redis.ttl(f"{key}:jobs") > 60


# ---------------------------------------------------------------------------
# Identifiers and aliases are bounded before they become Redis keys
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cohort_id",
    [
        pytest.param("3b045376-23b4-48e6-b3aa-12c0773d8137:jobs", id="reaches_the_member_set"),
        pytest.param("not-a-uuid", id="not_an_identifier"),
        pytest.param("*", id="glob"),
        pytest.param("x" * 5000, id="oversized"),
    ],
)
def test_a_cohort_id_that_is_not_one_this_service_minted_names_no_cohort(fake_redis, cohort_id: str) -> None:
    """An identifier is interpolated into a Redis key, so its shape is checked first.

    `<uuid>:jobs` is the cohort's own member Set. Reading it as a hash raises
    WRONGTYPE inside the client, which surfaces as a 500 -- the service reporting
    itself broken over a value the caller chose. Every value that is not one of
    the identifiers this service mints gets the one honest answer: it names no
    cohort.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
        cohort_id: The malformed identifier under test.
    """
    fake_redis.sadd("cohort:3b045376-23b4-48e6-b3aa-12c0773d8137:jobs", "job-a")

    with pytest.raises(LookupError, match="not found"):
        resolve_cohort(fake_redis, cohort_id=cohort_id, alias=None, passphrase=PASSPHRASE)


@pytest.mark.parametrize(
    "alias",
    [
        pytest.param("a" * (MAX_ALIAS_CHARS + 1), id="too_long"),
        pytest.param("family\na", id="newline"),
        pytest.param("family\x00a", id="null_byte"),
    ],
)
def test_an_alias_outside_the_documented_policy_is_refused(fake_redis, alias: str) -> None:
    """An alias becomes a Redis key too, so it is bounded on the way in.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
        alias: The alias under test.
    """
    with pytest.raises(ValueError, match="alias"):
        create_cohort_record(fake_redis, alias=alias, passphrase=PASSPHRASE, retention_seconds=RETENTION)

    assert fake_redis.keys(f"{ALIAS_KEY_PREFIX}*") == []
    assert fake_redis.keys(f"{COHORT_KEY_PREFIX}*") == []


def test_an_alias_at_the_limit_is_still_accepted(fake_redis) -> None:
    """The bound is a ceiling, not an off-by-one that rejects the longest legal alias.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    alias = "a" * MAX_ALIAS_CHARS

    created = create_cohort_record(fake_redis, alias=alias, passphrase=PASSPHRASE, retention_seconds=RETENTION)

    assert created["alias"] == alias
    assert fake_redis.get(alias_key(alias)) == created["cohort_id"]


def test_a_malformed_cohort_id_is_reported_as_not_found_rather_than_as_a_server_error(client) -> None:
    """End to end: the caller gets a constrained answer, not a 500.

    Args:
        client: TestClient fixture from conftest.
    """
    response = client.get(
        "/cohort-status/",
        params={"cohort_id": "3b045376-23b4-48e6-b3aa-12c0773d8137:jobs", "passphrase": PASSPHRASE},
    )

    assert response.status_code == 404, response.text
