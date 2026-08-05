"""Authorization for the cohort routes of the web service.

A cohort groups jobs whose results can be read, downloaded and re-analysed
together. The service has no accounts, no sessions and no tokens, so the pair
(cohort identifier, passphrase) is the only credential that exists. Everything
below pins the two properties that follow from that:

1. The server-generated ``cohort_id`` is the only thing that may select a
   cohort. The human-chosen ``alias`` is a label, verified against the resolved
   cohort when supplied, never used to find one.
2. A passphrase is required on reads as well as on writes, including for records
   that carry no stored hash.

The unit-level cases drive ``app.cohorts.resolve_cohort`` against ``fake_redis``
directly, because authorization is decided there and needs no app. The endpoint
cases drive the same rules through ``client`` to prove the guard is actually
wired into every route and that the failures map onto the right status codes.

The passphrase hashing every rule below depends on is pinned separately, in
``test_passphrase_hashing.py``.

``docker`` is put on ``sys.path`` by ``tests/unit/web/conftest.py``, which pytest
imports before this module, so ``app.cohorts`` is importable here.
"""

from pathlib import Path
from unittest.mock import patch

import pytest

pytestmark = pytest.mark.unit

from app.cohorts import (  # noqa: E402
    ALIAS_KEY_PREFIX,
    COHORT_KEY_PREFIX,
    cohort_job_ids,
    create_cohort_record,
    resolve_cohort,
)

RETENTION_SECONDS = 14 * 86400

# Any passphrase works; bcrypt is deliberately slow, so the tests reuse one.
GOOD_PASSPHRASE = "correct-horse-battery-staple"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _seed_unprotected_cohort(store, cohort_id: str, alias: str, job_ids=()) -> None:
    """Write a cohort record that carries no passphrase hash.

    This is exactly the shape the service used to store when a creator supplied
    an alias but no passphrase, so it stands in both for such a record already
    in Redis and for any future write that forgets the hash.

    Args:
        store: The Redis client backing the cohort database.
        cohort_id: Identifier to store the record under.
        alias: The cohort's human-chosen label.
        job_ids: Job identifiers to place in the cohort's member set.
    """
    store.hset(
        f"{COHORT_KEY_PREFIX}{cohort_id}",
        mapping={
            "alias": alias,
            "hashed_passphrase": "",
            "created_at": "2026-01-01T00:00:00+00:00",
        },
    )
    for job_id in job_ids:
        store.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", job_id)


def _create_cohort_via_api(client, alias: str, passphrase: str = GOOD_PASSPHRASE) -> str:
    """Create a passphrase-protected cohort through the public endpoint.

    Args:
        client: TestClient fixture from conftest.
        alias: The alias to give the cohort.
        passphrase: The passphrase protecting it.

    Returns:
        str: The identifier of the created cohort.
    """
    response = client.post("/create-cohort/", data={"alias": alias, "passphrase": passphrase})
    assert response.status_code == 200, response.text
    return response.json()["cohort_id"]


def _write_result_zip(tmp_path: Path, job_id: str) -> Path:
    """Place a stand-in result archive where the service looks for one.

    Args:
        tmp_path: The directory the `web_app` fixture configured as the output tree.
        job_id: The job the archive belongs to.

    Returns:
        Path: The archive that was written.
    """
    zip_path = tmp_path / "output" / f"{job_id}.zip"
    zip_path.write_bytes(b"PK\x05\x06" + b"\x00" * 18)
    return zip_path


# ---------------------------------------------------------------------------
# Unit level: resolve_cohort against fake_redis
# ---------------------------------------------------------------------------


def test_create_cohort_record_rejects_a_missing_passphrase(fake_redis) -> None:
    """A cohort cannot be created without the credential that will guard it.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    for passphrase in (None, "", "   "):
        with pytest.raises(ValueError):
            create_cohort_record(
                fake_redis,
                alias="study",
                passphrase=passphrase,
                retention_seconds=RETENTION_SECONDS,
            )

    assert fake_redis.keys("*") == []


def test_create_cohort_record_stores_a_hash_never_the_passphrase(fake_redis) -> None:
    """The passphrase is stored only as a bcrypt hash.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    record = create_cohort_record(
        fake_redis,
        alias="study",
        passphrase=GOOD_PASSPHRASE,
        retention_seconds=RETENTION_SECONDS,
    )

    stored = fake_redis.hgetall(f"{COHORT_KEY_PREFIX}{record['cohort_id']}")
    assert stored["hashed_passphrase"].startswith("$2")
    assert GOOD_PASSPHRASE not in stored["hashed_passphrase"]
    assert fake_redis.ttl(f"{COHORT_KEY_PREFIX}{record['cohort_id']}") > 0


def test_create_cohort_record_rejects_a_duplicate_alias(fake_redis) -> None:
    """Two cohorts cannot hold the same alias.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    first = create_cohort_record(
        fake_redis,
        alias="shared-study",
        passphrase=GOOD_PASSPHRASE,
        retention_seconds=RETENTION_SECONDS,
    )

    with pytest.raises(ValueError):
        create_cohort_record(
            fake_redis,
            alias="shared-study",
            passphrase="a-different-passphrase",
            retention_seconds=RETENTION_SECONDS,
        )

    # The loser of the race leaves nothing behind, and the winner still holds
    # the alias.
    assert fake_redis.get(f"{ALIAS_KEY_PREFIX}shared-study") == first["cohort_id"]
    assert len(fake_redis.keys(f"{COHORT_KEY_PREFIX}*")) == 1


def test_alias_is_claimed_with_a_single_atomic_command(fake_redis) -> None:
    """The alias claim is one conditional write, not a read followed by a write.

    A `SET ... NX` either creates the key or reports that someone else holds it.
    A separate "does it exist?" check followed by a write would let two
    simultaneous creations both pass the check and both believe they won.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    with patch.object(fake_redis, "set", wraps=fake_redis.set) as spy:
        create_cohort_record(
            fake_redis,
            alias="study",
            passphrase=GOOD_PASSPHRASE,
            retention_seconds=RETENTION_SECONDS,
        )

    assert spy.call_count == 1
    assert spy.call_args.kwargs.get("nx") is True

    # A claim held by nobody's cohort still blocks: the claim itself is the
    # authority, so it cannot be worked around by removing a cohort record.
    with pytest.raises(ValueError):
        create_cohort_record(
            fake_redis,
            alias="study",
            passphrase=GOOD_PASSPHRASE,
            retention_seconds=RETENTION_SECONDS,
        )


def test_resolve_cohort_requires_the_identifier(fake_redis) -> None:
    """An alias on its own never selects a cohort.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    _seed_unprotected_cohort(fake_redis, "cohort-1", "shared-study")

    with pytest.raises(ValueError):
        resolve_cohort(fake_redis, cohort_id=None, alias="shared-study", passphrase=None)

    with pytest.raises(ValueError):
        resolve_cohort(fake_redis, cohort_id=None, alias="shared-study", passphrase=GOOD_PASSPHRASE)


def test_resolve_cohort_requires_an_identifier_or_an_alias(fake_redis) -> None:
    """A request naming no cohort at all is malformed.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    with pytest.raises(ValueError):
        resolve_cohort(fake_redis, cohort_id=None, alias=None, passphrase=GOOD_PASSPHRASE)


def test_resolve_cohort_rejects_a_missing_passphrase(fake_redis) -> None:
    """Knowing the identifier is not on its own enough to read a cohort.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    record = create_cohort_record(
        fake_redis, alias="study", passphrase=GOOD_PASSPHRASE, retention_seconds=RETENTION_SECONDS
    )

    with pytest.raises(PermissionError):
        resolve_cohort(fake_redis, cohort_id=record["cohort_id"], alias=None, passphrase=None)

    with pytest.raises(PermissionError):
        resolve_cohort(fake_redis, cohort_id=record["cohort_id"], alias=None, passphrase="")


def test_resolve_cohort_rejects_a_wrong_passphrase(fake_redis) -> None:
    """A passphrase that does not verify is refused.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    record = create_cohort_record(
        fake_redis, alias="study", passphrase=GOOD_PASSPHRASE, retention_seconds=RETENTION_SECONDS
    )

    with pytest.raises(PermissionError):
        resolve_cohort(fake_redis, cohort_id=record["cohort_id"], alias=None, passphrase="not-it")


def test_resolve_cohort_reports_an_unknown_identifier_separately(fake_redis) -> None:
    """An identifier with no record behind it is a lookup failure, not a refusal.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    with pytest.raises(LookupError):
        resolve_cohort(fake_redis, cohort_id="no-such-cohort", alias=None, passphrase=GOOD_PASSPHRASE)


def test_resolve_cohort_denies_a_record_with_no_stored_hash(fake_redis) -> None:
    """A record carrying no hash is closed, not open.

    Treating an absent hash as "no passphrase needed" would turn every record
    written before the passphrase became mandatory into an unguarded one, which
    is the opposite of what this change is for.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    _seed_unprotected_cohort(fake_redis, "cohort-1", "shared-study")

    with pytest.raises(PermissionError):
        resolve_cohort(fake_redis, cohort_id="cohort-1", alias=None, passphrase=GOOD_PASSPHRASE)

    with pytest.raises(PermissionError):
        resolve_cohort(fake_redis, cohort_id="cohort-1", alias="shared-study", passphrase="anything")


def test_resolve_cohort_rejects_a_mismatched_alias(fake_redis) -> None:
    """A supplied alias must agree with the resolved cohort.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    record = create_cohort_record(
        fake_redis, alias="study", passphrase=GOOD_PASSPHRASE, retention_seconds=RETENTION_SECONDS
    )

    with pytest.raises(ValueError):
        resolve_cohort(fake_redis, cohort_id=record["cohort_id"], alias="some-other-study", passphrase=GOOD_PASSPHRASE)


def test_resolve_cohort_accepts_the_right_credentials(fake_redis) -> None:
    """The legitimate path keeps working, with and without the alias.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    record = create_cohort_record(
        fake_redis, alias="study", passphrase=GOOD_PASSPHRASE, retention_seconds=RETENTION_SECONDS
    )
    cohort_id = record["cohort_id"]
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")

    resolved = resolve_cohort(fake_redis, cohort_id=cohort_id, alias=None, passphrase=GOOD_PASSPHRASE)
    assert resolved["cohort_id"] == cohort_id
    assert resolved["alias"] == "study"
    assert cohort_job_ids(fake_redis, resolved["cohort_key"]) == ["job-a"]

    with_alias = resolve_cohort(fake_redis, cohort_id=cohort_id, alias="study", passphrase=GOOD_PASSPHRASE)
    assert with_alias["cohort_id"] == cohort_id


def test_resolve_cohort_never_returns_the_stored_hash(fake_redis) -> None:
    """The credential material stays inside the module that checks it.

    Args:
        fake_redis: In-process Redis stand-in from conftest.
    """
    record = create_cohort_record(
        fake_redis, alias="study", passphrase=GOOD_PASSPHRASE, retention_seconds=RETENTION_SECONDS
    )

    resolved = resolve_cohort(fake_redis, cohort_id=record["cohort_id"], alias=None, passphrase=GOOD_PASSPHRASE)

    assert "hashed_passphrase" not in resolved
    assert GOOD_PASSPHRASE not in repr(resolved)


# ---------------------------------------------------------------------------
# Endpoint level: /create-cohort/
# ---------------------------------------------------------------------------


def test_create_cohort_endpoint_requires_a_passphrase_field(client) -> None:
    """Creating a cohort without the passphrase field is refused by the schema.

    Args:
        client: TestClient fixture from conftest.
    """
    response = client.post("/create-cohort/", data={"alias": "study"})

    assert response.status_code == 422


def test_create_cohort_endpoint_rejects_a_blank_passphrase(client, fake_redis) -> None:
    """A present but empty passphrase is refused, and stores nothing.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    response = client.post("/create-cohort/", data={"alias": "study", "passphrase": "  "})

    assert response.status_code == 400
    assert fake_redis.keys(f"{COHORT_KEY_PREFIX}*") == []


def test_create_cohort_endpoint_rejects_a_duplicate_alias(client, fake_redis) -> None:
    """The second cohort to ask for an alias is refused it.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    _create_cohort_via_api(client, "shared-study")

    response = client.post("/create-cohort/", data={"alias": "shared-study", "passphrase": "another-passphrase"})

    assert response.status_code == 400
    assert len(fake_redis.keys(f"{COHORT_KEY_PREFIX}*")) == 1


def test_create_cohort_endpoint_still_creates_a_cohort(client, fake_redis) -> None:
    """The legitimate creation path is unchanged apart from the new requirement.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    response = client.post("/create-cohort/", data={"alias": "study", "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200
    body = response.json()
    assert body["alias"] == "study"
    assert fake_redis.hget(f"{COHORT_KEY_PREFIX}{body['cohort_id']}", "hashed_passphrase")


# ---------------------------------------------------------------------------
# Endpoint level: /cohort-status/
# ---------------------------------------------------------------------------


def test_cohort_status_denies_alias_only_access(client, fake_redis) -> None:
    """A guessed alias does not reveal a cohort's membership.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    _seed_unprotected_cohort(fake_redis, "cohort-1", "shared-study", job_ids=["job-a"])

    response = client.get("/cohort-status/", params={"alias": "shared-study"})

    assert response.status_code == 400
    assert "job-a" not in response.text


def test_cohort_status_denies_a_record_with_no_stored_hash(client, fake_redis) -> None:
    """Holding the identifier of an unguarded record is still not enough.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    _seed_unprotected_cohort(fake_redis, "cohort-1", "shared-study", job_ids=["job-a"])

    response = client.get("/cohort-status/", params={"cohort_id": "cohort-1"})

    assert response.status_code == 401
    assert "job-a" not in response.text


def test_cohort_status_requires_the_passphrase(client) -> None:
    """A protected cohort is not readable from its identifier alone.

    Args:
        client: TestClient fixture from conftest.
    """
    cohort_id = _create_cohort_via_api(client, "study")

    response = client.get("/cohort-status/", params={"cohort_id": cohort_id})

    assert response.status_code == 401


def test_cohort_status_rejects_a_wrong_passphrase(client) -> None:
    """A wrong passphrase is refused with the same code as a missing one.

    Args:
        client: TestClient fixture from conftest.
    """
    cohort_id = _create_cohort_via_api(client, "study")

    response = client.get("/cohort-status/", params={"cohort_id": cohort_id, "passphrase": "not-it"})

    assert response.status_code == 401


def test_cohort_status_reports_an_unknown_identifier_as_not_found(client) -> None:
    """An identifier with no record behind it is a 404.

    Args:
        client: TestClient fixture from conftest.
    """
    response = client.get("/cohort-status/", params={"cohort_id": "no-such-cohort", "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 404


def test_cohort_status_without_an_identifier_or_alias_is_a_bad_request(client) -> None:
    """Naming no cohort at all stays a 400.

    Args:
        client: TestClient fixture from conftest.
    """
    response = client.get("/cohort-status/", params={"passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 400


def test_cohort_status_works_with_the_right_passphrase(client, fake_redis) -> None:
    """The legitimate read path still returns the cohort's jobs.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    cohort_id = _create_cohort_via_api(client, "study")
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")

    response = client.get("/cohort-status/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200
    body = response.json()
    assert body["cohort_id"] == cohort_id
    assert body["alias"] == "study"
    assert [job["job_id"] for job in body["jobs"]] == ["job-a"]


# ---------------------------------------------------------------------------
# Endpoint level: /cohort-download/
# ---------------------------------------------------------------------------


def test_cohort_download_denies_alias_only_access(client, fake_redis, tmp_path: Path) -> None:
    """A guessed alias does not hand over every member's results.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the fixture configured as the output tree.
    """
    _seed_unprotected_cohort(fake_redis, "cohort-1", "shared-study", job_ids=["job-a"])
    _write_result_zip(tmp_path, "job-a")

    response = client.get("/cohort-download/", params={"alias": "shared-study"})

    assert response.status_code == 400


def test_cohort_download_requires_the_passphrase(client, fake_redis, tmp_path: Path) -> None:
    """The bulk download is not reachable from the identifier alone.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the fixture configured as the output tree.
    """
    cohort_id = _create_cohort_via_api(client, "study")
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")
    _write_result_zip(tmp_path, "job-a")

    response = client.get("/cohort-download/", params={"cohort_id": cohort_id})

    assert response.status_code == 401


def test_cohort_download_works_with_the_right_passphrase(client, fake_redis, tmp_path: Path) -> None:
    """The legitimate download path still serves the archive.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the fixture configured as the output tree.
    """
    cohort_id = _create_cohort_via_api(client, "study")
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")
    _write_result_zip(tmp_path, "job-a")

    response = client.get("/cohort-download/", params={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200
    assert response.headers["content-type"] == "application/zip"


# ---------------------------------------------------------------------------
# Endpoint level: /cohort-analysis/
# ---------------------------------------------------------------------------


def test_cohort_analysis_denies_alias_only_access(client, web_app, fake_redis, tmp_path: Path) -> None:
    """A guessed alias cannot start an analysis over someone else's cohort.

    Args:
        client: TestClient fixture from conftest.
        web_app: Patched `app.main` module, sharing this test's `tmp_path`.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the fixture configured as the output tree.
    """
    _seed_unprotected_cohort(fake_redis, "cohort-1", "shared-study", job_ids=["job-a"])
    _write_result_zip(tmp_path, "job-a")

    response = client.post("/cohort-analysis/", data={"alias": "shared-study"})

    assert response.status_code == 400
    web_app.run_cohort_analysis_job.delay.assert_not_called()


def test_cohort_analysis_requires_the_passphrase(client, web_app, fake_redis, tmp_path: Path) -> None:
    """An analysis is a write, and is refused without the credential.

    Args:
        client: TestClient fixture from conftest.
        web_app: Patched `app.main` module, sharing this test's `tmp_path`.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the fixture configured as the output tree.
    """
    cohort_id = _create_cohort_via_api(client, "study")
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")
    _write_result_zip(tmp_path, "job-a")

    response = client.post("/cohort-analysis/", data={"cohort_id": cohort_id})

    assert response.status_code == 401
    web_app.run_cohort_analysis_job.delay.assert_not_called()


def test_cohort_analysis_works_with_the_right_passphrase(client, web_app, fake_redis, tmp_path: Path) -> None:
    """The legitimate analysis path still enqueues the job.

    Args:
        client: TestClient fixture from conftest.
        web_app: Patched `app.main` module, sharing this test's `tmp_path`.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the fixture configured as the output tree.
    """
    cohort_id = _create_cohort_via_api(client, "study")
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", "job-a")
    _write_result_zip(tmp_path, "job-a")

    response = client.post("/cohort-analysis/", data={"cohort_id": cohort_id, "passphrase": GOOD_PASSPHRASE})

    assert response.status_code == 200
    web_app.run_cohort_analysis_job.delay.assert_called_once()


# ---------------------------------------------------------------------------
# Endpoint level: /run-job/, which also resolves a cohort
# ---------------------------------------------------------------------------


def test_run_job_denies_alias_only_cohort_assignment(client, fake_redis) -> None:
    """A job cannot be pushed into a cohort named only by its alias.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    _seed_unprotected_cohort(fake_redis, "cohort-1", "shared-study")

    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"x", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19", "alias": "shared-study"},
    )

    assert response.status_code == 400
    assert fake_redis.smembers(f"{COHORT_KEY_PREFIX}cohort-1:jobs") == set()


def test_run_job_requires_the_passphrase_for_cohort_assignment(client, fake_redis) -> None:
    """Adding a job to a cohort needs the same credential as reading one.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    cohort_id = _create_cohort_via_api(client, "study")

    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"x", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19", "cohort_id": cohort_id},
    )

    assert response.status_code == 401
    assert fake_redis.smembers(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs") == set()


def test_run_job_still_joins_a_cohort_with_the_right_passphrase(client, fake_redis) -> None:
    """The legitimate submission path still associates the job.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
    """
    cohort_id = _create_cohort_via_api(client, "study")

    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"x", "application/octet-stream")},
        data={
            "thread": "1",
            "reference_assembly": "hg19",
            "cohort_id": cohort_id,
            "passphrase": GOOD_PASSPHRASE,
        },
    )

    assert response.status_code == 200
    job_id = response.json()["job_id"]
    assert fake_redis.smembers(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs") == {job_id}
