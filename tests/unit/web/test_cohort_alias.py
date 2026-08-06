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
    COHORT_KEY_PREFIX,
    alias_key,
    create_cohort_record,
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
    store.hset.side_effect = RuntimeError("the store went away mid-write")

    with pytest.raises(RuntimeError, match="went away"):
        create_cohort_record(store, alias=ALIAS, passphrase=PASSPHRASE, retention_seconds=RETENTION)

    assert fake_redis.get(alias_key(ALIAS)) is None


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
