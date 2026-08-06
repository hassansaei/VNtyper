"""The passphrase hashing the web service's cohort authorization rests on.

A cohort passphrase is the only credential the service has. `hash_passphrase`
and `verify_passphrase` in `docker/app/utils.py` are what turn that credential
into something storable and what decide, later, whether a caller supplied it.
Every cohort rule in `test_cohort_authz.py` is decided by these two functions,
so they are pinned here on their own: if hashing does not work, nothing built on
top of it does either, and a failure needs to point at the primitive rather than
at twenty authorization tests at once.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.utils` is importable here.
"""

import pytest

pytestmark = pytest.mark.unit

from app.utils import MAX_PASSPHRASE_BYTES, hash_passphrase, verify_passphrase  # noqa: E402

PASSPHRASE = "correct-horse-battery-staple"


def test_a_passphrase_can_be_hashed_and_verified() -> None:
    """Hashing a passphrase and checking it back is a working operation."""
    hashed = hash_passphrase(PASSPHRASE)

    assert hashed.startswith("$2")
    assert PASSPHRASE not in hashed
    assert verify_passphrase(PASSPHRASE, hashed) is True
    assert verify_passphrase("not-it", hashed) is False


def test_hashing_the_same_passphrase_twice_gives_different_hashes() -> None:
    """Each hash carries its own salt."""
    assert hash_passphrase(PASSPHRASE) != hash_passphrase(PASSPHRASE)


def test_an_over_long_passphrase_is_refused_rather_than_truncated() -> None:
    """The algorithm's input limit is stated, not silently applied.

    bcrypt reads a bounded number of bytes. Accepting a longer passphrase and
    hashing a prefix of it would mean a passphrase that verifies against input
    the user never chose.
    """
    over_long = "p" * (MAX_PASSPHRASE_BYTES + 1)

    with pytest.raises(ValueError):
        hash_passphrase(over_long)

    at_limit = "p" * MAX_PASSPHRASE_BYTES
    hashed = hash_passphrase(at_limit)
    assert verify_passphrase(at_limit, hashed) is True
    assert verify_passphrase(over_long, hashed) is False


def test_verifying_against_an_unusable_stored_hash_returns_false() -> None:
    """A stored value that is not a hash fails the check instead of erroring.

    Raising here would turn an unopenable cohort into a 500 rather than a
    refusal, and the two must not be confusable.
    """
    for stored in ("", "not-a-bcrypt-hash", "$2b$12$too-short"):
        assert verify_passphrase(PASSPHRASE, stored) is False
