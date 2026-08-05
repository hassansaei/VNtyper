"""Cohort identity and authorization for the web service.

A cohort groups jobs whose results can be read, downloaded and re-analysed
together. The service has no accounts, no sessions and no tokens, so the pair
(cohort identifier, passphrase) is the only credential it has to work with.
Every cohort read and every cohort write is checked against that pair here,
in one place, rather than being restated at each endpoint.

Two rules follow from having a single credential, and both are enforced below:

1. The server-generated ``cohort_id`` is the only thing that may select a
   cohort. The ``alias`` is chosen by a person, so it is a label rather than an
   identifier: it is verified against the resolved cohort when one is supplied,
   and it is claimed exclusively at creation time, but it never selects a
   cohort on its own.
2. A passphrase is mandatory -- at creation, and on every subsequent read and
   write. A record that carries no stored hash cannot be opened at all; see
   ``resolve_cohort`` for why that direction is the safe one.

Failures are signalled with builtin exceptions used as a small vocabulary:
``ValueError`` for a malformed request, ``PermissionError`` for a refused one
and ``LookupError`` for a cohort that does not exist. The API layer maps those
onto 400, 401 and 404. Keeping the vocabulary builtin follows the project's
"no custom exception classes" convention, and keeping FastAPI out of this module
lets the rules be tested against a plain Redis client with no app in the way.

The Redis client is passed in rather than imported so that the caller decides
which database these keys live in, and so the tests can supply an in-process
stand-in.
"""

import logging
from datetime import datetime, timezone
from typing import Any
from uuid import uuid4

from .utils import hash_passphrase, verify_passphrase

logger = logging.getLogger(__name__)

# Cohort metadata lives in a hash per cohort.
COHORT_KEY_PREFIX = "cohort:"

# Each claimed alias gets one key, whose existence *is* the claim.
ALIAS_KEY_PREFIX = "cohort-alias:"


def cohort_key(cohort_id: str) -> str:
    """Build the Redis key holding a cohort's metadata.

    Args:
        cohort_id: The cohort's identifier.

    Returns:
        str: The key the cohort hash is stored under.
    """
    return f"{COHORT_KEY_PREFIX}{cohort_id}"


def alias_key(alias: str) -> str:
    """Build the Redis key that records ownership of an alias.

    Args:
        alias: The alias being claimed or looked up.

    Returns:
        str: The key the claim is stored under.
    """
    return f"{ALIAS_KEY_PREFIX}{alias}"


def _clean(value: str | None) -> str | None:
    """Reduce a form field to a meaningful value or nothing at all.

    An omitted field and a field holding only whitespace mean the same thing to
    every caller here, so they are collapsed to ``None`` before any check runs.

    Args:
        value: The raw field value.

    Returns:
        str | None: The trimmed value, or None if it carried no content.
    """
    return (value or "").strip() or None


def create_cohort_record(
    store: Any,
    *,
    alias: str | None,
    passphrase: str | None,
    retention_seconds: int,
) -> dict[str, str | None]:
    """Create a cohort, claiming its alias and storing its passphrase hash.

    The alias claim is made with a single ``SET ... NX EX`` command, which is
    atomic in Redis: it either creates the key or reports that another cohort
    already holds it. A separate "is this alias free?" read followed by a write
    would let two simultaneous creations both pass the read and both believe
    they had won the alias. The claim is made before the cohort hash is written,
    so a refused alias leaves nothing behind.

    Args:
        store: Redis client for the cohort database.
        alias: Optional human-chosen label for the cohort.
        passphrase: The credential that will guard the cohort. Required.
        retention_seconds: Lifetime of the cohort record and its alias claim.

    Returns:
        dict[str, str | None]: The new ``cohort_id`` and the ``alias`` it holds.

    Raises:
        ValueError: If no passphrase was supplied, or the alias is already taken.
    """
    alias = _clean(alias)
    # Normalised to `str` here, exactly as `resolve_cohort` does, so that the value
    # handed to `hash_passphrase` below is the one the caller sent. Only the
    # *emptiness* test is made on the trimmed form: `resolve_cohort` verifies the
    # raw passphrase, so trimming before hashing would lock out any passphrase with
    # leading or trailing whitespace.
    passphrase = passphrase or ""
    if not _clean(passphrase):
        msg = "A passphrase is required to create a cohort"
        logger.error(msg)
        raise ValueError(msg)

    cohort_id = str(uuid4())

    if alias is not None and not store.set(alias_key(alias), cohort_id, nx=True, ex=retention_seconds):
        msg = "Cohort alias is already in use"
        logger.error(msg)
        raise ValueError(msg)

    key = cohort_key(cohort_id)
    store.hset(
        key,
        mapping={
            "alias": alias or "",
            "hashed_passphrase": hash_passphrase(passphrase),
            "created_at": datetime.now(tz=timezone.utc).isoformat(),
        },
    )
    store.expire(key, retention_seconds)

    logger.info(f"Created cohort {cohort_id}")
    return {"cohort_id": cohort_id, "alias": alias}


def resolve_cohort(
    store: Any,
    *,
    cohort_id: str | None,
    alias: str | None,
    passphrase: str | None,
) -> dict[str, str]:
    """Identify a cohort and authorize the caller for it.

    Checks run in a fixed order so that a request carrying no credential is
    refused before it can reach the store at all:

    1. An identifier must be supplied. An alias alone is not one.
    2. A passphrase must be supplied.
    3. The cohort must exist.
    4. Any supplied alias must agree with the cohort that was found.
    5. The cohort must carry a stored hash, and the passphrase must verify
       against it.

    Step 5 refuses a record whose hash is absent. Records written before a
    passphrase was mandatory have an empty hash, and reading "no hash" as "no
    passphrase needed" would leave exactly those records open, which is the
    opposite of what the check is for. Such a record cannot be opened; cohort
    records carry a retention TTL, so the affected set empties itself.

    The comparison in step 5 is delegated to `utils.verify_passphrase`, which
    compares in constant time. No comparison is hand-written here.

    Args:
        store: Redis client for the cohort database.
        cohort_id: The cohort's identifier. Required.
        alias: Optional label, cross-checked against the resolved cohort.
        passphrase: The cohort's passphrase. Required.

    Returns:
        dict[str, str]: The resolved ``cohort_id``, its ``cohort_key`` and its
            stored ``alias``. The stored hash is deliberately not included.

    Raises:
        ValueError: If no identifier was supplied, or the alias disagrees with
            the cohort.
        PermissionError: If no passphrase was supplied, if the cohort carries no
            stored hash, or if the passphrase does not verify.
        LookupError: If no cohort exists for the identifier.
    """
    cohort_id = _clean(cohort_id)
    alias = _clean(alias)
    passphrase = passphrase or ""

    if cohort_id is None:
        msg = "Cohort identifier required"
        logger.error(msg)
        raise ValueError(msg)

    if not passphrase:
        msg = "Passphrase required for this cohort"
        logger.error(msg)
        raise PermissionError(msg)

    key = cohort_key(cohort_id)
    record = store.hgetall(key)
    if not record:
        msg = "Cohort ID not found"
        logger.error(msg)
        raise LookupError(msg)

    stored_alias = record.get("alias", "")
    if alias is not None and stored_alias != alias:
        msg = "Provided alias does not match the cohort's alias"
        logger.error(msg)
        raise ValueError(msg)

    stored_hash = record.get("hashed_passphrase", "")
    if not stored_hash:
        msg = "Cohort cannot be opened"
        logger.error(f"{msg}: cohort {cohort_id} carries no stored passphrase hash")
        raise PermissionError(msg)

    if not verify_passphrase(passphrase, stored_hash):
        msg = "Incorrect passphrase"
        logger.error(msg)
        raise PermissionError(msg)

    return {"cohort_id": cohort_id, "cohort_key": key, "alias": stored_alias}


def cohort_job_ids(store: Any, key: str) -> list[str]:
    """List the jobs belonging to a cohort.

    Args:
        store: Redis client for the cohort database.
        key: The cohort's metadata key, as returned by `resolve_cohort`.

    Returns:
        list[str]: The member job identifiers, empty if the cohort has none.
    """
    return list(store.smembers(f"{key}:jobs") or [])
