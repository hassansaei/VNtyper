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

# The two cohort read routes are GETs, so their credential used to have nowhere
# to travel but the query string -- part of the request line, and so written to
# every server and proxy access log on the path, kept in browser history and
# sent on in ``Referer`` headers. This header carries it instead. The query
# parameter still works and is documented as deprecated rather than removed:
# the web UI and the ``vntyper online`` CLI subcommand both use it. The three
# names below live here, beside ``preferred_passphrase``, so the rule and the
# words that describe it to a client author cannot drift apart.
COHORT_PASSPHRASE_HEADER = "X-Cohort-Passphrase"

PASSPHRASE_HEADER_DESCRIPTION = (
    "Passphrase protecting the cohort. Preferred over the `passphrase` query "
    "parameter, and used instead of it when both are supplied."
)

PASSPHRASE_QUERY_DESCRIPTION = (
    "Passphrase protecting the cohort. Deprecated: a query string is recorded "
    f"in access logs and browser history, so send the `{COHORT_PASSPHRASE_HEADER}` "
    "request header instead. Still accepted, and still required if the header is "
    "not sent."
)


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


def preferred_passphrase(header_value: str | None, query_value: str | None) -> str | None:
    """Choose which of the two ways in carried the caller's passphrase.

    The cohort read routes are GETs, so their passphrase used to travel in the
    query string -- which is part of the request line, and therefore ends up in
    server and proxy access logs, in browser history and in ``Referer`` headers.
    A request header carries it instead. The query parameter still works, and is
    marked deprecated rather than removed, because existing clients use it.

    The header wins when both are present. "Whichever verifies" would turn two
    credentials into two attempts at one, and a caller sending both is better
    told plainly which was used.

    Neither value is trimmed. Only the *emptiness* test uses the trimmed form, so
    that a header a proxy adds unconditionally does not blank out a credential
    the caller did send, while a passphrase with leading or trailing whitespace
    still verifies against the hash made from it.

    Args:
        header_value: The passphrase from the request header, if any.
        query_value: The passphrase from the query string, if any.

    Returns:
        str | None: The passphrase to authorize with, or None if neither way in
            carried one.
    """
    if _clean(header_value) is not None:
        return header_value
    return query_value


def _discard_partial_cohort(store: Any, key: str, claimed_alias: str | None) -> None:
    """Undo what a failed creation may have written.

    Called while an exception is propagating, so neither delete may replace it
    with one of its own: what the caller needs to see is why the creation failed,
    not why the tidying did.

    Args:
        store: Redis client for the cohort database.
        key: The cohort's metadata key.
        claimed_alias: The alias claimed for this cohort, if any.
    """
    try:
        store.delete(key)
    except Exception as exc:  # noqa: BLE001 - reported, never raised over the original
        logger.error(f"Could not remove the partial cohort record {key}: {exc}")
    if claimed_alias is None:
        return
    try:
        store.delete(alias_key(claimed_alias))
        logger.error(f"Released the claim on alias {claimed_alias}: its cohort record was not written")
    except Exception as exc:  # noqa: BLE001 - reported, never raised over the original
        logger.error(f"Could not release the claim on alias {claimed_alias}: {exc}")


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
    they had won the alias.

    A claim lasts a whole retention period, so it is only ever taken on behalf
    of a cohort that is going to exist. Everything the request can be refused
    for -- an absent passphrase, one the hashing scheme cannot store -- is
    decided *before* the claim, and the claim is released again if writing the
    record after it fails. Between those two, an alias is out of circulation
    only while it names something.

    Args:
        store: Redis client for the cohort database.
        alias: Optional human-chosen label for the cohort.
        passphrase: The credential that will guard the cohort. Required, and
            bounded by ``utils.MAX_PASSPHRASE_BYTES`` once encoded.
        retention_seconds: Lifetime of the cohort record and its alias claim.

    Returns:
        dict[str, str | None]: The new ``cohort_id`` and the ``alias`` it holds.

    Raises:
        ValueError: If no passphrase was supplied, if it is longer than the
            hashing scheme accepts, or if the alias is already taken.
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

    # Hashing is what can still refuse this passphrase -- it states and applies
    # the scheme's own input bound -- so it happens here, before anything is
    # written, rather than after the alias has been taken.
    hashed_passphrase = hash_passphrase(passphrase)

    cohort_id = str(uuid4())

    claimed_alias = None
    if alias is not None:
        if not store.set(alias_key(alias), cohort_id, nx=True, ex=retention_seconds):
            msg = "Cohort alias is already in use"
            logger.error(msg)
            raise ValueError(msg)
        claimed_alias = alias

    key = cohort_key(cohort_id)
    try:
        # The record and its expiry go in one round trip. Written separately, a
        # failure between them leaves a cohort hash with no TTL, holding an alias
        # and a passphrase hash: the identifier is minted in here and only reaches
        # the caller in the response, so a creation that raises never reveals it.
        # That record would be unreachable and permanent at the same time.
        pipeline = store.pipeline()
        pipeline.hset(
            key,
            mapping={
                "alias": alias or "",
                "hashed_passphrase": hashed_passphrase,
                "created_at": datetime.now(tz=timezone.utc).isoformat(),
            },
        )
        pipeline.expire(key, retention_seconds)
        pipeline.execute()
    except Exception:
        # Belt and braces behind the transaction: a client that applied the HSET
        # and then lost the connection before EXEC is indistinguishable from here,
        # so the key is deleted rather than assumed absent.
        _discard_partial_cohort(store, key, claimed_alias)
        raise

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


def extend_cohort_retention(store: Any, key: str, retention_seconds: int) -> None:
    """Push a cohort's whole record out by another retention period.

    A cohort is three keys, not one: its metadata hash, its member set, and the
    claim on its alias. They have to expire together. Extending only the first
    two lets an actively used cohort outlive the claim on its own name, after
    which a different cohort can take that name and the label stops identifying
    what it labels.

    The alias is read back from the record rather than passed in, so every
    caller extends the same three keys without having to know which alias -- if
    any -- the cohort holds.

    Args:
        store: Redis client for the cohort database.
        key: The cohort's metadata key, as returned by `resolve_cohort`.
        retention_seconds: The lifetime to give each of the three keys.
    """
    store.expire(key, retention_seconds)
    store.expire(f"{key}:jobs", retention_seconds)

    alias = store.hget(key, "alias")
    if alias:
        store.expire(alias_key(alias), retention_seconds)


def cohort_job_ids(store: Any, key: str) -> list[str]:
    """List the jobs belonging to a cohort.

    Args:
        store: Redis client for the cohort database.
        key: The cohort's metadata key, as returned by `resolve_cohort`.

    Returns:
        list[str]: The member job identifiers, empty if the cohort has none.
    """
    return list(store.smembers(f"{key}:jobs") or [])
