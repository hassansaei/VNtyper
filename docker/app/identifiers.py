"""The shape of the identifiers this service issues, and the check for it.

Every job identifier the API hands out is a version-4 UUID in its canonical
36-character form, minted by `uuid4()` at submission time. Callers hand those
identifiers back on three routes, and each of them uses the value directly: as a
key in the job database, and as part of a filename in the output directory. A
value that is not one of the service's own identifiers therefore cannot name a
job, but it can still name *something* -- another key in the same database,
another file in the same directory.

`is_job_id` is the one place that shape is stated. It is a pure predicate with
no Redis and no framework in it, so the routes can share it and it can be tested
on its own. The API layer decides what to do with a False answer; there is only
one sensible answer, and it is the same one an identifier that names no job
gets.

The pattern is anchored at both ends with `\\Z` rather than `$`, because `$` also
matches immediately before a trailing newline, and it accepts either hex case:
the service issues lower case, but an identifier that has been round-tripped
through a system that upper-cased it is still the same identifier.

Accepting either case is only half of that, and on its own it is worse than
rejecting upper case outright. Redis keys are bytes and the output directory is
a case-sensitive filesystem, so an upper-cased identifier used *unchanged*
validates and then matches nothing -- the caller gets 404, the answer reserved
for an identifier that names no job, for one that does. `canonical_id` closes
that: it is the validator, and it returns the identifier in the form the service
issued rather than the form the caller typed. Every caller uses its return
value; the `is_*` predicates are thin wrappers over it, kept because the shape
alone is worth stating and testing on its own.
"""

import logging
import re

logger = logging.getLogger(__name__)

# The canonical form of a UUID: five hyphen-separated groups of hex digits.
_JOB_ID_PATTERN = re.compile(r"[0-9a-fA-F]{8}(?:-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}\Z")


def canonical_id(value: str | None) -> str | None:
    """Validate an identifier and return it in the form this service issues.

    Args:
        value: The candidate identifier, as supplied by a caller. Untrusted.

    Returns:
        str | None: The identifier lower-cased, or None when the value is not one
            of this service's identifiers at all. The lower-casing is the whole
            point: the value goes on to be a Redis key and a path segment, both
            case-sensitive, so returning it as typed would validate an identifier
            and then fail to find it.
    """
    if not value or not _JOB_ID_PATTERN.match(value):
        return None
    return value.lower()


def is_job_id(value: str | None) -> bool:
    """Report whether a value has the form of an identifier this service issued.

    Args:
        value: The candidate identifier, as supplied by a caller. Untrusted.

    Returns:
        bool: True only for the canonical 36-character UUID form. An absent
            value is False, as is anything carrying a path separator, a glob, a
            control character or any other content the service never mints.
    """
    return canonical_id(value) is not None


def is_cohort_id(value: str | None) -> bool:
    """Report whether a value has the form of a cohort identifier this service issued.

    Cohort identifiers are minted by the same `uuid4()` call and used the same
    way -- interpolated into a Redis key. A value that is not one cannot name a
    cohort, but it can still name another key in the same database: the cohort's
    own member Set is `cohort:<id>:jobs`, so `<id>:jobs` reads a Set as a hash
    and the client raises WRONGTYPE, which the API can only report as its own
    failure.

    Args:
        value: The candidate identifier, as supplied by a caller. Untrusted.

    Returns:
        bool: True only for the canonical 36-character UUID form.
    """
    return is_job_id(value)
