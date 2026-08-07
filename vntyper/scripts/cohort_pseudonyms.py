"""
vntyper/scripts/cohort_pseudonyms.py

Module Purpose:
---------------
The name a cohort sample is reported under once ``--pseudonymize-samples`` is asked for.

The pseudonym replaces the sample's identity in every table, plot and export, and
``pseudonymization_table.tsv`` written beside the report is the only way back. That makes
two things load-bearing, and neither of them touches the filesystem: the digest has to be
wide enough that two patients cannot land on one name (#206), and the algorithm and width
- both read out of ``config["cohort"]["pseudonym"]`` - have to be refused *by name* when
they are unusable rather than silently substituted, because a silent substitution changes
every pseudonym in a report without saying so.

Split out of ``cohort_inputs.py``, which #205's identity work grew to 499 lines. That
module finds the samples of a cohort on disk and reads what the pipeline recorded for
each; naming one afterwards is a separate, purely arithmetic job, and this is the module
AGENTS.md rule 3 asks new pure logic to go into rather than back into the file it came
from.

Two things deliberately stayed behind. ``DiscoveredSample`` and
``discover_sample_directories`` are discovery, not pseudonymisation. So is
``duplicate_identity``: the collision it looks for is between two *identities*, it is
upstream of any digest, and no width chosen here can fix it - ``a/sample`` and
``b/sample`` are two patients with one basename (#206). Keeping it there also keeps this
module free of any dependency on the discovery record.

Nothing about the pseudonym changed in the move.
``tests/unit/test_cohort_pseudonyms.py`` characterises it,
``tests/unit/test_cohort_pseudonym_config.py`` states the specification its width answers
and exercises the validation of both settings, and
``tests/unit/test_cohort_identity.py`` states what no width can fix.
"""

import hashlib
import logging
from collections.abc import Mapping
from typing import Any

logger = logging.getLogger(__name__)

#: Digest used to build a cohort pseudonym. Overridable through
#: ``config["cohort"]["pseudonym"]``; declared here so a config that omits the key - which
#: ``--config-path`` produces, because it replaces rather than merges (AGENTS.md trap 2) -
#: does not raise.
DEFAULT_PSEUDONYM_ALGORITHM = "sha256"

#: Hex characters of the digest a pseudonym carries. Twelve is 48 bits: the birthday
#: probability of at least one collision, ``1 - exp(-n(n-1)/2**49)``, is 1.78e-9 at 1,000
#: samples and 1.78e-7 at 10,000. The previous value was five characters of MD5 - 20 bits,
#: which collides with probability ~37.9% at 1,000 samples (#206).
DEFAULT_PSEUDONYM_LENGTH = 12

#: Where the two settings live in the configuration, outermost first. Named once so the
#: warning a malformed level produces and the read itself cannot drift apart.
PSEUDONYM_CONFIG_PATH = ("cohort", "pseudonym")


#: Distinguishes "the key is not there" from "the key is there and holds null". The two
#: are the same to ``.get(key)`` and are not the same to an operator: an absent key is
#: every configuration written before this milestone, while a null one was typed.
_ABSENT = object()


def _mapping_at(container: Mapping[str, Any], key: str, path: str) -> Mapping[str, Any]:
    """Read one nested configuration level, tolerating anything JSON can put there.

    Args:
        container: The mapping to read ``key`` from.
        key: The key to read.
        path: The dotted name of ``key`` for the log message, e.g. ``cohort.pseudonym``.

    Returns:
        Mapping[str, Any]: The nested mapping, or an empty one when the key is absent or
        holds anything that is not a mapping. Only the second case is logged - at
        warning, naming the key - because falling back silently changes every pseudonym
        in the report, while an absent key is the ordinary consequence of AGENTS.md trap
        2 and says nothing about the operator's intent.
    """
    value = container.get(key, _ABSENT)
    if value is _ABSENT:
        return {}
    if isinstance(value, Mapping):
        return value
    shape = "null" if value is None else f"a {type(value).__name__}"
    logger.warning(
        f"Configuration key {path!r} is {shape} rather than a mapping; "
        f"the pseudonym defaults ({DEFAULT_PSEUDONYM_ALGORITHM}, {DEFAULT_PSEUDONYM_LENGTH}) are used instead."
    )
    return {}


def pseudonym_settings(config: Any) -> tuple[Any, Any]:
    """Read the digest algorithm and width out of a configuration.

    The two settings are configuration rather than code, and ``--config-path`` replaces
    the whole configuration rather than merging it (AGENTS.md trap 2) - so every level of
    the read has to survive a hand-written document. A ``.get(key, {})`` chain does not:
    it defends against a key being *absent* and not against it being present and null, and
    ``None.get`` is an ``AttributeError`` naming neither the key nor the file. That raised
    even for a run that had not asked for pseudonyms at all.

    Neither value is validated here. :func:`pseudonymized_sample_name` refuses an unusable
    algorithm or width by name, and it is the one place that can: it is also reached
    directly, with defaults, by callers that never see a configuration.

    Args:
        config: The loaded configuration, or anything at all - a non-mapping is reported
            and treated as empty.

    Returns:
        tuple[Any, Any]: The algorithm and the digest width, each falling back to its
        module default. The types are whatever the JSON held, deliberately: validating
        them here would duplicate the refusal that has to exist in the digest function
        anyway, and would report the failure from the wrong place.
    """
    if not isinstance(config, Mapping):
        logger.warning(
            f"The configuration is a {type(config).__name__} rather than a mapping; "
            f"the pseudonym defaults ({DEFAULT_PSEUDONYM_ALGORITHM}, {DEFAULT_PSEUDONYM_LENGTH}) are used instead."
        )
        return DEFAULT_PSEUDONYM_ALGORITHM, DEFAULT_PSEUDONYM_LENGTH

    outer, inner = PSEUDONYM_CONFIG_PATH
    section = _mapping_at(_mapping_at(config, outer, outer), inner, f"{outer}.{inner}")
    return (
        section.get("algorithm", DEFAULT_PSEUDONYM_ALGORITHM),
        section.get("digest_characters", DEFAULT_PSEUDONYM_LENGTH),
    )


def pseudonymized_sample_name(
    prefix: Any,
    original_sample: str,
    *,
    algorithm: str = DEFAULT_PSEUDONYM_ALGORITHM,
    length: int = DEFAULT_PSEUDONYM_LENGTH,
) -> str:
    """Build the pseudonym a sample is reported under.

    The pseudonym is the caller's prefix followed by the first ``length`` hex digits of
    the digest of the original name, so it is stable across runs and the pseudonymization
    table written beside the report stays meaningful.

    MD5 at five characters was the original scheme and is gone for two reasons: 20 bits
    collides at realistic cohort sizes (#206), and ``hashlib.md5()`` raises on a
    FIPS-enabled build unless it is called with ``usedforsecurity=False``.

    Args:
        prefix: The value ``--pseudonymize-samples`` supplied. Interpolated rather than
            concatenated, so a non-string does not raise.
        original_sample: The sample's identity.
        algorithm: A ``hashlib`` algorithm name.
        length: How many hex characters of the digest to keep.

    Returns:
        str: The pseudonym.

    Raises:
        ValueError: If ``algorithm`` is not a string, if it is a string that is not
            available in this interpreter's ``hashlib``, if it needs a digest length of its
            own (the SHAKE family), if the ``hashlib`` backend refuses it outright (a FIPS
            provider does that to a listed but non-approved digest, raising ``ValueError``
            from ``hashlib.new`` rather than the SHAKE family's ``TypeError``), or if
            ``length`` is not a positive integer no wider than the digest. Both settings
            come out of a JSON configuration, so both are checked; an unknown algorithm is
            refused by name rather than silently falling back, because a silent fallback
            changes every pseudonym in the report without saying so.
    """
    # The type test comes first because the membership test below is `in` against a `set`,
    # which hashes its left operand: a JSON list or object arrived as
    # `TypeError: unhashable type: 'list'` from inside the guard whose whole purpose is to
    # produce a `ValueError`, bypassing both the documented contract and the
    # `logger.error` + `raise` convention. `length` beside it was already ordered this way.
    if not isinstance(algorithm, str) or algorithm not in hashlib.algorithms_available:
        msg = f"Unknown pseudonym digest algorithm: {algorithm}"
        logger.error(msg)
        raise ValueError(msg)
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        msg = f"Pseudonym digest length must be a positive integer, got {length!r}"
        logger.error(msg)
        raise ValueError(msg)
    try:
        digest = hashlib.new(algorithm, original_sample.encode()).hexdigest()
    except TypeError as e:
        # shake_128 and shake_256 are in algorithms_available but take their output length
        # as an argument, so hexdigest() raises rather than returning a digest.
        msg = f"Pseudonym digest algorithm {algorithm} does not produce a fixed-length digest: {e}"
        logger.error(msg)
        raise ValueError(msg) from e
    except ValueError as e:
        # `algorithms_available` lists what this interpreter *knows about*, not what its
        # backend will compute: a FIPS-enforcing OpenSSL provider refuses a non-approved
        # digest at construction and hashlib re-raises that as ValueError. Translated
        # here so the message names the configured algorithm rather than quoting OpenSSL.
        msg = f"Pseudonym digest algorithm {algorithm} was refused by this interpreter's hashlib backend: {e}"
        logger.error(msg)
        raise ValueError(msg) from e
    if length > len(digest):
        msg = f"Pseudonym digest length {length} exceeds the {len(digest)} hex characters {algorithm} produces"
        logger.error(msg)
        raise ValueError(msg)
    return f"{prefix}{digest[:length]}"
