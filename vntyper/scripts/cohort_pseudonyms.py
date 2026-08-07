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
``tests/unit/test_cohort_pseudonyms.py`` characterises it and
``tests/unit/test_cohort_identity.py`` states the specification its width answers.
"""

import hashlib
import logging
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
        ValueError: If ``algorithm`` is not available in this interpreter's ``hashlib``,
            if it needs a digest length of its own (the SHAKE family), or if ``length`` is
            not a positive integer no wider than the digest. Both settings come out of a
            JSON configuration, so both are checked; an unknown algorithm is refused by
            name rather than silently falling back, because a silent fallback changes every
            pseudonym in the report without saying so.
    """
    if algorithm not in hashlib.algorithms_available:
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
    if length > len(digest):
        msg = f"Pseudonym digest length {length} exceeds the {len(digest)} hex characters {algorithm} produces"
        logger.error(msg)
        raise ValueError(msg)
    return f"{prefix}{digest[:length]}"
