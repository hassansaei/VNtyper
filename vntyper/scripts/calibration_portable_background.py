"""Pre-binding projection of native background parameters without cohort narrative."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from typing import NoReturn

logger = logging.getLogger(__name__)
PORTABLE_BACKGROUND_PROVENANCE = "VNtyper calibrated background; portable runtime parameters only"
_FIELDS = {"schema", "version", "provenance", "default_probability", "states"}


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _probability(value: object) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value < 1:
        _fail("portable background probabilities must be finite and strictly between zero and one")
    return value


def project_portable_background(value: object) -> dict[str, object]:
    """Replace only free narrative provenance before any candidate bytes are bound.

    Args:
        value: Closed native fitter background document with state probabilities.

    Returns:
        Fresh canonicalizable document preserving every numerical parameter.
        Native state-grammar compatibility remains the native loader's obligation.

    Raises:
        ValueError: If source fields or model parameter types/domains differ.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        _fail("portable background source fields differ from the closed native parameter contract")
    if value["schema"] != "advntr.frameshift.background" or type(value["version"]) is not int or value["version"] != 1:
        _fail("portable background source schema or version differs")
    if not isinstance(value["provenance"], str) or not value["provenance"].strip():
        _fail("portable background source provenance must be nonempty text")
    states = value["states"]
    if not isinstance(states, Mapping):
        _fail("portable background states must be an object")
    projected = {}
    for name, probability in states.items():
        if not isinstance(name, str) or not name or any(part != part.strip() or not part for part in name.split("&")):
            _fail("portable background state names must be nonempty exact lookup keys")
        projected[name] = _probability(probability)
    return {
        "schema": value["schema"],
        "version": value["version"],
        "provenance": PORTABLE_BACKGROUND_PROVENANCE,
        "default_probability": _probability(value["default_probability"]),
        "states": projected,
    }


def validate_portable_background(value: object) -> dict[str, object]:
    """Require an already-projected background without rewriting frozen payload bytes.

    Args:
        value: Parsed portable background whose original bytes are separately hashed.

    Returns:
        Validated numerical parameters and fixed public provenance.

    Raises:
        ValueError: If the source is malformed or still carries narrative provenance.
    """
    projected = project_portable_background(value)
    if value != projected:
        _fail("portable background provenance must be projected before candidate payload binding")
    return projected
