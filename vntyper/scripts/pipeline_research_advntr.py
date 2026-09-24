"""Explicit native adVNTR arguments for a derived research decision profile.

A research profile (``resolve_run_configuration(research_profile=...)``) carries the
calibrated adVNTR values in ``components.advntr.calibrated_calling`` but no approved
bundle, so the bundle preflight that normally renders the explicit native policy never
runs. This module renders that policy from the profile itself, under the fixed
CLI-representable capture semantics the calibration captures were produced under.
Only legacy mode is supported: exact mode needs an approved background.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import fields
from types import MappingProxyType
from typing import Final, NoReturn

from vntyper.modules.advntr.advntr_calibration_policy import (
    CapturePolicy,
    capture_policy_for_caller,
    decode_capture_policy,
)
from vntyper.modules.advntr.advntr_capture import calibrated_policy_argv
from vntyper.scripts.calibration_caller_policy import (
    ADVNTR_CALLER_POLICY_POINTERS,
    KESTREL_CALLER_POLICY_POINTERS,
    CallerPolicyValues,
    decode_caller_policy_values,
)
from vntyper.scripts.run_configuration import RunConfiguration

logger = logging.getLogger(__name__)

RESEARCH_CAPTURE_PARAMETERS: Final[Mapping[str, object]] = MappingProxyType(
    {
        "platform": "illumina",
        "frameshift_mode": True,
        "is_haploid": False,
        "caller_mode": "legacy",
        "minimum_read_length": None,
        "prune_reverse": False,
        "filter_adapter_readthrough": False,
        "minimum_read_match_ratio": None,
        "minimum_relative_ru_coverage": None,
        "use_reference_alignment": True,
        "fully_covered_ru_only": False,
        "maximum_error_rate": 0.05,
        "legacy_error_rate": 0.01,
        "mapq_cutoff": 0,
        "base_quality_cutoff": 20,
        "maximum_low_quality_fraction": 0.1,
        "enhanced_hmm": True,
        "trained_hmms": False,
    }
)

_CAPTURE_SCHEMA: Final = "advntr-runtime-capture-policy-v1"
_POLICY_SCHEMA: Final = "calibration-caller-policy-values-v1"
_MODE_POINTER: Final = "/components/advntr/calibrated_calling/mode"
#: Why a research profile cannot run exact-mode adVNTR.
RESEARCH_LEGACY_ONLY: Final[str] = (
    "research decision profiles support only legacy adVNTR calling; exact mode requires an approved calibration bundle"
)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def project_caller_policy(
    components: Mapping[str, object],
    pointers: Iterable[str],
    required_callers: Sequence[str],
) -> CallerPolicyValues:
    """Read calibrated caller pointers back from resolved decision components.

    Args:
        components: Resolved decision components keyed by component name.
        pointers: Exact ``/components/<name>/...`` JSON pointers to read.
        required_callers: Callers the rebuilt policy declares.

    Returns:
        The validated caller policy holding the values found at each pointer.

    Raises:
        ValueError: If a pointer is absent from the components or the values are invalid.
    """
    values: dict[str, object] = {}
    for pointer in pointers:
        parts = pointer.strip("/").split("/")
        node: object = components.get(parts[1])
        for part in parts[2:]:
            if not isinstance(node, Mapping) or part not in node:
                _fail(f"research decision profile does not expose {pointer} at runtime")
            node = node[part]
        values[pointer] = node
    return decode_caller_policy_values(
        {"schema_version": _POLICY_SCHEMA, "required_callers": list(required_callers), "values": values}
    )


def research_capture_policy(caller: CallerPolicyValues, threads: int) -> CapturePolicy:
    """The complete capture policy a research profile runs adVNTR under.

    Args:
        caller: A complete caller policy that includes adVNTR.
        threads: The native ``-t`` value.

    Returns:
        :data:`RESEARCH_CAPTURE_PARAMETERS` with ``threads``, with the caller's adVNTR
        values projected on (``capture_policy_for_caller``).

    Raises:
        ValueError: If the caller policy does not include adVNTR or is invalid.
    """
    baseline = decode_capture_policy(
        {"schema_version": _CAPTURE_SCHEMA, "parameters": {**RESEARCH_CAPTURE_PARAMETERS, "threads": threads}}
    )
    return capture_policy_for_caller(baseline, caller)


def research_capture_differences(capture: CapturePolicy, caller: CallerPolicyValues) -> tuple[str, ...]:
    """The capture fields a research profile could not reproduce at runtime.

    Args:
        capture: The capture policy some evidence was produced under.
        caller: The caller policy whose adVNTR values the evidence was produced with.

    Returns:
        The names of the fields that differ from :func:`research_capture_policy`, in
        ``CapturePolicy`` order; the thread count is not part of the contract.
    """
    expected = research_capture_policy(caller, capture.threads)
    return tuple(
        field.name
        for field in fields(CapturePolicy)
        if field.name not in {"threads", "sha256"} and getattr(capture, field.name) != getattr(expected, field.name)
    )


def research_policy_argv(advntr: Mapping[str, object], kestrel: Mapping[str, object], threads: int) -> tuple[str, ...]:
    """Render the native adVNTR arguments of a research profile's resolved components.

    This is the rendering :func:`research_advntr_policy_argv` performs at runtime, without a
    run configuration, so ``vntyper calibrate optimize`` can prove an exported profile runs.

    Args:
        advntr: The resolved ``advntr`` decision component, carrying ``calibrated_calling``.
        kestrel: The resolved ``kestrel`` decision component.
        threads: The native ``-t`` value.

    Returns:
        The explicit genotype policy arguments.

    Raises:
        ValueError: If the components select exact mode or their values are incomplete.
    """
    # calibrated_calling is only ever generated alongside the Kestrel pointers, so both callers are required.
    caller = project_caller_policy(
        {"advntr": advntr, "kestrel": kestrel},
        (*ADVNTR_CALLER_POLICY_POINTERS, *KESTREL_CALLER_POLICY_POINTERS),
        ("advntr", "kestrel"),
    )
    if caller.values[_MODE_POINTER] != "legacy":
        _fail(RESEARCH_LEGACY_ONLY)
    return tuple(calibrated_policy_argv(research_capture_policy(caller, threads), caller, None))


def research_advntr_policy_argv(configuration: RunConfiguration, threads: int) -> tuple[str, ...] | None:
    """Render explicit native adVNTR arguments for a research decision profile.

    Args:
        configuration: One resolved pipeline decision configuration.
        threads: The resolved native ``-t`` value for this run.

    Returns:
        The explicit genotype policy arguments, or None when the configuration is not
        a research profile with calibrated adVNTR calling (the packaged profile, a
        profile without adVNTR values, or an approved bundle, which has its own path).

    Raises:
        ValueError: If the profile selects exact mode or its values are incomplete.
    """
    if (
        configuration.caller_calibration is not None
        or configuration.decision_profile.profile_kind != "generated"
        or "calibrated_calling" not in configuration.advntr
    ):
        return None
    return research_policy_argv(configuration.advntr, configuration.kestrel, threads)
