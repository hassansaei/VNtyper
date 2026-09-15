"""Strict runtime sidecar for one selected calibrated adVNTR policy."""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import NoReturn

from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, caller_policy_values_document
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_PREFIX = "/components/advntr/calibrated_calling/"
_FIELDS = {
    "schema_version",
    "mode",
    "cutoff",
    "minimum_read_support",
    "rare_unit_fraction",
    "adapter_filter",
    "minimum_read_match_ratio",
    "prune_reverse",
    "background_sha256",
    "model_sha256",
    "capture_policy_sha256",
    "advntr_revision",
}
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_REVISION = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")


@dataclass(frozen=True)
class AdvntrRuntimePolicy:
    """Selected caller values plus immutable native runtime asset identities."""

    mode: str
    cutoff: float
    minimum_read_support: int
    rare_unit_fraction: float | None
    adapter_filter: bool
    minimum_read_match_ratio: float
    prune_reverse: bool
    background_sha256: str | None
    model_sha256: str
    capture_policy_sha256: str
    advntr_revision: str
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str, *, nullable: bool = False) -> str | None:
    if nullable and value is None:
        return None
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"adVNTR runtime {label} must be a lowercase SHA256 digest")
    return value


def _fraction(value: object, label: str, *, positive: bool) -> float:
    if not isinstance(value, float) or not math.isfinite(value) or not 0 <= value <= 1 or (positive and value == 0):
        _fail(f"adVNTR runtime {label} must be a finite float in its declared domain")
    return value


def decode_advntr_runtime_policy(value: object) -> AdvntrRuntimePolicy:
    """Decode the closed VNtyper runtime sidecar without opening bound assets.

    Args:
        value: Exact JSON-compatible sidecar object.

    Returns:
        Immutable selected values and native runtime commitments.

    Raises:
        ValueError: If fields, numeric domains, or conditional bindings differ.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        _fail("adVNTR runtime policy fields differ from the closed contract")
    if value["schema_version"] != "vntyper-advntr-calibrated-calling-v1":
        _fail("adVNTR runtime policy schema is unsupported")
    mode = value["mode"]
    if not isinstance(mode, str) or mode not in {"legacy", "exact"}:
        _fail("adVNTR runtime mode must be legacy or exact")
    support = value["minimum_read_support"]
    if isinstance(support, bool) or not isinstance(support, int) or support < 1:
        _fail("adVNTR runtime minimum_read_support must be a positive integer")
    rare = value["rare_unit_fraction"]
    if rare is not None:
        rare = _fraction(rare, "rare_unit_fraction", positive=True)
    adapter, prune = value["adapter_filter"], value["prune_reverse"]
    if not isinstance(adapter, bool) or not isinstance(prune, bool):
        _fail("adVNTR runtime adapter_filter and prune_reverse must be boolean")
    background = _digest(value["background_sha256"], "background", nullable=True)
    if mode == "exact" and background is None:
        _fail("adVNTR runtime exact mode requires a background")
    if mode == "legacy" and background is not None:
        _fail("adVNTR runtime legacy mode forbids a background")
    revision = value["advntr_revision"]
    if not isinstance(revision, str) or _REVISION.fullmatch(revision) is None:
        _fail("adVNTR runtime revision must be a full lowercase source revision")
    decoded = AdvntrRuntimePolicy(
        mode,
        _fraction(value["cutoff"], "cutoff", positive=True),
        support,
        rare,
        adapter,
        _fraction(value["minimum_read_match_ratio"], "minimum_read_match_ratio", positive=True),
        prune,
        background,
        str(_digest(value["model_sha256"], "model")),
        str(_digest(value["capture_policy_sha256"], "capture policy")),
        revision,
        "",
    )
    return replace(decoded, sha256=canonical_sha256(_document(decoded)))


def _document(policy: AdvntrRuntimePolicy) -> dict[str, object]:
    return {
        "schema_version": "vntyper-advntr-calibrated-calling-v1",
        "mode": policy.mode,
        "cutoff": policy.cutoff,
        "minimum_read_support": policy.minimum_read_support,
        "rare_unit_fraction": policy.rare_unit_fraction,
        "adapter_filter": policy.adapter_filter,
        "minimum_read_match_ratio": policy.minimum_read_match_ratio,
        "prune_reverse": policy.prune_reverse,
        "background_sha256": policy.background_sha256,
        "model_sha256": policy.model_sha256,
        "capture_policy_sha256": policy.capture_policy_sha256,
        "advntr_revision": policy.advntr_revision,
    }


def advntr_runtime_policy_document(policy: AdvntrRuntimePolicy) -> dict[str, object]:
    """Project a sidecar after revalidating all direct or replaced content."""
    if not isinstance(policy, AdvntrRuntimePolicy):
        _fail("adVNTR runtime policy projection requires AdvntrRuntimePolicy")
    document = _document(policy)
    if decode_advntr_runtime_policy(document) != policy:
        _fail("adVNTR runtime policy differs from its canonical content")
    return document


def build_advntr_runtime_policy(
    caller_policy: CallerPolicyValues,
    *,
    model_sha256: str,
    background_sha256: str | None,
    capture_policy_sha256: str,
    advntr_revision: str,
) -> AdvntrRuntimePolicy:
    """Bind the selected full caller policy to exact native runtime assets."""
    caller_policy_values_document(caller_policy)
    if "advntr" not in caller_policy.required_callers:
        _fail("adVNTR runtime policy requires adVNTR in selected callers")
    values = caller_policy.values
    document: dict[str, object] = {
        "schema_version": "vntyper-advntr-calibrated-calling-v1",
        **{
            name: values[f"{_PREFIX}{name}"]
            for name in (
                "mode",
                "cutoff",
                "minimum_read_support",
                "rare_unit_fraction",
                "adapter_filter",
                "minimum_read_match_ratio",
                "prune_reverse",
            )
        },
        "background_sha256": background_sha256,
        "model_sha256": model_sha256,
        "capture_policy_sha256": capture_policy_sha256,
        "advntr_revision": advntr_revision,
    }
    return decode_advntr_runtime_policy(document)


def validate_advntr_runtime_policy(
    runtime: AdvntrRuntimePolicy,
    caller_policy: CallerPolicyValues,
    *,
    background_raw_sha256: str | None,
) -> AdvntrRuntimePolicy:
    """Require exact selected values and the actually opened background digest."""
    document = advntr_runtime_policy_document(runtime)
    try:
        rebuilt = build_advntr_runtime_policy(
            caller_policy,
            model_sha256=runtime.model_sha256,
            background_sha256=runtime.background_sha256,
            capture_policy_sha256=runtime.capture_policy_sha256,
            advntr_revision=runtime.advntr_revision,
        )
    except ValueError as error:
        raise ValueError("adVNTR runtime values differ from the selected caller policy") from error
    if rebuilt != runtime or document != advntr_runtime_policy_document(rebuilt):
        _fail("adVNTR runtime values differ from the selected caller policy")
    if background_raw_sha256 != runtime.background_sha256:
        _fail("adVNTR runtime background bytes differ from their selected commitment")
    return runtime
