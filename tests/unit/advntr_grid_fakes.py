"""Synthetic native adVNTR replay grids for the cutoff-axis tests; the executable is never run.

This is a helper module, not a test module: it defines no tests and is imported by
``test_calibration_cutoff_advntr_axes`` and ``test_calibration_cutoff_optimize``.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path

from tests.unit.test_calibration_cutoff_advntr import _capture, _json_bytes, _policy
from vntyper.modules.advntr.advntr_calibration_policy import (
    ADVNTR_BACKGROUND_RECIPE_IDS,
    ADVNTR_CAPABILITIES,
    ADVNTR_CAPTURE_SCHEMA_VERSIONS,
    ADVNTR_POLICY_SCHEMA_VERSIONS,
    AdvntrCapabilities,
    advntr_canonical_sha256,
    decode_advntr_capabilities,
)
from vntyper.modules.advntr.advntr_replay import ReplayLocus
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_cutoff_advntr import (
    AdvntrCutoffGridResult,
    AdvntrCutoffPolicyResult,
    AdvntrCutoffSample,
)
from vntyper.scripts.calibration_cutoff_advntr_axes import AdvntrVisit, predicted_call
from vntyper.scripts.canonical_json import load_strict_json_object

#: First synthetic VNTR identifier; sample ``n`` (sorted roster order) uses ``FIRST_VNTR_ID + n``.
FIRST_VNTR_ID = 17

#: The producer every synthetic capture and locus result declares.
_PRODUCER = {"package_version": "2.4.0", "build_id": "a" * 64, "source_revision": "b" * 40}


def synthetic_capabilities(package_version: str = "2.4.0") -> AdvntrCapabilities:
    """A decodable installed-tool identity, so ``advntr_capabilities_document`` accepts it."""
    return decode_advntr_capabilities(
        {
            "schema_version": "advntr-capabilities-v1",
            "package_version": package_version,
            "build_id": _PRODUCER["build_id"],
            "source_revision": _PRODUCER["source_revision"],
            "capabilities": list(ADVNTR_CAPABILITIES),
            "capture_schema_versions": list(ADVNTR_CAPTURE_SCHEMA_VERSIONS),
            "policy_schema_versions": list(ADVNTR_POLICY_SCHEMA_VERSIONS),
            "background_recipe_ids": list(ADVNTR_BACKGROUND_RECIPE_IDS),
        }
    )


def visit_document(disposition: str, support: object, pvalue: object, *, statistic: bool = True) -> dict[str, object]:
    """One native-shaped decision visit (``statistic=False`` gives it a null statistic)."""
    return {
        "disposition": disposition,
        "plan": {"read_support": support},
        "statistic": {"pvalue": pvalue} if statistic else None,
    }


def _locus_raw(vntr_id: int, visits: Sequence[Mapping[str, object]], called: bool | None, record_sha256: str) -> bytes:
    """A complete upstream locus result document (every ``_LOCUS_RESULT_FIELDS`` key)."""
    document = {
        "schema_version": "advntr-frameshift-replay-result-v1",
        "vntr_id": vntr_id,
        "capture_record_sha256": record_sha256,
        "policy_sha256": "c" * 64,
        "capture_producer": dict(_PRODUCER),
        "capture_assets": {},
        "loaded_background_sha256": None,
        "baseline_parity": True,
        "decision_visits": list(visits),
        "calls": [{"state": "I22_2_G_LEN1"}] if called else [],
        "warnings": [],
        "capture_audit": {
            "attribution_outside_trials": ["synthetic audit failure"] if called is None else [],
            "calibrated_policy_domain_errors": [],
        },
    }
    return _json_bytes(document).rstrip(b"\n")


def _visit_documents(visits: Sequence[AdvntrVisit], cutoff: float, support: int) -> list[dict[str, object]]:
    """Native-shaped decision visits: unscored below the support, else called or cutoff."""
    documents: list[dict[str, object]] = []
    for visit in visits:
        if visit.read_support < support:
            documents.append(visit_document("insufficient-read-support", visit.read_support, None, statistic=False))
            continue
        called = visit.pvalue < Fraction(cutoff)
        documents.append(
            {
                "disposition": "called" if called else "cutoff",
                "plan": {"read_support": visit.read_support},
                "statistic": {"pvalue": float(visit.pvalue), "called": called},
            }
        )
    return documents


def _capture_record(path: Path) -> tuple[int, str]:
    """The one record of a synthetic capture: its VNTR ID and its canonical record digest."""
    lines = path.read_bytes().splitlines()
    if len(lines) != 1:
        raise AssertionError(f"synthetic capture {path.name} must hold exactly one record")
    document = load_strict_json_object(lines[0])
    locus = document["locus"]
    assert isinstance(locus, Mapping)
    return int(locus["vntr_id"]), advntr_canonical_sha256(document)


def advntr_grid_result(
    policies: Mapping[str, Mapping[str, bool | None] | None],
    *,
    visits: Mapping[str, Mapping[str, Sequence[AdvntrVisit] | None]] | None = None,
    thresholds: Mapping[str, tuple[float, int]] | None = None,
    baseline_policy_id: str | None = None,
    capture_paths: Mapping[str, Path] | None = None,
) -> AdvntrCutoffGridResult:
    """A synthetic native adVNTR grid result; the executable itself is never run.

    Args:
        policies: Policy ID to its explicit per-sample native calls (``None`` marks an
            unassessable sample), or ``None`` to compute every call with
            :func:`predicted_call` from that policy's ``visits`` and ``thresholds``.
        visits: Policy ID to per-sample scored visits written into each locus result's
            ``decision_visits`` (``None`` marks an unassessable sample); a policy without
            visits has empty decision visits.
        thresholds: Policy ID to its legacy ``(cutoff, minimum_read_support)``; required
            for a policy whose calls are computed, and it also sets the visit
            dispositions and the execution ID. Defaults to the baseline ``(0.001, 3)``.
        baseline_policy_id: The grid's baseline policy; defaults to the first policy.
        capture_paths: Optional sample key to a one-record capture file. When given, each
            locus carries that record's VNTR ID and its real canonical record digest, as a
            native replay would; otherwise VNTR ``17 + n`` in sorted roster order and a
            placeholder digest.

    Returns:
        A grid with one locus per sample and a decodable tool identity.
    """
    visits = visits or {}
    thresholds = thresholds or {}
    records = {key: _capture_record(path) for key, path in (capture_paths or {}).items()}
    policy_rows = []
    for policy_id, explicit in policies.items():
        cutoff, support = thresholds.get(policy_id, (0.001, 3))
        sample_visits = visits.get(policy_id, {})
        if explicit is None:
            if policy_id not in thresholds or policy_id not in visits:
                raise AssertionError(f"policy {policy_id} needs visits and thresholds to compute its calls")
            calls: Mapping[str, bool | None] = {
                key: predicted_call(None if items is None else tuple(items), cutoff, support)
                for key, items in sample_visits.items()
            }
        else:
            calls = explicit
        samples = []
        for index, key in enumerate(sorted(calls)):
            vntr_id, record_sha256 = records.get(key, (FIRST_VNTR_ID + index, "e" * 64))
            items = sample_visits.get(key, ())
            raw = _locus_raw(vntr_id, _visit_documents(items or (), cutoff, support), calls[key], record_sha256)
            locus = ReplayLocus(key, vntr_id, calls[key] is not None, raw, hashlib.sha256(raw).hexdigest())
            samples.append(AdvntrCutoffSample(key, calls[key] is not None, calls[key], (locus,)))
        execution = f"exec-{cutoff!r}-{support}" if policy_id in thresholds else f"exec-{policy_id}"
        policy_rows.append(AdvntrCutoffPolicyResult(policy_id, "c" * 64, execution, tuple(samples)))
    return AdvntrCutoffGridResult(
        Path("/nonexistent"),
        baseline_policy_id or next(iter(policies)),
        "b" * 64,
        synthetic_capabilities(),
        tuple(policy_rows),
        "d" * 64,
    )


def parity_capture(path: Path, vntr_id: int, called: bool | None, baseline: CallerPolicyValues | None = None) -> Path:
    """A one-record capture whose decision visit records the native call; ``None`` fails its audit."""
    _capture(path, vntr_id=vntr_id, baseline=baseline or _policy())
    document = json.loads(path.read_bytes())
    document["decision_visits"] = [
        {"statistic": None, "plan": None, "disposition": "insufficient-read-support"},
        {
            "disposition": "called" if called else "cutoff",
            "mean_coverage": 30.0,
            "plan": {"state": "I22_2_G_LEN1", "read_support": 5},
            "statistic": {"called": bool(called), "pvalue": 0.0004 if called else 0.4},
        },
    ]
    if called is None:
        document["warnings"] = [{"origin": "calibration-audit", "message": "synthetic"}]
    path.write_bytes(_json_bytes(document))
    return path
