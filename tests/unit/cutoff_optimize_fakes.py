"""The ``vntyper calibrate optimize`` test harness: synthetic cohorts, captures and native grids.

This is a helper module, not a test module: it defines no tests. It is imported by
``test_calibration_cutoff_optimize``, ``test_calibration_cutoff_optimize_advntr``,
``test_calibration_cutoff_document`` and ``test_pipeline_research_advntr``. Every capture
is synthetic and every sample name is invented, so nothing here can carry cohort identity.
"""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Callable, Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any, cast
from unittest.mock import patch

import pandas as pd

from tests.builders import kestrel_config
from tests.unit.advntr_grid_fakes import advntr_grid_result, parity_capture
from tests.unit.test_calibration_cutoff_kestrel import _native_negative, _native_production
from tests.unit.test_calibration_kestrel_replay import _capture, _policy, _raw
from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_caller_policy import (
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_cutoff_advntr_axes import AdvntrVisit, predicted_call
from vntyper.scripts.calibration_kestrel_capture import kestrel_capture_document
from vntyper.scripts.canonical_json import canonical_json_bytes

#: ``(alternate depth, active-region depth)`` pairs whose ratio is the Depth_Score.
DEPTHS: dict[str, tuple[int, int]] = {
    "0.014": (7, 500),
    "0.008": (4, 500),
    "0.006": (3, 500),
    "0.004": (2, 500),
    "0.001": (5, 5000),
    "0.0004": (2, 5000),
}


def _frame(scores: tuple[str, ...]) -> pd.DataFrame:
    """Build a raw Kestrel frame with one candidate row per requested Depth_Score."""
    if not scores:
        return _raw().iloc[0:0]
    rows = []
    for ordinal, score in enumerate(scores):
        alternate, region = DEPTHS[score]
        row = _raw(depth_alt=alternate, depth_region=region)
        row.loc[0, "POS"] = 67 + ordinal
        rows.append(row)
    return pd.concat(rows, ignore_index=True)


#: The adVNTR half of an adVNTR-bearing baseline: legacy calling at cutoff 0.001, support 3.
_ADVNTR_BASELINE: dict[str, object] = {
    "/components/advntr/calibrated_calling/mode": "legacy",
    "/components/advntr/calibrated_calling/cutoff": 0.001,
    "/components/advntr/calibrated_calling/minimum_read_support": 3,
    "/components/advntr/calibrated_calling/rare_unit_fraction": None,
    "/components/advntr/calibrated_calling/adapter_filter": False,
    "/components/advntr/calibrated_calling/minimum_read_match_ratio": 0.6,
    "/components/advntr/calibrated_calling/prune_reverse": False,
}
ADV_CUT = "/components/advntr/calibrated_calling/cutoff"
ADV_SUP = "/components/advntr/calibrated_calling/minimum_read_support"


def advntr_baseline_policy() -> CallerPolicyValues:
    """The standard Kestrel baseline plus the seven adVNTR legacy pointers."""
    kestrel = _policy(kestrel_config())
    document = caller_policy_values_document(kestrel)
    return decode_caller_policy_values(
        {**document, "required_callers": ["advntr", "kestrel"], "values": {**kestrel.values, **_ADVNTR_BASELINE}}
    )


def _write_capture(path: Path, scores: tuple[str, ...], policy: CallerPolicyValues | None = None) -> CallerPolicyValues:
    """Serialize one complete synthetic capture and return its baseline policy."""
    capture = _capture(_frame(scores), kestrel_config(), policy=policy)
    path.write_bytes(canonical_json_bytes(kestrel_capture_document(capture)))
    return capture.baseline_policy


def write_manifests(
    tmp_path: Path,
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]],
    *,
    natives: dict[str, str] | None = None,
    advntr: bool = False,
    advntr_baseline: bool = False,
    advntr_calls: Mapping[str, bool | None] | None = None,
    capture_parameters: Mapping[str, object] | None = None,
) -> tuple[Path, Path]:
    """Write the cohort TSV, the capture TSV and every capture the pair references.

    Args:
        tmp_path: Private directory the synthetic inputs are written into.
        cohort: Sample -> (genotype, Depth_Score labels, group id or None).
        natives: Optional sample -> ``"production"``, ``"negative"`` or ``"corrupt"``.
        advntr: Whether to declare an adVNTR capture column. Without ``advntr_calls`` the
            referenced files are placeholders the mocked native grid never opens.
        advntr_baseline: Whether the Kestrel captures declare an adVNTR-bearing baseline.
        advntr_calls: Sample -> the native baseline adVNTR call its one-record capture
            records, so adVNTR baseline parity is proven against real capture files.
        capture_parameters: Overrides of every adVNTR capture's capture-policy fields.

    Returns:
        The cohort manifest path and the capture association manifest path.
    """
    tmp_path.mkdir(mode=0o700, parents=True, exist_ok=True)
    cohort_path = tmp_path / "cohort.tsv"
    captures_path = tmp_path / "captures.tsv"
    cohort_lines = ["sample_id\tbam\tassembly\tgenotype\tgroup_id"]
    header = "sample_id\tkestrel_capture\tnative_kestrel" + ("\tadvntr_capture" if advntr else "")
    capture_lines = [header]
    for name, (genotype, scores, group) in cohort.items():
        cohort_lines.append(f"{name}\t{name}.bam\tGRCh38\t{genotype}\t{group or ''}")
        capture = tmp_path / f"{name}.capture.json"
        _write_capture(capture, scores, advntr_baseline_policy() if advntr_baseline else None)
        native = ""
        request = (natives or {}).get(name)
        if request is not None:
            native_path = tmp_path / f"{name}.native.tsv"
            if request == "negative":
                _native_negative(native_path)
            else:
                _native_production(
                    native_path, capture, changes={"Depth_Score": "0.99"} if request == "corrupt" else None
                )
            native = native_path.name
        row = f"{name}\t{capture.name}\t{native}"
        if advntr:
            locus = tmp_path / f"{name}.advntr.jsonl"
            if advntr_calls is None:
                locus.write_text("{}\n", encoding="utf-8")
            else:
                parity_capture(
                    locus, 17 + len(capture_lines), advntr_calls.get(name, False), capture_parameters=capture_parameters
                )
            row += f"\t{locus.name}"
        capture_lines.append(row)
    cohort_path.write_text("\n".join(cohort_lines) + "\n", encoding="utf-8")
    captures_path.write_text("\n".join(capture_lines) + "\n", encoding="utf-8")
    return cohort_path, captures_path


STANDARD_COHORT: dict[str, tuple[str, tuple[str, ...], str | None]] = {
    "specimen-alpha": ("positive", ("0.014",), "family-1"),
    "specimen-alpha-repeat": ("positive", ("0.014",), "family-1"),
    "specimen-bravo": ("positive", ("0.004",), None),
    "specimen-charlie": ("negative", ("0.0004",), None),
    "specimen-delta": ("negative", ("0.001",), None),
    "specimen-echo": ("positive", (), None),
    "specimen-foxtrot": ("unknown", ("0.014",), None),
}


def namespace(cohort_path: Path, captures_path: Path, **overrides: Any) -> argparse.Namespace:
    """Build the parsed-argument namespace the entry point reads."""
    values: dict[str, Any] = {
        "manifest": cohort_path,
        "captures": captures_path,
        "objective": "max-sensitivity-at-specificity",
        "min_sensitivity": None,
        "min_specificity": 1.0,
        "caller": "kestrel",
        "axes": None,
        "max_breakpoints": None,
        "folds": 3,
        "seed": 20260915,
        "workers": 1,
        "advntr_executable": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def run_optimize(
    tmp_path: Path,
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] | None = None,
    *,
    natives: dict[str, str] | None = None,
    **overrides: Any,
) -> tuple[bool, dict[str, Any], Path]:
    """Run the command through the real atomic adapter and read its report back."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = write_manifests(tmp_path, cohort or STANDARD_COHORT, natives=natives)
    output = tmp_path / "derived"
    args = namespace(cohort_path, captures_path, **overrides)
    successful = atomic_output(output, lambda staging: run_cutoff_optimization(args, staging))
    document = json.loads((output / "report.json").read_bytes())
    assert isinstance(document, dict)
    return successful, document, output


def visit(support: int, pvalue: float) -> AdvntrVisit:
    """One scored adVNTR visit with an exact p-value."""
    return AdvntrVisit(support, Fraction(pvalue))


#: Scored adVNTR visits every synthetic native replay reports (support 5, so the baseline
#: support 3 admits each). The decisive legacy p-value ``q`` per sample:
#:
#: ====================  =======  ==========
#: Sample                q        truth
#: ====================  =======  ==========
#: ``specimen-alpha``    0.0004   positive
#: ``specimen-bravo``    0.004    positive
#: ``specimen-delta``    0.006    negative
#: ``specimen-charlie``  0.02     negative
#: ``specimen-echo``     --       positive (no visit: never called)
#: ``specimen-foxtrot``  0.2      unknown
#: ====================  =======  ==========
PROBE_VISITS: dict[str, tuple[AdvntrVisit, ...] | None] = {
    "specimen-alpha": (visit(5, 0.0004),),
    "specimen-bravo": (visit(5, 0.004),),
    "specimen-charlie": (visit(5, 0.02),),
    "specimen-delta": (visit(5, 0.006),),
    "specimen-echo": (),
    "specimen-foxtrot": (visit(5, 0.2),),
}


def up(value: float) -> float:
    """The next float above ``value``: the smallest strict cutoff that calls a sample with p = ``value``."""
    return math.nextafter(value, math.inf)


GridOverride = Callable[[Any, Mapping[str, CallerPolicyValues]], Any]


def native_grid(
    seen: list[dict[str, Any]],
    *,
    visits: Mapping[str, tuple[AdvntrVisit, ...] | None] = PROBE_VISITS,
    override: GridOverride | None = None,
) -> Callable[..., Any]:
    """A native adVNTR grid fake: every policy's calls follow the legacy rule over ``visits``."""

    def grid(capture_paths: Mapping[str, Path], policies: Mapping[str, CallerPolicyValues], **kwargs: Any) -> Any:
        seen.append({"captures": dict(capture_paths), "policies": dict(policies), **kwargs})
        thresholds = {
            pid: (float(cast(float, p.values[ADV_CUT])), int(cast(int, p.values[ADV_SUP])))
            for pid, p in policies.items()
        }
        result = advntr_grid_result(
            dict.fromkeys(policies),
            visits=dict.fromkeys(policies, visits),
            thresholds=thresholds,
            baseline_policy_id=kwargs["baseline_policy_id"],
            capture_paths=capture_paths,
        )
        return result if override is None else override(result, policies)

    return grid


def run_advntr(
    tmp_path: Path,
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] | None = None,
    *,
    seen: list[dict[str, Any]] | None = None,
    visits: Mapping[str, tuple[AdvntrVisit, ...] | None] = PROBE_VISITS,
    override: GridOverride | None = None,
    capture_parameters: Mapping[str, object] | None = None,
    **overrides: Any,
) -> tuple[bool, dict[str, Any], Path]:
    """Run an adVNTR-bearing optimize with both native grid seams replaced by one fake.

    The captures carry an adVNTR-bearing baseline and one-record adVNTR captures whose
    recorded baseline call is the legacy rule's own, so adVNTR baseline parity runs for real.
    """
    from vntyper.scripts import calibration_cutoff_advntr_axes as axes_module
    from vntyper.scripts import calibration_cutoff_optimize as module

    calls = {key: predicted_call(items, 0.001, 3) for key, items in visits.items()}
    cohort_path, captures_path = write_manifests(
        tmp_path,
        cohort or STANDARD_COHORT,
        advntr=True,
        advntr_baseline=True,
        advntr_calls=calls,
        capture_parameters=capture_parameters,
    )
    args = namespace(cohort_path, captures_path, advntr_executable=tmp_path / "advntr", **overrides)
    grid = native_grid([] if seen is None else seen, visits=visits, override=override)
    output = tmp_path / "derived"
    with (
        patch.object(module, "evaluate_advntr_cutoff_grid", grid),
        patch.object(axes_module, "evaluate_advntr_cutoff_grid", grid),
    ):
        successful = atomic_output(output, lambda staging: module.run_cutoff_optimization(args, staging))
    document = json.loads((output / "report.json").read_bytes())
    assert isinstance(document, dict)
    return successful, document, output
