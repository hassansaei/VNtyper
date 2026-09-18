"""Subprocess adapter for the packaged adVNTR recipe-v1 background fitter."""

from __future__ import annotations

import hashlib
import logging
import math
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn

from vntyper.modules.advntr.advntr_calibration_policy import (
    AdvntrCapabilities,
    AdvntrToolPin,
    ProcessRunner,
    probe_advntr_capabilities,
)
from vntyper.scripts.calibration_secure_io import SecureDirectoryReader, read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)

_SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
_BACKGROUND_FIELDS = {"schema", "version", "provenance", "default_probability", "states"}


@dataclass(frozen=True)
class BackgroundFitPlan:
    """One immutable recipe-v1 fit invocation over controller-staged training data."""

    argv_prefix: tuple[str, ...]
    argv: tuple[str, ...]
    capture_root: Path
    labels_path: Path
    diagnostic_policy_path: Path
    output_directory: Path
    partition: str
    profile: str
    folds: int
    insert_lengths: int
    source_cohort: str
    design: str


@dataclass(frozen=True)
class BackgroundFitResult:
    """Validated upstream fit artifacts; probabilities remain upstream-owned bytes."""

    capabilities: AdvntrCapabilities
    background_sha256: str
    diagnostic_policy_sha256: str
    artifact_sha256: MappingProxyType[str, str]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _absolute(value: object, label: str) -> Path:
    if not isinstance(value, Path) or not value.is_absolute():
        _fail(f"{label} must be an absolute Path")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value or "\x00" in value or "\n" in value:
        _fail(f"{label} must be nonempty trimmed single-line text")
    return value


def _count(value: object, label: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _fail(f"{label} must be an integer >= {minimum}")
    return value


def fit_output_names(profile: str) -> tuple[str, ...]:
    """Return the exact packaged fit-background output inventory."""
    if not isinstance(profile, str) or _SAFE_ID.fullmatch(profile) is None:
        _fail("adVNTR fit profile must be one safe identifier")
    profiled = tuple(
        f"{profile}.{suffix}"
        for suffix in (
            "background.json",
            "build-report.json",
            "build-report.md",
            "cv.json",
            "falsification.json",
            "predictions.json",
            "sidecar.json",
            "states.tsv",
        )
    )
    return (*profiled, "loader-refusal-probe.json")


def build_background_fit_plan(
    argv_prefix: tuple[str, ...],
    *,
    capture_root: Path,
    labels_path: Path,
    diagnostic_policy_path: Path,
    output_directory: Path,
    partition: str,
    profile: str,
    folds: int,
    insert_lengths: int,
    source_cohort: str,
    design: str,
) -> BackgroundFitPlan:
    """Build exact shell-free argv for the frozen installed estimator."""
    if (
        not isinstance(argv_prefix, tuple)
        or not argv_prefix
        or any(not isinstance(token, str) or not token or "\x00" in token for token in argv_prefix)
    ):
        _fail("adVNTR executable prefix must be a nonempty tuple of command tokens")
    capture_root = _absolute(capture_root, "adVNTR fit capture root")
    labels_path = _absolute(labels_path, "adVNTR fit labels")
    diagnostic_policy_path = _absolute(diagnostic_policy_path, "adVNTR diagnostic policy")
    output_directory = _absolute(output_directory, "adVNTR fit output")
    if _SAFE_ID.fullmatch(_text(partition, "adVNTR fit partition")) is None:
        _fail("adVNTR fit partition must be one safe identifier")
    fit_output_names(profile)
    folds = _count(folds, "adVNTR fit folds", 2)
    insert_lengths = _count(insert_lengths, "adVNTR fit insertion lengths", 1)
    source_cohort = _text(source_cohort, "adVNTR fit source cohort")
    design = _text(design, "adVNTR fit design")
    argv = (
        *argv_prefix,
        "fit-background",
        "--capture-root",
        str(capture_root),
        "--labels",
        str(labels_path),
        "--partition",
        partition,
        "--out-dir",
        str(output_directory),
        "--profile",
        profile,
        "--profile-version",
        "1",
        "--folds",
        str(folds),
        "--insert-lengths",
        str(insert_lengths),
        "--sink-name",
        "calibration.jsonl",
        "--source-cohort",
        source_cohort,
        "--design",
        design,
        "--background-recipe",
        "recipe-v1",
        "--diagnostic-policy",
        str(diagnostic_policy_path),
    )
    return BackgroundFitPlan(
        argv_prefix,
        argv,
        capture_root,
        labels_path,
        diagnostic_policy_path,
        output_directory,
        partition,
        profile,
        folds,
        insert_lengths,
        source_cohort,
        design,
    )


def _require_plan(plan: BackgroundFitPlan) -> BackgroundFitPlan:
    if not isinstance(plan, BackgroundFitPlan):
        _fail("adVNTR background fit requires BackgroundFitPlan")
    rebuilt = build_background_fit_plan(
        plan.argv_prefix,
        capture_root=plan.capture_root,
        labels_path=plan.labels_path,
        diagnostic_policy_path=plan.diagnostic_policy_path,
        output_directory=plan.output_directory,
        partition=plan.partition,
        profile=plan.profile,
        folds=plan.folds,
        insert_lengths=plan.insert_lengths,
        source_cohort=plan.source_cohort,
        design=plan.design,
    )
    if rebuilt != plan:
        _fail("adVNTR background fit plan differs from its typed argv")
    return plan


def _caller_policy(raw: bytes) -> dict[str, object]:
    document = load_strict_json_object(raw)
    if set(document) != {"schema_version", "mode", "cutoff", "minimum_read_support"}:
        _fail("adVNTR diagnostic policy fields differ")
    cutoff = document["cutoff"]
    support = document["minimum_read_support"]
    if (
        document["schema_version"] != "advntr-frameshift-policy-v1"
        or document["mode"] not in {"legacy", "exact"}
        or not isinstance(cutoff, float)
        or not math.isfinite(cutoff)
        or not 0 < cutoff < 1
        or isinstance(support, bool)
        or not isinstance(support, int)
        or support < 1
    ):
        _fail("adVNTR diagnostic policy values are invalid")
    return document


def _probability(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value < 1:
        _fail(f"{label} must be a finite probability strictly between zero and one")
    return float(value)


def _validate_background(document: dict[str, object]) -> None:
    if (
        set(document) != _BACKGROUND_FIELDS
        or document["schema"] != "advntr.frameshift.background"
        or document["version"] != 1
    ):
        _fail("adVNTR background fields or schema are unsupported")
    _text(document["provenance"], "adVNTR background provenance")
    _probability(document["default_probability"], "adVNTR default background")
    states = document["states"]
    if not isinstance(states, dict):
        _fail("adVNTR background states must be an object")
    for name, value in states.items():
        _text(name, "adVNTR background state")
        _probability(value, "adVNTR state background")


def fit_background(
    plan: BackgroundFitPlan,
    pin: AdvntrToolPin,
    *,
    runner: ProcessRunner = subprocess.run,
) -> BackgroundFitResult:
    """Run the packaged estimator once and bind its exact complete output inventory."""
    plan = _require_plan(plan)
    capabilities = probe_advntr_capabilities(plan.argv_prefix, pin, runner=runner)
    if not plan.capture_root.is_dir() or not plan.labels_path.is_file():
        _fail("adVNTR fit staged capture root or labels are unavailable")
    if os.path.lexists(plan.output_directory):
        _fail("adVNTR fit output directory must be new")
    diagnostic_raw = read_regular_path(plan.diagnostic_policy_path)
    diagnostic = _caller_policy(diagnostic_raw)
    completed = runner(plan.argv, capture_output=True, text=True, check=False)
    if not isinstance(completed, subprocess.CompletedProcess) or completed.returncode != 0:
        raise RuntimeError("adVNTR background fit process failed")
    names = fit_output_names(plan.profile)
    with SecureDirectoryReader.open(plan.output_directory, set(names)) as reader:
        payload = reader.read_files(names)
    background_name = f"{plan.profile}.background.json"
    sidecar_name = f"{plan.profile}.sidecar.json"
    report_name = f"{plan.profile}.build-report.json"
    background_document = load_strict_json_object(payload[background_name])
    _validate_background(background_document)
    sidecar = load_strict_json_object(payload[sidecar_name])
    if (
        sidecar.get("schema") != "advntr-bench.frameshift.background.sidecar"
        or sidecar.get("version") != 1
        or sidecar.get("profile_name") != plan.profile
        or sidecar.get("partition") != plan.partition
        or sidecar.get("background_recipe_id") != "recipe-v1"
        or sidecar.get("diagnostic_policy") != diagnostic
        or sidecar.get("preregistration_overrides") != {}
    ):
        _fail("adVNTR background sidecar recipe or study bindings differ from the fit plan")
    capture_identity = sidecar.get("capture_identity")
    producer = capture_identity.get("producer") if isinstance(capture_identity, dict) else None
    expected_producer = {
        "package_version": capabilities.package_version,
        "build_id": capabilities.build_id,
        "source_revision": capabilities.source_revision,
    }
    if producer != expected_producer:
        _fail("adVNTR background capture producer differs from the capability preflight")
    report = load_strict_json_object(payload[report_name])
    if report.get("background_recipe_id") != "recipe-v1" or report.get("diagnostic_policy") != diagnostic:
        _fail("adVNTR background build report differs from the frozen recipe or diagnostic policy")
    loader_proof = sidecar.get("loader_proof")
    bad_key_probe = loader_proof.get("bad_key_probe") if isinstance(loader_proof, dict) else None
    if not isinstance(bad_key_probe, dict) or bad_key_probe.get("refused") is not True:
        _fail("adVNTR background sidecar lacks the packaged loader-refusal proof")
    load_strict_json_object(payload["loader-refusal-probe.json"])
    hashes = MappingProxyType({name: hashlib.sha256(payload[name]).hexdigest() for name in names})
    return BackgroundFitResult(
        capabilities,
        hashes[background_name],
        hashlib.sha256(diagnostic_raw).hexdigest(),
        hashes,
    )
