"""Declared sample-to-capture association for finite cutoff grid replays."""

from __future__ import annotations

import csv
import io
import logging
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_cohort_manifest import CohortSample
from vntyper.scripts.calibration_secure_io import read_regular_path

logger = logging.getLogger(__name__)

_CALLER_SETS: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {"kestrel": ("kestrel",), "advntr": ("advntr",), "both": ("kestrel", "advntr")}
)
_NATIVE_KESTREL = "native_kestrel"


@dataclass(frozen=True)
class CutoffCaptureInputs:
    """Immutable capture association; the mapping key is the declared sample identity.

    Capture payloads carry no internal sample identity, so this external association is
    the only statement of which evidence belongs to which sample. Identical capture
    bytes are therefore not a duplicate: two samples without surviving candidates
    legitimately serialize the same way.
    """

    kestrel: Mapping[str, Path]
    native_kestrel: Mapping[str, Path | None]
    advntr: Mapping[str, Path]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _require_regular(path: Path, label: str) -> None:
    try:
        read_regular_path(path)
    except (OSError, ValueError) as error:
        message = f"cutoff capture {label} is missing, unsafe, or unreadable"
        logger.error(message)
        raise ValueError(message) from error


def read_cutoff_captures(manifest: Path, samples: tuple[CohortSample, ...], *, caller: str) -> CutoffCaptureInputs:
    """Associate declared cohort samples with the capture files that replay them.

    Every declared identity is checked against the cohort roster, and every declared
    path against the paths already claimed, before any referenced capture is opened: a
    manifest that names a sample the cohort never declared is rejected without reading
    a byte of unaccounted evidence. Paths are interpreted relative to the manifest's own
    directory and returned absolute, so the association cannot silently follow the
    working directory of whoever runs the replay.

    One manifest can serve every caller selection. The columns belonging to a caller
    that is not selected (its ``<caller>_capture`` column, and ``native_kestrel`` when
    Kestrel is not selected) are accepted and ignored: their values may be empty, and
    they are never resolved, opened, checked for existence, counted in the duplicate
    path check or returned. Unknown, duplicated and missing required columns are still
    refused.

    Args:
        manifest: Local TSV with a required ``sample_id`` column, one
            ``<caller>_capture`` column per required caller, and an optional
            ``native_kestrel`` column naming an independently retained native final
            Kestrel TSV per sample; ignored columns of unselected callers are allowed.
        samples: Already validated cohort rows forming the declared roster.
        caller: Which capture columns are required: ``kestrel``, ``advntr`` or ``both``.

    Returns:
        Sample-sorted immutable capture mappings, with a native Kestrel entry for every
        Kestrel sample that is null when no native result was declared.

    Raises:
        ValueError: If the caller selector is unsupported, the manifest is unsafe,
            unreadable or malformed, an identity is outside the roster, an identity or a
            path is claimed twice, or a referenced capture is missing or unsafe.
    """
    callers = _CALLER_SETS.get(caller)
    if callers is None:
        _fail("cutoff capture caller must be kestrel, advntr, or both")
    if not isinstance(manifest, Path):
        _fail("cutoff capture manifest must be a Path")
    parent = manifest.resolve().parent
    capture_columns = tuple(f"{name}_capture" for name in callers)
    required = {"sample_id", *capture_columns}
    known = {"sample_id", _NATIVE_KESTREL, *(f"{name}_capture" for name in _CALLER_SETS["both"])}
    optional = known - required
    reads_native = "kestrel" in callers
    reader = csv.DictReader(io.StringIO(read_regular_path(manifest).decode("utf-8-sig")), delimiter="\t")
    fields = reader.fieldnames
    if (
        not fields
        or len(fields) != len(set(fields))
        or not set(fields) >= required
        or set(fields) - (required | optional)
    ):
        _fail("cutoff capture manifest columns must contain required and supported unique fields")
    roster = {sample.sample_id for sample in samples}
    rows: list[tuple[str, dict[str, Path], Path | None]] = []
    identities: set[str] = set()
    claimed: set[Path] = set()
    for raw in reader:
        if None in raw or any(value is None or value != value.strip() for value in raw.values()):
            _fail("cutoff capture manifest rows must match the header and contain trimmed values")
        if any(not raw[column] for column in required):
            _fail("cutoff capture required values must not be empty")
        identity = raw["sample_id"]
        if identity not in roster:
            _fail("cutoff capture manifest sample is absent from the declared cohort roster")
        captures = {column: (parent / raw[column]).resolve() for column in capture_columns}
        native = (parent / raw[_NATIVE_KESTREL]).resolve() if reads_native and raw.get(_NATIVE_KESTREL) else None
        declared = [*captures.values(), *(() if native is None else (native,))]
        if identity in identities or len(set(declared)) != len(declared) or not claimed.isdisjoint(declared):
            _fail("cutoff capture manifest declares a duplicate sample identity or duplicate capture path")
        identities.add(identity)
        claimed.update(declared)
        rows.append((identity, captures, native))
    if not rows:
        _fail("cutoff capture manifest must declare at least one sample")
    kestrel: dict[str, Path] = {}
    native_kestrel: dict[str, Path | None] = {}
    advntr: dict[str, Path] = {}
    for identity, captures, native in sorted(rows, key=lambda row: row[0]):
        for column, path in captures.items():
            _require_regular(path, f"{column} for {identity}")
        if native is not None:
            _require_regular(native, f"native Kestrel result for {identity}")
        if "kestrel_capture" in captures:
            kestrel[identity] = captures["kestrel_capture"]
            native_kestrel[identity] = native
        if "advntr_capture" in captures:
            advntr[identity] = captures["advntr_capture"]
    return CutoffCaptureInputs(MappingProxyType(kestrel), MappingProxyType(native_kestrel), MappingProxyType(advntr))


def primary_samples(rows: tuple[CohortSample, ...]) -> tuple[CohortSample, ...]:
    """Keep one representative per declared group so related samples cannot leak.

    Members of one ``group_id`` are biological duplicates: relatives, repeat libraries or
    the same individual sequenced twice. Scoring a cutoff on all of them would count one
    piece of biological evidence several times and let a group span a train/test split,
    so only the first-seen member of each group is retained. Group membership is never a
    predictor, and a group whose members disagree about the declared genotype is a
    manifest defect rather than a sample to silently drop.

    Args:
        rows: Already validated cohort rows in declaration order.

    Returns:
        The first-seen representative of each group, in declaration order.

    Raises:
        ValueError: If two members of one group declare conflicting genotypes.
    """
    representatives: dict[str, CohortSample] = {}
    for row in rows:
        first = representatives.setdefault(row.group_id, row)
        if first.genotype != row.genotype:
            _fail("cutoff primary samples reject conflicting genotypes inside one declared group")
    return tuple(representatives.values())
