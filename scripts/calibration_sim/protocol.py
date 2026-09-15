"""Frozen finite simulation designs; declarations confer no evidence authority."""

from __future__ import annotations

import hashlib
import json
import logging
import re
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NoReturn, cast

from .haplotypes import DiploidTruth, RepeatEdit, SimulatedHaplotype, build_haplotype, diploid_truth
from .reads import generate_read_pairs

logger = logging.getLogger(__name__)
_TOKEN = re.compile(r"[a-zA-Z0-9][a-zA-Z0-9_.-]{0,127}\Z")
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_ROLES = {"development-assessment", "training", "policy-selection", "validation", "locked-heldout"}
_FIELDS = {
    "schema_version",
    "purpose",
    "repeat_unit_bp",
    "maximum_haplotype_bp",
    "maximum_case_count",
    "maximum_generated_bases",
    "acceptance_protocol_sha256",
    "cases",
}
_CASE_FIELDS = {
    "case_id",
    "group_key",
    "backbone_family",
    "pair_family",
    "seed_family",
    "role",
    "primary",
    "strata",
    "scenario",
    "caller_positive",
    "truth_variant_ids",
    "haplotypes",
    "reads",
}
_READ_FIELDS = {
    "pair_count",
    "read_length",
    "fragment_length",
    "seed",
    "substitution_rate",
    "quality_score",
    "allele_copy_weights",
}


@dataclass(frozen=True)
class ReadDesign:
    """Finite read-generator arguments independent of production measurement."""

    pair_count: int
    read_length: int
    fragment_length: int
    seed: int
    substitution_rate: float
    quality_score: int
    allele_copy_weights: tuple[float, float]


@dataclass(frozen=True)
class SimulationCase:
    """One technical arm, explicit family identity and independent sequence truth."""

    case_id: str
    group_key: str
    backbone_family: str
    pair_family: str
    seed_family: str
    role: str
    primary: bool
    strata: tuple[str, ...]
    scenario: str
    caller_positive: bool
    truth_variant_ids: tuple[str, ...]
    haplotypes: tuple[SimulatedHaplotype, SimulatedHaplotype]
    truth: DiploidTruth
    reads: ReadDesign


@dataclass(frozen=True)
class SimulationProtocol:
    """Hashed declaration with derived budgets; never a scientific attestation."""

    purpose: str
    maximum_generated_bases: int
    cases: tuple[SimulationCase, ...]
    generated_bases: int
    pair_count: int
    independent_group_count: int
    canonical_json: bytes
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"simulation {name} must contain exactly its declared fields")
    return value


def _integer(value: object, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= 2**53 - 1:
        _fail(f"simulation {name} must be an integer within declared JSON-safe limits")
    return value


def _token(value: object, name: str) -> str:
    if not isinstance(value, str) or _TOKEN.fullmatch(value) is None:
        _fail(f"simulation {name} must be a safe nonempty identifier")
    return value


def _strings(value: object, name: str, *, empty: bool = False) -> tuple[str, ...]:
    if not isinstance(value, list) or (not value and not empty):
        _fail(f"simulation {name} must be a list with required members")
    result = tuple(_token(item, name) for item in value)
    if list(result) != sorted(set(result)):
        _fail(f"simulation {name} must have unique sorted members")
    return result


def _haplotype(value: object, unit_bp: int, maximum_bp: int) -> SimulatedHaplotype:
    raw = _object(value, {"repeat_units", "left_flank", "right_flank", "edits"}, "haplotype")
    units, edits = raw["repeat_units"], raw["edits"]
    if not isinstance(units, list) or not isinstance(edits, list):
        _fail("simulation haplotype units and edits must be lists")
    decoded_edits = []
    for edit in edits:
        row = _object(edit, {"repeat_index", "offset", "deleted_bases", "inserted_bases"}, "edit")
        decoded_edits.append(
            RepeatEdit(
                repeat_index=cast(int, row["repeat_index"]),
                offset=cast(int, row["offset"]),
                deleted_bases=cast(str, row["deleted_bases"]),
                inserted_bases=cast(str, row["inserted_bases"]),
            )
        )
    return build_haplotype(
        repeat_units=tuple(units),
        repeat_unit_bp=unit_bp,
        left_flank=cast(str, raw["left_flank"]),
        right_flank=cast(str, raw["right_flank"]),
        maximum_haplotype_bp=maximum_bp,
        edits=tuple(decoded_edits),
    )


def _reads(value: object, haplotypes: tuple[SimulatedHaplotype, SimulatedHaplotype], budget: int) -> ReadDesign:
    raw = _object(value, _READ_FIELDS, "read design")
    weights = raw["allele_copy_weights"]
    if not isinstance(weights, list) or len(weights) != 2:
        _fail("simulation allele copy weights must contain exactly two values")
    design = ReadDesign(
        pair_count=_integer(raw["pair_count"], "pair count"),
        read_length=_integer(raw["read_length"], "read length", 1),
        fragment_length=_integer(raw["fragment_length"], "fragment length", 1),
        seed=_integer(raw["seed"], "seed"),
        substitution_rate=cast(float, raw["substitution_rate"]),
        quality_score=_integer(raw["quality_score"], "quality score"),
        allele_copy_weights=cast(tuple[float, float], tuple(weights)),
    )
    # The generator validates eagerly; obtaining this iterator emits no reads.
    generate_read_pairs(
        haplotypes=haplotypes,
        maximum_generated_bases=budget,
        pair_count=design.pair_count,
        read_length=design.read_length,
        fragment_length=design.fragment_length,
        seed=design.seed,
        substitution_rate=design.substitution_rate,
        quality_score=design.quality_score,
        allele_copy_weights=design.allele_copy_weights,
    )
    return design


def _case(value: object, unit_bp: int, maximum_bp: int, budget: int) -> SimulationCase:
    raw = _object(value, _CASE_FIELDS, "case")
    tokens = {
        name: _token(raw[name], name)
        for name in (
            "case_id",
            "group_key",
            "backbone_family",
            "pair_family",
            "seed_family",
            "role",
            "scenario",
        )
    }
    if tokens["role"] not in _ROLES or tokens["scenario"] not in {"nominal", "stress"}:
        _fail("simulation role or scenario is not declared")
    if type(raw["primary"]) is not bool or type(raw["caller_positive"]) is not bool:
        _fail("simulation primary and caller truth must be exact booleans")
    variants = _strings(raw["truth_variant_ids"], "variant IDs", empty=True)
    if bool(variants) != raw["caller_positive"]:
        _fail("simulation caller truth must agree with its explicitly supplied variant identities")
    raw_haplotypes = raw["haplotypes"]
    if not isinstance(raw_haplotypes, list) or len(raw_haplotypes) != 2:
        _fail("simulation requires exactly two haplotypes")
    haplotypes = (
        _haplotype(raw_haplotypes[0], unit_bp, maximum_bp),
        _haplotype(raw_haplotypes[1], unit_bp, maximum_bp),
    )
    return SimulationCase(
        **tokens,
        primary=cast(bool, raw["primary"]),
        strata=_strings(raw["strata"], "strata"),
        caller_positive=cast(bool, raw["caller_positive"]),
        truth_variant_ids=variants,
        haplotypes=haplotypes,
        truth=diploid_truth(*haplotypes),
        reads=_reads(raw["reads"], haplotypes, budget),
    )


def _families(cases: tuple[SimulationCase, ...]) -> int:
    families: dict[tuple[str, object], tuple[str, str]] = {}
    primaries: Counter[str] = Counter()
    groups: set[str] = set()
    for case in cases:
        group = (case.group_key, case.role)
        groups.add(case.group_key)
        primaries[case.group_key] += int(case.primary)
        # Explicit genealogies remain required. These extra mechanical collisions
        # prevent a seed, ancestral haplotype, or paired backbone being renamed
        # into independence.
        haplotype_backbones = tuple(
            (item.repeat_unit_bp, item.repeat_units, item.left_flank, item.right_flank) for item in case.haplotypes
        )
        for key in (
            ("group", case.group_key),
            ("backbone", case.backbone_family),
            ("pair", case.pair_family),
            ("seed-family", case.seed_family),
            ("seed-value", case.reads.seed),
            ("backbone-bases", tuple(sorted(haplotype_backbones))),
            *(("haplotype-backbone", backbone) for backbone in haplotype_backbones),
        ):
            previous = families.setdefault(key, group)
            if previous != group:
                _fail("simulation shared families must remain in one leakage group and role")
    if any(primaries[group] != 1 for group in groups):
        _fail("simulation each independent group requires exactly one prespecified primary case")
    return len(groups)


def decode_simulation_protocol(value: object) -> SimulationProtocol:
    """Validate a finite design before generating reads or accessing outcomes.

    Args:
        value: Closed JSON-like declaration with sequence templates and read arms.

    Returns:
        Immutable design and exact resource totals. Study-design roles record
        intent only; external custody and power checks remain separate gates.

    Raises:
        ValueError: On malformed declarations, exceeded budgets, conflicting
            family assignments, or attempts to treat a development smoke as held out.
    """
    raw = _object(value, _FIELDS, "protocol")
    if raw["schema_version"] != "calibration-simulation-protocol-v1":
        _fail("simulation protocol schema is not supported")
    purpose = raw["purpose"]
    if not isinstance(purpose, str) or purpose not in {"development-smoke", "study-design"}:
        _fail("simulation purpose must declare development-smoke or study-design")
    unit_bp = _integer(raw["repeat_unit_bp"], "repeat unit length", 1)
    maximum_bp = _integer(raw["maximum_haplotype_bp"], "haplotype budget", 1)
    maximum_cases = _integer(raw["maximum_case_count"], "case budget", 1)
    budget = _integer(raw["maximum_generated_bases"], "generated base budget")
    protocols = _object(raw["acceptance_protocol_sha256"], {"callers", "length"}, "acceptance protocols")
    if any(not isinstance(digest, str) or _DIGEST.fullmatch(digest) is None for digest in protocols.values()):
        _fail("simulation acceptance protocols must be exact SHA256 identities")
    rows = raw["cases"]
    if not isinstance(rows, list) or not 1 <= len(rows) <= maximum_cases:
        _fail("simulation case list is empty or exceeds its declared resource budget")
    cases = tuple(_case(row, unit_bp, maximum_bp, budget) for row in rows)
    ids = [case.case_id for case in cases]
    if ids != sorted(set(ids)):
        _fail("simulation case IDs must be unique and sorted")
    independent_count = _families(cases)
    if purpose == "development-smoke" and any(case.role != "development-assessment" for case in cases):
        _fail("simulation development smoke cannot assign confirmation or fitting roles")
    total_bases = sum(case.reads.pair_count * 2 * case.reads.read_length for case in cases)
    if total_bases > budget:
        _fail("simulation total read bases exceed the declared resource budget")
    try:
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode()
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("simulation declaration must contain finite JSON data") from error
    return SimulationProtocol(
        purpose,
        budget,
        cases,
        total_bases,
        sum(case.reads.pair_count for case in cases),
        independent_count,
        encoded,
        hashlib.sha256(encoded).hexdigest(),
    )


def simulation_protocol_document(protocol: SimulationProtocol) -> dict[str, object]:
    """Project independently copied canonical JSON after verifying derived content.

    Args:
        protocol: Previously decoded immutable simulation declaration.

    Returns:
        Exact original finite protocol document, without derived claims.

    Raises:
        ValueError: If the typed object or any derived identity was altered.
    """
    if not isinstance(protocol, SimulationProtocol) or not isinstance(protocol.canonical_json, bytes):
        _fail("simulation protocol must be a typed canonical declaration")
    try:
        document = json.loads(protocol.canonical_json)
    except (ValueError, UnicodeError) as error:
        raise ValueError("simulation canonical protocol is invalid") from error
    if decode_simulation_protocol(document) != protocol:
        _fail("simulation protocol content or identity differs from its frozen declaration")
    return document
