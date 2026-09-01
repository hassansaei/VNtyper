"""Render the complete packaged decision profile from reviewed legacy inputs."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object  # noqa: E402
from vntyper.scripts.decision_profile_schema import (  # noqa: E402
    ValidationClass,
    flatten_decision_projection,
)

PACKAGE_ROOT = ROOT / "vntyper"
OUTPUT_ROOT = PACKAGE_ROOT / "profiles"

_CRITICAL_NUMERIC_METADATA: dict[str, tuple[str, str, bool]] = {
    "/components/kestrel/confidence_assignment/depth_score_thresholds/low": (
        "depth-score-ratio",
        "gte",
        True,
    ),
    "/components/kestrel/alt_filtering/gg_depth_score_threshold": ("depth-score-ratio", "gte", True),
    "/components/kestrel/confidence_assignment/depth_score_thresholds/high": (
        "depth-score-ratio",
        "gte",
        True,
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": (
        "alternate-kmer-path-depth",
        "lte",
        True,
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": (
        "alternate-kmer-path-depth",
        "gte",
        True,
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": (
        "alternate-kmer-path-depth",
        "gte-and-upper-exclusive",
        False,
    ),
    "/components/kestrel/confidence_assignment/var_active_region_threshold": (
        "active-region-kmer-depth",
        "lte",
        True,
    ),
    "/components/nomenclature/thresholds/bam_flank": ("base-pairs-per-side", "eq", True),
    "/components/nomenclature/thresholds/bam_thin_haplotype_record_support": (
        "resolved-haplotype-records",
        "lt",
        False,
    ),
    "/components/nomenclature/identity_reconciliation/kestrel_min_alternate_kmer_path_depth": (
        "alternate-kmer-path-depth",
        "gte",
        True,
    ),
    "/components/nomenclature/identity_reconciliation/advntr_min_sequencing_read_support": (
        "sequencing-reads",
        "gte",
        True,
    ),
}

_FIXED_PREFIXES = (
    "/components/advntr/artifact_evidence/",
    "/components/advntr/flagging_rules/Polymorphic_Call/",
    "/components/advntr/settings/",
    "/components/kestrel/confidence_assignment/",
    "/components/kestrel/selection/final_filter_columns/",
    "/components/kestrel/selection/frameshift/",
    "/components/nomenclature/",
)


def _load(path: Path) -> dict[str, Any]:
    return load_strict_json_object(path.read_bytes())


def build_decision_projection() -> dict[str, object]:
    """Build the exact governed projection from current reviewed behavior.

    Returns:
        Complete component mapping, including decisions historically held in
        JSON and the named hard-coded selection/reconciliation contracts.
    """
    kestrel = _load(PACKAGE_ROOT / "scripts" / "kestrel_config.json")
    advntr = _load(PACKAGE_ROOT / "modules" / "advntr" / "advntr_config.json")
    nomenclature = _load(PACKAGE_ROOT / "scripts" / "nomenclature_config.json")
    main = _load(PACKAGE_ROOT / "config.json")
    evidence_path = PACKAGE_ROOT / "modules" / "advntr" / "advntr_artifact_evidence.json"
    evidence = _load(evidence_path)
    evidence_digest = (evidence_path.with_suffix(".sha256")).read_text(encoding="ascii").strip()

    kestrel_projection = {
        key: kestrel[key]
        for key in (
            "confidence_assignment",
            "alt_filtering",
            "motif_filtering",
            "flagging_rules",
            "artifact_flags",
            "duplicate_flagging",
        )
    }
    kestrel_projection["selection"] = {
        "frameshift": {
            "modulus": 3,
            "insertion_remainder": 1,
            "deletion_remainder": 2,
        },
        "final_filter_columns": [
            "is_frameshift",
            "is_valid_frameshift",
            "depth_confidence_pass",
            "alt_filter_pass",
            "motif_filter_pass",
            "flag_filter_pass",
        ],
        "confidence_priority": {
            "High_Precision*": 3,
            "High_Precision": 2,
            "Low_Precision": 1,
            "Negative": 0,
        },
        "unflagged_value": "Not flagged",
        "sort_order": [
            {"column": "confidence_priority", "ascending": False},
            {"column": "is_unflagged", "ascending": False},
            {"column": "Depth_Score", "ascending": False},
            {"column": "haplo_count", "ascending": False},
            {"column": "POS", "ascending": True},
        ],
    }

    advntr_settings = advntr["advntr_settings"]
    assert isinstance(advntr_settings, dict)
    active_entries = [entry for entry in evidence["entries"] if entry["active"]]
    return {
        "kestrel": kestrel_projection,
        "advntr": {
            "settings": {
                key: advntr_settings[key] for key in ("max_frameshift", "frameshift_multiplier", "vid", "output_format")
            },
            "flagging_rules": advntr["flagging_rules"],
            "artifact_evidence": {
                "digest": evidence_digest,
                "assertion": evidence["assertion"],
                "exact_state_match_required": True,
                "active_states": [entry["state"] for entry in active_entries],
                "active_statuses": [entry["status"] for entry in active_entries],
                "matching_disposition": "identity-insufficient",
                "nonmatching_disposition": "admissible",
            },
        },
        "shark": {},
        "nomenclature": {
            "canonical_unit": nomenclature["canonical_unit"],
            "unit_length": nomenclature["unit_length"],
            "motifs": nomenclature["motifs"],
            "advntr": {
                "mappable_repeat_units": nomenclature["advntr"]["mappable_repeat_units"],
                "rotation_offset": nomenclature["advntr"]["rotation_offset"],
            },
            "sources": {"caller_of": nomenclature["sources"]["caller_of"]},
            "thresholds": {
                "bam_flank": nomenclature["thresholds"]["bam_flank"],
                "bam_thin_haplotype_record_support": nomenclature["thresholds"]["bam_thin_haplotype_record_support"],
            },
            "identity_reconciliation": {
                "kestrel_min_alternate_kmer_path_depth": nomenclature["thresholds"]["min_support_for_high_confidence"],
                "advntr_min_sequencing_read_support": nomenclature["thresholds"]["min_support_for_high_confidence"],
                "source_evidence_units": {
                    "kestrel_vcf": "alternate-kmer-path-depth",
                    "kestrel_bam": "resolved-haplotype-records",
                    "advntr": "sequencing-reads",
                },
                "admissible_disposition": "admissible",
                "independent_callers_required": 2,
                "tier_a_blocking_flags": [
                    "motif-context-diverges",
                    "sequence-undetermined",
                    "caller-disagreement",
                ],
            },
            "known_variants": nomenclature["known_variants"],
        },
        "cross_match": {
            "match_rule": main["cross_match"]["match_rule"],
            "required_advntr_evidence_disposition": "admissible",
        },
        "dominance": {
            "enabled": False,
            "minimum_record_count_margin": 0,
            "minimum_record_share": 0.0,
            "minimum_record_share_margin": 0.0,
            "xd_veto": "disabled",
            "abstain_on_inadmissible_advntr": False,
        },
    }


def _validation_class(pointer: str) -> ValidationClass:
    if pointer.startswith("/components/dominance/"):
        return ValidationClass.GENERATED_MUTABLE
    if pointer in _CRITICAL_NUMERIC_METADATA or pointer.startswith(_FIXED_PREFIXES):
        return ValidationClass.FIXED_SAFETY
    if pointer == "/components/kestrel/alt_filtering/gg_depth_score_threshold":
        return ValidationClass.FIXED_SAFETY
    if pointer in {
        "/components/cross_match/required_advntr_evidence_disposition",
        "/components/kestrel/flagging_rules/Low_Depth_Conserved_Motifs/all/1/operator",
        "/components/shark",
    }:
        return ValidationClass.FIXED_SAFETY
    return ValidationClass.EXPLICIT_CUSTOM


def _numeric_metadata(pointer: str) -> tuple[str, str, bool]:
    critical = _CRITICAL_NUMERIC_METADATA.get(pointer)
    if critical is not None:
        return critical
    if pointer.endswith("/Low_Depth_Conserved_Motifs/all/0/right/literal"):
        return "depth-score-ratio", "lt", False
    if pointer.endswith("/Low_Coverage/all/0/right/literal"):
        return "sequencing-reads", "lt", False
    if pointer.endswith("/position_threshold"):
        return "motif-pair-position", "split-at", True
    if pointer.endswith("/unit_length"):
        return "base-pairs", "eq", True
    if pointer.endswith("/rotation_offset"):
        return "base-pairs", "eq", True
    if pointer.endswith("/independent_callers_required"):
        return "independent-callers", "gte", True
    if "/confidence_priority/" in pointer:
        return "ordinal-rank", "eq", True
    if "/selection/frameshift/" in pointer:
        return "bases-modulo-frame", "eq", True
    if pointer.endswith("/max_frameshift"):
        return "accepted-series-terms", "eq", True
    if pointer.endswith("/frameshift_multiplier"):
        return "bases", "eq", True
    if pointer.endswith("/vid"):
        return "advntr-model-identifier", "eq", True
    if pointer.endswith("/minimum_record_count_margin"):
        return "resolved-haplotype-records", "gte", True
    if pointer.endswith(("/minimum_record_share", "/minimum_record_share_margin")):
        return "resolved-haplotype-record-fraction", "gte", True
    raise ValueError(f"numeric decision leaf lacks reviewed semantics: {pointer}")


def build_decision_profile(projection: dict[str, object]) -> dict[str, object]:
    """Build the packaged profile inventory from a complete projection.

    Args:
        projection: Complete reviewed decision projection.

    Returns:
        Canonical-profile-compatible object with one class per leaf.
    """
    inventory: dict[str, dict[str, object]] = {}
    for pointer, value in flatten_decision_projection(projection).items():
        field: dict[str, object] = {"class": _validation_class(pointer).value, "value": value}
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            unit, comparator, inclusive = _numeric_metadata(pointer)
            field.update({"unit": unit, "comparator": comparator, "inclusive": inclusive})
        inventory[pointer] = field
    return {
        "schema_version": 1,
        "profile_id": "vntyper-packaged-default",
        "profile_revision": "1",
        "profile_kind": "packaged",
        "generated_metadata": None,
        "inventory": inventory,
    }


def render(output_root: Path = OUTPUT_ROOT) -> None:
    """Write the canonical projection, profile, and digest sidecar.

    Args:
        output_root: Destination directory, overridable for round-trip tests.
    """
    projection = build_decision_projection()
    profile = build_decision_profile(projection)
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "decision_projection.json").write_bytes(canonical_json_bytes(projection))
    (output_root / "decision_profile.json").write_bytes(canonical_json_bytes(profile))
    (output_root / "decision_profile.sha256").write_text(f"{canonical_sha256(profile)}\n", encoding="ascii")


def main() -> None:
    """Render checked-in artifacts or verify that rerendering is byte-identical."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if checked-in artifacts differ from a rerender")
    args = parser.parse_args()
    if not args.check:
        render()
        return
    expected = {
        path.name: path.read_bytes()
        for path in OUTPUT_ROOT.iterdir()
        if path.name in {"decision_projection.json", "decision_profile.json", "decision_profile.sha256"}
    }
    with tempfile.TemporaryDirectory() as temporary:
        destination = Path(temporary)
        render(destination)
        actual = {path.name: path.read_bytes() for path in destination.iterdir()}
    if actual != expected:
        raise SystemExit("decision profile artifacts are stale; run scripts/render_decision_profile.py")


if __name__ == "__main__":
    main()
