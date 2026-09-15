"""Strict synthetic contracts for length feature provenance."""

from __future__ import annotations

from dataclasses import replace

import pytest

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_feature_provenance import (
    LENGTH_FEATURE_CONTEXT_SCHEMA_VERSION,
    LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION,
    DenominatorQc,
    LengthFeatureProvenance,
    decode_length_feature_context,
    decode_length_feature_provenance,
    encode_length_feature_context,
    encode_length_feature_provenance,
)

pytestmark = pytest.mark.unit


def context_raw(**changes: object) -> dict[str, object]:
    """Return one complete synthetic measurement context document."""
    value: dict[str, object] = {
        "schema_version": LENGTH_FEATURE_CONTEXT_SCHEMA_VERSION,
        "manifest_key": "synthetic-member-1",
        "input_sha256": "1" * 64,
        "assembly": "synthetic-build-v1",
        "assay_class": "synthetic-short-read",
        "input_scope": "regional",
        "original_contig": "synthetic-contig",
        "reference_fasta_sha256": "2" * 64,
        "annotation_sha256": "3" * 64,
        "aligner": {
            "name": "synthetic-aligner",
            "version": "1.0",
            "arguments_sha256": "4" * 64,
            "primary_secondary_marking": "primary-only",
        },
        "preprocessing_id": "synthetic-preprocessing-v1",
        "counting_policy": {
            "policy_id": "primary-mapq0-baseq0-overlap-count-v1",
            "samtools_revision": "samtools-synthetic-1.0",
            "minimum_mapping_quality": 0,
            "minimum_base_quality": 0,
            "excluded_alignment_flags": ["UNMAP", "SECONDARY", "QCFAIL", "DUP"],
            "supplementary_alignment_policy": "included-unless-excluded-by-another-flag",
            "overlap_policy": "count-overlapping-mates-independently",
            "base_counting_policy": "one-per-aligned-covered-base",
            "zero_coverage_policy": "emit-zero-for-every-queried-position",
            "queried_intervals": [{"start": 0, "end": 6}],
        },
    }
    value.update(changes)
    return value


def provenance_raw(**changes: object) -> dict[str, object]:
    """Return one complete synthetic bound provenance document."""
    value: dict[str, object] = {
        "schema_version": LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION,
        "measurement_context": context_raw(),
        "denominator_qc": {
            "evidence_kind": "read-pair-identity-qc-proxy",
            "supporting_fragment_counts": {"CORE": 2, "INVARIANT": 3, "ARRAY": 4, "BOTH_FLANKS": 2},
        },
    }
    value.update(changes)
    return value


def _context_object(name: str) -> dict[str, object]:
    value = context_raw()[name]
    assert isinstance(value, dict)
    return value


def test_context_encoding_binds_complete_measurement_identity_and_policy() -> None:
    context = decode_length_feature_context(context_raw())

    encoded = encode_length_feature_context(context)

    assert encoded == context_raw()
    assert context.sha256 == canonical_sha256(encoded)
    assert context.counting_policy_sha256 == canonical_sha256(encoded["counting_policy"])
    assert context.original_contig == "synthetic-contig"


@pytest.mark.parametrize(
    "mutation",
    [
        {"input_sha256": "9" * 64},
        {"original_contig": "alias-contig"},
        {"preprocessing_id": "synthetic-preprocessing-v2"},
        {"aligner": {**_context_object("aligner"), "version": "2.0"}},
        {
            "counting_policy": {
                **_context_object("counting_policy"),
                "samtools_revision": "samtools-synthetic-2.0",
            }
        },
    ],
)
def test_each_material_context_change_changes_its_digest(mutation: dict[str, object]) -> None:
    baseline = decode_length_feature_context(context_raw())
    changed = decode_length_feature_context(context_raw(**mutation))

    assert changed.sha256 != baseline.sha256


@pytest.mark.parametrize(
    ("path", "value"),
    [
        ("root", {**context_raw(), "unknown": "field"}),
        ("aligner", {**_context_object("aligner"), "unknown": "field"}),
        ("counting_policy", {**_context_object("counting_policy"), "unknown": "field"}),
    ],
)
def test_context_rejects_unknown_fields(path: str, value: object) -> None:
    raw = context_raw()
    if path == "root":
        raw = value  # type: ignore[assignment]
    else:
        raw[path] = value

    with pytest.raises(ValueError, match="fields"):
        decode_length_feature_context(raw)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("minimum_mapping_quality", True),
        ("minimum_mapping_quality", 1),
        ("minimum_base_quality", 0.0),
        ("queried_intervals", [{"start": 0, "end": 2}, {"start": 2, "end": 6}]),
        ("queried_intervals", [{"start": 0, "end": 4}, {"start": 3, "end": 6}]),
    ],
)
def test_counting_policy_rejects_boolean_thresholds_and_noncanonical_interval_unions(field: str, value: object) -> None:
    raw = context_raw()
    policy = dict(_context_object("counting_policy"))
    policy[field] = value
    raw["counting_policy"] = policy

    with pytest.raises(ValueError):
        decode_length_feature_context(raw)


def test_bound_provenance_round_trips_and_binds_both_flanks_count() -> None:
    provenance = decode_length_feature_provenance(provenance_raw())

    encoded = encode_length_feature_provenance(provenance)

    assert encoded == provenance_raw()
    assert provenance.sha256 == canonical_sha256(encoded)
    assert provenance.denominator_qc.supporting_fragment_counts["BOTH_FLANKS"] == 2


def test_provenance_rejects_unknown_or_invalid_denominator_qc() -> None:
    raw = provenance_raw()
    raw["denominator_qc"] = {
        "evidence_kind": "read-pair-identity-qc-proxy",
        "supporting_fragment_counts": {
            "CORE": 2,
            "INVARIANT": 3,
            "ARRAY": 4,
            "BOTH_FLANKS": True,
            "UNKNOWN": 1,
        },
    }
    with pytest.raises(ValueError):
        decode_length_feature_provenance(raw)


def test_depth_only_provenance_requires_all_denominator_counts_to_be_null() -> None:
    raw = provenance_raw()
    raw["denominator_qc"] = {
        "evidence_kind": "unavailable",
        "supporting_fragment_counts": {"CORE": 2, "INVARIANT": None, "ARRAY": None, "BOTH_FLANKS": None},
    }
    with pytest.raises(ValueError, match="unavailable"):
        decode_length_feature_provenance(raw)


def test_public_encoders_revalidate_replace_created_values_and_hashes() -> None:
    context = decode_length_feature_context(context_raw())
    with pytest.raises(ValueError, match="digest"):
        encode_length_feature_context(replace(context, sha256="f" * 64))
    with pytest.raises(ValueError, match="typed policy"):
        encode_length_feature_context(replace(context, aligner=True))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="digest"):
        encode_length_feature_context(replace(context, aligner=replace(context.aligner, version="2.0")))

    provenance = decode_length_feature_provenance(provenance_raw())
    bad_qc = DenominatorQc(
        evidence_kind="read-pair-identity-qc-proxy",
        supporting_fragment_counts={"CORE": 2, "INVARIANT": 3, "ARRAY": 4, "BOTH_FLANKS": 99},
    )
    with pytest.raises(ValueError, match="digest"):
        encode_length_feature_provenance(replace(provenance, denominator_qc=bad_qc))


def test_provenance_rejects_a_directly_constructed_wrong_typed_object() -> None:
    provenance = decode_length_feature_provenance(provenance_raw())
    malformed = LengthFeatureProvenance(
        schema_version=provenance.schema_version,
        measurement_context=provenance.measurement_context,
        denominator_qc=provenance.denominator_qc,
        sha256=True,  # type: ignore[arg-type]
    )
    with pytest.raises(ValueError, match="digest"):
        encode_length_feature_provenance(malformed)
