from __future__ import annotations

import math

import pytest

from vntyper.scripts.length_features import DepthPosition
from vntyper.scripts.length_standard_features import (
    STANDARD_FEATURE_ORDER,
    STANDARD_REFERENCE_LOCUS_SHA256,
    StandardReadSummary,
    decode_standard_length_measurement,
    encode_standard_length_measurement,
    extract_standard_length_features,
)

pytestmark = pytest.mark.unit


def _depth(*, core_depth: int = 20, invariant_depth: int = 10, flank_depth: int = 5, contig: str = "chr1"):
    rows = []
    for position in range(155188296, 155192429):
        if 155188726 <= position < 155191939:
            value = core_depth
        elif 155188486 <= position < 155188726 or 155191939 <= position < 155192239:
            value = invariant_depth
        else:
            value = flank_depth
        rows.append(DepthPosition(contig, position, value, () if value == 0 else (f"fragment-{position}",)))
    return tuple(rows)


def _reads() -> StandardReadSummary:
    return StandardReadSummary(
        eligible_read_count=4,
        mapq_zero_count=1,
        mapq_sum=120,
        soft_clipped_read_count=2,
        query_gc_bases=150,
        query_bases=300,
    )


def test_extract_standard_features_has_exact_order_and_formulas() -> None:
    measured = extract_standard_length_features(
        _depth(),
        _reads(),
        assembly="GRCh38",
        contig="chr1",
        reference_locus_sha256=STANDARD_REFERENCE_LOCUS_SHA256,
    )

    assert measured.feature_order == STANDARD_FEATURE_ORDER
    assert measured.values["A"] == 2.0
    assert measured.values["F"] == pytest.approx((3213 * 20 + 268 * 10) / 3481 / 5)
    assert measured.values["core_zero_fraction"] == 0
    assert measured.values["core_depth_cv"] == 0
    assert measured.values["core_bin_cv"] == 0
    assert measured.values["core_half_log_ratio"] == 0
    assert measured.values["invariant_end_log_ratio"] == 0
    assert measured.values["flank_end_log_ratio"] == 0
    assert measured.values["log_invariant_depth"] == pytest.approx(math.log1p(10))
    assert measured.values["mapq_zero_fraction"] == 0.25
    assert measured.values["mean_mapq"] == 30
    assert measured.values["soft_clipped_read_fraction"] == 0.5
    assert measured.values["query_sequence_gc_fraction"] == 0.5
    assert measured.reasons == ()
    assert measured.qc.left_flank_supporting_fragments == 190
    assert measured.qc.right_flank_supporting_fragments == 190
    with pytest.raises(TypeError):
        measured.values["A"] = 3  # type: ignore[index]


def test_standard_measurement_round_trip_revalidates_digest() -> None:
    measured = extract_standard_length_features(
        _depth(contig="1"),
        _reads(),
        assembly="GRCh38",
        contig="1",
        reference_locus_sha256=STANDARD_REFERENCE_LOCUS_SHA256,
    )
    document = encode_standard_length_measurement(measured)
    assert decode_standard_length_measurement(document) == measured

    document["values"]["A"] = 9  # type: ignore[index]
    with pytest.raises(ValueError, match="digest"):
        decode_standard_length_measurement(document)


def test_signed_log_ratio_round_trips() -> None:
    depths = list(_depth())
    for index, row in enumerate(depths):
        if 155188486 <= row.position_zero_based < 155188726:
            depths[index] = DepthPosition(row.contig, row.position_zero_based, 5, row.supporting_fragment_ids)
    measured = extract_standard_length_features(
        depths,
        _reads(),
        assembly="GRCh38",
        contig="chr1",
        reference_locus_sha256=STANDARD_REFERENCE_LOCUS_SHA256,
    )
    assert measured.values["invariant_end_log_ratio"] < 0  # type: ignore[operator]
    assert decode_standard_length_measurement(encode_standard_length_measurement(measured)) == measured


def test_zero_denominators_and_no_reads_are_explicitly_missing() -> None:
    measured = extract_standard_length_features(
        _depth(invariant_depth=0, flank_depth=0),
        StandardReadSummary(0, 0, 0, 0, 0, 0),
        assembly="GRCh38",
        contig="chr1",
        reference_locus_sha256=STANDARD_REFERENCE_LOCUS_SHA256,
    )

    assert measured.values["A"] is None
    assert measured.values["F"] is None
    assert measured.values["mean_mapq"] is None
    assert measured.reasons == (
        "zero_invariant_mean_depth",
        "zero_combined_flank_mean_depth",
        "no_eligible_reads",
    )


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "wrong-contig"])
def test_standard_feature_positions_must_be_exact(mutation: str) -> None:
    depths = list(_depth())
    if mutation == "missing":
        depths.pop()
    elif mutation == "duplicate":
        depths.append(depths[-1])
    else:
        depths[0] = DepthPosition("chr2", depths[0].position_zero_based, depths[0].depth, ("x",))
    with pytest.raises(ValueError):
        extract_standard_length_features(
            depths,
            _reads(),
            assembly="GRCh38",
            contig="chr1",
            reference_locus_sha256=STANDARD_REFERENCE_LOCUS_SHA256,
        )


def test_standard_read_summary_rejects_booleans_and_inconsistent_counts() -> None:
    with pytest.raises(ValueError):
        StandardReadSummary(True, 0, 0, 0, 0, 0)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        StandardReadSummary(1, 2, 0, 0, 0, 1)
