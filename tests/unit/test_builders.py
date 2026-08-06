"""Tests for the shared test-data builders.

The builders are an oracle used by every other test file, so they need their
own tests -- a silently wrong builder makes every consumer green for the wrong
reason.
"""

import pandas as pd
import pytest

pytestmark = pytest.mark.unit

from tests.builders import (  # noqa: E402
    STAGE_COLUMNS,
    bam_header,
    kestrel_config,
    kestrel_df,
    kestrel_row,
    kestrel_stage_frame,
)


def test_kestrel_row_defaults_are_a_coherent_call() -> None:
    row = kestrel_row()
    assert (
        row["Sample"] == f"Del:{row['Estimated_Depth_AlternateVariant']}:{row['Estimated_Depth_Variant_ActiveRegion']}"
    )
    assert len(row["ALT"]) != len(row["REF"]), "default row must be a frameshift"


def test_kestrel_row_overrides_flow_into_the_sample_field() -> None:
    row = kestrel_row(depth_alt=7, depth_region=700)
    assert row["Estimated_Depth_AlternateVariant"] == 7
    assert row["Sample"] == "Del:7:700"


def test_kestrel_df_preserves_row_order_and_dtypes() -> None:
    frame = kestrel_df(kestrel_row(pos=60), kestrel_row(pos=67))
    assert list(frame["POS"]) == [60, 67]
    assert frame["POS"].dtype.kind == "i"


def test_stage_frames_are_cumulative() -> None:
    """Each stage adds columns; nothing is ever removed."""
    order = ["raw", "scored", "confidence", "flagged", "final"]
    for earlier, later in zip(order, order[1:], strict=False):
        assert set(STAGE_COLUMNS[earlier]).issubset(set(STAGE_COLUMNS[later])), (
            f"{later} dropped columns from {earlier}"
        )


@pytest.mark.parametrize("stage", ["raw", "scored", "confidence", "flagged", "final"])
def test_stage_frame_has_exactly_its_declared_columns(stage: str) -> None:
    frame = kestrel_stage_frame(stage)
    assert tuple(frame.columns) == STAGE_COLUMNS[stage]


def test_stage_frame_rejects_an_unknown_stage() -> None:
    with pytest.raises(ValueError, match="Unknown stage"):
        kestrel_stage_frame("nonsense")


def test_stage_frame_rows_are_distinguishable() -> None:
    """A builder returning N identical rows is a trap for dedup tests."""
    frame = kestrel_stage_frame("raw", rows=3)
    assert len(frame) == 3
    assert list(frame["POS"]) == [67, 68, 69]
    assert not frame.duplicated().any()


def test_kestrel_config_returns_the_real_values_by_default() -> None:
    conf = kestrel_config()
    assert conf["confidence_assignment"]["depth_score_thresholds"]["low"] == 0.00469


def test_kestrel_config_applies_a_dotted_override_without_mutating_the_original() -> None:
    conf = kestrel_config(**{"confidence_assignment.depth_score_thresholds.low": 0.5})
    assert conf["confidence_assignment"]["depth_score_thresholds"]["low"] == 0.5
    assert kestrel_config()["confidence_assignment"]["depth_score_thresholds"]["low"] == 0.00469


@pytest.mark.parametrize(
    "convention,expected_first",
    [("ucsc", "chr1"), ("ensembl", "1"), ("ncbi", "NC_000001.10")],
)
def test_bam_header_uses_the_requested_naming_convention(convention: str, expected_first: str) -> None:
    header = bam_header(convention=convention, assembly="GRCh37")
    assert f"SN:{expected_first}\t" in header


@pytest.mark.parametrize("assembly,length", [("GRCh37", 249250621), ("GRCh38", 248956422)])
def test_bam_header_uses_the_correct_chr1_length(assembly: str, length: int) -> None:
    header = bam_header(assembly=assembly)
    assert f"LN:{length}" in header


def test_bam_header_rejects_an_unknown_assembly_and_convention() -> None:
    with pytest.raises(ValueError, match="Unknown assembly"):
        bam_header(assembly="hg19")
    with pytest.raises(ValueError, match="Unknown convention"):
        bam_header(convention="nonsense")


def test_the_scored_stage_frame_feeds_the_real_confidence_assignment() -> None:
    """STAGE_COLUMNS['scored'] must be the real input contract, not a guess.

    Every Wave 2 Kestrel test is built on this frame; if the column set is wrong
    the production function raises here rather than in six other files.
    """
    from vntyper.scripts.confidence_assignment import calculate_depth_score_and_assign_confidence

    out = calculate_depth_score_and_assign_confidence(kestrel_stage_frame("scored"), kestrel_config())
    assert set(STAGE_COLUMNS["confidence"]) - set(STAGE_COLUMNS["scored"]) <= set(out.columns)
    assert isinstance(out["Depth_Score"].iloc[0], float)
    assert out["Confidence"].iloc[0] == "High_Precision*"


def test_the_final_stage_frame_feeds_the_real_final_filter(tmp_path) -> None:
    """STAGE_COLUMNS['final'] must carry every boolean filter_final_dataframe ANDs."""
    from vntyper.scripts.kestrel_genotyping import filter_final_dataframe

    out = filter_final_dataframe(kestrel_stage_frame("final"), str(tmp_path))
    assert isinstance(out, pd.DataFrame)
    assert len(out) == 1, "the default row is a passing call and must survive the final filter"
