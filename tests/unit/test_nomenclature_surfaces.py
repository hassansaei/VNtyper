"""The nomenclature columns on every output surface.

The five columns are added in one pass so the result TSVs, the summary JSON built
from those same frames, and the HTML report all inherit them together.

Research use only.
"""

from pathlib import Path

import pandas as pd
import pytest

from tests.builders import STAGE_COLUMNS, kestrel_stage_frame
from vntyper.scripts.cohort_tables import ADVNTR_DISPLAY_COLUMNS as COHORT_ADVNTR
from vntyper.scripts.cohort_tables import KESTREL_DISPLAY_COLUMNS as COHORT_KESTREL
from vntyper.scripts.nomenclature_annotate import (
    NOMENCLATURE_COLUMNS,
    annotate_advntr_frame,
    annotate_kestrel_frame,
    reconcile_caller_outputs,
)
from vntyper.scripts.report_formatting import ADVNTR_DISPLAY_COLUMNS, KESTREL_DISPLAY_COLUMNS

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Stage contract
# ---------------------------------------------------------------------------


def test_the_named_stage_is_final_plus_the_five_columns() -> None:
    added = tuple(column for column in STAGE_COLUMNS["named"] if column not in STAGE_COLUMNS["final"])
    assert added == NOMENCLATURE_COLUMNS


def test_the_named_stage_keeps_the_subset_invariant() -> None:
    assert set(STAGE_COLUMNS["final"]).issubset(set(STAGE_COLUMNS["named"]))


def test_the_named_stage_builder_produces_real_names() -> None:
    frame = kestrel_stage_frame("named")
    assert tuple(frame.columns) == STAGE_COLUMNS["named"]
    assert frame.loc[0, "Nomenclature"], "the builder must produce a populated cell"


# ---------------------------------------------------------------------------
# Every surface carries the name and its tier
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "columns",
    [KESTREL_DISPLAY_COLUMNS, ADVNTR_DISPLAY_COLUMNS, COHORT_KESTREL, COHORT_ADVNTR],
    ids=["report-kestrel", "report-advntr", "cohort-kestrel", "cohort-advntr"],
)
def test_every_display_surface_carries_the_name_and_tier(columns) -> None:
    assert "Nomenclature" in columns
    assert "Nomenclature_Tier" in columns


# ---------------------------------------------------------------------------
# Negative runs
# ---------------------------------------------------------------------------


def test_a_negative_kestrel_run_keeps_its_ten_column_schema() -> None:
    """A negative run has no variant, so there is no name and no new columns.

    Its TSV is a separate 10-column schema whose first column is singular ``Motif``.
    Padding it with five permanently empty columns would widen a schema that
    deliberately differs, and the fields are nullable precisely so that they are
    never padded.
    """
    negative = pd.DataFrame(
        [
            {
                "Motif": "None",
                "Variant": "None",
                "POS": "None",
                "REF": "None",
                "ALT": "None",
                "Motif_sequence": "None",
                "Estimated_Depth_AlternateVariant": "None",
                "Estimated_Depth_Variant_ActiveRegion": "None",
                "Depth_Score": "None",
                "Confidence": "Negative",
            }
        ]
    )
    out = annotate_kestrel_frame(negative)
    assert list(out.columns) == list(negative.columns)
    assert len(out.columns) == 10


def test_a_negative_advntr_row_leaves_the_five_empty_not_placeheld() -> None:
    negative = pd.DataFrame([{"VID": "Negative", "Variant": "Not applicable", "NumberOfSupportingReads": "0"}])
    out = annotate_advntr_frame(negative)
    for column in NOMENCLATURE_COLUMNS:
        assert out.loc[0, column] == ""


def test_an_empty_frame_is_returned_untouched() -> None:
    empty = pd.DataFrame()
    assert annotate_kestrel_frame(empty).empty
    assert annotate_advntr_frame(empty).empty


# ---------------------------------------------------------------------------
# Nullability
# ---------------------------------------------------------------------------


def test_the_nullable_fields_are_empty_never_a_placeholder() -> None:
    """An unshiftable insertion has no ambiguity window and no tract."""
    frame = pd.DataFrame([{"Motifs": "X-X", "POS": 61, "REF": "T", "ALT": "TC"}])
    out = annotate_kestrel_frame(frame)
    assert out.loc[0, "Ambiguity_Interval"] == ""
    assert out.loc[0, "Repeat_Form"] == ""
    assert "N/A" not in set(out.iloc[0].astype(str))
    assert "Not applicable" not in set(out.iloc[0].astype(str))


def test_a_populated_window_and_tract_are_written() -> None:
    frame = pd.DataFrame([{"Motifs": "X-X", "POS": 67, "REF": "G", "ALT": "GG"}])
    out = annotate_kestrel_frame(frame)
    assert out.loc[0, "Ambiguity_Interval"] == "53_59"
    assert out.loc[0, "Repeat_Form"] == "53C[7]>53C[8]"


def test_a_below_tier_a_row_never_shows_a_bare_number() -> None:
    """One caller alone cannot reach tier A, so the cell must not read `59dupC`."""
    frame = pd.DataFrame([{"Motifs": "X-X", "POS": 67, "REF": "G", "ALT": "GG"}])
    out = annotate_kestrel_frame(frame)
    assert out.loc[0, "Nomenclature_Tier"] != "A"
    assert out.loc[0, "Nomenclature"] != "59dupC"


# ---------------------------------------------------------------------------
# Cross-caller reconciliation, on the files production actually writes
# ---------------------------------------------------------------------------


def _write_caller_outputs(tmp_path, advntr_state: str, support: str) -> tuple[Path, Path]:
    kestrel = tmp_path / "kestrel_result.tsv"
    kestrel.write_text(
        "## VNtyper Kestrel result\n"
        "Motifs\tPOS\tREF\tALT\tConfidence\tNomenclature\tNomenclature_Tier\t"
        "Nomenclature_Flags\tAmbiguity_Interval\tRepeat_Form\n"
        "X-X\t67\tG\tGG\tHigh_Precision\tduplication, position-ambiguous\tB\t\t53_59\t53C[7]>53C[8]\n"
    )
    advntr = tmp_path / "output_adVNTR_result.tsv"
    advntr.write_text(
        "VID\tVariant\tNumberOfSupportingReads\tNomenclature\tNomenclature_Tier\t"
        "Nomenclature_Flags\tAmbiguity_Interval\tRepeat_Form\n"
        f"25561\t{advntr_state}\t{support}\tduplication, position-ambiguous\tB\t\t53_59\t53C[7]>53C[8]\n"
    )
    return kestrel, advntr


def test_two_agreeing_callers_reach_tier_a_in_the_written_files(tmp_path) -> None:
    """Production must be able to produce a tier-A name, not just the API.

    Each caller names its own rows as it writes them, so neither stage can see the
    other; without a reconciliation step over the written results, no real run could
    ever emit tier A however well the two callers agreed.
    """
    kestrel, advntr = _write_caller_outputs(tmp_path, "I22_2_G_LEN1", "24")
    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature"] == "59dupC"
    assert written.loc[0, "Nomenclature_Tier"] == "A"

    mirrored = pd.read_csv(advntr, sep="\t", dtype=str)
    assert mirrored.loc[0, "Nomenclature"] == "59dupC", "both files must carry one verdict"


def test_disagreeing_callers_do_not_reach_tier_a_in_the_written_files(tmp_path) -> None:
    kestrel, advntr = _write_caller_outputs(tmp_path, "I23_2_C_LEN1", "24")
    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature_Tier"] != "A"
    assert written.loc[0, "Nomenclature"] != "59dupC"


def test_thin_advntr_support_does_not_reach_tier_a_in_the_written_files(tmp_path) -> None:
    kestrel, advntr = _write_caller_outputs(tmp_path, "I22_2_G_LEN1", "1")
    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature_Tier"] != "A"


def test_reconciliation_preserves_the_kestrel_header_lines(tmp_path) -> None:
    kestrel, advntr = _write_caller_outputs(tmp_path, "I22_2_G_LEN1", "24")
    reconcile_caller_outputs(kestrel, advntr)
    assert kestrel.read_text().startswith("## VNtyper Kestrel result")


def test_a_missing_advntr_file_is_not_an_error(tmp_path) -> None:
    kestrel, advntr = _write_caller_outputs(tmp_path, "I22_2_G_LEN1", "24")
    advntr.unlink()
    assert reconcile_caller_outputs(kestrel, advntr) is False
