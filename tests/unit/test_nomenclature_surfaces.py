"""The nomenclature columns on every output surface.

The five columns are added in one pass so the result TSVs, the summary JSON built
from those same frames, and the HTML report all inherit them together.

Research use only.
"""

from pathlib import Path
from unittest import mock

import pandas as pd
import pysam
import pytest

from tests.builders import STAGE_COLUMNS, kestrel_stage_frame
from vntyper.scripts import nomenclature
from vntyper.scripts.cohort_tables import ADVNTR_DISPLAY_COLUMNS as COHORT_ADVNTR
from vntyper.scripts.cohort_tables import KESTREL_DISPLAY_COLUMNS as COHORT_KESTREL
from vntyper.scripts.nomenclature import Nomenclature
from vntyper.scripts.nomenclature_annotate import (
    NOMENCLATURE_COLUMNS,
    annotate_advntr_frame,
    annotate_kestrel_frame,
    reconcile_caller_outputs,
)
from vntyper.scripts.report_formatting import ADVNTR_DISPLAY_COLUMNS, KESTREL_DISPLAY_COLUMNS
from vntyper.scripts.summary import parse_tsv

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Stage contract
# ---------------------------------------------------------------------------


def test_the_named_stage_is_final_plus_the_nomenclature_columns() -> None:
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


def test_a_below_tier_a_row_still_shows_its_name() -> None:
    """The name is emitted whenever one exists; the tier column carries the caveat.

    One caller alone cannot reach tier A, but withholding the name it computed
    loses information the reader can weigh against the tier themselves.
    """
    frame = pd.DataFrame([{"Motifs": "X-X", "POS": 67, "REF": "G", "ALT": "GG"}])
    out = annotate_kestrel_frame(frame)
    assert out.loc[0, "Nomenclature"] == "59dupC"
    assert out.loc[0, "Nomenclature_Tier"] == "B", "the lower confidence must still be reported"


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


def _write_disagreeing_outputs(tmp_path) -> tuple[Path, Path]:
    """Kestrel and adVNTR naming different alleles at the same locus.

    The real ``insG`` shape: Kestrel's VCF places the insertion one base 3' of where
    adVNTR places it.
    """
    kestrel = tmp_path / "kestrel_result.tsv"
    kestrel.write_text(
        "## VNtyper Kestrel result\n"
        "Motifs\tPOS\tREF\tALT\tConfidence\tEstimated_Depth_AlternateVariant\t"
        "Motif_fasta\tPOS_fasta\tNomenclature\tNomenclature_Tier\tNomenclature_Flags\t"
        "Ambiguity_Interval\tRepeat_Form\tNomenclature_Note\t"
        "Nomenclature_Kestrel\tNomenclature_adVNTR\n"
        "X-X\t61\tT\tTC\tHigh_Precision\t40\tX-X\t61\t59_60insG\tB\t\t\t\t"
        "representation of the caller's call, not a described variant; requires validation\t"
        "59_60insG\t\n"
    )
    advntr = tmp_path / "output_adVNTR_result.tsv"
    advntr.write_text(
        "VID\tVariant\tNumberOfSupportingReads\tNomenclature\tNomenclature_Tier\t"
        "Nomenclature_Flags\tAmbiguity_Interval\tRepeat_Form\tNomenclature_Note\t"
        "Nomenclature_Kestrel\tNomenclature_adVNTR\n"
        "25561\tI23_2_C_LEN1\t40\t58_59insG\tB\t\t\t\t\t\t58_59insG\n"
    )
    return kestrel, advntr


def test_each_caller_keeps_its_own_name_when_they_disagree(tmp_path) -> None:
    """Overwriting both files with one verdict destroyed the per-caller information.

    Which caller said what is the evidence a reader needs in order to weigh a
    disagreement at all; collapsing it to a single cell throws that away.
    """
    kestrel, advntr = _write_disagreeing_outputs(tmp_path)
    assert reconcile_caller_outputs(kestrel, advntr) is True

    for path, sep_kwargs in ((kestrel, {"comment": "#"}), (advntr, {})):
        written = pd.read_csv(path, sep="\t", dtype=str, **sep_kwargs)
        assert written.loc[0, "Nomenclature_Kestrel"] == "59_60insG"
        assert written.loc[0, "Nomenclature_adVNTR"] == "58_59insG"


def test_a_disagreement_is_declared_on_both_files(tmp_path) -> None:
    kestrel, advntr = _write_disagreeing_outputs(tmp_path)
    reconcile_caller_outputs(kestrel, advntr)

    for path, sep_kwargs in ((kestrel, {"comment": "#"}), (advntr, {})):
        written = pd.read_csv(path, sep="\t", dtype=str, **sep_kwargs)
        assert "caller-disagreement" in written.loc[0, "Nomenclature_Flags"]
        assert written.loc[0, "Nomenclature_Tier"] != "A"


def test_the_kestrel_file_is_complete_without_advntr_having_run(tmp_path) -> None:
    """adVNTR is optional; Kestrel is not.

    The primary surface must carry its own name and a full set of columns whether or
    not the optional module ever ran, so nothing downstream depends on it.
    """
    frame = kestrel_stage_frame("final")
    named = annotate_kestrel_frame(frame)
    for column in NOMENCLATURE_COLUMNS:
        assert column in named.columns
    assert named.loc[0, "Nomenclature_Kestrel"] == named.loc[0, "Nomenclature"]
    assert named.loc[0, "Nomenclature_adVNTR"] == "", "adVNTR did not run; it reported nothing"


def test_the_advntr_frame_records_its_own_name_as_the_advntr_column() -> None:
    frame = pd.DataFrame([{"VID": "25561", "Variant": "I22_2_G_LEN1", "NumberOfSupportingReads": "24"}])
    named = annotate_advntr_frame(frame)
    assert named.loc[0, "Nomenclature_adVNTR"] == "59dupC"
    assert named.loc[0, "Nomenclature_Kestrel"] == ""


def _write_pair_bam(path: Path, inserted_at: int, base: str, reads: int) -> None:
    """Write an indexed ``output.bam`` over the real ``X-X`` pair reference.

    Args:
        path: Destination ``output.bam``.
        inserted_at: 0-based reference position the insertion sits before.
        base: The inserted plus-strand base.
        reads: How many reads carry it.
    """
    pair = nomenclature.pair_sequence("X-X")
    assert pair is not None
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "X-X", "LN": len(pair)}]}
    trailing = len(pair) - inserted_at
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for index in range(reads):
            record = pysam.AlignedSegment(handle.header)
            record.query_name = f"read{index}"
            record.reference_id = 0
            record.reference_start = 0
            record.mapping_quality = 255
            record.cigarstring = f"{inserted_at}=1I{trailing}="
            record.query_sequence = pair[:inserted_at] + base + pair[inserted_at:]
            handle.write(record)
    pysam.index(str(path))  # type: ignore[attr-defined]


def test_the_cross_caller_stage_consults_the_reads(tmp_path) -> None:
    """The reads must reach the only step that can see both callers.

    This is the production path for the whole vote. The Kestrel VCF places the
    insertion at ``59_60insG``; adVNTR and the reads independently say ``58_59insG``,
    and two sources from two callers outvote one. Without the BAM the cross-caller
    step sees only the two callers disagreeing and the wrong VCF placement stands.
    """
    kestrel, advntr = _write_disagreeing_outputs(tmp_path)
    _write_pair_bam(tmp_path / "output.bam", inserted_at=62, base="C", reads=4)

    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature"] == "58_59insG", "the reads corroborated adVNTR; the VCF was outvoted"
    assert written.loc[0, "Nomenclature_Kestrel"] == "59_60insG", "Kestrel's own record is still reported"
    assert written.loc[0, "Nomenclature_adVNTR"] == "58_59insG"


def test_the_reads_alone_do_not_outvote_the_vcf(tmp_path) -> None:
    """Same reads, but adVNTR agreeing with Kestrel instead.

    The reads are Kestrel's own alignment, so on their own they are not a second
    caller and must not overturn the VCF -- the guard against the earlier policy that
    cost `dupA` 6 correct calls out of 10.
    """
    kestrel, advntr = _write_disagreeing_outputs(tmp_path)
    # I22_2_C_LEN1 names 59_60insG, the same allele the Kestrel VCF names.
    advntr.write_text(advntr.read_text().replace("I23_2_C_LEN1", "I22_2_C_LEN1"))
    _write_pair_bam(tmp_path / "output.bam", inserted_at=62, base="C", reads=4)

    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature"] == "59_60insG", "one caller's reads may not veto its own VCF"


def test_an_annotated_kestrel_result_survives_the_summary_parser(tmp_path) -> None:
    """The written result must round-trip into ``pipeline_summary.json``.

    The nomenclature columns are nullable by contract, and the last of them is empty
    on every run where adVNTR did not contribute. A row ending in an empty field is
    one field short once its line is stripped, so the summary parser discarded it as
    ragged -- and the pipeline then raised "Kestrel genotyping results not found" on
    a run that had in fact succeeded. Nothing about the row was wrong, which is why
    only an end-to-end run found it.
    """
    named = annotate_kestrel_frame(kestrel_stage_frame("final"))
    assert named.loc[0, "Nomenclature_adVNTR"] == "", "the last column must be the empty case"

    path = tmp_path / "kestrel_result.tsv"
    with path.open("w") as handle:
        handle.write("## VNtyper Kestrel result\n")
        named.to_csv(handle, sep="\t", index=False)

    parsed = parse_tsv(str(path))
    assert len(parsed["data"]) == len(named), "every written row must reach the summary"
    assert parsed["data"][0]["Nomenclature_Kestrel"] == named.loc[0, "Nomenclature"]


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


# ---------------------------------------------------------------------------
# The config is the source of truth, and it must ship
# ---------------------------------------------------------------------------


def test_the_nomenclature_config_ships_with_the_package() -> None:
    """Reference data lives in config, so a config that is not packaged is a bug.

    `MANIFEST.in` and `pyproject.toml` both list it, the way every other config in
    this package is listed; this asserts the file they name is actually there.
    """
    config = Path(nomenclature.__file__).with_name("nomenclature_config.json")
    assert config.is_file(), f"{config.name} must sit beside the module that loads it"

    root = Path(__file__).resolve().parents[2]
    assert "nomenclature_config.json" in (root / "MANIFEST.in").read_text()
    assert "nomenclature_config.json" in (root / "pyproject.toml").read_text()


def test_every_configured_table_is_populated() -> None:
    """A silently empty table would make every call undetermined rather than fail."""
    assert len(nomenclature.MOTIFS) == 34
    assert len(nomenclature.CANONICAL_UNIT) == nomenclature.UNIT_LENGTH
    assert nomenclature.MAPPABLE_RUS, "no adVNTR repeat unit is mappable"
    assert nomenclature.KNOWN_VARIANTS, "the described-variant check has nothing to check against"
    assert nomenclature.MIN_SUPPORT_FOR_TIER_A > 0


def test_no_reference_data_is_written_into_the_module() -> None:
    """Domain knowledge belongs in the config, not in the logic.

    Guards the rule directly: a future edit that pastes a motif sequence or a
    variant name back into the module fails here rather than being noticed later.
    """
    source = Path(nomenclature.__file__).read_text()
    for motif in list(nomenclature.MOTIFS.values())[:5]:
        assert motif not in source, "a motif sequence was written into the module"
    for name in nomenclature.KNOWN_VARIANTS:
        assert f'"{name}"' not in source, f"the variant {name} was written into the module"


def test_a_variant_can_be_added_by_editing_config_alone() -> None:
    """The described-variant list is data: extending it must not need a code change."""
    extended = dict(nomenclature.KNOWN_VARIANTS)
    extended["99dupT"] = "Fictional et al. 2099"
    call = Nomenclature(
        name="99dupT",
        event="duplication",
        unit="X",
        tier="B",
        flags=(),
        ambiguity=None,
        repeat_form=None,
        net_length=1,
        source="kestrel_vcf",
    )
    with mock.patch.object(nomenclature, "KNOWN_VARIANTS", extended):
        assert "Fictional" in nomenclature.confidence_note(call)
