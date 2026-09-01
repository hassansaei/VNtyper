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
from vntyper.scripts import nomenclature, nomenclature_annotate
from vntyper.scripts.cohort_tables import ADVNTR_DISPLAY_COLUMNS as COHORT_ADVNTR
from vntyper.scripts.cohort_tables import KESTREL_DISPLAY_COLUMNS as COHORT_KESTREL
from vntyper.scripts.molecular_identity import (
    MolecularIdentity,
    make_coding_edit,
    make_molecular_identity,
    serialize_molecular_identity,
)
from vntyper.scripts.nomenclature import (
    FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
    FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,
    Nomenclature,
    from_advntr,
    from_kestrel,
)
from vntyper.scripts.nomenclature_annotate import (
    NOMENCLATURE_COLUMNS,
    annotate_advntr_frame,
    annotate_kestrel_frame,
    reconcile_caller_outputs,
)
from vntyper.scripts.nomenclature_bam import BamRescuer
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
    # The depth column is present and populated: an agreement is only as strong as its
    # weakest source, and a source with no depth at all makes it unknown rather than
    # borrowing the other caller's.
    kestrel.write_text(
        "## VNtyper Kestrel result\n"
        "Motifs\tPOS\tREF\tALT\tConfidence\tEstimated_Depth_AlternateVariant\tNomenclature\t"
        "Nomenclature_Tier\tNomenclature_Flags\tAmbiguity_Interval\tRepeat_Form\n"
        "X-X\t67\tG\tGG\tHigh_Precision\t40\tduplication, position-ambiguous\tB\t\t53_59\t53C[7]>53C[8]\n"
    )
    advntr = tmp_path / "output_adVNTR_result.tsv"
    advntr.write_text(
        "VID\tVariant\tNumberOfSupportingReads\tNomenclature\tNomenclature_Tier\t"
        "Nomenclature_Flags\tAmbiguity_Interval\tRepeat_Form\n"
        f"25561\t{advntr_state}\t{support}\tduplication, position-ambiguous\tB\t\t53_59\t53C[7]>53C[8]\n"
    )
    return kestrel, advntr


def _write_identity_aware_outputs(
    tmp_path: Path,
    kestrel_identity: MolecularIdentity,
    *,
    advntr_flag: str = "Not flagged",
) -> tuple[Path, Path]:
    raw_key = '{"source":"kestrel","values":["X-X",67,"G","GG"]}'
    kestrel_row = {
        "Motifs": "X-X",
        "POS": "67",
        "REF": "G",
        "ALT": "GG",
        "Confidence": "High_Precision",
        "Estimated_Depth_AlternateVariant": "40",
        "Nomenclature": "59dupC",
        "Nomenclature_Tier": "B",
        "Nomenclature_Flags": "",
        "Ambiguity_Interval": "53_59",
        "Repeat_Form": "53C[7]>53C[8]",
        "__Identity_Raw_Representation_Key": raw_key,
        "__Identity_Molecular_Identity": serialize_molecular_identity(kestrel_identity),
        "__Identity_Translation_Status": "resolved",
        "__Identity_Translation_Failure": "absent",
        "__Identity_Context_Diverges": "false",
        "__Identity_Observation_Ordinal": "0",
        "__Identity_Selected_Raw_Representation_Key": raw_key,
        "__Identity_Equivalent_Representation_Count": "1",
        "__Identity_Hypothesis_Count": "1",
        "__Identity_Group_Blocking_Gates": "[]",
        "__Identity_Group_Flags": "[]",
        "__Identity_Selected_Observation_Ordinal": "0",
        "__Identity_Group_Context_Diverges": "false",
    }
    kestrel = tmp_path / "kestrel_result.tsv"
    pd.DataFrame([kestrel_row]).to_csv(kestrel, sep="\t", index=False)
    advntr = tmp_path / "output_adVNTR_result.tsv"
    pd.DataFrame(
        [
            {
                "VID": "25561",
                "Variant": "I22_2_G_LEN1",
                "NumberOfSupportingReads": "40",
                "MeanCoverage": "100",
                "RU": "2",
                "POS": "22",
                "Flag": advntr_flag,
                "Nomenclature": "59dupC",
                "Nomenclature_Tier": "B",
                "Nomenclature_Flags": "",
                "Ambiguity_Interval": "53_59",
                "Repeat_Form": "53C[7]>53C[8]",
            }
        ]
    ).to_csv(advntr, sep="\t", index=False)
    return kestrel, advntr


def test_production_reconciliation_uses_persisted_kestrel_identity_not_equal_display(tmp_path) -> None:
    dup_a = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))
    kestrel, advntr = _write_identity_aware_outputs(tmp_path, dup_a)

    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", dtype=str)
    assert written.loc[0, "Nomenclature"] == "59dupC"
    assert written.loc[0, "Nomenclature_Tier"] == "B"
    assert "caller-disagreement" in written.loc[0, "Nomenclature_Flags"]


def test_pr_a_polymorphic_advntr_row_remains_admissible_for_identity_agreement(tmp_path) -> None:
    dup_c = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
    kestrel, advntr = _write_identity_aware_outputs(tmp_path, dup_c, advntr_flag="Polymorphic_Call")

    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", dtype=str)
    assert written.loc[0, "Nomenclature_Tier"] == "A"


def test_production_identity_policy_is_resolved_from_checked_in_config_at_the_stage_boundary(tmp_path) -> None:
    dup_c = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
    kestrel, advntr = _write_identity_aware_outputs(tmp_path, dup_c)
    config = nomenclature.load_nomenclature_config()
    config["thresholds"] = {**config["thresholds"], "min_support_for_high_confidence": 41}

    with mock.patch.object(nomenclature_annotate, "load_nomenclature_config", return_value=config):
        assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", dtype=str)
    assert written.loc[0, "Nomenclature_Tier"] == "B"
    assert "low-kmer-path-support" in written.loc[0, "Nomenclature_Flags"]
    assert "low-read-support" in written.loc[0, "Nomenclature_Flags"]


def test_production_reconciliation_fails_closed_on_malformed_internal_identity_metadata(tmp_path) -> None:
    dup_c = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
    kestrel, advntr = _write_identity_aware_outputs(tmp_path, dup_c)
    frame = pd.read_csv(kestrel, sep="\t", dtype=str)
    frame.loc[0, "__Identity_Selected_Raw_Representation_Key"] = '{"source": "kestrel", "values": []}'
    frame.to_csv(kestrel, sep="\t", index=False)

    with pytest.raises(ValueError):
        reconcile_caller_outputs(kestrel, advntr)


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


def _write_pair_bam(
    path: Path,
    inserted_at: int,
    base: str,
    minimum_kmer_depths: tuple[int | None, ...],
) -> None:
    """Write an indexed ``output.bam`` over the real ``X-X`` pair reference.

    Args:
        path: Destination ``output.bam``.
        inserted_at: 0-based reference position the insertion sits before.
        base: The inserted plus-strand base.
        minimum_kmer_depths: One optional XD value per resolved haplotype record.
    """
    pair = nomenclature.pair_sequence("X-X")
    assert pair is not None
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "X-X", "LN": len(pair)}]}
    trailing = len(pair) - inserted_at
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for index, minimum_kmer_depth in enumerate(minimum_kmer_depths):
            record = pysam.AlignedSegment(handle.header)
            record.query_name = f"resolved_haplotype_{index}"
            record.reference_id = 0
            record.reference_start = 0
            record.mapping_quality = 255
            record.cigarstring = f"{inserted_at}=1I{trailing}="
            record.query_sequence = pair[:inserted_at] + base + pair[inserted_at:]
            if minimum_kmer_depth is not None:
                record.set_tag("XD", minimum_kmer_depth)
            handle.write(record)
    pysam.index(str(path))  # type: ignore[attr-defined]


def test_annotate_kestrel_frame_uses_record_support_not_xd_for_a_bam_candidate(tmp_path: Path) -> None:
    """The direct stage path must consult a candidate and keep XD observational."""
    _write_pair_bam(
        tmp_path / "output.bam",
        inserted_at=62,
        base="C",
        minimum_kmer_depths=(2_147_483_647, 2_147_483_647),
    )
    frame = pd.DataFrame(
        [
            {
                "Motifs": "X-X",
                "POS": 999,
                "REF": "G",
                "ALT": "GG",
                "Motif_fasta": "X-X",
                "POS_fasta": 61,
            }
        ]
    )
    observed: list[tuple[Nomenclature | None, int | None]] = []
    row_haplotype_call = nomenclature_annotate._row_haplotype_call

    def tracked_row_haplotype_call(row: pd.Series, rescuer: BamRescuer) -> tuple[Nomenclature | None, int | None]:
        result = row_haplotype_call(row, rescuer)
        observed.append(result)
        return result

    with mock.patch.object(
        nomenclature_annotate,
        "_row_haplotype_call",
        side_effect=tracked_row_haplotype_call,
    ):
        named = annotate_kestrel_frame(frame, tmp_path)

    assert len(observed) == 1, "the tier-C VCF call must remain BAM-eligible"
    assert observed[0][0] is not None
    assert observed[0][1] == 2, "the tuple must carry haplotype-record support, not either XD value"
    assert named.loc[0, "Nomenclature"] == "58_59insG"
    assert FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT in named.loc[0, "Nomenclature_Flags"]


def test_row_verdicts_keep_each_haplotype_call_paired_with_its_vcf_row() -> None:
    """Compacting a missing first BAM call must not shift row two onto row one."""
    first_vcf = from_kestrel("X-X", 67, "G", "GG")
    second_vcf = Nomenclature("55_56insA", "insertion", "X", "B", (), None, None, 1, "kestrel_vcf")
    second_bam = Nomenclature("55delinsAT", "delins", "X", "B", (), None, None, 1, "kestrel_bam")

    verdicts = nomenclature_annotate._row_verdicts([first_vcf, second_vcf], [None, second_bam])

    assert [call.name if call is not None else None for call in verdicts] == ["59dupC", "55delinsAT"]


def test_the_cross_caller_stage_consults_resolved_haplotype_records(tmp_path) -> None:
    """The haplotype records must reach the only step that sees both callers.

    This is the production path for the whole vote. The Kestrel VCF places the
    insertion at ``59_60insG``; adVNTR and the BAM records say ``58_59insG``,
    and two sources from two callers outvote one. Without the BAM the cross-caller
    step sees only the two callers disagreeing and the wrong VCF placement stands.
    """
    kestrel, advntr = _write_disagreeing_outputs(tmp_path)
    _write_pair_bam(tmp_path / "output.bam", inserted_at=62, base="C", minimum_kmer_depths=(5, 181, 7_416, 8_704))

    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature"] == "58_59insG", "the haplotype records corroborated adVNTR"
    assert written.loc[0, "Nomenclature_Kestrel"] == "59_60insG", "Kestrel's own record is still reported"
    assert written.loc[0, "Nomenclature_adVNTR"] == "58_59insG"
    assert FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT in written.loc[0, "Nomenclature_Flags"]
    assert written.loc[0, "Nomenclature_Tier"] == "B", "four records, not XD, set the support boundary"


def test_the_haplotype_records_alone_do_not_outvote_the_vcf(tmp_path) -> None:
    """Same haplotype records, but adVNTR agreeing with Kestrel instead.

    The records are Kestrel's own output, so on their own they are not a second
    caller and must not overturn the VCF -- the guard against the earlier policy that
    cost `dupA` 6 correct calls out of 10.
    """
    kestrel, advntr = _write_disagreeing_outputs(tmp_path)
    # I22_2_C_LEN1 names 59_60insG, the same allele the Kestrel VCF names.
    advntr.write_text(advntr.read_text().replace("I23_2_C_LEN1", "I22_2_C_LEN1"))
    _write_pair_bam(tmp_path / "output.bam", inserted_at=62, base="C", minimum_kmer_depths=(5, 181, 7_416, 8_704))

    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature"] == "59_60insG", "one caller's BAM records may not veto its own VCF"


@pytest.mark.parametrize(
    "minimum_kmer_depths",
    [
        (None, None, None, None),
        (0, 0, 0, 0),
        (5, 181, 7_416, 8_704),
        (2_147_483_647, 2_147_483_647, 2_147_483_647, 2_147_483_647),
    ],
)
def test_xd_cannot_change_production_nomenclature_decisions(
    tmp_path: Path,
    minimum_kmer_depths: tuple[int | None, ...],
) -> None:
    """Changing only XD cannot change any displayed decision or candidate fetch."""
    case_dir = tmp_path / str(minimum_kmer_depths[0])
    case_dir.mkdir()
    kestrel, advntr = _write_disagreeing_outputs(case_dir)
    _write_pair_bam(case_dir / "output.bam", 62, "C", minimum_kmer_depths)

    rescuers: list[BamRescuer] = []
    open_rescuer = nomenclature_annotate._open_rescuer

    def tracked_open(output_dir: str | Path | None) -> BamRescuer | None:
        rescuer = open_rescuer(output_dir)
        if rescuer is not None:
            rescuers.append(rescuer)
        return rescuer

    with mock.patch.object(nomenclature_annotate, "_open_rescuer", side_effect=tracked_open):
        assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    decision_columns = [
        "Nomenclature",
        "Nomenclature_Tier",
        "Nomenclature_Flags",
        "Nomenclature_Kestrel",
        "Nomenclature_adVNTR",
    ]
    assert written.loc[0, decision_columns].to_dict() == {
        "Nomenclature": "58_59insG",
        "Nomenclature_Tier": "B",
        "Nomenclature_Flags": "caller-disagreement;known-variant;low-haplotype-record-support",
        "Nomenclature_Kestrel": "59_60insG",
        "Nomenclature_adVNTR": "58_59insG",
    }
    assert len(rescuers) == 1
    assert rescuers[0].opens == 1
    assert rescuers[0].fetches == 1


def test_haplotype_support_is_the_record_count_not_xd(tmp_path: Path) -> None:
    """High XD must not replace four resolved records in the support mapping."""
    _write_pair_bam(
        tmp_path / "output.bam",
        62,
        "C",
        (2_147_483_647, 2_147_483_647, 2_147_483_647, 2_147_483_647),
    )
    row = pd.Series({"Motif_fasta": "X-X", "POS_fasta": 61})
    rescuer = BamRescuer(tmp_path / "output.bam")
    try:
        call, support = nomenclature_annotate._row_haplotype_call(row, rescuer)
    finally:
        rescuer.close()

    assert call is not None
    assert support == 4
    assert rescuer.opens == 1
    assert rescuer.fetches == 1


def test_non_candidate_never_opens_or_fetches_the_bam(tmp_path: Path) -> None:
    """A settled caller agreement must return before constructing a rescuer."""
    _write_pair_bam(tmp_path / "output.bam", 62, "C", (5, 181, 7_416, 8_704))
    supports: dict[str, int | None] = {"kestrel_vcf": 40, "advntr": 40}
    with mock.patch.object(nomenclature_annotate, "_open_rescuer", wraps=nomenclature_annotate._open_rescuer) as opened:
        calls = nomenclature_annotate._haplotype_calls(
            [pd.Series({"Motif_fasta": "X-X", "POS_fasta": 61})],
            tmp_path,
            [from_kestrel("X-X", 67, "G", "GG")],
            [from_advntr("I22_2_G_LEN1")[0]],
            supports,
        )

    assert calls == []
    opened.assert_not_called()


def test_candidate_rows_share_one_rescuer_and_keep_fetch_order(tmp_path: Path) -> None:
    """Recreating or reordering the rescuer would detach evidence from its row."""
    _write_pair_bam(tmp_path / "output.bam", 62, "C", (5, 181, 7_416, 8_704))
    rows = [
        pd.Series({"Motif_fasta": "X-X", "POS_fasta": 61}),
        pd.Series({"Motif_fasta": "X-X", "POS_fasta": 67}),
    ]
    supports: dict[str, int | None] = {"kestrel_vcf": 40, "advntr": 40}
    fetched_loci: list[tuple[str, int]] = []
    rescuers: list[BamRescuer] = []
    open_rescuer = nomenclature_annotate._open_rescuer

    def tracked_open(output_dir: str | Path | None) -> BamRescuer | None:
        rescuer = open_rescuer(output_dir)
        assert rescuer is not None
        original_rescue = rescuer.rescue

        def tracked_rescue(contig: str, position: int):
            fetched_loci.append((contig, position))
            return original_rescue(contig, position)

        rescuer.rescue = tracked_rescue  # type: ignore[method-assign]
        rescuers.append(rescuer)
        return rescuer

    with mock.patch.object(nomenclature_annotate, "_open_rescuer", side_effect=tracked_open):
        calls = nomenclature_annotate._haplotype_calls(
            rows,
            tmp_path,
            [from_kestrel("X-X", 61, "T", "TC")],
            [from_advntr("I23_2_C_LEN1")[0]],
            supports,
        )

    assert len(rescuers) == 1
    assert rescuers[0].opens == 1
    assert rescuers[0].fetches == 2
    assert fetched_loci == [("X-X", 61), ("X-X", 67)]
    assert len(calls) == 2
    assert supports["kestrel_bam"] == 4, "the production mapping must carry the record count, not XD"


def test_a_missing_allele_cell_does_not_become_a_name(tmp_path) -> None:
    """A blank ``REF`` must make the row untranslatable, not become DNA.

    Read with ``dtype=str`` a missing cell is NaN, and ``str(NaN)`` is the text
    ``nan`` -- which every emptiness check accepts as an allele. An absent ``REF``
    produced a confident ``52_53delGC`` out of nothing.
    """
    frame = kestrel_stage_frame("final")
    frame.loc[0, "REF"] = None

    named = annotate_kestrel_frame(frame)
    assert named.loc[0, "Nomenclature"] == ""
    assert named.loc[0, "Nomenclature_Kestrel"] == ""


def test_each_row_keeps_its_own_name_when_a_file_reports_two_variants(tmp_path) -> None:
    """Row identity survives reconciliation.

    Broadcasting one joined summary onto every row made a two-variant file say both
    names on both rows, so neither row described itself any more.
    """
    kestrel = tmp_path / "kestrel_result.tsv"
    kestrel.write_text(
        "## VNtyper Kestrel result\n"
        "Motifs\tPOS\tREF\tALT\tConfidence\tEstimated_Depth_AlternateVariant\tNomenclature\t"
        "Nomenclature_Tier\tNomenclature_Flags\tAmbiguity_Interval\tRepeat_Form\tNomenclature_Note\t"
        "Nomenclature_Kestrel\tNomenclature_adVNTR\n"
        "X-X\t67\tG\tGG\tHigh_Precision\t24\t59dupC\tB\t\t\t\t\t59dupC\t\n"
        "X-X\t62\tG\tGC\tHigh_Precision\t24\t58_59insG\tB\t\t\t\t\t58_59insG\t\n"
    )
    advntr = tmp_path / "output_adVNTR_result.tsv"
    advntr.write_text(
        "VID\tVariant\tNumberOfSupportingReads\tNomenclature\tNomenclature_Tier\t"
        "Nomenclature_Flags\tAmbiguity_Interval\tRepeat_Form\tNomenclature_Note\t"
        "Nomenclature_Kestrel\tNomenclature_adVNTR\n"
        "1\tI22_2_G_LEN1\t24\t59dupC\tB\t\t\t\t\t\t59dupC\n"
        "2\tI23_2_C_LEN1\t24\t58_59insG\tB\t\t\t\t\t\t58_59insG\n"
    )

    assert reconcile_caller_outputs(kestrel, advntr) is True

    written = pd.read_csv(kestrel, sep="\t", comment="#", dtype=str)
    assert written.loc[0, "Nomenclature_Kestrel"] == "59dupC"
    assert written.loc[1, "Nomenclature_Kestrel"] == "58_59insG", "row 2 must not inherit row 1's name"

    mirrored = pd.read_csv(advntr, sep="\t", dtype=str)
    assert mirrored.loc[0, "Nomenclature_adVNTR"] == "59dupC"
    assert mirrored.loc[1, "Nomenclature_adVNTR"] == "58_59insG"


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
