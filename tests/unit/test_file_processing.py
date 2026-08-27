"""Unit tests for file_processing.py -- the indel split that feeds Kestrel.

`filter_vcf` decides what counts as an indel and `filter_indel_vcf` decides
which side of the insertion/deletion split it lands on. Both were untested.
A row dropped here never reaches motif processing, scoring or the report.
"""

import io
import json
import logging

import pandas as pd
import pytest

from tests.builders import kestrel_config
from vntyper.scripts import file_processing
from vntyper.scripts.confidence_assignment import calculate_depth_score_and_assign_confidence
from vntyper.scripts.file_processing import filter_indel_vcf, filter_vcf
from vntyper.scripts.scoring import split_depth_and_calculate_frame_score

pytestmark = pytest.mark.unit

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"

# The classification specification, enumerated exhaustively rather than spot-checked.
#
# Each entry is (ref_len, alt_len, filter_vcf outcome, filter_indel_vcf destination):
#   "kept"     -- filter_vcf writes the row to the indel file (REF and ALT differ in length)
#   "dropped"  -- filter_vcf silently omits it (equal lengths: a substitution, not an indel)
#   "rejected" -- filter_vcf raises ValueError, because two multi-base alleles break the
#                 output contract of the pinned Kestrel 1.0.1 (see the guard in filter_vcf)
#
# The destination column records what filter_indel_vcf does when handed the row DIRECTLY,
# bypassing filter_vcf. That bypass is the only way to observe routing for the rows
# filter_vcf rejects, and filter_indel_vcf deliberately carries no contract guard of its
# own -- its input is a derived file, not Kestrel's output.
CLASSIFICATION = [
    (1, 1, "dropped", "deletion"),
    (1, 2, "kept", "insertion"),
    (1, 3, "kept", "insertion"),
    (2, 1, "kept", "deletion"),
    (2, 2, "rejected", "deletion"),
    (2, 3, "rejected", "insertion"),
    (3, 1, "kept", "deletion"),
    (3, 2, "rejected", "deletion"),
    (3, 3, "rejected", "deletion"),
    (1, 5, "kept", "insertion"),
    (5, 1, "kept", "deletion"),
    (2, 5, "rejected", "insertion"),
    (5, 2, "rejected", "deletion"),
    (4, 7, "rejected", "insertion"),
    (7, 4, "rejected", "deletion"),
    (4, 4, "rejected", "deletion"),
]


def _vcf(tmp_path, name, *rows):
    path = tmp_path / name
    path.write_text(
        HEADER
        + "".join(f"X-5\t{i}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:GDP:DP\t1:50:5000\n" for i, (ref, alt) in enumerate(rows))
    )
    return path


def _rows(path):
    return [line for line in path.read_text().splitlines() if not line.startswith(("##", "#CHROM"))]


def _alleles(ref_len, alt_len):
    """Distinct alleles of the requested lengths; the bases themselves never matter here."""
    return "A" * ref_len, "C" * alt_len


def test_the_header_is_carried_through_verbatim(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CC"))), str(out))
    assert out.read_text().startswith("##fileformat=VCFv4.2\n#CHROM")


def test_a_one_base_to_many_insertion_is_kept(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGGCA"))), str(out))
    assert len(_rows(out)) == 1


def test_a_many_to_one_base_deletion_is_kept(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("CGGCA", "C"))), str(out))
    assert len(_rows(out)) == 1


def test_a_single_base_substitution_is_dropped(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "G"))), str(out))
    assert _rows(out) == []


@pytest.mark.parametrize(("ref_len", "alt_len", "outcome", "destination"), CLASSIFICATION)
def test_filter_vcf_classifies_every_ref_alt_length_pair(tmp_path, ref_len, alt_len, outcome, destination):
    """SPECIFICATION: filter_vcf keeps a row exactly when REF and ALT differ in length.

    The test is the length *difference*, so an indel is recognised whichever allele
    carries the extra bases -- not only when one of them happens to be a single base.
    Equal-length rows are substitutions and stay dropped.

    Rows whose REF and ALT are both multi-base are rejected outright rather than
    classified: they break the output contract of the Kestrel build VNtyper pins. See
    the companion guard test for that reasoning. `destination` is unused here and is
    asserted by test_filter_indel_vcf_routes_every_ref_alt_length_pair, which walks the
    same table -- the two together are the specification.
    """
    ref, alt = _alleles(ref_len, alt_len)
    src = _vcf(tmp_path, "in.vcf", (ref, alt))
    out = tmp_path / "out.vcf"

    if outcome == "rejected":
        with pytest.raises(ValueError, match="Kestrel 1.0.1"):
            filter_vcf(str(src), str(out))
        return

    filter_vcf(str(src), str(out))
    assert len(_rows(out)) == (1 if outcome == "kept" else 0)


@pytest.mark.parametrize(("ref_len", "alt_len", "outcome", "destination"), CLASSIFICATION)
def test_filter_indel_vcf_routes_every_ref_alt_length_pair(tmp_path, ref_len, alt_len, outcome, destination):
    """SPECIFICATION: filter_indel_vcf routes by comparing the REF and ALT lengths.

    ALT longer than REF is an insertion; REF longer than or equal to ALT is a deletion,
    which preserves the `>` boundary the module docstring has always documented. There is
    no unconditional `else` catch-all any more, so a multi-base-REF insertion no longer
    falls into the deletion file.

    Rows are fed in directly, bypassing filter_vcf, so routing is observable even for the
    pairs filter_vcf rejects. `outcome` is unused here; it is asserted by the companion
    test_filter_vcf_classifies_every_ref_alt_length_pair over the same table.
    """
    ref, alt = _alleles(ref_len, alt_len)
    src = _vcf(tmp_path, "indel.vcf", (ref, alt))
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"

    filter_indel_vcf(str(src), str(ins), str(dele))

    landed = "insertion" if _rows(ins) else "deletion"
    assert landed == destination
    assert len(_rows(ins)) + len(_rows(dele)) == 1


def test_a_record_with_two_multi_base_alleles_is_rejected_as_off_contract(tmp_path, caplog):
    """SPECIFICATION: filter_vcf refuses a record that breaks Kestrel's output contract.

    filter_vcf is the only consumer of the raw VCF Kestrel writes, so it is where the
    upstream contract is asserted. The pinned Kestrel 1.0.1 in
    vntyper/dependencies/kestrel/kestrel.jar anchors every indel on a single reference
    base -- VariantInsertion.getVcfRef() and VariantDeletion.getVcfAlt() return
    Character.toString(char) on every one of their return paths -- so it emits only
    1-vs-1, 1-vs-N and N-vs-1 records. A record with two multi-base alleles therefore
    means the Kestrel build has changed.

    That guarantee is a property of this vendored jar and is enforced nowhere else in
    VNtyper, so the check fails loudly and names the version it is defending, rather than
    letting a jar swap change reported labels in silence.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.file_processing")
    src = _vcf(tmp_path, "in.vcf", ("AC", "ACGGG"))

    with pytest.raises(ValueError) as excinfo:
        filter_vcf(str(src), str(tmp_path / "out.vcf"))

    assert "Kestrel 1.0.1" in str(excinfo.value)
    assert "kestrel.jar" in str(excinfo.value)
    assert "Kestrel 1.0.1" in caplog.text


def test_a_line_with_too_few_columns_raises(tmp_path):
    """filter_vcf requires the pinned build's full 10-column record; a short line is not tolerated."""
    path = tmp_path / "short.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tC\n")
    with pytest.raises(ValueError):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_insertions_and_deletions_are_written_to_separate_files(tmp_path):
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    filter_indel_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGG"), ("CGG", "C"))), str(ins), str(dele))
    assert len(_rows(ins)) == 1
    assert len(_rows(dele)) == 1


def test_both_output_files_get_the_header(tmp_path):
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    filter_indel_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGG"))), str(ins), str(dele))
    assert dele.read_text().startswith("##fileformat")
    assert _rows(dele) == []


def test_equal_length_alleles_are_classified_as_a_deletion(tmp_path):
    """filter_indel_vcf's insertion test is strictly `alt > ref`.

    The docstring is explicit that a deletion is REF longer than *or equal
    to* ALT, so a same-length pair (which filter_vcf would normally have
    dropped as a substitution before this function ever sees it) falls to
    the deletion side. This pins the documented `>` boundary directly, since
    none of the other single-row tests in this file exercise ref_len == alt_len here.
    """
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    path = tmp_path / "indel.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tC\tG\t.\tPASS\t.\n")
    filter_indel_vcf(str(path), str(ins), str(dele))
    assert _rows(ins) == []
    assert len(_rows(dele)) == 1


def test_a_multi_base_ref_insertion_is_routed_to_the_insertion_file(tmp_path):
    """SPECIFICATION: a multi-base-REF insertion goes to the insertion file.

    REF="AC", ALT="ACG" adds one base, so it is an insertion by any reading -- nothing was
    deleted. It used to fail the old `len(ref) == 1 and len(alt) > 1` test and drop through
    an unconditional `else` into the deletion file.

    This is a latent defect, not a live one: the pinned Kestrel 1.0.1 cannot emit a record
    with two multi-base alleles, and filter_vcf now rejects one outright, so no such row
    reaches this function in the pipeline. Had one ever arrived, the damage would have been
    a wrong 'Variant' label -- "Deletion" for an insertion -- in kestrel_result.tsv, the
    cohort tables and the HTML report, and a corrupted 'Allele_Change' in the adVNTR
    cross-match, where compute_allele_change() returns the REF ("AC") instead of the
    inserted bases ("G"). The frameshift verdict is NOT affected either way:
    scoring.extract_frameshifts derives 'direction' and 'frameshift_amount' from the REF and
    ALT lengths themselves and never reads the 'Variant' label or the originating file.
    """
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    path = tmp_path / "indel.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tAC\tACG\t.\tPASS\t.\n")
    filter_indel_vcf(str(path), str(ins), str(dele))
    assert len(_rows(ins)) == 1
    assert _rows(dele) == []


def test_a_record_with_empty_alleles_is_refused_not_silently_dropped(tmp_path):
    """#223: the one malformed shape ``filter_vcf`` would otherwise discard without a trace.

    The indel test is ``len(ref) != len(alt)``, and two empty alleles compare equal -- so a
    truncated record is classified as a substitution and dropped. A file of them yields two
    valid-header-but-empty derived VCFs, passes every header check, and is reported as a
    confident negative.
    """
    src = tmp_path / "output.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "chr1\t1\t.\t\t\t.\t.\tDP=1\tGT:GDP:DP\t1:1:1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="must carry both alleles"):
        filter_vcf(str(src), str(tmp_path / "out_indel.vcf"))


# --------------------------------------------------------------------------------------
# C4 (2026-08-26 screen): the FORMAT / sample-field contract
# --------------------------------------------------------------------------------------
#
# scoring.split_depth_and_calculate_frame_score assigns the colon-split Sample column
# positionally and never reads FORMAT, which is dropped before scoring. The boundary
# therefore proves the exact configured shape and the numeric depth values before any
# row -- including a substitution -- can be discarded. All 2,418,255 measured raw VCF
# rows carry GT:GDP:DP with one sample column and finite scalar depths, so the expected
# output delta is zero.


def _vcf_line(fmt, sample, ref="C", alt="CC"):
    return HEADER + f"X-5\t60\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt}\t{sample}\n"


def _replace_file_processing_config(monkeypatch, config):
    monkeypatch.setattr(
        file_processing.pkg_resources,
        "open_text",
        lambda *_args, **_kwargs: io.StringIO(json.dumps(config)),
    )


def test_a_format_other_than_the_configured_one_is_rejected(tmp_path, caplog):
    """A reordered FORMAT must stop the run, not swap the two depths in silence."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.file_processing")
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:DP:GDP", "1:5000:50"))

    with pytest.raises(ValueError, match="GT:GDP:DP") as excinfo:
        filter_vcf(str(path), str(tmp_path / "out.vcf"))

    assert str(excinfo.value) in caplog.text


def test_the_format_expectation_is_required_config_not_a_hardcoded_fallback(tmp_path, monkeypatch):
    """Deleting the required contract must fail rather than silently restore a default."""
    _replace_file_processing_config(monkeypatch, {"file_processing": {"snv_length": 1}})
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", "1:50:5000"))

    with pytest.raises(KeyError, match="kestrel_format"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_a_changed_configured_format_controls_the_contract(tmp_path, monkeypatch):
    """The configured expectation is read, while the two depth positions remain GDP then DP."""
    config = {"file_processing": {"snv_length": 1, "kestrel_format": "DEL:GDP:DP"}}
    _replace_file_processing_config(monkeypatch, config)
    old_path = tmp_path / "old.vcf"
    old_path.write_text(_vcf_line("GT:GDP:DP", "1:50:5000"))

    with pytest.raises(ValueError, match="DEL:GDP:DP"):
        filter_vcf(str(old_path), str(tmp_path / "old-out.vcf"))

    changed_path = tmp_path / "changed.vcf"
    changed_path.write_text(_vcf_line("DEL:GDP:DP", "1:50:5000"))
    changed_out = tmp_path / "changed-out.vcf"
    filter_vcf(str(changed_path), str(changed_out))
    assert len(_rows(changed_out)) == 1


@pytest.mark.parametrize("configured_format", ["GT:DP:GDP", "GT:GDP:DP:EXTRA"])
def test_the_configured_format_must_preserve_the_positional_depth_layout(tmp_path, monkeypatch, configured_format):
    """Configuration cannot bless a layout the unchanged positional scorer would misread."""
    config = {"file_processing": {"snv_length": 1, "kestrel_format": configured_format}}
    _replace_file_processing_config(monkeypatch, config)
    path = tmp_path / "in.vcf"
    sample = "1:5000:50" if configured_format == "GT:DP:GDP" else "1:50:5000:unused"
    path.write_text(_vcf_line(configured_format, sample))

    with pytest.raises(ValueError, match="GDP and DP in positions 2 and 3"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


@pytest.mark.parametrize(
    "sample",
    [
        "1:.:5000",
        "1:text:5000",
        "1::5000",
        "1:50:.",
        "1:50:text",
        "1:50:",
        "1:1_0:5000",
        "1:50:5_000",
        "1:١٠:5000",
    ],
)
def test_a_nonnumeric_depth_is_rejected(tmp_path, sample):
    """GDP and DP must both be numeric before scoring can consume them."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", sample))

    with pytest.raises(ValueError, match="finite scalar numeric"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


@pytest.mark.parametrize(
    "sample",
    [
        "1:nan:5000",
        "1:inf:5000",
        "1:-inf:5000",
        "1:50:nan",
        "1:50:inf",
        "1:50:-inf",
    ],
)
def test_a_nonfinite_depth_is_rejected(tmp_path, sample):
    """NaN and infinities are numeric syntax but cannot represent read depths."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", sample))

    with pytest.raises(ValueError, match="finite scalar numeric"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


@pytest.mark.parametrize(
    ("sample", "invalid_depth"),
    [
        ("1:-1:5000", "GDP='-1'"),
        ("1:50:-5000", "DP='-5000'"),
        ("1:-50:-5000", "GDP='-50'"),
    ],
)
def test_a_negative_depth_is_rejected(tmp_path, caplog, sample, invalid_depth):
    """Negative GDP/DP values are impossible read counts and must not reach confidence scoring."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.file_processing")
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", sample))

    with pytest.raises(ValueError, match="cannot be negative") as excinfo:
        filter_vcf(str(path), str(tmp_path / "out.vcf"))

    assert invalid_depth in str(excinfo.value)
    assert str(excinfo.value) in caplog.text


@pytest.mark.parametrize("sample", ["1:0:5000", "1:50:0", "1:0:0"])
def test_a_zero_depth_remains_valid_at_the_nonnegative_boundary(tmp_path, sample):
    """Zero is a possible read count and must not be rejected with negative depths."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", sample))
    out = tmp_path / "out.vcf"

    filter_vcf(str(path), str(out))

    assert _rows(out)[0].endswith(f"\t{sample}")


@pytest.mark.parametrize("sample", ["1:3420,120:5000", "1:3420:5000,120"])
def test_a_multi_value_depth_is_rejected_not_scored_as_zero(tmp_path, sample):
    """A comma-joined GDP or DP is not one scalar depth and must fail closed."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", sample))

    with pytest.raises(ValueError, match="finite scalar numeric"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_a_sample_with_the_wrong_subfield_count_is_rejected(tmp_path):
    """A sample wider than FORMAT must not reach pandas as an opaque length mismatch."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", "1:50:5000:7"))

    with pytest.raises(ValueError, match="colon-separated"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_an_extra_sample_column_is_rejected(tmp_path):
    """The pinned producer has one sample; a second sample would otherwise be ignored."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", "1:50:5000").rstrip("\n") + "\t1:25:5000\n")

    with pytest.raises(ValueError, match="exactly 10"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_a_data_line_without_format_and_sample_columns_is_rejected(tmp_path):
    """Without FORMAT and the sample field, depths cannot be checked."""
    path = tmp_path / "in.vcf"
    path.write_text(HEADER + "X-5\t60\t.\tC\tCC\t.\tPASS\t.\n")

    with pytest.raises(ValueError, match="exactly 10"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_a_substitution_row_with_bad_depths_is_still_rejected(tmp_path):
    """The contract applies before equal-length alleles are discarded."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", "1:.:5000", ref="C", alt="G"))

    with pytest.raises(ValueError, match="finite scalar numeric"):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_a_conforming_raw_row_reaches_scoring_with_the_pinned_depth_values(tmp_path):
    """The accepted GDP/DP positions become the expected depth score and confidence."""
    path = tmp_path / "in.vcf"
    path.write_text(_vcf_line("GT:GDP:DP", "1:3420:5000"))
    out = tmp_path / "out.vcf"

    filter_vcf(str(path), str(out))
    fields = _rows(out)[0].split("\t")
    frame = pd.DataFrame({"REF": [fields[3]], "ALT": [fields[4]], "Sample": [fields[9]]})
    scored = split_depth_and_calculate_frame_score(frame)
    result = calculate_depth_score_and_assign_confidence(scored, kestrel_config()).iloc[0]

    assert result["Estimated_Depth_AlternateVariant"] == 3420
    assert result["Estimated_Depth_Variant_ActiveRegion"] == 5000
    assert result["Depth_Score"] == pytest.approx(0.684)
    assert result["Confidence"] == "High_Precision*"
