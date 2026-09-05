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
    bam_contigs,
    bam_contigs_mixed_conventions,
    bam_contigs_unknown_chr1_length,
    bam_contigs_without_chr1,
    bam_header,
    bam_header_malformed_chr1_length,
    bam_header_missing_chr1_length,
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
    assert conf["confidence_assignment"]["reporting_floor"] == 0.00469
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


# ---------------------------------------------------------------------------
# bam_contigs: the parsed shape production hands to assembly reconciliation
# ---------------------------------------------------------------------------


def test_bam_contigs_returns_the_shape_parse_contigs_from_header_returns() -> None:
    """`reconcile_assembly` consumes this exact shape: name (str), length (int)."""
    contigs = bam_contigs()
    assert isinstance(contigs, list)
    assert contigs, "builder produced no contigs; every assertion below would be vacuous"
    for contig in contigs:
        assert set(contig) == {"name", "length"}, f"unexpected keys: {sorted(contig)}"
        assert isinstance(contig["name"], str)
        assert isinstance(contig["length"], int)


@pytest.mark.parametrize("convention", ["ucsc", "ensembl", "ncbi"])
@pytest.mark.parametrize("assembly", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("n_contigs", [1, 3, 25])
def test_bam_header_and_bam_contigs_agree_through_the_production_parser(
    convention: str, assembly: str, n_contigs: int
) -> None:
    """The two builders must never drift; production's own parser is the judge."""
    from vntyper.scripts.fastq_bam_processing import parse_contigs_from_header

    header = bam_header(convention=convention, assembly=assembly, n_contigs=n_contigs)
    parsed = parse_contigs_from_header(header)

    assert len(parsed) == n_contigs
    assert parsed == bam_contigs(convention=convention, assembly=assembly, n_contigs=n_contigs)


def test_bam_header_renders_a_caller_supplied_contig_list() -> None:
    header = bam_header(contigs=[{"name": "chrZ", "length": 42}])
    assert "@SQ\tSN:chrZ\tLN:42" in header
    assert "chr1" not in header


# --- the defect: NCBI accessions are chromosome-specific, not one version -----


@pytest.mark.parametrize(
    "assembly,expected",
    [
        # Real GRCh37/GRCh38 RefSeq accessions. The old builder emitted ".10" for
        # every GRCh37 chromosome and ".11" for every GRCh38 one - accessions that
        # do not exist, so every consumer went green for the wrong reason.
        ("GRCh37", ["NC_000001.10", "NC_000002.11", "NC_000003.11", "NC_000004.11", "NC_000005.9"]),
        ("GRCh38", ["NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10"]),
    ],
)
def test_ncbi_accessions_are_chromosome_specific(assembly: str, expected: list[str]) -> None:
    names = [contig["name"] for contig in bam_contigs(convention="ncbi", assembly=assembly, n_contigs=5)]
    assert names == expected


def test_ncbi_accessions_do_not_drift_from_the_production_version_tables() -> None:
    """**Circular by construction, and kept anyway** -- as a drift detector, not an oracle.

    The builder calls ``_construct_ncbi_accession`` itself, so this can only fail if the
    builder stops delegating; it can never catch an error the two share. The five
    independently written accessions in
    :func:`test_ncbi_accessions_are_chromosome_specific` above are the real assertion, and
    they are what would catch a wrong version table. What this adds is coverage of
    chromosomes 6-25, where no hand-written expectation exists: if someone reimplements the
    builder's naming instead of delegating, the two diverge here.
    """
    from vntyper.scripts.chromosome_utils import _construct_ncbi_accession

    contigs = bam_contigs(convention="ncbi", assembly="GRCh37", n_contigs=25)
    assert len(contigs) == 25, "vacuity guard: the loop below would assert nothing"

    for number, contig in enumerate(contigs, start=1):
        assert contig["name"] == _construct_ncbi_accession(number, "hg19")


def test_ncbi_chr1_matches_the_accession_shipped_in_config() -> None:
    """config.json's `known_chromosome_naming` short-circuits chr1 in production."""
    import json
    from pathlib import Path

    config = json.loads((Path(__file__).resolve().parents[2] / "vntyper" / "config.json").read_text(encoding="utf-8"))
    naming = config["bam_processing"]["known_chromosome_naming"]
    for assembly in ("GRCh37", "GRCh38"):
        first = bam_contigs(convention="ncbi", assembly=assembly, n_contigs=1)[0]
        assert first["name"] == naming[assembly]["ncbi"]


@pytest.mark.parametrize(
    "convention,expected_tail",
    [("ucsc", ["chr22", "chrX", "chrY", "chrM"]), ("ensembl", ["22", "X", "Y", "MT"])],
)
def test_sex_and_mitochondrial_contigs_use_production_names(convention: str, expected_tail: list[str]) -> None:
    """No real BAM contains `chr23`; production never builds that name either."""
    names = [contig["name"] for contig in bam_contigs(convention=convention, n_contigs=25)]
    assert names[-4:] == expected_tail


def test_bam_contigs_rejects_an_unknown_assembly_and_convention() -> None:
    with pytest.raises(ValueError, match="Unknown assembly"):
        bam_contigs(assembly="hg19")
    with pytest.raises(ValueError, match="Unknown convention"):
        bam_contigs(convention="nonsense")


# --- fixtures the assembly-reconciliation work needs -------------------------


@pytest.mark.parametrize("convention", ["ucsc", "ensembl", "ncbi"])
def test_the_canonical_contigs_are_detected_as_their_assembly(convention: str) -> None:
    """The positive control the degenerate fixtures below are measured against."""
    from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

    for assembly in ("GRCh37", "GRCh38"):
        detected = detect_assembly_from_chr1_length(bam_contigs(convention=convention, assembly=assembly))
        # NCBI used to be pinned to None here: production searched only
        # "chr1"/"1", so an NCBI-named chr1 was invisible to it. That gap is
        # closed, so all three conventions now answer with the real build.
        assert detected == assembly


def test_unknown_chr1_length_fixture_matches_no_assembly() -> None:
    from vntyper.scripts.chromosome_utils import CHR1_LENGTHS, detect_assembly_from_chr1_length

    contigs = bam_contigs_unknown_chr1_length()
    chr1 = next(c for c in contigs if c["name"] == "chr1")
    assert chr1["length"] not in CHR1_LENGTHS.values()
    assert detect_assembly_from_chr1_length(contigs) is None


def test_without_chr1_fixture_really_has_no_chr1() -> None:
    from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

    for convention in ("ucsc", "ensembl", "ncbi"):
        contigs = bam_contigs_without_chr1(convention=convention)
        assert contigs, "fixture produced nothing to reason about"
        assert not any(c["name"].lower() in ("chr1", "1") for c in contigs)
        assert detect_assembly_from_chr1_length(contigs) is None


def test_missing_and_malformed_chr1_lengths_are_dropped_by_the_production_parser() -> None:
    """A header the parser silently discards is indistinguishable from no chr1."""
    from vntyper.scripts.fastq_bam_processing import parse_contigs_from_header

    for header in (bam_header_missing_chr1_length(), bam_header_malformed_chr1_length()):
        parsed = parse_contigs_from_header(header)
        assert parsed, "the whole header was discarded; the assertion below would be vacuous"
        assert not any(c["name"] == "chr1" for c in parsed), (
            "parse_contigs_from_header keeps only contigs with an integer length"
        )


def test_mixed_convention_fixture_defeats_convention_detection() -> None:
    from vntyper.scripts.chromosome_utils import detect_naming_convention

    contigs = bam_contigs_mixed_conventions()
    names = [c["name"] for c in contigs]
    assert len(names) >= 3, "too few contigs for a majority vote to mean anything"
    assert detect_naming_convention(names) == "unknown"
    # ...but chr1 is still UCSC-named and still carries a decisive length, so the
    # build is recoverable even when the convention is not.
    from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

    assert detect_assembly_from_chr1_length(contigs) == "GRCh37"


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
