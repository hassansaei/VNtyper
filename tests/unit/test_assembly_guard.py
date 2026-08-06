"""Tests for `assembly_guard.reconcile_assembly`.

`--reference-assembly` is taken on faith today. Declaring hg19 on an hg38
alignment slices a region ~30 kb away from the MUC1 VNTR, finds nothing, and
reports a confident negative. `reconcile_assembly` is the pure verdict that
makes that contradiction visible; the caller decides what to do about it.

Every expected value below is written out independently -- none is derived from
the function under test, and none is derived from a builder-computed value.
"""

import pytest

pytestmark = pytest.mark.unit

from tests.builders import (  # noqa: E402
    bam_contigs,
    bam_contigs_mixed_conventions,
    bam_contigs_unknown_chr1_length,
    bam_contigs_without_chr1,
    bam_header_malformed_chr1_length,
    bam_header_missing_chr1_length,
)
from vntyper.scripts.assembly_guard import (  # noqa: E402
    STATUS_AGREE,
    STATUS_CONFLICT,
    STATUS_MISMATCH,
    STATUS_UNDETERMINED,
    UNKNOWN,
    AssemblyVerdict,
    reconcile_assembly,
)

# Real chr1 lengths, written out rather than imported from production.
GRCH37_CHR1 = 249250621
GRCH38_CHR1 = 248956422


def parsed(header: str) -> list[dict]:
    """Run a builder header through production's own parser."""
    from vntyper.scripts.fastq_bam_processing import parse_contigs_from_header

    return parse_contigs_from_header(header)


# ---------------------------------------------------------------------------
# The verdict table
# ---------------------------------------------------------------------------

# (case id, declared, contigs, status, coordinate_system, declared_coordinate_system, naming_convention)
CASES = [
    # --- the declared build matches the alignment ---------------------------
    ("ucsc hg19 on GRCh37", "hg19", bam_contigs("ucsc", "GRCh37"), STATUS_AGREE, "GRCh37", "GRCh37", "ucsc"),
    ("ucsc hg38 on GRCh38", "hg38", bam_contigs("ucsc", "GRCh38"), STATUS_AGREE, "GRCh38", "GRCh38", "ucsc"),
    ("build name GRCh37", "GRCh37", bam_contigs("ncbi", "GRCh37"), STATUS_AGREE, "GRCh37", "GRCh37", "ncbi"),
    ("build name GRCh38", "GRCh38", bam_contigs("ncbi", "GRCh38"), STATUS_AGREE, "GRCh38", "GRCh38", "ncbi"),
    ("ncbi alias", "hg19_ncbi", bam_contigs("ncbi", "GRCh37"), STATUS_AGREE, "GRCh37", "GRCh37", "ncbi"),
    ("ensembl alias", "hg38_ensembl", bam_contigs("ensembl", "GRCh38"), STATUS_AGREE, "GRCh38", "GRCh38", "ensembl"),
    # The naming convention and the build are separate questions: a header
    # whose convention cannot be decided still has a decidable build.
    (
        "mixed convention, right build",
        "hg19",
        bam_contigs_mixed_conventions("GRCh37"),
        STATUS_AGREE,
        "GRCh37",
        "GRCh37",
        UNKNOWN,
    ),
    # --- the declared build contradicts the alignment -----------------------
    ("hg19 declared on GRCh38", "hg19", bam_contigs("ucsc", "GRCh38"), STATUS_MISMATCH, "GRCh38", "GRCh37", "ucsc"),
    (
        "hg38 declared on GRCh37",
        "hg38",
        bam_contigs("ensembl", "GRCh37"),
        STATUS_MISMATCH,
        "GRCh37",
        "GRCh38",
        "ensembl",
    ),
    ("GRCh38 declared on GRCh37", "GRCh38", bam_contigs("ncbi", "GRCh37"), STATUS_MISMATCH, "GRCh37", "GRCh38", "ncbi"),
    (
        "mixed convention, wrong build",
        "hg38",
        bam_contigs_mixed_conventions("GRCh37"),
        STATUS_MISMATCH,
        "GRCh37",
        "GRCh38",
        UNKNOWN,
    ),
    # --- nothing can be concluded ------------------------------------------
    ("no contigs at all", "hg19", [], STATUS_UNDETERMINED, UNKNOWN, "GRCh37", UNKNOWN),
    (
        "chr1 length matches neither build",
        "hg19",
        bam_contigs_unknown_chr1_length("ucsc"),
        STATUS_UNDETERMINED,
        UNKNOWN,
        "GRCh37",
        "ucsc",
    ),
    (
        "no chr1 in the header",
        "hg38",
        bam_contigs_without_chr1("ucsc", "GRCh38"),
        STATUS_UNDETERMINED,
        UNKNOWN,
        "GRCh38",
        "ucsc",
    ),
    (
        "chr1 has no LN field",
        "hg19",
        parsed(bam_header_missing_chr1_length()),
        STATUS_UNDETERMINED,
        UNKNOWN,
        "GRCh37",
        "ucsc",
    ),
    (
        "chr1 LN is not a number",
        "hg19",
        parsed(bam_header_malformed_chr1_length()),
        STATUS_UNDETERMINED,
        UNKNOWN,
        "GRCh37",
        "ucsc",
    ),
    (
        "declared assembly is not recognised",
        "b37",
        bam_contigs("ucsc", "GRCh37"),
        STATUS_UNDETERMINED,
        "GRCh37",
        UNKNOWN,
        "ucsc",
    ),
]


@pytest.mark.parametrize(
    "declared,contigs,status,coordinate_system,declared_coordinate_system,naming_convention",
    [case[1:] for case in CASES],
    ids=[case[0] for case in CASES],
)
def test_the_verdict_table(declared, contigs, status, coordinate_system, declared_coordinate_system, naming_convention):
    verdict = reconcile_assembly(declared, contigs)

    assert isinstance(verdict, AssemblyVerdict)
    assert verdict.status == status
    assert verdict.declared == declared
    assert verdict.coordinate_system == coordinate_system
    assert verdict.declared_coordinate_system == declared_coordinate_system
    assert verdict.naming_convention == naming_convention
    assert verdict.message


def test_the_table_covers_all_three_statuses():
    """Guards the table above from silently degenerating into one status."""
    statuses = {case[3] for case in CASES}
    assert statuses == {STATUS_AGREE, STATUS_MISMATCH, STATUS_UNDETERMINED}


# ---------------------------------------------------------------------------
# The mismatch message is the whole product: it is what a human reads
# ---------------------------------------------------------------------------


def test_a_mismatch_message_names_both_builds_and_the_declared_input():
    verdict = reconcile_assembly("hg19", bam_contigs("ucsc", "GRCh38"))

    assert verdict.status == STATUS_MISMATCH
    assert "hg19" in verdict.message
    assert "GRCh37" in verdict.message, "the declared build must be named"
    assert "GRCh38" in verdict.message, "the detected build must be named"


def test_a_mismatch_message_names_both_builds_in_the_other_direction():
    verdict = reconcile_assembly("hg38", bam_contigs("ucsc", "GRCh37"))

    assert verdict.status == STATUS_MISMATCH
    assert "hg38" in verdict.message
    assert "GRCh37" in verdict.message
    assert "GRCh38" in verdict.message


def test_an_agreeing_message_does_not_read_as_a_problem():
    verdict = reconcile_assembly("hg19", bam_contigs("ucsc", "GRCh37"))

    assert verdict.status == STATUS_AGREE
    assert "GRCh38" not in verdict.message


# ---------------------------------------------------------------------------
# chr1_length is the evidence behind the verdict
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "assembly,expected_length",
    [("GRCh37", GRCH37_CHR1), ("GRCh38", GRCH38_CHR1)],
)
@pytest.mark.parametrize("convention", ["ucsc", "ensembl", "ncbi"])
def test_the_verdict_carries_the_chr1_length_it_reasoned_from(convention, assembly, expected_length):
    verdict = reconcile_assembly(assembly, bam_contigs(convention, assembly))
    assert verdict.chr1_length == expected_length


def test_chr1_length_is_none_when_chr1_is_absent():
    verdict = reconcile_assembly("hg19", bam_contigs_without_chr1("ucsc"))
    assert verdict.chr1_length is None


def test_chr1_length_is_reported_even_when_it_matches_no_build():
    """The unrecognised length is the single most useful thing to log."""
    verdict = reconcile_assembly("hg19", bam_contigs_unknown_chr1_length("ucsc"))
    assert verdict.status == STATUS_UNDETERMINED
    assert verdict.chr1_length is not None
    assert verdict.chr1_length not in (GRCH37_CHR1, GRCH38_CHR1)
    assert str(verdict.chr1_length) in verdict.message.replace(",", "")


# ---------------------------------------------------------------------------
# The guard must be total: it is called on whatever a header happened to hold
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "contigs",
    [
        [],
        [{}],
        [{"length": 249250621}],
        [{"name": None, "length": 249250621}],
        [{"name": "chr1"}],
        [{"name": "chr1", "length": None}],
        [{"name": "chr1", "length": "249250621"}],
        [{"name": "chr1", "length": 249250621.0}],
        [None, {"name": "chr1", "length": 249250621}],
    ],
)
def test_degenerate_contig_lists_produce_a_verdict_rather_than_an_exception(contigs):
    verdict = reconcile_assembly("hg19", contigs)
    assert verdict.status in (STATUS_AGREE, STATUS_MISMATCH, STATUS_UNDETERMINED)


def test_a_usable_contig_still_decides_despite_junk_beside_it():
    """The junk entry above must be skipped, not swallow the real chr1."""
    verdict = reconcile_assembly("hg19", [None, {"name": "chr1", "length": GRCH37_CHR1}])
    assert verdict.status == STATUS_AGREE
    assert verdict.chr1_length == GRCH37_CHR1


# ---------------------------------------------------------------------------
# A header naming chromosome 1 twice must not let contig order decide
# ---------------------------------------------------------------------------


def hybrid_header_contigs(first: tuple[str, int], second: tuple[str, int]) -> list[dict]:
    """Two aliases of chromosome 1 plus one ordinary contig, in the given order."""
    return [
        {"name": first[0], "length": first[1]},
        {"name": second[0], "length": second[1]},
        {"name": "chr2", "length": 242193529},
    ]


CONFLICTING_ORDERS = [
    ("ensembl first", hybrid_header_contigs(("1", GRCH38_CHR1), ("chr1", GRCH37_CHR1))),
    ("ucsc first", hybrid_header_contigs(("chr1", GRCH37_CHR1), ("1", GRCH38_CHR1))),
    ("ncbi first", hybrid_header_contigs(("NC_000001.11", GRCH38_CHR1), ("chr1", GRCH37_CHR1))),
]


@pytest.mark.parametrize("label,contigs", CONFLICTING_ORDERS, ids=[c[0] for c in CONFLICTING_ORDERS])
@pytest.mark.parametrize("declared", ["hg19", "hg38", "GRCh37", "GRCh38"])
def test_conflicting_recognised_chr1_entries_conflict_whatever_the_order(label, contigs, declared):
    """
    A hybrid header carrying both `1` and `chr1` at two *recognised* build lengths is a
    `conflict`: the header asserts GRCh37 and GRCh38 at once, and no ordering of it makes
    that true. Taking the first match made the answer depend on the order the aligner
    happened to write the header, so that is still refused -- but the verdict is now fatal
    rather than `undetermined`.

    This reverses an earlier decision made in review, which had this case `undetermined`
    and therefore proceeding. It was wrong: `undetermined` means "the guard could not
    answer the question", and the caller responds by trusting the declared assembly
    unchecked. Here the header has answered, contradictorily. Whichever alias naming
    resolution then picks, the declared build's MUC1 coordinates get applied to a contig
    of the other build's length -- the wrong slice, and a plausible false negative, which
    is the exact failure class this guard exists to make visible. Contradictory evidence
    must not fail open.
    """
    verdict = reconcile_assembly(declared, contigs)

    assert verdict.status == STATUS_CONFLICT, label
    assert verdict.chr1_length is None, label
    assert verdict.coordinate_system == UNKNOWN, label


@pytest.mark.parametrize("declared", ["hg19", "hg38"])
@pytest.mark.parametrize(
    "aliases",
    [
        pytest.param((("chr1", GRCH37_CHR1), ("1", 12345)), id="recognised_then_unrecognised"),
        pytest.param((("1", 12345), ("chr1", GRCH38_CHR1)), id="unrecognised_then_recognised"),
        pytest.param((("chr1", 12345), ("1", 999)), id="two_unrecognised"),
    ],
)
def test_an_unrecognised_chr1_length_stays_undetermined(declared, aliases):
    """Only *mutually contradictory recognised* evidence is fatal.

    An unrecognised chr1 length is not a contradiction, it is an unanswered question: a
    patched, decoy-carrying or non-human reference can name chromosome 1 at a length this
    guard has no entry for. That must keep warning and proceeding exactly as before, or
    the fatal conflict above becomes a guard that rejects ordinary inputs.
    """
    verdict = reconcile_assembly(declared, hybrid_header_contigs(*aliases))

    assert verdict.status == STATUS_UNDETERMINED
    assert verdict.chr1_length is None


def test_the_verdict_is_identical_for_both_orderings_of_the_same_conflict():
    """The whole point: ordering is not evidence, so it must not move the verdict."""
    forward = reconcile_assembly("hg19", hybrid_header_contigs(("1", GRCH38_CHR1), ("chr1", GRCH37_CHR1)))
    reverse = reconcile_assembly("hg19", hybrid_header_contigs(("chr1", GRCH37_CHR1), ("1", GRCH38_CHR1)))

    assert forward.status == reverse.status
    assert forward.chr1_length is None and reverse.chr1_length is None
    assert forward.coordinate_system == reverse.coordinate_system == UNKNOWN
    assert forward.message == reverse.message


def test_a_conflicting_header_names_both_contigs_and_both_lengths_in_its_message():
    """This message becomes an exception the user sees, so it must show what contradicts."""
    verdict = reconcile_assembly("hg19", hybrid_header_contigs(("1", GRCH38_CHR1), ("chr1", GRCH37_CHR1)))

    plain = verdict.message.replace(",", "")
    assert str(GRCH37_CHR1) in plain
    assert str(GRCH38_CHR1) in plain
    assert "chr1" in verdict.message
    assert "GRCh37" in verdict.message and "GRCh38" in verdict.message
    assert "conflicting" in verdict.message


@pytest.mark.parametrize(
    "declared,length,expected",
    [("hg19", GRCH37_CHR1, STATUS_AGREE), ("hg38", GRCH37_CHR1, STATUS_MISMATCH)],
)
def test_chromosome_one_repeated_at_the_same_length_still_decides(declared, length, expected):
    """Duplicated aliases are only a problem when they *disagree*."""
    verdict = reconcile_assembly(declared, hybrid_header_contigs(("1", length), ("chr1", length)))

    assert verdict.status == expected
    assert verdict.chr1_length == length


@pytest.mark.parametrize("bad_length", [None, "249250621", 249250621.0, True])
def test_a_malformed_chr1_alias_no_longer_masks_a_well_formed_one(bad_length):
    """
    Previously the first chr1-named contig ended the search whatever its length was, so a
    header whose leading alias carried a junk length read as "no chr1" even though a usable
    one followed. Skipping the junk is a decision, not an accident of order.
    """
    verdict = reconcile_assembly(
        "hg19",
        [{"name": "1", "length": bad_length}, {"name": "chr1", "length": GRCH37_CHR1}],
    )

    assert verdict.status == STATUS_AGREE
    assert verdict.chr1_length == GRCH37_CHR1


def test_the_verdict_is_immutable():
    """F stores this in the pipeline summary; it must not be edited in place."""
    import dataclasses

    verdict = reconcile_assembly("hg19", bam_contigs("ucsc", "GRCh37"))
    with pytest.raises(dataclasses.FrozenInstanceError):
        verdict.status = STATUS_MISMATCH  # type: ignore[misc]


def test_reconcile_assembly_does_not_mutate_its_input():
    contigs = bam_contigs("ucsc", "GRCh38")
    before = [dict(contig) for contig in contigs]
    reconcile_assembly("hg19", contigs)
    assert contigs == before
