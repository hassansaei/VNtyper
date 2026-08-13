"""Cross-caller reconciliation and the confidence tier ladder.

Tier A is the only tier that emits a bare number, so it is the only tier that can
be confidently wrong. It therefore requires evidence no single caller can supply:
two independent sources agreeing after normalisation, a motif context that matches
canonical X, and read support at or above a threshold.

Research use only.
"""

import pytest

from vntyper.scripts.nomenclature import (
    Nomenclature,
    from_advntr,
    from_kestrel,
    reconcile,
    render,
)

pytestmark = pytest.mark.unit


def _kestrel_dupc() -> Nomenclature:
    return from_kestrel("X-X", 67, "G", "GG")


def _advntr_dupc() -> Nomenclature:
    (call,) = from_advntr("I22_2_G_LEN1")
    return call


# ---------------------------------------------------------------------------
# Promotion to tier A
# ---------------------------------------------------------------------------


def test_two_agreeing_callers_with_support_reach_tier_a() -> None:
    merged = reconcile(_kestrel_dupc(), _advntr_dupc(), support=24)
    assert merged.tier == "A"
    assert merged.name == "59dupC"


def test_one_caller_alone_never_reaches_tier_a() -> None:
    """However clean it looks. Kestrel's insG family looks clean and is wrong."""
    assert reconcile(_kestrel_dupc(), support=500).tier != "A"
    assert reconcile(_advntr_dupc(), support=500).tier != "A"


def test_thin_support_blocks_tier_a() -> None:
    """A 1-read consensus is reported honestly, never silently accepted."""
    merged = reconcile(_kestrel_dupc(), _advntr_dupc(), support=1)
    assert merged.tier != "A"
    assert "low-read-support" in merged.flags


def test_absent_support_blocks_tier_a() -> None:
    """Unknown support is not the same as sufficient support."""
    assert reconcile(_kestrel_dupc(), _advntr_dupc()).tier != "A"


def test_a_divergent_motif_context_blocks_tier_a() -> None:
    """S-C projects onto X; the projection is not itself evidence."""
    merged = reconcile(from_kestrel("S-C", 67, "G", "GG"), _advntr_dupc(), support=40)
    assert merged.tier != "A"


# ---------------------------------------------------------------------------
# Disagreement
# ---------------------------------------------------------------------------


def test_callers_naming_different_alleles_never_reach_tier_a() -> None:
    kestrel = _kestrel_dupc()
    (advntr,) = from_advntr("I23_2_C_LEN1")
    merged = reconcile(kestrel, advntr, support=40)
    assert merged.tier != "A"
    assert "caller-disagreement" in merged.flags


def test_disagreement_on_event_class_falls_to_tier_c() -> None:
    """If the two callers do not even agree what kind of event it is, no number."""
    kestrel = _kestrel_dupc()
    (advntr,) = from_advntr("D17_2&D18_2&D19_2&D20_2&D21_2")
    merged = reconcile(kestrel, advntr, support=40)
    assert merged.tier == "C"
    assert merged.name is None


def test_reconciling_nothing_yields_an_undetermined_call() -> None:
    merged = reconcile()
    assert merged.tier == "C"
    assert merged.name is None


def test_duplicate_rows_from_one_caller_are_not_two_sources() -> None:
    """Independence is counted over the calls that carry the name, not all calls.

    An unnamed call -- an adVNTR state in an unmappable repeat unit, say -- must not
    donate a second `source` to two duplicate rows from one caller. Counting every
    supplied call let a single caller's placement be promoted as though two had
    agreed on it, and Kestrel's insG placements are exactly the ones that are wrong.
    """
    kestrel = from_kestrel("X-X", 61, "T", "TC")
    (unnamed,) = from_advntr("I69_3_G_LEN1")
    assert unnamed.name is None

    merged = reconcile(kestrel, kestrel, unnamed, support=10)
    assert merged.tier != "A"
    assert render(merged) != kestrel.name


def test_support_must_belong_to_the_agreeing_evidence() -> None:
    """A well-covered but unrelated observation cannot lend depth to a thin agreement."""
    kestrel = from_kestrel("X-X", 62, "G", "GC")
    (advntr,) = from_advntr("I23_2_C_LEN1")
    (unrelated,) = from_advntr("I69_3_G_LEN1")

    thin = reconcile(kestrel, advntr, unrelated, supports={"kestrel_vcf": 40, "advntr": 1})
    assert thin.tier != "A", "a 1-read adVNTR agreement must not borrow depth"
    assert "low-read-support" in thin.flags

    real = reconcile(kestrel, advntr, supports={"kestrel_vcf": 40, "advntr": 24})
    assert real.tier == "A"
    assert real.name == "58_59insG"


def test_an_agreement_is_only_as_strong_as_its_weakest_source() -> None:
    merged = reconcile(
        from_kestrel("X-X", 62, "G", "GC"),
        from_advntr("I23_2_C_LEN1")[0],
        supports={"kestrel_vcf": 500, "advntr": 2},
    )
    assert merged.tier != "A"


def test_several_simultaneous_events_never_reach_tier_a() -> None:
    """A locus reporting three events at once is not one simple allele.

    This is the real ins25bp sample pair_4094: Kestrel truncates the 25 bp
    insertion to 1 bp and adVNTR's first state agrees with that truncation, so
    passing only the first state promoted a wrong `30dupC` to tier A. Passing
    every event the caller reported is what makes the locus's complexity visible.
    """
    kestrel = from_kestrel("5C-M", 31, "G", "GG")
    events = [call for state in ("I51_2_G_LEN1", "I69_3_G_LEN1", "D49_3&I49_3_C_LEN5") for call in from_advntr(state)]
    merged = reconcile(kestrel, *events, support=42)
    assert merged.tier != "A"
    assert merged.name is None


# ---------------------------------------------------------------------------
# Rendering: tier decides what may be shown
# ---------------------------------------------------------------------------


def test_tier_a_renders_the_bare_name() -> None:
    merged = reconcile(_kestrel_dupc(), _advntr_dupc(), support=24)
    assert render(merged) == "59dupC"


def test_tier_b_renders_the_event_but_never_a_bare_number() -> None:
    """Spec §3.3: tier B emits the event and interval, no bare number."""
    call = reconcile(_kestrel_dupc(), support=40)
    text = render(call)
    assert "59dupC" not in text
    assert "duplication" in text


def test_a_locus_with_nothing_to_say_does_not_claim_a_frameshift() -> None:
    """`frameshift +0` is not a frameshift.

    It is what a locus with no usable call reduces to, and printing it states
    something untrue about a locus we know nothing about. On the benchmark 11 of
    the 200 samples rendered this.
    """
    assert render(reconcile()) == "allele undetermined"
    assert "frameshift" not in render(reconcile())


def test_tier_c_renders_a_frameshift_statement_with_no_number() -> None:
    undetermined = Nomenclature(
        name=None,
        event="insertion",
        unit=None,
        tier="C",
        flags=("allele-unrepresentable-in-vcf",),
        ambiguity=None,
        repeat_form=None,
        net_length=1,
        source="kestrel_vcf",
    )
    text = render(undetermined)
    assert text == "frameshift +1, allele undetermined"
    assert not any(character.isdigit() and character != "1" for character in text.replace("+1", ""))
