"""Cross-caller reconciliation and the confidence tier ladder.

Tier A is the only tier that emits a bare number, so it is the only tier that can
be confidently wrong. It therefore requires evidence no single caller can supply:
two independent sources agreeing after normalisation, a motif context that matches
canonical X, and source-specific support at or above a threshold.

Research use only.
"""

from dataclasses import replace

import pytest

from vntyper.scripts.identity_reconciliation import (
    IdentityReconciliationObservation,
    IdentityReconciliationPolicy,
)
from vntyper.scripts.molecular_identity import (
    EvidenceDisposition,
    EvidenceDispositionValue,
    IdentityTranslation,
    MolecularIdentity,
    make_coding_edit,
    make_molecular_identity,
)
from vntyper.scripts.nomenclature import (
    CALLER_OF,
    FLAG_CALLER_DISAGREEMENT,
    FLAG_LOW_EVIDENCE_SUPPORT,
    FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
    FLAG_LOW_KMER_PATH_SUPPORT,
    FLAG_LOW_READ_SUPPORT,
    FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,
    KNOWN_VARIANTS,
    MIN_SUPPORT_FOR_TIER_A,
    Nomenclature,
    from_advntr,
    from_kestrel,
    reconcile,
    render,
)
from vntyper.scripts.nomenclature_bam import refine

pytestmark = pytest.mark.unit

IDENTITY_POLICY = IdentityReconciliationPolicy(5, 5)
_DUPC_IDENTITY = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
_DUPA_IDENTITY = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))


def _identity_observation(
    identity: MolecularIdentity | None,
    call: Nomenclature,
    *,
    kmer_depth: int | None = None,
    advntr_reads: int | None = None,
    disposition: EvidenceDispositionValue = "admissible",
    presentation_call_index: int | None = None,
) -> IdentityReconciliationObservation:
    translation = (
        IdentityTranslation(identity, "resolved", None, False)
        if identity is not None
        else IdentityTranslation(None, "unresolved", "missing-motif-context", False)
    )
    return IdentityReconciliationObservation(
        translation=translation,
        display_name=call.name,
        source=call.source,
        event=call.event,
        net_length=call.net_length,
        flags=frozenset(call.flags),
        disposition=EvidenceDisposition(disposition),
        known_variant_match=call.name in KNOWN_VARIANTS,
        kestrel_alternate_kmer_path_depth=kmer_depth,
        advntr_sequencing_read_support=advntr_reads,
        presentation_call_index=presentation_call_index,
    )


def test_typed_adapter_rejects_equal_display_names_with_different_identities() -> None:
    kestrel = _kestrel_dupc()
    advntr = _advntr_dupc()

    merged = reconcile(
        kestrel,
        advntr,
        identity_observations=(
            _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40),
            _identity_observation(_DUPA_IDENTITY, advntr, advntr_reads=40),
        ),
        identity_policy=IDENTITY_POLICY,
    )

    assert merged.name == "59dupC", "the packaged legacy display row remains primary"
    assert merged.tier == "B"
    assert FLAG_CALLER_DISAGREEMENT in merged.flags


def test_typed_adapter_agrees_by_identity_without_replacing_the_legacy_display_row() -> None:
    kestrel = Nomenclature("legacy-display", "duplication", "X", "B", (), None, None, 1, "kestrel_vcf")
    advntr = Nomenclature("other-display", "duplication", "X", "B", (), None, None, 1, "advntr")

    merged = reconcile(
        kestrel,
        advntr,
        identity_observations=(
            _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40),
            _identity_observation(_DUPC_IDENTITY, advntr, advntr_reads=40),
        ),
        identity_policy=IDENTITY_POLICY,
    )

    assert merged.name == "legacy-display"
    assert merged.source == "reconciled"
    assert FLAG_CALLER_DISAGREEMENT not in merged.flags


def test_typed_adapter_rejects_call_observation_cardinality_drift() -> None:
    kestrel = _kestrel_dupc()
    with pytest.raises(ValueError, match="one identity observation per call"):
        reconcile(kestrel, identity_observations=(), identity_policy=IDENTITY_POLICY)


def test_explicit_binding_rejects_a_duplicate_even_when_every_call_index_is_present() -> None:
    kestrel = _kestrel_dupc()
    advntr = _advntr_dupc()
    with pytest.raises(ValueError, match="exactly one identity observation per call"):
        reconcile(
            kestrel,
            advntr,
            identity_observations=(
                _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40, presentation_call_index=0),
                _identity_observation(_DUPC_IDENTITY, advntr, advntr_reads=40, presentation_call_index=1),
                _identity_observation(_DUPC_IDENTITY, advntr, advntr_reads=40, presentation_call_index=1),
            ),
            identity_policy=IDENTITY_POLICY,
        )


def test_explicit_binding_rejects_a_missing_call_index() -> None:
    kestrel = _kestrel_dupc()
    advntr = _advntr_dupc()
    with pytest.raises(ValueError, match="exactly one identity observation per call"):
        reconcile(
            kestrel,
            advntr,
            identity_observations=(
                _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40, presentation_call_index=0),
            ),
            identity_policy=IDENTITY_POLICY,
        )


def test_explicit_binding_rejects_an_out_of_range_extra_call_index() -> None:
    kestrel = _kestrel_dupc()
    with pytest.raises(ValueError, match="exactly one identity observation per call"):
        reconcile(
            kestrel,
            identity_observations=(
                _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40, presentation_call_index=0),
                _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40, presentation_call_index=1),
            ),
            identity_policy=IDENTITY_POLICY,
        )


def test_explicit_binding_accepts_reordered_observations_and_preserves_call_order() -> None:
    kestrel = _kestrel_dupc()
    advntr = Nomenclature("same-identity-other-display", "duplication", "X", "B", (), None, None, 1, "advntr")

    merged = reconcile(
        kestrel,
        advntr,
        identity_observations=(
            _identity_observation(_DUPC_IDENTITY, advntr, advntr_reads=40, presentation_call_index=1),
            _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40, presentation_call_index=0),
        ),
        identity_policy=IDENTITY_POLICY,
    )

    assert merged.name == "59dupC"


def test_explicit_binding_rejects_a_display_bearing_unbound_observation() -> None:
    kestrel = _kestrel_dupc()
    malformed_unbound = replace(
        _identity_observation(_DUPC_IDENTITY, _advntr_dupc(), advntr_reads=40),
        advntr_mean_coverage=100,
    )

    with pytest.raises(ValueError, match="unbound observation must be a complete adVNTR row"):
        reconcile(
            kestrel,
            identity_observations=(
                _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40, presentation_call_index=0),
                malformed_unbound,
            ),
            identity_policy=IDENTITY_POLICY,
        )


def test_identity_insufficient_different_event_cannot_withdraw_the_admissible_kestrel_presentation() -> None:
    kestrel = _kestrel_dupc()
    advntr = Nomenclature("58_59insG", "insertion", "X", "B", (), None, None, 1, "advntr")

    merged = reconcile(
        kestrel,
        advntr,
        identity_observations=(
            _identity_observation(_DUPC_IDENTITY, kestrel, kmer_depth=40),
            _identity_observation(
                _DUPA_IDENTITY,
                advntr,
                advntr_reads=400,
                disposition="identity-insufficient",
            ),
        ),
        identity_policy=IDENTITY_POLICY,
    )

    assert merged.name == "59dupC"
    assert merged.event == "duplication"
    assert merged.tier == "B"
    assert FLAG_CALLER_DISAGREEMENT not in merged.flags


def test_unresolved_vcf_abstention_preserves_exact_legacy_low_support_flags() -> None:
    kestrel = _kestrel_dupc()
    advntr = _advntr_dupc()
    supports = {"kestrel_vcf": 2, "advntr": 4}

    legacy = reconcile(kestrel, advntr, supports=supports)
    typed = reconcile(
        kestrel,
        advntr,
        supports=supports,
        identity_observations=(
            _identity_observation(None, kestrel, kmer_depth=2),
            _identity_observation(_DUPC_IDENTITY, advntr, advntr_reads=4),
        ),
        identity_policy=IDENTITY_POLICY,
    )

    assert typed.flags == legacy.flags
    assert {FLAG_LOW_KMER_PATH_SUPPORT, FLAG_LOW_READ_SUPPORT}.issubset(typed.flags)


def test_unbound_bam_advntr_abstention_preserves_exact_legacy_low_support_flags() -> None:
    kestrel = from_kestrel("X-X", 61, "T", "TC")
    bam = replace(
        _from_haplotypes("58_59insG", "insertion"),
        flags=(FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,),
    )
    advntr = from_advntr("I23_2_C_LEN1")[0]
    supports = {"kestrel_vcf": 40, "kestrel_bam": 2, "advntr": 4}

    legacy = reconcile(kestrel, bam, advntr, supports=supports)
    typed = reconcile(
        kestrel,
        bam,
        advntr,
        supports=supports,
        identity_observations=(
            _identity_observation(_DUPA_IDENTITY, kestrel, kmer_depth=40),
            _identity_observation(None, bam),
            _identity_observation(_DUPC_IDENTITY, advntr, advntr_reads=4),
        ),
        identity_policy=IDENTITY_POLICY,
    )

    assert typed.flags == legacy.flags
    assert {
        FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
        FLAG_LOW_READ_SUPPORT,
        FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,
    }.issubset(typed.flags)


def test_future_source_abstention_preserves_generic_legacy_low_support_flag() -> None:
    future = Nomenclature("59dupC", "duplication", "X", "B", (), None, None, 1, "future_caller")
    advntr = _advntr_dupc()
    supports = {"future_caller": 2, "advntr": 4}

    legacy = reconcile(future, advntr, supports=supports)
    typed = reconcile(
        future,
        advntr,
        supports=supports,
        identity_observations=(
            _identity_observation(None, future),
            _identity_observation(_DUPC_IDENTITY, advntr, advntr_reads=4),
        ),
        identity_policy=IDENTITY_POLICY,
    )

    assert typed.flags == legacy.flags
    assert FLAG_LOW_EVIDENCE_SUPPORT in typed.flags


def _kestrel_dupc() -> Nomenclature:
    return from_kestrel("X-X", 67, "G", "GG")


def _advntr_dupc() -> Nomenclature:
    (call,) = from_advntr("I22_2_G_LEN1")
    return call


def _from_haplotypes(name: str, event: str) -> Nomenclature:
    """A call recovered from Kestrel's resolved ``output.bam`` haplotypes."""
    return Nomenclature(
        name=name,
        event=event,
        unit="X",
        tier="B",
        flags=(),
        ambiguity=None,
        repeat_form=None,
        net_length=1,
        source="kestrel_bam",
    )


# ---------------------------------------------------------------------------
# Independent-caller majority
# ---------------------------------------------------------------------------


def test_two_independent_callers_outvote_a_single_dissenting_caller() -> None:
    """The real ``insG`` case, where the VCF alone is wrong.

    Kestrel places the whole ``insG`` family one base 3' of truth. adVNTR and the
    resolved haplotypes say ``58_59insG``; the Kestrel VCF alone says ``59_60insG``.
    Two sources from two different callers outvote one, so the name that two
    independent lines of evidence support is the one reported.
    """
    merged = reconcile(
        from_kestrel("X-X", 61, "T", "TC"),
        from_advntr("I23_2_C_LEN1")[0],
        _from_haplotypes("58_59insG", "insertion"),
    )
    assert merged.name == "58_59insG"


def test_kestrel_agreeing_with_its_own_alignment_is_not_a_majority() -> None:
    """``kestrel_bam`` is Kestrel's alignment, not a second opinion on it.

    Counting the two as independent would let one caller outvote the other simply by
    agreeing with itself, which is not evidence about the allele.
    """
    merged = reconcile(
        from_kestrel("X-X", 61, "T", "TC"),
        _from_haplotypes("59_60insG", "insertion"),
        from_advntr("I23_2_C_LEN1")[0],
    )
    assert merged.name == "59_60insG", "the VCF stands; nothing outvoted it"
    assert "caller-disagreement" in merged.flags


def test_kestrel_agreeing_with_its_own_alignment_never_reaches_tier_a() -> None:
    """The false-confidence path that adding BAM records as a third source opens.

    ``59dupC`` is a described variant with ample support and no divergent context,
    so every other tier-A condition is met: only the independence of the two sources
    stands between this and a confidently stated name resting on one caller.
    """
    merged = reconcile(
        _kestrel_dupc(),
        _from_haplotypes("59dupC", "duplication"),
        supports={"kestrel_vcf": 40, "kestrel_bam": 40},
    )
    assert merged.tier != "A"


def test_the_haplotype_records_corroborating_a_second_caller_reach_tier_a() -> None:
    """Independence is the test, not which source the name came from."""
    merged = reconcile(
        _advntr_dupc(),
        _from_haplotypes("59dupC", "duplication"),
        supports={"advntr": 40, "kestrel_bam": 40},
    )
    assert merged.tier == "A"


def test_unknown_depth_from_one_source_blocks_tier_a() -> None:
    """Unknown is not sufficient, and one caller's depth is not the other's evidence.

    Dropping the unknown and taking the minimum of what remained let a well-covered
    Kestrel record supply the depth for an adVNTR agreement whose own depth column was
    blank -- a tier-A name resting on one source's reads while claiming two.
    """
    merged = reconcile(
        _kestrel_dupc(),
        _advntr_dupc(),
        supports={"kestrel_vcf": 40, "advntr": None},
    )
    assert merged.tier != "A"


def test_two_corroborated_alleles_are_a_conflict_rather_than_a_vote() -> None:
    """Only an unambiguous majority decides.

    adVNTR may report several events at one locus, so it can corroborate two
    different alleles at once -- here it backs the VCF and BAM alleles. Both
    then clear the independence bar, and picking whichever happened to be seen first
    would settle a real conflict by input order.
    """
    merged = reconcile(
        from_kestrel("X-X", 61, "T", "TC"),  # 59_60insG, backed by the VCF alone
        _from_haplotypes("58_59insG", "insertion"),
        _from_haplotypes("57_58insG", "insertion"),
        from_advntr("I23_2_C_LEN1")[0],  # 58_59insG -- corroborates the first BAM call
        from_advntr("I24_2_C_LEN1")[0],  # 57_58insG -- corroborates the second
    )
    assert merged.tier != "A"
    assert "caller-disagreement" in merged.flags
    assert merged.name == "59_60insG", "two corroborated alleles must not be settled by input order"


def test_support_is_taken_from_the_sources_backing_the_chosen_name() -> None:
    """A dissenting source's depth is not the agreement's depth.

    The two callers backing ``58_59insG`` are both deep; the Kestrel VCF dissents on
    one unit of k-mer-path support. Letting the dissent set support would withhold tier A from an
    agreement that is in fact well covered.
    """
    merged = reconcile(
        from_kestrel("X-X", 61, "T", "TC"),
        from_advntr("I23_2_C_LEN1")[0],
        _from_haplotypes("58_59insG", "insertion"),
        supports={"kestrel_vcf": 1, "advntr": 40, "kestrel_bam": 40},
    )
    assert merged.name == "58_59insG"
    assert not {
        FLAG_LOW_EVIDENCE_SUPPORT,
        FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
        FLAG_LOW_KMER_PATH_SUPPORT,
        FLAG_LOW_READ_SUPPORT,
    }.intersection(merged.flags)


@pytest.mark.parametrize(
    ("low_source", "other_source", "expected_flag"),
    [
        ("kestrel_bam", "advntr", FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT),
        ("kestrel_vcf", "advntr", FLAG_LOW_KMER_PATH_SUPPORT),
        ("advntr", "kestrel_vcf", FLAG_LOW_READ_SUPPORT),
        ("future_caller", "kestrel_vcf", FLAG_LOW_EVIDENCE_SUPPORT),
    ],
)
def test_low_support_flags_name_each_backing_sources_evidence_unit(
    low_source: str,
    other_source: str,
    expected_flag: str,
) -> None:
    calls = [
        Nomenclature("59dupC", "duplication", "X", "B", (), None, None, 1, low_source),
        Nomenclature("59dupC", "duplication", "X", "B", (), None, None, 1, other_source),
    ]

    merged = reconcile(*calls, supports={low_source: 4, other_source: 40})

    assert merged.name == "59dupC"
    assert render(merged) == "59dupC"
    assert merged.tier == "B"
    assert expected_flag in merged.flags
    assert {
        FLAG_LOW_EVIDENCE_SUPPORT,
        FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
        FLAG_LOW_KMER_PATH_SUPPORT,
        FLAG_LOW_READ_SUPPORT,
    }.intersection(merged.flags) == {expected_flag}


def test_all_low_known_backing_units_are_reported_truthfully() -> None:
    merged = reconcile(
        _kestrel_dupc(),
        _advntr_dupc(),
        supports={"kestrel_vcf": 2, "advntr": 4},
    )

    assert merged.name == "59dupC"
    assert {FLAG_LOW_KMER_PATH_SUPPORT, FLAG_LOW_READ_SUPPORT}.issubset(merged.flags)
    assert FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT not in merged.flags
    assert FLAG_LOW_EVIDENCE_SUPPORT not in merged.flags


def test_a_low_non_backing_source_cannot_contribute_a_flag() -> None:
    merged = reconcile(
        from_kestrel("X-X", 61, "T", "TC"),
        from_advntr("I23_2_C_LEN1")[0],
        _from_haplotypes("58_59insG", "insertion"),
        supports={"kestrel_vcf": 1, "advntr": 40, "kestrel_bam": 40},
    )

    assert merged.name == "58_59insG"
    assert not {
        FLAG_LOW_EVIDENCE_SUPPORT,
        FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
        FLAG_LOW_KMER_PATH_SUPPORT,
        FLAG_LOW_READ_SUPPORT,
    }.intersection(merged.flags)


@pytest.mark.parametrize(
    "calls",
    [
        (_kestrel_dupc(), _advntr_dupc()),
        (_advntr_dupc(), _kestrel_dupc()),
    ],
)
def test_one_unknown_backing_value_suppresses_every_low_support_token(calls: tuple[Nomenclature, ...]) -> None:
    merged = reconcile(*calls, supports={"kestrel_vcf": None, "advntr": 2})

    assert merged.tier == "B"
    assert not {
        FLAG_LOW_EVIDENCE_SUPPORT,
        FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
        FLAG_LOW_KMER_PATH_SUPPORT,
        FLAG_LOW_READ_SUPPORT,
    }.intersection(merged.flags)


def test_all_known_backing_values_enable_the_truthful_low_support_token() -> None:
    merged = reconcile(_kestrel_dupc(), _advntr_dupc(), supports={"kestrel_vcf": 40, "advntr": 2})

    assert merged.tier == "B"
    assert FLAG_LOW_READ_SUPPORT in merged.flags
    assert FLAG_LOW_KMER_PATH_SUPPORT not in merged.flags


def test_scalar_support_retains_its_legacy_read_support_contract() -> None:
    merged = reconcile(_kestrel_dupc(), _advntr_dupc(), support=4)

    assert merged.tier == "B"
    assert FLAG_LOW_READ_SUPPORT in merged.flags
    assert FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT not in merged.flags


def test_bam_support_of_four_is_not_thin_but_remains_below_the_tier_a_threshold() -> None:
    merged = reconcile(
        _from_haplotypes("59dupC", "duplication"),
        _advntr_dupc(),
        supports={"kestrel_bam": 4, "advntr": 40},
    )

    assert MIN_SUPPORT_FOR_TIER_A == 5
    assert FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT not in merged.flags
    assert FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT in merged.flags
    assert merged.tier == "B"


def test_thin_haplotype_records_cannot_settle_a_caller_disagreement() -> None:
    """Renaming the thin token must not bypass the established refine veto."""
    disagreement = reconcile(
        _kestrel_dupc(),
        from_advntr("D17_2&D18_2&D19_2&D20_2&D21_2")[0],
        supports={"kestrel_vcf": 40, "advntr": 40},
    )
    thin = Nomenclature(
        name="58_59insG",
        event="insertion",
        unit="X",
        tier="B",
        flags=(FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,),
        ambiguity=None,
        repeat_form=None,
        net_length=1,
        source="kestrel_bam",
    )

    assert disagreement.name is None
    assert refine(disagreement, thin) is disagreement


def test_support_at_five_is_not_low_and_preserves_tier_a() -> None:
    merged = reconcile(
        _from_haplotypes("59dupC", "duplication"),
        _advntr_dupc(),
        supports={"kestrel_bam": 5, "advntr": 40},
    )

    assert FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT not in merged.flags
    assert merged.tier == "A"


def test_kestrel_vcf_and_bam_remain_one_caller() -> None:
    assert CALLER_OF["kestrel_vcf"] == CALLER_OF["kestrel_bam"] == "kestrel"
    merged = reconcile(
        _kestrel_dupc(),
        _from_haplotypes("59dupC", "duplication"),
        supports={"kestrel_vcf": 40, "kestrel_bam": 40},
    )

    assert merged.tier == "B"


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
    assert merged.tier != "A", "one caller's placement must not be promoted as two"


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


def test_a_name_is_shown_whenever_one_could_be_computed() -> None:
    """The tier is a confidence indicator, not a gate on the name.

    Suppressing the number below tier A discarded most of the useful output: on
    the 200-sample benchmark 129 samples had the correct name computed and only 46
    displayed it. Withholding a name that was right 83 times, to avoid showing one
    that was wrong 0 times, is a bad trade -- and it denies the reader information
    they can weigh themselves against the tier and the flags.
    """
    call = reconcile(_kestrel_dupc(), support=40)
    assert call.tier != "A"
    assert render(call) == "59dupC"


def test_the_tier_still_reports_the_lower_confidence() -> None:
    """Showing the name must not quietly upgrade how much it is trusted."""
    call = reconcile(_kestrel_dupc(), support=40)
    assert call.tier == "B"


def test_a_locus_with_nothing_to_say_does_not_claim_a_frameshift() -> None:
    """`frameshift +0` is not a frameshift.

    It is what a locus with no usable call reduces to, and printing it states
    something untrue about a locus we know nothing about. On the benchmark 11 of
    the 200 samples rendered this.
    """
    assert render(reconcile()) == "allele undetermined"
    assert "frameshift" not in render(reconcile())


def test_tier_c_still_carries_no_position_because_it_has_no_name() -> None:
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
