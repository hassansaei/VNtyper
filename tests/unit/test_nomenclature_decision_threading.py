"""Resolved nomenclature and cross-match policy reaches every decision family."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts import nomenclature
from vntyper.scripts.cross_match import cross_match_variants
from vntyper.scripts.identity_reconciliation import IdentityReconciliationPolicy
from vntyper.scripts.nomenclature_annotate import annotate_kestrel_frame
from vntyper.scripts.nomenclature_bam import BamConsensus, BamRescuer
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def test_profile_projects_one_typed_nomenclature_configuration() -> None:
    run = resolve_run_configuration()
    builder = getattr(nomenclature, "decision_config_from_component", None)

    assert callable(builder)
    decision = builder(run.nomenclature)
    assert decision.canonical_unit == run.nomenclature["canonical_unit"]
    assert decision.bam_flank == 8
    assert decision.bam_thin_haplotype_record_support == 3


@pytest.mark.parametrize(
    "change, expected",
    [
        (lambda component: component.pop("canonical_unit"), "nomenclature component fields differ"),
        (lambda component: component.__setitem__("motifs", "not-a-mapping"), "nomenclature motifs must be a mapping"),
    ],
)
def test_nomenclature_configuration_rejects_malformed_components(change, expected: str) -> None:
    component = dict(resolve_run_configuration().nomenclature)
    change(component)
    builder = getattr(nomenclature, "decision_config_from_component", None)

    assert callable(builder)
    with pytest.raises(ValueError, match=expected):
        builder(component)


def test_identity_policy_retains_the_two_value_five_names_and_units() -> None:
    component = resolve_run_configuration().nomenclature["identity_reconciliation"]

    policy = IdentityReconciliationPolicy.from_component(component)

    assert policy.kestrel_min_alternate_kmer_path_depth == 5
    assert policy.advntr_min_sequencing_read_support == 5
    assert policy.source_evidence_units["kestrel_vcf"] == "alternate-kmer-path-depth"
    assert policy.source_evidence_units["advntr"] == "sequencing-reads"
    assert policy.source_evidence_units["kestrel_vcf"] != policy.source_evidence_units["advntr"]


def test_identity_policy_rejects_an_incomplete_profile_component() -> None:
    component = dict(resolve_run_configuration().nomenclature["identity_reconciliation"])
    component.pop("source_evidence_units")

    with pytest.raises(ValueError, match="identity reconciliation component fields differ"):
        IdentityReconciliationPolicy.from_component(component)


def test_custom_nomenclature_context_cannot_fall_back_before_an_early_return() -> None:
    with pytest.raises(ValueError, match="custom Nomenclature run context requires an explicit resolved component"):
        annotate_kestrel_frame(pd.DataFrame(), custom_context_active=True)


def test_bam_consensus_uses_its_resolved_record_threshold() -> None:
    thin = BamConsensus(
        kind="insertion",
        start=10,
        ref_span=0,
        inserted=1,
        supporting_haplotype_records=3,
        fetched_haplotype_records=3,
        distinct_edit_count=1,
        supporting_record_minimum_kmer_depths=(1, 1, 1),
        thin_haplotype_record_support=4,
    )
    sufficient = BamConsensus(
        kind="insertion",
        start=10,
        ref_span=0,
        inserted=1,
        supporting_haplotype_records=4,
        fetched_haplotype_records=4,
        distinct_edit_count=1,
        supporting_record_minimum_kmer_depths=(1, 1, 1, 1),
        thin_haplotype_record_support=4,
    )

    assert thin.is_thin is True
    assert sufficient.is_thin is False


def test_bam_rescuer_retains_the_resolved_flank_and_thin_threshold(tmp_path: Path) -> None:
    rescuer = BamRescuer(tmp_path / "output.bam", flank=17, thin_haplotype_record_support=4)

    assert rescuer.flank == 17
    assert rescuer.thin_haplotype_record_support == 4


def test_cross_match_uses_the_resolved_required_advntr_disposition() -> None:
    component = dict(resolve_run_configuration().cross_match)
    component["required_advntr_evidence_disposition"] = "identity-insufficient"
    kestrel = [{"REF": "C", "ALT": "CC", "POS": 67}]
    advntr = [
        {
            "REF": "C",
            "ALT": "CC",
            "POS": 1,
            "Evidence_Disposition": "identity-insufficient",
        }
    ]

    result = cross_match_variants(
        kestrel,
        advntr,
        resolved_component=component,
        custom_context_active=True,
    )

    assert result["overall_match"] == "Yes"


def test_custom_cross_match_context_cannot_fall_back_to_package_policy() -> None:
    with pytest.raises(ValueError, match="custom cross-match run context requires an explicit resolved component"):
        cross_match_variants([], [], custom_context_active=True)


def test_cross_match_rejects_an_unsupported_resolved_disposition() -> None:
    component = dict(resolve_run_configuration().cross_match)
    component["required_advntr_evidence_disposition"] = "unknown"

    with pytest.raises(ValueError, match="required_advntr_evidence_disposition is unsupported"):
        cross_match_variants([], [], resolved_component=component)
