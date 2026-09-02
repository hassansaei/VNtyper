"""Validated whole-locus BAM evidence universe."""

import pytest

from vntyper.scripts.molecular_identity import make_coding_edit, make_molecular_identity
from vntyper.scripts.nomenclature_bam_evidence import BamIdentityEvidence, BamLocusEvidence
from vntyper.scripts.nomenclature_bam_replay import BamReplayArtifact, BamReplayLocus
from vntyper.scripts.nomenclature_bam_universe import validated_whole_locus_bam_evidence

pytestmark = pytest.mark.unit

_A = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))
_C = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))


def _artifact(identity=_A) -> BamReplayArtifact:
    record = BamIdentityEvidence((identity,), ((0,),), 10)
    return BamReplayArtifact((BamReplayLocus((0,), "observed", BamLocusEvidence((record,), 1, {identity: 1})),))


def test_validated_universe_preserves_one_unweighted_record() -> None:
    evidence = validated_whole_locus_bam_evidence(_artifact(), {0: _A})

    assert evidence is not None
    assert evidence.eligible_record_count == 1
    assert evidence.counts == {_A: 1}


def test_validated_universe_rejects_an_identity_bound_to_the_wrong_ordinal() -> None:
    with pytest.raises(ValueError, match="binding.*authoritative"):
        validated_whole_locus_bam_evidence(_artifact(_C), {0: _A})


def test_validated_universe_rejects_incomplete_ordinal_custody() -> None:
    with pytest.raises(ValueError, match="ordinals.*authoritative"):
        validated_whole_locus_bam_evidence(_artifact(), {0: _A, 1: _C})


def test_validated_universe_keeps_unavailable_bam_missing() -> None:
    artifact = BamReplayArtifact((BamReplayLocus((0,), "unavailable", None),))

    assert validated_whole_locus_bam_evidence(artifact, {0: _A}) is None
