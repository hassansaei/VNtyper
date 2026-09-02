"""Validated whole-locus projection of retained BAM identity evidence."""

from __future__ import annotations

from collections.abc import Mapping

from vntyper.scripts.molecular_identity import MolecularIdentity
from vntyper.scripts.nomenclature_bam_evidence import BamLocusEvidence
from vntyper.scripts.nomenclature_bam_replay import BamReplayArtifact


def validated_whole_locus_bam_evidence(
    artifact: BamReplayArtifact,
    authoritative_identities: Mapping[int, MolecularIdentity | None],
) -> BamLocusEvidence | None:
    """Validate identity bindings and combine producer-deduplicated replay loci.

    Args:
        artifact: Strict retained BAM replay artifact.
        authoritative_identities: Captured candidate identity keyed by ordinal.

    Returns:
        All observed records in replay order, or ``None`` when BAM was never observed.

    Raises:
        ValueError: If artifact, mapping, ordinal custody, or an identity binding differs.
    """
    if not isinstance(artifact, BamReplayArtifact):
        raise ValueError("whole-locus BAM evidence requires a BamReplayArtifact")
    if not isinstance(authoritative_identities, Mapping) or any(
        isinstance(ordinal, bool)
        or not isinstance(ordinal, int)
        or ordinal < 0
        or (identity is not None and not isinstance(identity, MolecularIdentity))
        for ordinal, identity in authoritative_identities.items()
    ):
        raise ValueError("authoritative BAM candidate identities must map non-negative ordinals to identities")
    expected = set(authoritative_identities)
    retained = {ordinal for locus in artifact.loci for ordinal in locus.candidate_observation_ordinals}
    if retained != expected:
        raise ValueError("BAM replay candidate ordinals differ from authoritative candidate identities")
    records = []
    observed = False
    for locus in artifact.loci:
        if locus.state != "observed" or locus.evidence is None:
            continue
        observed = True
        for record in locus.evidence.records:
            for identity, ordinals in zip(record.identities, record.candidate_observation_ordinals, strict=True):
                if any(authoritative_identities.get(ordinal) != identity for ordinal in ordinals):
                    raise ValueError(
                        "BAM identity binding differs from its authoritative captured pre-result candidate identity"
                    )
            records.append(record)
    if not observed:
        return None
    retained_records = tuple(records)
    counts = {
        identity: sum(identity in record.identities for record in retained_records)
        for record in retained_records
        for identity in record.identities
    }
    return BamLocusEvidence(retained_records, len(retained_records), counts)
