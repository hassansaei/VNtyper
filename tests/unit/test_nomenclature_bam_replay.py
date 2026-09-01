"""Closed codec and atomic sidecar for Kestrel BAM replay evidence."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, cast
from unittest import mock

import pytest

from vntyper.scripts.molecular_identity import make_coding_edit, make_molecular_identity, serialize_molecular_identity
from vntyper.scripts.nomenclature_bam_evidence import BamIdentityEvidence, BamLocusEvidence
from vntyper.scripts.nomenclature_bam_replay import (
    BAM_REPLAY_FILENAME,
    BAM_REPLAY_SCHEMA_VERSION,
    BamReplayArtifact,
    BamReplayLocus,
    decode_bam_replay_artifact,
    encode_bam_replay_artifact,
    invalidate_bam_replay_artifact,
    read_bam_replay_artifact,
    write_bam_replay_artifact,
)

pytestmark = pytest.mark.unit

_DUPC = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
_DUPA = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))


def _winner_runner_up() -> BamLocusEvidence:
    records = (
        BamIdentityEvidence((_DUPC,), ((2,),), 5),
        BamIdentityEvidence((_DUPC, _DUPA), ((2,), (7,)), None),
        BamIdentityEvidence((_DUPA,), ((7,),), 2_147_483_647),
        BamIdentityEvidence((), (), 0),
        BamIdentityEvidence((_DUPC,), ((2,),), 181),
    )
    return BamLocusEvidence(records, 5, {_DUPC: 3, _DUPA: 2})


def _tie() -> BamLocusEvidence:
    records = _winner_runner_up().records[:-1]
    return BamLocusEvidence(records, 4, {_DUPC: 2, _DUPA: 2})


def _artifact() -> BamReplayArtifact:
    return BamReplayArtifact(
        (
            BamReplayLocus(2, "observed", _winner_runner_up()),
            BamReplayLocus(7, "observed", _tie()),
            BamReplayLocus(11, "observed", BamLocusEvidence((), 0, {})),
            BamReplayLocus(13, "unavailable", None),
            BamReplayLocus(17, "not-consulted", None),
        )
    )


def test_replay_codec_uses_closed_canonical_primitive_schema() -> None:
    encoded = encode_bam_replay_artifact(_artifact())

    assert encoded == {
        "schema_version": "bam-identity-replay-v1",
        "loci": [
            {
                "observation_ordinal": 2,
                "state": "observed",
                "evidence": {
                    "eligible_record_count": 5,
                    "records": [
                        {
                            "identities": [serialize_molecular_identity(_DUPC)],
                            "candidate_observation_ordinals": [[2]],
                            "minimum_kmer_depth": 5,
                        },
                        {
                            "identities": [
                                serialize_molecular_identity(_DUPC),
                                serialize_molecular_identity(_DUPA),
                            ],
                            "candidate_observation_ordinals": [[2], [7]],
                            "minimum_kmer_depth": None,
                        },
                        {
                            "identities": [serialize_molecular_identity(_DUPA)],
                            "candidate_observation_ordinals": [[7]],
                            "minimum_kmer_depth": 2_147_483_647,
                        },
                        {
                            "identities": [],
                            "candidate_observation_ordinals": [],
                            "minimum_kmer_depth": 0,
                        },
                        {
                            "identities": [serialize_molecular_identity(_DUPC)],
                            "candidate_observation_ordinals": [[2]],
                            "minimum_kmer_depth": 181,
                        },
                    ],
                    "counts": [
                        {"identity": serialize_molecular_identity(_DUPC), "record_count": 3},
                        {"identity": serialize_molecular_identity(_DUPA), "record_count": 2},
                    ],
                },
            },
            {
                "observation_ordinal": 7,
                "state": "observed",
                "evidence": {
                    "eligible_record_count": 4,
                    "records": [
                        {
                            "identities": [serialize_molecular_identity(_DUPC)],
                            "candidate_observation_ordinals": [[2]],
                            "minimum_kmer_depth": 5,
                        },
                        {
                            "identities": [
                                serialize_molecular_identity(_DUPC),
                                serialize_molecular_identity(_DUPA),
                            ],
                            "candidate_observation_ordinals": [[2], [7]],
                            "minimum_kmer_depth": None,
                        },
                        {
                            "identities": [serialize_molecular_identity(_DUPA)],
                            "candidate_observation_ordinals": [[7]],
                            "minimum_kmer_depth": 2_147_483_647,
                        },
                        {
                            "identities": [],
                            "candidate_observation_ordinals": [],
                            "minimum_kmer_depth": 0,
                        },
                    ],
                    "counts": [
                        {"identity": serialize_molecular_identity(_DUPC), "record_count": 2},
                        {"identity": serialize_molecular_identity(_DUPA), "record_count": 2},
                    ],
                },
            },
            {
                "observation_ordinal": 11,
                "state": "observed",
                "evidence": {"eligible_record_count": 0, "records": [], "counts": []},
            },
            {"observation_ordinal": 13, "state": "unavailable", "evidence": None},
            {"observation_ordinal": 17, "state": "not-consulted", "evidence": None},
        ],
    }
    assert BAM_REPLAY_SCHEMA_VERSION == "bam-identity-replay-v1"
    assert BAM_REPLAY_FILENAME == "bam_identity_replay.v1.json"
    json.dumps(encoded)


def test_replay_codec_round_trips_empty_winner_runner_up_tie_and_states() -> None:
    artifact = _artifact()

    decoded = decode_bam_replay_artifact(encode_bam_replay_artifact(artifact))

    assert decoded == artifact
    assert decoded.loci[0].evidence is not None
    assert decoded.loci[0].evidence.record_identity_sets == (
        (_DUPC,),
        (_DUPC, _DUPA),
        (_DUPA,),
        (),
        (_DUPC,),
    )
    assert decoded.loci[0].evidence.winning_identity == _DUPC
    assert decoded.loci[0].evidence.counts[_DUPA] == 2
    assert decoded.loci[1].evidence is not None
    assert decoded.loci[1].evidence.winning_identity is None
    assert decoded.loci[2].evidence == BamLocusEvidence((), 0, {})


@pytest.mark.parametrize(
    "factory",
    [
        lambda: BamReplayLocus(True, "not-consulted", None),
        lambda: BamReplayLocus(0, "unknown", None),  # type: ignore[arg-type]
        lambda: BamReplayLocus(0, "observed", None),
        lambda: BamReplayLocus(0, "unavailable", BamLocusEvidence((), 0, {})),
        lambda: BamReplayArtifact([]),  # type: ignore[arg-type]
        lambda: BamReplayArtifact((BamReplayLocus(2, "unavailable", None), BamReplayLocus(1, "unavailable", None))),
        lambda: BamReplayArtifact((BamReplayLocus(2, "unavailable", None), BamReplayLocus(2, "unavailable", None))),
    ],
)
def test_replay_values_reject_invalid_state_pairing_and_ordinals(factory) -> None:
    with pytest.raises(ValueError):
        factory()


@pytest.mark.parametrize(
    "mutate",
    [
        lambda value: value.update({"extra": True}),
        lambda value: value["loci"][0].update({"extra": True}),
        lambda value: value["loci"][0]["evidence"].update({"extra": True}),
        lambda value: value["loci"][0]["evidence"]["records"][0].update({"extra": True}),
        lambda value: value["loci"][0]["evidence"]["counts"][0].update({"extra": True}),
        lambda value: value["loci"][0]["evidence"]["counts"][0].update({"record_count": 4}),
        lambda value: value.update({"schema_version": "bam-identity-replay-v2"}),
    ],
)
def test_replay_decoder_rejects_extra_keys_versions_and_inconsistent_derived_counts(mutate) -> None:
    encoded = cast(dict[str, Any], encode_bam_replay_artifact(_artifact()))
    mutate(encoded)

    with pytest.raises(ValueError):
        decode_bam_replay_artifact(encoded)


def test_replay_decoder_rejects_noncanonical_identity_and_count_order() -> None:
    encoded = cast(dict[str, Any], encode_bam_replay_artifact(_artifact()))
    evidence = encoded["loci"][0]["evidence"]
    evidence["records"][0]["identities"][0] += "|"
    with pytest.raises(ValueError):
        decode_bam_replay_artifact(encoded)

    reordered = cast(dict[str, Any], encode_bam_replay_artifact(_artifact()))
    reordered["loci"][0]["evidence"]["counts"].reverse()
    with pytest.raises(ValueError, match="canonical"):
        decode_bam_replay_artifact(reordered)


def test_atomic_replay_write_round_trips_and_uses_sibling_replace(tmp_path: Path) -> None:
    artifact = _artifact()
    real_replace = os.replace
    with mock.patch("vntyper.scripts.nomenclature_bam_replay.os.replace", wraps=real_replace) as replaced:
        path = write_bam_replay_artifact(tmp_path, artifact)

    assert path == tmp_path / BAM_REPLAY_FILENAME
    assert read_bam_replay_artifact(tmp_path) == artifact
    source, destination = replaced.call_args.args
    assert Path(source).parent == tmp_path
    assert Path(source).name == BAM_REPLAY_FILENAME + ".tmp"
    assert Path(destination) == path
    assert not Path(source).exists()


def test_replay_reader_fails_closed_on_missing_invalid_and_duplicate_keys(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_bam_replay_artifact(tmp_path)

    path = tmp_path / BAM_REPLAY_FILENAME
    path.write_text("not json", encoding="utf-8")
    with pytest.raises(ValueError, match="valid JSON"):
        read_bam_replay_artifact(tmp_path)

    path.write_text(
        '{"schema_version":"bam-identity-replay-v1","schema_version":"duplicate","loci":[]}', encoding="utf-8"
    )
    with pytest.raises(ValueError, match="duplicate"):
        read_bam_replay_artifact(tmp_path)


def test_failed_atomic_install_preserves_the_previous_replay_artifact(tmp_path: Path) -> None:
    path = write_bam_replay_artifact(tmp_path, _artifact())
    previous = path.read_bytes()

    with (
        mock.patch("vntyper.scripts.nomenclature_bam_replay.os.replace", side_effect=OSError("install failed")),
        pytest.raises(OSError, match="install failed"),
    ):
        write_bam_replay_artifact(tmp_path, BamReplayArtifact(()))

    assert path.read_bytes() == previous


def test_replay_invalidation_removes_only_the_exact_stale_sidecar(tmp_path: Path) -> None:
    path = write_bam_replay_artifact(tmp_path, _artifact())
    neighbour = tmp_path / "keep.json"
    neighbour.write_text("keep", encoding="utf-8")

    invalidate_bam_replay_artifact(tmp_path)
    invalidate_bam_replay_artifact(tmp_path)

    assert not path.exists()
    assert neighbour.read_text(encoding="utf-8") == "keep"
