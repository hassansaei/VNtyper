"""A supported pipeline run can emit replayable Kestrel calibration evidence.

Without this, calibrating a new cohort would mean patching the pipeline at runtime,
which is neither reproducible nor reviewable. These tests pin the writer's contract:
it emits exactly one complete capture per run, it never writes a partial file, and an
empty candidate population is recorded as evidence of a completed run rather than as
a missing file.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from tests.unit.test_calibration_caller_policy import policy_document
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
from vntyper.scripts.calibration_kestrel_capture import KESTREL_RAW_COLUMNS, decode_kestrel_capture
from vntyper.scripts.identity_candidates import translation_component_from_config
from vntyper.scripts.kestrel_genotyping import _resolve_selection
from vntyper.scripts.nomenclature import nomenclature_config
from vntyper.scripts.pipeline_kestrel_capture import (
    CaptureAssets,
    build_capture_policy,
    kestrel_capture_writer,
)

pytestmark = pytest.mark.unit

_SHA = "3" * 64


def _baseline():
    document = policy_document(include_advntr=False)
    values = document["values"]
    assert isinstance(values, dict)
    values["/components/kestrel/confidence_assignment/reporting_floor"] = 0.00469
    return decode_caller_policy_values(document)


def _assets(tmp_path: Path) -> CaptureAssets:
    reference = tmp_path / "reference.fa"
    motifs = tmp_path / "motifs.fa"
    jar = tmp_path / "kestrel.jar"
    reference.write_bytes(b">X-5\nSEQ1SEQ2\n")
    motifs.write_bytes(b">X\nSEQ1\n>5\nSEQ2\n")
    jar.write_bytes(b"invented-kestrel-jar")
    return CaptureAssets(
        reference_file=reference,
        motif_reference_file=motifs,
        kestrel_jar=jar,
        baseline_policy=_baseline(),
        decision_profile_sha256=_SHA,
        capture_policy=build_capture_policy(kmer_sizes=(20,), threads=2, input_scope="bam-fast-mode"),
    )


def _motifs() -> pd.DataFrame:
    return pd.DataFrame({"Motif": ["5", "X"], "Motif_sequence": ["SEQ2", "SEQ1"]})


def _observe(writer, frame: pd.DataFrame) -> None:
    writer(
        frame,
        _motifs(),
        kestrel_config(),
        _resolve_selection(kestrel_config()),
        translation_component_from_config(nomenclature_config),
    )


def test_writer_emits_one_decodable_capture_of_the_complete_raw_population(tmp_path: Path) -> None:
    destination = tmp_path / "kestrel_capture.json"
    writer = kestrel_capture_writer(destination, assets=_assets(tmp_path))
    frame = kestrel_stage_frame("raw", rows=3)

    _observe(writer, frame)

    document = json.loads(destination.read_text())
    capture = decode_kestrel_capture(document)
    assert len(capture.rows) == 3
    assert document["source_row_count"] == 3
    assert [row.source_row_ordinal for row in capture.rows] == [0, 1, 2]
    assert destination.stat().st_mode & 0o777 == 0o600


def test_a_second_observation_in_one_run_is_refused(tmp_path: Path) -> None:
    destination = tmp_path / "kestrel_capture.json"
    writer = kestrel_capture_writer(destination, assets=_assets(tmp_path))
    frame = kestrel_stage_frame("raw")
    _observe(writer, frame)

    with pytest.raises(ValueError, match="once"):
        _observe(writer, frame)


def test_an_empty_population_is_captured_as_a_completed_run(tmp_path: Path) -> None:
    """A run that found no indel still produced evidence: a complete empty population.

    Treating the absent file as "no evidence" would silently drop every true negative
    with no candidates from a calibration cohort.
    """
    destination = tmp_path / "kestrel_capture.json"
    writer = kestrel_capture_writer(destination, assets=_assets(tmp_path))

    writer.write_empty(_motifs(), kestrel_config())

    document = json.loads(destination.read_text())
    capture = decode_kestrel_capture(document)
    assert capture.rows == ()
    assert document["source_row_count"] == 0


def test_empty_capture_is_refused_once_a_population_was_already_written(tmp_path: Path) -> None:
    destination = tmp_path / "kestrel_capture.json"
    writer = kestrel_capture_writer(destination, assets=_assets(tmp_path))
    _observe(writer, kestrel_stage_frame("raw"))

    with pytest.raises(ValueError, match="once"):
        writer.write_empty(_motifs(), kestrel_config())


def test_raw_columns_match_the_capture_contract(tmp_path: Path) -> None:
    destination = tmp_path / "kestrel_capture.json"
    writer = kestrel_capture_writer(destination, assets=_assets(tmp_path))

    _observe(writer, kestrel_stage_frame("raw"))

    document = json.loads(destination.read_text())
    assert document["raw_columns"] == list(KESTREL_RAW_COLUMNS)
    assert document["population"] == "complete-prefilter-v1"


def test_capture_policy_records_the_actual_run_parameters() -> None:
    policy = build_capture_policy(kmer_sizes=(20, 25), threads=4, input_scope="fastq")

    assert policy["kmer_sizes"] == [20, 25]
    assert policy["threads"] == 4
    assert policy["input_scope"] == "fastq"
    assert policy["schema_version"] == "vntyper-kestrel-recruitment-v1"
    assert policy["vntyper_version"]


@pytest.mark.parametrize(
    ("kmer_sizes", "threads", "input_scope"),
    [((), 2, "fastq"), ((20,), 0, "fastq"), ((20,), 2, ""), ((20, 20), 2, "fastq"), ((0,), 2, "fastq")],
)
def test_capture_policy_refuses_unusable_run_parameters(kmer_sizes, threads, input_scope) -> None:
    with pytest.raises(ValueError):
        build_capture_policy(kmer_sizes=kmer_sizes, threads=threads, input_scope=input_scope)


def test_a_failed_write_leaves_no_partial_capture(tmp_path: Path) -> None:
    destination = tmp_path / "nested" / "kestrel_capture.json"
    writer = kestrel_capture_writer(destination, assets=_assets(tmp_path))

    with pytest.raises(ValueError, match="motif"):
        writer(
            kestrel_stage_frame("raw"),
            pd.DataFrame({"wrong": [1]}),
            kestrel_config(),
            _resolve_selection(kestrel_config()),
            translation_component_from_config(nomenclature_config),
        )

    assert not destination.exists()


def test_missing_asset_files_are_reported_before_any_capture_is_built(tmp_path: Path) -> None:
    assets = _assets(tmp_path)
    assets.kestrel_jar.unlink()
    destination = tmp_path / "kestrel_capture.json"
    writer = kestrel_capture_writer(destination, assets=assets)

    with pytest.raises(ValueError, match="cannot read"):
        _observe(writer, kestrel_stage_frame("raw"))

    assert not destination.exists()


def test_a_non_writable_destination_is_reported_and_leaves_nothing(tmp_path: Path) -> None:
    """An unwritable destination must fail loudly rather than lose the evidence."""
    blocked = tmp_path / "blocked"
    blocked.mkdir()
    blocked.chmod(0o500)
    try:
        writer = kestrel_capture_writer(blocked / "kestrel_capture.json", assets=_assets(tmp_path))
        with pytest.raises(ValueError, match="cannot write"):
            _observe(writer, kestrel_stage_frame("raw"))
    finally:
        blocked.chmod(0o700)
    assert list(blocked.iterdir()) == []


def test_the_writer_factory_refuses_unvalidated_assets(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="CaptureAssets"):
        kestrel_capture_writer(tmp_path / "capture.json", assets={"reference_file": "x"})
