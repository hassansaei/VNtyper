"""Unit tests for the CRAM fixture generator.

The generator's whole value is its refusal to record a conversion it has not proved
lossless, so that refusal is what these tests pin. They use a stub samtools rather than
the real binary: the losslessness of ``no_ref=1`` is a property of htslib, measured once
against the real cohort, not something a unit test should re-litigate on every run.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.make_cram_fixtures import (
    CRAM_WRITE_OPTIONS,
    Fixture,
    LossyConversionError,
    Summary,
    derive_cram,
    discover_source_bams,
    write_manifest,
)

pytestmark = pytest.mark.unit


def _touch(path: Path, content: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


class _FakeSamtools:
    """Records the commands asked of it and replays scripted digests.

    ``digests`` is consumed in call order by ``samtools view <file>``: the first entry
    answers for the BAM, the second for the CRAM. Handing it two different values is how
    a test simulates a lossy conversion.
    """

    def __init__(self, digests: list[str]) -> None:
        self.digests = list(digests)
        self.commands: list[list[str]] = []

    def __call__(self, argv: list[str], **kwargs: object) -> object:
        self.commands.append(argv)

        class _Result:
            returncode = 0
            stderr = ""
            stdout = ""

        result = _Result()
        if argv[1] == "view" and "-C" in argv:
            # The conversion itself: create the output file the caller expects.
            Path(argv[argv.index("-o") + 1]).write_text("cram bytes")
        elif argv[1] == "view" and "-c" in argv:
            result.stdout = "7\n"
        return result


def test_discovery_finds_every_bam_and_excludes_the_fixture_root(tmp_path: Path) -> None:
    data_root = tmp_path / "data"
    fixture_root = data_root / "cram"
    _touch(data_root / "a.bam")
    _touch(data_root / "remapped" / "b.bam")
    # A previously derived fixture directory must not become a source for the next run.
    _touch(fixture_root / "a.bam")

    found = discover_source_bams(data_root, fixture_root)

    assert found == [data_root / "a.bam", data_root / "remapped" / "b.bam"]


def test_a_lossy_conversion_raises_rather_than_being_recorded(tmp_path: Path, monkeypatch) -> None:
    """The generator must never emit a fixture whose reads differ from its source.

    This is the property the whole script exists for: a CRAM that silently dropped or
    altered records would make every BAM-versus-CRAM equivalence claim meaningless.
    """
    data_root = tmp_path / "data"
    fixture_root = data_root / "cram"
    bam = _touch(data_root / "sample.bam")

    fake = _FakeSamtools([])
    monkeypatch.setattr("scripts.make_cram_fixtures.subprocess.run", fake)
    # Source and derived digests disagree - exactly what a lossy encoder would produce.
    monkeypatch.setattr(
        "scripts.make_cram_fixtures._record_digest",
        lambda _samtools, alignment: ("aaaa", 7) if alignment == bam else ("bbbb", 7),
    )

    with pytest.raises(LossyConversionError, match="is not lossless"):
        derive_cram("samtools", bam, data_root, fixture_root)


def test_a_faithful_conversion_is_recorded_with_its_evidence(tmp_path: Path, monkeypatch) -> None:
    data_root = tmp_path / "data"
    fixture_root = data_root / "cram"
    bam = _touch(data_root / "remapped" / "sample.bam", "bam bytes")

    fake = _FakeSamtools([])
    monkeypatch.setattr("scripts.make_cram_fixtures.subprocess.run", fake)
    monkeypatch.setattr("scripts.make_cram_fixtures._record_digest", lambda _s, _a: ("deadbeef", 7))

    fixture = derive_cram("samtools", bam, data_root, fixture_root)

    assert fixture.cram == fixture_root / "remapped" / "sample.cram"
    assert fixture.records == 7
    assert fixture.record_digest == "deadbeef"
    # The digest is the evidence; it must be the one shared by both sides.
    assert fixture.unmapped_pairs == 7


def test_the_encoding_options_are_passed_to_samtools(tmp_path: Path, monkeypatch) -> None:
    """``no_ref=1`` is load-bearing, not incidental.

    The cohort's BAM headers carry no ``M5`` tags, so a reference-compressed CRAM could
    not be decoded by a pipeline that passes no ``-T``. If this option is ever dropped the
    fixtures become undecodable by the very code they exist to exercise.
    """
    data_root = tmp_path / "data"
    bam = _touch(data_root / "sample.bam")

    fake = _FakeSamtools([])
    monkeypatch.setattr("scripts.make_cram_fixtures.subprocess.run", fake)
    monkeypatch.setattr("scripts.make_cram_fixtures._record_digest", lambda _s, _a: ("d", 1))

    derive_cram("samtools", bam, data_root, data_root / "cram")

    convert = next(c for c in fake.commands if "-C" in c)
    # Assert the option by name rather than by iterating CRAM_WRITE_OPTIONS: looping over
    # the constant would pass vacuously if the constant were ever emptied, which is the
    # exact regression this test exists to catch.
    assert "no_ref=1" in CRAM_WRITE_OPTIONS, "dropping no_ref makes every fixture need a reference"
    assert "--output-fmt-option" in convert
    assert "no_ref=1" in convert, "no_ref=1 must reach samtools or the fixture needs a reference"


def test_the_manifest_records_what_was_derived_and_what_was_skipped(tmp_path: Path) -> None:
    summary = Summary(
        fixtures=[
            Fixture(
                source_bam=Path("tests/data/s.bam"),
                cram=Path("tests/data/cram/s.cram"),
                records=7,
                unmapped_pairs=2,
                record_digest="abc",
                source_bytes=100,
                cram_bytes=70,
            )
        ],
        skipped=[(Path("tests/data/broken.bam"), "truncated")],
    )
    manifest = tmp_path / "manifest.json"

    write_manifest(summary, manifest)

    payload = json.loads(manifest.read_text())
    assert payload["encoding"] == list(CRAM_WRITE_OPTIONS)
    assert payload["fixtures"][0]["record_digest"] == "abc"
    # A skipped BAM must stay visible: a silently shorter fixture set would weaken the
    # equivalence claim without anyone noticing.
    assert payload["skipped"] == [{"bam": "tests/data/broken.bam", "reason": "truncated"}]
