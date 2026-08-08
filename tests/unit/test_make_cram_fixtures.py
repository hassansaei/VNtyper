"""Unit tests for the CRAM fixture generator.

The generator's whole value is its refusal to record a conversion it has not proved
lossless, so that refusal is what these tests pin. They use a stub samtools rather than
the real binary: the losslessness of ``no_ref=1`` is a property of htslib, measured once
against the real cohort, not something a unit test should re-litigate on every run.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pysam
import pytest

import scripts.make_cram_fixtures as cram_fixtures
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

pysam_any: Any = pysam
REPO_ROOT = Path(__file__).resolve().parents[2]


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


def test_the_direct_script_entry_point_loads_its_sibling_selection_module_without_pythonpath(tmp_path: Path) -> None:
    """A direct invocation cannot rely on importing the repository's ``scripts`` package."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "make_cram_fixtures.py"), "--help"],
        cwd=REPO_ROOT,
        env={"PATH": os.defpath},
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--all" in result.stdout


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


def test_the_default_deriver_uses_only_bams_declared_in_the_data_config(tmp_path: Path, monkeypatch) -> None:
    """An incidental BAM must not silently enlarge the normal fixture set."""
    data_root = tmp_path / "data"
    declared = _touch(data_root / "declared.bam")
    incidental = _touch(data_root / "incidental.bam")
    config = tmp_path / "test_data_config.json"
    config.write_text(json.dumps({"file_resources": [{"local_path": "tests/data", "filename": declared.name}]}))
    derived: list[Path] = []
    monkeypatch.setattr(
        "scripts.make_cram_fixtures.derive_cram",
        lambda _s, bam, _d, _f: derived.append(bam),
    )

    cram_fixtures.build_fixtures("samtools", data_root, data_root / "cram", data_config=config, include_all=False)

    assert derived == [declared]
    assert incidental not in derived


def test_the_all_switch_derives_every_discovered_bam(tmp_path: Path, monkeypatch) -> None:
    """The gate can intentionally request all cohort BAMs, including new discoveries."""
    data_root = tmp_path / "data"
    declared = _touch(data_root / "declared.bam")
    incidental = _touch(data_root / "incidental.bam")
    config = tmp_path / "test_data_config.json"
    config.write_text(json.dumps({"file_resources": [{"local_path": "tests/data", "filename": declared.name}]}))
    derived: list[Path] = []
    monkeypatch.setattr(
        "scripts.make_cram_fixtures.derive_cram",
        lambda _s, bam, _d, _f: derived.append(bam),
    )

    cram_fixtures.build_fixtures("samtools", data_root, data_root / "cram", data_config=config, include_all=True)

    assert derived == [declared, incidental]


def test_a_known_task9_derived_fixture_declaration_is_validated_without_changing_cram_selection(
    tmp_path: Path, monkeypatch
) -> None:
    """Task9's schema is recognized, but its single-end builder stays independently owned."""
    data_root = tmp_path / "data"
    declared = _touch(data_root / "declared.bam")
    config = tmp_path / "test_data_config.json"
    config.write_text(
        json.dumps(
            {
                "file_resources": [{"local_path": "tests/data", "filename": declared.name}],
                "derived_fixtures": [
                    {
                        "name": "declared-single-end",
                        "kind": "single_end_bam",
                        "source_bam": "tests/data/declared.bam",
                        "output_bam": "tests/data/derived/declared-single-end.bam",
                    }
                ],
            }
        )
    )
    derived: list[Path] = []
    monkeypatch.setattr("scripts.make_cram_fixtures.derive_cram", lambda _s, bam, _d, _f: derived.append(bam))

    cram_fixtures.build_fixtures("samtools", data_root, data_root / "cram", data_config=config)

    assert derived == [declared]


def test_an_unknown_declared_fixture_kind_is_refused_instead_of_being_ignored(tmp_path: Path) -> None:
    """A new derived-fixture kind must get an explicit dispatch decision."""
    data_root = tmp_path / "data"
    config = tmp_path / "test_data_config.json"
    config.write_text(json.dumps({"derived_fixtures": [{"name": "unknown", "kind": "new_kind"}]}))

    with pytest.raises(ValueError, match="unsupported derived fixture kind 'new_kind'"):
        cram_fixtures.build_fixtures("samtools", data_root, data_root / "cram", data_config=config)


def test_all_mode_still_refuses_an_unknown_declared_fixture_kind(tmp_path: Path) -> None:
    """Selection breadth must not bypass declaration validation."""
    data_root = tmp_path / "data"
    config = tmp_path / "test_data_config.json"
    config.write_text(json.dumps({"derived_fixtures": [{"name": "unknown", "kind": "new_kind"}]}))

    with pytest.raises(ValueError, match="unsupported derived fixture kind 'new_kind'"):
        cram_fixtures.build_fixtures("samtools", data_root, data_root / "cram", data_config=config, include_all=True)


def test_a_malformed_known_derived_fixture_declaration_is_refused(tmp_path: Path) -> None:
    """Recognizing Task9's kind includes checking the fields its later dispatcher needs."""
    data_root = tmp_path / "data"
    config = tmp_path / "test_data_config.json"
    config.write_text(
        json.dumps({"derived_fixtures": [{"name": "single", "kind": "single_end_bam", "source_bam": "x.bam"}]})
    )

    with pytest.raises(ValueError, match="single_end_bam.*output_bam"):
        cram_fixtures.build_fixtures("samtools", data_root, data_root / "cram", data_config=config)


def test_reference_dependent_fixture_has_a_local_ur_target_that_can_be_removed(tmp_path: Path, monkeypatch) -> None:
    """A missing-reference test is valid only after its header ``UR:`` target is gone."""
    monkeypatch.setenv("PATH", str(tmp_path / "no-samtools"))
    fixture = cram_fixtures.build_reference_dependent_fixture(tmp_path)
    with pysam.AlignmentFile(str(fixture.cram), "rc", reference_filename=str(fixture.reference)) as alignment:
        sequence = alignment.header.to_dict()["SQ"][0]

    assert sequence["M5"]
    assert Path(sequence["UR"]) == fixture.reference
    missing = fixture.reference.with_name("reference-is-missing.fa")
    fixture.reference.rename(missing)
    local_only_ref_path = str(tmp_path / "no-reference-cache" / "%2s" / "%2s" / "%s")
    real_view = pysam_any.view
    monkeypatch.setenv("REF_PATH", local_only_ref_path)

    def guarded_view(*args: str) -> str:
        assert os.environ.get("REF_PATH") == local_only_ref_path
        return real_view(*args)

    monkeypatch.setattr(pysam_any, "view", guarded_view)
    with pytest.raises(pysam.SamtoolsError):
        pysam_any.view("-h", str(fixture.cram))


def test_placed_flag12_fixture_proves_idxstats_column_four_requires_the_stream_scan(
    tmp_path: Path, monkeypatch
) -> None:
    """The indexed ``'*'`` query loses precisely the placed flag-12 records."""
    monkeypatch.setenv("PATH", str(tmp_path / "no-samtools"))
    fixture = cram_fixtures.build_placed_flag12_fixture(tmp_path)

    assert int(pysam_any.view("-c", "-f", "12", str(fixture.cram))) == 130
    assert int(pysam_any.view("-c", "-f", "12", str(fixture.cram), "*")) == 80
    fields = pysam_any.idxstats(str(fixture.cram)).splitlines()[0].split("\t")
    assert fields == ["chr1", "20000", "600", "50"]


def test_the_deriver_command_also_builds_the_purpose_specific_cram_fixtures(tmp_path: Path, monkeypatch) -> None:
    """``make cram-fixtures`` must make the #209 and A-SCAN-1 fixtures available."""
    data_root = tmp_path / "data"
    data_root.mkdir()
    calls: list[Path] = []
    selections: list[bool] = []

    def fake_build(*_args, **kwargs) -> Summary:
        selections.append(kwargs["include_all"])
        return Summary()

    monkeypatch.setattr("scripts.make_cram_fixtures.build_fixtures", fake_build)
    monkeypatch.setattr("scripts.make_cram_fixtures.build_reference_dependent_fixture", lambda root: calls.append(root))
    monkeypatch.setattr("scripts.make_cram_fixtures.build_placed_flag12_fixture", lambda root: calls.append(root))

    exit_code = cram_fixtures.main(["--data-root", str(data_root), "--fixture-root", str(tmp_path / "cram")])
    all_exit_code = cram_fixtures.main(
        ["--data-root", str(data_root), "--fixture-root", str(tmp_path / "all-cram"), "--all"]
    )

    assert exit_code == 1
    assert all_exit_code == 1
    assert selections == [False, True]
    assert calls == [tmp_path / "cram", tmp_path / "cram", tmp_path / "all-cram", tmp_path / "all-cram"]
