"""Unit tests for the CRAM fixture generator.

The generator's whole value is its refusal to record a conversion it has not proved
lossless, so that refusal is what these tests pin. They use a stub samtools rather than
the real binary: the losslessness of ``no_ref=1`` is a property of htslib, measured once
against the real cohort, not something a unit test should re-litigate on every run.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any
from unittest import mock

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


def _stub_reference_validation(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cram_fixtures, "validate_reference_fasta", mock.Mock(return_value=mock.sentinel.authority))


def _paired_source_bam(path: Path) -> Path:
    """Write a minimal indexed-coordinate BAM suitable for single-end derivation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as output:
        read = pysam.AlignedSegment()
        read.query_name = "pair-1"
        read.query_sequence = "ACGT"
        read.query_qualities = pysam.qualitystring_to_array("IIII")
        read.flag = 99
        read.reference_id = 0
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigarstring = "4M"
        read.next_reference_id = 0
        read.next_reference_start = 200
        read.template_length = 104
        output.write(read)
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


def test_discovery_returns_an_empty_list_when_no_bams_exist(tmp_path: Path) -> None:
    assert discover_source_bams(tmp_path / "data", tmp_path / "data/cram") == []


def test_run_failure_includes_stderr(monkeypatch: pytest.MonkeyPatch) -> None:
    completed = subprocess.CompletedProcess(["samtools"], 2, stdout="", stderr="broken index")
    monkeypatch.setattr(cram_fixtures.subprocess, "run", mock.Mock(return_value=completed))
    with pytest.raises(RuntimeError, match="samtools exited 2: broken index"):
        cram_fixtures._run(["samtools"])


def test_run_returns_stdout_on_success(monkeypatch: pytest.MonkeyPatch) -> None:
    completed = subprocess.CompletedProcess(["samtools"], 0, stdout="record-count\n", stderr="")
    monkeypatch.setattr(cram_fixtures.subprocess, "run", mock.Mock(return_value=completed))
    assert cram_fixtures._run(["samtools"]) == "record-count\n"


def test_run_to_file_reports_binary_command_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    completed = subprocess.CompletedProcess(["samtools"], 2, stdout=None, stderr=b"broken reheader")
    monkeypatch.setattr(cram_fixtures.subprocess, "run", mock.Mock(return_value=completed))

    with pytest.raises(RuntimeError, match="samtools exited 2: broken reheader"):
        cram_fixtures._run_to_file(["samtools"], tmp_path / "prepared.bam")


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


def test_the_direct_module_dispatches_its_sibling_single_end_builder_without_pythonpath(tmp_path: Path) -> None:
    """Direct-module mode must import and run the sibling builder without package imports."""
    repository_root = tmp_path / "repo"
    data_root = repository_root / "tests" / "data"
    source = _paired_source_bam(data_root / "source.bam")
    output = data_root / "derived" / "source-single.bam"
    config = repository_root / "tests" / "test_data_config.json"
    config.parent.mkdir(parents=True, exist_ok=True)
    config.write_text(
        json.dumps(
            {
                "derived_fixtures": [
                    {
                        "name": "source-single",
                        "kind": "single_end_bam",
                        "source_bam": str(source.relative_to(repository_root)),
                        "output_bam": str(output.relative_to(repository_root)),
                    }
                ]
            }
        )
    )

    dispatch = (
        "import sys; from pathlib import Path; "
        f"sys.path.insert(0, {str(REPO_ROOT / 'scripts')!r}); "
        "import make_cram_fixtures as fixtures; "
        "fixtures._derive_declared_single_end_fixtures(Path('tests/test_data_config.json'), Path('.'))"
    )
    result = subprocess.run(
        [sys.executable, "-c", dispatch],
        cwd=repository_root,
        env={"PATH": os.defpath},
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "ModuleNotFoundError" not in result.stderr
    assert output.is_file()
    assert output.with_suffix(".bam.bai").is_file()


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
    assert fixture.unmapped_reads == 7
    count_command = next(command for command in fake.commands if "-c" in command)
    assert count_command[count_command.index("-f") + 1] == "4"


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
                unmapped_reads=2,
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


def test_reference_compressed_derivation_requires_an_existing_reference(tmp_path: Path) -> None:
    """Removing the selected FASTA must fail before a misleading CRAM is emitted."""
    source = _touch(tmp_path / "tests/data/example_b178_hg19_subset.bam")
    missing_reference = tmp_path / "reference/alignment/chr1.hg19.fa"

    with pytest.raises(FileNotFoundError, match="Reference FASTA does not exist"):
        cram_fixtures.derive_reference_compressed_cram(
            "samtools",
            source,
            missing_reference,
            tmp_path / "tests/data/cram",
        )


def test_wrong_reference_is_rejected_before_htslib_can_fall_back(tmp_path: Path) -> None:
    """A wrong FASTA must not become a valid-looking embedded/non-reference CRAM."""
    source = _touch(tmp_path / "tests/data/example_b178_hg19_subset.bam")
    wrong_reference = _touch(tmp_path / "reference/wrong.fa", ">chr1\nA\n")

    with pytest.raises(ValueError, match="Reference FASTA SHA-256 mismatch"):
        cram_fixtures.derive_reference_compressed_cram(
            "samtools", source, wrong_reference, tmp_path / "tests/data/cram"
        )


def test_reference_compressed_manifest_pins_encoding_source_and_reference_digests(tmp_path: Path) -> None:
    """The manifest must identify both the decoded source stream and exact FASTA bytes."""
    fixture = cram_fixtures.ReferenceCompressedFixture(
        source_bam=Path("tests/data/example_b178_hg19_subset.bam"),
        cram=Path("tests/data/cram/reference-compressed/example_b178_hg19_subset.cram"),
        reference=Path("reference/alignment/chr1.hg19.fa"),
        records=34_214,
        source_record_digest="37428a5a9c95b4791f063c663dc4973548796dcbe9342b870fd4e8425a2d8cc6",
        decoded_record_digest="37428a5a9c95b4791f063c663dc4973548796dcbe9342b870fd4e8425a2d8cc6",
        indexed_region="chr1:155160500-155162000",
        indexed_region_records=13_868,
        source_indexed_region_digest="region-digest",
        decoded_indexed_region_digest="region-digest",
        reference_sha256="0c19925c13b1312f0cbdc2b804f62da260345589b8f9e8ad655abfb5d6e99338",
        source_bytes=3_893_446,
        cram_bytes=2_613_506,
    )
    manifest = tmp_path / "manifest.json"

    write_manifest(Summary(), manifest, reference_compressed=fixture)

    entry = json.loads(manifest.read_text(encoding="utf-8"))["reference_compressed"]
    assert entry == {
        "encoding": "reference-compressed",
        "source_bam": "tests/data/example_b178_hg19_subset.bam",
        "cram": "tests/data/cram/reference-compressed/example_b178_hg19_subset.cram",
        "reference_fasta": "reference/alignment/chr1.hg19.fa",
        "records": 34_214,
        "source_record_digest": "37428a5a9c95b4791f063c663dc4973548796dcbe9342b870fd4e8425a2d8cc6",
        "decoded_record_digest": "37428a5a9c95b4791f063c663dc4973548796dcbe9342b870fd4e8425a2d8cc6",
        "indexed_region": "chr1:155160500-155162000",
        "indexed_region_records": 13_868,
        "source_indexed_region_digest": "region-digest",
        "decoded_indexed_region_digest": "region-digest",
        "reference_sha256": "0c19925c13b1312f0cbdc2b804f62da260345589b8f9e8ad655abfb5d6e99338",
        "source_bytes": 3_893_446,
        "cram_bytes": 2_613_506,
    }


def test_manifest_creation_creates_its_parent_directory(tmp_path: Path) -> None:
    manifest = tmp_path / "nested" / "reports" / "manifest.json"
    write_manifest(Summary(), manifest)
    assert json.loads(manifest.read_text(encoding="utf-8"))["fixtures"] == []


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


def test_select_source_bams_distinguishes_declared_and_all_modes(tmp_path: Path) -> None:
    data_root = tmp_path / "tests/data"
    declared = _touch(data_root / "declared.bam")
    incidental = _touch(data_root / "incidental.bam")
    config = tmp_path / "tests/test_data_config.json"
    config.write_text(
        json.dumps({"file_resources": [{"local_path": "tests/data", "filename": declared.name}]}),
        encoding="utf-8",
    )
    discovered = [declared, incidental]

    assert cram_fixtures._select_source_bams(
        discovered,
        data_config=config,
        data_root=data_root,
        include_all=False,
    ) == [declared]
    assert (
        cram_fixtures._select_source_bams(
            discovered,
            data_config=config,
            data_root=data_root,
            include_all=True,
        )
        == discovered
    )


def test_build_fixtures_retains_successes_when_a_later_conversion_is_skipped(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_root = tmp_path / "tests/data"
    good = _touch(data_root / "a-good.bam", "good")
    bad = _touch(data_root / "b-bad.bam", "bad")
    config = tmp_path / "tests/test_data_config.json"
    config.write_text(
        json.dumps(
            {
                "file_resources": [
                    {"local_path": "tests/data", "filename": good.name},
                    {"local_path": "tests/data", "filename": bad.name},
                ]
            }
        ),
        encoding="utf-8",
    )
    good_fixture = Fixture(good, data_root / "cram/a-good.cram", 7, 2, "digest", 4, 3)

    def derive(_samtools: str, bam: Path, _data_root: Path, _fixture_root: Path) -> Fixture:
        if bam == bad:
            raise RuntimeError("conversion failed")
        return good_fixture

    monkeypatch.setattr(cram_fixtures, "derive_cram", derive)
    summary = cram_fixtures.build_fixtures(
        "samtools",
        data_root,
        data_root / "cram",
        data_config=config,
    )

    assert summary.fixtures == [good_fixture]
    assert summary.skipped == [(bad, "conversion failed")]


def test_a_known_task9_declaration_is_dispatched_without_changing_cram_selection(tmp_path: Path, monkeypatch) -> None:
    """Task9 dispatch is additional to, rather than a replacement for, CRAM selection."""
    repository_root = tmp_path / "repo"
    data_root = repository_root / "tests" / "data"
    declared = _paired_source_bam(data_root / "declared.bam")
    output = data_root / "derived" / "declared-single-end.bam"
    config = repository_root / "tests" / "test_data_config.json"
    config.parent.mkdir(parents=True, exist_ok=True)
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

    cram_fixtures.build_fixtures(
        "samtools",
        data_root,
        data_root / "cram",
        data_config=config,
        repository_root=repository_root,
    )

    assert derived == [declared]
    assert output.is_file()


@pytest.mark.parametrize("include_all", [False, True])
def test_build_fixtures_dispatches_declared_single_end_bams_without_all_bypassing_them(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    include_all: bool,
) -> None:
    """Task9 declarations are materialized by the same command that validates them."""
    repository_root = tmp_path / "repo"
    data_root = repository_root / "tests" / "data"
    source = _paired_source_bam(data_root / "source.bam")
    output = data_root / "derived" / "source-single.bam"
    config = repository_root / "tests" / "test_data_config.json"
    config.parent.mkdir(parents=True, exist_ok=True)
    config.write_text(
        json.dumps(
            {
                "derived_fixtures": [
                    {
                        "name": "source-single",
                        "kind": "single_end_bam",
                        "source_bam": str(source.relative_to(repository_root)),
                        "output_bam": str(output.relative_to(repository_root)),
                    }
                ]
            }
        )
    )
    derived: list[Path] = []
    monkeypatch.setattr(cram_fixtures, "derive_cram", lambda _s, bam, _d, _f: derived.append(bam))

    cram_fixtures.build_fixtures(
        "samtools",
        data_root,
        data_root / "cram",
        data_config=config,
        include_all=include_all,
        repository_root=repository_root,
    )

    assert output.is_file()
    assert output.with_suffix(".bam.bai").is_file()
    with pysam.AlignmentFile(str(output), "rb") as alignment:
        records = list(alignment.fetch(until_eof=True))
    assert len(records) == 1
    assert records[0].is_paired is False
    assert derived == ([source] if include_all else [])


def test_the_ordinary_command_materializes_declared_single_end_outputs_relative_to_the_repo(
    tmp_path: Path,
    monkeypatch,
) -> None:
    repository_root = tmp_path / "repo"
    data_root = repository_root / "tests" / "data"
    source = _paired_source_bam(data_root / "source.bam")
    output = data_root / "derived" / "source-single.bam"
    config = repository_root / "tests" / "test_data_config.json"
    config.parent.mkdir(parents=True, exist_ok=True)
    config.write_text(
        json.dumps(
            {
                "file_resources": [{"local_path": "tests/data", "filename": source.name}],
                "derived_fixtures": [
                    {
                        "name": "source-single",
                        "kind": "single_end_bam",
                        "source_bam": str(source.relative_to(repository_root)),
                        "output_bam": str(output.relative_to(repository_root)),
                    }
                ],
            }
        )
    )

    def fake_derive_cram(_samtools: str, bam: Path, _data_root: Path, fixture_root: Path) -> Fixture:
        return Fixture(bam, fixture_root / "source.cram", 1, 0, "digest", 1, 1)

    monkeypatch.chdir(repository_root)
    monkeypatch.setattr(cram_fixtures, "derive_cram", fake_derive_cram)
    _stub_reference_validation(monkeypatch)
    monkeypatch.setattr(
        cram_fixtures,
        "derive_reference_compressed_cram",
        lambda *_args, **_kwargs: mock.Mock(as_manifest_entry=lambda: {"encoding": "reference-compressed"}),
    )
    monkeypatch.setattr(cram_fixtures, "build_reference_dependent_fixture", lambda _root: None)
    monkeypatch.setattr(cram_fixtures, "build_placed_flag12_fixture", lambda _root: None)
    monkeypatch.setattr(cram_fixtures, "build_indexed_safe_fixture", lambda _root: None)

    exit_code = cram_fixtures.main(
        [
            "--data-root",
            "tests/data",
            "--fixture-root",
            "tests/data/cram",
            "--data-config",
            "tests/test_data_config.json",
        ]
    )

    assert exit_code == 0
    assert output.is_file()
    assert output.with_suffix(".bam.bai").is_file()


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


def test_indexed_safe_fixture_has_nonempty_identical_indexed_and_stream_record_sets(
    tmp_path: Path, monkeypatch
) -> None:
    """HIGH2: column-four zero must authorize a real, nonempty indexed extraction."""
    monkeypatch.setenv("PATH", str(tmp_path / "no-samtools"))
    fixture = cram_fixtures.build_indexed_safe_fixture(tmp_path)

    idxstats = [line.split("\t") for line in pysam_any.idxstats(str(fixture.cram)).splitlines()]
    stream_names = sorted(
        record.partition("\t")[0] for record in pysam_any.view("-f", "4", str(fixture.cram)).splitlines()
    )
    indexed_names = sorted(record.partition("\t")[0] for record in pysam_any.view(str(fixture.cram), "*").splitlines())
    digest = hashlib.sha256("".join(f"{name}\n" for name in stream_names).encode()).hexdigest()

    assert idxstats == [["chr1", "20000", "20", "0"], ["*", "0", "0", "20"]]
    assert len(stream_names) == 20
    assert indexed_names == stream_names
    assert digest == "16a0efa7785630c3d80716d9a386ddaa24f4933b5671f4ecd221b42a8dffe740"


def test_indexed_safe_fixture_mapped_records_are_matching_read_pairs(tmp_path: Path, monkeypatch) -> None:
    """Mapped records must not add a single-end FASTQ beside the unmapped pairs."""
    monkeypatch.setenv("PATH", str(tmp_path / "no-samtools"))
    fixture = cram_fixtures.build_indexed_safe_fixture(tmp_path)

    mapped_records = [
        (fields[0], int(fields[1]))
        for record in pysam_any.view(str(fixture.cram), "chr1").splitlines()
        if (fields := record.split("\t"))
    ]

    assert [name for name, _flag in mapped_records[::2]] == [f"mapped-{number}" for number in range(10)]
    assert mapped_records[::2] == [(f"mapped-{number}", 65) for number in range(10)]
    assert mapped_records[1::2] == [(f"mapped-{number}", 129) for number in range(10)]


def test_the_deriver_command_also_builds_the_purpose_specific_cram_fixtures(tmp_path: Path, monkeypatch) -> None:
    """``make cram-fixtures`` must make the #209, A-SCAN-1 and A-178-2 fixtures available."""
    data_root = tmp_path / "data"
    data_root.mkdir()
    calls: list[Path] = []
    selections: list[bool] = []

    def fake_build(*_args, **kwargs) -> Summary:
        selections.append(kwargs["include_all"])
        return Summary()

    def record_reference_fixture(
        _samtools: str, _source: Path, _reference: Path, root: Path, **_kwargs: object
    ) -> mock.Mock:
        calls.append(root)
        return mock.Mock(as_manifest_entry=lambda: {"encoding": "reference-compressed"})

    monkeypatch.setattr("scripts.make_cram_fixtures.build_fixtures", fake_build)
    _stub_reference_validation(monkeypatch)
    monkeypatch.setattr(
        "scripts.make_cram_fixtures.derive_reference_compressed_cram",
        record_reference_fixture,
    )
    monkeypatch.setattr("scripts.make_cram_fixtures.build_reference_dependent_fixture", lambda root: calls.append(root))
    monkeypatch.setattr("scripts.make_cram_fixtures.build_placed_flag12_fixture", lambda root: calls.append(root))
    monkeypatch.setattr("scripts.make_cram_fixtures.build_indexed_safe_fixture", lambda root: calls.append(root))

    exit_code = cram_fixtures.main(["--data-root", str(data_root), "--fixture-root", str(tmp_path / "cram")])
    all_exit_code = cram_fixtures.main(
        ["--data-root", str(data_root), "--fixture-root", str(tmp_path / "all-cram"), "--all"]
    )

    assert exit_code == 1
    assert all_exit_code == 1
    assert selections == [False, True]
    assert calls == [
        tmp_path / "cram",
        tmp_path / "cram",
        tmp_path / "cram",
        tmp_path / "cram",
        tmp_path / "all-cram",
        tmp_path / "all-cram",
        tmp_path / "all-cram",
        tmp_path / "all-cram",
    ]


def test_the_deriver_command_fails_when_any_declared_cram_was_skipped(tmp_path: Path, monkeypatch) -> None:
    """A partial fixture tree must not make the ordinary Make target report success."""
    data_root = tmp_path / "data"
    data_root.mkdir()
    source = data_root / "good.bam"
    source.write_bytes(b"BAM")
    fixture = cram_fixtures.Fixture(
        source_bam=source,
        cram=tmp_path / "cram" / "good.cram",
        records=1,
        unmapped_reads=0,
        record_digest="digest",
        source_bytes=3,
        cram_bytes=2,
    )
    summary = Summary(fixtures=[fixture], skipped=[(data_root / "broken.bam", "truncated")])

    monkeypatch.setattr("scripts.make_cram_fixtures.build_fixtures", lambda *_args, **_kwargs: summary)
    _stub_reference_validation(monkeypatch)
    monkeypatch.setattr(
        "scripts.make_cram_fixtures.derive_reference_compressed_cram",
        lambda *_args, **_kwargs: mock.Mock(as_manifest_entry=lambda: {"encoding": "reference-compressed"}),
    )
    monkeypatch.setattr("scripts.make_cram_fixtures.build_reference_dependent_fixture", lambda _root: None)
    monkeypatch.setattr("scripts.make_cram_fixtures.build_placed_flag12_fixture", lambda _root: None)
    monkeypatch.setattr("scripts.make_cram_fixtures.build_indexed_safe_fixture", lambda _root: None)

    exit_code = cram_fixtures.main(["--data-root", str(data_root), "--fixture-root", str(tmp_path / "cram")])

    assert exit_code == 1


def test_lossy_build_aborts_before_a_manifest_can_record_the_failed_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_root = tmp_path / "data"
    data_root.mkdir()
    config = tmp_path / "test_data_config.json"
    config.write_text("{}", encoding="utf-8")
    manifest = tmp_path / "manifest.json"

    def fail_build(*args: object, **kwargs: object) -> Summary:
        del args, kwargs
        raise LossyConversionError("not lossless")

    monkeypatch.setattr(cram_fixtures, "build_fixtures", fail_build)
    _stub_reference_validation(monkeypatch)
    with pytest.raises(LossyConversionError, match="not lossless"):
        cram_fixtures.main(
            [
                "--data-root",
                str(data_root),
                "--data-config",
                str(config),
                "--manifest",
                str(manifest),
            ]
        )
    assert not manifest.exists()


@pytest.mark.parametrize("missing", ["data", "config"])
def test_main_returns_one_when_required_input_is_missing(tmp_path: Path, missing: str) -> None:
    data_root = tmp_path / "data"
    config = tmp_path / "test_data_config.json"
    if missing != "data":
        data_root.mkdir()
    if missing != "config":
        config.write_text("{}", encoding="utf-8")

    assert cram_fixtures.main(["--data-root", str(data_root), "--data-config", str(config)]) == 1


def test_the_reference_contract_purpose_fixtures_are_registered_with_portable_paths() -> None:
    payload = json.loads((REPO_ROOT / "tests" / "test_data_config.json").read_text())
    fixtures = payload["purpose_fixtures"]["cram_reference_contract"]

    assert fixtures == {
        "reference_dependent_cram": "tests/data/cram/reference-dependent/reference-dependent.cram",
        "reference_fasta": "tests/data/cram/reference-dependent/reference.fa",
        "no_ref_cram": "tests/data/cram/placed-flag12/placed-flag12.cram",
        "indexed_safe_cram": "tests/data/cram/indexed-safe/indexed-safe.cram",
    }
    assert all(not Path(path).is_absolute() and ".." not in Path(path).parts for path in fixtures.values())
