"""What goes into an immutable reference release, and what is kept out of it.

`refs-v1` is built once, published as a GitHub release and never edited. A bundle that
omits a file, ships one twice, or records a digest that does not describe the bytes it
names cannot be patched - it forces a new release version, a full rebuild and new
digests in every pin that references it.

The load-bearing property tested here is therefore not "the seven assets exist". It is
that **every file under the reference tree is assigned to exactly one asset or matches
an explicit exclusion, and anything else fails the build by name**. Two omissions -
`filter_config.json` never being staged, and `vntr_db_advntr.zip` being left unassigned -
were caught by hand before this rule existed; the rule is what makes catching them
mechanical.

The tests build real tar archives over real (tiny) files and re-hash what comes back out
of them, rather than asserting against mocks: a manifest that agrees with a mock but not
with the archive is exactly the failure that cannot be fixed after publication.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import sys
import tarfile
from pathlib import Path
from typing import Any

import pytest

from vntyper.scripts import reference_bundle, verify_seeds
from vntyper.scripts.reference_bundle import safe_extract

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import bundle_release  # noqa: E402

#: Captured before any test patches the module attribute, so the one test that wants the
#: real subprocess probe gets it back rather than the pinned stand-in.
_REAL_CAPTURE = bundle_release._capture

#: Physical reference id -> the extracted chromosome FASTA `--from-source` leaves behind.
#: NCBI ships `.fna`, everyone else `.fa`; both spellings have to survive the grouper.
GENOMES = {
    "hg19": "alignment/chr1.hg19.fa",
    "hg38": "alignment/chr1.hg38.fa",
    "GRCh37": "alignment/chr1.GRCh37.fna",
    "GRCh38": "alignment/chr1.GRCh38.fna",
    "hg19_ensembl": "alignment/chr1.hg19_ensembl.fa",
    "hg38_ensembl": "alignment/chr1.hg38_ensembl.fa",
}

#: `samtools faidx` plus the five files `bwa index` writes beside a FASTA.
SIDECARS = (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa")

#: Every MUC1-scoped FASTA in the tree. The first three are the "three MUC1 FASTAs";
#: the last two are the SHARK region cuts.
MUC1_FASTAS = (
    "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
    "MUC1_motifs_Rev_com.fa",
    "code-adVNTR_RUs.fa",
    "muc1_region_hg19.fa",
    "muc1_region_hg38.fa",
)

#: Outputs `run_derivations` produces and verifies against a committed digest.
DERIVED = (
    "muc1_region_hg19.fa",
    "muc1_region_hg38.fa",
    "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
)

#: The four non-derivable artefacts the data repository commits.
SEEDS = ("MUC1_motifs_Rev_com.fa", "code-adVNTR_RUs.fa", "vntr_db_advntr.zip", "filter_config.json")

TAG = "refs-v1"

EXPECTED_ASSETS = {
    "vntyper-references-refs-v1-ucsc-hg19.tar.gz",
    "vntyper-references-refs-v1-ucsc-hg38.tar.gz",
    "vntyper-references-refs-v1-ncbi-GRCh37.tar.gz",
    "vntyper-references-refs-v1-ncbi-GRCh38.tar.gz",
    "vntyper-references-refs-v1-ensembl-hg19.tar.gz",
    "vntyper-references-refs-v1-ensembl-hg38.tar.gz",
    "vntyper-references-refs-v1-muc1.tar.gz",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _tree_files(root: Path) -> set[str]:
    return {str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}


@pytest.fixture
def refs(tmp_path: Path) -> Path:
    """A miniature of the tree `install-references --from-source` actually produces.

    The inventory is taken from a real six-genome run: the extracted FASTA and its six
    index sidecars per assembly, the downloaded `.gz` beside them, five MUC1 FASTAs with
    their `.fai`, the two adVNTR databases and the zip they came out of, the staged
    `filter_config.json`, and the installer's own log.
    """
    root = tmp_path / "refs"
    for ref_id, fasta in GENOMES.items():
        payload = f">chr1 {ref_id}\nACGTACGT\n".encode()
        _write(root / fasta, payload)
        _write(root / f"{fasta}.gz", gzip.compress(payload))
        for suffix in SIDECARS:
            _write(root / f"{fasta}{suffix}", f"{ref_id}{suffix} index bytes".encode())
    for name in MUC1_FASTAS:
        _write(root / name, f">{name}\nACGTACGT\n".encode())
        _write(root / f"{name}.fai", f"{name}\t8\t{len(name) + 2}\t8\t9\n".encode())
    _write(root / "vntr_db_advntr" / "hg19_muc1.db", b"hg19 advntr database")
    _write(root / "vntr_db_advntr" / "hg38_muc1.db", b"hg38 advntr database")
    _write(root / "vntr_db_advntr.zip", b"PK the advntr database archive")
    _write(root / "filter_config.json", b'{"1": ["2"]}')
    _write(root / "install_references.log", b"2026-08-11 [INFO] Logging initialized.\n")
    return root


def _spec_document(refs_root: Path) -> dict[str, Any]:
    return {
        "release_tag": TAG,
        "source_commit": "b" * 40,
        "bwa_version": "0.7.18-r1243",
        "samtools_version": "1.20",
        "sources": {
            ref_id: {
                "url": f"https://example.invalid/{ref_id}/chr1.gz",
                "source_sha256": _sha256(refs_root / f"{fasta}.gz"),
            }
            for ref_id, fasta in GENOMES.items()
        },
        "seeds": {name: {"sha256": _sha256(refs_root / name)} for name in SEEDS},
        "derivations": [
            {
                "output": name,
                "command": f"samtools faidx <source> > {name}",
                "expected_sha256": _sha256(refs_root / name),
            }
            for name in DERIVED
        ],
    }


@pytest.fixture
def spec(tmp_path: Path, refs: Path) -> Path:
    path = tmp_path / "refs-v1.json"
    path.write_text(json.dumps(_spec_document(refs), indent=2), encoding="utf-8")
    return path


@pytest.fixture(autouse=True)
def _pinned_toolchain(monkeypatch: pytest.MonkeyPatch) -> None:
    """Report the versions the spec pins, so no test depends on bwa or samtools existing.

    Unit tests must not shell out to the aligner stack; `_capture` is the single seam
    where `bundle_release` does.
    """
    outputs = {
        "bwa": "Program: bwa (alignment via Burrows-Wheeler transform)\nVersion: 0.7.18-r1243\n",
        "samtools": "samtools 1.20\nUsing htslib 1.20\n",
    }
    monkeypatch.setattr(bundle_release, "_capture", lambda argv: outputs.get(argv[0], ""))


def _argv(refs_root: Path, spec_path: Path, out: Path, *extra: str) -> list[str]:
    return [
        "--refs",
        str(refs_root),
        "--spec",
        str(spec_path),
        "--tag",
        TAG,
        "--out",
        str(out),
        "--data-sha",
        "a" * 40,
        "--builder-sha",
        "b" * 40,
        *extra,
    ]


@pytest.fixture
def built(tmp_path: Path, refs: Path, spec: Path) -> Path:
    """Run a complete build and hand back the populated `dist/` directory."""
    out = tmp_path / "dist"
    assert bundle_release.main(_argv(refs, spec, out)) == 0
    return out


class TestAssetGrouping:
    def test_the_seven_asset_names_are_the_frozen_map(self, refs: Path) -> None:
        """The genome asset names are source-prefixed, not the physical id."""
        assets = bundle_release.plan_assets(refs, TAG, list(GENOMES))
        assert {asset.name for asset in assets} == EXPECTED_ASSETS

    def test_the_tag_is_taken_from_the_argument_not_hardcoded(self, refs: Path) -> None:
        assets = bundle_release.plan_assets(refs, "refs-v2", list(GENOMES))
        assert "vntyper-references-refs-v2-muc1.tar.gz" in {asset.name for asset in assets}

    def test_a_genome_asset_carries_its_fasta_and_all_six_sidecars(self, refs: Path) -> None:
        (asset,) = [a for a in bundle_release.plan_assets(refs, TAG, ["hg38"]) if a.reference_id == "hg38"]
        assert set(asset.members) == {"alignment/chr1.hg38.fa"} | {
            f"alignment/chr1.hg38.fa{suffix}" for suffix in SIDECARS
        }

    def test_the_ncbi_assemblies_keep_their_fna_spelling(self, refs: Path) -> None:
        (asset,) = [a for a in bundle_release.plan_assets(refs, TAG, ["GRCh37"]) if a.reference_id == "GRCh37"]
        assert "alignment/chr1.GRCh37.fna" in asset.members

    def test_the_muc1_asset_carries_every_muc1_fasta_its_index_and_both_databases(self, refs: Path) -> None:
        (asset,) = [a for a in bundle_release.plan_assets(refs, TAG, []) if a.reference_id is None]
        expected = {name for fasta in MUC1_FASTAS for name in (fasta, f"{fasta}.fai")}
        expected |= {"vntr_db_advntr/hg19_muc1.db", "vntr_db_advntr/hg38_muc1.db"}
        assert set(asset.members) == expected


class TestEveryFileIsAccountedFor:
    def test_the_assets_and_the_exclusions_together_cover_the_whole_tree(self, refs: Path) -> None:
        assets = bundle_release.plan_assets(refs, TAG, list(GENOMES))
        excluded = bundle_release.assign_every_file(_tree_files(refs), assets)
        assigned = {member for asset in assets for member in asset.members}
        assert assigned | set(excluded) == _tree_files(refs)

    def test_no_file_is_carried_by_two_assets(self, refs: Path) -> None:
        assets = bundle_release.plan_assets(refs, TAG, list(GENOMES))
        members = [member for asset in assets for member in asset.members]
        assert len(members) == len(set(members))

    def test_an_unexpected_file_fails_the_build_by_name(self, refs: Path) -> None:
        _write(refs / "alignment" / "chr1.hg38.fa.bwt.old", b"leftover")
        assets = bundle_release.plan_assets(refs, TAG, list(GENOMES))
        with pytest.raises(ValueError, match=r"chr1\.hg38\.fa\.bwt\.old"):
            bundle_release.assign_every_file(_tree_files(refs), assets)

    def test_an_unexpected_file_stops_the_whole_build(self, tmp_path: Path, refs: Path, spec: Path) -> None:
        _write(refs / "surprise.txt", b"not in any asset")
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 1
        assert not list(out.glob("*.tar.gz"))

    @pytest.mark.parametrize("missing", ["alignment/chr1.hg38.fa.fai", "alignment/chr1.hg38.fa.sa"])
    def test_a_genome_missing_an_index_file_fails_by_name(self, refs: Path, missing: str) -> None:
        (refs / missing).unlink()
        with pytest.raises(ValueError, match=Path(missing).name.replace(".", r"\.")):
            bundle_release.plan_assets(refs, TAG, list(GENOMES))

    def test_a_genome_with_no_extracted_fasta_fails_by_reference_id(self, refs: Path) -> None:
        for path in (refs / "alignment").glob("chr1.hg19.fa*"):
            if path.name != "chr1.hg19.fa.gz":
                path.unlink()
        with pytest.raises(ValueError, match="hg19"):
            bundle_release.plan_assets(refs, TAG, list(GENOMES))

    def test_a_missing_muc1_index_fails_by_name(self, refs: Path) -> None:
        (refs / "code-adVNTR_RUs.fa.fai").unlink()
        with pytest.raises(ValueError, match=r"code-adVNTR_RUs\.fa\.fai"):
            bundle_release.plan_assets(refs, TAG, list(GENOMES))

    @pytest.mark.parametrize(
        "relative",
        [
            "install_references.log",
            "md5_checksums.txt",
            "filter_config.json",
            "vntr_db_advntr.zip",
            "alignment/chr1.hg38.fa.gz",
        ],
    )
    def test_the_exclusion_list_is_explicit_and_reasoned(self, relative: str) -> None:
        reason = bundle_release.exclusion_reason(relative)
        assert reason and len(reason) > 10

    def test_nothing_else_is_excluded_by_accident(self) -> None:
        assert bundle_release.exclusion_reason("alignment/chr1.hg38.fa") is None
        assert bundle_release.exclusion_reason("muc1_region_hg19.fa") is None

    def test_a_symlink_in_the_tree_stops_the_build(self, refs: Path) -> None:
        """An archived link records a digest that describes something else."""
        (refs / "alignment" / "chr1.link.fa").symlink_to(refs / "alignment" / "chr1.hg38.fa")
        with pytest.raises(ValueError, match=r"chr1\.link\.fa"):
            bundle_release.inventory(refs)

    def test_an_unknown_reference_id_is_rejected(self, refs: Path) -> None:
        with pytest.raises(ValueError, match="chm13"):
            bundle_release.plan_assets(refs, TAG, ["chm13"])

    def test_two_assets_cannot_claim_the_same_file(self, refs: Path) -> None:
        assets = bundle_release.plan_assets(refs, TAG, ["hg38", "hg38"])
        with pytest.raises(ValueError, match="claimed by both"):
            bundle_release.assign_every_file(_tree_files(refs), assets)

    def test_an_asset_member_absent_from_the_tree_is_named(self, refs: Path) -> None:
        assets = bundle_release.plan_assets(refs, TAG, ["hg38"])
        (refs / "alignment" / "chr1.hg38.fa.pac").unlink()
        with pytest.raises(ValueError, match=r"chr1\.hg38\.fa\.pac"):
            bundle_release.assign_every_file(_tree_files(refs), assets)


class TestArchiveLayout:
    def test_every_asset_extracts_straight_back_into_the_reference_tree_layout(
        self, tmp_path: Path, refs: Path, built: Path
    ) -> None:
        """Member paths are tree-relative, so extraction reproduces the layout directly.

        `safe_extract` rejects absolute paths and `..` components, so this also proves
        the installer will accept every archive.
        """
        restored = tmp_path / "restored"
        for archive in sorted(built.glob("*.tar.gz")):
            safe_extract(archive, restored)

        rebuilt = _tree_files(restored) - {"release-manifest.json", "BUILD_INFO.json"}
        expected = _tree_files(refs) - {path for path in _tree_files(refs) if bundle_release.exclusion_reason(path)}
        assert rebuilt == expected
        for relative in sorted(rebuilt):
            assert (restored / relative).read_bytes() == (refs / relative).read_bytes()

    def test_no_member_path_is_prefixed_with_the_tree_root(self, built: Path) -> None:
        for archive in sorted(built.glob("*.tar.gz")):
            with tarfile.open(archive, "r:gz") as tar:
                names = tar.getnames()
            assert not [name for name in names if name.startswith(("/", "refs/", "./"))], archive.name

    def test_every_tarball_carries_its_own_manifest_and_build_info(self, built: Path) -> None:
        """PR-2 verifies extracted files against a manifest under the archive's digest.

        A loose metadata file beside the assets has no committed digest to check it
        against, so an actor who can replace an asset can replace the metadata too.
        """
        archives = sorted(built.glob("*.tar.gz"))
        assert len(archives) == len(EXPECTED_ASSETS)
        for archive in archives:
            with tarfile.open(archive, "r:gz") as tar:
                names = set(tar.getnames())
            assert {"release-manifest.json", "BUILD_INFO.json"} <= names, archive.name

    def test_the_in_archive_manifest_describes_that_archive_alone(self, built: Path) -> None:
        archive = built / "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        with tarfile.open(archive, "r:gz") as tar:
            member = tar.extractfile("release-manifest.json")
            assert member is not None
            manifest = json.load(member)
            names = set(tar.getnames())
        assert manifest["asset"] == archive.name
        assert manifest["reference_id"] == "hg19"
        assert {entry["path"] for entry in manifest["files"]} == names - set(manifest["metadata"])

    def test_the_in_archive_manifest_digests_describe_the_extracted_bytes(self, tmp_path: Path, built: Path) -> None:
        restored = tmp_path / "verify"
        archive = built / "vntyper-references-refs-v1-muc1.tar.gz"
        safe_extract(archive, restored)
        manifest = json.loads((restored / "release-manifest.json").read_text(encoding="utf-8"))
        assert manifest["files"]
        for entry in manifest["files"]:
            extracted = restored / entry["path"]
            assert extracted.exists(), entry["path"]
            assert _sha256(extracted) == entry["sha256"], entry["path"]
            assert extracted.stat().st_size == entry["size"], entry["path"]

    def test_the_in_archive_build_info_names_the_versions_the_installer_compares(self, built: Path) -> None:
        with tarfile.open(built / "vntyper-references-refs-v1-muc1.tar.gz", "r:gz") as tar:
            member = tar.extractfile("BUILD_INFO.json")
            assert member is not None
            info = json.load(member)
        assert info["bwa_version"] == "0.7.18-r1243"
        assert info["samtools_version"] == "1.20"
        assert info["builder_commit"] == "b" * 40
        assert info["data_commit"] == "a" * 40
        assert info["release_tag"] == TAG

    def test_the_source_gz_never_reaches_an_asset(self, built: Path) -> None:
        for archive in sorted(built.glob("*.tar.gz")):
            with tarfile.open(archive, "r:gz") as tar:
                assert not [name for name in tar.getnames() if name.endswith(".gz")], archive.name

    def test_the_installer_log_and_the_build_input_never_reach_an_asset(self, built: Path) -> None:
        shipped: set[str] = set()
        for archive in sorted(built.glob("*.tar.gz")):
            with tarfile.open(archive, "r:gz") as tar:
                shipped |= set(tar.getnames())
        assert "install_references.log" not in shipped
        assert "filter_config.json" not in shipped
        assert "vntr_db_advntr.zip" not in shipped


class TestReleaseLevelOutputs:
    def test_the_release_publishes_the_reviewer_facing_files(self, built: Path) -> None:
        for name in ("SHA256SUMS", "release-manifest.json", "BUILD_INFO.json", "verification-report.json"):
            assert (built / name).exists(), name

    def test_sha256sums_describes_every_asset_and_nothing_else(self, built: Path) -> None:
        lines = (built / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        recorded = {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in lines}
        assert set(recorded) == EXPECTED_ASSETS
        for name, digest in recorded.items():
            assert _sha256(built / name) == digest

    def test_the_top_level_manifest_covers_every_asset_and_every_exclusion(self, refs: Path, built: Path) -> None:
        manifest = json.loads((built / "release-manifest.json").read_text(encoding="utf-8"))
        assert {asset["name"] for asset in manifest["assets"]} == EXPECTED_ASSETS
        assigned = {entry["path"] for asset in manifest["assets"] for entry in asset["files"]}
        excluded = {entry["path"] for entry in manifest["excluded"]}
        assert assigned | excluded == _tree_files(refs)
        assert not assigned & excluded
        assert all(entry["reason"] for entry in manifest["excluded"])

    def test_the_top_level_manifest_records_each_asset_digest_and_size(self, built: Path) -> None:
        manifest = json.loads((built / "release-manifest.json").read_text(encoding="utf-8"))
        for asset in manifest["assets"]:
            archive = built / asset["name"]
            assert asset["sha256"] == _sha256(archive)
            assert asset["size"] == archive.stat().st_size

    def test_provenance_reaches_the_manifest(self, built: Path) -> None:
        manifest = json.loads((built / "release-manifest.json").read_text(encoding="utf-8"))
        entries = {entry["path"]: entry for asset in manifest["assets"] for entry in asset["files"]}
        assert entries["alignment/chr1.hg38.fa"]["source_url"] == "https://example.invalid/hg38/chr1.gz"
        assert "bwa index" in entries["alignment/chr1.hg38.fa.bwt"]["produced_by"]
        assert "samtools faidx" in entries["muc1_region_hg19.fa"]["produced_by"]
        assert "seed" in entries["MUC1_motifs_Rev_com.fa"]["produced_by"]
        assert "vntr_db_advntr.zip" in entries["vntr_db_advntr/hg19_muc1.db"]["produced_by"]

    def test_the_verification_report_records_every_digest_it_checked(self, built: Path) -> None:
        report = json.loads((built / "verification-report.json").read_text(encoding="utf-8"))
        assert report["ok"] is True
        checked = {(check["check"], check["name"]) for check in report["checks"]}
        for name in DERIVED:
            assert ("derivation", name) in checked
        for name in SEEDS:
            assert ("seed", name) in checked
        assert ("source-archive", "alignment/chr1.hg19.fa.gz") in checked
        assert all("verdict" in check for check in report["checks"])


class TestFailClosed:
    def _rewrite_spec(self, spec_path: Path, mutate: Any) -> None:
        document = json.loads(spec_path.read_text(encoding="utf-8"))
        mutate(document)
        spec_path.write_text(json.dumps(document), encoding="utf-8")

    def test_a_derivation_that_does_not_match_the_spec_digest_stops_the_build(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["derivations"][0]["expected_sha256"] = "0" * 64

        self._rewrite_spec(spec, mutate)
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 1
        assert not list(out.glob("*.tar.gz"))

    def test_a_failed_check_still_leaves_a_verification_report_behind(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["seeds"]["code-adVNTR_RUs.fa"]["sha256"] = "0" * 64

        self._rewrite_spec(spec, mutate)
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 1
        report = json.loads((out / "verification-report.json").read_text(encoding="utf-8"))
        assert report["ok"] is False
        failed = [check for check in report["checks"] if not check["ok"]]
        assert [check["name"] for check in failed] == ["code-adVNTR_RUs.fa"]

    def test_a_source_archive_that_does_not_match_the_spec_stops_the_build(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["sources"]["GRCh38"]["source_sha256"] = "0" * 64

        self._rewrite_spec(spec, mutate)
        assert bundle_release.main(_argv(refs, spec, tmp_path / "dist")) == 1

    def test_a_spec_that_does_not_name_a_source_fails_by_reference_id(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        def mutate(document: dict[str, Any]) -> None:
            del document["sources"]["hg19_ensembl"]

        self._rewrite_spec(spec, mutate)
        assert bundle_release.main(_argv(refs, spec, tmp_path / "dist")) == 1

    def test_a_spec_written_for_another_tag_stops_the_build(self, tmp_path: Path, refs: Path, spec: Path) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["release_tag"] = "refs-v9"

        self._rewrite_spec(spec, mutate)
        assert bundle_release.main(_argv(refs, spec, tmp_path / "dist")) == 1

    def test_a_toolchain_that_differs_from_the_spec_is_reported_but_not_fatal(
        self, tmp_path: Path, refs: Path, spec: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The digest assertions all passed; the bundle is a draft a human reviews."""
        monkeypatch.setattr(bundle_release, "_capture", lambda argv: "Version: 0.7.17-r1188\nsamtools 1.19\n")
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 0
        report = json.loads((out / "verification-report.json").read_text(encoding="utf-8"))
        toolchain = [check for check in report["checks"] if check["check"] == "toolchain"]
        assert {check["name"] for check in toolchain} == {"bwa_version", "samtools_version"}
        assert all(not check["ok"] and not check["fatal"] for check in toolchain)
        assert report["ok"] is True

    def test_a_seed_the_spec_names_but_the_build_did_not_stage_is_noted_not_failed(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        """`verify_seeds.py` already checked it in the seeds directory."""

        def mutate(document: dict[str, Any]) -> None:
            document["seeds"]["never_staged.fa"] = {"sha256": "0" * 64}

        self._rewrite_spec(spec, mutate)
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 0
        report = json.loads((out / "verification-report.json").read_text(encoding="utf-8"))
        (record,) = [check for check in report["checks"] if check["name"] == "never_staged.fa"]
        assert record["ok"] and not record["fatal"]

    def test_a_derivation_the_spec_names_but_the_build_did_not_produce_stops_it(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["derivations"].append({"output": "muc1_region_chm13.fa", "expected_sha256": "0" * 64})

        self._rewrite_spec(spec, mutate)
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 1
        report = json.loads((out / "verification-report.json").read_text(encoding="utf-8"))
        (record,) = [check for check in report["checks"] if check["name"] == "muc1_region_chm13.fa"]
        assert not record["ok"] and record["fatal"]

    def test_a_source_with_no_digest_and_a_discarded_download_are_both_only_noted(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        """The installer verified the download before decompressing it either way."""

        def mutate(document: dict[str, Any]) -> None:
            document["sources"]["hg19"].pop("source_sha256")

        self._rewrite_spec(spec, mutate)
        (refs / "alignment" / "chr1.hg38.fa.gz").unlink()
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 0
        report = json.loads((out / "verification-report.json").read_text(encoding="utf-8"))
        noted = {
            check["name"]: check
            for check in report["checks"]
            if check["check"] == "source-archive" and not check["fatal"]
        }
        assert set(noted) == {"alignment/chr1.hg19.fa.gz", "alignment/chr1.hg38.fa.gz"}
        assert all(check["ok"] for check in noted.values())

    def test_a_spec_that_pins_no_toolchain_versions_records_no_toolchain_check(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document.pop("bwa_version")
            document.pop("samtools_version")

        self._rewrite_spec(spec, mutate)
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 0
        report = json.loads((out / "verification-report.json").read_text(encoding="utf-8"))
        assert not [check for check in report["checks"] if check["check"] == "toolchain"]

    def test_a_spec_with_no_derivations_still_builds_and_says_so_in_the_manifest(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        """The derivations were already asserted against the committed install config."""

        def mutate(document: dict[str, Any]) -> None:
            document["derivations"] = []

        self._rewrite_spec(spec, mutate)
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 0
        manifest = json.loads((out / "release-manifest.json").read_text(encoding="utf-8"))
        entries = {entry["path"]: entry for asset in manifest["assets"] for entry in asset["files"]}
        assert entries["muc1_region_hg19.fa"]["produced_by"] == "install-references --from-source"
        report = json.loads((out / "verification-report.json").read_text(encoding="utf-8"))
        assert not [check for check in report["checks"] if check["check"] == "derivation"]

    def test_a_spec_with_no_seeds_block_is_rejected_before_anything_is_built(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        def mutate(document: dict[str, Any]) -> None:
            document["seeds"] = {}

        self._rewrite_spec(spec, mutate)
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out)) == 1
        assert not list(out.glob("*.tar.gz")) if out.exists() else True


class TestSpecPreflight:
    def test_check_spec_only_passes_a_well_formed_spec_without_reading_the_tree(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        """The workflow runs this before three hours of downloading and indexing."""
        empty = tmp_path / "not-built-yet"
        empty.mkdir()
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(empty, spec, out, "--check-spec-only")) == 0
        assert not out.exists()

    def test_check_spec_only_rejects_a_spec_missing_a_source(self, tmp_path: Path, refs: Path, spec: Path) -> None:
        document = json.loads(spec.read_text(encoding="utf-8"))
        del document["sources"]["GRCh37"]
        spec.write_text(json.dumps(document), encoding="utf-8")
        empty = tmp_path / "not-built-yet"
        empty.mkdir()
        assert bundle_release.main(_argv(empty, spec, tmp_path / "dist", "--check-spec-only")) == 1

    def test_check_spec_only_rejects_a_source_with_no_url(self, tmp_path: Path, refs: Path, spec: Path) -> None:
        document = json.loads(spec.read_text(encoding="utf-8"))
        document["sources"]["GRCh37"].pop("url")
        spec.write_text(json.dumps(document), encoding="utf-8")
        empty = tmp_path / "not-built-yet"
        empty.mkdir()
        assert bundle_release.main(_argv(empty, spec, tmp_path / "dist", "--check-spec-only")) == 1

    def test_a_spec_that_is_not_json_fails_with_its_path(self, tmp_path: Path) -> None:
        broken = tmp_path / "broken.json"
        broken.write_text("{ not json", encoding="utf-8")
        with pytest.raises(ValueError, match="broken.json"):
            bundle_release.load_spec(broken)

    def test_a_spec_that_cannot_be_read_fails_with_its_path(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="absent.json"):
            bundle_release.load_spec(tmp_path / "absent.json")

    def test_a_spec_that_is_not_an_object_is_rejected(self, tmp_path: Path) -> None:
        listy = tmp_path / "listy.json"
        listy.write_text("[]", encoding="utf-8")
        with pytest.raises(ValueError, match="must be a JSON object"):
            bundle_release.load_spec(listy)

    def test_a_seed_entry_with_no_digest_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="MUC1_motifs_Rev_com.fa"):
            bundle_release.spec_seed_digests({"seeds": {"MUC1_motifs_Rev_com.fa": {}}})

    def test_a_spec_with_no_sources_block_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="no 'sources' block"):
            bundle_release.spec_source({}, "hg19")

    @pytest.mark.parametrize(
        "derivations",
        ["not-a-list", [{"output": "muc1_region_hg19.fa"}], [{"expected_sha256": "0" * 64}]],
    )
    def test_a_malformed_derivations_block_is_rejected(self, refs: Path, derivations: Any) -> None:
        document = _spec_document(refs)
        document["derivations"] = derivations
        with pytest.raises(ValueError, match="derivation"):
            bundle_release.validate_spec(document, TAG, [])


class TestPruning:
    def test_pruning_removes_each_asset_only_after_its_archive_is_written(
        self, tmp_path: Path, refs: Path, spec: Path
    ) -> None:
        """Peak disk is the reason: six indexed genomes plus their tarballs is ~7.5 GB."""
        out = tmp_path / "dist"
        assert bundle_release.main(_argv(refs, spec, out, "--prune")) == 0
        assert len(list(out.glob("*.tar.gz"))) == len(EXPECTED_ASSETS)
        assert not (refs / "alignment" / "chr1.hg38.fa").exists()
        assert not (refs / "alignment" / "chr1.hg38.fa.gz").exists()
        assert not (refs / "muc1_region_hg19.fa").exists()
        # The archives are still complete and still verify against SHA256SUMS.
        lines = (out / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        for line in lines:
            digest, name = line.split("  ", 1)
            assert _sha256(out / name) == digest

    def test_without_pruning_the_tree_survives_the_build(self, refs: Path, built: Path) -> None:
        assert (refs / "alignment" / "chr1.hg38.fa").exists()
        assert (refs / "alignment" / "chr1.hg38.fa.gz").exists()


class TestToolchainCapture:
    def test_a_missing_binary_reports_no_version_rather_than_crashing(self) -> None:
        assert _REAL_CAPTURE(["definitely-not-a-real-binary-9f2a"]) == ""

    def test_a_real_probe_returns_both_streams(self) -> None:
        """bwa prints its banner to stderr and exits 1, so neither may be dropped."""
        assert "Python" in _REAL_CAPTURE([sys.executable, "--version"])

    def test_versions_are_parsed_out_of_the_usual_banners(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(
            bundle_release,
            "_capture",
            lambda argv: {
                "bwa": "Program: bwa\nVersion: 0.7.18-r1243\nContact: ...\n",
                "samtools": "samtools 1.20\nUsing htslib 1.20\n",
            }[argv[0]],
        )
        assert bundle_release.observed_toolchain() == {
            "bwa_version": "0.7.18-r1243",
            "samtools_version": "1.20",
        }

    def test_an_unrecognisable_banner_yields_none(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(bundle_release, "_capture", lambda argv: "")
        assert bundle_release.observed_toolchain() == {"bwa_version": None, "samtools_version": None}


class TestWorkflowAgreement:
    """The workflow is the only caller, and two of its arguments are load-bearing."""

    WORKFLOW = REPO_ROOT / ".github" / "workflows" / "build-reference-bundles.yml"

    def _workflow(self) -> dict[str, Any]:
        yaml = pytest.importorskip("yaml")
        if not self.WORKFLOW.exists():
            pytest.skip("no GitHub Actions workflows present in this tree")
        return dict(yaml.safe_load(self.WORKFLOW.read_text(encoding="utf-8")))

    def _steps(self) -> list[dict[str, Any]]:
        return list(self._workflow()["jobs"]["build"]["steps"])

    def _run_text(self) -> str:
        return "\n".join(str(step.get("run", "")) for step in self._steps())

    def test_the_staging_step_copies_the_filter_config(self) -> None:
        """The literal derivation reads it; without it the build dies after every
        genome has been downloaded and indexed."""
        assert "filter_config.json" in self._run_text()

    def test_the_provenance_only_seed_is_not_staged(self) -> None:
        """`generate_vntr_reference.py` is superseded by `merge_pairwise_motifs`.
        Copying it would leave an unassigned file for the grouper to trip over."""
        assert "generate_vntr_reference.py" not in self._run_text()

    def test_all_six_genomes_are_selected(self) -> None:
        """A short list makes an out-of-scope derivation skip, silently producing a
        bundle with no `muc1_region_hg38.fa` in it."""
        run = self._run_text()
        assert "--references hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl" in run

    def test_the_spec_is_checked_before_the_genomes_are_downloaded(self) -> None:
        run_texts = [str(step.get("run", "")) for step in self._steps()]
        preflight = [index for index, text in enumerate(run_texts) if "--check-spec-only" in text]
        install = [index for index, text in enumerate(run_texts) if "install-references" in text]
        assert preflight and install
        assert min(preflight) < min(install)

    def test_the_release_is_published_as_a_draft(self) -> None:
        assert "--draft" in self._run_text()

    def test_untrusted_inputs_are_passed_through_the_environment(self) -> None:
        """`${{ inputs.* }}` interpolated straight into a `run:` block is a script
        injection; every caller-supplied value is bound to an env var instead."""
        for step in self._steps():
            run = str(step.get("run", ""))
            assert "${{ inputs." not in run, step.get("name")


class TestDigestAgreement:
    """`bundle_release.sha256_of` is a copy, and a copy that drifts is a wrong digest."""

    def test_the_copy_agrees_with_reference_bundle_on_every_shape_of_input(self, tmp_path: Path) -> None:
        for index, payload in enumerate([b"", b"\x00", b">chr1\nACGT\n", b"x" * (1024 * 1024 + 7)]):
            sample = tmp_path / f"sample{index}"
            sample.write_bytes(payload)
            assert bundle_release.sha256_of(sample) == reference_bundle.sha256_of(sample)
            assert bundle_release.sha256_of(sample) == _sha256(sample)

    def test_the_script_imports_nothing_from_vntyper(self) -> None:
        """It is run by path from a second checkout, where `import vntyper` may not resolve."""
        source = (REPO_ROOT / "scripts" / "bundle_release.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        imported = {
            node.module
            for node in tree.body
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("vntyper")
        }
        assert not imported


class TestVerifySeeds:
    """The seeds are checked against the spec before three hours of work begins."""

    def _seed_dir(self, tmp_path: Path) -> Path:
        seeds = tmp_path / "seeds"
        seeds.mkdir()
        for name in SEEDS:
            _write(seeds / name, f"contents of {name}".encode())
        # The data repository also keeps files the spec does not name.
        _write(seeds / "generate_vntr_reference.py", b"# provenance only\n")
        return seeds

    def _spec(self, tmp_path: Path, seeds_dir: Path, **overrides: Any) -> Path:
        document: dict[str, Any] = {
            "release_tag": TAG,
            "seeds": {name: {"sha256": _sha256(seeds_dir / name)} for name in SEEDS},
        }
        document.update(overrides)
        path = tmp_path / "spec.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def test_a_matching_seed_set_exits_zero(self, tmp_path: Path) -> None:
        seeds = self._seed_dir(tmp_path)
        spec = self._spec(tmp_path, seeds)
        assert verify_seeds.main(["--spec", str(spec), "--seeds", str(seeds)]) == 0

    def test_a_seed_that_does_not_match_exits_non_zero(self, tmp_path: Path) -> None:
        seeds = self._seed_dir(tmp_path)
        spec = self._spec(tmp_path, seeds)
        (seeds / "filter_config.json").write_bytes(b"tampered")
        assert verify_seeds.main(["--spec", str(spec), "--seeds", str(seeds)]) == 1

    def test_a_missing_seed_is_named(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        seeds = self._seed_dir(tmp_path)
        spec = self._spec(tmp_path, seeds)
        (seeds / "vntr_db_advntr.zip").unlink()
        assert verify_seeds.main(["--spec", str(spec), "--seeds", str(seeds)]) == 1
        assert "vntr_db_advntr.zip" in caplog.text

    def test_an_empty_seeds_block_is_rejected(self, tmp_path: Path) -> None:
        seeds = self._seed_dir(tmp_path)
        spec = self._spec(tmp_path, seeds, seeds={})
        assert verify_seeds.main(["--spec", str(spec), "--seeds", str(seeds)]) == 1

    def test_a_seed_entry_with_no_digest_is_rejected(self, tmp_path: Path) -> None:
        seeds = self._seed_dir(tmp_path)
        spec = self._spec(tmp_path, seeds)
        document = json.loads(spec.read_text(encoding="utf-8"))
        document["seeds"]["MUC1_motifs_Rev_com.fa"] = {}
        spec.write_text(json.dumps(document), encoding="utf-8")
        assert verify_seeds.main(["--spec", str(spec), "--seeds", str(seeds)]) == 1

    def test_a_bare_string_digest_is_accepted(self, tmp_path: Path) -> None:
        seeds = self._seed_dir(tmp_path)
        spec = self._spec(tmp_path, seeds)
        document = json.loads(spec.read_text(encoding="utf-8"))
        document["seeds"] = {name: entry["sha256"] for name, entry in document["seeds"].items()}
        spec.write_text(json.dumps(document), encoding="utf-8")
        assert verify_seeds.main(["--spec", str(spec), "--seeds", str(seeds)]) == 0

    def test_a_spec_that_cannot_be_read_exits_non_zero(self, tmp_path: Path) -> None:
        seeds = self._seed_dir(tmp_path)
        assert verify_seeds.main(["--spec", str(tmp_path / "absent.json"), "--seeds", str(seeds)]) == 1

    @pytest.mark.parametrize("payload", ["{ not json", "[]"])
    def test_a_spec_that_is_not_a_json_object_exits_non_zero(self, tmp_path: Path, payload: str) -> None:
        seeds = self._seed_dir(tmp_path)
        broken = tmp_path / "broken.json"
        broken.write_text(payload, encoding="utf-8")
        assert verify_seeds.main(["--spec", str(broken), "--seeds", str(seeds)]) == 1

    def test_a_seed_name_that_escapes_the_seed_directory_is_rejected(self, tmp_path: Path) -> None:
        """A spec is committed data from another repository; treat its keys as paths."""
        seeds = self._seed_dir(tmp_path)
        spec = self._spec(tmp_path, seeds, seeds={"../outside.fa": {"sha256": "0" * 64}})
        assert verify_seeds.main(["--spec", str(spec), "--seeds", str(seeds)]) == 1
