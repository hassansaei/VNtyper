"""`install-references`'s default path: fetch the published, checksummed release bundle.

Milestone 5 replaces the six-third-party-host download with a fetch of one immutable
GitHub release (`berntpopp/vntyper-data`, tag `refs-v1`). `install_from_bundle` is the
function that does it: for every selected genome plus the always-installed common
asset, it downloads `https://github.com/{repository}/releases/download/{release_tag}/
{asset}`, verifies the bytes against the SHA-256 committed in
`install_references_config.json` *before* extracting a single byte (Task 3 review
carry-forward), extracts through `safe_extract`, and activates through
`staged_install` - so a failure partway through a multi-asset install leaves any
previously installed tree untouched rather than half-replaced.

Every asset also carries `release-manifest.json` and `BUILD_INFO.json` *inside* the
archive (never fetched separately - a loose metadata file has no committed digest to
check it against, and trusting one would reopen exactly the hole `asset_sha256`
closes). This file's fixtures build real `.tar.gz` archives with both members inside,
the way `scripts/bundle_release.py` does, rather than asserting against mocks: a
manifest that agrees with a mock but not with the archive is exactly the failure this
installer exists to catch.

Network access is faked throughout via `download_file`, which every scenario
monkeypatches to copy from a `tmp_path`-local "release" directory instead of reaching
the network - see AGENTS.md's unit-tier rule (`tmp_path` + `unittest.mock`, no network).
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import subprocess
import tarfile
from pathlib import Path

import pytest

from vntyper.scripts import install_references

pytestmark = pytest.mark.unit


# =============================================================================
# Fixture builders: real tar.gz archives, real digests, no mocks standing in for bytes.
# =============================================================================


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _write_bundle_asset(
    release_dir: Path,
    name: str,
    files: dict[str, bytes],
    *,
    reference_id: str | None,
    bwa_version: str | None = "9.9.9",
) -> str:
    """Build one `.tar.gz` the way `scripts/bundle_release.py` does: the named files
    plus `release-manifest.json` and `BUILD_INFO.json` written *inside* the archive.

    Args:
        release_dir: Directory standing in for the published release; the fake
            `download_file` below serves assets from here.
        name: Asset file name, e.g. ``vntyper-references-refs-v1-ucsc-hg19.tar.gz``.
        files: Tree-relative path -> payload, for every non-metadata member.
        reference_id: `release-manifest.json`'s `reference_id` field.
        bwa_version: `BUILD_INFO.json`'s `bwa_version` field. None omits the key
            entirely, matching a build where the toolchain probe found nothing.

    Returns:
        str: The archive's own SHA-256 - the value a fixture config's `asset_sha256`
        (or `common_asset_sha256`) must carry for the happy path to verify.
    """
    manifest_files = [
        {"path": path, "size": len(payload), "sha256": _sha256_bytes(payload)} for path, payload in files.items()
    ]
    manifest = json.dumps({"asset": name, "reference_id": reference_id, "files": manifest_files}).encode("utf-8")
    build_info_doc: dict[str, object] = {"samtools_version": "1.20"}
    if bwa_version is not None:
        build_info_doc["bwa_version"] = bwa_version
    build_info = json.dumps(build_info_doc).encode("utf-8")

    release_dir.mkdir(parents=True, exist_ok=True)
    archive = release_dir / name
    with tarfile.open(archive, "w:gz") as tar:
        for path, payload in {**files, "release-manifest.json": manifest, "BUILD_INFO.json": build_info}.items():
            info = tarfile.TarInfo(path)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))

    return _sha256_file(archive)


def _write_raw_asset(release_dir: Path, name: str, members: dict[str, bytes]) -> str:
    """Build a `.tar.gz` from exactly the given members - no metadata auto-injected.

    Used to construct malformed or incomplete assets: one missing `release-manifest.json`
    entirely, one whose metadata is not valid JSON, or one whose manifest names a file
    that is not actually a member.

    Args:
        release_dir: Directory standing in for the published release.
        name: Asset file name.
        members: Every member to write, metadata included where the test wants it.

    Returns:
        str: The archive's own SHA-256.
    """
    release_dir.mkdir(parents=True, exist_ok=True)
    archive = release_dir / name
    with tarfile.open(archive, "w:gz") as tar:
        for path, payload in members.items():
            info = tarfile.TarInfo(path)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return _sha256_file(archive)


def _genome_files(ref_id: str, fasta: bytes = b">chr1\nACGTACGTAC\n") -> dict[str, bytes]:
    """The seven members every genome asset carries: FASTA, `.fai`, five BWA sidecars."""
    prefix = f"alignment/chr1.{ref_id}.fa"
    return {
        prefix: fasta,
        f"{prefix}.fai": b"chr1\t10\t6\t10\t11\n",
        f"{prefix}.amb": b"amb-bytes",
        f"{prefix}.ann": b"ann-bytes",
        f"{prefix}.bwt": b"bwt-bytes",
        f"{prefix}.pac": b"pac-bytes",
        f"{prefix}.sa": b"sa-bytes",
    }


def _common_files() -> dict[str, bytes]:
    """A lean stand-in for the muc1 asset's twelve members: enough to prove the
    common asset installs regardless of `--references`, without the derivation outputs
    this module's own logic never reads."""
    return {
        "MUC1_motifs_Rev_com.fa": b">motif\nACGT\n",
        "code-adVNTR_RUs.fa": b">ru\nACGT\n",
        "vntr_db_advntr/hg19_muc1.db": b"hg19-db-bytes",
        "vntr_db_advntr/hg38_muc1.db": b"hg38-db-bytes",
    }


def _genome_entry(ref_id: str, asset: str, asset_sha256: str, *, kind: str = "bwa") -> dict:
    return {
        "kind": kind,
        "installed_path": f"alignment/chr1.{ref_id}.fa",
        "index_command": "bwa index {path}",
        "asset": asset,
        "asset_sha256": asset_sha256,
    }


def _install_config(*, ucsc: dict, bundle: dict, ncbi: dict | None = None, ensembl: dict | None = None) -> dict:
    return {
        "bundle": bundle,
        "ucsc_references": ucsc,
        "ncbi_references": ncbi or {},
        "ensembl_references": ensembl or {},
    }


def _bundle_pointer(*, common_asset: str, common_asset_sha256: str, repository: str = "acme/vntyper-data") -> dict:
    return {
        "repository": repository,
        "release_tag": "refs-v1",
        "common_asset": common_asset,
        "common_asset_sha256": common_asset_sha256,
    }


def _fake_download_from(release_dir: Path):
    """A `download_file` stand-in that copies from `release_dir` instead of the network.

    Mirrors the real function's contract closely enough for `install_from_bundle` to
    behave identically under it: a missing asset raises `SystemExit`, matching what
    `download_file` does when `urlretrieve` fails.

    Args:
        release_dir: Directory the fake "release" assets live in.

    Returns:
        Callable[[str, Path], None]: A `download_file(url, dest_path)` stand-in.
    """
    import shutil

    def _download(url: str, dest_path: Path) -> None:
        source = release_dir / Path(url).name
        if not source.is_file():
            raise SystemExit(1)
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(source, dest_path)

    return _download


COMMON_ASSET_NAME = "vntyper-references-refs-v1-muc1.tar.gz"


# =============================================================================
# Happy path
# =============================================================================


class TestInstallFromBundleHappyPath:
    def test_the_selected_genome_is_installed_an_unselected_one_is_not_and_the_common_asset_always_is(
        self, tmp_path, monkeypatch
    ):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"

        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19")
        hg38_asset = "vntyper-references-refs-v1-ucsc-hg38.tar.gz"
        _write_bundle_asset(release_dir, hg38_asset, _genome_files("hg38"), reference_id="hg38")
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)

        config = _install_config(
            ucsc={
                "hg19": _genome_entry("hg19", hg19_asset, hg19_sha256),
                "hg38": _genome_entry("hg38", hg38_asset, "0" * 64),  # wrong on purpose: must never be fetched
            },
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))
        monkeypatch.setattr(install_references, "_local_bwa_version", lambda: "9.9.9")

        install_references.install_from_bundle(config, output_dir, ["hg19"])

        assert (output_dir / "alignment" / "chr1.hg19.fa").read_bytes() == b">chr1\nACGTACGTAC\n"
        for ext in (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"):
            assert (output_dir / "alignment" / f"chr1.hg19.fa{ext}").exists()

        assert not (output_dir / "alignment" / "chr1.hg38.fa").exists(), "an unselected genome must not be fetched"

        # "--references hg19 -> the common asset is still installed" (brief bullet).
        assert (output_dir / "MUC1_motifs_Rev_com.fa").read_bytes() == b">motif\nACGT\n"
        assert (output_dir / "code-adVNTR_RUs.fa").exists()
        assert (output_dir / "vntr_db_advntr" / "hg19_muc1.db").read_bytes() == b"hg19-db-bytes"
        assert (output_dir / "vntr_db_advntr" / "hg38_muc1.db").read_bytes() == b"hg38-db-bytes"

        # The two metadata members are read and then removed - they describe one asset
        # each and must not leak into, or be mistaken for, the installed reference tree.
        assert not (output_dir / "release-manifest.json").exists()
        assert not (output_dir / "BUILD_INFO.json").exists()


# =============================================================================
# asset_sha256 mismatch
# =============================================================================


class TestInstallFromBundleDigestMismatch:
    def test_a_mismatch_names_the_asset_and_both_digests_and_installs_nothing(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        _write_bundle_asset(release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19")
        actual_sha256 = _sha256_file(release_dir / hg19_asset)
        wrong_digest = "f" * 64
        assert wrong_digest != actual_sha256

        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, wrong_digest)},
            bundle=_bundle_pointer(common_asset="unused.tar.gz", common_asset_sha256="1" * 64),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))

        with pytest.raises(ValueError) as excinfo:
            install_references.install_from_bundle(config, output_dir, ["hg19"])

        message = str(excinfo.value)
        assert hg19_asset in message
        assert wrong_digest in message
        assert actual_sha256 in message
        assert not output_dir.exists()

    def test_the_common_asset_is_checked_before_extraction_too(self, tmp_path, monkeypatch):
        """Carry-forward from the Task 3 review: `verify_sha256` runs before
        `safe_extract` for *every* asset, common included, not only the first one a run
        happens to process."""
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19")
        _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)
        actual_common_sha256 = _sha256_file(release_dir / COMMON_ASSET_NAME)
        wrong_digest = "e" * 64
        assert wrong_digest != actual_common_sha256

        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=wrong_digest),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))

        with pytest.raises(ValueError, match=COMMON_ASSET_NAME):
            install_references.install_from_bundle(config, output_dir, ["hg19"])

        # The genome asset that verified fine must not have been left activated either:
        # the whole run is one staged install.
        assert not output_dir.exists()


# =============================================================================
# Malformed in-archive metadata
#
# `release-manifest.json` and `BUILD_INFO.json` live inside an already
# `asset_sha256`-verified archive, so nothing outside this module can corrupt them - but
# a hand-assembled or future-format asset could still omit one, ship invalid JSON, or
# describe a file that was not actually extracted. Every case below still verifies its
# whole-archive digest against a real `asset_sha256`; only the *contents* are malformed.
# =============================================================================


class TestInstallFromBundleMalformedMetadata:
    def test_an_asset_missing_release_manifest_entirely_is_refused(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        members = {**_genome_files("hg19"), "BUILD_INFO.json": json.dumps({"bwa_version": "9.9.9"}).encode("utf-8")}
        hg19_sha256 = _write_raw_asset(release_dir, hg19_asset, members)
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)

        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))

        with pytest.raises(ValueError, match="release-manifest.json"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])
        assert not output_dir.exists()

    def test_an_asset_whose_metadata_is_not_valid_json_is_refused(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        members = {
            **_genome_files("hg19"),
            "release-manifest.json": b"{not valid json",
            "BUILD_INFO.json": b"{}",
        }
        hg19_sha256 = _write_raw_asset(release_dir, hg19_asset, members)
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)

        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))

        with pytest.raises(ValueError, match="not valid JSON"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])
        assert not output_dir.exists()

    def test_a_manifest_naming_a_file_extraction_did_not_produce_is_refused(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        bogus_manifest = json.dumps(
            {
                "asset": hg19_asset,
                "reference_id": "hg19",
                "files": [{"path": "alignment/chr1.hg19.fa.missing", "size": 1, "sha256": "a" * 64}],
            }
        ).encode("utf-8")
        members = {**_genome_files("hg19"), "release-manifest.json": bogus_manifest, "BUILD_INFO.json": b"{}"}
        hg19_sha256 = _write_raw_asset(release_dir, hg19_asset, members)
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)

        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))

        with pytest.raises(ValueError, match="chr1.hg19.fa.missing"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])
        assert not output_dir.exists()


# =============================================================================
# BUILD_INFO.json bwa version mismatch
# =============================================================================


class TestInstallFromBundleToolchainMismatch:
    def test_a_bwa_version_mismatch_logs_a_warning_and_reindexes(self, tmp_path, monkeypatch, caplog):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(
            release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19", bwa_version="0.7.17"
        )
        common_sha256 = _write_bundle_asset(
            release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None, bwa_version="0.7.17"
        )
        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))
        monkeypatch.setattr(install_references, "_local_bwa_version", lambda: "0.7.18-r1243")

        reindex_calls: list[list[str]] = []

        def _fake_run(argv, **kwargs):
            reindex_calls.append(list(argv))
            return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

        monkeypatch.setattr(install_references.subprocess, "run", _fake_run)

        with caplog.at_level(logging.WARNING, logger="vntyper.scripts.install_references"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])

        warnings = [m for m in caplog.messages if "bwa" in m.lower() and "0.7.17" in m]
        assert warnings, "the bundled and local bwa versions must both be named in the warning"
        assert any("0.7.18-r1243" in m for m in warnings)

        assert len(reindex_calls) == 1
        assert reindex_calls[0][:2] == ["bwa", "index"]
        # Re-indexed in the staging tree, before activation - not against `output_dir`
        # itself, which does not hold the new bytes until the whole install succeeds.
        assert reindex_calls[0][-1].endswith("/alignment/chr1.hg19.fa")
        assert reindex_calls[0][-1] != str(output_dir / "alignment" / "chr1.hg19.fa")

    def test_a_matching_bwa_version_does_not_reindex(self, tmp_path, monkeypatch, caplog):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(
            release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19", bwa_version="9.9.9"
        )
        common_sha256 = _write_bundle_asset(
            release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None, bwa_version="9.9.9"
        )
        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))
        monkeypatch.setattr(install_references, "_local_bwa_version", lambda: "9.9.9")
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda *a, **k: pytest.fail("bwa must not be re-run when the versions already match"),
        )

        with caplog.at_level(logging.WARNING, logger="vntyper.scripts.install_references"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])

        assert not [m for m in caplog.messages if "re-index" in m.lower()]

    def test_a_non_bwa_kind_never_triggers_the_bwa_version_check(self, tmp_path, monkeypatch, caplog):
        """Confirms `kind` is actually consumed, not merely written: a hypothetical
        future non-bwa entry must not be re-indexed with `bwa index`, however the local
        and bundled bwa versions compare."""
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(
            release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19", bwa_version="0.7.17"
        )
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)
        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256, kind="none")},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))
        monkeypatch.setattr(install_references, "_local_bwa_version", lambda: "0.7.18-r1243")
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda *a, **k: pytest.fail("a non-bwa-kind entry must never be re-indexed with bwa"),
        )

        with caplog.at_level(logging.WARNING, logger="vntyper.scripts.install_references"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])

        assert not [m for m in caplog.messages if "re-index" in m.lower()]
        assert (output_dir / "alignment" / "chr1.hg19.fa").exists()

    def test_bwa_absent_locally_is_treated_as_unknown_not_as_a_mismatch(self, tmp_path, monkeypatch):
        """`_local_bwa_version` returns None when bwa is not on PATH (an OSError from
        `subprocess.run`). That must not be treated as a positive mismatch - the whole
        point is to react to a *confirmed* difference, not to the absence of a bwa to
        compare against."""
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(
            release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19", bwa_version="0.7.17"
        )
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)
        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda argv, **k: (
                (_ for _ in ()).throw(OSError("bwa not found"))
                if argv == ["bwa"]
                else pytest.fail(f"unexpected subprocess.run call: {argv}")
            ),
        )

        install_references.install_from_bundle(config, output_dir, ["hg19"])

        assert (output_dir / "alignment" / "chr1.hg19.fa").exists()

    def test_a_mismatch_with_no_index_command_configured_is_refused(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(
            release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19", bwa_version="0.7.17"
        )
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)
        entry = _genome_entry("hg19", hg19_asset, hg19_sha256)
        del entry["index_command"]
        config = _install_config(
            ucsc={"hg19": entry},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))
        monkeypatch.setattr(install_references, "_local_bwa_version", lambda: "0.7.18-r1243")

        with pytest.raises(ValueError, match="index_command"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])
        assert not output_dir.exists()

    def test_a_failed_reindex_raises_and_the_failure_message_survives(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(
            release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19", bwa_version="0.7.17"
        )
        common_sha256 = _write_bundle_asset(release_dir, COMMON_ASSET_NAME, _common_files(), reference_id=None)
        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256=common_sha256),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))
        monkeypatch.setattr(install_references, "_local_bwa_version", lambda: "0.7.18-r1243")
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda argv, **k: subprocess.CompletedProcess(argv, 1, stdout="", stderr="disk full"),
        )

        with pytest.raises(RuntimeError, match="disk full"):
            install_references.install_from_bundle(config, output_dir, ["hg19"])
        assert not output_dir.exists()


# =============================================================================
# Download failure
# =============================================================================


class TestInstallFromBundleDownloadFailure:
    def test_a_failed_download_names_the_repository_tag_and_asset(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        release_dir.mkdir()  # exists, but the asset itself is never written: a 404.
        output_dir = tmp_path / "refs"
        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"

        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, "a" * 64)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256="b" * 64),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))

        with pytest.raises(RuntimeError) as excinfo:
            install_references.install_from_bundle(config, output_dir, ["hg19"])

        message = str(excinfo.value)
        assert "acme/vntyper-data" in message
        assert "refs-v1" in message
        assert hg19_asset in message
        assert not output_dir.exists()


# =============================================================================
# Interruption mid-install
# =============================================================================


class TestInstallFromBundleInterruption:
    def test_an_interruption_after_the_first_asset_leaves_the_previous_tree_intact(self, tmp_path, monkeypatch):
        release_dir = tmp_path / "release"
        output_dir = tmp_path / "refs"
        output_dir.mkdir()
        (output_dir / "alignment").mkdir()
        (output_dir / "alignment" / "chr1.hg19.fa").write_text("previous install", encoding="utf-8")

        hg19_asset = "vntyper-references-refs-v1-ucsc-hg19.tar.gz"
        hg19_sha256 = _write_bundle_asset(release_dir, hg19_asset, _genome_files("hg19"), reference_id="hg19")
        # The common asset is never written to release_dir: its fetch 404s, simulating an
        # interruption after the first (genome) asset has already verified and extracted.

        config = _install_config(
            ucsc={"hg19": _genome_entry("hg19", hg19_asset, hg19_sha256)},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256="c" * 64),
        )
        monkeypatch.setattr(install_references, "download_file", _fake_download_from(release_dir))

        with pytest.raises(RuntimeError):
            install_references.install_from_bundle(config, output_dir, ["hg19"])

        assert (output_dir / "alignment" / "chr1.hg19.fa").read_text(encoding="utf-8") == "previous install"
        assert not (output_dir / "MUC1_motifs_Rev_com.fa").exists()
        assert not list(tmp_path.glob(f".{output_dir.name}.staging.*"))
        assert not list(tmp_path.glob(f".{output_dir.name}.previous.*"))


# =============================================================================
# Configuration errors
# =============================================================================


class TestInstallFromBundleConfiguration:
    def test_a_missing_bundle_section_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="bundle"):
            install_references.install_from_bundle({"ucsc_references": {}}, tmp_path / "refs", ["hg19"])

    def test_a_selected_genome_with_no_asset_configured_is_refused(self, tmp_path):
        config = _install_config(
            ucsc={"hg19": {"kind": "bwa", "installed_path": "alignment/chr1.hg19.fa"}},
            bundle=_bundle_pointer(common_asset=COMMON_ASSET_NAME, common_asset_sha256="a" * 64),
        )
        with pytest.raises(ValueError, match="hg19"):
            install_references.install_from_bundle(config, tmp_path / "refs", ["hg19"])
