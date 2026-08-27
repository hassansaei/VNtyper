"""Seed archives must stay inside staging and record only extracted files."""

import hashlib
import io
import json
import tarfile
import zipfile

import pytest

from vntyper.scripts import install_references

pytestmark = pytest.mark.unit


def _zip_bytes(name: str = "vntr_db_advntr_v2/hg19_muc1.db") -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(name, "hg19 db")
    return buffer.getvalue()


def _tar_bytes(name: str = "vntr_db_advntr_v2/hg19_muc1.db") -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        payload = b"hg19 db"
        member = tarfile.TarInfo(name)
        member.size = len(payload)
        member.mode = 0o644
        archive.addfile(member, io.BytesIO(payload))
    return buffer.getvalue()


def _config(payload: bytes, *, suffix: str, extract_to: str = "database") -> dict[str, dict[str, str]]:
    return {
        "vntr_db_advntr": {
            "url": f"https://seed.invalid/vntr_db_advntr{suffix}",
            "target_path": f"vntr_db_advntr{suffix}",
            "source_sha256": hashlib.sha256(payload).hexdigest(),
            "extract_to": extract_to,
        }
    }


def _download_payload(monkeypatch: pytest.MonkeyPatch, payload: bytes) -> None:
    def _download(_url, destination):
        destination.write_bytes(payload)
        return True

    monkeypatch.setattr(install_references, "download_file", _download)


def test_the_zip_branch_extracts_and_records_regular_members(tmp_path, monkeypatch):
    payload = _zip_bytes()
    _download_payload(monkeypatch, payload)

    install_references.process_vntyper_references(_config(payload, suffix=".zip"), tmp_path, "bwa", True, {})

    assert (tmp_path / "database/vntr_db_advntr_v2/hg19_muc1.db").read_bytes() == b"hg19 db"
    provenance = json.loads((tmp_path / "install_provenance.json").read_text(encoding="utf-8"))
    assert "database/vntr_db_advntr_v2/hg19_muc1.db" in provenance


def test_the_zip_branch_refuses_parent_traversal_without_recording_or_writing_it(tmp_path, monkeypatch):
    payload = _zip_bytes("../escaped.txt")
    _download_payload(monkeypatch, payload)

    with pytest.raises(RuntimeError, match="Failed to extract") as excinfo:
        install_references.process_vntyper_references(_config(payload, suffix=".zip"), tmp_path, "bwa", True, {})

    assert isinstance(excinfo.value.__cause__, ValueError)
    assert not (tmp_path / "escaped.txt").exists()
    provenance = json.loads((tmp_path / "install_provenance.json").read_text(encoding="utf-8"))
    assert "escaped.txt" not in provenance


def test_the_zip_branch_refuses_a_preexisting_symlink_escape(tmp_path, monkeypatch):
    payload = _zip_bytes("data/file.txt")
    _download_payload(monkeypatch, payload)
    outside = tmp_path / "outside"
    outside.mkdir()
    extract_dir = tmp_path / "database"
    extract_dir.mkdir()
    (extract_dir / "data").symlink_to(outside, target_is_directory=True)

    with pytest.raises(RuntimeError, match="Failed to extract") as excinfo:
        install_references.process_vntyper_references(_config(payload, suffix=".zip"), tmp_path, "bwa", True, {})

    assert isinstance(excinfo.value.__cause__, ValueError)
    assert not (outside / "file.txt").exists()


def test_the_zip_branch_refuses_a_symlinked_extract_target(tmp_path, monkeypatch):
    payload = _zip_bytes()
    _download_payload(monkeypatch, payload)
    outside = tmp_path / "outside"
    outside.mkdir()
    (tmp_path / "database").symlink_to(outside, target_is_directory=True)

    with pytest.raises(RuntimeError, match="Failed to extract") as excinfo:
        install_references.process_vntyper_references(_config(payload, suffix=".zip"), tmp_path, "bwa", True, {})

    assert isinstance(excinfo.value.__cause__, ValueError)
    assert not any(outside.iterdir())


def test_a_refused_symlink_escape_never_activates_the_staged_tree(tmp_path, monkeypatch):
    payload = _zip_bytes()
    _download_payload(monkeypatch, payload)
    output_dir = tmp_path / "refs"
    output_dir.mkdir()
    (output_dir / "sentinel.txt").write_text("previous install", encoding="utf-8")
    outside = tmp_path / "outside"
    outside.mkdir()
    (output_dir / "database").symlink_to(outside, target_is_directory=True)
    config = {"vntyper_references": _config(payload, suffix=".zip")}

    with pytest.raises(RuntimeError, match="Failed to extract"):
        install_references._install_references(
            install_config=config,
            output_dir=output_dir,
            config_path=None,
            skip_indexing=True,
            index_threads=1,
            aligners_to_use=None,
            found_refs=set(),
            from_source=True,
            release_spec_path=None,
            aligner_config={},
        )

    assert (output_dir / "sentinel.txt").read_text(encoding="utf-8") == "previous install"
    assert (output_dir / "database").is_symlink()
    assert not any(outside.iterdir())
    assert not (output_dir / "vntr_db_advntr.zip").exists()
    assert not (output_dir / "install_provenance.json").exists()


def test_the_tar_branch_routes_through_containment_and_wraps_refusal(tmp_path, monkeypatch):
    payload = _tar_bytes("../escaped.txt")
    _download_payload(monkeypatch, payload)

    with pytest.raises(RuntimeError, match="Failed to extract") as excinfo:
        install_references.process_vntyper_references(_config(payload, suffix=".tar.gz"), tmp_path, "bwa", True, {})

    assert isinstance(excinfo.value.__cause__, ValueError)
    assert not (tmp_path / "escaped.txt").exists()


def test_the_tar_branch_records_the_regular_members_returned_by_safe_extraction(tmp_path, monkeypatch):
    payload = _tar_bytes()
    _download_payload(monkeypatch, payload)

    install_references.process_vntyper_references(
        _config(payload, suffix=".tar.gz", extract_to="."), tmp_path, "bwa", True, {}
    )

    provenance = json.loads((tmp_path / "install_provenance.json").read_text(encoding="utf-8"))
    assert "vntr_db_advntr_v2/hg19_muc1.db" in provenance
    assert "vntr_db_advntr_v2" not in provenance
