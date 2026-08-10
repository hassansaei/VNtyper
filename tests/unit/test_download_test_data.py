import hashlib
import json
import sys
import zipfile
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace

import pytest

import scripts.download_test_data as dtd

pytestmark = pytest.mark.unit


def test_compute_md5_hashes_file_in_chunks(tmp_path: Path) -> None:
    path = tmp_path / "payload"
    path.write_bytes(b"abc" * 30000)
    assert dtd.compute_md5(path) == hashlib.md5(path.read_bytes()).hexdigest()


@pytest.mark.parametrize(
    ("members", "expected", "absent"),
    [
        ({**{f"data/f{i}.txt": str(i) for i in range(10)}, "README": "r"}, "f0.txt", "data"),
        ({"a/x": "x", "b/y": "y", "root": "r"}, "a/x", "not-stripped"),
    ],
)
def test_extract_archive_handles_dominant_and_mixed_layouts(
    tmp_path: Path, members: dict[str, str], expected: str, absent: str
) -> None:
    archive = tmp_path / "data.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        for name, value in members.items():
            handle.writestr(name, value)
    output = tmp_path / "out"
    dtd.extract_archive(archive, output)
    if absent == "data":
        assert (output / expected).read_text() == "0"
        assert not (output / "data").exists()
    else:
        assert (output / expected).read_text() == "x"
        assert (output / "b/y").read_text() == "y"


@pytest.mark.parametrize("dominant", [True, False])
def test_extract_archive_preserves_safe_empty_directories_in_every_layout(tmp_path: Path, dominant: bool) -> None:
    archive = tmp_path / "directories.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        if dominant:
            handle.writestr("data/empty/nested/", "")
            for index in range(10):
                handle.writestr(f"data/file-{index}", str(index))
            handle.writestr("README", "metadata")
        else:
            handle.writestr("a/empty/", "")
            handle.writestr("a/file", "a")
            handle.writestr("b/file", "b")
            handle.writestr("root", "root")

    output = tmp_path / "out"
    dtd.extract_archive(archive, output)

    expected = output / ("empty/nested" if dominant else "a/empty")
    assert expected.is_dir()


@pytest.mark.parametrize("member", ["data/../../escaped", "data/link/payload"])
def test_extract_archive_rejects_members_that_escape_destination(tmp_path: Path, member: str) -> None:
    archive = tmp_path / "data.zip"
    output = tmp_path / "out"
    output.mkdir()
    if member == "data/link/payload":
        (output / "link").symlink_to(tmp_path / "outside", target_is_directory=True)
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(member, "owned")
    with pytest.raises(ValueError, match="archive member escapes extraction directory"):
        dtd.extract_archive(archive, output)
    assert not (tmp_path / "escaped").exists()
    assert not (tmp_path / "outside/payload").exists()


def test_mixed_layout_rejects_preexisting_symlink_escape(tmp_path: Path) -> None:
    archive = tmp_path / "mixed.zip"
    output = tmp_path / "out"
    output.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (output / "link").symlink_to(outside, target_is_directory=True)
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("link/payload", "owned")
        handle.writestr("other/value", "ordinary")
    with pytest.raises(ValueError, match="archive member escapes extraction directory"):
        dtd.extract_archive(archive, output)
    assert not (outside / "payload").exists()


@pytest.mark.parametrize("relative", [Path(), Path("/absolute")])
def test_confined_archive_target_rejects_empty_and_absolute_paths(tmp_path: Path, relative: Path) -> None:
    with pytest.raises(ValueError, match="archive member escapes extraction directory"):
        dtd._confined_archive_target(tmp_path, relative)


def test_verify_reports_missing_and_md5_mismatch(tmp_path: Path) -> None:
    data = tmp_path / "data"
    data.mkdir()
    (data / "bad.txt").write_text("bad")
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "file_resources": [
                    {"local_path": "tests/data", "filename": "missing.txt", "md5sum": "0"},
                    {"local_path": "tests/data", "filename": "bad.txt", "md5sum": "deadbeef"},
                ]
            }
        )
    )
    ok, errors = dtd.verify_test_data(config, data)
    assert ok is False
    assert errors[0] == f"Missing: {data / 'missing.txt'}"
    assert errors[1].startswith(f"MD5 mismatch: {data / 'bad.txt'}")


def test_download_rejects_a_short_response(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    response = SimpleNamespace(
        headers={"content-length": "4"},
        raise_for_status=lambda: None,
        iter_content=lambda chunk_size: iter([b"abc"]),
    )
    monkeypatch.setitem(sys.modules, "requests", SimpleNamespace(get=lambda *a, **k: response))
    with pytest.raises(RuntimeError, match="Incomplete download: 3/4 bytes"):
        dtd.download_file_requests("https://example.invalid/data.zip", tmp_path / "out")


def test_verify_only_returns_existing_verification_status(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"url": "https://example.invalid/x"}}))
    monkeypatch.setattr(dtd, "__file__", str(script))
    monkeypatch.setattr(dtd, "verify_test_data", lambda *a, **k: (False, ["missing"]))
    assert dtd.main(["--verify-only"]) == 1


def test_verify_skip_md5_accepts_wrong_checksum(tmp_path: Path) -> None:
    data = tmp_path / "data"
    data.mkdir()
    (data / "present.txt").write_text("payload")
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "file_resources": [
                    {"local_path": "tests/data", "filename": "present.txt", "md5sum": "deliberately-wrong"}
                ]
            }
        )
    )
    assert dtd.verify_test_data(config, data, skip_md5=True) == (True, [])


def test_extract_archive_propagates_bad_zip(tmp_path: Path) -> None:
    archive = tmp_path / "bad.zip"
    archive.write_bytes(b"not a zip archive")
    with pytest.raises(zipfile.BadZipFile):
        dtd.extract_archive(archive, tmp_path / "out")


def test_main_returns_one_when_config_is_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    monkeypatch.setattr(dtd, "__file__", str(script))
    assert dtd.main([]) == 1


def test_main_returns_one_when_archive_config_is_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text("{}")
    monkeypatch.setattr(dtd, "__file__", str(script))
    assert dtd.main([]) == 1


@pytest.mark.parametrize("contents", ["{broken", "[]"])
def test_main_returns_one_and_diagnoses_invalid_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture, contents: str
) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(contents, encoding="utf-8")
    monkeypatch.setattr(dtd, "__file__", str(script))

    with caplog.at_level("ERROR", logger=dtd.__name__):
        assert dtd.main([]) == 1
    assert f"Could not load config file {config}:" in caplog.text


def test_main_returns_one_and_diagnoses_missing_archive_url(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"md5sum": "unused"}}), encoding="utf-8")
    monkeypatch.setattr(dtd, "__file__", str(script))

    with caplog.at_level("ERROR", logger=dtd.__name__):
        assert dtd.main(["--force"]) == 1
    assert "Archive URL is missing or invalid in test_data_config.json" in caplog.text


def test_main_returns_one_and_diagnoses_temporary_archive_creation_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"url": "https://example.invalid/data.zip"}}), encoding="utf-8")
    monkeypatch.setattr(dtd, "__file__", str(script))

    def fail_temporary_file(**kwargs: object) -> None:
        del kwargs
        raise OSError("temporary storage unavailable")

    monkeypatch.setattr(dtd.tempfile, "NamedTemporaryFile", fail_temporary_file)
    with caplog.at_level("ERROR", logger=dtd.__name__):
        assert dtd.main(["--force"]) == 1
    assert "Could not create temporary archive: temporary storage unavailable" in caplog.text


def test_main_returns_one_and_diagnoses_temporary_archive_cleanup_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"url": "https://example.invalid/data.zip"}}), encoding="utf-8")
    temporary_archive = tmp_path / "named-temporary.zip"
    temporary_archive.touch()
    original_unlink = Path.unlink

    def fail_cleanup(path: Path, missing_ok: bool = False) -> None:
        if path == temporary_archive:
            raise OSError("cleanup denied")
        original_unlink(path, missing_ok=missing_ok)

    monkeypatch.setattr(dtd, "__file__", str(script))
    monkeypatch.setattr(
        dtd.tempfile,
        "NamedTemporaryFile",
        lambda **kwargs: nullcontext(SimpleNamespace(name=str(temporary_archive))),
    )
    monkeypatch.setattr(dtd, "download_file_requests", lambda *args, **kwargs: None)
    monkeypatch.setattr(dtd, "extract_archive", lambda *args, **kwargs: None)
    monkeypatch.setattr(dtd, "verify_test_data", lambda *args, **kwargs: (True, []))
    monkeypatch.setattr(Path, "unlink", fail_cleanup)

    with caplog.at_level("ERROR", logger=dtd.__name__):
        assert dtd.main(["--force"]) == 1
    assert f"Could not remove temporary archive {temporary_archive}: cleanup denied" in caplog.text


def test_force_download_succeeds_and_removes_temporary_archive(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"url": "https://example.invalid/data.zip"}}))
    temporary_archive = tmp_path / "named-temporary.zip"
    calls: list[tuple[str, Path]] = []

    def fake_download(url: str, destination: Path, timeout: int = 300) -> None:
        del timeout
        calls.append((url, destination))
        destination.write_bytes(b"archive")

    def fake_extract(archive: Path, output: Path) -> None:
        calls.append(("extract", archive))
        output.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(dtd, "__file__", str(script))
    monkeypatch.setattr(
        dtd.tempfile,
        "NamedTemporaryFile",
        lambda **kwargs: nullcontext(SimpleNamespace(name=str(temporary_archive))),
    )
    monkeypatch.setattr(dtd, "download_file_requests", fake_download)
    monkeypatch.setattr(dtd, "extract_archive", fake_extract)
    monkeypatch.setattr(dtd, "verify_test_data", lambda *a, **k: (True, []))

    assert dtd.main(["--force"]) == 0
    assert calls == [
        ("https://example.invalid/data.zip", temporary_archive),
        ("extract", temporary_archive),
    ]
    assert not temporary_archive.exists()


def test_post_extraction_verification_failure_removes_temporary_archive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"url": "https://example.invalid/data.zip"}}))
    temporary_archive = tmp_path / "named-temporary.zip"

    def fake_download(url: str, destination: Path, timeout: int = 300) -> None:
        del url, timeout
        destination.write_bytes(b"archive")

    monkeypatch.setattr(dtd, "__file__", str(script))
    monkeypatch.setattr(
        dtd.tempfile,
        "NamedTemporaryFile",
        lambda **kwargs: nullcontext(SimpleNamespace(name=str(temporary_archive))),
    )
    monkeypatch.setattr(dtd, "download_file_requests", fake_download)
    monkeypatch.setattr(dtd, "extract_archive", lambda *a, **k: None)
    monkeypatch.setattr(dtd, "verify_test_data", lambda *a, **k: (False, ["missing"]))

    assert dtd.main(["--force"]) == 1
    assert not temporary_archive.exists()


def test_download_exception_returns_one_and_removes_temporary_archive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"url": "https://example.invalid/data.zip"}}))
    temporary_archive = tmp_path / "named-temporary.zip"
    temporary_archive.touch()

    def fail_download(url: str, destination: Path, timeout: int = 300) -> None:
        del url, destination, timeout
        raise RuntimeError("download failed")

    monkeypatch.setattr(dtd, "__file__", str(script))
    monkeypatch.setattr(
        dtd.tempfile,
        "NamedTemporaryFile",
        lambda **kwargs: nullcontext(SimpleNamespace(name=str(temporary_archive))),
    )
    monkeypatch.setattr(dtd, "download_file_requests", fail_download)

    assert dtd.main(["--force"]) == 1
    assert not temporary_archive.exists()
