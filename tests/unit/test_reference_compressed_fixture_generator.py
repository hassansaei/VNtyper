"""Unit tests for reference-compressed CRAM generation and artifact publication."""

from __future__ import annotations

import hashlib
import io
import subprocess
import sys
from pathlib import Path
from unittest import mock

import pytest

import scripts.make_cram_fixtures as generator
from scripts.cram_reference_contract import (
    REGISTERED_B178_INDEXED_REGION_DIGEST,
    REGISTERED_B178_INDEXED_REGION_RECORDS,
    LossyConversionError,
)

pytestmark = pytest.mark.unit


def _touch(path: Path, content: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


class _ReferenceSamtools:
    """Materialize deterministic staged artifacts while recording every command."""

    def __init__(self, *, failure: str | None = None, create_index: bool = True) -> None:
        self.failure = failure
        self.create_index = create_index
        self.commands: list[list[str]] = []

    def __call__(self, argv: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        self.commands.append(argv)
        if argv[1:4] == ["view", "-H", "--no-PG"]:
            operation = "header"
        elif argv[1] == "reheader":
            operation = "reheader"
        elif "-C" in argv:
            operation = "encode"
        else:
            operation = "index"
        if operation == self.failure:
            return subprocess.CompletedProcess(argv, 1, stdout="", stderr=f"{operation} failed")
        if operation == "header":
            header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:249250621\n@SQ\tSN:chr2\tLN:243199373\n"
            return subprocess.CompletedProcess(argv, 0, stdout=header, stderr="")
        if operation == "reheader":
            output = kwargs["stdout"]
            assert hasattr(output, "name")
            Path(output.name).write_bytes(b"reheadered bam")
        elif operation == "encode":
            Path(argv[argv.index("-o") + 1]).write_bytes(b"reference-compressed cram")
        elif operation == "index" and self.create_index:
            Path(f"{argv[-1]}.crai").write_bytes(b"reference-compressed crai")
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")


class _DigestProcess:
    def __init__(self, stdout: str, *, returncode: int = 0) -> None:
        self.stdout = io.StringIO(stdout)
        self.returncode = returncode

    def __enter__(self) -> _DigestProcess:
        return self

    def __exit__(self, *args: object) -> None:
        del args

    def wait(self) -> None:
        return None


def _registered_digest(
    _samtools: str,
    _alignment: Path,
    _explicit_reference: Path | None = None,
    *,
    region: str | None = None,
) -> tuple[str, int]:
    if region is not None:
        return REGISTERED_B178_INDEXED_REGION_DIGEST, REGISTERED_B178_INDEXED_REGION_RECORDS
    return "normalized-record-digest", 34_214


def test_normalized_digest_builds_the_explicit_indexed_query(monkeypatch: pytest.MonkeyPatch) -> None:
    record = "read-1\t0\tchr1\t155160500\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:4\n"
    popen = mock.Mock(return_value=_DigestProcess(record))
    monkeypatch.setattr(generator.subprocess, "Popen", popen)
    reference = Path("reference/alignment/chr1.hg19.fa")

    digest, count = generator._normalized_record_digest(
        "samtools", Path("sample.cram"), reference, region=generator.REFERENCE_VALIDATION_REGION
    )

    normalized = b"read-1\t0\tchr1\t155160500\t60\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:4\tNM:i:0\n"
    assert (digest, count) == (hashlib.sha256(normalized).hexdigest(), 1)
    assert popen.call_args.args[0] == [
        "samtools",
        "view",
        "-T",
        str(reference),
        "sample.cram",
        generator.REFERENCE_VALIDATION_REGION,
    ]


def test_normalized_digest_reports_decode_stderr(monkeypatch: pytest.MonkeyPatch) -> None:
    def popen(*_args: object, **kwargs: object) -> _DigestProcess:
        stderr_file = kwargs["stderr"]
        assert hasattr(stderr_file, "write")
        stderr_file.write("reference unavailable")
        return _DigestProcess("partial\n", returncode=3)

    monkeypatch.setattr(generator.subprocess, "Popen", popen)

    with pytest.raises(RuntimeError, match="exited 3: reference unavailable"):
        generator._normalized_record_digest("samtools", Path("sample.cram"), Path("reference.fa"))


@pytest.mark.timeout(5)
def test_stream_digest_does_not_deadlock_when_stderr_exceeds_pipe_capacity() -> None:
    script = (
        "import signal,sys; signal.alarm(3); "
        "sys.stderr.write('diagnostic-' * 200000); sys.stderr.flush(); "
        "sys.stdout.write('record\\n')"
    )

    digest, count = generator._stream_record_digest([sys.executable, "-c", script], lambda line: line.encode())

    assert (digest, count) == (hashlib.sha256(b"record\n").hexdigest(), 1)


@pytest.mark.parametrize("missing", ["absent", "wrong"])
def test_reference_identity_is_checked_before_encoding(tmp_path: Path, missing: str) -> None:
    source = _touch(tmp_path / "source.bam")
    reference = tmp_path / "reference.fa"
    if missing == "wrong":
        _touch(reference, ">chr1\nA\n")
        expected = "SHA-256 mismatch"
    else:
        expected = "does not exist"

    with pytest.raises((FileNotFoundError, ValueError), match=expected):
        generator.derive_reference_compressed_cram("samtools", source, reference, tmp_path / "cram")


def test_derivation_stages_then_publishes_explicit_reference_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _touch(tmp_path / "source.bam", "source bam")
    reference = _touch(tmp_path / "reference.fa", ">chr1\nACGT\n")
    root = tmp_path / "cram"
    fake = _ReferenceSamtools()
    digest = mock.Mock(side_effect=_registered_digest)
    monkeypatch.setattr(generator.subprocess, "run", fake)
    monkeypatch.setattr(generator, "_normalized_record_digest", digest)

    fixture = generator.derive_reference_compressed_cram(
        "samtools",
        source,
        reference,
        root,
        expected_reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
    )

    final_cram = root / "reference-compressed/source.cram"
    staged_cram = Path(fake.commands[2][fake.commands[2].index("-o") + 1])
    assert staged_cram != final_cram
    assert fake.commands[2][0:6] == ["samtools", "view", "-C", "-T", str(reference), "--no-PG"]
    assert fake.commands[3] == ["samtools", "index", str(staged_cram)]
    assert fixture.cram == final_cram
    assert final_cram.read_bytes() == b"reference-compressed cram"
    assert Path(f"{final_cram}.crai").read_bytes() == b"reference-compressed crai"
    assert fixture.indexed_region_records == REGISTERED_B178_INDEXED_REGION_RECORDS
    assert fixture.source_indexed_region_digest == REGISTERED_B178_INDEXED_REGION_DIGEST


def test_index_failure_aborts_before_decode_or_publication(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = _touch(tmp_path / "source.bam")
    reference = _touch(tmp_path / "reference.fa", ">chr1\nACGT\n")
    fake = _ReferenceSamtools(failure="index")
    digest = mock.Mock()
    monkeypatch.setattr(generator.subprocess, "run", fake)
    monkeypatch.setattr(generator, "_normalized_record_digest", digest)

    with pytest.raises(RuntimeError, match="index failed"):
        generator.derive_reference_compressed_cram(
            "samtools",
            source,
            reference,
            tmp_path / "cram",
            expected_reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        )

    digest.assert_not_called()
    assert not (tmp_path / "cram/reference-compressed/source.cram").exists()


def test_missing_crai_aborts_before_decode_or_publication(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A zero-exit index command cannot validate or preserve a stale final CRAI."""
    source = _touch(tmp_path / "source.bam")
    reference = _touch(tmp_path / "reference.fa", ">chr1\nACGT\n")
    root = tmp_path / "cram"
    final_crai = _touch(root / "reference-compressed/source.cram.crai", "stale crai")
    fake = _ReferenceSamtools(create_index=False)
    digest = mock.Mock()
    monkeypatch.setattr(generator.subprocess, "run", fake)
    monkeypatch.setattr(generator, "_normalized_record_digest", digest)

    with pytest.raises(RuntimeError, match="did not create expected CRAI"):
        generator.derive_reference_compressed_cram(
            "samtools",
            source,
            reference,
            root,
            expected_reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        )

    digest.assert_not_called()
    assert final_crai.read_bytes() == b"stale crai"
    assert not (root / "reference-compressed/source.cram").exists()


@pytest.mark.parametrize("failure", ["full", "indexed"])
def test_lossy_evidence_aborts_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    source = _touch(tmp_path / "source.bam")
    reference = _touch(tmp_path / "reference.fa", ">chr1\nACGT\n")
    fake = _ReferenceSamtools()
    calls = 0

    def digest(*args: object, **kwargs: object) -> tuple[str, int]:
        nonlocal calls
        calls += 1
        if failure == "full" and calls == 2:
            return "decoded-full", 34_214
        if kwargs.get("region") is not None:
            if failure == "indexed" and calls == 4:
                return "decoded-region", 4_975
            return REGISTERED_B178_INDEXED_REGION_DIGEST, REGISTERED_B178_INDEXED_REGION_RECORDS
        return "full", 34_214

    monkeypatch.setattr(generator.subprocess, "run", fake)
    monkeypatch.setattr(generator, "_normalized_record_digest", digest)

    with pytest.raises(LossyConversionError, match="not lossless"):
        generator.derive_reference_compressed_cram(
            "samtools",
            source,
            reference,
            tmp_path / "cram",
            expected_reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        )

    assert not (tmp_path / "cram/reference-compressed/source.cram").exists()


def test_rerun_replaces_mutated_cram_and_crai_with_identical_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _touch(tmp_path / "source.bam", "source bam")
    reference = _touch(tmp_path / "reference.fa", ">chr1\nACGT\n")
    root = tmp_path / "cram"
    final_cram = _touch(root / "reference-compressed/source.cram", "stale cram")
    final_crai = _touch(Path(f"{final_cram}.crai"), "stale crai")
    fake = _ReferenceSamtools()
    monkeypatch.setattr(generator.subprocess, "run", fake)
    monkeypatch.setattr(generator, "_normalized_record_digest", _registered_digest)
    reference_sha256 = hashlib.sha256(reference.read_bytes()).hexdigest()

    first = generator.derive_reference_compressed_cram(
        "samtools", source, reference, root, expected_reference_sha256=reference_sha256
    )
    first_bytes = (final_cram.read_bytes(), final_crai.read_bytes())
    final_cram.write_bytes(b"mutated stale cram")
    final_crai.write_bytes(b"mutated stale crai")
    second = generator.derive_reference_compressed_cram(
        "samtools", source, reference, root, expected_reference_sha256=reference_sha256
    )

    assert second == first
    assert (final_cram.read_bytes(), final_crai.read_bytes()) == first_bytes
    assert first_bytes == (b"reference-compressed cram", b"reference-compressed crai")
