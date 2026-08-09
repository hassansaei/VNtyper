"""Private CRAM M5 cache ownership and completeness tests."""

from __future__ import annotations

import gzip
import hashlib
import os
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from vntyper.scripts.reference_cache_binding import (
    PrivateReferenceCache,
    cache_entry_path,
    probe_private_reference_cache,
)
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution, restore_reference_resolution
from vntyper.scripts.reference_uri_policy import LocalHeaderReference

pytestmark = pytest.mark.unit


def _digest(sequence: str) -> str:
    return hashlib.md5(sequence.encode()).hexdigest()


def _cache_entry(root: Path, digest: str) -> Path:
    return root / digest[:2] / digest[2:4] / digest[4:]


def test_compressed_header_reference_closes_its_underlying_descriptor(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa.gz"
    with gzip.open(reference, "wb") as handle:
        handle.write(b">chr1\nACGT\n")
    cache = object.__new__(PrivateReferenceCache)
    before = set(os.listdir("/proc/self/fd"))

    with cache._open_fasta(str(reference)) as handle:
        assert handle.read() == b">chr1\nACGT\n"

    assert set(os.listdir("/proc/self/fd")) == before


def test_cache_entry_path_consumes_digest_progressively() -> None:
    digest = "0123456789abcdef0123456789abcdef"

    assert cache_entry_path("/cache/%2s/%2s/%s", digest) == Path("/cache/01/23/456789abcdef0123456789abcdef")


@pytest.mark.parametrize("pattern", ["/cache/no-token", "/cache/%40s", "/cache/%2s"])
def test_cache_entry_path_rejects_templates_that_do_not_consume_the_digest(pattern: str) -> None:
    with pytest.raises(ValueError, match="Invalid local CRAM reference cache template"):
        cache_entry_path(pattern, "0123456789abcdef0123456789abcdef")


def test_second_local_ref_path_entry_can_supply_a_verified_read_only_cache(tmp_path: Path) -> None:
    sequence = "ACGT" * 25
    digest = _digest(sequence)
    source_root = tmp_path / "operator-cache"
    source = _cache_entry(source_root, digest)
    source.parent.mkdir(parents=True)
    source.write_text(sequence, encoding="ascii")
    before = source.read_bytes()
    source.chmod(0o444)
    source.parent.chmod(0o555)
    (tmp_path / "run").mkdir()
    binding: PrivateReferenceCache | None = None

    try:
        binding = PrivateReferenceCache(
            (("chr1", digest),),
            f"invalid-template:{source_root}/%2s/%2s/%s",
            (),
            tmp_path / "run",
            "input",
        )
        private_entry = cache_entry_path(binding.pattern, digest)
        assert private_entry.read_bytes() == before
        assert source.read_bytes() == before
    finally:
        if binding is not None:
            binding.close()
        source.parent.chmod(0o755)
        source.chmod(0o644)

    assert not (tmp_path / "run" / ".input_reference_cache").exists()


def test_local_cache_entry_precedes_remote_ambient_fallback_in_the_retained_environment(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sequence = "ACGT" * 25
    digest = _digest(sequence)
    source_root = tmp_path / "operator-cache"
    source = _cache_entry(source_root, digest)
    source.parent.mkdir(parents=True)
    source.write_text(sequence, encoding="ascii")
    ambient = f"{source_root}/%2s/%2s/%s:https://refget.example/%s"
    monkeypatch.setenv("REF_PATH", ambient)
    observations: list[tuple[str, str]] = []

    previous = pin_reference_resolution({"cram": {"allow_ambient_reference_resolution": True}})
    binding: PrivateReferenceCache | None = None

    def observe_environment(_path: str | None, _position: int) -> tuple[bool, str]:
        observations.append((os.environ["REF_CACHE"], os.environ["REF_PATH"]))
        return True, ""

    try:
        binding, reason = probe_private_reference_cache(
            header_contigs=("chr1",),
            header_m5s=(("chr1", digest),),
            header_references=(LocalHeaderReference("chr1", digest, str(tmp_path / "stale.fa")),),
            config={"cram": {"allow_ambient_reference_resolution": True}},
            has_remote_header_reference=False,
            output_dir=str(tmp_path),
            output_name="input",
            position=4,
            probe=observe_environment,
        )

        assert reason is None
        assert binding is not None
        assert observations == [(binding.pattern, "https://refget.example/%s")]
    finally:
        if binding is not None:
            binding.close()
        restore_reference_resolution(previous)

    assert source.read_text(encoding="ascii") == sequence
    assert not (tmp_path / ".input_reference_cache").exists()


def test_default_mode_binds_all_header_cache_entries_against_later_path_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sequences = {"chr1": "AAAA", "chr2": "CCCC"}
    source_root = tmp_path / "operator-cache"
    for sequence in sequences.values():
        source = _cache_entry(source_root, _digest(sequence))
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(sequence, encoding="ascii")
    configured = f"{source_root}/%2s/%2s/%s"
    monkeypatch.setenv("REF_PATH", configured)

    previous = pin_reference_resolution({"cram": {"local_ref_path": configured}})
    binding: PrivateReferenceCache | None = None

    def replace_second_source(_path: str | None, _position: int) -> tuple[bool, str]:
        private_pattern = os.environ["REF_CACHE"]
        assert os.environ["REF_PATH"] == private_pattern
        second = _cache_entry(source_root, _digest(sequences["chr2"]))
        replacement = tmp_path / "replacement"
        replacement.write_text("GGGG", encoding="ascii")
        os.replace(replacement, second)
        return True, ""

    try:
        binding, reason = probe_private_reference_cache(
            header_contigs=("chr1", "chr2"),
            header_m5s=tuple((contig, _digest(sequence)) for contig, sequence in sequences.items()),
            header_references=(LocalHeaderReference("chr1", _digest("AAAA"), str(tmp_path / "chr1.fa")),),
            config={"cram": {"local_ref_path": configured}},
            has_remote_header_reference=False,
            output_dir=str(tmp_path),
            output_name="input",
            position=4,
            probe=replace_second_source,
        )

        assert reason is None
        assert binding is not None
        assert cache_entry_path(binding.pattern, _digest("CCCC")).read_text(encoding="ascii") == "CCCC"
    finally:
        if binding is not None:
            binding.close()
        restore_reference_resolution(previous)

    assert not (tmp_path / ".input_reference_cache").exists()


def test_ambient_remote_header_contig_does_not_require_a_local_cache_entry(tmp_path: Path) -> None:
    sequence = "AAAA"
    digest = _digest(sequence)
    local_reference = tmp_path / "chr1.fa"
    local_reference.write_text(f">chr1\n{sequence}\n", encoding="ascii")

    with patch("vntyper.scripts.reference_cache_binding.activate_private_reference_cache"):
        binding, reason = probe_private_reference_cache(
            header_contigs=("chr1", "chr2"),
            header_m5s=(("chr1", digest), ("chr2", _digest("CCCC"))),
            header_references=(LocalHeaderReference("chr1", digest, str(local_reference)),),
            config={"cram": {"allow_ambient_reference_resolution": True}},
            has_remote_header_reference=True,
            output_dir=str(tmp_path),
            output_name="input",
            position=4,
            probe=Mock(return_value=(True, "")),
        )

    assert reason is None
    assert binding is not None
    binding.close()


def test_digest_read_failure_cleans_the_immediately_tracked_entry(tmp_path: Path) -> None:
    sequence = "ACGT"
    digest = _digest(sequence)
    source_root = tmp_path / "operator-cache"
    source = _cache_entry(source_root, digest)
    source.parent.mkdir(parents=True)
    source.write_text(sequence, encoding="ascii")
    (tmp_path / "run").mkdir()

    with (
        patch("vntyper.scripts.reference_cache_binding._digest_file", side_effect=OSError("read failed")),
        pytest.raises(OSError, match="read failed"),
    ):
        PrivateReferenceCache(
            (("chr1", digest),),
            f"{source_root}/%2s/%2s/%s",
            (),
            tmp_path / "run",
            "input",
        )

    assert not (tmp_path / "run" / ".input_reference_cache").exists()


def test_digest_mismatch_is_unresolved_and_leaves_no_private_namespace(tmp_path: Path) -> None:
    required_digest = _digest("AAAA")
    source_root = tmp_path / "operator-cache"
    source = _cache_entry(source_root, required_digest)
    source.parent.mkdir(parents=True)
    source.write_text("CCCC", encoding="ascii")
    (tmp_path / "run").mkdir()

    with pytest.raises(ValueError, match="could not resolve contig=chr1"):
        PrivateReferenceCache(
            (("chr1", required_digest),),
            f"{source_root}/%2s/%2s/%s",
            (),
            tmp_path / "run",
            "input",
        )

    assert source.read_text(encoding="ascii") == "CCCC"
    assert not (tmp_path / "run" / ".input_reference_cache").exists()


def test_separate_local_fasta_files_populate_every_required_header_digest(tmp_path: Path) -> None:
    sequences = {"chr1": "ACGT" * 25, "chr2": "TGCA" * 25}
    references: list[LocalHeaderReference] = []
    before: dict[Path, bytes] = {}
    for contig, sequence in sequences.items():
        path = tmp_path / "operator" / f"{contig}.fa"
        path.parent.mkdir(exist_ok=True)
        path.write_text(f">{contig}\n{sequence}\n", encoding="ascii")
        before[path] = path.read_bytes()
        references.append(LocalHeaderReference(contig, _digest(sequence), str(path)))
    (tmp_path / "run").mkdir()

    binding = PrivateReferenceCache(
        tuple((contig, _digest(sequence)) for contig, sequence in sequences.items()),
        str(tmp_path / "empty-cache" / "%2s" / "%2s" / "%s"),
        references,
        tmp_path / "run",
        "input",
    )
    try:
        for sequence in sequences.values():
            assert cache_entry_path(binding.pattern, _digest(sequence)).read_text(encoding="ascii") == sequence
        assert {path: path.read_bytes() for path in before} == before
    finally:
        binding.close()

    assert not (tmp_path / "run" / ".input_reference_cache").exists()


def test_private_probe_requires_m5_for_every_header_contig_without_running_a_command(tmp_path: Path) -> None:
    probe = Mock()

    binding, reason = probe_private_reference_cache(
        header_contigs=("chr1", "chr2"),
        header_m5s=(("chr1", _digest("AAAA")),),
        header_references=(LocalHeaderReference("chr2", None, str(tmp_path / "chr2.fa")),),
        config={"cram": {"local_ref_path": str(tmp_path / "%2s" / "%2s" / "%s")}},
        has_remote_header_reference=False,
        output_dir=str(tmp_path),
        output_name="input",
        position=4,
        probe=probe,
    )

    assert binding is None
    assert reason == "header M5 is missing for contig chr2"
    probe.assert_not_called()


@pytest.mark.parametrize("probe_result", [(True, ""), (False, "decode failed")], ids=["success", "failure"])
def test_private_probe_activation_retains_only_a_successful_complete_cache(
    tmp_path: Path,
    probe_result: tuple[bool, str],
) -> None:
    sequence = "AAAA"
    digest = _digest(sequence)
    source_root = tmp_path / "operator-cache"
    source = _cache_entry(source_root, digest)
    source.parent.mkdir(parents=True)
    source.write_text(sequence, encoding="ascii")
    probe = Mock(return_value=probe_result)

    with patch("vntyper.scripts.reference_cache_binding.activate_private_reference_cache") as activate:
        binding, reason = probe_private_reference_cache(
            header_contigs=("chr1",),
            header_m5s=(("chr1", digest),),
            header_references=(LocalHeaderReference("chr1", digest, str(tmp_path / "missing.fa")),),
            config={"cram": {"local_ref_path": f"{source_root}/%2s/%2s/%s"}},
            has_remote_header_reference=False,
            output_dir=str(tmp_path),
            output_name="input",
            position=4,
            probe=probe,
        )

    probe.assert_called_once_with(None, 4)
    expected_pattern = str(tmp_path / ".input_reference_cache" / "%2s" / "%2s" / "%s")
    activate.assert_called_once_with(expected_pattern, remote_reference_path=None)
    if probe_result[0]:
        assert binding is not None
        assert reason is None
        binding.close()
    else:
        assert binding is None
        assert reason == "decode failed"
    assert not (tmp_path / ".input_reference_cache").exists()


def test_reserved_private_cache_namespace_collision_is_fatal(tmp_path: Path) -> None:
    reserved = tmp_path / ".input_reference_cache"
    reserved.mkdir()

    with pytest.raises(RuntimeError, match="reserved CRAM reference cache directory already exists"):
        PrivateReferenceCache((), str(tmp_path / "%s"), (), tmp_path, "input")

    assert reserved.is_dir()
