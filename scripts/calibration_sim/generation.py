"""Atomic deterministic sequence/read bundles for calibration development."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import re
import stat
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, NoReturn

from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object

from .protocol import SimulationCase, SimulationProtocol, decode_simulation_protocol, simulation_protocol_document
from .reads import generate_read_pairs

logger = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_SOURCES = {"generation.py", "haplotypes.py", "protocol.py", "reads.py"}
_CASE_FILES = ("haplotypes.fa", "origins.tsv", "reads_R1.fastq", "reads_R2.fastq", "truth.json")
_MANIFEST_FIELDS = {
    "schema_version",
    "evidence_status",
    "protocol_sha256",
    "generator",
    "pair_count",
    "generated_bases",
    "independent_group_count",
    "files",
}


@dataclass(frozen=True)
class GeneratedSimulation:
    """Installed input bundle and externally pinnable exact manifest identity."""

    output: Path
    manifest_sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _json(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + "\n"
    ).encode()


def _private_file(path: Path) -> BinaryIO:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600)
    try:
        return os.fdopen(descriptor, "wb")
    except BaseException:
        os.close(descriptor)
        raise


def _write(path: Path, data: bytes) -> None:
    with _private_file(path) as stream:
        stream.write(data)


def _truth(case: SimulationCase) -> dict[str, object]:
    return {
        "schema_version": "calibration-simulation-truth-v1",
        "case_id": case.case_id,
        "group_key": case.group_key,
        "role": case.role,
        "primary": case.primary,
        "strata": list(case.strata),
        "scenario": case.scenario,
        "backbone_family": case.backbone_family,
        "pair_family": case.pair_family,
        "seed_family": case.seed_family,
        "caller_positive": case.caller_positive,
        "truth_variant_ids": list(case.truth_variant_ids),
        "allele_repeat_counts": list(case.truth.allele_repeat_counts),
        "total_repeat_count": case.truth.total_repeat_count,
        "haplotype_sha256": list(case.truth.haplotype_sha256),
        "repeat_unit_bp": case.haplotypes[0].repeat_unit_bp,
    }


def _write_case(root: Path, case: SimulationCase, budget: int) -> None:
    folder = root / case.case_id
    folder.mkdir(mode=0o700)
    _write(folder / "truth.json", _json(_truth(case)))
    fasta = "".join(f">allele-{index + 1}\n{haplotype.sequence}\n" for index, haplotype in enumerate(case.haplotypes))
    _write(folder / "haplotypes.fa", fasta.encode("ascii"))
    with ExitStack() as stack:
        first = stack.enter_context(_private_file(folder / "reads_R1.fastq"))
        second = stack.enter_context(_private_file(folder / "reads_R2.fastq"))
        origins = stack.enter_context(_private_file(folder / "origins.tsv"))
        origins.write(b"read\tallele_index\tfragment_start\tfragment_end\tread1_substitutions\tread2_substitutions\n")
        reads = case.reads
        for pair in generate_read_pairs(
            haplotypes=case.haplotypes,
            pair_count=reads.pair_count,
            read_length=reads.read_length,
            fragment_length=reads.fragment_length,
            seed=reads.seed,
            substitution_rate=reads.substitution_rate,
            quality_score=reads.quality_score,
            allele_copy_weights=reads.allele_copy_weights,
            maximum_generated_bases=budget,
        ):
            first.write(f"@{pair.name}/1\n{pair.read1}\n+\n{pair.qualities}\n".encode("ascii"))
            second.write(f"@{pair.name}/2\n{pair.read2}\n+\n{pair.qualities}\n".encode("ascii"))
            origins.write(
                f"{pair.name}\t{pair.allele_index}\t{pair.fragment_start}\t{pair.fragment_end}\t"
                f"{pair.read1_substitutions}\t{pair.read2_substitutions}\n".encode("ascii")
            )


def _hash_child(parent: int, name: str) -> tuple[int, str]:
    descriptor = os.open(name, os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | os.O_CLOEXEC, dir_fd=parent)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            _fail("simulation inventory contains a nonregular file")
        digest = hashlib.sha256()
        observed = 0
        while True:
            data = os.read(descriptor, 1024 * 1024)
            if not data:
                break
            digest.update(data)
            observed += len(data)
        after = os.fstat(descriptor)
        if observed != before.st_size or any(
            getattr(before, key) != getattr(after, key)
            for key in ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        ):
            _fail("simulation inventory file changed while hashing")
        return observed, digest.hexdigest()
    finally:
        os.close(descriptor)


def _inventory(root: Path, expected_directories: set[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    observed_directories = set()
    if root.is_symlink() or not root.is_dir():
        _fail("simulation inventory root must be a nonsymlink directory")
    try:
        for folder, directories, files, descriptor in os.fwalk(root, follow_symlinks=False):
            for name in directories:
                relative_directory = (Path(folder) / name).relative_to(root).as_posix()
                if relative_directory not in expected_directories:
                    _fail("simulation directory inventory contains an undeclared directory")
                observed_directories.add(relative_directory)
                if not stat.S_ISDIR(os.stat(name, dir_fd=descriptor, follow_symlinks=False).st_mode):
                    _fail("simulation inventory contains a symlink directory")
            for name in files:
                relative = (Path(folder) / name).relative_to(root).as_posix()
                size, digest = _hash_child(descriptor, name)
                rows.append({"path": relative, "size": size, "sha256": digest})
    except OSError as error:
        raise ValueError("simulation inventory is unreadable, changed, or contains a symlink") from error
    if observed_directories != expected_directories:
        _fail("simulation directory inventory is incomplete")
    return sorted(rows, key=lambda row: str(row["path"]))


def _generator_identity() -> dict[str, object]:
    root = Path(__file__).parent
    return {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "source_sha256": {
            name: hashlib.sha256(read_regular_path(root / name)).hexdigest() for name in sorted(_SOURCES)
        },
    }


def _validate_generator(value: object) -> None:
    if not isinstance(value, dict) or set(value) != {"python_implementation", "python_version", "source_sha256"}:
        _fail("simulation generator identity fields differ from the closed contract")
    if any(
        not isinstance(value[key], str) or not value[key].strip() for key in ("python_implementation", "python_version")
    ):
        _fail("simulation generator runtime identity must contain nonempty text")
    sources = value["source_sha256"]
    if not isinstance(sources, dict) or set(sources) != _SOURCES:
        _fail("simulation generator source inventory differs")
    if any(not isinstance(digest, str) or _DIGEST.fullmatch(digest) is None for digest in sources.values()):
        _fail("simulation generator source identities must be lowercase SHA256 digests")


def generate_simulation_bundle(protocol: SimulationProtocol, output: Path) -> GeneratedSimulation:
    """Write finite paired reads, source truth, and hashed files atomically.

    Args:
        protocol: Validated design with frozen role/family/resource declarations.
        output: New output directory; existing paths are always refused.

    Returns:
        Exact manifest digest for independent verification. It identifies only
        generated inputs, not caller performance, scientific power or custody.

    Raises:
        ValueError: For altered protocols, occupied outputs or inconsistent files.
        RuntimeError: If atomic publication is unavailable or generation fails.
        OSError: If private output files cannot be created or written.
    """
    document = simulation_protocol_document(protocol)
    manifest_digest = ""

    def produce(staging: Path) -> bool:
        nonlocal manifest_digest
        _write(staging / "protocol.json", _json(document))
        for case in protocol.cases:
            _write_case(staging, case, protocol.maximum_generated_bases)
        manifest = {
            "schema_version": "calibration-simulation-inputs-v1",
            "evidence_status": "generated-inputs-only",
            "protocol_sha256": protocol.sha256,
            "generator": _generator_identity(),
            "pair_count": protocol.pair_count,
            "generated_bases": protocol.generated_bases,
            "independent_group_count": protocol.independent_group_count,
            "files": _inventory(staging, {case.case_id for case in protocol.cases}),
        }
        encoded = _json(manifest)
        manifest_digest = hashlib.sha256(encoded).hexdigest()
        _write(staging / "manifest.json", encoded)
        verify_generated_bundle(staging, expected_manifest_sha256=manifest_digest)
        return True

    atomic_output(output, produce)
    return GeneratedSimulation(output, manifest_digest)


def verify_generated_bundle(root: Path, *, expected_manifest_sha256: str) -> SimulationProtocol:
    """Verify a complete input inventory against an externally pinned manifest.

    Args:
        root: Previously generated input directory.
        expected_manifest_sha256: Exact manifest digest recorded outside the bundle.

    Returns:
        Revalidated declared simulation protocol. This verifies bytes and layout;
        an independent outcome oracle must separately test algorithm correctness.

    Raises:
        ValueError: On changed bindings, malformed declarations or incomplete files.
    """
    if not isinstance(expected_manifest_sha256, str) or _DIGEST.fullmatch(expected_manifest_sha256) is None:
        _fail("simulation expected manifest identity must be a lowercase SHA256")
    if not isinstance(root, Path) or root.is_symlink() or not root.is_dir():
        _fail("simulation inventory root must be a nonsymlink directory")
    encoded = read_regular_path(root / "manifest.json")
    if hashlib.sha256(encoded).hexdigest() != expected_manifest_sha256:
        _fail("simulation manifest differs from its externally recorded identity")
    try:
        manifest = load_strict_json_object(encoded)
        protocol = decode_simulation_protocol(load_strict_json_object(read_regular_path(root / "protocol.json")))
    except (ValueError, UnicodeError) as error:
        raise ValueError("simulation manifest or protocol is malformed") from error
    if not isinstance(manifest, dict) or set(manifest) != _MANIFEST_FIELDS:
        _fail("simulation manifest fields differ from the closed contract")
    _validate_generator(manifest["generator"])
    for field in ("pair_count", "generated_bases", "independent_group_count"):
        if type(manifest[field]) is not int or manifest[field] < 0:
            _fail("simulation manifest resource counts must be nonnegative integers")
    if (
        manifest["schema_version"] != "calibration-simulation-inputs-v1"
        or manifest["evidence_status"] != "generated-inputs-only"
        or manifest["protocol_sha256"] != protocol.sha256
        or manifest["pair_count"] != protocol.pair_count
        or manifest["generated_bases"] != protocol.generated_bases
        or manifest["independent_group_count"] != protocol.independent_group_count
    ):
        _fail("simulation manifest differs from its declared protocol")
    rows = _inventory(root, {case.case_id for case in protocol.cases})
    expected_names = {"protocol.json"} | {f"{case.case_id}/{name}" for case in protocol.cases for name in _CASE_FILES}
    payload = [row for row in rows if row["path"] != "manifest.json"]
    if {row["path"] for row in payload} != expected_names or payload != manifest["files"]:
        _fail("simulation inventory differs from its exact file manifest")
    own = [row for row in rows if row["path"] == "manifest.json"]
    if len(own) != 1 or own[0]["sha256"] != expected_manifest_sha256:
        _fail("simulation manifest changed during inventory verification")
    return protocol
