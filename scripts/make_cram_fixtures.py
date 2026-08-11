#!/usr/bin/env python3
"""Derive lossless CRAM fixtures from the BAM test cohort.

Most fixtures use ``no_ref=1`` because cohort BAM headers lack M5 tags. The registered
b178 fixture exercises explicit reference compression and a pinned indexed-region oracle;
purpose-built fixtures cover unavailable-reference and placed-unmapped failure modes.
Outputs live under ignored ``tests/data`` because the cohort ships separately and
BAM/CRAM equivalence is stronger evidence than a committed standalone CRAM.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import logging
import subprocess
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import pysam

if __package__:
    from .cram_reference_contract import (
        HG19_CHR1_REFERENCE_SHA256,
        REFERENCE_VALIDATION_REGION,
        LossyConversionError,
        ReferenceCompressedFixture,
        ValidatedReferenceSnapshot,
        header_with_hg19_m5,
        normalize_sam_record,
        snapshot_reference_fasta,
        validate_registered_b178_index_evidence,
    )
else:
    from cram_reference_contract import (  # type: ignore[no-redef]
        HG19_CHR1_REFERENCE_SHA256,
        REFERENCE_VALIDATION_REGION,
        LossyConversionError,
        ReferenceCompressedFixture,
        ValidatedReferenceSnapshot,
        header_with_hg19_m5,
        normalize_sam_record,
        snapshot_reference_fasta,
        validate_registered_b178_index_evidence,
    )

logger = logging.getLogger("cram_fixtures")

#: Default encoding is reference-independent; see the module docstring.
CRAM_WRITE_OPTIONS = ("no_ref=1",)

DEFAULT_FIXTURE_ROOT = Path("tests/data/cram")
DEFAULT_DATA_CONFIG = Path("tests/test_data_config.json")
DEFAULT_REFERENCE_COMPRESSED_SOURCE = Path("tests/data/example_b178_hg19_subset.bam")
DEFAULT_REFERENCE_COMPRESSED_FASTA = Path("reference/alignment/chr1.hg19.fa")

pysam_any: Any = pysam


@dataclass
class Fixture:
    """One derived CRAM and the evidence that it is faithful to its source."""

    source_bam: Path
    cram: Path
    records: int
    unmapped_reads: int
    #: Digest of the decoded record stream, identical for the BAM and the CRAM.
    record_digest: str
    source_bytes: int
    cram_bytes: int

    def as_manifest_entry(self) -> dict[str, object]:
        return {
            "source_bam": str(self.source_bam),
            "cram": str(self.cram),
            "records": self.records,
            "unmapped_reads": self.unmapped_reads,
            "record_digest": self.record_digest,
            "source_bytes": self.source_bytes,
            "cram_bytes": self.cram_bytes,
        }


@dataclass
class Summary:
    """What a run produced, for the caller and for the attestation."""

    fixtures: list[Fixture] = field(default_factory=list)
    skipped: list[tuple[Path, str]] = field(default_factory=list)

    @property
    def total_source_bytes(self) -> int:
        return sum(f.source_bytes for f in self.fixtures)

    @property
    def total_cram_bytes(self) -> int:
        return sum(f.cram_bytes for f in self.fixtures)


@dataclass(frozen=True)
class ReferenceDependentFixture:
    """A reference-compressed CRAM and the local reference named by its header."""

    cram: Path
    reference: Path


@dataclass(frozen=True)
class PlacedFlag12Fixture:
    """A CRAM whose placed flag-12 records make indexed extraction lossy."""

    cram: Path


@dataclass(frozen=True)
class IndexedSafeFixture:
    """A CRAM with nonempty unplaced reads and no placed-unmapped records."""

    cram: Path


def _run(argv: list[str]) -> str:
    """Run a samtools command, returning stdout and raising with stderr on failure."""
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{' '.join(argv)} exited {completed.returncode}: {completed.stderr.strip()}")
    return completed.stdout


def _run_to_file(argv: list[str], output: Path) -> None:
    """Run a command whose binary stdout is the requested output file."""
    with output.open("wb") as handle:
        completed = subprocess.run(argv, stdout=handle, stderr=subprocess.PIPE, check=False)
    if completed.returncode != 0:
        stderr = completed.stderr.decode() if isinstance(completed.stderr, bytes) else completed.stderr
        raise RuntimeError(f"{' '.join(argv)} exited {completed.returncode}: {stderr.strip()}")


def _record_digest(samtools: str, alignment: Path) -> tuple[str, int]:
    """Digest every decoded record of an alignment, with no reference supplied.

    Reading the CRAM with no ``-T`` is the point: it is what the pipeline does. If the
    fixture needed a reference this call would fail rather than quietly succeed.
    """
    return _stream_record_digest([samtools, "view", str(alignment)], lambda line: line.encode())


def _stream_record_digest(argv: list[str], normalize: Callable[[str], bytes]) -> tuple[str, int]:
    """Digest command stdout while draining arbitrary stderr without a pipe deadlock."""
    digest = hashlib.sha256()
    count = 0
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as stderr_file:
        with subprocess.Popen(
            argv,
            stdout=subprocess.PIPE,
            stderr=stderr_file,
            text=True,
            bufsize=1024 * 1024,
        ) as proc:
            assert proc.stdout is not None
            for line in proc.stdout:
                digest.update(normalize(line))
                count += 1
            proc.wait()
        stderr_file.seek(0)
        stderr = stderr_file.read()
    if proc.returncode != 0:
        raise RuntimeError(f"{' '.join(argv)} exited {proc.returncode}: {stderr.strip()}")
    return digest.hexdigest(), count


def _normalized_record_digest(
    samtools: str,
    alignment: Path,
    explicit_reference: Path | None = None,
    *,
    region: str | None = None,
) -> tuple[str, int]:
    """Digest decoded SAM records while ignoring optional-tag ordering.

    Args:
        samtools: The samtools executable.
        alignment: BAM or CRAM whose decoded records are digested.
        explicit_reference: FASTA passed with ``-T`` for CRAM decoding.
        region: Optional indexed region appended to the samtools query.
    """
    argv = [samtools, "view"]
    if explicit_reference is not None:
        argv.extend(["-T", str(explicit_reference)])
    argv.append(str(alignment))
    if region is not None:
        argv.append(region)

    return _stream_record_digest(argv, normalize_sam_record)


def discover_source_bams(data_root: Path, fixture_root: Path) -> list[Path]:
    """Every BAM in the cohort, excluding anything already under the fixture root."""
    bams = sorted(p for p in data_root.rglob("*.bam") if fixture_root not in p.parents)
    logger.info("discovered %d source BAMs under %s", len(bams), data_root)
    return bams


def _select_source_bams(discovered: list[Path], *, data_config: Path, data_root: Path, include_all: bool) -> list[Path]:
    """Import the selection policy in package and direct-script execution modes."""
    if __package__:
        from .cram_fixture_selection import select_source_bams as package_selector

        return cast(
            list[Path],
            package_selector(discovered, data_config=data_config, data_root=data_root, include_all=include_all),
        )
    from cram_fixture_selection import select_source_bams as direct_selector

    return cast(
        list[Path],
        direct_selector(discovered, data_config=data_config, data_root=data_root, include_all=include_all),
    )


def _derive_declared_single_end_fixtures(data_config: Path, repository_root: Path) -> None:
    """Dispatch every validated Task9 fixture declaration.

    Args:
        data_config: Manifest containing top-level ``derived_fixtures`` entries.
        repository_root: Root against which their portable paths are resolved.
    """
    payload = json.loads(data_config.read_text())
    declarations = payload.get("derived_fixtures", [])
    module_name = f"{__package__}.single_end_fixture" if __package__ else "single_end_fixture"
    fixture_builder = importlib.import_module(module_name)

    for declaration in declarations:
        spec = fixture_builder.parse_single_end_fixture(declaration, root=repository_root)
        fixture_builder.derive_single_end_bam(spec)


def derive_cram(samtools: str, bam: Path, data_root: Path, fixture_root: Path) -> Fixture:
    """Convert one BAM to CRAM and verify the conversion lost nothing.

    Raises:
        LossyConversionError: if the CRAM's decoded records differ from the BAM's.
    """
    relative = bam.relative_to(data_root)
    cram = (fixture_root / relative).with_suffix(".cram")
    cram.parent.mkdir(parents=True, exist_ok=True)

    options: list[str] = []
    for option in CRAM_WRITE_OPTIONS:
        options += ["--output-fmt-option", option]
    _run([samtools, "view", "-C", *options, "-o", str(cram), str(bam)])
    _run([samtools, "index", str(cram)])

    bam_digest, bam_records = _record_digest(samtools, bam)
    cram_digest, cram_records = _record_digest(samtools, cram)
    if bam_digest != cram_digest:
        raise LossyConversionError(
            f"{bam} -> {cram} is not lossless: {bam_records} source records digest "
            f"{bam_digest[:16]}, {cram_records} decoded records digest {cram_digest[:16]}"
        )

    # The pipeline recovers every flag-4 read, including single-end and placed-unmapped
    # records, so the fixture evidence must use the same inclusive predicate.
    unmapped = int(_run([samtools, "view", "-c", "-f", "4", str(cram)]).strip())

    logger.info(
        "%s -> %s  (%d records, %d unmapped reads, %.1f%% of BAM size)",
        relative,
        cram.relative_to(fixture_root),
        bam_records,
        unmapped,
        100.0 * cram.stat().st_size / bam.stat().st_size,
    )
    return Fixture(
        source_bam=bam,
        cram=cram,
        records=bam_records,
        unmapped_reads=unmapped,
        record_digest=bam_digest,
        source_bytes=bam.stat().st_size,
        cram_bytes=cram.stat().st_size,
    )


def derive_reference_compressed_cram(
    samtools: str,
    source_bam: Path,
    reference: Path,
    fixture_root: Path,
    *,
    expected_reference_sha256: str = HG19_CHR1_REFERENCE_SHA256,
    validated_reference: ValidatedReferenceSnapshot | None = None,
) -> ReferenceCompressedFixture:
    """Derive and prove a real-read CRAM decoded with an explicit FASTA.

    Args:
        samtools: The samtools executable.
        source_bam: Registered b178 BAM used as the record source.
        reference: Exact chr1 hg19 FASTA used for encoding and decoding.
        fixture_root: Root beneath which the fixture is regenerated.
        expected_reference_sha256: Required identity of the reference bytes.
        validated_reference: Optional immutable snapshot from an earlier validation.

    Returns:
        The derived fixture and stable digest evidence.

    Raises:
        FileNotFoundError: If the explicit reference does not exist.
        ValueError: If the supplied reference is not the pinned chr1 hg19 FASTA.
        LossyConversionError: If decoded records or record counts differ.
        RuntimeError: If encoding, indexing, or decoding fails.
    """
    if validated_reference is None:
        with snapshot_reference_fasta(reference, expected_reference_sha256) as snapshot:
            return derive_reference_compressed_cram(
                samtools,
                source_bam,
                reference,
                fixture_root,
                expected_reference_sha256=expected_reference_sha256,
                validated_reference=snapshot,
            )
    reference_for_tools = validated_reference.tool_path_for(reference, expected_reference_sha256)

    cram = fixture_root / "reference-compressed" / f"{source_bam.stem}.cram"
    crai = Path(f"{cram}.crai")
    cram.parent.mkdir(parents=True, exist_ok=True)
    source_header = _run([samtools, "view", "-H", "--no-PG", str(source_bam)])
    prepared_header = header_with_hg19_m5(source_header, validated_reference.source_uri)
    with tempfile.TemporaryDirectory(prefix="reference-compressed-", dir=cram.parent) as temp_dir:
        temp_root = Path(temp_dir)
        header_path = temp_root / "header.sam"
        prepared_bam = temp_root / source_bam.name
        encoded_cram = temp_root / f"encoded-{cram.name}"
        staged_cram = temp_root / cram.name
        staged_crai = Path(f"{staged_cram}.crai")
        header_path.write_text(prepared_header, encoding="utf-8")
        _run_to_file([samtools, "reheader", "-P", str(header_path), str(source_bam)], prepared_bam)
        _run(
            [
                samtools,
                "view",
                "-C",
                "-T",
                str(reference_for_tools),
                "--no-PG",
                "-o",
                str(encoded_cram),
                str(prepared_bam),
            ]
        )
        _run_to_file([samtools, "reheader", "-P", str(header_path), str(encoded_cram)], staged_cram)
        _run([samtools, "index", str(staged_cram)])
        if not staged_crai.is_file():
            raise RuntimeError(f"samtools index did not create expected CRAI: {staged_crai}")

        source_digest, source_records = _normalized_record_digest(samtools, source_bam)
        decoded_digest, decoded_records = _normalized_record_digest(samtools, staged_cram, reference_for_tools)
        if source_digest != decoded_digest or source_records != decoded_records:
            raise LossyConversionError(
                f"{source_bam} -> {cram} reference-compressed CRAM is not lossless: "
                f"{source_records} source records digest {source_digest[:16]}, "
                f"{decoded_records} decoded records digest {decoded_digest[:16]}"
            )

        source_region_digest, source_region_records = _normalized_record_digest(
            samtools, source_bam, region=REFERENCE_VALIDATION_REGION
        )
        decoded_region_digest, decoded_region_records = _normalized_record_digest(
            samtools,
            staged_cram,
            reference_for_tools,
            region=REFERENCE_VALIDATION_REGION,
        )
        validate_registered_b178_index_evidence(
            source_region_digest,
            source_region_records,
            decoded_region_digest,
            decoded_region_records,
        )
        staged_crai.replace(crai)
        staged_cram.replace(cram)

    return ReferenceCompressedFixture(
        source_bam=source_bam,
        cram=cram,
        reference=reference,
        records=source_records,
        source_record_digest=source_digest,
        decoded_record_digest=decoded_digest,
        indexed_region=REFERENCE_VALIDATION_REGION,
        indexed_region_records=source_region_records,
        source_indexed_region_digest=source_region_digest,
        decoded_indexed_region_digest=decoded_region_digest,
        reference_sha256=validated_reference.sha256,
        source_bytes=source_bam.stat().st_size,
        cram_bytes=cram.stat().st_size,
    )


def build_fixtures(
    samtools: str,
    data_root: Path,
    fixture_root: Path,
    limit: int | None = None,
    *,
    data_config: Path = DEFAULT_DATA_CONFIG,
    include_all: bool = False,
    repository_root: Path = Path("."),
) -> Summary:
    """Derive verified CRAMs for declared BAMs, or every BAM when requested.

    Args:
        samtools: The samtools executable.
        data_root: The source ``tests/data`` root.
        fixture_root: Destination root for derived CRAMs.
        limit: Optional smoke-test cap after selection.
        data_config: Manifest defining the normal fixture set.
        include_all: Derive every discovered source BAM for the golden cohort.
        repository_root: Root used to resolve portable derived-fixture paths.

    Returns:
        The derivation summary.
    """
    summary = Summary()
    bams = _select_source_bams(
        discover_source_bams(data_root, fixture_root),
        data_config=data_config,
        data_root=data_root,
        include_all=include_all,
    )
    _derive_declared_single_end_fixtures(data_config, repository_root)
    if limit is not None:
        bams = bams[:limit]
    for bam in bams:
        try:
            summary.fixtures.append(derive_cram(samtools, bam, data_root, fixture_root))
        except LossyConversionError:
            # A lossy conversion is never acceptable: the whole point of the fixture is
            # that BAM and CRAM carry the same reads. Fail the run rather than record it.
            raise
        except RuntimeError as exc:
            logger.error("skipping %s: %s", bam, exc)
            summary.skipped.append((bam, str(exc)))
    return summary


def _write_reference(path: Path, *, length: int) -> str:
    """Write and index a deterministic single-contig FASTA, returning its sequence."""
    sequence = "A" * length
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">chr1\n{sequence}\n")
    pysam_any.faidx(str(path))
    return sequence


def _segment(name: str, *, flag: int, reference_id: int, reference_start: int) -> pysam.AlignedSegment:
    """Build one short, deterministic alignment record for a purpose-built fixture."""
    record = pysam.AlignedSegment()
    record.query_name = name
    record.flag = flag
    record.reference_id = reference_id
    record.reference_start = reference_start
    record.mapping_quality = 60
    record.cigarstring = "100M" if reference_id >= 0 else None
    record.query_sequence = "A" * 100
    record.query_qualities = pysam.qualitystring_to_array("I" * 100)
    return record


def build_reference_dependent_fixture(fixture_root: Path) -> ReferenceDependentFixture:
    """Create a CRAM that needs the copied local FASTA named by its ``UR:`` header.

    The 10 kb/50-read fixture deliberately stores a reference-compressed payload. Tests
    exercising a missing reference must rename or remove ``reference`` first: htslib may
    otherwise satisfy a no-``-T`` decode from the header's still-resolvable ``UR:`` path.

    Args:
        fixture_root: Destination directory for the purpose-built fixture.
    Returns:
        The CRAM and the copied reference which its header names.
    """
    root = fixture_root / "reference-dependent"
    reference = root / "reference.fa"
    sequence = _write_reference(reference, length=10_000)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": len(sequence), "M5": hashlib.md5(sequence.encode()).hexdigest(), "UR": str(reference)}
        ],
    }
    cram = root / "reference-dependent.cram"
    with pysam.AlignmentFile(str(cram), "wc", header=header, reference_filename=str(reference)) as output:
        for number in range(50):
            output.write(_segment(f"reference-{number}", flag=0, reference_id=0, reference_start=100 + number * 90))
    pysam_any.index(str(cram))
    return ReferenceDependentFixture(cram=cram, reference=reference)


def build_placed_flag12_fixture(fixture_root: Path) -> PlacedFlag12Fixture:
    """Create the 130-read fixture proving ``idxstats`` column four detects loss.

    It holds 600 mapped records, 25 placed flag-12 pairs and 40 unplaced flag-12 pairs.
    ``samtools view -f 12 <cram> '*'`` therefore sees 80 records while a full scan sees
    130, and ``idxstats`` reports the missing 50 in its fourth ``chr1`` column.

    Args:
        fixture_root: Destination directory for the purpose-built fixture.
    Returns:
        The indexed CRAM fixture.
    """
    root = fixture_root / "placed-flag12"
    root.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 20_000}]}
    cram = root / "placed-flag12.cram"
    with pysam_any.AlignmentFile(
        str(cram), "wc", header=header, format_options=[option.encode() for option in CRAM_WRITE_OPTIONS]
    ) as output:
        for number in range(600):
            output.write(_segment(f"mapped-{number}", flag=0, reference_id=0, reference_start=number * 20))
        for number in range(50):
            output.write(
                _segment(f"placed-{number // 2}", flag=12, reference_id=0, reference_start=13_000 + number * 10)
            )
        for number in range(80):
            output.write(_segment(f"unplaced-{number // 2}", flag=12, reference_id=-1, reference_start=-1))
    pysam_any.index(str(cram))
    return PlacedFlag12Fixture(cram=cram)


def build_indexed_safe_fixture(fixture_root: Path) -> IndexedSafeFixture:
    """Create a nonempty CRAM on which indexed and stream extraction are equivalent.

    Args:
        fixture_root: Destination directory for the purpose-built fixture.

    Returns:
        The indexed CRAM fixture.
    """
    root = fixture_root / "indexed-safe"
    root.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 20_000}]}
    cram = root / "indexed-safe.cram"
    with pysam_any.AlignmentFile(
        str(cram), "wc", header=header, format_options=[option.encode() for option in CRAM_WRITE_OPTIONS]
    ) as output:
        for number in range(10):
            output.write(_segment(f"mapped-{number}", flag=65, reference_id=0, reference_start=number * 200))
            output.write(_segment(f"mapped-{number}", flag=129, reference_id=0, reference_start=number * 200 + 100))
        for number in range(10):
            output.write(_segment(f"unplaced-{number}", flag=77, reference_id=-1, reference_start=-1))
            output.write(_segment(f"unplaced-{number}", flag=141, reference_id=-1, reference_start=-1))
    pysam_any.index(str(cram))
    return IndexedSafeFixture(cram=cram)


def write_manifest(
    summary: Summary,
    manifest_path: Path,
    *,
    reference_compressed: ReferenceCompressedFixture | None = None,
) -> None:
    """Record what was derived, so a gate run can cite it rather than re-derive it."""
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, object] = {
        "encoding": list(CRAM_WRITE_OPTIONS),
        "verified": "every fixture's decoded record stream digests identically to its source BAM",
        "fixtures": [f.as_manifest_entry() for f in summary.fixtures],
        "skipped": [{"bam": str(p), "reason": r} for p, r in summary.skipped],
    }
    if reference_compressed is not None:
        payload["reference_compressed"] = reference_compressed.as_manifest_entry()
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    logger.info("wrote manifest with %d fixtures to %s", len(summary.fixtures), manifest_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path, default=Path("tests/data"))
    parser.add_argument("--fixture-root", type=Path, default=DEFAULT_FIXTURE_ROOT)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--limit", type=int, default=None, help="derive only the first N, for a smoke run")
    parser.add_argument("--data-config", type=Path, default=DEFAULT_DATA_CONFIG)
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REFERENCE_COMPRESSED_FASTA, help="hg19 FASTA")
    parser.add_argument("--all", action="store_true", help="derive every discovered BAM, for the golden cohort")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(name)s: %(message)s")

    if not args.data_root.is_dir():
        logger.error("no cohort at %s; run `make download-test-data` first", args.data_root)
        return 1
    if not args.data_config.is_file():
        logger.error("test-data config not found: %s", args.data_config)
        return 1
    if not args.reference_fasta.is_file():
        reference_root = args.reference_fasta.parent.parent
        logger.error(
            "required pinned hg19 reference FASTA not found: %s; "
            "run `vntyper install-references -d %s --references hg19` before deriving CRAM fixtures",
            args.reference_fasta,
            reference_root,
        )
        return 1

    with snapshot_reference_fasta(args.reference_fasta) as validated_reference:
        summary = build_fixtures(
            args.samtools,
            args.data_root,
            args.fixture_root,
            args.limit,
            data_config=args.data_config,
            include_all=args.all,
        )
        reference_compressed = derive_reference_compressed_cram(
            args.samtools,
            DEFAULT_REFERENCE_COMPRESSED_SOURCE,
            args.reference_fasta,
            args.fixture_root,
            validated_reference=validated_reference,
        )
        build_reference_dependent_fixture(args.fixture_root)
        build_placed_flag12_fixture(args.fixture_root)
        build_indexed_safe_fixture(args.fixture_root)
        write_manifest(
            summary,
            args.manifest or args.fixture_root / "manifest.json",
            reference_compressed=reference_compressed,
        )

        if not summary.fixtures or summary.skipped:
            logger.error(
                "fixture derivation incomplete: %d verified, %d skipped",
                len(summary.fixtures),
                len(summary.skipped),
            )
            return 1
        logger.info(
            "derived %d verified CRAM fixtures: %.1f MiB from %.1f MiB of BAM (%.0f%%), %d skipped",
            len(summary.fixtures),
            summary.total_cram_bytes / 1024 / 1024,
            summary.total_source_bytes / 1024 / 1024,
            100.0 * summary.total_cram_bytes / summary.total_source_bytes,
            len(summary.skipped),
        )
        return 0


if __name__ == "__main__":
    sys.exit(main())
