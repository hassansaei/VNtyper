#!/usr/bin/env python3
"""Derive CRAM fixtures from the BAM test cohort, and prove the derivation is lossless.

Why this exists
---------------
VNtyper accepts CRAM but no CRAM has ever been exercised by a test (#188), so the CRAM
branch of ``process_bam_to_fastq`` -- including its whole-file unmapped-read scan and the
merge that follows it -- has never run against real data. Every claim about the CRAM path
has therefore been an argument from source reading.

Why the fixtures are derived rather than committed
--------------------------------------------------
``/tests/data/`` is git-ignored; the cohort ships as a Zenodo archive. A committed CRAM
would need a new archive, new checksums, and a second copy of every read. Deriving each
CRAM from the BAM beside it costs about a second per sample and buys a much stronger
assertion than a standalone fixture ever could:

    the same sample, as BAM and as CRAM, must produce the same genotype.

That is a real equivalence test. A shipped CRAM fixture could only be compared against a
recorded expectation, which is exactly the kind of evidence this repository has learned
not to trust.

Why ``no_ref=1``
----------------
The cohort's BAM headers carry no ``M5`` tags, so htslib cannot resolve a reference by
digest: no ``REF_PATH``/``REF_CACHE`` lookup and no refget request can succeed. A
conventionally reference-compressed CRAM would therefore need an explicit ``-T`` at read
time, and ``process_bam_to_fastq`` passes ``cram_ref_option = ""`` unconditionally -- so
such a CRAM could not be decoded by the pipeline at all. ``no_ref=1`` stores sequences
verbatim, needs no reference to write or read, and round-trips byte-for-byte.

``embed_ref=2`` was measured and rejected: htslib cannot derive a consensus reference from
a region subset, warns, falls back to non-ref mode part-way through a container, and the
round trip is **not** lossless.

What this script does NOT establish
-----------------------------------
A ``no_ref`` CRAM exercises the container format, the CRAM decoder, ``.crai`` indexing and
the unmapped-read scan. It does **not** exercise reference resolution, because it needs
none. The separate failure mode of a CRAM whose reference is unavailable -- the ordinary
externally-referenced CRAM a diagnostic lab would send -- is a different fixture and a
different test; see ``build_reference_dependent_fixture``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import logging
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pysam

logger = logging.getLogger("cram_fixtures")

#: Encoding options for the derived fixtures. See the module docstring for why this is
#: ``no_ref`` rather than a reference-compressed or ``embed_ref`` CRAM.
CRAM_WRITE_OPTIONS = ("no_ref=1",)

#: Where derived fixtures are written, mirroring the source layout underneath.
DEFAULT_FIXTURE_ROOT = Path("tests/data/cram")
DEFAULT_DATA_CONFIG = Path("tests/test_data_config.json")

# Pysam's typed interface exposes ``AlignmentFile`` but not its htslib-command wrappers.
# Purpose-built fixtures use the latter deliberately so the unit tier never needs PATH's
# external samtools binary.
pysam_any: Any = pysam


class LossyConversionError(RuntimeError):
    """A derived CRAM did not round-trip to the records of its source BAM."""


@dataclass
class Fixture:
    """One derived CRAM and the evidence that it is faithful to its source."""

    source_bam: Path
    cram: Path
    records: int
    unmapped_pairs: int
    #: Digest of the decoded record stream, identical for the BAM and the CRAM.
    record_digest: str
    source_bytes: int
    cram_bytes: int

    def as_manifest_entry(self) -> dict[str, object]:
        return {
            "source_bam": str(self.source_bam),
            "cram": str(self.cram),
            "records": self.records,
            "unmapped_pairs": self.unmapped_pairs,
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


def _run(argv: list[str]) -> str:
    """Run a samtools command, returning stdout and raising with stderr on failure."""
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{' '.join(argv)} exited {completed.returncode}: {completed.stderr.strip()}")
    return completed.stdout


def _record_digest(samtools: str, alignment: Path) -> tuple[str, int]:
    """Digest every decoded record of an alignment, with no reference supplied.

    Reading the CRAM with no ``-T`` is the point: it is what the pipeline does. If the
    fixture needed a reference this call would fail rather than quietly succeed.
    """
    digest = hashlib.sha256()
    count = 0
    with subprocess.Popen(
        [samtools, "view", str(alignment)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            digest.update(line.encode())
            count += 1
        stderr = proc.stderr.read() if proc.stderr else ""
        # Wait before reading returncode: draining stdout does not reap the child, and an
        # unreaped child reports returncode None, which would compare unequal to 0 and
        # turn every success into a spurious failure.
        proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools view {alignment} exited {proc.returncode}: {stderr.strip()}")
    return digest.hexdigest(), count


def discover_source_bams(data_root: Path, fixture_root: Path) -> list[Path]:
    """Every BAM in the cohort, excluding anything already under the fixture root."""
    bams = sorted(p for p in data_root.rglob("*.bam") if fixture_root not in p.parents)
    logger.info("discovered %d source BAMs under %s", len(bams), data_root)
    return bams


def _select_source_bams(discovered: list[Path], *, data_config: Path, data_root: Path, include_all: bool) -> list[Path]:
    """Import the selection policy in package and direct-script execution modes."""
    if __package__:
        from .cram_fixture_selection import select_source_bams as package_selector

        return package_selector(discovered, data_config=data_config, data_root=data_root, include_all=include_all)
    from cram_fixture_selection import select_source_bams as direct_selector

    return direct_selector(discovered, data_config=data_config, data_root=data_root, include_all=include_all)


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

    # -f 12 is the flag combination the pipeline's unmapped-read scan uses; record it so a
    # gate case can assert the CRAM presents the same unmapped set the BAM does.
    unmapped = int(_run([samtools, "view", "-c", "-f", "12", str(cram)]).strip())

    logger.info(
        "%s -> %s  (%d records, %d unmapped pairs, %.1f%% of BAM size)",
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
        unmapped_pairs=unmapped,
        record_digest=bam_digest,
        source_bytes=bam.stat().st_size,
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


def write_manifest(summary: Summary, manifest_path: Path) -> None:
    """Record what was derived, so a gate run can cite it rather than re-derive it."""
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "encoding": list(CRAM_WRITE_OPTIONS),
        "verified": "every fixture's decoded record stream digests identically to its source BAM",
        "fixtures": [f.as_manifest_entry() for f in summary.fixtures],
        "skipped": [{"bam": str(p), "reason": r} for p, r in summary.skipped],
    }
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
    parser.add_argument("--all", action="store_true", help="derive every discovered BAM, for the golden cohort")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(name)s: %(message)s")

    if not args.data_root.is_dir():
        logger.error("no cohort at %s; run `make download-test-data` first", args.data_root)
        return 2

    summary = build_fixtures(
        args.samtools,
        args.data_root,
        args.fixture_root,
        args.limit,
        data_config=args.data_config,
        include_all=args.all,
    )
    build_reference_dependent_fixture(args.fixture_root)
    build_placed_flag12_fixture(args.fixture_root)
    write_manifest(summary, args.manifest or args.fixture_root / "manifest.json")

    if not summary.fixtures:
        logger.error("no fixtures derived")
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
