"""Generated-alignment parity tests for the strict length depth adapter."""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation
from vntyper.scripts.length_depth_io import read_length_depth
from vntyper.scripts.length_feature_provenance import LengthFeatureContext, decode_length_feature_context

pytestmark = pytest.mark.integration


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _annotation(contig: str, length: int, reference_sha256: str) -> LengthAnnotation:
    if length == 20:
        core = [{"start": 2, "end": 10}]
        invariant = [{"start": 10, "end": 18}]
        array_interval = {"start": 2, "end": 18}
        left = {"start": 0, "end": 2}
        right = {"start": 18, "end": 20}
    else:
        core = [{"start": 1, "end": 2}]
        invariant = [{"start": 2, "end": 4}]
        array_interval = {"start": 1, "end": 4}
        left = {"start": 0, "end": 1}
        right = {"start": 4, "end": 5}
    return decode_length_annotation(
        {
            "schema_version": "length-annotation-v1",
            "assembly": "synthetic-build-v1",
            "contig": contig,
            "accepted_contig_aliases": [],
            "reference_fasta_sha256": reference_sha256,
            "coordinate_system": "zero-based-half-open",
            "boundary_definition": "complete-core-plus-invariant-units-v1",
            "repeat_unit_bp": 1,
            "regions": {
                "CORE": core,
                "INVARIANT": invariant,
                "ARRAY": array_interval,
                "LEFT_FLANK": left,
                "RIGHT_FLANK": right,
            },
            "array_boundary_geometry": {"array_only_bp": 0, "target_only_bp": 0},
            "target_boundary_conversion_sha256": None,
            "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
            "annotation_provenance": "generated-synthetic-integration-contract",
            "annotation_version": "synthetic-1",
        }
    )


def _samtools_revision(samtools_path: Path) -> str:
    completed = subprocess.run(
        [str(samtools_path), "--version"], capture_output=True, check=True, text=True, shell=False
    )
    lines = completed.stdout.splitlines()
    return f"samtools={lines[0].split()[1]};htslib={lines[1].split()[2]}"


def _context(
    bam_path: Path,
    annotation: LengthAnnotation,
    samtools_revision: str,
    length: int,
) -> LengthFeatureContext:
    return decode_length_feature_context(
        {
            "schema_version": "length-feature-measurement-context-v1",
            "manifest_key": f"generated-{annotation.contig}",
            "input_sha256": _sha256(bam_path),
            "assembly": annotation.assembly,
            "assay_class": "synthetic-short-read",
            "input_scope": "regional",
            "original_contig": annotation.contig,
            "reference_fasta_sha256": annotation.reference_fasta_sha256,
            "annotation_sha256": annotation.sha256,
            "aligner": {
                "name": "synthetic-aligner",
                "version": "1.0",
                "arguments_sha256": "a" * 64,
                "primary_secondary_marking": "primary-plus-supplementary",
            },
            "fragment_reader": {
                "name": "pysam",
                "version": pysam.__version__,
                "htslib_version": pysam.version.__htslib_version__,
                "alignment_semantics": "explicit-filtered-aligned-pairs-v1",
            },
            "preprocessing_id": "generated-synthetic-v1",
            "counting_policy": {
                "policy_id": "primary-mapq0-baseq0-overlap-count-v1",
                "samtools_revision": samtools_revision,
                "minimum_mapping_quality": 0,
                "minimum_base_quality": 0,
                "excluded_alignment_flags": ["UNMAP", "SECONDARY", "QCFAIL", "DUP"],
                "supplementary_alignment_policy": "included-unless-excluded-by-another-flag",
                "overlap_policy": "count-overlapping-mates-independently",
                "base_counting_policy": "one-per-aligned-covered-base",
                "zero_coverage_policy": "emit-zero-for-every-queried-position",
                "queried_intervals": [{"start": 0, "end": length}],
            },
        }
    )


def _record(
    header: pysam.AlignmentHeader,
    name: str,
    start: int,
    cigar: tuple[tuple[int, int], ...],
    *,
    flag: int = 0,
    mapping_quality: int = 60,
    read_group: str = "rg-1",
    mate_start: int = -1,
) -> pysam.AlignedSegment:
    query_length = sum(length for operation, length in cigar if operation in {0, 1, 4, 7, 8})
    record = pysam.AlignedSegment(header)
    record.query_name = name
    record.query_sequence = "A" * query_length
    record.flag = flag
    record.reference_id = 0 if start >= 0 else -1
    record.reference_start = start
    record.mapping_quality = mapping_quality
    record.cigartuples = list(cigar)
    record.query_qualities = pysam.qualitystring_to_array("!" * query_length)
    record.next_reference_id = 0 if mate_start >= 0 else -1
    record.next_reference_start = mate_start
    record.set_tag("RG", read_group)
    return record


@pytest.fixture
def generated_alignment(tmp_path: Path) -> tuple[Path, Path, Path, str]:
    """Generate a coordinate-sorted indexed BAM and its local reference."""
    reference = tmp_path / "reference with spaces.fa"
    reference.write_text(">synthetic-contig\n" + "A" * 20 + "\n>unused-contig\nAAAAA\n", encoding="ascii")
    pysam.faidx(str(reference))  # type: ignore[attr-defined]
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": "synthetic-contig", "LN": 20, "AS": "synthetic-build-v1"},
                {"SN": "unused-contig", "LN": 5, "AS": "synthetic-build-v1"},
            ],
            "RG": [{"ID": "rg-1"}, {"ID": "rg-2"}],
        }
    )
    bam = tmp_path / "alignment with spaces.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
        for record in (
            _record(header, "mapq-and-baseq-zero", 0, ((0, 4),), mapping_quality=0),
            _record(header, "secondary", 0, ((0, 4),), flag=0x100),
            _record(header, "supplementary", 0, ((0, 4),), flag=0x800),
            _record(header, "qc-fail", 0, ((0, 4),), flag=0x200),
            _record(header, "duplicate", 0, ((0, 4),), flag=0x400),
            _record(header, "deletion", 4, ((0, 2), (2, 2), (0, 2))),
            _record(header, "reference-skip", 4, ((0, 2), (3, 2), (0, 2))),
            _record(header, "overlapping-pair", 10, ((0, 5),), flag=99, mate_start=12),
            _record(header, "overlapping-pair", 12, ((0, 5),), flag=147, mate_start=10),
            _record(header, "overlapping-pair", 15, ((0, 2),), read_group="rg-2"),
            _record(header, "unmapped", -1, ((0, 4),), flag=0x4),
        ):
            output.write(record)
    pysam.index(str(bam))  # type: ignore[attr-defined]
    samtools = Path(sys.executable).parent / "samtools"
    assert samtools.is_file()
    return bam, reference, samtools, _samtools_revision(samtools)


def test_generated_bam_matches_external_depth_for_flags_cigar_overlap_and_zero_policy(
    generated_alignment: tuple[Path, Path, Path, str],
) -> None:
    bam, reference, samtools, revision = generated_alignment
    annotation = _annotation("synthetic-contig", 20, _sha256(reference))
    context = _context(bam, annotation, revision, 20)

    depths = read_length_depth(bam, reference, annotation, context, samtools)

    assert [item.depth for item in depths] == [2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0]
    assert len(depths[12].supporting_fragment_ids or ()) == 1
    assert len(depths[15].supporting_fragment_ids or ()) == 2


def test_samtools_a_missing_unused_contig_positions_fails_instead_of_inventing_zeros(
    generated_alignment: tuple[Path, Path, Path, str],
) -> None:
    bam, reference, samtools, revision = generated_alignment
    annotation = _annotation("unused-contig", 5, _sha256(reference))
    context = _context(bam, annotation, revision, 5)

    with pytest.raises(ValueError, match="missing an expected position"):
        read_length_depth(bam, reference, annotation, context, samtools)


def test_generated_cram_preserves_bam_depth_and_fragment_identity_with_pinned_reference(
    generated_alignment: tuple[Path, Path, Path, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    bam, reference, samtools, revision = generated_alignment
    cram = bam.with_suffix(".cram")
    with (
        pysam.AlignmentFile(str(bam), "rb") as source,
        pysam.AlignmentFile(str(cram), "wc", header=source.header, reference_filename=str(reference)) as output,
    ):
        for record in source.fetch(until_eof=True):
            output.write(record)
    pysam.index(str(cram))  # type: ignore[attr-defined]
    annotation = _annotation("synthetic-contig", 20, _sha256(reference))
    bam_depths = read_length_depth(bam, reference, annotation, _context(bam, annotation, revision, 20), samtools)
    # A hostile ambient resolver must not affect either decoding engine.
    monkeypatch.setenv("REF_PATH", "https://synthetic.invalid/reference/%s")
    cram_depths = read_length_depth(cram, reference, annotation, _context(cram, annotation, revision, 20), samtools)

    assert cram_depths == bam_depths
    assert [item.depth for item in cram_depths] == [2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0]
    assert len(cram_depths[12].supporting_fragment_ids or ()) == 1
    assert len(cram_depths[15].supporting_fragment_ids or ()) == 2
    assert os.environ["REF_PATH"] == "https://synthetic.invalid/reference/%s"
