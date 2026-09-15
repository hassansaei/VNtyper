"""Pure adapter tests for strict length depth I/O."""

from __future__ import annotations

import hashlib
import subprocess
import tempfile
from collections.abc import Callable
from pathlib import Path
from unittest import mock

import pysam
import pytest

from vntyper.scripts import length_depth_io
from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation
from vntyper.scripts.length_depth_io import (
    AlignmentEvidence,
    FragmentPositionEvidence,
    build_samtools_depth_argv,
    parse_samtools_depth,
    read_length_depth,
)
from vntyper.scripts.length_feature_provenance import (
    LengthFeatureContext,
    QueryInterval,
    decode_length_feature_context,
)

pytestmark = pytest.mark.unit


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _annotation(reference_sha256: str = "a" * 64) -> LengthAnnotation:
    return decode_length_annotation(
        {
            "schema_version": "length-annotation-v1",
            "assembly": "synthetic-build-v1",
            "contig": "synthetic-contig",
            "accepted_contig_aliases": ["synthetic-alias"],
            "reference_fasta_sha256": reference_sha256,
            "coordinate_system": "zero-based-half-open",
            "boundary_definition": "complete-core-plus-invariant-units-v1",
            "repeat_unit_bp": 1,
            "regions": {
                "CORE": [{"start": 1, "end": 3}],
                "INVARIANT": [{"start": 3, "end": 5}],
                "ARRAY": {"start": 1, "end": 5},
                "LEFT_FLANK": {"start": 0, "end": 1},
                "RIGHT_FLANK": {"start": 5, "end": 6},
            },
            "array_boundary_geometry": {"array_only_bp": 0, "target_only_bp": 0},
            "target_boundary_conversion_sha256": None,
            "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
            "annotation_provenance": "generated-synthetic-test-contract",
            "annotation_version": "synthetic-1",
        }
    )


def _context(input_sha256: str, annotation: LengthAnnotation, **changes: object) -> LengthFeatureContext:
    raw: dict[str, object] = {
        "schema_version": "length-feature-measurement-context-v1",
        "manifest_key": "synthetic-member-1",
        "input_sha256": input_sha256,
        "assembly": annotation.assembly,
        "assay_class": "synthetic-short-read",
        "input_scope": "regional",
        "original_contig": annotation.contig,
        "reference_fasta_sha256": annotation.reference_fasta_sha256,
        "annotation_sha256": annotation.sha256,
        "aligner": {
            "name": "synthetic-aligner",
            "version": "1.0",
            "arguments_sha256": "b" * 64,
            "primary_secondary_marking": "primary-plus-supplementary",
        },
        "fragment_reader": {
            "name": "pysam",
            "version": pysam.__version__,
            "htslib_version": pysam.version.__htslib_version__,
            "alignment_semantics": "explicit-filtered-aligned-pairs-v1",
        },
        "preprocessing_id": "synthetic-preprocessing-v1",
        "counting_policy": {
            "policy_id": "primary-mapq0-baseq0-overlap-count-v1",
            "samtools_revision": "samtools=1.20;htslib=1.23",
            "minimum_mapping_quality": 0,
            "minimum_base_quality": 0,
            "excluded_alignment_flags": ["UNMAP", "SECONDARY", "QCFAIL", "DUP"],
            "supplementary_alignment_policy": "included-unless-excluded-by-another-flag",
            "overlap_policy": "count-overlapping-mates-independently",
            "base_counting_policy": "one-per-aligned-covered-base",
            "zero_coverage_policy": "emit-zero-for-every-queried-position",
            "queried_intervals": [{"start": 0, "end": 6}],
        },
    }
    raw.update(changes)
    return decode_length_feature_context(raw)


def _paths(tmp_path: Path) -> tuple[Path, Path, Path, LengthAnnotation, LengthFeatureContext]:
    input_path = tmp_path / "input alignment.bam"
    reference_path = tmp_path / "reference genome.fa"
    samtools_path = tmp_path / "pinned tools" / "samtools"
    input_path.write_bytes(b"synthetic alignment bytes")
    reference_path.write_bytes(b">synthetic-contig\nAAAAAA\n")
    samtools_path.parent.mkdir()
    samtools_path.write_bytes(b"synthetic executable")
    samtools_path.chmod(0o700)
    annotation = _annotation(_sha256(reference_path))
    return input_path, reference_path, samtools_path, annotation, _context(_sha256(input_path), annotation)


def _alignment_evidence() -> AlignmentEvidence:
    return AlignmentEvidence(
        contig_length=6,
        assembly="synthetic-build-v1",
        positions={
            position: FragmentPositionEvidence(
                base_contributions=depth,
                fragment_ids=frozenset() if depth == 0 else frozenset({f"fragment-{position}"}),
            )
            for position, depth in enumerate((5, 0, 20, 5, 5, 5))
        },
    )


def _completed(argv: list[str] | tuple[str, ...], stdout: str, returncode: int = 0) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(argv, returncode, stdout, "synthetic stderr" if returncode else "")


def test_depth_argv_is_shell_free_and_preserves_paths_with_spaces(tmp_path: Path) -> None:
    argv = build_samtools_depth_argv(
        tmp_path / "tool dir" / "samtools",
        tmp_path / "input dir" / "input.bam",
        tmp_path / "reference dir" / "reference.fa",
        tmp_path / "private dir" / "union.bed",
    )

    assert argv == (
        str(tmp_path / "tool dir" / "samtools"),
        "depth",
        "-a",
        "-q",
        "0",
        "-Q",
        "0",
        "-b",
        str(tmp_path / "private dir" / "union.bed"),
        "--reference",
        str(tmp_path / "reference dir" / "reference.fa"),
        str(tmp_path / "input dir" / "input.bam"),
    )
    assert all("'" not in argument and '"' not in argument for argument in argv)


def test_depth_parser_converts_one_based_positions_and_attaches_fragment_evidence() -> None:
    support = {
        0: FragmentPositionEvidence(2, frozenset({"fragment-a"})),
        1: FragmentPositionEvidence(0, frozenset()),
    }

    depths = parse_samtools_depth(
        "synthetic-contig\t1\t2\nsynthetic-contig\t2\t0\n", "synthetic-contig", {0, 1}, support
    )

    assert [(item.position_zero_based, item.depth, item.supporting_fragment_ids) for item in depths] == [
        (0, 2, ("fragment-a",)),
        (1, 0, ()),
    ]


@pytest.mark.parametrize(
    ("stdout", "expected", "message"),
    [
        ("", {0}, "missing"),
        ("synthetic-contig\t1\t0\nsynthetic-contig\t1\t0\n", {0}, "duplicate"),
        ("synthetic-contig\t2\t0\n", {0}, "out-of-range"),
        ("wrong-contig\t1\t0\n", {0}, "contig"),
        ("synthetic-contig\t0\t0\n", {0}, "one-based"),
        ("synthetic-contig\tone\t0\n", {0}, "integer"),
        ("synthetic-contig\t+1\t0\n", {0}, "canonical"),
        ("synthetic-contig\t01\t0\n", {0}, "canonical"),
        ("synthetic-contig\t1\t 0\n", {0}, "canonical"),
        ("synthetic-contig\t1\t-1\n", {0}, "non-negative"),
        ("synthetic-contig\t1\n", {0}, "three tab-separated"),
        ("synthetic-contig\t1\t0\textra\n", {0}, "three tab-separated"),
    ],
)
def test_depth_parser_rejects_missing_duplicate_out_of_range_and_malformed_rows(
    stdout: str, expected: set[int], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        parse_samtools_depth(stdout, "synthetic-contig", expected, {})


def test_depth_parser_fails_when_fragment_reader_base_contributions_differ() -> None:
    support = {0: FragmentPositionEvidence(1, frozenset({"fragment-a"}))}

    with pytest.raises(ValueError, match="base contribution mismatch"):
        parse_samtools_depth("synthetic-contig\t1\t2\n", "synthetic-contig", {0}, support)


@pytest.mark.parametrize(
    ("constructor", "message"),
    [
        (lambda: FragmentPositionEvidence(True, frozenset()), "non-negative integer"),
        (lambda: FragmentPositionEvidence(0, {"fragment"}), "frozenset"),  # type: ignore[arg-type]
        (lambda: FragmentPositionEvidence(0, frozenset({" fragment"})), "trimmed"),
        (lambda: FragmentPositionEvidence(0, frozenset({"fragment"})), "cannot exceed"),
        (lambda: AlignmentEvidence(True, None, {}), "positive integer"),
        (lambda: AlignmentEvidence(1, " build", {}), "trimmed"),
        (lambda: AlignmentEvidence(1, None, {True: FragmentPositionEvidence(0, frozenset())}), "positions"),
        (lambda: AlignmentEvidence(1, None, {0: object()}), "FragmentPositionEvidence"),  # type: ignore[dict-item]
    ],
)
def test_evidence_types_reject_bool_and_malformed_values(constructor: Callable[[], object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        constructor()


def test_public_evidence_boundaries_revalidate_corrupted_typed_values() -> None:
    corrupted = FragmentPositionEvidence(1, frozenset({"fragment"}))
    object.__setattr__(corrupted, "base_contributions", True)

    with pytest.raises(ValueError, match="non-negative integer"):
        AlignmentEvidence(1, None, {0: corrupted})
    with pytest.raises(ValueError, match="non-negative integer"):
        parse_samtools_depth("synthetic-contig\t1\t1\n", "synthetic-contig", {0}, {0: corrupted})


@pytest.mark.parametrize(
    ("stdout", "contig", "positions", "evidence", "message"),
    [
        (b"", "synthetic-contig", {0}, {}, "must be text"),
        ("", " synthetic-contig", {0}, {}, "trimmed"),
        ("", "synthetic-contig", {True}, {}, "non-negative integers"),
        ("", "synthetic-contig", frozenset({0}), {}, "must be a set"),
        ("", "synthetic-contig", set(), [], "position mapping"),
        ("synthetic-contig\t1\t0\n", "synthetic-contig", {0}, {}, "exact expected"),
        ("synthetic-contig\t1\t0\n", "synthetic-contig", {0}, {0: object()}, "FragmentPositionEvidence"),
    ],
)
def test_depth_parser_revalidates_all_public_inputs(
    stdout: object,
    contig: str,
    positions: object,
    evidence: object,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        parse_samtools_depth(stdout, contig, positions, evidence)  # type: ignore[arg-type]


def test_fragment_reader_scans_explicit_aligned_bases_and_deduplicates_fragments() -> None:
    class Header:
        def to_dict(self) -> dict[str, object]:
            return {"SQ": [{"SN": "synthetic-contig", "LN": 3, "AS": "synthetic-build-v1"}]}

    class Record:
        def __init__(
            self,
            name: str,
            flag: int,
            pairs: list[tuple[int | None, int | None]],
            *,
            read_group: str | None,
            qualities: list[int] | None,
        ) -> None:
            self.query_name = name
            self.flag = flag
            self._pairs = pairs
            self._read_group = read_group
            self.query_qualities = qualities

        def get_aligned_pairs(self, *, matches_only: bool) -> list[tuple[int | None, int | None]]:
            assert matches_only is False
            return self._pairs

        def has_tag(self, tag: str) -> bool:
            assert tag == "RG"
            return self._read_group is not None

        def get_tag(self, tag: str) -> str:
            assert tag == "RG"
            assert self._read_group is not None
            return self._read_group

    included_mate_one = Record(
        "same-pair",
        0,
        [(0, 0), (None, 1), (1, None), (1, 2), (2, 3)],
        read_group="rg-1",
        qualities=[0, 0, 0],
    )
    included_mate_two = Record("same-pair", 0, [(0, 0)], read_group="rg-1", qualities=None)
    no_read_group = Record("other", 0x800, [(0, 1)], read_group=None, qualities=[0])
    excluded = Record("excluded", 0x100, [(0, 0)], read_group="rg-1", qualities=[0])

    class Alignment:
        header = Header()

        def __enter__(self) -> Alignment:
            return self

        def __exit__(self, *args: object) -> None:
            return None

        def fetch(self, contig: str, start: int, end: int) -> list[Record]:
            assert (contig, start, end) == ("synthetic-contig", 0, 3)
            return [included_mate_one, included_mate_two, no_read_group, excluded]

    with mock.patch("vntyper.scripts.length_depth_io.pysam.AlignmentFile", return_value=Alignment()) as opener:
        evidence = length_depth_io._read_alignment_evidence(
            Path("/synthetic/input.bam"),
            Path("/synthetic/reference.fa"),
            "synthetic-contig",
            (QueryInterval(0, 3),),
            {0, 1, 2},
            is_cram=False,
        )

    opener.assert_called_once_with("/synthetic/input.bam", "rb", reference_filename="/synthetic/reference.fa")
    assert evidence.contig_length == 3
    assert evidence.assembly == "synthetic-build-v1"
    assert evidence.positions[0].base_contributions == 2
    assert len(evidence.positions[0].fragment_ids) == 1
    assert evidence.positions[1].base_contributions == 1
    assert len(evidence.positions[1].fragment_ids) == 1
    assert evidence.positions[2].base_contributions == 1


def test_fragment_reader_wraps_invalid_header_or_record_content() -> None:
    class Header:
        def to_dict(self) -> dict[str, object]:
            return {"SQ": []}

    alignment = mock.MagicMock()
    alignment.__enter__.return_value = alignment
    alignment.header = Header()
    with (
        mock.patch("vntyper.scripts.length_depth_io.pysam.AlignmentFile", return_value=alignment),
        pytest.raises(RuntimeError, match="failed to read alignment records"),
    ):
        length_depth_io._read_alignment_evidence(
            Path("/synthetic/input.cram"),
            Path("/synthetic/reference.fa"),
            "synthetic-contig",
            (QueryInterval(0, 1),),
            {0},
            is_cram=True,
        )


def test_read_length_depth_validates_then_runs_exact_union_and_removes_private_tempdir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path, reference_path, samtools_path, annotation, context = _paths(tmp_path)
    temporary_parent = tmp_path / "temporary parent"
    temporary_parent.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(temporary_parent))
    observed_bed: list[str] = []
    observed_depth_argv: list[str] = []

    def run(argv: list[str] | tuple[str, ...], **kwargs: object) -> subprocess.CompletedProcess[str]:
        assert kwargs == {"capture_output": True, "check": False, "text": True, "env": None, "shell": False}
        if argv[1:] == ("--version",):
            return _completed(argv, "samtools 1.20\nUsing htslib 1.23\n")
        observed_depth_argv.extend(argv)
        bed_path = Path(argv[argv.index("-b") + 1])
        observed_bed.append(bed_path.read_text(encoding="utf-8"))
        return _completed(
            argv,
            "".join(
                f"synthetic-contig\t{position + 1}\t{depth}\n" for position, depth in enumerate((5, 0, 20, 5, 5, 5))
            ),
        )

    with (
        mock.patch("vntyper.scripts.length_depth_io.subprocess.run", side_effect=run),
        mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence", return_value=_alignment_evidence()),
    ):
        depths = read_length_depth(input_path, reference_path, annotation, context, samtools_path)

    assert [item.depth for item in depths] == [5, 0, 20, 5, 5, 5]
    assert observed_bed == ["synthetic-contig\t0\t6\n"]
    assert observed_depth_argv == list(
        build_samtools_depth_argv(samtools_path, input_path, reference_path, Path(observed_depth_argv[8]))
    )
    assert list(temporary_parent.iterdir()) == []


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"assembly": "wrong-build"}, "assembly mismatch"),
        ({"original_contig": "wrong-contig"}, "contig mismatch"),
        ({"annotation_sha256": "f" * 64}, "annotation digest mismatch"),
    ],
)
def test_contract_mismatch_fails_before_any_tool_or_alignment_io(
    tmp_path: Path, change: dict[str, object], message: str
) -> None:
    input_path, reference_path, samtools_path, annotation, context = _paths(tmp_path)
    context = _context(_sha256(input_path), annotation, **change)

    with (
        mock.patch("vntyper.scripts.length_depth_io.subprocess.run") as run,
        mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence") as alignment_read,
        pytest.raises(ValueError, match=message),
    ):
        read_length_depth(input_path, reference_path, annotation, context, samtools_path)
    run.assert_not_called()
    alignment_read.assert_not_called()


@pytest.mark.parametrize(("which", "message"), [("input", "input digest"), ("reference", "reference digest")])
def test_file_digest_mismatch_fails_before_tool_or_alignment_io(tmp_path: Path, which: str, message: str) -> None:
    input_path, reference_path, samtools_path, annotation, context = _paths(tmp_path)
    if which == "input":
        input_path.write_bytes(b"changed")
    else:
        reference_path.write_bytes(b"changed")

    with (
        mock.patch("vntyper.scripts.length_depth_io.subprocess.run") as run,
        mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence") as alignment_read,
        pytest.raises(ValueError, match=message),
    ):
        read_length_depth(input_path, reference_path, annotation, context, samtools_path)
    run.assert_not_called()
    alignment_read.assert_not_called()


@pytest.mark.parametrize(
    ("reference_text", "message"),
    [
        (">different-contig\nAAAAAA\n", "missing the original input contig"),
        (">synthetic-contig\nAAAAA\n", "does not cover"),
        (">synthetic-contig\nAAA\n>synthetic-contig\nAAA\n", "duplicate contig"),
        ("AAAAAA\n", "malformed"),
    ],
)
def test_reference_fasta_identity_and_length_fail_before_tool_io(
    tmp_path: Path, reference_text: str, message: str
) -> None:
    input_path, reference_path, samtools_path, _, _ = _paths(tmp_path)
    reference_path.write_text(reference_text, encoding="ascii")
    annotation = _annotation(_sha256(reference_path))
    context = _context(_sha256(input_path), annotation)

    with (
        mock.patch("vntyper.scripts.length_depth_io.subprocess.run") as run,
        mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence") as alignment_read,
        pytest.raises(ValueError, match=message),
    ):
        read_length_depth(input_path, reference_path, annotation, context, samtools_path)
    run.assert_not_called()
    alignment_read.assert_not_called()


def test_header_contig_length_and_assembly_are_validated_before_depth_process(tmp_path: Path) -> None:
    input_path, reference_path, samtools_path, annotation, context = _paths(tmp_path)

    for evidence, message in (
        (AlignmentEvidence(5, "synthetic-build-v1", {}), "contig length"),
        (AlignmentEvidence(6, "wrong-build", {}), "header assembly"),
    ):
        with (
            mock.patch(
                "vntyper.scripts.length_depth_io.subprocess.run",
                return_value=_completed((str(samtools_path), "--version"), "samtools 1.20\nUsing htslib 1.23\n"),
            ) as run,
            mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence", return_value=evidence),
            pytest.raises(ValueError, match=message),
        ):
            read_length_depth(input_path, reference_path, annotation, context, samtools_path)
        assert run.call_count == 1
        assert run.call_args.args[0][1:] == ("--version",)


def test_cram_uses_reference_environment_seam_and_restores_on_depth_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path, reference_path, samtools_path, annotation, context = _paths(tmp_path)
    input_path.write_bytes(b"CRAM" + b" synthetic alignment bytes")
    context = _context(_sha256(input_path), annotation)
    temporary_parent = tmp_path / "temporary parent"
    temporary_parent.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(temporary_parent))

    def run(argv: list[str] | tuple[str, ...], **kwargs: object) -> subprocess.CompletedProcess[str]:
        if argv[1:] == ("--version",):
            return _completed(argv, "samtools 1.20\nUsing htslib 1.23\n")
        return _completed(argv, "", returncode=1)

    with (
        mock.patch("vntyper.scripts.length_depth_io.subprocess.run", side_effect=run),
        mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence", return_value=_alignment_evidence()),
        mock.patch("vntyper.scripts.length_depth_io.pin_reference_resolution", return_value="old") as pin,
        mock.patch("vntyper.scripts.length_depth_io.restore_reference_resolution") as restore,
        pytest.raises(RuntimeError, match="samtools depth failed"),
    ):
        read_length_depth(input_path, reference_path, annotation, context, samtools_path)

    pin.assert_called_once_with(
        {"cram": {"allow_ambient_reference_resolution": False, "local_ref_path": str(reference_path)}}
    )
    restore.assert_called_once_with("old")
    assert list(temporary_parent.iterdir()) == []


@pytest.mark.parametrize(
    ("version_output", "fragment_reader", "message"),
    [
        ("samtools 1.21\nUsing htslib 1.23\n", None, "samtools revision"),
        ("samtools 1.20\n", None, "truncated"),
        (
            "samtools 1.20\nUsing htslib 1.23\n",
            {
                "name": "pysam",
                "version": "999.0",
                "htslib_version": pysam.version.__htslib_version__,
                "alignment_semantics": "explicit-filtered-aligned-pairs-v1",
            },
            "pysam version",
        ),
        (
            "samtools 1.20\nUsing htslib 1.23\n",
            {
                "name": "pysam",
                "version": pysam.__version__,
                "htslib_version": "999.0",
                "alignment_semantics": "explicit-filtered-aligned-pairs-v1",
            },
            "pysam htslib version",
        ),
    ],
)
def test_external_and_fragment_reader_versions_must_match_before_alignment_io(
    tmp_path: Path,
    version_output: str,
    fragment_reader: dict[str, object] | None,
    message: str,
) -> None:
    input_path, reference_path, samtools_path, annotation, context = _paths(tmp_path)
    if fragment_reader is not None:
        context = _context(_sha256(input_path), annotation, fragment_reader=fragment_reader)

    with (
        mock.patch(
            "vntyper.scripts.length_depth_io.subprocess.run",
            return_value=_completed((str(samtools_path), "--version"), version_output),
        ),
        mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence") as alignment_read,
        pytest.raises((ValueError, RuntimeError), match=message),
    ):
        read_length_depth(input_path, reference_path, annotation, context, samtools_path)
    alignment_read.assert_not_called()


def test_input_change_during_measurement_is_detected_and_private_tempdir_is_removed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path, reference_path, samtools_path, annotation, context = _paths(tmp_path)
    temporary_parent = tmp_path / "temporary parent"
    temporary_parent.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(temporary_parent))

    def run(argv: list[str] | tuple[str, ...], **kwargs: object) -> subprocess.CompletedProcess[str]:
        if argv[1:] == ("--version",):
            return _completed(argv, "samtools 1.20\nUsing htslib 1.23\n")
        input_path.write_bytes(b"same size changed bytes!!")
        return _completed(
            argv,
            "".join(
                f"synthetic-contig\t{position + 1}\t{depth}\n" for position, depth in enumerate((5, 0, 20, 5, 5, 5))
            ),
        )

    with (
        mock.patch("vntyper.scripts.length_depth_io.subprocess.run", side_effect=run),
        mock.patch("vntyper.scripts.length_depth_io._read_alignment_evidence", return_value=_alignment_evidence()),
        pytest.raises(RuntimeError, match="changed during"),
    ):
        read_length_depth(input_path, reference_path, annotation, context, samtools_path)

    assert list(temporary_parent.iterdir()) == []
