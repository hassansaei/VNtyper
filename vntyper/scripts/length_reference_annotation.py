"""Reviewed GRCh38 measurement geometry derived from public reference evidence."""

from __future__ import annotations

from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation

GRCH38_REFERENCE_FASTA_SHA256 = "04ee6db2e94ccc4daddc168453189d8c17a01454ba67ac0443a73a0401408ee0"
_TARGET_START = 155188486
_TARGET_END = 155192239
_PRE_REPEAT_SEQUENCE = (
    "AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG"
    "AGTATGACCAGCAGCGTACTCTCCAGCCACAGCCCCGGTTCAGGCTCCTCCACCACTCAG"
    "GGACAGGATGTCACTCTGGCCCCGGCCACGGAACCAGCTTCAGGTTCAGCTGCCACCTGG"
    "GGACAGGATGTCACCTCGGTCCCAGTCACCAGGCCAGCCCTGGGCTCCACCACCCCGCCA"
    "GCCCACGATGTCACCTCAGCCCCGGACAACAAGCCAGCCCCGGGCTCCACCGCCCCCCCA"
)
_AFTER_REPEAT_SEQUENCE = (
    "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCGGGCCCCGGGCTCCACCCCGGCCCCG"
    "GGCTCCACCGCCCCCCCAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCG"
    "GGCTCCACCGCCCCCCCAGCCCATGGTGTCACCTCGGCCCCGGACAACAGGCCCGCCTTG"
    "GGCTCCACCGCCCCTCCAGTCCACAATGTCACCTCGGCCTCAGGCTCTGCATCAGGCTCA"
)


def grch38_length_annotation() -> LengthAnnotation:
    """Return the reviewed public GRCh38 length measurement definition.

    The primary method defines five pre-repeats and four after-repeats.
    Their exact public sequences locate the two invariant intervals on the
    negative-strand GRCh38 reference. The intervening CORE contains sequence
    indels, so its bp length does not establish an integral reference count.

    Returns:
        Immutable annotation bound to the pinned public chromosome FASTA.
    """
    return decode_length_annotation(
        {
            "schema_version": "length-annotation-v1",
            "assembly": "GRCh38",
            "contig": "chr1",
            "accepted_contig_aliases": ["1"],
            "reference_fasta_sha256": GRCH38_REFERENCE_FASTA_SHA256,
            "coordinate_system": "zero-based-half-open",
            "boundary_definition": "complete-core-plus-invariant-units-v1",
            "repeat_unit_bp": 60,
            "regions": {
                "CORE": [{"start": 155188726, "end": 155191939}],
                "INVARIANT": [
                    {"start": 155188486, "end": 155188726},
                    {"start": 155191939, "end": 155192239},
                ],
                "ARRAY": {"start": 155188529, "end": 155192010},
                "LEFT_FLANK": {"start": 155188296, "end": 155188486},
                "RIGHT_FLANK": {"start": 155192239, "end": 155192429},
            },
            "array_boundary_geometry": {"array_only_bp": 0, "target_only_bp": 272},
            "target_boundary_conversion_sha256": None,
            "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
            "annotation_provenance": (
                "PMCID:PMC12458345;public-repeat-repository:"
                "pristanna/muc1repeats@74a8bab867ef984798991dcd9c848faef4371cc9;"
                "VNtyper-array-derivation:b476789a"
            ),
            "annotation_version": "grch38-public-reference-terminal-witness-v1",
        }
    )


def verify_grch38_length_reference(annotation: LengthAnnotation, target_forward_sequence: str) -> None:
    """Verify exact public terminal-unit witnesses in the target reference slice.

    Args:
        annotation: Annotation returned by :func:`grch38_length_annotation`.
        target_forward_sequence: Forward-strand sequence for the exact annotation
            target span ``chr1:155188487-155192239`` (one-based closed).

    Raises:
        ValueError: If the annotation, sequence type/length, alphabet, or either
            terminal witness differs from the reviewed reference definition.
    """
    if annotation != grch38_length_annotation():
        raise ValueError("length reference verifier requires the approved GRCh38 annotation")
    if not isinstance(target_forward_sequence, str):
        raise ValueError("length reference target sequence must be DNA text")
    sequence = target_forward_sequence.upper()
    if len(sequence) != _TARGET_END - _TARGET_START:
        raise ValueError("length reference target sequence length differs from the annotation")
    if any(base not in "ACGT" for base in sequence):
        raise ValueError("length reference target sequence must contain only A, C, G, and T")
    after_forward = _reverse_complement(_AFTER_REPEAT_SEQUENCE)
    pre_forward = _reverse_complement(_PRE_REPEAT_SEQUENCE)
    if not sequence.startswith(after_forward) or not sequence.endswith(pre_forward):
        raise ValueError("length reference terminal sequence witness mismatch")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
