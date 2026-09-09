"""Pure decision logic for unrepresentable and representation-limited variant alleles.

Certain variant classes (such as complex deletion-insertions like 54_56delinsAT or
large/complex insertions like ins25bp) cannot be faithfully represented within
Kestrel's atomic VCF variant model (1:1, 1:N, N:1). When such an allele is
encountered, the pipeline abstains from emitting a fabricated or misleading
positional name, emits a 'representation-limited' disposition, retains 100% frameshift
detection sensitivity, and tracks the disposition separately from erroneous calls.

Research use only.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Iterable

DISPOSITION_REPRESENTATION_LIMITED = "representation-limited"
FLAG_ALLELE_UNREPRESENTABLE = "allele-unrepresentable-in-vcf"
NOTE_REPRESENTATION_LIMITED = (
    "Positional name withheld: allele molecular class cannot be represented in VCF; requires validation"
)


def is_unrepresentable_molecular_class(
    event: str | None,
    ref_span: int = 0,
    inserted: int = 0,
    flags: Iterable[str] = (),
) -> bool:
    """Determine whether an edit's molecular class cannot be represented in VCF.

    Args:
        event: Variant event string (e.g. 'delins', 'insertion', 'deletion').
        ref_span: Span of reference bases deleted or replaced.
        inserted: Number of bases inserted.
        flags: Existing flags associated with the call.

    Returns:
        bool: True if the molecular class is unrepresentable.
    """
    return bool(FLAG_ALLELE_UNREPRESENTABLE in flags or event == "delins" or (ref_span > 0 and inserted > 0))


def format_representation_limited_name(net_length: int) -> str:
    """Format the display string for a representation-limited call.

    Args:
        net_length: Net length change of the allele (positive for net insertion,
            negative for net deletion, zero for in-frame substitution).

    Returns:
        str: Non-empty description that satisfies downstream table parsers while
        avoiding digit-initial positional naming (so it is not treated as an HGVS coordinate).
    """
    if net_length == 0:
        return DISPOSITION_REPRESENTATION_LIMITED
    sign = "+" if net_length > 0 else "-"
    kind = "in-frame" if net_length % 3 == 0 else "frameshift"
    return f"{kind} {sign}{abs(net_length)}, {DISPOSITION_REPRESENTATION_LIMITED}"


def format_masthead_representation_limited(net_length: int) -> str:
    """Format the patient-facing/clinician masthead string for representation-limited results.

    Args:
        net_length: Net length change.

    Returns:
        str: Clear clinical statement distinguishing frameshift detection from positional naming.
    """
    if net_length == 0:
        return f"Variant detected (positional naming withheld: {DISPOSITION_REPRESENTATION_LIMITED})"
    sign = "+" if net_length > 0 else "-"
    kind = "In-frame variant" if net_length % 3 == 0 else "Frameshift"
    return f"{kind} detected ({sign}{abs(net_length)} bp, positional naming withheld: {DISPOSITION_REPRESENTATION_LIMITED})"
