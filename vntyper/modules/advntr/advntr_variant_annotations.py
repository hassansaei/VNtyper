"""Pure parsing of repeat-unit annotations from adVNTR variant state strings."""

import re
from collections.abc import Iterable

INSERTION_PATTERN: re.Pattern[str] = re.compile(r"^I(\d+)_([0-9]+)_([ACGT])_LEN(\d+)$")
DELETION_PATTERN: re.Pattern[str] = re.compile(r"^D(\d+)_([0-9]+)$")


def derive_ru_and_pos(variant_values: Iterable[object]) -> tuple[list[str], list[str]]:
    """Derive repeat-unit and position annotations from adVNTR state strings.

    Repeat-unit identity and position are encoded directly in each state part, so this
    parsing needs no repeat-unit FASTA. Compound states retain one comma-joined value per
    part, in input order. Unparseable and non-string values contribute ``"."`` to both
    annotations rather than aborting result processing.

    Args:
        variant_values: Variant state values, each optionally containing ``&``-joined
            parts.

    Returns:
        A pair containing the repeat-unit annotations and position annotations.
    """
    ru_annotations: list[str] = []
    pos_annotations: list[str] = []

    for variant in variant_values:
        if not isinstance(variant, str):
            ru_annotations.append(".")
            pos_annotations.append(".")
            continue

        ru_parts: list[str] = []
        pos_parts: list[str] = []
        for part in variant.split("&"):
            match = INSERTION_PATTERN.match(part.strip()) or DELETION_PATTERN.match(part.strip())
            if match is None:
                ru_parts.append(".")
                pos_parts.append(".")
                continue
            ru_parts.append(match.group(2))
            pos_parts.append(str(int(match.group(1))))

        ru_annotations.append(",".join(ru_parts))
        pos_annotations.append(",".join(pos_parts))

    return ru_annotations, pos_annotations
