"""
kestrel_vcf_contract.py

Whether a VCF Kestrel produced can be parsed into variant records at all (#223).

Split out of ``file_processing.py``, which is about *transforming* VCFs -- both of its
public functions read one file and write another. Deciding whether a file is readable in
the first place is a different responsibility, it is the only pure decision on this path,
and it is the one a future VCF contract check belongs beside. ``file_processing`` keeps the
per-record allele-shape contract, which is about a record rather than a file.

**Why this exists.** ``run_kestrel`` used to treat "the VCF path exists" as success. A
Kestrel invocation that exits 0 and writes a file nothing can parse therefore satisfied the
completion check, both derived frames came back empty, and ``output_empty_result`` emitted
the ``Negative`` placeholder -- exit 0, a normal-looking report and no ERROR record
anywhere. The severe case is not the zero-byte file in the issue title: ``filter_vcf``
copies data lines through whether or not a header is present, so a VCF that lost its header
while *carrying real indel records* had those records silently discarded and a positive
was reported as a negative.

**The contract is ordered and structural, not a presence check.**
``variant_parsing.read_vcf_without_comments`` takes any line beginning ``#CHROM`` as its
column definition without validating the columns, and collects a data line only once that
line has been seen. So two files pass a naive check and still parse to nothing: one whose
``#CHROM`` line is truncated, and one whose header arrives after some of its records.

**A valid header with no records is usable.** That is a genuine empty result and reporting
it as a negative is correct. Keeping those two apart -- "ran and found nothing" versus
"produced nothing readable" -- is the entire point of the module.

Functions:
    describe_unusable_vcf: A VCF path to a reason it cannot be parsed, or ``None``
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

#: The eight mandatory VCF columns, in order, as VCFv4.2 fixes them.
#:
#: A ``#CHROM`` line carrying fewer than these is truncated, not a header. It matters because
#: ``variant_parsing.read_vcf_without_comments`` accepts *any* line beginning ``#CHROM`` as its
#: column definition, so a bare ``#CHROM`` sets a one-element column list, collects no records,
#: and yields an empty frame - which two steps later is the ``Negative`` placeholder (#223).
_MANDATORY_VCF_COLUMNS: tuple[str, ...] = (
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
)


def describe_unusable_vcf(vcf_path: str | Path) -> str | None:
    """
    Say why a VCF cannot be parsed into variant records, or ``None`` if it can.

    This is the second half of this module's job as the trust boundary for Kestrel's output.
    :func:`_assert_kestrel_allele_contract` checks that a *record* has the shape the pinned
    build emits; this checks that the *file* has a header those records can be read against.

    Three things make a VCF unusable here, and each is a distinct silent-negative path:

    1. **It cannot be read at all.** ``read_vcf_without_comments`` converts any read failure
       into an empty frame, which is indistinguishable from "no variants found".
    2. **It has no ``#CHROM`` line, or a truncated one.** The parser takes any line beginning
       ``#CHROM`` as its column definition without checking the columns, so a bare ``#CHROM``
       is accepted and yields an empty frame.
    3. **A data line precedes the ``#CHROM`` line.** The contract is *ordered*, not merely
       present: the parser collects a data line only once the header has been seen, so records
       arriving before it are discarded in silence and a file whose header is late parses to
       nothing while passing any presence check.

    Why this matters: :func:`filter_vcf` copies data lines through whether or not a header is
    present, so a headerless VCF yields headerless derived files, and
    ``kestrel_genotyping.process_kestrel_output`` renders two empty frames as the ``Negative``
    placeholder. A Kestrel run that exits 0 after losing its header would otherwise turn real
    indel records into a confident negative genotype with no ERROR logged at all (#223).

    A valid header with **no records** is usable. That is a genuine empty result, and reporting
    it as a negative is correct - the distinction between "ran and found nothing" and "produced
    nothing readable" is the entire point.

    Args:
        vcf_path (str | pathlib.Path): Path to the VCF to inspect.

    Returns:
        str | None: ``None`` when the file can be parsed into records. Otherwise a lowercase
        clause naming the problem, phrased to read after "but " in a caller's message.
    """
    records_before_header = 0
    try:
        with open(vcf_path) as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    columns = tuple(line.rstrip("\n").split("\t"))
                    if columns[: len(_MANDATORY_VCF_COLUMNS)] != _MANDATORY_VCF_COLUMNS:
                        return (
                            f"its #CHROM line is truncated - it carries {len(columns)} column(s) rather "
                            f"than the {len(_MANDATORY_VCF_COLUMNS)} mandatory ones, so no record can be "
                            "read against it"
                        )
                    if records_before_header:
                        return (
                            f"{records_before_header} data line(s) appear before the #CHROM header, and "
                            "every one of them is discarded when the file is parsed"
                        )
                    return None
                if line.startswith("##") or not line.strip():
                    continue
                records_before_header += 1
    except (OSError, UnicodeDecodeError) as exc:
        # UnicodeDecodeError is caught alongside OSError deliberately. This function's contract
        # is that it *returns a reason*; letting a decode error escape would hand the caller an
        # exception it does not expect, from a function whose whole job is to classify
        # unreadable output. Corrupt bytes where a VCF should be is exactly the condition being
        # detected, not an unexpected one (#223).
        return f"it could not be read ({exc})"

    if records_before_header:
        return (
            f"it has no #CHROM header line, so all {records_before_header} data line(s) it carries are "
            "discarded when the file is parsed"
        )
    return "it has no #CHROM header line"
