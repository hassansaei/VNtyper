"""Input guards ``run_pipeline`` applies before it commits to a coordinate system.

``assembly_guard.reconcile_assembly`` computes a verdict and never raises: it is pure
by design, so that the *policy* - what a mismatch should do to a run - could land and
be reverted on its own. This module is that policy, extracted rather than added to
``pipeline.py``, which is already over the 650-line limit (AGENTS.md rule 3).

What it guards against: the MUC1 VNTR sits at chr1:155,158,000-155,163,000 in GRCh37
and chr1:155,184,000-155,194,000 in GRCh38, about 30 kb apart. Declaring the wrong
``--reference-assembly`` slices a region that does not contain the VNTR. Kestrel then
sees no supporting reads, ``kestrel_result.tsv`` is empty, and the report's configured
**negative default** applies. The run succeeds, exits 0, and says the patient does not
carry the variant.

Three rules, and the third is the one that keeps the guard safe to add:

1. A **decided disagreement** between two known builds stops the run.
2. **Undetermined** - unreadable header, unknown declared name, no chr1, a chr1 length
   matching neither build - warns and continues. It is neither a pass nor a failure,
   and treating it as a failure would reject inputs that are perfectly fine.
3. The header read itself is **non-fatal here**. ``extract_bam_header`` runs samtools
   with ``check=True``, so a CRAM whose reference cannot be resolved raises
   ``CalledProcessError`` - which is neither ``KeyError`` nor ``ValueError``, and so is
   not caught by ``get_region_string_with_fallback``. A guard that turned an unreadable
   header into a crash would be worse than the defect it closes.

FASTQ input is deliberately **not** guarded. A FASTQ has no header until it has been
aligned, and the header it then has describes the reference that was indexed for BWA,
not the sample. A mismatch there is a real misconfiguration with the same
false-negative consequence, but the remediation is different - so enforcing it with
this message would send users to fix an input file that is not wrong.

Functions:
    read_alignment_header: Read a BAM/CRAM header, or None if it cannot be read.
    enforce_declared_assembly: Apply the three rules above to a header.
"""

from __future__ import annotations

import logging

from vntyper.scripts.assembly_guard import (
    STATUS_MISMATCH,
    STATUS_UNDETERMINED,
    AssemblyVerdict,
    reconcile_assembly,
)
from vntyper.scripts.fastq_bam_processing import extract_bam_header, parse_contigs_from_header

logger = logging.getLogger(__name__)


def read_alignment_header(input_file: str, config: dict) -> str | None:
    """Read a BAM or CRAM header, returning None when samtools cannot.

    Args:
        input_file: Path to the BAM or CRAM file.
        config: Configuration mapping, for ``config["tools"]["samtools"]``.

    Returns:
        str | None: The header text, or None if it could not be read.
    """
    try:
        return extract_bam_header(input_file, config)
    except Exception as exc:
        logger.warning(f"Could not read the alignment header of {input_file} for assembly reconciliation: {exc}")
        return None


def enforce_declared_assembly(reference_assembly: str, header: str | None) -> AssemblyVerdict:
    """Stop the run when the declared assembly and the header disagree.

    Args:
        reference_assembly: The assembly the caller declared, verbatim, as
            ``run_pipeline`` received it.
        header: The alignment header text, or None if it could not be read.

    Returns:
        AssemblyVerdict: The verdict, so a caller can record it.

    Raises:
        ValueError: If the declared assembly and the header name different builds.
            The message names the declared string, both builds, the chr1 length that
            decided it, and the ``--reference-assembly`` value that would be right.
    """
    contigs = parse_contigs_from_header(header) if header else []
    verdict = reconcile_assembly(reference_assembly, contigs)

    if verdict.status == STATUS_MISMATCH:
        logger.error(verdict.message)
        raise ValueError(verdict.message)

    if verdict.status == STATUS_UNDETERMINED:
        logger.warning(verdict.message)
    else:
        logger.info(verdict.message)

    return verdict
