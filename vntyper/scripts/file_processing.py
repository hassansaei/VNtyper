# vntyper/scripts/file_processing.py

import importlib.resources as pkg_resources
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def _assert_kestrel_allele_contract(ref, alt, snv_length, source_path):
    """
    Reject a VCF record that breaks the allele-shape contract of the pinned Kestrel build.

    VNtyper pins Kestrel 1.0.1 (``vntyper/dependencies/kestrel/kestrel.jar``, vendored and
    invoked by ``kestrel_genotyping.run_kestrel``). That build anchors every indel on a single
    reference base: ``VariantInsertion.getVcfRef()`` and ``VariantDeletion.getVcfAlt()`` return
    ``Character.toString(char)`` on every one of their return paths, and ``VariantType`` has
    only ``SNP``, ``INSERTION`` and ``DELETION``. It can therefore emit only three record
    shapes -- 1-vs-1, 1-vs-N and N-vs-1 -- and never two multi-base alleles.

    Nothing else in VNtyper enforces that, so a jar swap would otherwise change reported
    variant labels in silence. This check names the version it defends so that whoever
    replaces the jar is told what they broke.

    Args:
        ref (str): Reference allele from the record.
        alt (str): Alternate allele from the record.
        snv_length (int): Length of a single-base allele, from config.json.
        source_path (str): Path of the VCF being read, for the error message.

    Raises:
        ValueError: If both alleles are longer than ``snv_length``.
    """
    if len(ref) <= snv_length or len(alt) <= snv_length:
        return

    msg = (
        f"Off-contract VCF record in {source_path}: REF={ref!r} and ALT={alt!r} are both longer "
        f"than {snv_length} base(s). VNtyper pins Kestrel 1.0.1 "
        "(vntyper/dependencies/kestrel/kestrel.jar), which anchors every indel on a single "
        "reference base and so emits only 1-vs-1, 1-vs-N or N-vs-1 records. Two multi-base "
        "alleles mean the Kestrel build has changed; re-verify the insertion/deletion split and "
        "the calibrated thresholds before relaxing this check."
    )
    logger.error(msg)
    raise ValueError(msg)


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
    except OSError as exc:
        return f"it could not be read ({exc})"

    if records_before_header:
        return (
            f"it has no #CHROM header line, so all {records_before_header} data line(s) it carries are "
            "discarded when the file is parsed"
        )
    return "it has no #CHROM header line"


def filter_vcf(input_path, output_path):
    """
    Filter a VCF file to extract indels (insertions and deletions) and write them to a new file.

    A row is kept when the reference allele (REF) and the alternate allele (ALT) differ in
    length. Equal-length rows are substitutions and are dropped. The test is the length
    *difference*, so an indel is recognised whichever allele carries the extra bases -- not
    only when one of them happens to be a single base.

    This function is the trust boundary for Kestrel's output: it is the sole consumer of the
    raw VCF Kestrel writes (see ``kestrel_genotyping.process_kestrel_output``). Every record is
    checked against the pinned build's allele-shape contract before it is classified; see
    :func:`_assert_kestrel_allele_contract`.

    Args:
        input_path (str): Path to the input VCF file.
        output_path (str): Path to the output VCF file containing only indels.

    Raises:
        ValueError: If a data line has fewer than five tab-separated fields, or if a record
            carries two multi-base alleles and so breaks the pinned Kestrel output contract.
    """
    with pkg_resources.open_text("vntyper", "config.json") as f:
        config_data = json.load(f)
    snv_length = config_data.get("file_processing", {}).get("snv_length", 1)

    with open(input_path) as vcf_file, open(output_path, "w") as indel_file:
        for line in vcf_file:
            if line.startswith(("##", "#CHROM")):
                indel_file.write(line)
            else:
                _, _, _, ref, alt, *_ = line.split("\t")
                _assert_kestrel_allele_contract(ref, alt, snv_length, input_path)
                if len(ref) != len(alt):
                    indel_file.write(line)


def filter_indel_vcf(indel_vcf, output_ins, output_del):
    """
    Further filter the indel VCF file to separate insertions and deletions into two separate files.

    Insertions and deletions are differentiated by comparing the length of the reference allele
    (REF) with that of the alternate allele (ALT):
        - An insertion occurs when the ALT allele is longer than the REF allele.
        - A deletion occurs when the REF allele is longer than or equal to the ALT allele.

    The routing is that length comparison and nothing else, so a multi-base-REF insertion such
    as ``REF="AC", ALT="ACG"`` lands in the insertion file. It previously failed an insertion
    test of ``len(ref) == 1 and len(alt) > 1`` and fell through an unconditional ``else`` into
    the deletion file. Equal-length rows stay on the deletion side, which preserves the
    documented ``>`` boundary above.

    That mis-routing was latent rather than live: the pinned Kestrel 1.0.1 cannot emit a record
    with two multi-base alleles, and :func:`filter_vcf` now rejects one outright, so no such row
    reaches this function in the pipeline. Had one arrived, the consequence would have been a
    wrong ``Variant`` label ("Deletion" for an insertion) in ``kestrel_result.tsv``, the cohort
    tables and the HTML report, and a corrupted ``Allele_Change`` in the adVNTR cross-match,
    where ``cross_match.compute_allele_change`` would return the REF ("AC") instead of the
    inserted bases ("G"). The frameshift verdict is unaffected either way:
    ``scoring.extract_frameshifts`` derives ``direction`` and ``frameshift_amount`` from the REF
    and ALT lengths themselves and never reads the ``Variant`` label or the originating file.

    No allele-shape guard is applied here. This function's input is :func:`filter_vcf`'s output,
    a derived file rather than Kestrel's, so asserting a Kestrel-version contract on it would
    misattribute the failure. The contract is checked once, at the boundary, in
    :func:`filter_vcf`.

    Args:
        indel_vcf (str): Path to the input VCF file containing indels.
        output_ins (str): Path to the output VCF file for insertions.
        output_del (str): Path to the output VCF file for deletions.

    Raises:
        ValueError: If a data line has fewer than five tab-separated fields.
    """
    with open(indel_vcf) as vcf_file, open(output_ins, "w") as insertion_file, open(output_del, "w") as deletion_file:
        for line in vcf_file:
            if line.startswith(("##", "#CHROM")):
                insertion_file.write(line)
                deletion_file.write(line)
            else:
                _, _, _, ref, alt, *_ = line.split("\t")

                if len(alt) > len(ref):
                    insertion_file.write(line)
                else:
                    deletion_file.write(line)
