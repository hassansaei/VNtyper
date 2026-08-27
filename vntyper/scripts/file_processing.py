# vntyper/scripts/file_processing.py

import importlib.resources as pkg_resources
import json
import logging
import math
import re

logger = logging.getLogger(__name__)

_SCALAR_DEPTH_PATTERN = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


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


def _assert_empty_allele_contract(ref, alt, source_path):
    """
    Reject a VCF record with an empty REF or ALT allele.

    Every VCF record has both alleles. A record missing one is truncated, and it is the one
    malformed shape this function's caller would otherwise **discard in silence**: the indel
    test is ``len(ref) != len(alt)``, and two empty alleles are equal, so such a record is
    classified as a substitution and dropped. A Kestrel VCF whose records were all truncated
    that way would produce two valid-header-but-empty derived files, pass every header check,
    and be reported as a confident ``Negative`` (#223).

    Args:
        ref (str): Reference allele from the record.
        alt (str): Alternate allele from the record.
        source_path (str): Path of the VCF being read, for the error message.

    Raises:
        ValueError: If either allele is empty.
    """
    if ref and alt:
        return

    msg = (
        f"Truncated VCF record in {source_path}: REF={ref!r} and ALT={alt!r}, and a record must "
        "carry both alleles. Two empty alleles compare equal, so this record would be classified "
        "as a substitution and dropped without a trace -- and a file of them would be reported as "
        "a negative genotype rather than as the malformed output it is."
    )
    logger.error(msg)
    raise ValueError(msg)


def _assert_kestrel_format_contract(
    format_field: str,
    sample_field: str,
    expected_format: str,
    source_path: str,
) -> None:
    """Reject a Kestrel sample field that the positional scorer cannot consume safely.

    Scoring splits the sample value into exactly three positional columns and does not
    retain or inspect ``FORMAT``. The configured contract may change the genotype label,
    but GDP and DP must remain the second and third fields respectively. Both depths must
    be finite, non-negative scalar numeric values so malformed data cannot silently alter
    confidence assignment.

    Args:
        format_field: The record's FORMAT column.
        sample_field: The record's single sample column.
        expected_format: Required configured FORMAT expectation.
        source_path: VCF path included in diagnostics.

    Raises:
        ValueError: If the configured layout is incompatible with positional scoring, the
            record's FORMAT or width differs, or GDP/DP is not finite, non-negative scalar
            numeric data.
    """
    expected_fields = expected_format.split(":")
    if len(expected_fields) != 3 or expected_fields[1:] != ["GDP", "DP"]:
        msg = (
            f"Off-contract Kestrel FORMAT configuration for {source_path}: "
            f"file_processing.kestrel_format is {expected_format!r}, but the positional scorer "
            "requires exactly three fields with GDP and DP in positions 2 and 3. Re-verify and "
            "update the scorer before configuring another depth layout."
        )
        logger.error(msg)
        raise ValueError(msg)

    if format_field != expected_format:
        msg = (
            f"Off-contract VCF record in {source_path}: FORMAT is {format_field!r}, but "
            f"config.json file_processing.kestrel_format requires {expected_format!r}. A changed "
            "FORMAT can alter the positional depth mapping, so the run stops rather than mis-scoring."
        )
        logger.error(msg)
        raise ValueError(msg)

    sample_fields = sample_field.split(":")
    if len(sample_fields) != len(expected_fields):
        msg = (
            f"Off-contract VCF record in {source_path}: sample field {sample_field!r} has "
            f"{len(sample_fields)} colon-separated values, but FORMAT {expected_format!r} declares "
            f"{len(expected_fields)}. The positional depth mapping cannot be trusted."
        )
        logger.error(msg)
        raise ValueError(msg)

    for field_name, value in zip(expected_fields[1:], sample_fields[1:], strict=True):
        try:
            numeric_value = float(value)
        except ValueError:
            numeric_value = math.nan
        if _SCALAR_DEPTH_PATTERN.fullmatch(value) is None or not math.isfinite(numeric_value):
            msg = (
                f"Off-contract VCF record in {source_path}: {field_name}={value!r} in sample field "
                f"{sample_field!r} is not a finite scalar numeric depth. Malformed depths can coerce "
                "to zero during scoring, so the run stops rather than reporting a silent Negative."
            )
            logger.error(msg)
            raise ValueError(msg)
        if numeric_value < 0:
            msg = (
                f"Off-contract VCF record in {source_path}: {field_name}={value!r} in sample field "
                f"{sample_field!r} cannot be negative because depths are read counts. Negative depths "
                "can alter confidence assignment, so the run stops rather than reporting a silent genotype."
            )
            logger.error(msg)
            raise ValueError(msg)


def filter_vcf(input_path, output_path):
    """
    Filter a VCF file to extract indels (insertions and deletions) and write them to a new file.

    A row is kept when the reference allele (REF) and the alternate allele (ALT) differ in
    length. Equal-length rows are substitutions and are dropped. The test is the length
    *difference*, so an indel is recognised whichever allele carries the extra bases -- not
    only when one of them happens to be a single base.

    This function is the trust boundary for Kestrel's output: it is the sole consumer of the
    raw VCF Kestrel writes (see ``kestrel_genotyping.process_kestrel_output``). Every record is
    checked against the pinned build's allele-shape and FORMAT/sample contracts before it
    is classified; see :func:`_assert_kestrel_allele_contract` and
    :func:`_assert_kestrel_format_contract`.

    Args:
        input_path (str): Path to the input VCF file.
        output_path (str): Path to the output VCF file containing only indels.

    Raises:
        ValueError: If a record does not have exactly ten tab-separated fields, carries an
            empty allele or two multi-base alleles, or violates the configured FORMAT/sample
            depth contract.
    """
    with pkg_resources.open_text("vntyper", "config.json") as f:
        config_data = json.load(f)
    snv_length = config_data.get("file_processing", {}).get("snv_length", 1)
    kestrel_format = config_data["file_processing"]["kestrel_format"]

    with open(input_path) as vcf_file, open(output_path, "w") as indel_file:
        for line in vcf_file:
            if line.startswith(("##", "#CHROM")):
                indel_file.write(line)
            else:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 10:
                    msg = (
                        f"Off-contract VCF record in {input_path}: found {len(fields)} tab-separated "
                        "fields, but the pinned Kestrel output has exactly 10 (eight mandatory VCF "
                        "columns, FORMAT, and one sample column). The run stops rather than ignoring "
                        "or mis-scoring sample data."
                    )
                    logger.error(msg)
                    raise ValueError(msg)
                ref, alt = fields[3], fields[4]
                _assert_empty_allele_contract(ref, alt, input_path)
                _assert_kestrel_allele_contract(ref, alt, snv_length, input_path)
                _assert_kestrel_format_contract(fields[8], fields[9], kestrel_format, input_path)
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
