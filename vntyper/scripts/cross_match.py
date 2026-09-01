"""
vntyper/scripts/cross_match.py

This module compares variant calls from Kestrel and adVNTR outputs.
It performs a configurable allele match by subtracting the REF from ALT for insertions/duplications
(or vice versa for deletions) and comparing the resulting allele change along with the variant type.
The matching rule is defined in the main configuration (under the "cross_match" key) and can be customized.
If a match is found according to the configured logic, the variant is considered concordant.
The module returns a summary dictionary that includes both all the individual comparisons
and an overall flag indicating if at least one combination has matched.

This module is designed to be used as a helper within the pipeline.
It accepts already‑parsed results (e.g. from a pipeline summary) without re‑parsing TSV files.
"""

import csv
import logging
from collections.abc import Mapping
from typing import NoReturn

from vntyper.scripts.comparator_rules import adapt_legacy_rule, evaluate_rule, validate_rule

# These names are matched by exact string comparison against what pipeline.py
# records. A typo does not fail - it silently drops the cross-match section
# (AGENTS.md trap 5), so they are named, never spelled out.
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_KESTREL

logger = logging.getLogger(__name__)

_LEGACY_MATCH_LOGIC = (
    "Kestrel_Allele_Change == Advntr_Allele_Change and Kestrel_Variant_Type.lower() == Advntr_Variant_Type.lower()"
)

DEFAULT_MATCH_RULE: dict[str, object] = {
    "all": [
        {
            "left": {"column": "Kestrel_Allele_Change"},
            "operator": "eq",
            "right": {"column": "Advntr_Allele_Change"},
        },
        {
            "left": {"column": "Kestrel_Variant_Type"},
            "operator": "casefold_eq",
            "right": {"column": "Advntr_Variant_Type"},
        },
    ]
}

CROSS_MATCH_COLUMNS: frozenset[str] = frozenset(
    {
        "Kestrel_POS",
        "Kestrel_REF",
        "Kestrel_ALT",
        "Kestrel_Allele_Change",
        "Kestrel_Variant_Type",
        "Advntr_POS",
        "Advntr_REF",
        "Advntr_ALT",
        "Advntr_Allele_Change",
        "Advntr_Variant_Type",
    }
)
_EVIDENCE_DISPOSITIONS = frozenset({"admissible", "identity-insufficient"})


def _configuration_error(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _configured_match_rule(config: object) -> object:
    """Resolve the structured rule or the exact historical compatibility value."""
    if config is None:
        return DEFAULT_MATCH_RULE
    if not isinstance(config, Mapping):
        _configuration_error("cross_match configuration must be a mapping")
    if "cross_match" not in config:
        return DEFAULT_MATCH_RULE
    configured = config["cross_match"]
    if not isinstance(configured, Mapping):
        _configuration_error("cross_match configuration block must be a mapping")
    unsupported_keys = sorted(set(configured) - {"match_rule", "match_logic"}, key=repr)
    if unsupported_keys:
        rendered_keys = ", ".join(repr(key) for key in unsupported_keys)
        _configuration_error(f"cross_match configuration block contains unsupported keys: {rendered_keys}")
    has_rule = "match_rule" in configured
    has_logic = "match_logic" in configured
    if has_rule and has_logic:
        _configuration_error("cross_match configuration cannot contain both 'match_rule' and 'match_logic'")
    if has_rule:
        return configured["match_rule"]
    if has_logic:
        legacy_logic = configured["match_logic"]
        if not isinstance(legacy_logic, str):
            _configuration_error("cross_match.match_logic must contain the exact historical string")
        return adapt_legacy_rule(
            legacy_logic,
            exact_rules={_LEGACY_MATCH_LOGIC: DEFAULT_MATCH_RULE},
            context="cross_match.match_logic",
        )
    return DEFAULT_MATCH_RULE


def _advntr_evidence_disposition(record: Mapping[str, object]) -> str:
    """Read the closed additive disposition, defaulting only legacy rows."""
    value = record.get("Evidence_Disposition", "admissible")
    if value not in _EVIDENCE_DISPOSITIONS:
        message = f"adVNTR record has unsupported Evidence_Disposition: {value}"
        logger.error(message)
        raise ValueError(message)
    return str(value)


def determine_variant_type(ref, alt):
    """
    Determine the variant type based on the lengths of REF and ALT.

    Args:
        ref (str): Reference allele.
        alt (str): Alternate allele.

    Returns:
        str: "Insertion", "Deletion", or "Other".
    """
    ref = str(ref)
    alt = str(alt)
    if len(alt) > len(ref):
        return "Insertion"
    elif len(alt) < len(ref):
        return "Deletion"
    else:
        return "Other"


def compute_allele_change(ref, alt, variant_type):
    """
    Compute the allele change string based on variant type.

    For insertions/duplications, the allele change is ALT with REF removed from its beginning.
    For deletions, it is REF with ALT removed from its beginning.

    Args:
        ref (str): Reference allele.
        alt (str): Alternate allele.
        variant_type (str): Variant type ("Insertion" or "Deletion").

    Returns:
        str: The allele change.
    """
    ref = str(ref)
    alt = str(alt)
    if variant_type.lower() in ["insertion", "duplication"]:
        if alt.startswith(ref):
            return alt[len(ref) :]
        return alt
    elif variant_type.lower() == "deletion":
        if ref.startswith(alt):
            return ref[len(alt) :]
        return ref
    return ""


def cross_match_variants(kestrel_records, advntr_records, config=None):
    """
    Cross-match variants from Kestrel and adVNTR outputs using a configurable rule.

    For each combination of Kestrel and adVNTR record, the function computes the allele change
    and variant type (if not already set) and then evaluates the matching rule.
    The rule is obtained from ``cross_match.match_rule``. The exact expression
    shipped before issue #286 is accepted under ``match_logic`` for migration.

    Args:
        kestrel_records (list of dict): Kestrel genotyping records.
        advntr_records (list of dict): adVNTR genotyping records.
        config (dict, optional): Main configuration dictionary.
            If provided, the matching rule is read from config["cross_match"]["match_rule"].

    Returns:
        dict: A dictionary with keys:
            "matches" - list of individual comparison records,
            "overall_match" - the STRING "Yes" if at least one combination matched,
            else "No". Never a boolean; downstream code matches on the string.

    Raises:
        ValueError: If the configured rule is malformed or cannot be evaluated.
            Validation happens before any input record is modified.
    """
    configured_rule = _configured_match_rule(config)
    compiled_rule = validate_rule(
        configured_rule,
        allowed_columns=CROSS_MATCH_COLUMNS,
        context="cross_match.match_rule",
    )

    results = []
    overall = False
    advntr_dispositions = [_advntr_evidence_disposition(record) for record in advntr_records]

    # Precompute allele change for Kestrel records.
    for k in kestrel_records:
        k_variant = k.get("Variant", "").strip() or determine_variant_type(k.get("REF", ""), k.get("ALT", ""))
        k["Variant_Type"] = k_variant
        k["Allele_Change"] = compute_allele_change(k.get("REF", ""), k.get("ALT", ""), k_variant)

    # Precompute allele change for adVNTR records.
    for a in advntr_records:
        a_variant = determine_variant_type(a.get("REF", ""), a.get("ALT", ""))
        a["Variant_Type"] = a_variant
        a["Allele_Change"] = compute_allele_change(a.get("REF", ""), a.get("ALT", ""), a_variant)

    # Evaluate each combination.
    for k in kestrel_records:
        for a, evidence_disposition in zip(advntr_records, advntr_dispositions, strict=True):
            result = {
                "Kestrel_POS": k.get("POS", ""),
                "Kestrel_REF": k.get("REF", ""),
                "Kestrel_ALT": k.get("ALT", ""),
                "Kestrel_Allele_Change": k["Allele_Change"],
                "Kestrel_Variant_Type": k["Variant_Type"],
                "Advntr_POS": a.get("POS", ""),
                "Advntr_REF": a.get("REF", ""),
                "Advntr_ALT": a.get("ALT", ""),
                "Advntr_Allele_Change": a["Allele_Change"],
                "Advntr_Variant_Type": a["Variant_Type"],
            }
            structural_match = evaluate_rule(compiled_rule, result, context="cross_match.match_rule")
            match = structural_match and evidence_disposition == "admissible"
            result["Match"] = "Yes" if match else "No"
            if match:
                overall = True
            results.append(result)

    overall_match = "Yes" if overall else "No"
    return {"matches": results, "overall_match": overall_match}


def write_results_tsv(results, output_path):
    """
    Write cross-match results to a TSV file.

    Args:
        results (list of dict): List of individual cross-match records.
        output_path (str or Path): File path to write the TSV.
    """
    if not results:
        logger.info("No results to write.")
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def extract_results_from_pipeline_summary(summary):
    """
    Extract Kestrel and adVNTR genotyping results from a pipeline summary dictionary.

    This function searches the summary's steps for those named "Kestrel Genotyping" and
    "adVNTR Genotyping" and returns their parsed results.

    Args:
        summary (dict): Pipeline summary dictionary.

    Returns:
        tuple: (kestrel_records, advntr_records) where each is a list of dictionaries.
               Returns (None, None) if not found.
    """
    kestrel_records = None
    advntr_records = None
    for step in summary.get("steps", []):
        if step.get("step") == STEP_KESTREL:
            kestrel_records = step.get("parsed_result", {}).get("data", [])
        elif step.get("step") == STEP_ADVNTR:
            advntr_records = step.get("parsed_result", {}).get("data", [])
    return kestrel_records, advntr_records
