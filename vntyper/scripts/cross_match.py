#!/usr/bin/env python3
"""
vntyper/scripts/cross_match.py

This module compares variant calls from Kestrel and adVNTR outputs.
It performs a configurable allele match by subtracting the REF from ALT for insertions/duplications
(or vice versa for deletions) and comparing the resulting allele change along with the variant type.
The matching logic is defined in the main configuration (under the "cross_match" key) and can be customized.
If a match is found according to the configured logic, the variant is considered concordant.
The module returns a summary dictionary that includes both all the individual comparisons
and an overall flag indicating if at least one combination has matched.

This module is designed to be used as a helper within the pipeline.
It accepts already‑parsed results (e.g. from a pipeline summary) without re‑parsing TSV files.
"""

import csv
import logging

DEFAULT_MATCH_LOGIC = (
    "Kestrel_Allele_Change == Advntr_Allele_Change and Kestrel_Variant_Type.lower() == Advntr_Variant_Type.lower()"
)


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
    Cross-match variants from Kestrel and adVNTR outputs using configurable logic.

    For each combination of Kestrel and adVNTR record, the function computes the allele change
    and variant type (if not already set) and then evaluates the matching logic.
    The matching logic is obtained from the configuration dictionary under the key
    "cross_match" -> "match_logic". If not provided, the default logic is used.

    Args:
        kestrel_records (list of dict): Kestrel genotyping records.
        advntr_records (list of dict): adVNTR genotyping records.
        config (dict, optional): Main configuration dictionary.
            If provided, the matching logic is read from config["cross_match"]["match_logic"].

    Returns:
        dict: A dictionary with keys:
            "matches" - list of individual comparison records,
            "overall_match" - "Yes" if at least one combination matched, else "No".
    """
    if config is not None:
        match_logic = config.get("cross_match", {}).get("match_logic", DEFAULT_MATCH_LOGIC)
    else:
        match_logic = DEFAULT_MATCH_LOGIC

    results = []
    overall = False

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
        for a in advntr_records:
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
            try:
                # Evaluate the matching condition in a restricted namespace.
                match = bool(eval(match_logic, {"__builtins__": {}}, result))
            except Exception as e:
                logging.error(f"Error evaluating match logic: {e}")
                match = False
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
        logging.info("No results to write.")
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
        if step.get("step") == "Kestrel Genotyping":
            kestrel_records = step.get("parsed_result", {}).get("data", [])
        elif step.get("step") == "adVNTR Genotyping":
            advntr_records = step.get("parsed_result", {}).get("data", [])
    return kestrel_records, advntr_records
