"""
flagging.py

Module for applying configurable flagging rules to a DataFrame output from a tool.
This module provides functionality to add a flag column based on a set of logical rules
defined in a JSON configuration. Rules use the validated comparator data schema from
``comparator_rules`` and are compiled against the DataFrame's columns before any row is
processed.

Example flagging rule:
{
    "Low_Depth": {
        "all": [
            {
                "left": {"column": "Depth_Score"},
                "operator": "lt",
                "right": {"literal": 0.2}
            }
        ]
    }
}
"""

from __future__ import annotations

import logging
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from typing import NoReturn

import pandas as pd

from vntyper.scripts.comparator_rules import CompiledRule, adapt_legacy_rule, evaluate_rule, validate_rule

logger = logging.getLogger(__name__)

_RESERVED_FLAG_NAMES = frozenset({"Not flagged", "Not applicable"})
KESTREL_FLAG_COLUMNS = frozenset(
    {
        "Motifs",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Sample",
        "Motif_sequence",
        "Del",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "ref_len",
        "alt_len",
        "Frame_Score",
        "is_frameshift",
        "direction",
        "frameshift_amount",
        "is_valid_frameshift",
        "Depth_Score",
        "Confidence",
        "depth_confidence_pass",
        "haplo_count",
        "alt_filter_pass",
        "motif_filter_pass",
        "Motif",
        "Motif_fasta",
        "POS_fasta",
    }
)
ADVNTR_FLAG_COLUMNS = frozenset({"VID", "Variant", "NumberOfSupportingReads", "MeanCoverage", "Pvalue", "RU", "POS"})
_LEGACY_FLAG_RULES: dict[str, dict[str, object]] = {
    "False_Positive_4bp_Insertion": {
        "(REF == 'C') and (ALT == 'CGGCA')": {
            "all": [
                {"left": {"column": "REF"}, "operator": "eq", "right": {"literal": "C"}},
                {"left": {"column": "ALT"}, "operator": "eq", "right": {"literal": "CGGCA"}},
            ]
        }
    },
    "Low_Depth_Conserved_Motifs": {
        "(Depth_Score < 0.4) and (Motif in ['1', '2', '3', '4', '6', '7', '8', '9'])": {
            "all": [
                {"left": {"column": "Depth_Score"}, "operator": "lt", "right": {"literal": 0.4}},
                {
                    "left": {"column": "Motif"},
                    "operator": "in",
                    "right": {"literal": ["1", "2", "3", "4", "6", "7", "8", "9"]},
                },
            ]
        }
    },
    "Low_Coverage": {
        "NumberOfSupportingReads < 10": {
            "all": [
                {
                    "left": {"column": "NumberOfSupportingReads"},
                    "operator": "lt",
                    "right": {"literal": 10},
                }
            ]
        }
    },
    "Repeat_Unit_7": {"RU == '7'": {"all": [{"left": {"column": "RU"}, "operator": "eq", "right": {"literal": "7"}}]}},
    "Polymorphic_Call": {
        "Variant in ['I10_2_A_LEN1', 'D8_2&D9_2&I9_2_A_LEN9', 'D2_2&I2_2_C_LEN5', "
        "'I39_2_A_LEN4', 'I52_2_A_LEN7', 'I45_2_A_LEN4', 'D45_2&I45_2_A_LEN2', "
        "'D14_2&I14_2_G_LEN14', 'D58_2&D59_2', 'I60_2_A_LEN10', 'I14_2_G_LEN16', "
        "'I18_2_T_LEN1', 'I21_2_G_LEN4', 'D29_2&I29_2_A_LEN2', 'D8_2&I8_2_A_LEN20', "
        "'D20_2&D21_2', 'D21_2&D22_2', 'I14_2_A_LEN1', 'I11_2_G_LEN1', 'I26_7_A_LEN25', "
        "'D17_2&D18_2&D19_2&D20_2&D21_2', 'I14_2_C_LEN4', 'I23_6_G_LEN1', 'I21_2_T_LEN1']": {
            "all": [
                {
                    "left": {"column": "Variant"},
                    "operator": "in",
                    "right": {
                        "literal": [
                            "I10_2_A_LEN1",
                            "D8_2&D9_2&I9_2_A_LEN9",
                            "D2_2&I2_2_C_LEN5",
                            "I39_2_A_LEN4",
                            "I52_2_A_LEN7",
                            "I45_2_A_LEN4",
                            "D45_2&I45_2_A_LEN2",
                            "D14_2&I14_2_G_LEN14",
                            "D58_2&D59_2",
                            "I60_2_A_LEN10",
                            "I14_2_G_LEN16",
                            "I18_2_T_LEN1",
                            "I21_2_G_LEN4",
                            "D29_2&I29_2_A_LEN2",
                            "D8_2&I8_2_A_LEN20",
                            "D20_2&D21_2",
                            "D21_2&D22_2",
                            "I14_2_A_LEN1",
                            "I11_2_G_LEN1",
                            "I26_7_A_LEN25",
                            "D17_2&D18_2&D19_2&D20_2&D21_2",
                            "I14_2_C_LEN4",
                            "I23_6_G_LEN1",
                            "I21_2_T_LEN1",
                        ]
                    },
                }
            ]
        }
    },
}


@dataclass(frozen=True)
class CompiledFlagRules:
    """Immutable flag names and comparator trees validated for one consumer schema."""

    rules: tuple[tuple[str, CompiledRule], ...]


def _invalid_flagging(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def compile_flag_rules(flag_rules: object, columns: Collection[str]) -> CompiledFlagRules:
    """Validate and compile a complete ordered flag-rule mapping.

    Args:
        flag_rules: Untrusted flag-name-to-rule mapping from configuration.
        columns: Columns available in the consumer DataFrame.

    Returns:
        Immutable ordered flag names and compiled rules.

    Raises:
        ValueError: If the rule mapping, a flag name, or any rule is invalid.
    """
    if not isinstance(flag_rules, Mapping):
        _invalid_flagging("flagging_rules must be a mapping of flag names to structured rules")

    compiled: list[tuple[str, CompiledRule]] = []
    for flag, configured in flag_rules.items():
        if not isinstance(flag, str) or not flag:
            _invalid_flagging("flagging_rules flag name must be a non-empty string")
        if "," in flag:
            _invalid_flagging(f"flagging_rules flag name {flag!r} must not contain a comma")
        if flag in _RESERVED_FLAG_NAMES:
            _invalid_flagging(f"flagging_rules flag name {flag!r} is reserved")
        context = f"flagging_rules.{flag}"
        migrated = adapt_legacy_rule(configured, exact_rules=_LEGACY_FLAG_RULES.get(flag, {}), context=context)
        compiled.append((flag, validate_rule(migrated, allowed_columns=columns, context=context)))
    return CompiledFlagRules(rules=tuple(compiled))


def add_flags(
    df: pd.DataFrame,
    flag_rules: object,
    duplicates_config: dict | None = None,
) -> pd.DataFrame:
    """
    Applies flagging rules to the DataFrame and adds a 'Flag' column with the matched flags.

    The complete rule set is validated against the DataFrame columns before the frame is
    copied or any row is processed. For each row, each compiled structured rule is then
    evaluated. If a rule evaluates to True, the flag name is added to that row's flag list.
    If multiple flags apply, they are concatenated with a comma. If no flag is applied, the
    'Flag' column will be set to 'Not flagged'.

    Additionally, if duplicates_config is provided, rows that appear to be duplicates
    (according to grouping and sorting logic) will be flagged as configured.

    Args:
        df (pd.DataFrame): The input DataFrame containing tool output.
        flag_rules: Either an untrusted mapping from flag names to structured comparator
            rules, or a :class:`CompiledFlagRules` validated by the consumer before any
            sample-dependent early return.
        duplicates_config (dict, optional): Configuration for marking potential duplicates.
            Example structure:
            {
                "enabled": true,
                "flag_name": "Potential_Duplicate",
                "group_by": ["REF", "ALT"],
                "sort_by": [
                    {"column": "Depth_Score", "ascending": false},
                    {"column": "Motif", "ascending": true},
                    {"column": "Position", "ascending": true}
                ]
            }

    Returns:
        pd.DataFrame: A copy of the input DataFrame with an added 'Flag' column.
    """
    compiled = flag_rules if isinstance(flag_rules, CompiledFlagRules) else compile_flag_rules(flag_rules, df.columns)

    # Create a copy to avoid modifying the original DataFrame
    df_copy = df.copy()
    logger.debug("Created a copy of the DataFrame for flag processing.")

    # Initialize a list to store flags for each row
    flags: list[list[str]] = [[] for _ in range(len(df_copy))]
    logger.debug("Initialized flags list for each row.")

    # Evaluate each flag rule
    for flag, compiled_rule in compiled.rules:
        context = f"flagging_rules.{flag}"
        logger.debug(f"Evaluating validated flag rule '{flag}'.")
        mask = df_copy.apply(
            lambda row, rule=compiled_rule, rule_context=context: evaluate_rule(
                rule, row.to_dict(), context=rule_context
            ),
            axis=1,
        )
        matching_count = mask.sum()
        logger.debug(f"Flag rule '{flag}' matched {matching_count} rows.")
        for i, condition_met in enumerate(mask):
            if condition_met:
                flags[i].append(flag)
                logger.debug(f"Row {i} meets condition for flag '{flag}'.")

    # Create the 'Flag' column as a comma-separated string of flags for each row,
    # or 'Not flagged' if no flags were applied.
    df_copy["Flag"] = [", ".join(flag_list) if flag_list else "Not flagged" for flag_list in flags]
    logger.debug("Added 'Flag' column to DataFrame with flag values.")

    # Mark potential duplicates if configured
    if duplicates_config and duplicates_config.get("enabled", False):
        logger.debug(f"Duplicates config detected: {duplicates_config}")
        group_cols = duplicates_config.get("group_by", ["REF", "ALT"])
        sort_info = duplicates_config.get("sort_by", [])
        if not sort_info:
            # Fallback if sort_by not provided
            sort_cols = ["Depth_Score"]
            sort_ascending = [False]
        else:
            sort_cols = [item["column"] for item in sort_info]
            sort_ascending = [item["ascending"] for item in sort_info]

        logger.debug(
            "Calling mark_potential_duplicates with "
            f"group_cols={group_cols}, sort_cols={sort_cols}, sort_ascending={sort_ascending}"
        )
        df_copy = mark_potential_duplicates(
            df_copy,
            group_cols=group_cols,
            sort_cols=sort_cols,
            sort_ascending=sort_ascending,
            duplicate_flag_name=duplicates_config.get("flag_name", "Potential_Duplicate"),
        )
    else:
        logger.debug("No duplicates_config or 'enabled' is False; skipping duplicate flagging.")

    return df_copy


def add_artifact_gate(df: pd.DataFrame, artifact_flags: Sequence[str]) -> pd.DataFrame:
    """Mark rows carrying a declared artifact flag as failing the final filter.

    A flag is an annotation; most flags mean "this call warrants scrutiny". A few name a
    known technical artifact instead - a row that is not a candidate variant at all - and
    those must not reach ``kestrel_result.tsv``, where a consumer would read them as a
    call. Before #174, an artifact-only sample produced one flagged row that
    ``report_config.json`` mapped to ``High_Precision_flagged``, which ``is_finding``
    treats as a positive.

    Which flags are artifacts is **configuration** (``kestrel_config.json``'s
    ``artifact_flags``), never a name written into this function. Emptying that list
    restores the pre-#174 behaviour without a code change.

    This stage marks; it does not filter. ``filter_final_dataframe`` drops the rows, and
    ``kestrel_pre_result.tsv`` keeps every one of them with ``flag_filter_pass=False`` so
    the evidence for a sample is never destroyed.

    Args:
        df (pd.DataFrame): The variant frame. A missing ``Flag`` column is legitimate -
            a negative run carries none - and yields an all-True gate.
        artifact_flags (Sequence[str]): Flag names that disqualify a row.

    Returns:
        pd.DataFrame: A copy with a boolean ``flag_filter_pass`` column, present on every
        row including when the frame is empty. ``filter_final_dataframe`` raises on a
        missing gate column, so unconditional presence is the contract.

    Note:
        ``Flag`` is a comma-joined string, so membership is tested per element rather than
        by substring: an artifact named ``X`` must not exclude a flag named ``XY``. Values
        that are not strings - ``None``, ``NaN``, ``pd.NA``, all accepted inputs today -
        carry no declared artifact and therefore pass. ``select_single_best_variant``
        still deprioritises them.
    """
    # A bare string is a valid `Sequence[str]`, and `set("Foo")` is a set of *characters* -
    # so writing a plain string instead of a one-element list in kestrel_config.json, which
    # is valid JSON and an easy slip, would silently degrade the gate to matching single
    # letters and let every artifact through. Refuse it loudly: the whole point of #174 is
    # that a known artifact must not be reported as a call, and a config typo must not
    # quietly undo that.
    if isinstance(artifact_flags, str):
        msg = (
            f"artifact_flags must be a list of flag names, not the string {artifact_flags!r}. "
            "A string would be iterated character by character, silently disabling the "
            "artifact gate. See issue #174."
        )
        logger.error(msg)
        raise ValueError(msg)

    result = df.copy()
    artifacts = {str(flag).strip() for flag in artifact_flags}

    # `Flag` is joined with ", " and split on ",", so a declared name containing a comma
    # could never match any element and would be silently inert.
    commas = sorted(flag for flag in artifacts if "," in flag)
    if commas:
        msg = (
            f"artifact_flags entries must not contain a comma: {commas}. The 'Flag' column is "
            "comma-joined, so such a name can never match and the gate would silently do "
            "nothing. See issue #174."
        )
        logger.error(msg)
        raise ValueError(msg)

    if "Flag" not in result.columns or not artifacts:
        logger.debug("No artifact flags to apply (or no 'Flag' column); every row passes the gate.")
        result["flag_filter_pass"] = pd.Series(True, index=result.index, dtype=bool)
        return result

    def _is_clean(value: object) -> bool:
        """Return True when this ``Flag`` value names no declared artifact.

        Args:
            value (object): One ``Flag`` cell. Not necessarily a string.

        Returns:
            bool: False only when a declared artifact is one of the comma-separated
            elements of a string value.
        """
        if not isinstance(value, str):
            return True
        return not artifacts.intersection(part.strip() for part in value.split(","))

    result["flag_filter_pass"] = result["Flag"].map(_is_clean).astype(bool)
    excluded = int((~result["flag_filter_pass"]).sum())
    logger.debug(f"Artifact gate applied with {sorted(artifacts)}; {excluded} row(s) marked for exclusion.")
    return result


def mark_potential_duplicates(
    df: pd.DataFrame,
    group_cols: list,
    sort_cols: list,
    sort_ascending: list,
    duplicate_flag_name: str = "Potential_Duplicate",
) -> pd.DataFrame:
    """
    Mark potential duplicates in the DataFrame by grouping and sorting.

    This function identifies duplicates based on a grouping of columns and a specified
    sorting order. The first row in each group (after sorting) is considered the primary
    record, and all subsequent rows in the same group are flagged as potential duplicates.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        group_cols (list): A list of column names to group by (e.g. ["REF", "ALT"]).
        sort_cols (list): A list of column names to sort by in order of priority.
        sort_ascending (list): A list of booleans indicating ascending/descending order
                               for each corresponding column in sort_cols.
        duplicate_flag_name (str, optional): The flag name to assign for duplicates.
                                             Defaults to "Potential_Duplicate".

    Returns:
        pd.DataFrame: A copy of the input DataFrame with duplicate rows flagged in the 'Flag' column.
    """
    df_copy = df.copy()
    logger.debug(
        f"Marking potential duplicates with group_cols={group_cols}, "
        f"sort_cols={sort_cols}, sort_ascending={sort_ascending}, "
        f"duplicate_flag_name={duplicate_flag_name}"
    )

    # Keep track of original index to restore ordering after grouping/sorting
    df_copy["__original_index"] = df_copy.index
    logger.debug(f"DataFrame shape before sorting: {df_copy.shape}")

    # Sort by the specified columns and order.
    # kind="stable" because sort_by is now a single key by default (#197): pandas'
    # default quicksort does not guarantee tied rows keep their input order, so which
    # tied row keeps the first position -- and therefore which row is left unflagged
    # -- would be arbitrary between runs. This is an implementation detail beyond the
    # #197 decision, which only concerns the sort *key*, not sort stability.
    df_copy.sort_values(by=sort_cols, ascending=sort_ascending, inplace=True, kind="stable")
    logger.debug(f"DataFrame shape after sorting: {df_copy.shape}")

    # Mark duplicates within each group (the first row in each group has count=0)
    df_copy["__dup_indicator"] = df_copy.groupby(group_cols).cumcount()
    df_copy["__is_duplicate"] = df_copy["__dup_indicator"] > 0

    # Log how many duplicates we have
    num_duplicates = df_copy["__is_duplicate"].sum()
    logger.debug(f"Number of duplicates found: {num_duplicates}")

    # Prepare a list for the new flags
    new_flags: list[list[str]] = []
    for _idx, row in df_copy.iterrows():
        if row["__is_duplicate"]:
            new_flags.append([duplicate_flag_name])
        else:
            new_flags.append([])

    logger.debug("Combining existing flags with new duplicate flags.")
    combined_flags = []
    for i, row_tuple in enumerate(df_copy.itertuples(index=False)):
        # itertuples returns a named tuple - access Flag attribute safely
        # Use getattr to avoid mypy errors with dynamic pandas named tuple attributes
        existing_flag = row_tuple.Flag
        dup_flag_list = new_flags[i]

        if existing_flag == "Not flagged":
            # If not flagged and is a duplicate, set the new flag
            combined_flags.append(", ".join(dup_flag_list) if dup_flag_list else "Not flagged")
        else:
            # If already flagged, append the new flag if it's a duplicate
            if dup_flag_list:
                combined_flags.append(existing_flag + ", " + ", ".join(dup_flag_list))
            else:
                combined_flags.append(existing_flag)

    df_copy["Flag"] = combined_flags

    # Restore original ordering
    logger.debug("Restoring original row order.")
    df_copy.sort_values(by="__original_index", inplace=True)
    df_copy.drop(columns=["__original_index", "__dup_indicator", "__is_duplicate"], inplace=True)

    logger.debug(f"Done marking potential duplicates. Final shape: {df_copy.shape}")
    return df_copy
