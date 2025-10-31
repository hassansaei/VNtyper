#!/usr/bin/env python3
"""
flagging.py

Module for applying configurable flagging rules to a DataFrame output from a tool.
This module provides functionality to add a flag column based on a set of logical rules
defined in a JSON configuration. The flagging rules are specified as key-value pairs,
where the key is the flag name to be added if the condition is met, and the value is a
Python logical expression that is evaluated for each row. The condition is evaluated in a
context that includes the row's column values as variables and a helper function
`regex_match` for regex matching. Standard Python "in" operations are also supported.

Example flagging rule:
{
    "High_Depth": "Depth_Score >= 0.2 and regex_match('^D', Motif)",
    "Low_Depth": "Depth_Score < 0.2 or Motif in ['X', 'Y']"
}
"""

import logging
import re
from typing import Optional

import pandas as pd


def regex_match(pattern, value):
    """
    Helper function to perform regex matching.

    Args:
        pattern (str): Regular expression pattern.
        value (any): The value to match against (converted to string).

    Returns:
        bool: True if the pattern matches the string representation of value, else False.
    """
    try:
        return bool(re.search(pattern, str(value)))
    except Exception as e:
        logging.error(f"Error in regex_match with pattern {pattern} and value {value}: {e}")
        return False


def evaluate_condition(row, condition):
    """
    Evaluates a condition in the context of a DataFrame row.

    The local namespace is populated with the row's data and the helper function regex_match.
    If a column referenced in the condition is missing, a warning is logged and the condition
    evaluates to False for that row.

    Args:
        row (pd.Series): A row from the DataFrame.
        condition (str): The condition to evaluate.

    Returns:
        bool: The boolean result of evaluating the condition for this row.
    """
    local_vars = dict(row)
    local_vars["regex_match"] = regex_match
    try:
        result = eval(condition, {"__builtins__": {}}, local_vars)
        return bool(result)
    except NameError as ne:
        logging.warning(f"NameError while evaluating condition '{condition}' with row {row.to_dict()}: {ne}")
        return False
    except Exception as e:
        logging.error(f"Error evaluating condition '{condition}' with row {row.to_dict()}: {e}")
        return False


def add_flags(df: pd.DataFrame, flag_rules: dict, duplicates_config: Optional[dict] = None) -> pd.DataFrame:
    """
    Applies flagging rules to the DataFrame and adds a 'Flag' column with the matched flags.

    For each row in the DataFrame, each flag rule is evaluated by applying the condition
    to the row's values (and using the helper function regex_match if needed). If a rule's
    condition evaluates to True for a given row, the flag name is added to that row's flag list.
    If multiple flags apply, they are concatenated with a comma. If no flag is applied, the
    'Flag' column will be set to 'Not flagged'.

    Additionally, if duplicates_config is provided, rows that appear to be duplicates
    (according to grouping and sorting logic) will be flagged as configured.

    Args:
        df (pd.DataFrame): The input DataFrame containing tool output.
        flag_rules (dict): A dictionary where keys are flag names and values are
                           Python logical expressions (as strings) to be evaluated.
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
    # Create a copy to avoid modifying the original DataFrame
    df_copy = df.copy()
    logging.debug("Created a copy of the DataFrame for flag processing.")

    # Initialize a list to store flags for each row
    flags: list[list[str]] = [[] for _ in range(len(df_copy))]
    logging.debug("Initialized flags list for each row.")

    # Evaluate each flag rule
    for flag, condition in flag_rules.items():
        logging.debug(f"Evaluating flag rule '{flag}': {condition}")
        mask = df_copy.apply(lambda row, cond=condition: evaluate_condition(row, cond), axis=1)
        matching_count = mask.sum()
        logging.debug(f"Flag rule '{flag}' matched {matching_count} rows.")
        for i, condition_met in enumerate(mask):
            if condition_met:
                flags[i].append(flag)
                logging.debug(f"Row {i} meets condition for flag '{flag}'.")

    # Create the 'Flag' column as a comma-separated string of flags for each row,
    # or 'Not flagged' if no flags were applied.
    df_copy["Flag"] = [", ".join(flag_list) if flag_list else "Not flagged" for flag_list in flags]
    logging.debug("Added 'Flag' column to DataFrame with flag values.")

    # Mark potential duplicates if configured
    if duplicates_config and duplicates_config.get("enabled", False):
        logging.debug(f"Duplicates config detected: {duplicates_config}")
        group_cols = duplicates_config.get("group_by", ["REF", "ALT"])
        sort_info = duplicates_config.get("sort_by", [])
        if not sort_info:
            # Fallback if sort_by not provided
            sort_cols = ["Depth_Score"]
            sort_ascending = [False]
        else:
            sort_cols = [item["column"] for item in sort_info]
            sort_ascending = [item["ascending"] for item in sort_info]

        logging.debug(
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
        logging.debug("No duplicates_config or 'enabled' is False; skipping duplicate flagging.")

    return df_copy


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
    logging.debug(
        f"Marking potential duplicates with group_cols={group_cols}, "
        f"sort_cols={sort_cols}, sort_ascending={sort_ascending}, "
        f"duplicate_flag_name={duplicate_flag_name}"
    )

    # Keep track of original index to restore ordering after grouping/sorting
    df_copy["__original_index"] = df_copy.index
    logging.debug(f"DataFrame shape before sorting: {df_copy.shape}")

    # Sort by the specified columns and order
    df_copy.sort_values(by=sort_cols, ascending=sort_ascending, inplace=True)
    logging.debug(f"DataFrame shape after sorting: {df_copy.shape}")

    # Mark duplicates within each group (the first row in each group has count=0)
    df_copy["__dup_indicator"] = df_copy.groupby(group_cols).cumcount()
    df_copy["__is_duplicate"] = df_copy["__dup_indicator"] > 0

    # Log how many duplicates we have
    num_duplicates = df_copy["__is_duplicate"].sum()
    logging.debug(f"Number of duplicates found: {num_duplicates}")

    # Prepare a list for the new flags
    new_flags = []
    for _idx, row in df_copy.iterrows():
        if row["__is_duplicate"]:
            new_flags.append([duplicate_flag_name])
        else:
            new_flags.append([])

    logging.debug("Combining existing flags with new duplicate flags.")
    combined_flags = []
    for i, row_tuple in enumerate(df_copy.itertuples(index=False)):
        # itertuples returns a named tuple - access Flag attribute safely
        # Use getattr to avoid mypy errors with dynamic pandas named tuple attributes
        existing_flag = getattr(row_tuple, "Flag")  # noqa: B009
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
    logging.debug("Restoring original row order.")
    df_copy.sort_values(by="__original_index", inplace=True)
    df_copy.drop(columns=["__original_index", "__dup_indicator", "__is_duplicate"], inplace=True)

    logging.debug(f"Done marking potential duplicates. Final shape: {df_copy.shape}")
    return df_copy
