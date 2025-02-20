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
        logging.error(
            f"Error in regex_match with pattern {pattern} and value {value}: {e}"
        )
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
        logging.warning(
            f"NameError while evaluating condition '{condition}' with row {row.to_dict()}: {ne}"
        )
        return False
    except Exception as e:
        logging.error(
            f"Error evaluating condition '{condition}' with row {row.to_dict()}: {e}"
        )
        return False


def add_flags(df: pd.DataFrame, flag_rules: dict) -> pd.DataFrame:
    """
    Applies flagging rules to the DataFrame and adds a 'Flag' column with the matched flags.

    For each row in the DataFrame, each flag rule is evaluated by applying the condition
    to the row's values (and using the helper function regex_match if needed). If a rule's
    condition evaluates to True for a given row, the flag name is added to that row's flag list.
    If multiple flags apply, they are concatenated with a comma. If no flag is applied, the
    'Flag' column will be set to None.

    Args:
        df (pd.DataFrame): The input DataFrame containing tool output.
        flag_rules (dict): A dictionary where keys are flag names and values are
                           Python logical expressions (as strings) to be evaluated.

    Returns:
        pd.DataFrame: A copy of the input DataFrame with an added 'Flag' column.
    """
    # Create a copy to avoid modifying the original DataFrame
    df_copy = df.copy()
    logging.debug("Created a copy of the DataFrame for flag processing.")

    # Initialize a list to store flags for each row
    flags = [[] for _ in range(len(df_copy))]
    logging.debug("Initialized flags list for each row.")

    # Iterate over each flag rule and evaluate its condition for each row
    for flag, condition in flag_rules.items():
        logging.debug(f"Evaluating flag rule '{flag}': {condition}")
        # Evaluate condition row-wise using apply with our custom evaluator
        mask = df_copy.apply(lambda row: evaluate_condition(row, condition), axis=1)
        matching_count = mask.sum()
        logging.debug(f"Flag rule '{flag}' matched {matching_count} rows.")
        # Append the flag to rows where the condition is True
        for i, condition_met in enumerate(mask):
            if condition_met:
                flags[i].append(flag)
                logging.debug(f"Row {i} meets condition for flag '{flag}'.")

    # Create the 'Flag' column as a comma-separated string of flags for each row,
    # or None if no flags were applied.
    df_copy["Flag"] = [
        ", ".join(flag_list) if flag_list else "Not flagged" for flag_list in flags
    ]
    logging.debug("Added 'Flag' column to DataFrame with flag values.")
    return df_copy
