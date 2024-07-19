"""
utils.py

This module contains utility functions for the report generation.
"""

import pandas as pd
import random
from typing import Tuple
import natsort

def convert_to_space_separated_string(array):
    """
    Converts a list to a space-separated string.

    Args:
        array (list): List of strings.

    Returns:
        str: Space-separated string.
    """
    try:
        import ast

        # Convert array to list and extract the string
        string_repr = array.tolist()[0]

        # Evaluate the string to convert it to an actual list
        list_repr = ast.literal_eval(string_repr)

        # Join the elements of the list into a space-separated string
        return " ".join(list_repr)
    except:
        return array

def split_text(text):
    """
    Helper function to split text on '-' characters.

    Args:
        text (str): The text to be split.

    Returns:
        str: The text with '-' characters replaced by '-\n'.
    """
    return text.replace("-", "-\n")

def get_target_outliers(df):
    """
    Identifies outliers in target coverage.

    Args:
        df (pd.DataFrame): DataFrame containing target coverage data.

    Returns:
        pd.DataFrame: DataFrame containing outliers.
    """
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=natsort.natsorted(df["chrom"].unique()), ordered=True
    )
    df = df.sort_values("chrom")
    Q1 = df["coverage"].quantile(0.25)
    Q3 = df["coverage"].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1 * IQR
    upper_bound = Q3 + 1 * IQR
    outliers = df[(df["coverage"] < lower_bound) | (df["coverage"] > upper_bound)]
    return outliers
