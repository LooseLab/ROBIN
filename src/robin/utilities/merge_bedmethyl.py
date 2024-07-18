"""
Helper functions to sort and merge bedmethyl files.
"""

import pandas as pd
import csv
import logging
from typing import List, Dict, Any

# Configure logging
logger = logging.getLogger(__name__)

def configure_logging(level: int = logging.INFO) -> None:
    """
    Configure the logging for this module.

    Args:
        level (int): The logging level to set. Defaults to logging.INFO.
    """
    logger.setLevel(level)

    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

def save_bedmethyl(result_df: pd.DataFrame, output_file: str) -> None:
    """
    Save a bedmethyl dataframe to a file.

    This code is only compatible with modkit >= 0.3.0

    Args:
        result_df (pd.DataFrame): The dataframe containing bedmethyl data.
        output_file (str): The path to the output file.

    Raises:
        IOError: If there's an error writing to the file.
    """
    try:
        result_df.to_csv(
            output_file,
            sep="\t",
            header=None,
            index=False,
            quoting=csv.QUOTE_NONNUMERIC,
            quotechar='"',
            escapechar="\\",
        )
        logger.info(f"Saved bedmethyl data to {output_file}")
    except IOError as e:
        logger.error(f"Error saving bedmethyl data to {output_file}: {e}")
        raise

def collapse_bedmethyl(concat_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse and aggregate bedmethyl dataframe.

    This function aggregates the bedmethyl data by specified columns and calculates the fraction.

    Args:
        concat_df (pd.DataFrame): The concatenated dataframe to collapse.

    Returns:
        pd.DataFrame: The collapsed and aggregated dataframe.
    """
    try:
        # Hack strand for aggregation
        concat_df["strand"] = concat_df["strand"].astype(str).replace({"+": ".", "-": "."})

        groupby_columns: List[str] = [
            "chrom", "start_pos", "end_pos", "mod", "strand",
            "start_pos2", "end_pos2", "colour"
        ]

        agg_funcs: Dict[str, str] = {
            "score": "sum",
            "Nvalid": "sum",
            "Nmod": "sum",
            "Ncanon": "sum",
            "Nother": "sum",
            "Ndel": "sum",
            "Nfail": "sum",
            "Ndiff": "sum",
            "Nnocall": "sum",
        }

        grouped = concat_df.groupby(groupby_columns, as_index=False, observed=True)
        result_df = grouped.agg(agg_funcs).reset_index()

        result_df["fraction"] = result_df["Nmod"] / result_df["Nvalid"] * 100

        column_order: List[str] = [
            "chrom", "start_pos", "end_pos", "mod", "score", "strand",
            "start_pos2", "end_pos2", "colour", "Nvalid", "fraction", "Nmod",
            "Ncanon", "Nother", "Ndel", "Nfail", "Ndiff", "Nnocall",
        ]

        merged_df = result_df[column_order].sort_values(by=["chrom", "start_pos"])
        logger.info("Successfully collapsed and aggregated the bedmethyl data.")
        return merged_df
    except Exception as e:
        logger.error(f"Error collapsing bedmethyl data: {e}")
        raise

def merge_bedmethyl(dfA: pd.DataFrame, dfB: pd.DataFrame) -> pd.DataFrame:
    """
    Merge two bedmethyl dataframes into a single dataframe.

    Given two bedmethyl dataframes, merge them into a single dataframe.

    Args:
        dfA (pd.DataFrame): The first bedmethyl dataframe.
        dfB (pd.DataFrame): The second bedmethyl dataframe.

    Returns:
        pd.DataFrame: The merged dataframe.

    Raises:
        ValueError: If the input dataframes are empty or have incompatible schemas.
    """
    try:
        if dfA.empty or dfB.empty:
            raise ValueError("One or both input dataframes are empty.")

        if not set(dfA.columns) == set(dfB.columns):
            raise ValueError("Input dataframes have incompatible schemas.")

        concat_df = pd.concat([dfA, dfB], ignore_index=True)
        merged_df = collapse_bedmethyl(concat_df)
        logger.info("Successfully merged two bedmethyl dataframes.")
        return merged_df
    except Exception as e:
        logger.error(f"Error merging bedmethyl dataframes: {e}")
        raise

if __name__ == "__main__":
    configure_logging(level=logging.DEBUG)

    # Example usage
    try:
        # Assuming you have two sample DataFrames dfA and dfB
        # dfA = pd.read_csv("sample_A.bedmethyl", sep="\t", header=None)
        # dfB = pd.read_csv("sample_B.bedmethyl", sep="\t", header=None)

        # merged_df = merge_bedmethyl(dfA, dfB)
        # save_bedmethyl(merged_df, "merged_output.bedmethyl")

        logger.info("Bedmethyl processing completed successfully.")
    except Exception as e:
        logger.error(f"An error occurred during bedmethyl processing: {e}")