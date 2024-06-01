"""
Helper functions to sort and merge bedmethyl files.
"""

import pandas as pd
import csv
import os
import logging

# Configure logging
logger = logging.getLogger(__name__)


def save_bedmethyl(result_df: pd.DataFrame, output_file: str) -> None:
    """
    Save a bedmethyl dataframe to a file.

    Given a bedmethyl dataframe, save it to a file with specific formatting.

    Args:
        result_df (pd.DataFrame): The dataframe containing bedmethyl data.
        output_file (str): The path to the output file.
    """
    # Combine columns up to the 9th as tab-separated values
    tab_sep_cols = result_df.iloc[:, :9].astype(str).apply("\t".join, axis=1)

    # Combine remaining columns as space-separated values
    space_sep_cols = result_df.iloc[:, 9:].astype(str).apply(" ".join, axis=1)

    # Combine both sets of columns into a single DataFrame
    combined_df = pd.DataFrame({"tab_sep": tab_sep_cols, "space_sep": space_sep_cols})

    # Save the combined data to a file
    file_path = output_file

    combined_df.to_csv(
        file_path,
        sep="\t",
        header=None,
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        quotechar='"',
        escapechar='\\'
    )
    logger.info(f"Saved bedmethyl data to {file_path}")

    os.system(f"sed -i.bak 's/\t/\t/g' {file_path}")
    logger.info(f"Post-processed the file with sed: {file_path}")


def collapse_bedmethyl(concat_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse and aggregate bedmethyl dataframe.

    This function aggregates the bedmethyl data by specified columns and calculates the fraction.

    Args:
        concat_df (pd.DataFrame): The concatenated dataframe to collapse.

    Returns:
        pd.DataFrame: The collapsed and aggregated dataframe.
    """
    # Hack strand for aggregation
    concat_df["strand"] = concat_df["strand"].astype(str).replace({"+": ".", "-": "."})
    grouped = concat_df.groupby(
        [
            "chrom",
            "start_pos",
            "end_pos",
            "mod",
            "strand",
            "start_pos2",
            "end_pos2",
            "colour",
        ],
        as_index=False,
        observed=True,
    )
    agg_funcs = {
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

    result_df = grouped.agg(agg_funcs).reset_index()

    result_df["fraction"] = result_df["Nmod"] / result_df["Nvalid"] * 100

    column_order = [
        "chrom",
        "start_pos",
        "end_pos",
        "mod",
        "score",
        "strand",
        "start_pos2",
        "end_pos2",
        "colour",
        "Nvalid",
        "fraction",
        "Nmod",
        "Ncanon",
        "Nother",
        "Ndel",
        "Nfail",
        "Ndiff",
        "Nnocall",
    ]
    merged_df = result_df[column_order].sort_values(by=["chrom", "start_pos"])
    logger.info("Collapsed and aggregated the bedmethyl data.")
    return merged_df


def merge_bedmethyl(dfA: pd.DataFrame, dfB: pd.DataFrame) -> pd.DataFrame:
    """
    Merge two bedmethyl dataframes into a single dataframe.

    Given two bedmethyl dataframes, merge them into a single dataframe.

    Args:
        dfA (pd.DataFrame): The first bedmethyl dataframe.
        dfB (pd.DataFrame): The second bedmethyl dataframe.

    Returns:
        pd.DataFrame: The merged dataframe.
    """
    concat_df = pd.concat([dfA, dfB])
    merged_df = collapse_bedmethyl(concat_df)
    logger.info("Merged two bedmethyl dataframes.")
    return merged_df
