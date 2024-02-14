"""
Helper functions to sort and merge bedmethyl files.
"""

import pandas as pd
import csv
import os


def save_bedmethyl(result_df, output_file):
    """
    Save a bedmethyl dataframe to a file.
    Given a bedmethyl dataframe, save it to a file.
    """
    # Combine columns up to the 11th as tab-separated values
    tab_sep_cols = result_df.iloc[:, :9].astype(str).apply("    ".join, axis=1)

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
        quoting=csv.QUOTE_NONE,
        quotechar="",
    )
    os.system(f"sed -i.bak 's/    /\t/g' {file_path}")


def collapse_bedmethyl(concat_df):
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
    merged_df = result_df[column_order]
    # merged_df[['start_pos2', 'end_pos2','score','Nvalid','Nmod', 'Ncanon', 'Nother', 'Ndel', 'Nfail', 'Ndiff', 'Nnocall']] = merged_df[['start_pos2', 'end_pos2','score','Nvalid','Nmod', 'Ncanon', 'Nother', 'Ndel', 'Nfail', 'Ndiff', 'Nnocall']].astype(int)
    merged_df.sort_values(by=["chrom", "start_pos"], inplace=True)
    return merged_df


def merge_bedmethyl(dfA, dfB):
    """
    Merge bedmethyl files into a single dataframe.
    Given one or more bedmethyl files, merge them into a single dataframe.
    """

    concat_df = pd.concat([dfA, dfB])

    merged_df = collapse_bedmethyl(concat_df)

    return merged_df
