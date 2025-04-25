"""
Helper functions to sort and merge bedmethyl files.
"""

import pandas as pd
import csv
import logging
from typing import List, Dict

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
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
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
        concat_df["strand"] = (
            concat_df["strand"].astype(str).replace({"+": ".", "-": "."})
        )

        groupby_columns: List[str] = [
            "chrom",
            "start_pos",
            "end_pos",
            "mod",
            "strand",
            "start_pos2",
            "end_pos2",
            "colour",
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


def parquet_to_bed(parquet_file: str, output_bed: str) -> None:
    """
    Convert a parquet file containing methylation data to BED format.
    
    The function reads a parquet file with methylation data and converts it to
    standard BED format. The input parquet file should have the following columns:
    - chrom: chromosome name
    - chromStart: start position
    - chromEnd: end position
    - mod_code: modification code
    - score_bed: score
    - strand: strand (+ or -)
    - thickStart: thick start position
    - thickEnd: thick end position
    - color: color value
    - valid_cov: valid coverage
    - percent_modified: percentage modified
    - n_mod: number of modifications
    - n_canonical: number of canonical bases
    - n_othermod: number of other modifications
    - n_delete: number of deletions
    - n_fail: number of failed calls
    - n_diff: number of different modifications
    - n_nocall: number of no-calls

    Args:
        parquet_file (str): Path to the input parquet file
        output_bed (str): Path to the output BED file

    Returns:
        None: Writes the BED file to the specified output path
    """
    try:
        # Read the parquet file
        df = pd.read_parquet(parquet_file)
        
        # Ensure all required columns are present
        required_columns = [
            "chrom", "chromStart", "chromEnd", "mod_code", "score_bed",
            "strand", "thickStart", "thickEnd", "color", "valid_cov",
            "percent_modified", "n_mod", "n_canonical", "n_othermod",
            "n_delete", "n_fail", "n_diff", "n_nocall"
        ]
        
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        # Create BED format DataFrame with correct column order
        bed_df = df[required_columns].copy()
        
        # Convert numeric columns to appropriate types
        numeric_columns = ["chromStart", "chromEnd", "thickStart", "thickEnd", "valid_cov"]
        for col in numeric_columns:
            bed_df[col] = pd.to_numeric(bed_df[col], errors='coerce').astype(int)
        
        # Ensure strand is properly formatted
        bed_df['strand'] = bed_df['strand'].astype(str).map({'.': '+', '-': '-'})
        
        # Sort by chromosome and start position
        bed_df = bed_df.sort_values(['chrom', 'chromStart'])
        
        # Write to BED file
        bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
        logging.info(f"Successfully converted {parquet_file} to BED format: {output_bed}")
        
    except Exception as e:
        logging.error(f"Error converting parquet to BED: {e}")
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
