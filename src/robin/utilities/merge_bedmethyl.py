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


def load_minimal_modkit_data(parquet_path):
    """
    Load only essential bedmethyl columns for optimized storage and processing.

    This function loads only the minimal columns required by the classifiers:
    - chrom: Chromosome name (required by all classifiers)
    - chromStart: Start position (required by all classifiers)
    - percent_modified: Primary methylation data (required by all classifiers)
    - mod_code: Modification code (required by Sturgeon for filtering)
    - strand: Strand information (required for proper aggregation)

    The function handles both old (18-column) and new (8-column optimized) formats.

    Args:
        parquet_path (str): Path to the parquet file containing bedmethyl data

    Returns:
        pd.DataFrame: DataFrame with only essential columns, sorted by chrom and chromStart
    """
    import time

    for attempt in range(5):  # Retry up to 5 times
        try:
            merged_modkit_df = pd.read_parquet(parquet_path)
            logger.debug("Successfully read the Parquet file.")
            break
        except Exception as e:
            logger.debug(f"Attempt {attempt+1}: File not ready ({e}). Retrying...")
            time.sleep(10)
    else:
        logger.debug("Failed to read Parquet file after multiple attempts.")
        return None

    # Debug: Check what columns are actually in the parquet file
    logger.info(f"Parquet file columns: {list(merged_modkit_df.columns)}")
    logger.info(f"Parquet file shape: {merged_modkit_df.shape}")

    # Check if this is the new optimized format (8 columns) or old format (18 columns)
    is_optimized_format = (
        len(merged_modkit_df.columns) <= 10
    )  # 8 essential columns + potential extras

    if is_optimized_format:
        logger.info("Detected optimized format - using essential columns directly")
        # Optimized format already has the essential columns we need
        essential_columns = [
            "chrom",
            "chromStart",
            "percent_modified",
            "mod_code",
            "strand",
        ]

        # Verify all essential columns are present
        missing_columns = [
            col for col in essential_columns if col not in merged_modkit_df.columns
        ]
        if missing_columns:
            logger.error(
                f"Missing essential columns in optimized format: {missing_columns}"
            )
            raise KeyError(f"Missing required columns: {missing_columns}")

        # Select only the essential columns
        df = merged_modkit_df[essential_columns].copy()

    else:
        logger.info("Detected legacy format - extracting essential columns")
        # Legacy format - use the original column mapping logic
        # Define minimal columns needed by classifiers with possible variations
        column_mappings = {
            "chrom": ["chrom", "chr", "chromosome"],
            "chromStart": ["chromStart", "start", "start_pos", "pos"],
            "percent_modified": [
                "percent_modified",
                "methylation_percent",
                "mod_percent",
                "score",
            ],
            "mod_code": ["mod_code", "mod", "modification_code"],
            "strand": ["strand", "strand_info"],
        }

        # Map available columns to our expected names
        mapped_columns = {}
        for expected_name, possible_names in column_mappings.items():
            found = False
            for possible_name in possible_names:
                if possible_name in merged_modkit_df.columns:
                    mapped_columns[expected_name] = possible_name
                    found = True
                    break
            if not found:
                logger.error(
                    f"Could not find column for {expected_name}. Available columns: {list(merged_modkit_df.columns)}"
                )
                raise KeyError(f"Missing required column: {expected_name}")

        # Create a new DataFrame with the expected column names
        df = merged_modkit_df[
            [
                mapped_columns[col]
                for col in [
                    "chrom",
                    "chromStart",
                    "percent_modified",
                    "mod_code",
                    "strand",
                ]
            ]
        ].copy()
        df.columns = ["chrom", "chromStart", "percent_modified", "mod_code", "strand"]

    logger.info(
        f"Loaded minimal bedmethyl data with {len(df)} rows and {len(df.columns)} columns"
    )
    return df.sort_values(by=["chrom", "chromStart"])


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


def collapse_minimal_bedmethyl(concat_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse and aggregate minimal bedmethyl dataframe.

    This function aggregates the minimal bedmethyl data by genomic position
    and calculates the fraction of methylation. Optimized for the minimal
    column set used by classifiers.

    Args:
        concat_df (pd.DataFrame): The concatenated dataframe with minimal columns
                                 (chrom, chromStart, percent_modified, mod_code, strand)

    Returns:
        pd.DataFrame: The collapsed and aggregated dataframe with fraction calculation
    """
    try:
        # Hack strand for aggregation (treat all strands the same for methylation analysis)
        concat_df["strand"] = (
            concat_df["strand"].astype(str).replace({"+": ".", "-": "."})
        )

        # Group by genomic position and modification code
        groupby_columns: List[str] = [
            "chrom",
            "chromStart",
            "mod_code",
            "strand",
        ]

        # Aggregate percent_modified by taking the mean (weighted average)
        agg_funcs: Dict[str, str] = {
            "percent_modified": "mean",
        }

        grouped = concat_df.groupby(groupby_columns, as_index=False, observed=True)
        result_df = grouped.agg(agg_funcs).reset_index()

        # Rename for consistency with existing code
        result_df = result_df.rename(
            columns={"chromStart": "start_pos", "percent_modified": "fraction"}
        )

        # Add end_pos column (same as start_pos for single-base positions)
        result_df["end_pos"] = result_df["start_pos"] + 1

        # Sort by chromosome and position
        merged_df = result_df.sort_values(by=["chrom", "start_pos"])
        logger.info("Successfully collapsed and aggregated the minimal bedmethyl data.")
        return merged_df
    except Exception as e:
        logger.error(f"Error collapsing minimal bedmethyl data: {e}")
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
    Convert a parquet file to BED format.

    This function handles both old (18-column) and new (8-column optimized) formats.
    For optimized format, it reconstructs the missing columns before conversion.

    Expected columns in the parquet file:
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

        # Check if this is the optimized format (8 columns) or full format (18 columns)
        is_optimized_format = (
            len(df.columns) <= 10
        )  # 8 essential columns + potential extras

        if is_optimized_format:
            logging.info(
                "Detected optimized format - reconstructing full format for BED conversion"
            )
            # Reconstruct missing columns for BED format
            df["chromEnd"] = df["chromStart"] + 1
            df["score_bed"] = df["percent_modified"]
            df["thickStart"] = df["chromStart"]
            df["thickEnd"] = df["chromEnd"]
            df["color"] = "255,0,0"
            df["n_othermod"] = 0
            df["n_delete"] = 0
            df["n_fail"] = 0
            df["n_diff"] = 0
            df["n_nocall"] = 0

        # Ensure all required columns are present
        required_columns = [
            "chrom",
            "chromStart",
            "chromEnd",
            "mod_code",
            "score_bed",
            "strand",
            "thickStart",
            "thickEnd",
            "color",
            "valid_cov",
            "percent_modified",
            "n_mod",
            "n_canonical",
            "n_othermod",
            "n_delete",
            "n_fail",
            "n_diff",
            "n_nocall",
        ]

        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")

        # Create BED format DataFrame with correct column order
        bed_df = df[required_columns].copy()

        # Convert numeric columns to appropriate types
        numeric_columns = [
            "chromStart",
            "chromEnd",
            "thickStart",
            "thickEnd",
            "valid_cov",
        ]
        for col in numeric_columns:
            bed_df[col] = pd.to_numeric(bed_df[col], errors="coerce").astype(int)

        # Ensure strand is properly formatted
        bed_df["strand"] = bed_df["strand"].astype(str).map({".": "+", "-": "-"})

        # Sort by chromosome and start position
        bed_df = bed_df.sort_values(["chrom", "chromStart"])

        # Write to BED file
        bed_df.to_csv(output_bed, sep="\t", header=False, index=False)
        logging.info(
            f"Successfully converted {parquet_file} to BED format: {output_bed}"
        )

    except Exception as e:
        logging.error(f"Error converting parquet to BED: {e}")
        raise


def reconstruct_full_bedmethyl_data(minimal_df: pd.DataFrame) -> pd.DataFrame:
    """
    Reconstruct full bedmethyl data structure from minimal data for RandomForest compatibility.

    This function takes minimal bedmethyl data and reconstructs the full 18-column
    structure expected by the RandomForest R script. Missing columns are filled with
    reasonable defaults or calculated values.

    The function handles both old (18-column) and new (8-column optimized) formats.

    Args:
        minimal_df (pd.DataFrame): Minimal bedmethyl data with columns:
                                  chrom, chromStart, percent_modified, mod_code, strand

    Returns:
        pd.DataFrame: Full bedmethyl data with all 18 columns
    """
    try:
        # Check if this is already in full format
        if len(minimal_df.columns) >= 15:  # Full format has 18 columns
            logger.info("Data already in full format, returning as-is")
            return minimal_df

        # Create a copy of the minimal data
        full_df = minimal_df.copy()

        # Add missing columns with reasonable defaults
        full_df["chromEnd"] = full_df["chromStart"] + 1
        full_df["score_bed"] = full_df[
            "percent_modified"
        ]  # Use methylation percentage as score
        full_df["thickStart"] = full_df["chromStart"]
        full_df["thickEnd"] = full_df["chromEnd"]
        full_df["color"] = "0,0,0"  # Default black color

        # Calculate coverage and modification counts based on methylation percentage
        # Assume a reasonable coverage value and calculate modifications
        full_df["valid_cov"] = 10  # Default coverage of 10 reads
        full_df["n_mod"] = (
            (full_df["percent_modified"] * full_df["valid_cov"] / 100)
            .round()
            .astype(int)
        )
        full_df["n_canonical"] = full_df["valid_cov"] - full_df["n_mod"]

        # Set other modification counts to 0 (not used by RandomForest)
        full_df["n_othermod"] = 0
        full_df["n_delete"] = 0
        full_df["n_fail"] = 0
        full_df["n_diff"] = 0
        full_df["n_nocall"] = 0

        # Rename columns to match expected format
        full_df = full_df.rename(
            columns={
                "chromStart": "start_pos",
                "chromEnd": "end_pos",
                "mod_code": "mod",
                "thickStart": "start_pos2",
                "thickEnd": "end_pos2",
                "color": "colour",
                "n_canonical": "Ncanon",
                "n_delete": "Ndel",
                "n_diff": "Ndiff",
                "n_fail": "Nfail",
                "n_mod": "Nmod",
                "n_nocall": "Nnocall",
                "n_othermod": "Nother",
                "valid_cov": "Nvalid",
                "percent_modified": "score",
            }
        )

        # Ensure all required columns are present
        required_columns = [
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
            "Nmod",
            "Ncanon",
            "Nother",
            "Ndel",
            "Nfail",
            "Ndiff",
            "Nnocall",
        ]

        for col in required_columns:
            if col not in full_df.columns:
                full_df[col] = 0  # Default to 0 for missing columns

        logger.info(f"Reconstructed full bedmethyl data with {len(full_df)} rows")
        return full_df[required_columns]

    except Exception as e:
        logger.error(f"Error reconstructing full bedmethyl data: {e}")
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
