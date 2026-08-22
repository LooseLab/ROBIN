"""
Helper functions to sort and merge bedmethyl files. Requires Python 3.12+.
"""

from __future__ import annotations

import sys

if sys.version_info < (3, 12):
    raise RuntimeError("robin merge_bedmethyl utilities require Python 3.12 or newer")

import csv
import gc
import logging
import os
from copy import deepcopy
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# Sturgeon-related imports (must be installed)
from sturgeon.utils import read_probes_file

# from sturgeon.utils import read_probes_file, map_methyl_calls_to_probes_chr

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


def load_modkit_data(parquet_path):
    """
    Load minimal bedmethyl data for RandomForest analysis.

    This function loads only the essential columns needed for RandomForest classification:
    - chrom: Chromosome name
    - chromStart: Start position
    - percent_modified: Primary methylation data
    - mod_code: Modification code
    - strand: Strand information

    Args:
        parquet_path (str): Path to the parquet file containing bedmethyl data

    Returns:
        pd.DataFrame: DataFrame with minimal columns, sorted by chrom and chromStart
    """
    if load_minimal_modkit_data is None:
        raise ImportError(
            "load_minimal_modkit_data function not available from robin package"
        )
    return load_minimal_modkit_data(parquet_path)


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


def modkit_pileup_file_to_bed(
    input_data,  # Accepts either a DataFrame or file path
    output_file: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code: str = "m",
) -> pd.DataFrame:
    """Processes a modkit pileup file or DataFrame and maps methylation data to probes."""

    # Initialize variables to None for cleanup
    modkit_df = None
    probes_df = None
    probes_methyl_df = None
    calls_per_probe = None
    temp_file = None
    result_df = None

    try:
        # Check if input_data is a file path or a DataFrame
        if isinstance(input_data, str):
            # Read from file
            modkit_df = pd.read_csv(input_data, delim_whitespace=True, header=None)
        elif isinstance(input_data, pd.DataFrame):
            # Use DataFrame directly
            modkit_df = input_data.copy()
        else:
            raise ValueError(
                "input_data must be either a file path (str) or a DataFrame"
            )

        # For minimal data, we expect only the essential columns
        # Support various optimized formats: 5, 8, 9, or 10 columns
        if isinstance(input_data, pd.DataFrame) and len(modkit_df.columns) <= 10:
            # Minimal/optimized data format (5-10 columns)
            # Columns should already be named correctly
            expected_columns = [
                "chrom",
                "chromStart",
                "percent_modified",
                "mod_code",
                "strand",
            ]
            # Check if all expected columns exist in the DataFrame
            has_expected_cols = all(
                col in modkit_df.columns for col in expected_columns
            )

            if has_expected_cols:
                # Data has the essential columns, just filter and select them
                modkit_df = modkit_df[expected_columns].copy()
                modkit_df = modkit_df[modkit_df["mod_code"] == fivemc_code]
            elif len(modkit_df.columns) == 5:
                # Exactly 5 unnamed columns - assign names
                modkit_df.columns = expected_columns
                modkit_df = modkit_df[modkit_df["mod_code"] == fivemc_code]
            else:
                # Optimized format with more columns - try to extract the essential ones
                modkit_df = modkit_df[expected_columns].copy()
                modkit_df = modkit_df[modkit_df["mod_code"] == fivemc_code]
        else:
            # Full data format - use original logic
            # Define column names
            column_names = [
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

            # Validate number of columns
            if modkit_df.shape[1] != len(column_names):
                raise AssertionError(
                    f"Invalid modkit pileup file. Expected {len(column_names)} columns, got {modkit_df.shape[1]}."
                )

            # Assign column names
            modkit_df.columns = column_names

            # Filter by modification code
            modkit_df = modkit_df[modkit_df["mod_code"] == fivemc_code]

            # Drop unnecessary columns
            modkit_df.drop(
                columns=[
                    "mod_code",
                    "thickStart",
                    "thickEnd",
                    "color",
                    "valid_cov",
                    "n_mod",
                    "n_canonical",
                    "n_othermod",
                    "n_delete",
                    "n_fail",
                    "n_diff",
                    "n_nocall",
                ],
                inplace=True,
            )

        # Rename and normalize score column
        modkit_df = modkit_df.rename(
            columns={
                "chrom": "chr",
                "chromStart": "reference_pos",
                "percent_modified": "score",
            }
        )
        modkit_df["score"] /= 100  # Convert from percentage to decimal fraction

        # Remove invalid positions
        modkit_df = modkit_df[
            (modkit_df["reference_pos"] != -1) & (modkit_df["chr"] != ".")
        ]

        # Load probes file
        # The probes file has a header and uses whitespace (spaces) as delimiter
        probes_df = pd.read_csv(probes_file, sep=r"\s+", header=0)

        # Rename columns to expected names if needed
        if "ID_REF" in probes_df.columns:
            probes_df = probes_df.rename(columns={"ID_REF": "probe_name"})

        # Keep only the columns we need
        probes_df = probes_df[["chr", "start", "end", "probe_name"]].copy()

        # Ensure chromosome names match
        probes_df["chr"] = probes_df["chr"].astype(str).str.removeprefix("chr")
        modkit_df["chr"] = modkit_df["chr"].astype(str).str.removeprefix("chr")

        # Get unique chromosomes
        chromosomes = np.unique(probes_df["chr"].astype(str))

        # Copy probes data for methylation processing
        probes_methyl_df = deepcopy(probes_df)

        # Initialize methylation count columns
        probes_methyl_df["methylation_calls"] = 0
        probes_methyl_df["unmethylation_calls"] = 0
        probes_methyl_df["total_calls"] = 0

        # Process each chromosome
        calls_per_probe = []
        for chrom in chromosomes:
            chrom_str = str(chrom)  # Ensure correct format

            # Efficient filtering
            probe_mask = probes_methyl_df["chr"] == chrom_str
            methyl_mask = modkit_df["chr"] == chrom_str

            if probe_mask.sum() == 0 or methyl_mask.sum() == 0:
                continue  # Skip if no relevant data

            calls_per_probe_chr = map_methyl_calls_to_probes_chr(
                probes_df=probes_methyl_df.loc[probe_mask].copy(),
                methyl_calls_per_read=modkit_df.loc[methyl_mask].copy(),
                margin=margin,
                neg_threshold=neg_threshold,
                pos_threshold=pos_threshold,
            )

            # Rename 'probe_name' to 'probe_id' for consistency
            if "probe_name" in calls_per_probe_chr.columns:
                calls_per_probe_chr = calls_per_probe_chr.rename(
                    columns={"probe_name": "probe_id"}
                )

            calls_per_probe.append(calls_per_probe_chr)

            calls = calls_per_probe_chr["total_calls"].sum()
            logging.debug(
                f"Found {calls} methylation array sites on chromosome {chrom_str}"
            )

        # Ensure we have data before concatenation
        if not calls_per_probe:
            logging.warning("No methylation data found matching the probes.")
            return pd.DataFrame()

        # Merge all results efficiently
        calls_per_probe = pd.concat(calls_per_probe, ignore_index=True)

        # Save intermediate output
        temp_file = output_file + ".tmp"
        calls_per_probe.to_csv(temp_file, header=True, index=False, sep="\t")

        # Rename columns for the final output format
        calls_per_probe.rename(
            columns={
                "chr": "chrom",
                "start": "chromStart",
                "end": "chromEnd",
                "ID_REF": "probe_id",
                "methylation_calls": "methylation_call",
            },
            inplace=True,
        )

        # Filter out rows with zero total_calls
        calls_per_probe = calls_per_probe[calls_per_probe["total_calls"] > 0]

        # Select final output columns
        calls_per_probe = calls_per_probe[
            ["chrom", "chromStart", "chromEnd", "methylation_call", "probe_id"]
        ]

        # Save final processed file
        calls_per_probe.to_csv(output_file, header=True, index=False, sep="\t")

        # Store result for return
        result_df = calls_per_probe.copy()

        # Rename 'probe_name' to 'probe_id' for Sturgeon compatibility
        if "probe_name" in result_df.columns:
            result_df = result_df.rename(columns={"probe_name": "probe_id"})

    finally:
        # Clean up large DataFrames that are no longer needed
        if modkit_df is not None:
            del modkit_df
        if probes_df is not None:
            del probes_df
        if probes_methyl_df is not None:
            del probes_methyl_df
        if calls_per_probe is not None:
            del calls_per_probe

        # Clean up temporary file
        if temp_file is not None:
            try:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
                    logging.debug(f"Cleaned up temporary file: {temp_file}")
            except OSError as e:
                logging.warning(f"Failed to clean up temporary file {temp_file}: {e}")

        # Force garbage collection
        gc.collect()

    return result_df


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


def map_methyl_calls_to_probes_chr(
    probes_df: pd.DataFrame,
    methyl_calls_per_read: pd.DataFrame,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
) -> pd.DataFrame:
    """Maps calls per read to probe locations in a chromosome using NumPy for performance."""

    # Convert Pandas DataFrames to NumPy arrays for performance
    probes_start = probes_df["start"].to_numpy()
    methyl_pos = methyl_calls_per_read["reference_pos"].to_numpy()
    scores = methyl_calls_per_read["score"].to_numpy()

    # Define search ranges
    starts = probes_start - margin
    ends = starts + 2 * margin + 1

    # Vectorized binary search
    s = np.searchsorted(methyl_pos, starts, side="left")
    n = np.searchsorted(methyl_pos, ends, side="right")

    # Filter where matches exist
    valid_idx = s != n
    s, n, valid_idx = s[valid_idx], n[valid_idx], np.nonzero(valid_idx)[0]

    # Initialize call counters
    methylation_calls = np.zeros(len(probes_df), dtype=int)
    unmethylation_calls = np.zeros(len(probes_df), dtype=int)

    # Vectorized processing
    for idx, (ss, nn) in enumerate(zip(s, n)):
        current_scores = scores[ss:nn]
        bin_scores = np.zeros_like(current_scores)
        bin_scores[current_scores > pos_threshold] = 1
        bin_scores[current_scores < neg_threshold] = -1

        if len(bin_scores[bin_scores != 0]) > 0:
            final_score = int(np.median(bin_scores[bin_scores != 0]))
            if final_score == 1:
                methylation_calls[valid_idx[idx]] += 1
            elif final_score == -1:
                unmethylation_calls[valid_idx[idx]] += 1

    # Assign back to DataFrame
    probes_df["methylation_calls"] = methylation_calls
    probes_df["unmethylation_calls"] = unmethylation_calls
    probes_df["total_calls"] = methylation_calls + unmethylation_calls

    return probes_df
