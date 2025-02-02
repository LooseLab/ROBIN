"""
A module for analyzing and visualizing Sturgeon methylation data in real-time.

This module provides functionality for processing BAM files to analyze methylation patterns
using the Sturgeon tool suite. It includes capabilities for:

- Real-time BAM file processing
- Methylation data extraction using modkit
- Sturgeon prediction and visualization
- Time series analysis of methylation patterns

The module integrates with the NiceGUI framework for interactive visualization
and uses temporary files for efficient data processing.

Dependencies
-----------
- pysam
- pandas
- nicegui
- sturgeon
- modkit

Notes
-----
The module requires proper configuration of input/output directories and
assumes the presence of necessary model files for Sturgeon predictions.
"""

from robin.subpages.base_analysis import BaseAnalysis
import os
import sys
import tempfile
import time
import shutil
import pandas as pd
from nicegui import ui, app, run
from robin import theme
import pysam
from robin import models
from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)
import click
from pathlib import Path
from typing import Optional, List, Tuple
import logging

import json
import zipfile
from copy import deepcopy


import numpy as np


import sturgeon

# Sturgeon-related imports (must be installed)
from sturgeon.utils import validate_model_file, get_model_path, read_probes_file
from sturgeon.prediction import load_model, predict_sample
from sturgeon.plot import plot_prediction

from robin import models

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


def load_modkit_data(parquet_path):
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

    # Keep only the original 18 columns
    return merged_modkit_df[column_names]


def modkit_pileup_file_to_bed(
    input_data,  # Accepts either a DataFrame or file path
    output_file: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code: str = "C",
) -> pd.DataFrame:
    """Processes a modkit pileup file or DataFrame and maps methylation data to probes."""

    # Check if input_data is a file path or a DataFrame
    if isinstance(input_data, str):
        # Read from file
        modkit_df = pd.read_csv(input_data, delim_whitespace=True, header=None)
    elif isinstance(input_data, pd.DataFrame):
        # Use DataFrame directly
        modkit_df = input_data.copy()
    else:
        raise ValueError("input_data must be either a file path (str) or a DataFrame")

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

    # Rename and normalize score column
    modkit_df = modkit_df.rename(
        columns={
            "chrom": "chr",
            "chromStart": "reference_pos",
            "percent_modified": "score",
        }
    )
    modkit_df["score"] /= 100  # Convert from percentage to decimal fraction

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

    # Remove invalid positions
    modkit_df = modkit_df[
        (modkit_df["reference_pos"] != -1) & (modkit_df["chr"] != ".")
    ]

    # Load probes file
    probes_df = read_probes_file(probes_file)

    # Ensure chromosome names match
    probes_df["chr"] = probes_df["chr"].astype(str)  # Make sure probes are strings
    modkit_df["chr"] = (
        modkit_df["chr"].astype(str).str.replace("^chr", "", regex=True)
    )  # Remove "chr" prefix

    # Print to verify
    # print("Normalized Chromosomes in probes:", np.unique(probes_df['chr']))
    # print("Normalized Chromosomes in modkit:", np.unique(modkit_df['chr']))

    # Copy probes data for methylation processing
    probes_methyl_df = deepcopy(probes_df)

    # Get unique chromosomes
    chromosomes = np.unique(probes_df["chr"].astype(str))

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
    calls_per_probe.to_csv(output_file + ".tmp", header=True, index=False, sep="\t")

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

    return calls_per_probe


def predict_sample_from_dataframe(
    bed_df: pd.DataFrame,
    # model_file: str,
    # output_dir: str,
    sample_name: str = "sample",
    plot_results: bool = False,
):
    """
    Runs `predict_sample` using a DataFrame instead of a BED file.

    Parameters:
        bed_df (pd.DataFrame): The BED file data as a DataFrame.
        model_file (str): Path to the trained model file.
        #output_dir (str): Directory to save results.
        #sample_name (str): Name identifier for output files.
        plot_results (bool): Whether to generate a plot.

    Returns:
        pd.DataFrame: The prediction results.
    """
    # Ensure output directory exists
    # os.makedirs(output_dir, exist_ok=True)
    modelfile = os.path.join(
        os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
    )

    # Validate and load model
    model_path = get_model_path(modelfile)
    if not validate_model_file(modelfile):
        raise ValueError(f"Invalid model file: {modelfile}")

    logging.info(f"Loading model from {modelfile}...")
    inference_session, probes_df, decoding_dict, temperatures, merge_dict = load_model(
        model_path
    )

    # Save DataFrame as a temporary BED file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".bed") as tmp_bed:
        temp_bed_file = tmp_bed.name
        bed_df.to_csv(temp_bed_file, sep="\t", index=False, header=True)

    # Run prediction directly on the DataFrame
    logging.info("Running prediction on the provided DataFrame...")
    prediction_df = predict_sample(
        inference_session=inference_session,
        bed_file=temp_bed_file,
        decoding_dict=deepcopy(decoding_dict),
        probes_df=probes_df,
        temperatures=temperatures,
        merge_dict=merge_dict,
    )

    # Save results as CSV
    # output_csv = os.path.join(output_dir, f"{sample_name}_prediction.csv")
    # logging.info(f"Saving predictions to {output_csv}")
    # prediction_df.to_csv(output_csv, index=False)
    """
    # Plot results if requested
    if plot_results:
        output_pdf = os.path.join(output_dir, f"{sample_name}_prediction.pdf")
        logging.info(f"Plotting results to {output_pdf}")

        with zipfile.ZipFile(model_path, 'r') as zipf:
            try:
                color_dict = json.load(zipf.open('colors.json'))
            except FileNotFoundError:
                color_dict = None  # No color dictionary in model

        plot_prediction(prediction_df=prediction_df, color_dict=color_dict, output_file=output_pdf)
    """
    return prediction_df


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


def run_probes_methyl_calls(merged_output_file, bed_output_file):
    """
    Convert merged methylation calls to BED format.

    Parameters
    ----------
    merged_output_file : str
        Path to the input file containing merged methylation calls
    bed_output_file : str
        Path where the output BED file should be written

    Notes
    -----
    This function is a wrapper around sturgeon's probes_methyl_calls_to_bed function.
    """
    probes_methyl_calls_to_bed(merged_output_file, bed_output_file)


def run_sturgeon_merge_probes(calls_per_probe_file, merged_output_file):
    """
    Merge multiple probe methylation call files.

    Parameters
    ----------
    calls_per_probe_file : str
        Path to the input file containing methylation calls per probe
    merged_output_file : str
        Path where the merged output should be written

    Notes
    -----
    This function merges multiple probe methylation call files into a single output file
    using sturgeon's merge_probes_methyl_calls function.
    """
    merge_probes_methyl_calls(
        [calls_per_probe_file, merged_output_file],
        merged_output_file,
    )


def pysam_cat(tempbam, tomerge):
    """
    Concatenate multiple BAM files using pysam.

    Parameters
    ----------
    tempbam : str
        Path where the concatenated BAM file should be written
    tomerge : list
        List of BAM file paths to merge

    Notes
    -----
    This function uses pysam's cat functionality to merge multiple BAM files
    into a single output file.
    """
    pysam.cat("-o", tempbam, *tomerge)


def run_modkit(file, temp, threads):
    """
    Run modkit to extract methylation data from a BAM file.

    Parameters
    ----------
    file : str
        Path to input BAM file
    temp : str
        Path for temporary output
    threads : int
        Number of threads to use for processing

    Notes
    -----
    This function adapts its command based on the installed modkit version.
    For versions >= 0.4, it uses 'extract full', otherwise just 'extract'.
    """
    try:
        # Get modkit version
        import subprocess

        version_output = subprocess.check_output(
            ["modkit", "--version"], text=True
        ).strip()
        version = version_output.split()[-1]  # Gets '0.4.1' from 'mod_kit 0.4.1'

        # Parse version number
        major, minor, *_ = version.split(".")
        version_num = float(f"{major}.{minor}")

        # Choose appropriate command based on version
        extract_cmd = "extract full" if version_num >= 0.4 else "extract"

        os.system(
            f"modkit {extract_cmd} --ignore h -t {threads} {file} {temp} "
            f"--force --suppress-progress >/dev/null 2>&1"
        )
    except Exception as e:
        logger.error(f"Error in run_modkit: {e}")
        pass


def run_sturgeon_predict(bedDir, dataDir, modelfile):
    """
    Run Sturgeon prediction on methylation data.

    Parameters
    ----------
    bedDir : str
        Directory containing BED format methylation data
    dataDir : str
        Directory where prediction output will be saved
    modelfile : str
        Path to the Sturgeon model file to use for predictions

    Notes
    -----
    This function executes the Sturgeon predict command with the specified
    model file and suppresses terminal output.
    """
    os.system(
        f"sturgeon predict -i {bedDir} -o {dataDir} "
        f"--model-files {modelfile} >/dev/null 2>&1"
    )


def run_sturgeon_inputtobed(temp, temp2):
    """
    Convert Sturgeon input format to BED format.

    Parameters
    ----------
    temp : str
        Path to input file in modkit format
    temp2 : str
        Directory where BED format output will be saved

    Notes
    -----
    This function converts modkit format methylation data to BED format
    using hg38 as the reference genome. Errors are caught and printed.
    """
    try:
        os.system(
            f"sturgeon inputtobed -i {temp} -o {temp2} -s modkit "
            f"--reference-genome hg38 >/dev/null 2>&1"
        )
    except Exception as e:
        logger.error(f"Error in run_sturgeon_inputtobed: {e}")
        pass


class Sturgeon_object(BaseAnalysis):
    """
    A class for processing and visualizing Sturgeon methylation analysis results.

    This class extends BaseAnalysis to provide specialized functionality for
    Sturgeon methylation analysis. It handles real-time processing of BAM files,
    methylation data extraction, and visualization of results.

    Parameters
    ----------
    *args
        Variable length argument list passed to BaseAnalysis
    **kwargs
        Arbitrary keyword arguments passed to BaseAnalysis

    Attributes
    ----------
    sturgeon_df_store : dict
        Store for Sturgeon dataframes
    threshold : float
        Confidence threshold for predictions (default: 0.05)
    first_run : dict
        Tracks first run status for each sample
    modelfile : str
        Path to the Sturgeon model file
    dataDir : dict
        Dictionary of temporary directories for data
    bedDir : dict
        Dictionary of temporary directories for BED files
    st_num_probes : dict
        Dictionary tracking number of probes per sample
    """

    def __init__(self, *args, **kwargs):
        self.sturgeon_df_store = {}
        self.threshold = 0.05
        self.first_run = {}
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        self.dataDir = {}
        self.bedDir = {}
        self.st_num_probes = {}
        reference_genome = "hg38"
        self.probes_file = os.path.join(
            os.path.dirname(sturgeon.__file__),
            "include/static",
            "probes_{}.bed".format(reference_genome),
        )
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        """
        Set up the user interface components for Sturgeon visualization.

        This method creates and arranges the UI cards and grids for displaying
        Sturgeon analysis results. It sets up:
        - Main card with dark theme
        - Grid layout for charts
        - Sturgeon classification chart
        - Time series chart
        - Initial classification label

        Notes
        -----
        The layout is responsive and adjusts based on the MENU_BREAKPOINT value.
        """
        self.card = ui.card().classes("dark:bg-black w-full p-2")
        with self.card:
            with ui.grid(columns=8).classes("w-full h-auto gap-2"):
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-3 max-[{self.MENU_BREAKPOINT}px]:col-span-8 dark:bg-black shadow-lg rounded-lg p-2"
                ):
                    self.create_sturgeon_chart("Sturgeon")
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-5 max-[{self.MENU_BREAKPOINT}px]:col-span-8 dark:bg-black shadow-lg rounded-lg p-2"
                ):
                    self.create_sturgeon_time_chart("Sturgeon Time Series")
        if self.summary:
            with self.summary:
                ui.label(f"Sturgeon classification: Unknown")
        if self.browse:
            self.show_previous_data()
        else:
            ui.timer(5, lambda: self.show_previous_data())

    def show_previous_data(self):
        """
        Load and display previously generated Sturgeon analysis results.
        """
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        if self.check_file_time(os.path.join(output, "sturgeon_scores.csv")):
            self.sturgeon_df_store = pd.read_csv(
                os.path.join(os.path.join(output, "sturgeon_scores.csv")),
                index_col=0,
            )
            columns_greater_than_threshold = (
                self.sturgeon_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = self.sturgeon_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()
            self.update_sturgeon_time_chart(self.sturgeon_df_store.drop(columns=result))
            lastrow = self.sturgeon_df_store.iloc[-1].drop("number_probes")
            lastrow_plot = lastrow.sort_values(ascending=False).head(10)
            lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)

            # Update summary with new card
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    classification_text = (
                        f"Sturgeon classification: {lastrow_plot_top.index[0]}"
                    )
                    self.create_summary_card(
                        classification_text=classification_text,
                        confidence_value=lastrow_plot_top.values[0],
                        features_found=int(
                            self.sturgeon_df_store.iloc[-1]["number_probes"]
                        ),
                    )

            self.update_sturgeon_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values),
                "All",
                self.sturgeon_df_store.iloc[-1]["number_probes"],
            )

    async def process_bam(self, bamfile: List[Tuple[str, float]]) -> None:
        """
        Process BAM files and perform Sturgeon analysis.

        This method handles the complete workflow of processing BAM files:
        1. Merging BAM files
        2. Extracting methylation data using modkit
        3. Converting to BED format
        4. Running Sturgeon predictions
        5. Updating visualizations

        Parameters
        ----------
        bamfile : List[Tuple[str, float]]
            List of tuples containing BAM file paths and their timestamps

        Notes
        -----
        The method processes files in batches and updates progress trackers.
        Results are stored in temporary directories and visualized in real-time.
        """
        sampleID = self.sampleID
        # Initialize directories for each sampleID if not already present
        if sampleID not in self.dataDir.keys():
            self.dataDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
            self.bedDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
        tomerge = []
        latest_file = 0
        if (
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bam_count"
            ]
            > 0
        ):
            if latest_file:
                currenttime = latest_file * 1000
            else:
                currenttime = time.time() * 1000

            parquet_path = os.path.join(
                self.check_and_create_folder(self.output, sampleID),
                f"{sampleID}.parquet",
            )

            try:
                if self.check_file_time(parquet_path):
                    tomerge_length_file = os.path.join(
                        self.check_and_create_folder(self.output, sampleID),
                        "tomerge_length.txt",
                    )
                    with open(tomerge_length_file, "r") as f:
                        tomerge_length = int(f.readline().strip().split(": ")[1])

                    merged_modkit_df = await run.cpu_bound(
                        load_modkit_data, parquet_path
                    )

                    temp_pileup = tempfile.NamedTemporaryFile(
                        dir=self.check_and_create_folder(self.output, sampleID)
                    )

                    # Pass the cleaned data
                    result_df = await run.cpu_bound(
                        modkit_pileup_file_to_bed,
                        merged_modkit_df,
                        temp_pileup.name,
                        self.probes_file,
                    )

                    diagnosis = predict_sample_from_dataframe(result_df)
                    self.st_num_probes[sampleID] = diagnosis.iloc[-1]["number_probes"]
                    # lastrow = mydf.iloc[-1].drop("number_probes")
                    mydf_to_save = diagnosis
                    mydf_to_save["timestamp"] = currenttime

                    if sampleID not in self.sturgeon_df_store:
                        self.sturgeon_df_store[sampleID] = pd.DataFrame()

                    # Exclude empty or all-NA entries before concatenation
                    if not mydf_to_save.dropna(how="all").empty:
                        self.sturgeon_df_store[sampleID] = pd.concat(
                            [
                                self.sturgeon_df_store[sampleID],
                                mydf_to_save.set_index("timestamp"),
                            ]
                        )
                        self.sturgeon_df_store[sampleID].to_csv(
                            os.path.join(
                                self.check_and_create_folder(self.output, sampleID),
                                "sturgeon_scores.csv",
                            )
                        )
                    app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                        "bam_processed"
                    ] = tomerge_length

                    # app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                    #    "bams_in_processing"
                # ] = 0
                # app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                #    "bam_count"
                # ] -= tomerge_length
                # print(app.storage.general[self.mainuuid][sampleID][self.name]["counters"])
            except Exception as e:
                logger.error(f"Error in process_bam (sturgeon): {e}")

        self.running = False

    def create_sturgeon_chart(self, title):
        """
        Create a bar chart for displaying Sturgeon classification results.
        """
        self.echart2 = self.create_chart(title)
        self.echart2.options.update(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": title,
                    "left": "center",
                    "top": 10,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "normal",
                        "color": "#000000",
                    },
                },
                "tooltip": {
                    "trigger": "axis",
                    "axisPointer": {"type": "shadow"},
                    "formatter": "{b}: {c}%",
                    "textStyle": {"fontSize": 14},
                },
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "5%",
                    "top": "25%",
                    "containLabel": True,
                },
                "xAxis": {
                    "type": "value",
                    "min": 0,
                    "max": 100,
                    "interval": 20,
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{value}%",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "yAxis": {
                    "type": "category",
                    "inverse": True,
                    "data": [],
                    "axisLabel": {
                        "fontSize": 12,
                        "width": 250,
                        "overflow": "break",
                        "interval": 0,
                        "align": "right",
                        "color": "#666666",
                    },
                },
                "series": [
                    {
                        "type": "bar",
                        "name": "Confidence",
                        "barMaxWidth": "60%",
                        "itemStyle": {"color": "#007AFF", "borderRadius": [0, 4, 4, 0]},
                        "label": {
                            "show": True,
                            "position": "right",
                            "formatter": "{c}%",
                            "fontSize": 12,
                            "color": "#666666",
                        },
                        "data": [],
                    }
                ],
            }
        )

    def create_sturgeon_time_chart(self, title):
        """
        Create a time series chart for Sturgeon results.
        """
        self.sturgeon_time_chart = self.create_time_chart(title)
        self.sturgeon_time_chart.options.update(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": title,
                    "left": "center",
                    "top": 5,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "normal",
                        "color": "#000000",
                    },
                    "padding": [0, 0, 20, 0],
                },
                "tooltip": {
                    "trigger": "axis",
                    "axisPointer": {"type": "line"},
                    "textStyle": {"fontSize": 14},
                },
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "5%",
                    "top": "25%",
                    "containLabel": True,
                },
                "legend": {
                    "type": "scroll",
                    "orient": "horizontal",
                    "top": 45,
                    "width": "90%",
                    "left": "center",
                    "textStyle": {"fontSize": 12, "color": "#666666"},
                    "pageButtonPosition": "end",
                    "pageButtonGap": 5,
                    "pageButtonItemGap": 5,
                    "pageIconColor": "#666666",
                    "pageIconInactiveColor": "#aaa",
                    "pageIconSize": 12,
                    "pageTextStyle": {"color": "#666666"},
                    "itemGap": 25,
                    "itemWidth": 14,
                    "itemHeight": 14,
                    "selectedMode": True,
                },
                "xAxis": {
                    "type": "time",
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{HH}:{mm}",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "yAxis": {
                    "type": "value",
                    "min": 0,
                    "max": 100,
                    "interval": 20,
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{value}%",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
            }
        )

    def update_sturgeon_plot(self, x, y, count, st_num_probes):
        """
        Update the Sturgeon visualization plot with new data.

        Parameters
        ----------
        x : List[str]
            List of tumor types
        y : List[float]
            Confidence scores for each tumor type
        count : str
            Number of BAM files processed
        st_num_probes : int
            Number of probes found
        """
        # Convert values to percentages and format
        formatted_values = [float(f"{val * 100:.1f}") for val in y]

        # Sort the data in descending order
        sorted_indices = sorted(
            range(len(formatted_values)),
            key=lambda k: formatted_values[k],
            reverse=True,
        )
        sorted_values = [formatted_values[i] for i in sorted_indices]
        sorted_labels = [x[i] for i in sorted_indices]

        # Create descriptive title with key information
        title_text = (
            f"Sturgeon Analysis Results\n"
            f"{count} samples processed • {int(st_num_probes)} probes found"
        )

        self.echart2.options["title"]["text"] = title_text
        self.echart2.options["yAxis"]["data"] = sorted_labels
        self.echart2.options["series"][0].update(
            {
                "data": sorted_values,
                "itemStyle": {"color": "#007AFF", "borderRadius": [0, 4, 4, 0]},
            }
        )
        self.echart2.update()

    def update_sturgeon_time_chart(self, datadf):
        """
        Update the time series chart with new data.

        Parameters
        ----------
        datadf : pd.DataFrame
            DataFrame containing time series data for visualization
        """
        self.sturgeon_time_chart.options["series"] = []

        # iOS color palette for multiple series
        colors = [
            "#007AFF",  # Blue
            "#34C759",  # Green
            "#FF9500",  # Orange
            "#FF2D55",  # Red
            "#5856D6",  # Purple
            "#FF3B30",  # Red-Orange
            "#5AC8FA",  # Light Blue
            "#4CD964",  # Light Green
        ]

        for idx, (series, data) in enumerate(datadf.to_dict().items()):
            if series != "number_probes":
                # Convert values to percentages
                data_list = [
                    [key, float(f"{value * 100:.1f}")] for key, value in data.items()
                ]
                self.sturgeon_time_chart.options["series"].append(
                    {
                        "name": series,
                        "type": "line",
                        "smooth": True,
                        "animation": False,
                        "symbolSize": 6,
                        "emphasis": {
                            "focus": "series",
                            "itemStyle": {"borderWidth": 2},
                        },
                        "endLabel": {
                            "show": True,
                            "formatter": "{a}: {c}%",
                            "distance": 10,
                            "fontSize": 12,
                        },
                        "lineStyle": {"width": 2, "color": colors[idx % len(colors)]},
                        "itemStyle": {"color": colors[idx % len(colors)]},
                        "data": data_list,
                    }
                )

        # Update chart title with summary
        latest_data = datadf.iloc[-1].drop(
            "number_probes"
        )  # Remove number_probes from latest data
        max_confidence = latest_data.max() * 100  # Convert to percentage
        max_type = latest_data.idxmax()
        self.sturgeon_time_chart.options["title"]["text"] = (
            f"Classification Confidence Over Time\n"
            f"Current highest confidence: {max_type} ({max_confidence:.1f}%)"
        )

        self.sturgeon_time_chart.update()


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
):
    """
    Test function to run the Sturgeon analysis interface.

    Parameters
    ----------
    port : int
        Port number for the web interface
    threads : int
        Number of processing threads
    watchfolder : str
        Directory to watch for new BAM files
    output : str
        Directory for output files
    reload : bool, optional
        Whether to enable auto-reload (default: False)
    browse : bool, optional
        Whether to run in browse mode (default: False)

    Notes
    -----
    This function sets up the web interface and starts processing BAM files
    either from a watch folder or in browse mode.
    """
    my_connection = None
    with theme.frame("Sturgeon Rapid CNS Diagnostic.", my_connection):
        TestObject = Sturgeon_object(threads, output, progress=True, batch=True)
    if not browse:
        path = watchfolder
        searchdirectory = os.fsencode(path)
        for root, d_names, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
    else:
        logger.info("Browse mode not implemented.")
        TestObject.progress_trackers.visible = False
        # TestObject.show_previous_data(output)
    ui.run(port=port, reload=False)


@click.command()
@click.option(
    "--port",
    default=12345,
    help="Port for GUI",
)
@click.option("--threads", default=4, help="Number of threads available.")
@click.argument(
    "watchfolder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.argument(
    "output",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.option(
    "--browse",
    is_flag=True,
    show_default=True,
    default=False,
    help="Browse Historic Data.",
)
def mainrun(port, threads, watchfolder, output, browse):
    """
    Main entry point for the Sturgeon analysis application.

    This function handles command line arguments and launches the application
    either in watch mode or browse mode.

    Parameters
    ----------
    port : int
        Port number for the web interface
    threads : int
        Number of processing threads
    watchfolder : Path
        Directory to watch for new BAM files
    output : Path
        Directory for output files
    browse : bool
        Whether to run in browse mode

    Notes
    -----
    In watch mode, both watchfolder and output are required.
    In browse mode, only output is required.
    """
    if browse:
        # Handle the case when --browse is set
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            # simtime=simtime,
            watchfolder=None,
            output=watchfolder,
            # sequencing_summary=sequencing_summary,
            # showerrors=showerrors,
            browse=browse,
            # exclude=exclude,
        )
        # Your logic for browse mode
    else:
        # Handle the case when --browse is not set
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            # simtime=simtime,
            watchfolder=watchfolder,
            output=output,
            # sequencing_summary=sequencing_summary,
            # showerrors=showerrors,
            browse=browse,
            # exclude=exclude,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    mainrun()
