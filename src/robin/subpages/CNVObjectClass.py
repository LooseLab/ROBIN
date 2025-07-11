"""
Copy Number Variation (CNV) Analysis Module
=========================================

This module provides functionality for analyzing copy number variations (CNVs) from BAM files
generated during sequencing runs. It implements real-time CNV detection and visualization
following Apple Human Interface Guidelines.

Key Features
-----------
- Real-time CNV detection from BAM files
- Interactive visualization of CNV data
- Genetic sex estimation (XX/XY)
- Break point detection and analysis
- Time series tracking of CNV changes

Classes
-------
Result
    Simple container class for CNV results.

CNVAnalysis
    Main analysis class inheriting from BaseAnalysis, providing CNV-specific functionality.

Functions
---------
run_ruptures(r_cnv, penalty_value, bin_width)
    Detects change points in CNV data using Kernel Change Point Detection.

iterate_bam(bamfile, threads, mapq_filter, copy_numbers, log_level)
    Processes BAM files to extract CNV data.

moving_average(data, n)
    Calculates moving averages for smoothing CNV data.

pad_arrays(arr1, arr2, pad_value)
    Utility function to pad arrays to equal lengths.

Visualization
------------
The module implements several visualization components following Apple HIG:

- CNV scatter plots with interactive zooming
- Time series charts for tracking changes
- Summary cards with key metrics
- Chromosome and gene selection controls

Dependencies
-----------
- cnv_from_bam.iterate_bam_file
- robin.subpages.base_analysis.BaseAnalysis
- natsort
- pandas
- logging
- numpy
- nicegui (ui, app, run)
- ruptures

Example Usage
-----------
.. code-block:: python

    from robin.subpages.CNV_object import CNVAnalysis

    # Initialize analysis
    cnv_analysis = CNVAnalysis(
        threads=4,
        output_dir="output/",
        target_panel="rCNS2"
    )

    # Process BAM files
    cnv_analysis.add_bam("sample.bam")

    # Run the UI
    cnv_analysis.setup_ui()

Notes
-----
The module follows Apple Human Interface Guidelines for:
- Color usage and accessibility
- Typography and spacing
- Interactive elements
- Visual feedback
- Information hierarchy

Authors
-------
Matt Loose
"""

import pysam
from cnv_from_bam import iterate_bam_file
from robin.subpages.base_analysis import BaseAnalysis
from robin.utilities.break_point_detector import CNVChangeDetectorTracker
from robin.utilities.bed_file import MasterBedTree, BedTree
from robin import theme, resources
import pandas as pd
import logging
import numpy as np
import os
import sys
from nicegui import ui, app, run, background_tasks
import click
from pathlib import Path
import pickle
import ruptures as rpt
from typing import Optional, Tuple, BinaryIO, List
import re

from collections import Counter
from scipy.ndimage import uniform_filter1d

import time

from robin.core.state import state, ProcessState
import shutil

os.environ["CI"] = "1"
# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


class Result:
    """
    A class to store CNV results.
    """

    def __init__(self, cnv_dict: dict) -> None:
        self.cnv = cnv_dict


def run_ruptures(
    r_cnv: list, penalty_value: int, bin_width: int
) -> List[Tuple[int, int]]:
    """
    Detect change points in a time series using the Kernel Change Point Detection algorithm from the ruptures library.

    Args:
        r_cnv: A list containing the time series data.
        penalty_value (int): The penalty value for the change point detection algorithm.
        bin_width (int): The width of the bins used to scale the detected change points.

    Returns:
        List[Tuple[float, float]]: A list of tuples where each tuple represents a detected change point range
                                    as (start, end) with respect to the bin width.
    """
    # Initialize the Kernel Change Point Detection algorithm with the Radial Basis Function (RBF) kernel
    # x_coords = range(0, (len(r_cnv) + 1))
    # signal = np.array(list(zip(x_coords, r_cnv)))
    algo_c = rpt.KernelCPD(kernel="rbf").fit(np.array(r_cnv))

    # Predict the change points using the provided penalty value
    ruptures_result = algo_c.predict(pen=penalty_value)

    # Compute the ranges around each change point
    return [
        (cp * bin_width - (bin_width), cp * bin_width + (bin_width))
        for cp in ruptures_result
    ]


def reduce_list(lst: List, max_length: int = 1000) -> List:
    """
    Reduce the length of a list to a specified maximum length by subsampling.

    Args:
        lst (List): The list to reduce.
        max_length (int): The maximum length of the list.

    Returns:
        List: The reduced list.
    """
    while len(lst) > max_length:
        lst = lst[::2]
    return lst


class CNV_Difference:
    """
    A class to store CNV difference data.
    """

    def __init__(self, *args, **kwargs) -> None:
        self.cnv = {}


def moving_average(data: np.ndarray, n: int = 3) -> np.ndarray:
    """
    Calculate the moving average of a given data array using scipy's uniform_filter1d.

    Args:
        data (np.ndarray): The data array.
        n (int): The window size for the moving average.

    Returns:
        np.ndarray: The array of moving averages.
    """
    return uniform_filter1d(data, size=n, mode="nearest")


def moving_average_orig(data: np.ndarray, n: int = 3) -> np.ndarray:
    """
    Calculate the moving average of a given data array.

    Args:
        data (np.ndarray): The data array.
        n (int): The window size for the moving average.

    Returns:
        np.ndarray: The array of moving averages.
    """
    return np.convolve(data, np.ones(n) / n, mode="same")


def pad_arrays(
    arr1: np.ndarray, arr2: np.ndarray, pad_value: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pad two arrays to the same length with a specified value.

    Args:
        arr1 (np.ndarray): The first array.
        arr2 (np.ndarray): The second array.
        pad_value (int): The value to pad with.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The padded arrays.
    """
    len1, len2 = len(arr1), len(arr2)
    max_len = max(len1, len2)

    if len1 < max_len:
        arr1 = np.pad(
            arr1, (0, max_len - len1), mode="constant", constant_values=pad_value
        )
    if len2 < max_len:
        arr2 = np.pad(
            arr2, (0, max_len - len2), mode="constant", constant_values=pad_value
        )

    return arr1, arr2


def pad_arrays_orig(
    arr1: np.ndarray, arr2: np.ndarray, pad_value: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pad two arrays to the same length with a specified value.

    Args:
        arr1 (np.ndarray): The first array.
        arr2 (np.ndarray): The second array.
        pad_value (int): The value to pad with.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The padded arrays.
    """
    len_diff = abs(len(arr1) - len(arr2))
    if len(arr1) < len(arr2):
        arr1 = np.pad(arr1, (0, len_diff), mode="constant", constant_values=pad_value)
    elif len(arr1) > len(arr2):
        arr2 = np.pad(arr2, (0, len_diff), mode="constant", constant_values=pad_value)
    return arr1, arr2


def get_data(output: str) -> Tuple[Result, Result, Result, dict]:
    """
    Get data from the output directory.
    """
    result = Result(
        np.load(os.path.join(output, "CNV.npy"), allow_pickle="TRUE").item()
    )
    result2 = Result(
        np.load(os.path.join(output, "CNV2.npy"), allow_pickle="TRUE").item()
    )
    result3 = Result(
        np.load(os.path.join(output, "CNV3.npy"), allow_pickle="TRUE").item()
    )
    cnv_dict = np.load(os.path.join(output, "CNV_dict.npy"), allow_pickle="TRUE").item()
    return result, result2, result3, cnv_dict

def iterate_bam(
    bamfile,
    _threads: int,
    mapq_filter: int,
    copy_numbers: dict,
    log_level: int,
    ref_cnv_dict,
) -> Tuple[dict, int, float, dict]:
    """
    Iterate over a BAM file and return CNV data and associated metrics.

    Args:
        bamfile (BinaryIO): The BAM file to process.
        _threads (int): Number of threads to use.
        mapq_filter (int): MAPQ filter value.
        copy_numbers (dict): Dictionary to store copy number data.
        log_level (int): Logging level.

    Returns:
        Tuple[dict, int, float, dict]: CNV data, bin width, variance, and updated copy numbers.
    """
    result = iterate_bam_file(
        bamfile,
        _threads=_threads,
        mapq_filter=mapq_filter,
        copy_numbers=copy_numbers,
        log_level=log_level,
    )

    result2 = iterate_bam_file(
        bamfile,
        _threads=_threads,
        mapq_filter=mapq_filter,
        copy_numbers=ref_cnv_dict,
        log_level=log_level,
        bin_width=result.bin_width,
    )

    return (
        result.cnv,
        result.bin_width,
        result.variance,
        copy_numbers,
        result.genome_length,
        result2.cnv,
    )


class CNVAnalysis(BaseAnalysis):
    """
    Class for analyzing copy number variations (CNVs) from BAM files.

    Inherits from `BaseAnalysis` and provides specific methods for CNV analysis.
    """

    def __init__(
        self,
        *args,
        target_panel: Optional[str] = None,
        reference_file: Optional[str] = None,
        bed_file: Optional[str] = None,
        readfish_toml: Optional[Path] = None,
        # This needs to be shared between different tools but is most likely a memory hog
        master_bed_tree: Optional[MasterBedTree] = None,
        **kwargs,
    ) -> None:
        # self.file_list = []
        super().__init__(*args, **kwargs)
        # Remove state tracking for CNV Analysis
        # state.start_process("CNV Analysis", ProcessType.BATCH)
        state.set_process_state("CNV Analysis", ProcessState.WAITING_FOR_DATA)
        self.target_panel = target_panel
        self.reference_file = reference_file
        self.bed_file = bed_file
        self.readfish_toml = readfish_toml
        # Define dtype for memmap - using numpy dtype
        self.dtype = np.dtype([("name", "U10"), ("start", "i8"), ("end", "i8")])
        #self.cnv_dict = {"bin_width": 0, "variance": 0}
        self.result = None
        self.result3 = CNV_Difference()
        self.ref_result = None
        self.working_dir = None
        self.chromosome_events = {}
        self.sex_estimate = None
        self.timer1 = None
        self.bin_width = 1_000_000
        self.total_reads = 0
        self.current_file = None
        self.timest = 0
        self.chrom_filter = "All"
        # Color mode: "chromosome" for coloring by chromosome, "value" for red/blue based on values
        self.color_mode = "chromosome"

        # Initialize data array related attributes
        self.local_data_array = None

        # self.len_tracker = defaultdict(lambda: 0)
        self.map_tracker = Counter()
        self.load_data = False
        with open(
            os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "HG01280_control_new.pkl",
            ),
            "rb",
        ) as f:
            self.ref_cnv_dict = pickle.load(f)
        if self.target_panel == "rCNS2":
            self.gene_bed_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "rCNS2_panel_name_uniq.bed",
            )
        elif self.target_panel == "AML":
            self.gene_bed_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "AML_panel_name_uniq.bed",
            )
        self.gene_bed = pd.read_table(
            self.gene_bed_file,
            names=["chrom", "start_pos", "end_pos", "gene"],
            header=None,
            sep="\s+",
        )
        self.centromeres_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "cenSatRegions.bed",
        )
        self.centromere_bed = pd.read_csv(
            self.centromeres_file,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start_pos", "end_pos", "name"],
            header=None,
            sep="\s+",
        )

        # Add cytobands file loading
        self.cytobands_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "cytoBand.txt",
        )
        self.cytobands_bed = pd.read_csv(
            self.cytobands_file,
            names=["chrom", "start_pos", "end_pos", "name", "gieStain"],
            header=None,
            sep="\s+",
        )

        self.master_bed_tree = master_bed_tree
        # super().__init__(*args, **kwargs)

        # Add state directory path
        self.state_dir = None  # Will be set in do_cnv_work

        # Add target_table as instance variable
        self.target_table = None
        self.target_table_placeholder = None
        self.last_bed_check = 0  # Track when we last checked for new BED files

        # Add y_axis_log attribute
        self.y_axis_log = False  # Default to linear scale
        # Add show_breakpoints attribute
        self.show_breakpoints = True  # Default to showing breakpoints

        

    def _get_state_dir(self) -> Path:
        """Get the directory for storing CNV detector state."""
        if not self.state_dir:
            self.state_dir = (
                Path(self.check_and_create_folder(self.output, self.sampleID))
                / "cnv_detector_state"
            )
            self.state_dir.mkdir(exist_ok=True)
        return self.state_dir

    def _cleanup_state(self) -> None:
        """Clean up temporary files but preserve CNV detector state."""
        # Don't delete the state directory as it contains the CNV detector state
        # that needs to persist between processing calls
        pass

    def _get_update_cnv_dict_path(self) -> str:
        """Get the path for storing the update_cnv_dict."""
        return os.path.join(
            self.check_and_create_folder(self.output, self.sampleID),
            "update_cnv_dict.pkl",
        )

    def _load_update_cnv_dict(self) -> dict:
        """Load the update_cnv_dict from file."""
        update_cnv_dict_path = self._get_update_cnv_dict_path()
        if os.path.exists(update_cnv_dict_path):
            try:
                with open(update_cnv_dict_path, "rb") as f:
                    update_cnv_dict = pickle.load(f)
                logger.info("Successfully loaded update_cnv_dict from file")
                return update_cnv_dict
            except Exception as e:
                logger.warning(f"Failed to load update_cnv_dict: {e}")
                return {}
        return {}

    def _save_update_cnv_dict(self, update_cnv_dict: dict) -> None:
        """Save the update_cnv_dict to file."""
        try:
            update_cnv_dict_path = self._get_update_cnv_dict_path()
            with open(update_cnv_dict_path, "wb") as f:
                pickle.dump(update_cnv_dict, f)
            logger.info("Successfully saved update_cnv_dict to file")
        except Exception as e:
            logger.error(f"Failed to save update_cnv_dict: {e}")
            raise

    def get_update_cnv_dict(self) -> dict:
        """Get the current update_cnv_dict data from file."""
        return self._load_update_cnv_dict()



    def calculate_chromosome_stats(self, result, ref_result):
        """Calculate chromosome-wide statistics and baselines.

        Args:
            result: CNV result object for sample
            ref_result: CNV result object for reference

        Returns:
            Dictionary of chromosome statistics including means, baselines, and thresholds
        """
        stats = {}
        autosome_means = []

        # Calculate normalized values and stats for each chromosome
        for chrom in result.cnv.keys():
            if chrom != "chrM" and chrom in ref_result:
                # Calculate normalized CNV values
                sample_avg = moving_average(result.cnv[chrom])
                ref_avg = moving_average(ref_result[chrom])

                # Pad arrays if needed
                max_len = max(len(sample_avg), len(ref_avg))
                if len(sample_avg) < max_len:
                    sample_avg = np.pad(sample_avg, (0, max_len - len(sample_avg)))
                if len(ref_avg) < max_len:
                    ref_avg = np.pad(ref_avg, (0, max_len - len(ref_avg)))

                # Calculate normalized CNV
                normalized_cnv = sample_avg - ref_avg

                # Calculate basic statistics
                chr_mean = np.mean(normalized_cnv)
                chr_std = np.std(normalized_cnv)

                # Store autosome means for global statistics
                if chrom.startswith("chr") and chrom[3:].isdigit():
                    autosome_means.append(chr_mean)

                # Set baseline and thresholds based on chromosome and sex
                if chrom == "chrX":
                    if self.XYestimate == "XX":  # Female
                        baseline = 1.0  # Expected +1 relative to male control
                    else:  # Male
                        baseline = 0.0  # Expected same as male control
                elif chrom == "chrY":
                    if self.XYestimate == "XY":  # Male
                        baseline = 0.0
                    else:  # Female
                        baseline = -1.0  # Expected absence
                else:  # Autosomes
                    baseline = 0.0

                stats[chrom] = {
                    "mean": chr_mean,
                    "std": chr_std,
                    "baseline": baseline,
                    "normalized_cnv": normalized_cnv,
                }

        # Calculate global autosome statistics
        global_mean = np.mean(autosome_means)
        global_std = np.std(autosome_means)

        # Store global stats
        stats["global"] = {"mean": global_mean, "std": global_std}

        self.chromosome_stats = stats
        return stats

    def estimate_XY(self) -> None:
        """
        Estimate genetic sex (XX or XY) based on CNV data.
        """
        X = round(np.average([i for i in self.result3.cnv["chrX"] if i != 0]), 2)
        Y = round(np.average([i for i in self.result3.cnv["chrY"] if i != 0]), 2)
        if X >= 0.1 and Y <= -0.1:
            self.sex_estimate = "Female"
        elif X <= 0.1 and Y >= -0.2:
            self.sex_estimate = "Male"
        elif X >= 0.1 and Y >= -0.1:
            self.sex_estimate = "Male (query X/Y copy number changes)"
        elif X > 0.1 and Y > 0.1:
            self.sex_estimate = "Male (query X/Y copy number changes)"
        elif X < 0.1 and Y < -0.2:
            self.sex_estimate = "Unknown (Query XY copy number changes)"
        else:
            self.sex_estimate = "Unknown"

        # print(X,Y, self.sex_estimate)
        with open(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "XYestimate.pkl",
            ),
            "wb",
        ) as file:
            pickle.dump(self.sex_estimate, file)

    async def process_bam(self, bamfile: BinaryIO, timestamp: float) -> None:
        """
        Process a BAM file to extract CNV data.

        Args:
            bamfile (BinaryIO): The BAM file to process.
            timestamp (float): The timestamp indicating when the file was generated.
        """

        state.set_process_state("CNV Analysis", ProcessState.RUNNING)

        try:

            # Set up state directory
            state_dir = self._get_state_dir()

            # Initialize CNV detector only if it doesn't exist
            
            CNVchangedetector = CNVChangeDetectorTracker(
                base_proportion=0.02, temp_dir=self.output
            )
            
            # Load previous state if it exists
            if (state_dir / "tracker_metadata.pkl").exists():
                try:
                    CNVchangedetector.load_state(str(state_dir))
                    logger.info("Successfully loaded CNV detector state")
                except Exception as e:
                    logger.warning(f"Failed to load CNV detector state: {e}")                    

            # Initialize data structures
            update_cnv_dict = self._load_update_cnv_dict()
            if self.sampleID not in update_cnv_dict.keys():
                update_cnv_dict[self.sampleID] = {}

            # Set up bed tree if needed
            if self.master_bed_tree[self.sampleID] is None:
                self.master_bed_tree.add_bed_tree(
                    sample_id=self.sampleID,
                    preserve_original_tree=True,
                    reference_file=f"{self.reference_file}.fai",
                )
                NewBed = self.master_bed_tree.bed_trees[self.sampleID]

                NewBed.load_from_file(self.bed_file)
            else:
                NewBed = self.master_bed_tree.bed_trees[self.sampleID]

            # Set up BAM file and data array
            bamdata = pysam.AlignmentFile(bamfile, "rb")
            data_array_path = os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "cnv_data_array.npy",
            )

            # Initialize or load data array using standard file operations
            if not os.path.exists(data_array_path):
                # Initialize with empty array
                DATA_ARRAY = np.array([], dtype=self.dtype)
                data_array_size = 0
            else:
                # Load existing data
                try:
                    DATA_ARRAY = np.load(data_array_path, allow_pickle=True)
                    data_array_size = len(DATA_ARRAY)
                except Exception as e:
                    logger.warning(f"Failed to load existing data array: {e}")
                    DATA_ARRAY = np.array([], dtype=self.dtype)
                    data_array_size = 0

            # Update map tracker
            self.map_tracker.update(
                Counter(
                    {
                        stat.contig: stat.mapped
                        for stat in bamdata.get_index_statistics()
                    }
                )
            )

            (
                r_cnv,
                r_bin,
                r_var,
                updated_copy_numbers,
                genome_length,
                r2_cnv,
            ) = await run.cpu_bound(
                iterate_bam,
                bamfile,
                self.threads,
                60,
                update_cnv_dict[self.sampleID],
                int(logging.ERROR),
                self.ref_cnv_dict,
            )
            
            # Update the local copy with the returned data
            update_cnv_dict[self.sampleID] = updated_copy_numbers
            
            cnv_dict = {}
            cnv_dict["bin_width"] = r_bin
            cnv_dict["variance"] = r_var

            # Process breakpoints and update CNV detector
            cnvupdate = False
            for key in r_cnv.keys():
                local_data_array = np.empty(0, dtype=self.dtype)
                if key != "chrM" and re.match(r"^chr(\d+|X|Y)$", key):
                    moving_avg_data1 = await run.cpu_bound(moving_average, r_cnv[key])
                    moving_avg_data2 = await run.cpu_bound(moving_average, r2_cnv[key])
                    moving_avg_data1, moving_avg_data2 = await run.cpu_bound(
                        pad_arrays, moving_avg_data1, moving_avg_data2
                    )
                    self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2

                    if len(r_cnv[key]) > 3 and self.map_tracker[key] > 2000:
                        paired_changepoints = await run.cpu_bound(
                            run_ruptures,
                            r_cnv[key],
                            5,  # penalty_value
                            cnv_dict["bin_width"],
                        )

                        # Check if any change points were detected by the ruptures algorithm
                        if len(paired_changepoints) > 0:
                            # Calculate approximate chromosome length based on number of bins and bin width
                            approx_chrom_length = len(r_cnv[key]) * r_bin
                            # Define padding around centromere regions to avoid false positives
                            padding = 2_500_000  # 2.5 Mb padding

                            # Process each detected change point (start, end) pair
                            for start, end in paired_changepoints:
                                # Filter out invalid or problematic regions:
                                # 1. Skip if the change point spans the entire chromosome (likely false positive)
                                # 2. Skip if start position is negative (invalid coordinate)
                                # 3. Skip if the change point overlaps with centromere regions (known problematic areas)
                                if (
                                    start < approx_chrom_length < end
                                    or start < 0
                                    or (
                                        start
                                        < self.centromere_bed[
                                            self.centromere_bed["chrom"].eq(key)
                                        ]["end_pos"].max()
                                        + padding
                                        and end
                                        > self.centromere_bed[
                                            self.centromere_bed["chrom"].eq(key)
                                        ]["start_pos"].min()
                                        - padding
                                    )
                                ):
                                    continue

                                # Create a structured array entry for this breakpoint
                                # dtype format: (chromosome_name, start_position, end_position)
                                item = np.array([(key, start, end)], dtype=self.dtype)

                                # Append to local array (current processing batch)
                                local_data_array = np.append(
                                    local_data_array, item
                                )
                                # Append to persistent data array (all breakpoints for this sample)
                                DATA_ARRAY = np.append(DATA_ARRAY, item)

                                # Save the updated data array to file immediately
                                try:
                                    np.save(data_array_path, DATA_ARRAY)
                                except Exception as e:
                                    logger.warning(f"Failed to save data array: {e}")

                            # Reset the mapping tracker for this chromosome to prevent reprocessing
                            self.map_tracker[key] = 0
                            # Flag that we have CNV updates to save
                            cnvupdate = True

                            # Extract all breakpoints for this chromosome from persistent storage
                            breakpoints = DATA_ARRAY[
                                DATA_ARRAY["name"] == key
                            ]
                            # Extract breakpoints from current processing batch
                            local_breakpoints = local_data_array[
                                local_data_array["name"] == key
                            ]

                            # Only proceed if we have valid breakpoints to process
                            if len(breakpoints) > 0:
                                try:
                                    # Add the detected breakpoints to the CNV change detector
                                    # This updates the internal tracking and generates BED targets
                                    CNVchangedetector.add_breakpoints(
                                        key,  # chromosome name
                                        breakpoints,  # all historical breakpoints for this chromosome
                                        local_breakpoints,  # newly detected breakpoints
                                        approx_chrom_length,  # chromosome length estimate
                                        r_bin,  # current bin width
                                    )
                                except Exception as e:
                                    logger.error(
                                        f"Error adding breakpoints for {key}: {e}"
                                    )
                                    # Re-raise the exception to ensure proper error handling upstream
                                    raise

            self.estimate_XY()

            # Save state if we have updates
            if cnvupdate:
                try:

                    # Save CNV detector state
                    CNVchangedetector.save_state(str(state_dir))
                    logger.info("Successfully saved CNV detector state")

                    # Generate and save BED content
                    bedcontent = ""
                    bedcontent2 = ""
                    for chrom in r_cnv.keys():
                        tempbedcontent = CNVchangedetector.get_bed_targets(chrom)
                        if len(tempbedcontent) > 0:
                            bedcontent += tempbedcontent + "\n"

                        tempbedcontent2 = (
                            CNVchangedetector.get_bed_targets_breakpoints(chrom)
                        )
                        if len(tempbedcontent2) > 0:
                            bedcontent2 += tempbedcontent2 + "\n"

                    if len(bedcontent2) > 0:
                        NewBed.load_from_string(
                            bedcontent2,
                            merge=False,
                            write_files=True,
                            output_location=str(
                                self.check_and_create_folder(self.output, self.sampleID)
                            ),
                            source_type="CNV",
                        )

                    # Save other analysis results
                    np.save(
                        os.path.join(
                            self.check_and_create_folder(self.output, self.sampleID),
                            "CNV.npy",
                        ),
                        r_cnv,
                    )
                    np.save(
                        os.path.join(
                            self.check_and_create_folder(self.output, self.sampleID),
                            "CNV2.npy",
                        ),
                        r2_cnv,
                    )
                    np.save(
                        os.path.join(
                            self.check_and_create_folder(self.output, self.sampleID),
                            "CNV3.npy",
                        ),
                        self.result3.cnv,
                    )
                    np.save(
                        os.path.join(
                            self.check_and_create_folder(self.output, self.sampleID),
                            "CNV_dict.npy",
                        ),
                        cnv_dict,
                    )

                except Exception as e:
                    logger.error(f"Error saving CNV detector state: {e}")
                    raise

            # Always save the update_cnv_dict since it may have been modified by iterate_bam
            try:
                self._save_update_cnv_dict(update_cnv_dict)
            except Exception as e:
                logger.error(f"Error saving update_cnv_dict: {e}")
                raise

        finally:
            # Clean up temporary files but preserve CNV detector
            self._cleanup_state()
            self.running = False
            state.set_process_state("CNV Analysis", ProcessState.WAITING_FOR_DATA)

    def get_latest_bed_file(self) -> Optional[str]:
        """Get the path to the latest bed file in the output directory."""
        if not hasattr(self, "output") or not hasattr(self, "sampleID"):
            # print("No output or sampleID found")
            return None

        bed_dir = os.path.join(self.output, "bed_files")
        if not os.path.exists(bed_dir):
            # print(f"Bed directory does not exist: {bed_dir}")
            return None

        bed_files = [
            f
            for f in os.listdir(bed_dir)
            if f.startswith("new_file_") and f.endswith(".bed")
        ]
        if not bed_files:
            # print("No bed files found")
            return None

        # Sort by the numeric part of the filename to get the latest
        latest_file = sorted(
            bed_files, key=lambda x: int(x.split("_")[2].split(".")[0])
        )[-1]
        # print(f"Latest bed file: {latest_file}")
        return os.path.join(bed_dir, latest_file)

    def check_and_update_from_bed_file(self) -> None:
        """Check for new BED files and update the table if needed."""
        current_time = time.time()
        # Only check every 5 seconds to avoid excessive file system access
        if current_time - self.last_bed_check < 5:
            return

        self.last_bed_check = current_time
        latest_bed = self.get_latest_bed_file()

        if not latest_bed or not os.path.exists(latest_bed):
            return
        try:
            # Read the BED file
            bed_data = pd.read_csv(
                latest_bed,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "name", "score", "strand"],
            )

            # Prepare the rows data
            table_rows = []
            for _, row in self.gene_bed.iterrows():
                target_status = "Original"
                target_source = "Panel"

                # Check for overlaps with this gene
                overlaps = bed_data[
                    (bed_data["chrom"] == row.chrom)
                    & (
                        (
                            (bed_data["start"] >= row.start_pos)
                            & (bed_data["start"] <= row.end_pos)
                        )
                        | (
                            (bed_data["end"] >= row.start_pos)
                            & (bed_data["end"] <= row.end_pos)
                        )
                        | (
                            (bed_data["start"] <= row.start_pos)
                            & (bed_data["end"] >= row.end_pos)
                        )
                    )
                ]

                if not overlaps.empty:
                    target_status = "Active"
                    source_name = (
                        overlaps.iloc[0]["name"] if "name" in overlaps.columns else ""
                    )
                    if source_name == "CNV_detected":
                        target_source = "CNV_detected"
                    elif source_name:
                        target_source = source_name

                table_rows.append(
                    {
                        "chrom": str(row.chrom),
                        "gene": str(row.gene),
                        "start_pos": int(row.start_pos),
                        "end_pos": int(row.end_pos),
                        "size": int(row.end_pos - row.start_pos),
                        "status": target_status,
                        "source": target_source,
                    }
                )

            # Update the table if it exists
            if self.target_table:
                # print("Updating table with new rows")
                # Clear existing rows first
                self.target_table.rows = []
                ui.update(self.target_table)
                # Add new rows
                self.target_table.rows = table_rows
                ui.update(self.target_table)
                # print("Table update complete")
            # else:
            # print("No target table found to update")

        except Exception as e:
            logger.error(f"Error reading bed file {latest_bed}: {e}")
            # print(f"Error updating table: {e}")
            if self.target_table:
                self.target_table.rows = []
                ui.update(self.target_table)

    async def stop_analysis(self):
        """Stop the CNV analysis."""
        state.set_process_state("CNV Analysis", ProcessState.STOPPING)
        state.stop_process("CNV Analysis")
        await super().stop_analysis()

    def __del__(self):
        """Cleanup when the object is destroyed."""
        try:
            self._cleanup_state()
        except Exception as e:
            logger.error(f"Error during CNVAnalysis cleanup: {e}")


def test_me(
    port: int,
    threads: int,
    watchfolder: Optional[str],
    output: str,
    reload: bool = False,
    browse: bool = False,
    target_panel: str = "rCNS2",
) -> None:
    """
    Helper function to run the CNV analysis application for testing purposes.

    Args:
        port (int): The port to serve the app on.
        threads (int): Number of threads to use.
        watchfolder (Optional[str]): The directory containing BAM files to process.
        output (str): The directory to store output files.
        reload (bool): Flag to indicate if the app should reload on changes.
        browse (bool): Flag to enable browsing of results.
        target_panel (str): The target gene panel for analysis.
    """
    app.add_static_files("/fonts", str(Path(__file__).parent / "../fonts"))
    with theme.frame("Copy Number Variation Testing."):
        TestObject = CNVAnalysis(
            threads,
            output,
            progress=True,
            target_panel=target_panel,
        )
    if not browse:
        searchdirectory = os.fsencode(watchfolder)
        for root, _, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
    else:
        TestObject.progress_trackers.visible = False
        TestObject.show_previous_data(output)
    ui.run(port=port, reload=reload)


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
@click.option(
    "--target_panel",
    "-t",
    default="rCNS2",
    help="Select analysis gene panel from one of these options. Default is rCNS2",
    type=click.Choice(
        ["rCNS2", "AML"],
        case_sensitive=True,
    ),
)
def main(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    browse: bool,
    target_panel: str,
) -> None:
    """
    Entry point for the command-line interface to start the CNV analysis application.

    Args:
        port (int): The port to serve the app on.
        threads (int): Number of threads to use.
        watchfolder (str): The directory containing BAM files to process.
        output (str): The directory to store output files.
        browse (bool): Flag to enable browsing of results.
        target_panel (str): The target gene panel for analysis.
    """
    if browse:
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=None,
            output=watchfolder,
            browse=browse,
        )
    else:
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=watchfolder,
            output=output,
            browse=browse,
            target_panel=target_panel,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()
