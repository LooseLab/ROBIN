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
from robin.utilities.bed_file import BedTree
import natsort
from robin import theme, resources
import pandas as pd
import logging
import numpy as np
import os
import sys
from nicegui import ui, app, run
import click
from pathlib import Path
import pickle
import ruptures as rpt
from typing import Optional, Tuple, List, BinaryIO
import re

from collections import Counter
from scipy.ndimage import uniform_filter1d

import math
import time

os.environ["CI"] = "1"
# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


def filter_and_find_max(data, window_size=50, z_threshold=2):
    """
    Removes outliers based on a rolling window of local statistics (mean and standard deviation).
    Uses the z-score within each window to classify outliers. Contiguous regions of outliers
    are retained if they consistently deviate in the window. Returns the maximum value
    from the filtered data.

    Args:
        data (np.ndarray): Array of points to process.
        window_size (int): Size of the rolling window for local statistics.
        z_threshold (float): The z-score threshold to classify a point as an outlier.

    Returns:
        float: The maximum value after filtering.
    """
    # Ensure data is a numpy array
    data = np.array(data)

    # Step 1: Convert to pandas Series for rolling calculations
    data_series = pd.Series(data)

    # Step 2: Calculate rolling mean and standard deviation
    rolling_mean = data_series.rolling(window=window_size, center=True).mean()
    rolling_std = data_series.rolling(window=window_size, center=True).std()

    # Fill NaN values for edge cases where rolling window cannot be applied
    rolling_mean.fillna(data_series.mean(), inplace=True)
    rolling_std.fillna(data_series.std(), inplace=True)

    # Step 3: Calculate z-scores based on rolling statistics
    z_scores = np.abs((data_series - rolling_mean) / rolling_std)

    # Step 4: Identify points where the z-score exceeds the threshold
    outlier_mask = z_scores > z_threshold

    # Step 5: Filter out outliers
    filtered_data = data[~outlier_mask]

    # Step 6: Calculate and return the maximum value from the filtered data
    return np.max(filtered_data) if filtered_data.size > 0 else None


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


def iterate_bam(
    bamfile, _threads: int, mapq_filter: int, copy_numbers: dict, log_level: int
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

    return (
        result.cnv,
        result.bin_width,
        result.variance,
        copy_numbers,
        result.genome_length,
    )


def iterate_bam_bin(
    bamfile: BinaryIO,
    _threads: int,
    mapq_filter: int,
    copy_numbers: dict,
    log_level: int,
    bin_width: int,
) -> Tuple[dict, int, float, dict]:
    """
    Iterate over a BAM file with specified bin width and return CNV data and associated metrics.

    Args:
        bamfile (BinaryIO): The BAM file to process.
        _threads (int): Number of threads to use.
        mapq_filter (int): MAPQ filter value.
        copy_numbers (dict): Dictionary to store copy number data.
        log_level (int): Logging level.
        bin_width (int): Bin width for CNV calculation.

    Returns:
        Tuple[dict, int, float, dict]: CNV data, bin width, variance, and updated copy numbers.
    """
    result = iterate_bam_file(
        bamfile,
        _threads=_threads,
        mapq_filter=mapq_filter,
        copy_numbers=copy_numbers,
        log_level=log_level,
        bin_width=bin_width,
    )
    return result.cnv, result.bin_width, result.variance, copy_numbers


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
        NewBed: Optional[BedTree] = None,
        **kwargs,
    ) -> None:
        # self.file_list = []
        self.cnv_dict = {"bin_width": 0, "variance": 0}
        self.update_cnv_dict = {}
        self.result = None
        self.result2 = None
        self.result3 = CNV_Difference()
        self.reference_file = reference_file
        self.bed_file = bed_file
        self.readfish_toml = readfish_toml
        self.XYestimate = "Unknown"
        self.dtype = [("name", "U10"), ("start", "i8"), ("end", "i8")]
        self.DATA_ARRAY = np.empty(0, dtype=self.dtype)

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
        self.target_panel = target_panel
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
        self.NewBed = NewBed
        super().__init__(*args, **kwargs)
        # Only initialize BedTree if reference file is provided

        self.CNVchangedetector = CNVChangeDetectorTracker(base_proportion=0.02)
        # Add target_table as instance variable
        self.target_table = None
        self.target_table_placeholder = None
        self.last_bed_check = 0  # Track when we last checked for new BED files

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

    def detect_chromosome_events(self, z_score_threshold=3.0):
        """Detect significant chromosome-wide CNV events.

        Args:
            z_score_threshold: Number of standard deviations for significance

        Returns:
            List of chromosome-wide CNV events
        """
        if not self.chromosome_stats:
            raise ValueError("Must call calculate_chromosome_stats first")

        events = []
        global_stats = self.chromosome_stats["global"]

        for chrom, stats in self.chromosome_stats.items():
            if chrom == "global":
                continue

            mean_cnv = stats["mean"]
            baseline = stats["baseline"]

            if chrom == "chrX":
                if self.XYestimate == "XY":  # Male
                    if mean_cnv > 0.3:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "GAIN",
                                "description": "Potential XXY",
                            }
                        )
                    elif mean_cnv < -0.3:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "LOSS",
                                "description": "X chromosome loss",
                            }
                        )
                elif self.XYestimate == "XX":  # Female
                    if mean_cnv >= 1.0:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "GAIN",
                                "description": "Potential trisomy X",
                            }
                        )
                    elif mean_cnv < -0.75:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "LOSS",
                                "description": "X chromosome loss",
                            }
                        )
            elif chrom == "chrY":
                if self.XYestimate == "XY":  # Male
                    if mean_cnv > 0.5:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "GAIN",
                                "description": "Y chromosome gain",
                            }
                        )
                    elif mean_cnv < -0.5:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "LOSS",
                                "description": "Y chromosome loss",
                            }
                        )
                elif mean_cnv > -0.2 and self.XYestimate == "XX":
                    events.append(
                        {
                            "chromosome": chrom,
                            "mean_cnv": mean_cnv,
                            "type": "PRESENT",
                            "description": "Y chromosome material detected",
                        }
                    )
            else:  # Autosomes
                adjusted_mean = mean_cnv - baseline
                z_score = (adjusted_mean - global_stats["mean"]) / global_stats["std"]

                if z_score > z_score_threshold:
                    events.append(
                        {
                            "chromosome": chrom,
                            "mean_cnv": mean_cnv,
                            "type": "GAIN",
                            "description": f"Whole chromosome gain (z-score: {z_score:.2f})",
                        }
                    )
                elif z_score < -z_score_threshold:
                    events.append(
                        {
                            "chromosome": chrom,
                            "mean_cnv": mean_cnv,
                            "type": "LOSS",
                            "description": f"Whole chromosome loss (z-score: {z_score:.2f})",
                        }
                    )

        return events

    def estimate_XY(self) -> None:
        """
        Estimate genetic sex (XX or XY) based on CNV data.
        """
        X = round(np.average([i for i in self.result3.cnv["chrX"] if i != 0]), 2)
        Y = round(np.average([i for i in self.result3.cnv["chrY"] if i != 0]), 2)
        if X >= 0.1 and Y <= 0.1:
            self.XYestimate = "XX"
        elif X <= 0.1 and Y >= -0.2:
            self.XYestimate = "XY"
        else:
            self.XYestimate = "Unknown"
        with open(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "XYestimate.pkl",
            ),
            "wb",
        ) as file:
            pickle.dump(self.XYestimate, file)

    async def process_bam(self, bamfile: BinaryIO, timestamp: float) -> None:
        """
        Process a BAM file to extract CNV data.

        Args:
            bamfile (BinaryIO): The BAM file to process.
            timestamp (float): The timestamp indicating when the file was generated.
        """
        # self.file_list.append(bamfile)
        # try:
        await self.do_cnv_work(bamfile)
        # except Exception as e:
        #    logger.error(e)
        #    logger.error("line 313")

    async def do_cnv_work(self, bamfile: BinaryIO) -> None:
        """
        Perform CNV analysis on a BAM file.

        Args:
            bamfile (BinaryIO): The BAM file to process.
        """
        # main_start_time = time.time()
        if self.sampleID not in self.update_cnv_dict.keys():
            self.update_cnv_dict[self.sampleID] = {}

        bamdata = pysam.AlignmentFile(bamfile, "rb")

        self.map_tracker.update(
            Counter(
                {stat.contig: stat.mapped for stat in bamdata.get_index_statistics()}
            )
        )

        r_cnv, r_bin, r_var, self.update_cnv_dict[self.sampleID], genome_length = (
            await run.cpu_bound(
                iterate_bam,
                bamfile,
                self.threads,
                60,
                self.update_cnv_dict[self.sampleID],
                int(logging.ERROR),
            )
        )

        self.cnv_dict["bin_width"] = r_bin
        self.cnv_dict["variance"] = r_var

        r2_cnv, r2_bin, r2_var, self.ref_cnv_dict = await run.cpu_bound(
            iterate_bam_bin,
            bamfile,
            self.threads,
            60,
            self.ref_cnv_dict,
            int(logging.ERROR),
            bin_width=self.cnv_dict["bin_width"],
        )

        if self.load_data:
            if os.path.exists(os.path.join(self.output, self.sampleID, "ruptures.npy")):
                self.CNVchangedetector.coordinates = np.load(
                    os.path.join(self.output, self.sampleID, "ruptures.npy"),
                    allow_pickle="TRUE",
                ).item()
        cnvupdate = False
        for key in r_cnv.keys():
            self.local_data_array = np.empty(0, dtype=self.dtype)
            if key != "chrM" and re.match(r"^chr(\d+|X|Y)$", key):
                moving_avg_data1 = await run.cpu_bound(moving_average, r_cnv[key])
                moving_avg_data2 = await run.cpu_bound(moving_average, r2_cnv[key])
                moving_avg_data1, moving_avg_data2 = await run.cpu_bound(
                    pad_arrays, moving_avg_data1, moving_avg_data2
                )
                self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2
                approx_chrom_length = 0
                if len(r_cnv[key]) > 3:
                    penalty_value = 3
                    if (
                        self.map_tracker[key] > 2000
                    ):  # Make sure we have this number of new reads at least on a chromosome before updating.
                        paired_changepoints = await run.cpu_bound(
                            run_ruptures,
                            r_cnv[key],
                            penalty_value,
                            self.cnv_dict["bin_width"],
                        )
                        if len(paired_changepoints) > 0:
                            approx_chrom_length = len(r_cnv[key]) * r_bin
                            padding = 2_500_000
                            for start, end in paired_changepoints:
                                if (
                                    start < approx_chrom_length < end
                                ):  # Captures change points that overlap with the very end of the chromosome
                                    pass
                                elif (
                                    start < 0
                                ):  # Captures change points that overlap with the very start of the chromosome
                                    pass
                                elif (
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
                                ):  # This ignores any event that may be occuring within a centromere.
                                    pass
                                else:
                                    item = np.array(
                                        [(key, start, end)], dtype=self.dtype
                                    )
                                    self.local_data_array = np.append(
                                        self.local_data_array, item
                                    )
                                    self.DATA_ARRAY = np.append(self.DATA_ARRAY, item)
                        self.map_tracker[key] = 0
                        cnvupdate = True
                        breakpoints = self.DATA_ARRAY[self.DATA_ARRAY["name"] == key]
                        local_breakpoints = self.local_data_array[
                            self.local_data_array["name"] == key
                        ]
                        if len(breakpoints) > 0:
                            try:
                                self.CNVchangedetector.add_breakpoints(
                                    key,
                                    breakpoints,
                                    local_breakpoints,
                                    approx_chrom_length,
                                    r_bin,
                                )
                            except Exception as e:
                                raise (e)

        self.estimate_XY()

        if cnvupdate:
            bedcontent = ""
            bedcontent2 = ""
            self.load_data = True
            for chrom in r_cnv.keys():
                tempbedcontent = self.CNVchangedetector.get_bed_targets(chrom)
                if len(tempbedcontent) > 0:
                    bedcontent += tempbedcontent
                    bedcontent += "\n"

                tempbedcontent2 = self.CNVchangedetector.get_bed_targets_breakpoints(
                    chrom
                )
                if len(tempbedcontent2) > 0:
                    bedcontent2 += tempbedcontent2
                    bedcontent2 += "\n"

            if len(bedcontent2) > 0:
                self.NewBed.load_from_string(
                    bedcontent2,
                    merge=False,
                    write_files=True,
                    output_location=os.path.join(
                        self.check_and_create_folder(self.output, self.sampleID)
                    ),
                    source_type="CNV",
                )
                # Update the target table after loading new BedTree data
                # self.update_target_table()

            np.save(
                os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID),
                    "ruptures.npy",
                ),
                self.CNVchangedetector.coordinates,
            )

        np.save(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID), "CNV.npy"
            ),
            r_cnv,
        )
        np.save(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID), "CNV_dict.npy"
            ),
            self.cnv_dict,
        )

        self.running = False

    def update_plots(self, gene_target: Optional[str] = None) -> None:
        """Update CNV plots with new data and annotations."""
        if not gene_target:
            # Reset zoom when switching to "All" chromosomes view
            if self.chrom_select.value == "All":
                self.scatter_echart.options["dataZoom"][0].update(
                    {"startValue": None, "endValue": None}
                )
                self.difference_scatter_echart.options["dataZoom"][0].update(
                    {"startValue": None, "endValue": None}
                )

            self._update_cnv_plot(
                plot_to_update=self.scatter_echart, result=self.result, title="CNV"
            )
            self._update_cnv_plot(
                plot_to_update=self.reference_scatter_echart,
                result=self.result2,
                title="Reference CNV",
            )
            self._update_cnv_plot(
                plot_to_update=self.difference_scatter_echart,
                result=self.result3,
                title="Difference CNV",
                min="dataMin",
            )
        else:
            self._update_cnv_plot(
                plot_to_update=self.scatter_echart,
                result=self.result,
                gene_target=gene_target,
                title="CNV",
            )
            self._update_cnv_plot(
                plot_to_update=self.reference_scatter_echart,
                result=self.result2,
                gene_target=gene_target,
                title="Reference CNV",
            )
            self._update_cnv_plot(
                plot_to_update=self.difference_scatter_echart,
                result=self.result3,
                gene_target=gene_target,
                title="Difference CNV",
                min="dataMin",
            )

    def _handle_table_row_click(self, e, plot_to_update):
        """Handle click events on table rows by zooming the plot to the selected region."""
        # Only handle clicks when in single chromosome view
        if self.chrom_filter == "All":
            return

        row = e.args["row"]
        start_pos = row["start_pos"]
        end_pos = row["end_pos"]

        # Calculate padding (10% of the region width)
        padding = (end_pos - start_pos) * 0.1

        # Update the plot's x-axis zoom to focus on the selected region
        plot_to_update.options["dataZoom"][0].update(
            {"startValue": max(0, start_pos - padding), "endValue": end_pos + padding}
        )

        # Highlight the selected region more prominently
        for series in plot_to_update.options["series"]:
            if series.get("name") == "cytobands_highlight":
                for area in series["markArea"]["data"]:
                    if area[0]["xAxis"] == start_pos and area[1]["xAxis"] == end_pos:
                        # Temporarily increase opacity for the selected region
                        area[0]["itemStyle"]["opacity"] = 0.5
                    else:
                        area[0]["itemStyle"]["opacity"] = 0.15

        ui.update(plot_to_update)

    def setup_ui(self) -> None:
        """
        Set up the user interface for CNV analysis.
        """

        self.display_row = ui.row().style("width: 100")
        if self.summary:
            with self.summary:
                with ui.card().classes("w-full p-4 mb-4"):
                    with ui.row().classes("w-full items-center justify-between"):
                        # Left side - CNV Status
                        with ui.column().classes("gap-2"):
                            ui.label("Copy Number Variation Analysis").classes(
                                "text-lg font-medium"
                            )
                            with ui.row().classes("items-center gap-2"):
                                ui.label("Status: Awaiting Data").classes(
                                    "text-gray-600"
                                )
                                ui.label("--").classes(
                                    "px-2 py-1 rounded bg-gray-100 text-gray-600"
                                )

                        # Right side - Analysis metrics
                        with ui.column().classes("gap-2 text-right"):
                            ui.label("Analysis Details").classes("font-medium")
                            ui.label("Genetic Sex: --").classes("text-gray-600")
                            ui.label("CNV Events: --").classes("text-gray-600")

                    # Bottom row - Information
                    with ui.row().classes(
                        "w-full mt-4 text-sm text-gray-500 justify-center"
                    ):
                        ui.label(
                            "Copy number analysis across genome with breakpoint detection"
                        )
        with self.display_row:
            ui.label("Copy Number Variation").classes(
                "text-sky-600 dark:text-white"
            ).style("font-size: 150%; font-weight: 300").tailwind(
                "drop-shadow", "font-bold"
            )
            with ui.expansion(
                "Methods",
                caption="A description of the methods used to identify genome-wide CNV events",
            ).classes("w-full"):
                ui.restructured_text(
                    """
                    Copy Number Variation (CNV) analysis is performed through the following steps:

                    1. **Bin-based Coverage Analysis**
                       Reads from BAM files are counted in fixed-width bins across the genome, with 
                       configurable bin sizes and mapping quality filters.

                    2. **Reference Normalization**
                       Sample coverage is compared against a control dataset to identify relative 
                       copy number changes and reduce systematic biases.

                    3. **Change Point Detection**
                       The Kernel Change Point Detection algorithm identifies significant shifts in 
                       copy number profiles, using:
                       * Ruptures library with RBF kernel
                       * Adaptive penalty values
                       * Minimum segment size filtering

                    The analysis maintains strict criteria:
                    * Minimum mapping quality > 60
                    * Bin-level variance tracking
                    * Centromere region masking
                    * Moving average smoothing

                    Each CNV event includes:
                    * Precise genomic coordinates
                    * Mean copy number value
                    * Event classification (Gain/Loss)
                    * Affected genes and cytobands
                    * Statistical confidence metrics

                    Results are displayed with interactive visualizations including:
                    * Chromosome-wide CNV plots
                    * Gene-level zoom capability
                    * Time series tracking of changes
                    * Comprehensive tabular summaries
                """
                ).style("font-size: 100%; font-weight: 300")
        with ui.row():
            self.chrom_select = ui.select(
                options={"All": "All"},
                on_change=self.update_plots,
                label="Select Chromosome",
                value="All",
            ).style("width: 150px")
            self.gene_select = ui.select(
                options={"All": "All"},
                on_change=lambda e: (
                    self.update_plots()
                    if e.value == "All"
                    else self.update_plots(gene_target=e.value)
                ),
                label="Select Gene",
                value="All",
            ).style("width: 150px")
            ui.label().bind_text_from(
                self.cnv_dict, "bin_width", backward=lambda n: f"Bin Width: {n:,}"
            )
            ui.label().bind_text_from(
                self.cnv_dict, "variance", backward=lambda n: f"Variance: {round(n, 3)}"
            )

        self.scatter_echart = self.generate_chart(title="CNV Scatter Plot")
        self.difference_scatter_echart = self.generate_chart(
            title="Difference Plot", initmin=-2, initmax=2
        )

        # Add CNV Analysis Results Table
        with ui.card().classes("w-full p-4"):
            with ui.row().classes("w-full items-center justify-between mb-4"):
                ui.label("CNV Analysis Results").classes("text-lg font-medium")

            # Create the table first
            self.cnv_table = (
                ui.table(
                    columns=[
                        {
                            "name": "chromosome",
                            "label": "Chromosome",
                            "field": "chrom",
                            "align": "left",
                            "sortable": True,
                        },
                        {
                            "name": "start",
                            "label": "Start",
                            "field": "start_pos",
                            "align": "right",
                            ":format": "val => Number(val).toLocaleString()",
                            "sortable": True,
                        },
                        {
                            "name": "end",
                            "label": "End",
                            "field": "end_pos",
                            "align": "right",
                            ":format": "val => Number(val).toLocaleString()",
                            "sortable": True,
                        },
                        {
                            "name": "length",
                            "label": "Length",
                            "field": "length",
                            "align": "right",
                            ":format": "val => Number(val).toLocaleString()",
                            "sortable": True,
                        },
                        {
                            "name": "cytoband",
                            "label": "Cytoband",
                            "field": "name",
                            "align": "left",
                            "sortable": True,
                        },
                        {
                            "name": "mean_cnv",
                            "label": "Mean CNV",
                            "field": "mean_cnv",
                            "align": "right",
                            ":format": "val => Number(val).toFixed(3)",
                            "sortable": True,
                        },
                        {
                            "name": "state",
                            "label": "State",
                            "field": "cnv_state",
                            "align": "center",
                            "sortable": True,
                        },
                        {
                            "name": "genes",
                            "label": "Genes",
                            "field": "genes",
                            "align": "left",
                            "sortable": True,
                        },
                    ],
                    rows=[],
                    row_key="name",
                    pagination=25,
                )
                .classes("w-full")
                .props("dense")
            )

            # Now add the search input after the table is created
            ui.input("Search by chromosome, cytoband or genes...").bind_value(
                self.cnv_table, "filter"
            ).classes("w-64")

            # Add slot for conditional formatting of the CNV state and row styling
            self.cnv_table.add_slot(
                "body",
                """
                <q-tr :props="props" @click="$parent.$emit('row-click', props)">
                    <q-td key="chromosome" :props="props">
                        {{ props.row.chrom }}
                    </q-td>
                    <q-td key="start" :props="props" class="text-right">
                        {{ Number(props.row.start_pos).toLocaleString() }}
                    </q-td>
                    <q-td key="end" :props="props" class="text-right">
                        {{ Number(props.row.end_pos).toLocaleString() }}
                    </q-td>
                    <q-td key="length" :props="props" class="text-right">
                        {{ Number(props.row.length).toLocaleString() }}
                    </q-td>
                    <q-td key="cytoband" :props="props">
                        {{ props.row.name }}
                    </q-td>
                    <q-td key="mean_cnv" :props="props" class="text-right">
                        {{ Number(props.row.mean_cnv).toFixed(3) }}
                    </q-td>
                    <q-td key="state" :props="props">
                        <q-badge :color="props.row.cnv_state === 'GAIN' ? 'positive' : props.row.cnv_state === 'LOSS' ? 'negative' : 'grey'"
                                 :label="props.row.cnv_state"/>
                    </q-td>
                    <q-td key="genes" :props="props">
                        {{ props.row.genes ? props.row.genes.join(', ') : '' }}
                    </q-td>
                </q-tr>
            """,
            )

            # Add event handler for row clicks
            self.cnv_table.on(
                "row-click",
                lambda e: self._handle_table_row_click(e, self.scatter_echart),
            )

        with ui.expansion("See Reference DataSet", icon="loupe").classes("w-full"):
            self.reference_scatter_echart = self.generate_chart(
                title="Reference CNV Scatter Plot"
            )

        with ui.card().classes("w-full"):
            with ui.column().classes("w-full"):
                ui.label("New Target Information").classes("text-lg font-medium")

                with ui.card().classes("w-full"):
                    ui.label(f"Panel: {self.target_panel or 'Not specified'}").classes(
                        "text-gray-600"
                    )
                    # if self.gene_bed_file:
                    #    ui.label(f"Total Targets: {len(self.gene_bed)}").classes("text-gray-600")

                    # Create the table at class level if it doesn't exist
                    if self.target_table is None:
                        # print("Creating new target table")
                        self.target_table = (
                            ui.table(
                                columns=[
                                    {
                                        "name": "chrom",
                                        "label": "Chromosome",
                                        "field": "chrom",
                                        "sortable": True,
                                        "align": "left",
                                    },
                                    {
                                        "name": "gene",
                                        "label": "Gene",
                                        "field": "gene",
                                        "sortable": True,
                                        "align": "left",
                                    },
                                    {
                                        "name": "start",
                                        "label": "Start",
                                        "field": "start_pos",
                                        "sortable": True,
                                        "align": "right",
                                        ":format": "val => Number(val).toLocaleString()",
                                    },
                                    {
                                        "name": "end",
                                        "label": "End",
                                        "field": "end_pos",
                                        "sortable": True,
                                        "align": "right",
                                        ":format": "val => Number(val).toLocaleString()",
                                    },
                                    {
                                        "name": "size",
                                        "label": "Size (bp)",
                                        "field": "size",
                                        "sortable": True,
                                        "align": "right",
                                        ":format": "val => Number(val).toLocaleString()",
                                    },
                                    {
                                        "name": "status",
                                        "label": "Status",
                                        "field": "status",
                                        "sortable": True,
                                        "align": "center",
                                    },
                                    {
                                        "name": "source",
                                        "label": "Source",
                                        "field": "source",
                                        "sortable": True,
                                        "align": "left",
                                    },
                                ],
                                rows=[],
                                row_key="gene",
                                pagination=25,
                            )
                            .classes("w-full")
                            .props("dense rows-per-page-options=[10,25,50,0] filter")
                        )

                        # Add slot for conditional formatting
                        self.target_table.add_slot(
                            "body",
                            """
                            <q-tr :props="props">
                                <q-td key="chrom" :props="props">{{ props.row.chrom }}</q-td>
                                <q-td key="gene" :props="props">{{ props.row.gene }}</q-td>
                                <q-td key="start" :props="props" class="text-right">
                                    {{ Number(props.row.start_pos).toLocaleString() }}
                                </q-td>
                                <q-td key="end" :props="props" class="text-right">
                                    {{ Number(props.row.end_pos).toLocaleString() }}
                                </q-td>
                                <q-td key="size" :props="props" class="text-right">
                                    {{ Number(props.row.size).toLocaleString() }}
                                </q-td>
                                <q-td key="status" :props="props">
                                    <q-badge :color="props.row.status === 'Active' ? 'positive' : 'grey'"
                                                :label="props.row.status"/>
                                </q-td>
                                <q-td key="source" :props="props">
                                    <q-badge :color="props.row.source === 'Panel' ? 'primary' : 
                                                    props.row.source === 'CNV_detected' ? 'warning' : 
                                                    'info'"
                                                :label="props.row.source"/>
                                </q-td>
                            </q-tr>
                            """,
                        )

                        # Add search input after table is created
                        with ui.row().classes("w-full my-2"):
                            ui.input("Search targets...").bind_value_to(
                                self.target_table, "filter"
                            )

                        # print("Target table created successfully")

        with ui.card().classes("w-full"):
            ui.label("Proportion Over Time Information")
            self.create_proportion_time_chart2("Proportions over time - Genome Wide.")
            self.create_proportion_time_chart("Proportions over time.")

        if self.browse:
            ui.timer(0.1, lambda: self.show_previous_data(), once=True)
        else:
            ui.timer(15, lambda: self.show_previous_data())

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
        print("running here")
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

    def create_proportion_time_chart(self, title: str) -> None:
        """
        Creates the NanoDX time series chart.

        Args:
            title (str): Title of the chart.
        """
        self.proportion_time_chart = self.create_time_chart(title)

    def update_proportion_time_chart(self, datadf: pd.DataFrame) -> None:
        """
        Updates the CNV time series chart with new data.

        Parameters
        ----------
        datadf : pd.DataFrame
            DataFrame containing the time series data to plot.

        Notes
        -----
        Handles empty datasets gracefully by displaying an appropriate message
        and maintaining chart visibility.
        """
        try:
            if datadf.empty:
                logger.warning("No data available for time series chart")
                self.proportion_time_chart.options["title"][
                    "text"
                ] = "No Data Available"
                self.proportion_time_chart.options["series"] = []
                self.proportion_time_chart.update()
                return

            self.proportion_time_chart.options["series"] = []

            # iOS color palette for consistent styling
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
                    data_list = [
                        [key, value] for key, value in data.items() if pd.notnull(value)
                    ]
                    if data_list:  # Only add series if it has valid data points
                        self.proportion_time_chart.options["series"].append(
                            {
                                "animation": False,
                                "type": "line",
                                "smooth": True,
                                "name": series,
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
                                "lineStyle": {
                                    "width": 2,
                                    "color": colors[idx % len(colors)],
                                },
                                "itemStyle": {"color": colors[idx % len(colors)]},
                                "data": data_list,
                            }
                        )

            # Update title with latest data summary if available
            if not datadf.empty:
                latest_data = datadf.iloc[-1]
                if not latest_data.empty and not latest_data.isna().all():
                    max_value = latest_data.max()
                    max_type = latest_data.idxmax()
                    self.proportion_time_chart.options["title"]["text"] = (
                        f"Proportions Over Time\n"
                        f"Current highest: {max_type} ({max_value:.1f}%)"
                    )

            self.proportion_time_chart.update()

        except Exception as e:
            logger.error(f"Error updating time chart: {str(e)}", exc_info=True)
            self.proportion_time_chart.options["title"]["text"] = "Error Updating Chart"
            self.proportion_time_chart.options["series"] = []
            self.proportion_time_chart.update()

    def create_proportion_time_chart2(self, title: str) -> None:
        """
        Creates the genome-wide CNV time series chart.

        Parameters
        ----------
        title : str
            Title of the chart.

        Notes
        -----
        Creates a time series chart specifically for genome-wide data,
        using consistent styling with other charts in the application.
        """
        self.proportion_time_chart2 = self.create_time_chart(title)

    def update_proportion_time_chart2(self, datadf: pd.DataFrame) -> None:
        """
        Updates the genome-wide CNV time series chart with new data.

        Parameters
        ----------
        datadf : pd.DataFrame
            DataFrame containing the genome-wide time series data to plot.

        Notes
        -----
        Handles empty datasets gracefully by displaying an appropriate message
        and maintaining chart visibility.
        """
        try:
            if datadf.empty:
                logger.warning("No genome-wide data available for time series chart")
                self.proportion_time_chart2.options["title"][
                    "text"
                ] = "No Data Available"
                self.proportion_time_chart2.options["series"] = []
                self.proportion_time_chart2.update()
                return

            self.proportion_time_chart2.options["series"] = []

            # iOS color palette for consistent styling
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
                    data_list = [
                        [key, value] for key, value in data.items() if pd.notnull(value)
                    ]
                    if data_list:  # Only add series if it has valid data points
                        self.proportion_time_chart2.options["series"].append(
                            {
                                "animation": False,
                                "type": "line",
                                "smooth": True,
                                "name": series,
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
                                "lineStyle": {
                                    "width": 2,
                                    "color": colors[idx % len(colors)],
                                },
                                "itemStyle": {"color": colors[idx % len(colors)]},
                                "data": data_list,
                            }
                        )

            # Update title with latest data summary if available
            if not datadf.empty:
                latest_data = datadf.iloc[-1]
                if not latest_data.empty and not latest_data.isna().all():
                    max_value = latest_data.max()
                    max_type = latest_data.idxmax()
                    self.proportion_time_chart2.options["title"]["text"] = (
                        f"Genome-wide Proportions Over Time\n"
                        f"Current highest: {max_type} ({max_value:.1f}%)"
                    )

            self.proportion_time_chart2.update()

        except Exception as e:
            logger.error(
                f"Error updating genome-wide time chart: {str(e)}", exc_info=True
            )
            self.proportion_time_chart2.options["title"][
                "text"
            ] = "Error Updating Chart"
            self.proportion_time_chart2.options["series"] = []
            self.proportion_time_chart2.update()

    def generate_chart(
        self,
        title: Optional[str] = None,
        initmax: Optional[int] = None,
        initmin: int = 0,
        type: str = "value",
    ) -> ui.echart:
        """
        Generate an ECharts object for displaying CNV scatter plots.
        """
        return (
            ui.echart(
                {
                    "backgroundColor": "transparent",
                    "textStyle": {
                        "fontFamily": "SF Pro Text, -apple-system, BlinkMacSystemFont, Helvetica, Arial, sans-serif",
                        "fontSize": 12,
                    },
                    "animation": False,
                    "grid": {
                        "top": "25%",
                        "bottom": "20%",
                        "left": "5%",
                        "right": "10%",
                        "containLabel": True,
                    },
                    "title": {
                        "text": f"{title}",
                        "left": "center",
                        "top": 20,
                        "textStyle": {
                            "fontSize": 16,
                            "fontWeight": "500",
                            "color": "#1D1D1F",
                        },
                    },
                    "toolbox": {
                        "show": True,
                        "right": 20,
                        "feature": {"saveAsImage": {"title": "Save Chart"}},
                    },
                    "tooltip": {
                        "trigger": "axis",
                        "backgroundColor": "rgba(255, 255, 255, 0.9)",
                        "borderColor": "#E5E5EA",
                        "textStyle": {"color": "#1D1D1F"},
                        "axisPointer": {"type": "cross"},
                        ":formatter": "params => `Position: ${params[0].value[0]}<br/>Value: ${params[0].value[1].toFixed(2)}`",
                    },
                    "xAxis": {
                        "type": f"{type}",
                        "max": "dataMax",
                        "splitLine": {
                            "show": True,
                            "lineStyle": {"type": "dashed", "color": "#E5E5EA"},
                        },
                        "axisLabel": {"color": "#86868B"},
                    },
                    "yAxis": [
                        {
                            "type": "value",
                            "logBase": 2,
                            "name": "Ploidy",
                            "position": "left",
                            "nameTextStyle": {
                                "color": "#86868B",
                                "fontSize": 12,
                                "padding": [0, 30, 0, 0],
                            },
                            "axisLabel": {"color": "#86868B"},
                            "splitLine": {
                                "show": True,
                                "lineStyle": {"type": "dashed", "color": "#E5E5EA"},
                            },
                        },
                        {
                            "type": "value",
                            "name": "Candidate Break Points",
                            "position": "right",
                            "nameTextStyle": {
                                "color": "#007AFF",
                                "fontSize": 12,
                                "padding": [0, 0, 0, 30],
                            },
                            "axisLine": {"lineStyle": {"color": "#007AFF"}},
                            "axisLabel": {"color": "#007AFF"},
                            "splitLine": {"show": False},
                        },
                    ],
                    "dataZoom": [
                        {
                            "type": "slider",
                            "yAxisIndex": "none",
                            "xAxisIndex": "0",
                            "filterMode": "none",
                            "height": 20,
                            "bottom": 35,
                            "borderColor": "#E5E5EA",
                            "backgroundColor": "#F5F5F7",
                            "fillerColor": "rgba(0, 122, 255, 0.2)",
                            "handleStyle": {"color": "#007AFF"},
                        },
                        {
                            "type": "slider",
                            "yAxisIndex": "0",
                            "xAxisIndex": "none",
                            "filterMode": "none",
                            "width": 20,
                            "right": 50,
                            "startValue": initmin,
                            "endValue": initmax,
                            "borderColor": "#E5E5EA",
                            "backgroundColor": "#F5F5F7",
                            "fillerColor": "rgba(0, 122, 255, 0.2)",
                            "handleStyle": {"color": "#007AFF"},
                        },
                    ],
                    "series": [
                        {
                            "type": "scatter",
                            "symbolSize": 4,
                            "itemStyle": {"color": "#007AFF", "opacity": 0.7},
                            "emphasis": {
                                "itemStyle": {
                                    "color": "#007AFF",
                                    "opacity": 1,
                                    "borderColor": "#FFFFFF",
                                    "borderWidth": 2,
                                }
                            },
                            "data": [],
                        },
                        {
                            "type": "scatter",
                            "name": "centromeres_highlight",
                            "data": [],  # Empty data because this series is just for highlighting
                            "symbolSize": 3,
                            "markArea": {
                                "itemStyle": {
                                    "color": "rgba(135, 206, 250, 0.4)"  # Light blue color
                                },
                                "data": [],
                            },
                        },
                        {
                            "type": "scatter",
                            "name": "cytobands_highlight",
                            "data": [],  # Empty data because this series is just for highlighting
                            "symbolSize": 3,
                            "markArea": {
                                "itemStyle": {
                                    "color": "rgba(200, 200, 200, 0.4)"  # Light gray color for default
                                },
                                "data": [],
                            },
                            "markLine": {"symbol": "none", "data": []},
                        },
                    ],
                }
            )
            .style("height: 450px; margin: 20px 0;")
            .classes("border rounded-lg shadow-sm")
        )

    def _update_cnv_plot(
        self,
        plot_to_update: ui.echart,
        result: Optional[Result] = None,
        gene_target: Optional[str] = None,
        title: Optional[str] = None,
        min: Optional[float] = 0,
        ui_mode: bool = True,
    ) -> None:
        """
        Update CNV plots with new data and annotations.

        Args:
            plot_to_update (ui.echart): The ECharts object to update.
            result (Optional[Result]): The CNV result data.
            gene_target (Optional[str]): Target gene for updating plots.
            title (Optional[str]): Title of the plot.
            min (Optional[float]): Minimum value for the x-axis.
            ui_mode (bool): Flag to indicate if in UI mode.
        """

        # Check if there is CNV result data available, either passed to the function or stored in the instance
        if result or self.result:
            # Initialize variables for data accumulation and dropdown options
            total = 0  # Accumulates the total length of the CNV data
            valueslist = {
                "All": "All"
            }  # Dictionary to store chromosome dropdown options
            genevalueslist = {"All": "All"}  # Dictionary to store gene dropdown options

            # Try to set the chromosome filter from the UI selection, default to "All" if not available
            try:
                self.chrom_filter = self.chrom_select.value
            except AttributeError:
                self.chrom_filter = "All"

            # Set initial min and max values for the x-axis
            min = min  # This can be modified later based on gene targeting
            max = "dataMax"  # Placeholder to automatically determine the maximum value

            # If a specific gene is targeted, adjust the plot to focus on the region around that gene
            if gene_target:
                # Retrieve the start and end positions and chromosome of the targeted gene
                start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                end_pos = self.gene_bed.iloc[int(gene_target)].end_pos
                chrom = self.gene_bed.iloc[int(gene_target)].chrom

                # Loop through the sorted CNV data to find the chromosome of the targeted gene
                for counter, contig in enumerate(
                    natsort.natsorted(result.cnv), start=1
                ):
                    valueslist[counter] = (
                        contig  # Add chromosome to the dropdown options
                    )
                    if contig == chrom:  # Stop when the target chromosome is found
                        break

                # Set the chromosome filter to the index of the target chromosome
                self.chrom_filter = counter

                # Adjust the x-axis min and max to focus on the region around the gene target
                min = start_pos - 10 * self.cnv_dict["bin_width"]
                max = end_pos + 10 * self.cnv_dict["bin_width"]

                # Further adjust the axis limits if the gene region is large
                if start_pos - min > 2_000_000:
                    min = start_pos - 2_000_000
                if max - end_pos > 2_000_000:
                    max = end_pos + 2_000_000

                if min < 0:
                    min = 0  # Ensure the minimum x-axis value is not negative

            # If all chromosomes are selected, prepare the plot to display all chromosomes
            if self.chrom_filter == "All":
                counter = 0  # Reset counter for chromosome loop

                # Update the plot title to reflect that all chromosomes are being shown
                plot_to_update.options["title"]["text"] = f"{title} - All Chromosomes"
                plot_to_update.options["series"] = []

                # Collect CNV analysis for all chromosomes
                all_cytoband_analysis = []
                for contig, cnv in natsort.natsorted(result.cnv.items()):
                    # Skip non-standard chromosomes
                    if contig == "chrM" or not re.match(r"^chr(\d+|X|Y)$", contig):
                        continue

                    # Get CNV analysis for this chromosome
                    chromosome_analysis = self.analyze_cytoband_cnv(
                        self.result3.cnv, contig
                    )
                    if not chromosome_analysis.empty:
                        all_cytoband_analysis.append(chromosome_analysis)

                # Update the table with combined analysis results
                try:
                    combined_analysis = (
                        pd.concat(all_cytoband_analysis)
                        if all_cytoband_analysis
                        else pd.DataFrame()
                    )
                    if (
                        combined_analysis.empty
                        and self.cnv_dict["bin_width"] > 10_000_000
                    ):
                        # Create a single row DataFrame with a message
                        message_df = pd.DataFrame(
                            [
                                {
                                    "chrom": "",
                                    "start_pos": 0,
                                    "end_pos": 0,
                                    "name": "More data needed for CNV analysis",
                                    "mean_cnv": 0,
                                    "cnv_state": f'Current bin width: {self.cnv_dict["bin_width"]:,}bp (need ≤ 10,000,000bp)',
                                    "length": 0,
                                }
                            ]
                        )
                        self.cnv_table.rows = message_df.to_dict("records")
                    else:
                        self.cnv_table.rows = combined_analysis.to_dict("records")
                    ui.update(self.cnv_table)
                except AttributeError:
                    pass  # Ignore if the table component is not available

                # Loop through the sorted CNV data and plot each chromosome
                for contig, cnv in natsort.natsorted(result.cnv.items()):
                    # Skip non-standard chromosomes (e.g., mitochondrial DNA or unrecognized chromosomes)
                    if contig == "chrM" or not re.match(r"^chr(\d+|X|Y)$", contig):
                        continue

                    counter += 1  # Increment chromosome counter
                    valueslist[counter] = contig  # Add chromosome to dropdown options

                    # Prepare data for plotting, adjusting positions by bin width
                    data = list(
                        zip(
                            (np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"],
                            cnv,
                        )
                    )

                    data = reduce_list(
                        data
                    )  # Reduce the data points for more efficient plotting
                    total += len(
                        cnv
                    )  # Update total length with current chromosome data length

                    # Set the x-axis min and max values
                    plot_to_update.options["xAxis"]["max"] = max
                    plot_to_update.options["xAxis"]["min"] = min

                    # Append the current chromosome data as a scatter plot series
                    plot_to_update.options["series"].append(
                        {
                            "type": "scatter",
                            "name": contig,
                            "data": data,
                            "symbolSize": 5,  # Size of the scatter plot points
                            "markLine": {
                                "symbol": "none",
                                "data": [
                                    {
                                        "lineStyle": {"width": 1},
                                        "label": {"formatter": contig},
                                        "name": contig,
                                        "xAxis": (
                                            (total - len(cnv) / 2)
                                            * self.cnv_dict["bin_width"]
                                        ),
                                    },
                                    {
                                        "lineStyle": {"width": 2},
                                        "label": {"normal": {"show": False}},
                                        "xAxis": ((total) * self.cnv_dict["bin_width"]),
                                    },
                                ],
                            },
                        }
                    )

                    # Add gene information to the dropdown options for the current chromosome
                    for index, gene in self.gene_bed[
                        self.gene_bed["chrom"] == contig
                    ].iterrows():
                        genevalueslist[index] = f"{gene.chrom} - {gene.gene}"

            # If a specific chromosome is selected, plot only that chromosome
            else:
                # Update dropdown options for chromosomes
                for counter, contig in enumerate(
                    natsort.natsorted(result.cnv), start=1
                ):
                    valueslist[counter] = contig

                # Define the main chromosomes to be included in the plot
                main_chromosomes = [
                    "chr1",
                    "chr2",
                    "chr3",
                    "chr4",
                    "chr5",
                    "chr6",
                    "chr7",
                    "chr8",
                    "chr9",
                    "chr10",
                    "chr11",
                    "chr12",
                    "chr13",
                    "chr14",
                    "chr15",
                    "chr16",
                    "chr17",
                    "chr18",
                    "chr19",
                    "chr20",
                    "chr21",
                    "chr22",
                    "chrX",
                    "chrY",
                ]

                # Filter CNV data to include only the main chromosomes
                filtered_data = {
                    k: v for k, v in result.cnv.items() if k in main_chromosomes
                }

                # Select the data for the filtered chromosome based on the user's selection
                contig, cnv = natsort.natsorted(filtered_data.items())[
                    int(self.chrom_filter) - 1
                ]

                # Prepare the data for plotting
                data = list(
                    zip((np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"], cnv)
                )

                ymax = math.ceil(filter_and_find_max(np.array(cnv)))

                if not gene_target:
                    min = min
                    max = "dataMax"
                else:
                    # If the gene target is "All", adjust axis limits to show the entire chromosome
                    if gene_target == "All":
                        min = min
                        max = "dataMax"
                    else:
                        # For a specific gene, adjust the axis limits to zoom in on the gene
                        start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                        end_pos = self.gene_bed.iloc[int(gene_target)].end_pos

                        chrom = self.gene_bed.iloc[int(gene_target)].chrom
                        counter = 0
                        for counter, contig in enumerate(
                            natsort.natsorted(result.cnv), start=1
                        ):
                            if contig == chrom:
                                self.chrom_filter = counter
                                break

                        # Set axis limits to focus on the gene, with a margin of 10 times the bin width
                        min = start_pos - 10 * self.cnv_dict["bin_width"]
                        max = end_pos + 10 * self.cnv_dict["bin_width"]

                        # Further adjust the axis limits if the gene region is large
                        if start_pos - min > 2_000_000:
                            min = start_pos - 2_000_000
                        if max - end_pos > 2_000_000:
                            max = end_pos + 2_000_000

                        if min < 0:
                            min = 0  # Ensure the minimum x-axis value is not negative

                # Update the plot title to reflect the selected chromosome
                if "Difference" in title:
                    plot_to_update.options["title"][
                        "text"
                    ] = f"Copy Number Variation (Relative Difference) - {contig}"
                else:
                    plot_to_update.options["title"][
                        "text"
                    ] = f"Copy Number Variation (Absolute) - {contig}"

                plot_to_update.options["xAxis"]["max"] = max
                plot_to_update.options["xAxis"]["min"] = min
                plot_to_update.options["dataZoom"][1]["endValue"] = ymax

                # Now initialize the series after we have contig defined
                plot_to_update.options["series"] = [
                    {
                        "type": "scatter",
                        "name": contig,
                        "data": data,
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(255, 173, 177, 0.4)"},
                            "data": [],
                        },
                    },
                    {
                        "type": "scatter",
                        "name": "centromeres_highlight",
                        "data": [],
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},
                            "data": [],
                        },
                    },
                    {
                        "type": "scatter",
                        "name": "cytobands_highlight",
                        "data": [],
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(200, 200, 200, 0.4)"},
                            "data": [],
                        },
                        "markLine": {"symbol": "none", "data": []},
                    },
                ]

                # Add gene information to the dropdown options and highlight gene regions in the plot
                for index, gene in self.gene_bed[
                    self.gene_bed["chrom"] == contig
                ].iterrows():
                    genevalueslist[index] = f"{gene.chrom} - {gene.gene}"

                # Highlight the regions of the genes in the plot with shaded areas
                for _, row in self.gene_bed[
                    self.gene_bed["chrom"] == contig
                ].iterrows():
                    plot_to_update.options["series"][0]["markArea"]["data"].append(
                        [
                            {
                                "name": row["gene"],
                                "xAxis": row["start_pos"],
                                "label": {
                                    "position": "insideTop",
                                    "distance": 0,  # Reduced from 25 to 15
                                    "color": "#000000",
                                    "fontSize": 12,
                                    "show": True,
                                    "emphasis": {
                                        "show": True,
                                        "distance": 0,  # Reduced from 40 to 25
                                    },
                                },
                            },
                            {
                                "xAxis": row["end_pos"],
                            },
                        ]
                    )

                # Highlight the centromeres with shaded area:
                for _, row in self.centromere_bed[
                    self.centromere_bed["chrom"] == contig
                ].iterrows():
                    plot_to_update.options["series"][1]["markArea"]["data"].append(
                        [
                            {
                                "name": row["name"],
                                "xAxis": row["start_pos"],
                                "label": {"position": "insideBottom", "distance": 10},
                            },
                            {
                                "xAxis": row["end_pos"],
                            },
                        ]
                    )

                # Add cytoband highlighting with CNV state colors
                for _, row in self.cytobands_bed[
                    self.cytobands_bed["chrom"] == contig
                ].iterrows():
                    # Get CNV state for this region from analysis
                    region_analysis = self.analyze_cytoband_cnv(
                        self.result3.cnv, contig
                    )
                    matching_region = region_analysis[
                        (region_analysis["start_pos"] <= row["start_pos"])
                        & (region_analysis["end_pos"] >= row["end_pos"])
                    ]

                    # Set color based on CNV state
                    if not matching_region.empty:
                        cnv_state = matching_region.iloc[0]["cnv_state"]
                        if cnv_state == "GAIN":
                            color = "rgba(52, 199, 89, 0.15)"  # Green for gains
                        elif cnv_state == "LOSS":
                            color = "rgba(255, 45, 85, 0.15)"  # Red for losses
                        else:
                            color = "rgba(0, 0, 0, 0.02)"  # Very subtle gray for normal
                    else:
                        color = "rgba(0, 0, 0, 0.02)"  # Default subtle gray

                    # Add colored regions to the plot
                    plot_to_update.options["series"][2]["markArea"]["data"].append(
                        [
                            {
                                "name": row[
                                    "name"
                                ],  # Just the cytoband name without state
                                "xAxis": row["start_pos"],
                                "itemStyle": {"color": color},
                                "label": {
                                    "position": "insideTop",
                                    "distance": 25,
                                    "color": cnv_state == "GAIN"
                                    and "#34C759"  # Green for gains
                                    or cnv_state == "LOSS"
                                    and "#E0162B"  # Brighter, more saturated red for losses
                                    or "#000000",  # Black for normal
                                    "fontWeight": cnv_state in ["GAIN", "LOSS"]
                                    and "bold"
                                    or "normal",
                                    "show": True,
                                },
                            },
                            {
                                "xAxis": row["end_pos"],
                            },
                        ]
                    )

                    # Add horizontal reference lines for CNV thresholds only on the difference plot
                    if "Difference" in plot_to_update.options["title"]["text"]:
                        plot_to_update.options["series"][2]["markLine"]["data"].extend(
                            [
                                {
                                    "name": "Gain Threshold",
                                    "yAxis": 0.5,
                                    "lineStyle": {
                                        "color": "rgba(52, 199, 89, 0.5)",
                                        "type": "dashed",
                                    },
                                    "label": {
                                        "show": True,
                                        "formatter": "Gain Threshold",
                                        "position": "insideEndTop",
                                    },
                                },
                                {
                                    "name": "Loss Threshold",
                                    "yAxis": -0.5,
                                    "lineStyle": {
                                        "color": "rgba(255, 45, 85, 0.5)",
                                        "type": "dashed",
                                    },
                                    "label": {
                                        "show": True,
                                        "formatter": "Loss Threshold",
                                        "position": "insideEndBottom",
                                    },
                                },
                            ]
                        )

            # Update the chromosome dropdown options in the UI
            try:
                self.chrom_select.set_options(valueslist)
            except AttributeError:
                pass  # Ignore if the UI component is not available

            # Update the gene dropdown options in the UI
            try:
                self.gene_select.set_options(genevalueslist)
            except AttributeError:
                pass  # Ignore if the UI component is not available

            # If in UI mode, update the plot in the UI
            if ui_mode:
                ui.update(plot_to_update)
            else:
                return (
                    plot_to_update  # If not in UI mode, return the updated plot object
                )

            # If a specific chromosome is selected, update the CNV analysis table
            if self.chrom_filter != "All":
                contig = valueslist[int(self.chrom_filter)]
                # Get CNV analysis for the current chromosome
                cytoband_analysis = self.analyze_cytoband_cnv(self.result3.cnv, contig)
                # Update the table with the analysis results
                try:
                    if (
                        cytoband_analysis.empty
                        and self.cnv_dict["bin_width"] > 10_000_000
                    ):
                        # Create a single row DataFrame with a message
                        message_df = pd.DataFrame(
                            [
                                {
                                    "chrom": contig,
                                    "start_pos": 0,
                                    "end_pos": 0,
                                    "name": "More data needed for CNV analysis",
                                    "mean_cnv": 0,
                                    "cnv_state": f'Current bin width: {self.cnv_dict["bin_width"]:,}bp (need ≤ 10,000,000bp)',
                                    "length": 0,
                                }
                            ]
                        )
                        self.cnv_table.rows = message_df.to_dict("records")
                    else:
                        self.cnv_table.rows = cytoband_analysis.to_dict("records")
                    ui.update(self.cnv_table)
                except AttributeError:
                    pass  # Ignore if the table component is not available

            # Add click event handler for the table rows
            self.cnv_table.on(
                "row-click", lambda e: self._handle_table_row_click(e, plot_to_update)
            )

    def create_summary_card(
        self, xy_estimate: str, bin_width: int, variance: float
    ) -> None:
        """
        Create a summary card displaying key CNV analysis information.

        Creates a visually consistent summary card that matches other analysis modules,
        displaying genetic sex estimation and key CNV metrics.

        Parameters
        ----------
        xy_estimate : str
            Estimated genetic sex (XX, XY, or Unknown)
        bin_width : int
            Current bin width used in the analysis
        variance : float
            Current variance value from the analysis
        """
        with self.summary:
            self.summary.clear()
            with ui.card().classes("w-full p-4 mb-4"):
                with ui.row().classes("w-full items-center justify-between"):
                    # Left side - CNV Status
                    with ui.column().classes("gap-2"):
                        ui.label("Copy Number Variation Analysis").classes(
                            "text-lg font-medium"
                        )
                        with ui.row().classes("items-center gap-2"):
                            if xy_estimate != "Unknown":
                                status_color = (
                                    "text-blue-600"
                                    if xy_estimate == "XY"
                                    else "text-pink-600"
                                )
                                status_bg = (
                                    "bg-blue-100"
                                    if xy_estimate == "XY"
                                    else "bg-pink-100"
                                )
                                if xy_estimate == "XY":
                                    ui.icon("man").classes("text-4xl text-blue-500")
                                else:
                                    ui.icon("woman").classes("text-4xl text-pink-500")
                                ui.label(f"Genetic Sex: {xy_estimate}").classes(
                                    f"{status_color} font-medium"
                                )
                            else:
                                ui.label("Status: Awaiting Data").classes(
                                    "text-gray-600"
                                )
                                ui.label("--").classes(
                                    "px-2 py-1 rounded bg-gray-100 text-gray-600"
                                )

                    # Right side - Analysis metrics
                    with ui.column().classes("gap-2 text-right"):
                        ui.label("Analysis Details").classes("font-medium")
                        ui.label(f"Bin Width: {bin_width:,}").classes("text-gray-600")
                        ui.label(f"Variance: {variance:.3f}").classes("text-gray-600")

                        # Calculate gene counts for gains and losses
                        total_gained_genes = set()
                        total_lost_genes = set()

                        if hasattr(self, "result3") and self.result3.cnv:
                            main_chromosomes = [f"chr{i}" for i in range(1, 23)] + [
                                "chrX",
                                "chrY",
                            ]
                            for chrom in main_chromosomes:
                                if chrom in self.result3.cnv:
                                    analysis = self.analyze_cytoband_cnv(
                                        self.result3.cnv, chrom
                                    )
                                    if not analysis.empty:
                                        # Get genes in gained regions
                                        gained = analysis[
                                            analysis["cnv_state"] == "GAIN"
                                        ]
                                        for _, row in gained.iterrows():
                                            if row["genes"]:
                                                total_gained_genes.update(row["genes"])

                                        # Get genes in lost regions
                                        lost = analysis[analysis["cnv_state"] == "LOSS"]
                                        for _, row in lost.iterrows():
                                            if row["genes"]:
                                                total_lost_genes.update(row["genes"])

                        # Display gene counts in the right column
                        with ui.row().classes("gap-2 justify-end mt-2"):
                            with ui.card().classes("py-1 px-2 bg-green-50 rounded"):
                                ui.label(f"Gained: {len(total_gained_genes)}").classes(
                                    "text-sm text-green-800"
                                )
                            with ui.card().classes("py-1 px-2 bg-red-50 rounded"):
                                ui.label(f"Lost: {len(total_lost_genes)}").classes(
                                    "text-sm text-red-800"
                                )

                # Bottom row - Information
                with ui.row().classes(
                    "w-full mt-4 text-sm text-gray-500 justify-center"
                ):
                    ui.label(
                        "Copy number analysis across genome with breakpoint detection"
                    )

    async def show_previous_data(self) -> None:
        """
        Load and display previously computed CNV data.

        Loads saved CNV analysis results and updates the visualization components
        with the data. Creates a consistent visual presentation that matches
        other analysis modules.

        Notes
        -----
        - Loads CNV data from saved numpy files
        - Updates scatter plots and time series charts
        - Creates summary cards with key metrics
        - Maintains visual consistency with other modules
        """
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        if self.check_file_time(os.path.join(output, "CNV.npy")):
            result = np.load(
                os.path.join(output, "CNV.npy"), allow_pickle="TRUE"
            ).item()
            self.result = Result(result)
            cnv_dict = np.load(
                os.path.join(output, "CNV_dict.npy"), allow_pickle=True
            ).item()
            self.cnv_dict["bin_width"] = cnv_dict["bin_width"]
            self.cnv_dict["variance"] = cnv_dict["variance"]

            r2_cnv, r2_bin, r2_var, self.ref_cnv_dict = await run.cpu_bound(
                iterate_bam_bin,
                None,
                1,
                60,
                self.ref_cnv_dict,
                int(logging.ERROR),
                bin_width=self.cnv_dict["bin_width"],
            )

            self.result2 = Result(r2_cnv)

            for key in self.result.cnv.keys():
                if key != "chrM" and re.match(r"^chr(\d+|X|Y)$", key):
                    moving_avg_data1 = moving_average(self.result.cnv[key])
                    moving_avg_data2 = moving_average(r2_cnv[key])
                    moving_avg_data1, moving_avg_data2 = pad_arrays(
                        moving_avg_data1, moving_avg_data2
                    )
                    self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2

            self.CNVResults = {}
            if self.check_file_time(os.path.join(output, "ruptures.npy")):
                self.CNVResults = np.load(
                    os.path.join(output, "ruptures.npy"), allow_pickle="TRUE"
                ).item()

            self.update_plots()

            bedcontent = ""
            local_update = False
            for chrom in natsort.natsorted(self.CNVResults):
                if len(self.CNVResults[chrom]["bed_data"]) > 0:
                    bedcontent += f'{self.CNVResults[chrom]["bed_data_breakpoints"]}\n'
                    local_update = True

            self.NewBed.load_from_string(bedcontent, merge=False, source_type="CNV")

            if self.check_file_time(os.path.join(output, "bedranges.csv")):
                self.proportions_df_store = pd.read_csv(
                    os.path.join(os.path.join(output, "bedranges.csv")),
                    index_col=0,
                )

                df_filtered = self.proportions_df_store[
                    self.proportions_df_store["chromosome"] != "Genome-wide"
                ]
                pivot_df = df_filtered.pivot(
                    index="timestamp", columns="chromosome", values="proportion"
                )

                df_filtered2 = self.proportions_df_store[
                    self.proportions_df_store["chromosome"] == "Genome-wide"
                ]
                pivot_df2 = df_filtered2.pivot(
                    index="timestamp", columns="chromosome", values="proportion"
                )

                self.update_proportion_time_chart(pivot_df)
                self.update_proportion_time_chart2(pivot_df2)

            self.update_target_table()

            if self.summary:
                with self.summary:
                    self.summary.clear()
                    with open(os.path.join(output, "XYestimate.pkl"), "rb") as file:
                        xy_estimate = pickle.load(file)
                        self.create_summary_card(
                            xy_estimate=xy_estimate,
                            bin_width=self.cnv_dict["bin_width"],
                            variance=self.cnv_dict["variance"],
                        )

    def create_time_chart(self, title: str) -> ui.echart:
        """
        Create a time series chart for CNV data visualization.
        """
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

        return (
            ui.echart(
                {
                    "backgroundColor": "transparent",
                    "textStyle": {
                        "fontFamily": "SF Pro Text, -apple-system, BlinkMacSystemFont, Helvetica, Arial, sans-serif",
                        "fontSize": 12,
                    },
                    "animation": False,
                    "grid": {
                        "top": "25%",
                        "bottom": "20%",
                        "left": "5%",
                        "right": "10%",
                        "containLabel": True,
                    },
                    "title": {
                        "text": title,
                        "left": "center",
                        "top": 20,
                        "textStyle": {
                            "fontSize": 16,
                            "fontWeight": "500",
                            "color": "#1D1D1F",
                        },
                    },
                    "tooltip": {
                        "trigger": "axis",
                        "backgroundColor": "rgba(255, 255, 255, 0.9)",
                        "borderColor": "#E5E5EA",
                        "textStyle": {"color": "#1D1D1F"},
                        "axisPointer": {"type": "cross"},
                        ":formatter": "params => `Position: ${params[0].value[0]}<br/>Value: ${params[0].value[1].toFixed(2)}`",
                    },
                    "legend": {
                        "type": "scroll",
                        "top": 50,
                        "textStyle": {"color": "#86868B"},
                    },
                    "toolbox": {
                        "show": True,
                        "right": 20,
                        "feature": {"saveAsImage": {"title": "Save Chart"}},
                    },
                    "xAxis": {
                        "type": "time",
                        "splitLine": {
                            "show": True,
                            "lineStyle": {"type": "dashed", "color": "#E5E5EA"},
                        },
                        "axisLabel": {"color": "#86868B"},
                    },
                    "yAxis": {
                        "type": "value",
                        "name": "Proportion",
                        "nameTextStyle": {
                            "color": "#86868B",
                            "fontSize": 12,
                            "padding": [0, 30, 0, 0],
                        },
                        "axisLabel": {"color": "#86868B", "formatter": "{value}%"},
                        "splitLine": {
                            "show": True,
                            "lineStyle": {"type": "dashed", "color": "#E5E5EA"},
                        },
                    },
                    "color": colors,
                    "series": [],
                }
            )
            .style("height: 450px; margin: 20px 0;")
            .classes("border rounded-lg shadow-sm")
        )

    def analyze_cytoband_cnv(self, cnv_data: dict, chromosome: str) -> pd.DataFrame:
        """
        Analyze CNV values within each cytoband to detect duplications and deletions.
        Uses dynamic thresholds based on data variation for more robust detection.

        Args:
            cnv_data (dict): Dictionary containing CNV values
            chromosome (str): Chromosome to analyze

        Returns:
            pd.DataFrame: DataFrame containing merged cytoband CNV analysis results
        """
        logger.debug(f"\n{'='*50}")
        logger.debug(f"Starting CNV analysis for {chromosome}")
        logger.debug(f"CNV data keys: {list(cnv_data.keys())}")
        logger.debug(f"Bin width: {self.cnv_dict['bin_width']}")

        # Check if bin width is small enough for accurate CNV calling
        if self.cnv_dict["bin_width"] > 10_000_000:
            logger.debug("Resolution insufficient for CNV calling")
            return pd.DataFrame()

        bin_width = self.cnv_dict["bin_width"]
        chromosome_cytobands = self.cytobands_bed[
            self.cytobands_bed["chrom"] == chromosome
        ].copy()
        logger.debug(
            f"Number of cytobands for {chromosome}: {len(chromosome_cytobands)}"
        )

        # Initialize results storage
        merged_cytobands = []
        whole_chr_event = False
        whole_chr_state = "NORMAL"

        # First, analyze the whole chromosome for potential aneuploidy
        if chromosome in cnv_data:
            logger.debug(
                f"\nAnalyzing chromosome {chromosome} for whole chromosome events:"
            )

            # Calculate mean CNV for the entire chromosome, excluding centromeric regions
            centromere = self.centromere_bed[self.centromere_bed["chrom"] == chromosome]
            if not centromere.empty:
                cent_start_bin = int(centromere["start_pos"].iloc[0] / bin_width)
                cent_end_bin = int(centromere["end_pos"].iloc[0] / bin_width)
                chr_cnv = np.concatenate(
                    [
                        cnv_data[chromosome][:cent_start_bin],
                        cnv_data[chromosome][cent_end_bin:],
                    ]
                )
                logger.debug(
                    f"Excluded centromere region: {cent_start_bin}-{cent_end_bin}"
                )
            else:
                chr_cnv = cnv_data[chromosome]

            # Calculate chromosome-wide statistics
            chr_mean = np.mean(chr_cnv)
            chr_std = np.std(chr_cnv)
            logger.debug(f"Chromosome-wide mean: {chr_mean:.3f}, std: {chr_std:.3f}")

            # Calculate SD of chromosome means for whole chromosome event detection
            chromosome_means = []
            for chrom in cnv_data:
                if chrom.startswith("chr") and chrom[3:].isdigit():  # Only autosomes
                    chrom_data = cnv_data[chrom]
                    # Exclude centromere regions if present
                    cent = self.centromere_bed[self.centromere_bed["chrom"] == chrom]
                    if not cent.empty:
                        cent_start = int(cent["start_pos"].iloc[0] / bin_width)
                        cent_end = int(cent["end_pos"].iloc[0] / bin_width)
                        chrom_data = np.concatenate(
                            [chrom_data[:cent_start], chrom_data[cent_end:]]
                        )
                    chromosome_means.append(np.mean(chrom_data))

            means_std = np.std(chromosome_means)
            means_mean = np.mean(chromosome_means)
            logger.debug(
                f"Mean of chromosome means: {means_mean:.3f}, std of means: {means_std:.3f}"
            )

            # Base thresholds on standard deviations from the mean
            if chromosome.startswith("chr") and chromosome[3:].isdigit():  # Autosomes
                # For whole chromosome events, use SD of means with 70% confidence
                gain_threshold = means_mean + (1.0 * means_std)  # ~70% confidence
                loss_threshold = means_mean - (1.0 * means_std)
                # For focal events, use chromosome-specific SD
                cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
                cytoband_loss_threshold = chr_mean - (1.0 * chr_std)
            elif chromosome == "chrX":
                if self.XYestimate == "XY":  # Male
                    gain_threshold = means_mean + (1.0 * means_std)
                    loss_threshold = means_mean - (1.0 * means_std)
                    cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
                    cytoband_loss_threshold = chr_mean - (1.0 * chr_std)
                else:  # Female
                    gain_threshold = means_mean + (1.0 * means_std)
                    loss_threshold = means_mean - (1.0 * means_std)
                    cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
                    cytoband_loss_threshold = chr_mean - (1.0 * chr_std)
            elif chromosome == "chrY":
                if self.XYestimate == "XY":  # Male
                    gain_threshold = means_mean + (1.0 * means_std)
                    loss_threshold = means_mean - (1.0 * means_std)
                    cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
                    cytoband_loss_threshold = chr_mean - (1.0 * chr_std)
                else:  # Female
                    # For Y in females, use slightly more extreme thresholds
                    gain_threshold = means_mean + (1.2 * means_std)
                    loss_threshold = means_mean - (1.2 * means_std)
                    cytoband_gain_threshold = chr_mean + (1.2 * chr_std)
                    cytoband_loss_threshold = chr_mean - (1.2 * chr_std)

            logger.debug(
                f"Thresholds - Whole chr gain: {gain_threshold:.3f}, loss: {loss_threshold:.3f}"
            )
            logger.debug(
                f"Thresholds - Cytoband gain: {cytoband_gain_threshold:.3f}, loss: {cytoband_loss_threshold:.3f}"
            )

            # Calculate proportion of bins supporting gain/loss using thresholds
            bins_above_gain = np.sum(chr_cnv > gain_threshold) / len(chr_cnv)
            bins_below_loss = np.sum(chr_cnv < loss_threshold) / len(chr_cnv)

            logger.debug(
                f"Proportion of bins - Above gain: {bins_above_gain:.3f}, Below loss: {bins_below_loss:.3f}"
            )

            # Detect whole chromosome events with proportion threshold
            min_proportion = 0.7  # Require at least 50% of bins to support the event
            if bins_above_gain > min_proportion:
                whole_chr_event = True
                whole_chr_state = "GAIN"
                logger.debug(f"WHOLE CHROMOSOME EVENT DETECTED: {chromosome} GAIN")
            elif bins_below_loss > min_proportion:
                whole_chr_event = True
                whole_chr_state = "LOSS"
                logger.debug(f"WHOLE CHROMOSOME EVENT DETECTED: {chromosome} LOSS")

            # If whole chromosome event detected, add it to results
            if whole_chr_event:
                genes_in_chr = self.gene_bed[self.gene_bed["chrom"] == chromosome][
                    "gene"
                ].tolist()

                merged_cytobands.append(
                    {
                        "chrom": chromosome,
                        "start_pos": chromosome_cytobands["start_pos"].min(),
                        "end_pos": chromosome_cytobands["end_pos"].max(),
                        "name": f"{chromosome} WHOLE CHROMOSOME {whole_chr_state}",
                        "mean_cnv": chr_mean,
                        "cnv_state": whole_chr_state,
                        "length": chromosome_cytobands["end_pos"].max()
                        - chromosome_cytobands["start_pos"].min(),
                        "genes": genes_in_chr,
                    }
                )

            # Now analyze individual cytobands regardless of whole chromosome event
            current_group = None

            for _, cytoband in chromosome_cytobands.iterrows():
                start_bin = int(cytoband["start_pos"] / bin_width)
                end_bin = int(cytoband["end_pos"] / bin_width)

                if start_bin < len(cnv_data[chromosome]):
                    region_cnv = cnv_data[chromosome][start_bin : end_bin + 1]
                    mean_cnv = np.mean(region_cnv) if len(region_cnv) > 0 else 0

                    # Determine cytoband state relative to whole chromosome state
                    if whole_chr_event:
                        # For whole chromosome events, use more stringent thresholds relative to chromosome mean
                        if mean_cnv > chr_mean + (2.0 * chr_std):  # More stringent threshold for additional gains
                            state = "HIGH_GAIN"  # Report significant additional gains
                        elif mean_cnv < chr_mean - (2.0 * chr_std):  # More stringent threshold for additional losses
                            state = "DEEP_LOSS"  # Report significant additional losses
                        else:
                            state = whole_chr_state  # Maintain the whole chromosome state
                    else:
                        # Normal threshold-based state determination
                        if mean_cnv > cytoband_gain_threshold:
                            state = "GAIN"
                        elif mean_cnv < cytoband_loss_threshold:
                            state = "LOSS"
                        else:
                            state = "NORMAL"
                else:
                    mean_cnv = 0
                    state = "NO_DATA"

                # Group cytobands with same state
                if current_group is None:
                    current_group = {
                        "chrom": cytoband["chrom"],
                        "start_pos": cytoband["start_pos"],
                        "end_pos": cytoband["end_pos"],
                        "name": cytoband["name"],
                        "mean_cnv": [mean_cnv],
                        "cnv_state": state,
                        "bands": [cytoband["name"]],
                        "length": cytoband["end_pos"] - cytoband["start_pos"],
                        "genes": [],
                    }
                elif state == current_group["cnv_state"]:
                    current_group["end_pos"] = cytoband["end_pos"]
                    current_group["mean_cnv"].append(mean_cnv)
                    current_group["bands"].append(cytoband["name"])
                else:
                    # Process current group
                    if current_group["cnv_state"] in [
                        "GAIN",
                        "LOSS",
                        "HIGH_GAIN",
                        "DEEP_LOSS",
                    ]:
                        genes_in_region = self.gene_bed[
                            (self.gene_bed["chrom"] == current_group["chrom"])
                            & (self.gene_bed["start_pos"] <= current_group["end_pos"])
                            & (self.gene_bed["end_pos"] >= current_group["start_pos"])
                        ]["gene"].tolist()
                        current_group["genes"] = genes_in_region

                    current_group["name"] = (
                        f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}"
                    )
                    current_group["mean_cnv"] = np.mean(current_group["mean_cnv"])
                    current_group["length"] = (
                        current_group["end_pos"] - current_group["start_pos"]
                    )

                    # Add all significant changes, including those on chromosomes with whole chromosome events
                    if current_group["cnv_state"] in ["GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"]:
                        merged_cytobands.append(current_group)

                    # Start new group
                    current_group = {
                        "chrom": cytoband["chrom"],
                        "start_pos": cytoband["start_pos"],
                        "end_pos": cytoband["end_pos"],
                        "name": cytoband["name"],
                        "mean_cnv": [mean_cnv],
                        "cnv_state": state,
                        "bands": [cytoband["name"]],
                        "length": cytoband["end_pos"] - cytoband["start_pos"],
                        "genes": [],
                    }

            # Process the last group
            if current_group is not None:
                if current_group["cnv_state"] in [
                    "GAIN",
                    "LOSS",
                    "HIGH_GAIN",
                    "DEEP_LOSS",
                ]:
                    genes_in_region = self.gene_bed[
                        (self.gene_bed["chrom"] == current_group["chrom"])
                        & (self.gene_bed["start_pos"] <= current_group["end_pos"])
                        & (self.gene_bed["end_pos"] >= current_group["start_pos"])
                    ]["gene"].tolist()
                    current_group["genes"] = genes_in_region

                current_group["name"] = (
                    f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}"
                )
                current_group["mean_cnv"] = np.mean(current_group["mean_cnv"])
                current_group["length"] = (
                    current_group["end_pos"] - current_group["start_pos"]
                )

                # Add all significant changes, including those on chromosomes with whole chromosome events
                if current_group["cnv_state"] in ["GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"]:
                    merged_cytobands.append(current_group)

        # Convert results to DataFrame and sort
        merged_df = pd.DataFrame(merged_cytobands)
        if not merged_df.empty:
            merged_df = merged_df.sort_values("start_pos")

        return merged_df

    def get_cytoband_cnv_summary(self, chromosome: str) -> str:
        """
        Generate a summary of CNV states for cytobands in a chromosome.

        Args:
            chromosome (str): Chromosome to summarize

        Returns:
            str: Summary string of gains and losses
        """
        if not hasattr(self, "result3") or not self.result3.cnv:
            return "No CNV data available"

        cytoband_analysis = self.analyze_cytoband_cnv(self.result3.cnv, chromosome)

        # Filter for gains and losses
        gains = cytoband_analysis[cytoband_analysis["cnv_state"] == "GAIN"]
        losses = cytoband_analysis[cytoband_analysis["cnv_state"] == "LOSS"]

        summary = []
        if not gains.empty:
            gain_bands = [
                f"{row['name']} ({row['mean_cnv']:.2f})" for _, row in gains.iterrows()
            ]
            summary.append(f"Gains: {', '.join(gain_bands)}")

        if not losses.empty:
            loss_bands = [
                f"{row['name']} ({row['mean_cnv']:.2f})" for _, row in losses.iterrows()
            ]
            summary.append(f"Losses: {', '.join(loss_bands)}")

        return "\n".join(summary) if summary else "No significant CNV changes detected"

    def update_target_table(self) -> None:
        """Update the target panel information table with current bed file data."""
        latest_bed = self.get_latest_bed_file()

        try:
            if not latest_bed:
                if self.target_table:
                    message_df = pd.DataFrame(
                        [
                            {
                                "chrom": "",
                                "gene": "No target data available",
                                "start_pos": 0,
                                "end_pos": 0,
                                "size": 0,
                                "status": "Inactive",
                                "source": "None",
                            }
                        ]
                    )
                    self.target_table.rows = message_df.to_dict("records")
                    ui.update(self.target_table)
                return

            # Read the BED file
            bed_data = pd.read_csv(
                latest_bed,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "name", "score", "strand"],
            )

            # Prepare rows data
            table_rows = []
            for _, bed_row in bed_data.iterrows():
                # Find overlapping genes from gene_bed
                overlapping_genes = self.gene_bed[
                    (self.gene_bed["chrom"] == bed_row["chrom"])
                    & (
                        (
                            (self.gene_bed["start_pos"] >= bed_row["start"])
                            & (self.gene_bed["start_pos"] <= bed_row["end"])
                        )
                        | (
                            (self.gene_bed["end_pos"] >= bed_row["start"])
                            & (self.gene_bed["end_pos"] <= bed_row["end"])
                        )
                        | (
                            (self.gene_bed["start_pos"] <= bed_row["start"])
                            & (self.gene_bed["end_pos"] >= bed_row["end"])
                        )
                    )
                ]

                # Get list of overlapping gene names
                gene_names = (
                    overlapping_genes["gene"].tolist()
                    if not overlapping_genes.empty
                    else ["Unknown"]
                )
                gene_label = ", ".join(gene_names)

                table_rows.append(
                    {
                        "chrom": str(bed_row["chrom"]),
                        "gene": gene_label,
                        "start_pos": int(bed_row["start"]),
                        "end_pos": int(bed_row["end"]),
                        "size": int(bed_row["end"] - bed_row["start"]),
                        "status": "Active",
                        "source": str(bed_row["name"]),
                    }
                )

            # Convert to DataFrame for sorting
            target_data = pd.DataFrame(table_rows)
            if not target_data.empty:
                # Sort by chromosome and start position
                target_data = target_data.sort_values(["chrom", "start_pos"])

            try:
                if self.target_table:
                    self.target_table.rows = target_data.to_dict("records")
                    ui.update(self.target_table)
            except Exception as e:
                logger.error(f"Table update error: {e}")

        except Exception as e:
            logger.error(f"Error processing bed file {latest_bed}: {e}")
            if self.target_table:
                message_df = pd.DataFrame(
                    [
                        {
                            "chrom": "",
                            "gene": f"Error loading target data: {str(e)}",
                            "start_pos": 0,
                            "end_pos": 0,
                            "size": 0,
                            "status": "Error",
                            "source": "Error",
                        }
                    ]
                )
                self.target_table.rows = message_df.to_dict("records")
                ui.update(self.target_table)


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
