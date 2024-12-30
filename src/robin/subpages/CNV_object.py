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
        (cp  * bin_width - (bin_width), cp * bin_width + (bin_width))
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

    def __init__(self, *args, target_panel: Optional[str] = None, reference_file: Optional[str] = None, bed_file: Optional[str] = None, readfish_toml: Optional[Path] = None, **kwargs) -> None:
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
        super().__init__(*args, **kwargs)
        self.NewBed = BedTree(preserve_original_tree=True, reference_file=f"{self.reference_file}.fai", readfish_toml=self.readfish_toml)
        self.NewBed.load_from_file(self.bed_file)
        self.CNVchangedetector = CNVChangeDetectorTracker(base_proportion=0.02)


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
                    penalty_value = 5
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
                                    self.local_data_array = np.append(self.local_data_array, item)
                                    self.DATA_ARRAY = np.append(self.DATA_ARRAY, item)
                        self.map_tracker[key] = 0
                        cnvupdate = True
                        breakpoints = self.DATA_ARRAY[self.DATA_ARRAY["name"] == key]
                        local_breakpoints = self.local_data_array[self.local_data_array["name"] == key]
                        if len(breakpoints) > 0:
                            try:
                                self.CNVchangedetector.add_breakpoints(
                                    key, breakpoints, local_breakpoints, approx_chrom_length, r_bin
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
                if len(tempbedcontent)>0:
                    bedcontent += tempbedcontent
                    bedcontent += "\n"

                tempbedcontent2 = self.CNVchangedetector.get_bed_targets_breakpoints(chrom)
                if len(tempbedcontent2)>0:
                    bedcontent2 += tempbedcontent2
                    bedcontent2 += "\n"

            #if len(bedcontent)>0:
            #    print ("bedcontent")
            #    print(f"{bedcontent}")
            #if len(bedcontent2)>0:
            #    print ("bedcontent2")
            #    print(f"{bedcontent2}")
            if len(bedcontent2)>0:
                #print(bedcontent2)
                self.NewBed.load_from_string(bedcontent2, merge=False, write_files=True, output_location=os.path.join(
                              self.check_and_create_folder(self.output, self.sampleID))
                                             )

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
        """
        Update CNV plots with new data and annotations.

        Args:
            gene_target (Optional[str]): Target gene for updating plots.
        """
        if not gene_target:
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

    def setup_ui(self) -> None:
        """
        Set up the user interface for CNV analysis.
        """

        self.display_row = ui.row().style("width: 100")
        if self.summary:
            with self.summary:
                ui.label("No CNV data available.")
        with self.display_row:
            ui.label("Copy Number Variation").classes('text-sky-600 dark:text-white').style(
                "font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
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
                self.cnv_dict, "bin_width", backward=lambda n: f"Bin Width: {n}"
            )
            ui.label().bind_text_from(
                self.cnv_dict, "variance", backward=lambda n: f"Variance: {round(n, 3)}"
            )

        self.scatter_echart = self.generate_chart(title="CNV Scatter Plot")
        self.difference_scatter_echart = self.generate_chart(
            title="Difference Plot", initmin=-2, initmax=2
        )

        with ui.expansion("See Reference DataSet", icon="loupe").classes("w-full"):
            self.reference_scatter_echart = self.generate_chart(
                title="Reference CNV Scatter Plot"
            )



        with ui.card().classes("w-full"):

            ui.label("New Target Information - intentionally blank.")

            # self.orig_tree = ui.tree(
            #     [self.NewBed.tree_data]
            #     , label_key="id",
            # )
            # self.orig_tree.add_slot(
            #     "default-body",
            #     """
            #     <div v-if="props.node.description">
            #             <span class="text-weight-bold">{{ props.node.description }}</span>
            #         </div>
            #     <div v-if="props.node.range_sum">
            #             <span class="text-weight-bold">{{props.node.range_sum}} bases</span>
            #     </div>
            #     <div v-if="props.node.count">
            #             <span class="text-weight-bold">{{ props.node.count }} targets</span>
            #     </div>
            #     <div v-if="props.node.proportion">
            #             <span class="text-weight-bold">{{ props.node.proportion }} proportion</span>
            #     </div>
            #     <div v-if="props.node.chromosome_length">
            #             <span class="text-weight-bold">{{ props.node.chromosome_length }} chromosome length</span>
            #     </div>
            #
            # """,
            # )


        with ui.card().classes("w-full"):
            ui.label("Proportion Over Time Information")
            self.create_proportion_time_chart2("Proportions over time - Genome Wide.")
            self.create_proportion_time_chart("Proportions over time.")


        if self.browse:
            ui.timer(0.1, lambda: self.show_previous_data(), once=True)
        else:
            ui.timer(15, lambda: self.show_previous_data())

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
                self.proportion_time_chart.options["title"]["text"] = "No Data Available"
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
                    data_list = [[key, value] for key, value in data.items() if pd.notnull(value)]
                    if data_list:  # Only add series if it has valid data points
                        self.proportion_time_chart.options["series"].append({
                            "animation": False,
                            "type": "line",
                            "smooth": True,
                            "name": series,
                            "emphasis": {
                                "focus": "series",
                                "itemStyle": {
                                    "borderWidth": 2
                                }
                            },
                            "endLabel": {
                                "show": True,
                                "formatter": "{a}: {c}%",
                                "distance": 10,
                                "fontSize": 12
                            },
                            "lineStyle": {
                                "width": 2,
                                "color": colors[idx % len(colors)]
                            },
                            "itemStyle": {
                                "color": colors[idx % len(colors)]
                            },
                            "data": data_list
                        })

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

    def generate_chart(
        self,
        title: Optional[str] = None,
        initmax: Optional[int] = None,
        initmin: int = 0,
        type: str = "value",
    ) -> ui.echart:
        """
        Generate an ECharts object for displaying CNV scatter plots.

        Creates a chart with consistent styling following Apple HIG guidelines and matching
        other analysis modules. The chart includes interactive features like zooming and
        saving, while maintaining visual consistency across the application.

        Parameters
        ----------
        title : str, optional
            Title of the chart. Will be displayed prominently at the top.
        initmax : int, optional
            Initial maximum value for the y-axis range.
        initmin : int, default=0
            Initial minimum value for the y-axis range.
        type : str, default="value"
            Type of x-axis to use.

        Returns
        -------
        ui.echart
            Configured ECharts object ready for CNV data visualization.

        Notes
        -----
        The chart follows Apple HIG principles:
        - Uses consistent colors from the iOS palette
        - Maintains clear visual hierarchy
        - Provides appropriate contrast for accessibility
        - Includes interactive elements for data exploration
        """
        return (
            ui.echart(
                {
                    "backgroundColor": "transparent",
                    "textStyle": {
                        "fontFamily": "SF Pro Text, -apple-system, BlinkMacSystemFont, Helvetica, Arial, sans-serif",
                        "fontSize": 12
                    },
                    "animation": False,
                    "grid": {
                        "top": 80,
                        "bottom": 90,
                        "left": 80,
                        "right": 80,
                        "containLabel": True
                    },
                    "title": {
                        "text": f"{title}",
                        "left": "center",
                        "top": 20,
                        "textStyle": {
                            "fontSize": 16,
                            "fontWeight": "500",
                            "color": "#1D1D1F"
                        }
                    },
                    "toolbox": {
                        "show": True,
                        "right": 20,
                        "feature": {
                            "saveAsImage": {
                                "title": "Save Chart"
                            }
                        }
                    },
                    "tooltip": {
                        "trigger": "axis",
                        "backgroundColor": "rgba(255, 255, 255, 0.9)",
                        "borderColor": "#E5E5EA",
                        "textStyle": {
                            "color": "#1D1D1F"
                        }
                    },
                    "xAxis": {
                        "type": f"{type}",
                        "max": "dataMax",
                        "splitLine": {
                            "show": True,
                            "lineStyle": {
                                "type": "dashed",
                                "color": "#E5E5EA"
                            }
                        },
                        "axisLabel": {
                            "color": "#86868B"
                        }
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
                                "padding": [0, 30, 0, 0]
                            },
                            "axisLabel": {
                                "color": "#86868B"
                            },
                            "splitLine": {
                                "show": True,
                                "lineStyle": {
                                    "type": "dashed",
                                    "color": "#E5E5EA"
                                }
                            }
                        },
                        {
                            "type": "value",
                            "name": "Candidate Break Points",
                            "position": "right",
                            "nameTextStyle": {
                                "color": "#007AFF",
                                "fontSize": 12,
                                "padding": [0, 0, 0, 30]
                            },
                            "axisLine": {
                                "lineStyle": {
                                    "color": "#007AFF"
                                }
                            },
                            "axisLabel": {
                                "color": "#007AFF"
                            },
                            "splitLine": {
                                "show": False
                            }
                        }
                    ],
                    "dataZoom": [
                        {
                            "type": "slider",
                            "yAxisIndex": "none",
                            "xAxisIndex": "0",
                            "filterMode": "none",
                            "height": 20,
                            "bottom": 50,
                            "borderColor": "#E5E5EA",
                            "backgroundColor": "#F5F5F7",
                            "fillerColor": "rgba(0, 122, 255, 0.2)",
                            "handleStyle": {
                                "color": "#007AFF"
                            }
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
                            "handleStyle": {
                                "color": "#007AFF"
                            }
                        }
                    ],
                    "series": [
                        {
                            "type": "scatter",
                            "symbolSize": 4,
                            "itemStyle": {
                                "color": "#007AFF",
                                "opacity": 0.7
                            },
                            "emphasis": {
                                "itemStyle": {
                                    "color": "#007AFF",
                                    "opacity": 1,
                                    "borderColor": "#FFFFFF",
                                    "borderWidth": 2
                                }
                            },
                            "data": []
                        }
                    ]
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

            # If all chromosomes are selected, prepare the plot to display all chromosomes
            if self.chrom_filter == "All":
                counter = 0  # Reset counter for chromosome loop

                # Update the plot title to reflect that all chromosomes are being shown
                plot_to_update.options["title"]["text"] = f"{title} - All Chromosomes"
                plot_to_update.options["series"] = (
                    []
                )  # Clear any existing series in the plot

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
                plot_to_update.options["series"] = (
                    []
                )  # Clear existing series in the plot

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
                plot_to_update.options["title"][
                    "text"
                ] = f"Copy Number Variation - {contig}"

                plot_to_update.options["xAxis"]["max"] = max
                plot_to_update.options["xAxis"]["min"] = min
                #print(plot_to_update.options["dataZoom"][1])
                #plot_to_update.options["dataZoom"][1]["startValue"] = 0
                plot_to_update.options["dataZoom"][1]["endValue"] = ymax
                plot_to_update.options["series"].append(
                    {
                        "type": "scatter",
                        "name": contig,
                        "data": data,
                        "symbolSize": 3,  # Smaller point size for this view
                        "markArea": {
                            "itemStyle": {"color": "rgba(255, 173, 177, 0.4)"},
                            "data": [],
                        },
                    }
                )

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
                            },
                            {
                                "xAxis": row["end_pos"],
                            },
                        ]
                    )

                plot_to_update.options["series"].append(
                    {
                        "type": "scatter",
                        "name": contig + "_highlight",
                        "data": [],  # Empty data because this series is just for highlighting
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {
                                "color": "rgba(135, 206, 250, 0.4)"
                            },  # Light blue color
                            "data": [],
                        },
                    }
                )

                # Highlight the centromeres with another shaded area:
                for _, row in self.centromere_bed[
                    self.centromere_bed["chrom"] == contig
                ].iterrows():
                    plot_to_update.options["series"][1]["markArea"]["data"].append(
                        [
                            {
                                "name": row["name"],
                                "xAxis": row["start_pos"],
                            },
                            {
                                "xAxis": row["end_pos"],
                            },
                        ]
                    )

                if contig in self.CNVResults.keys():
                    coordinates = self.CNVResults[contig]["positions"]
                    threshold = self.CNVResults[contig]["threshold"]
                    plot_to_update.options["series"].append(
                        {
                            "type": "line",
                            "name": "Line Plot",
                            "data": coordinates,  # Your line plot data (e.g., [(x1, y1), (x2, y2), ...])
                            "yAxisIndex": 1,  # Specify to use the secondary y-axis (index 1)
                            "lineStyle": {
                                "color": "blue"
                            },  # Color the line to match the y-axis
                            "symbol": "none",  # Optionally remove symbols from the line plot
                            "markLine": {
                                "symbol": "none",
                                "data": [
                                    {
                                        "lineStyle": {"width": 2},
                                        "label": {"formatter": "Selection Threshold"},
                                        "name": "Threshold",
                                        "color": "green",
                                        "yAxis": (threshold),
                                    },
                                ],
                            },
                        }
                    )
                    starts = self.CNVResults[contig]["start_positions"]
                    ends = self.CNVResults[contig]["end_positions"]
                    # (starts, ends) = detector.get_breakpoints()
                    # (starts, ends) = [],[]
                    for start in starts:
                        plot_to_update.options["series"][2]["markLine"]["data"].append(
                            {
                                "lineStyle": {"width": 2},
                                "label": {"formatter": "s"},
                                "name": "Target",
                                "color": "green",
                                "xAxis": (start),
                            },
                        )

                    for end in ends:
                        plot_to_update.options["series"][2]["markLine"]["data"].append(
                            {
                                "lineStyle": {"width": 2},
                                "label": {"formatter": "e"},
                                "name": "Target",
                                "color": "red",
                                "xAxis": (end),
                            },
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

    def create_summary_card(self, xy_estimate: str, bin_width: int, variance: float) -> None:
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

        Notes
        -----
        The summary card follows Apple HIG guidelines for:
        - Typography and spacing
        - Color usage and contrast
        - Information hierarchy
        - Visual feedback
        """
        with ui.card().classes("w-full p-4 bg-white dark:bg-gray-800"):
            with ui.row().classes("items-center gap-4"):
                # Genetic Sex Icon and Information
                if xy_estimate != "Unknown":
                    if xy_estimate == "XY":
                        ui.icon("man").classes("text-4xl text-blue-500")
                    else:
                        ui.icon("woman").classes("text-4xl text-pink-500")
                    ui.label(f"Estimated Genetic Sex: {xy_estimate}").classes(
                        "text-lg font-medium text-gray-900 dark:text-white"
                    )
                
                # Analysis Metrics
                with ui.column().classes("gap-2"):
                    ui.label(f"Bin Width: {bin_width:,}").classes(
                        "text-sm text-gray-600 dark:text-gray-300"
                    )
                    ui.label(f"Variance: {variance:.3f}").classes(
                        "text-sm text-gray-600 dark:text-gray-300"
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

            self.NewBed.load_from_string(bedcontent, merge=False)

            if self.check_file_time(os.path.join(output, "bedranges.csv")):
                self.proportions_df_store = pd.read_csv(
                    os.path.join(os.path.join(output, "bedranges.csv")),
                    index_col=0,
                )

                df_filtered = self.proportions_df_store[self.proportions_df_store['chromosome'] != 'Genome-wide']
                pivot_df = df_filtered.pivot(index='timestamp', columns='chromosome', values='proportion')

                df_filtered2 = self.proportions_df_store[self.proportions_df_store['chromosome'] == 'Genome-wide']
                pivot_df2 = df_filtered2.pivot(index='timestamp', columns='chromosome', values='proportion')

                self.update_proportion_time_chart(pivot_df)
                self.update_proportion_time_chart2(pivot_df2)

            if self.summary:
                with self.summary:
                    self.summary.clear()
                    with open(os.path.join(output, "XYestimate.pkl"), "rb") as file:
                        xy_estimate = pickle.load(file)
                        self.create_summary_card(
                            xy_estimate=xy_estimate,
                            bin_width=self.cnv_dict["bin_width"],
                            variance=self.cnv_dict["variance"]
                        )


    def create_time_chart(self, title: str) -> ui.echart:
        """
        Create a time series chart for CNV data visualization.

        Creates a chart specifically designed for time-based CNV data display,
        following Apple HIG guidelines and maintaining consistency with other
        analysis modules.

        Parameters
        ----------
        title : str
            Title of the chart to be displayed.

        Returns
        -------
        ui.echart
            Configured ECharts object for time series visualization.

        Notes
        -----
        The chart includes:
        - iOS-style color palette for multiple series
        - Interactive legend with emphasis states
        - Smooth animations for data updates
        - Consistent typography and spacing
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

        return ui.echart({
            "backgroundColor": "transparent",
            "textStyle": {
                "fontFamily": "SF Pro Text, -apple-system, BlinkMacSystemFont, Helvetica, Arial, sans-serif",
                "fontSize": 12
            },
            "animation": False,
            "grid": {
                "top": 80,
                "bottom": 90,
                "left": 80,
                "right": 80,
                "containLabel": True
            },
            "title": {
                "text": title,
                "left": "center",
                "top": 20,
                "textStyle": {
                    "fontSize": 16,
                    "fontWeight": "500",
                    "color": "#1D1D1F"
                }
            },
            "tooltip": {
                "trigger": "axis",
                "backgroundColor": "rgba(255, 255, 255, 0.9)",
                "borderColor": "#E5E5EA",
                "textStyle": {
                    "color": "#1D1D1F"
                }
            },
            "legend": {
                "type": "scroll",
                "top": 50,
                "textStyle": {
                    "color": "#86868B"
                }
            },
            "toolbox": {
                "show": True,
                "right": 20,
                "feature": {
                    "saveAsImage": {
                        "title": "Save Chart"
                    }
                }
            },
            "xAxis": {
                "type": "time",
                "splitLine": {
                    "show": True,
                    "lineStyle": {
                        "type": "dashed",
                        "color": "#E5E5EA"
                    }
                },
                "axisLabel": {
                    "color": "#86868B"
                }
            },
            "yAxis": {
                "type": "value",
                "name": "Proportion",
                "nameTextStyle": {
                    "color": "#86868B",
                    "fontSize": 12,
                    "padding": [0, 30, 0, 0]
                },
                "axisLabel": {
                    "color": "#86868B",
                    "formatter": "{value}%"
                },
                "splitLine": {
                    "show": True,
                    "lineStyle": {
                        "type": "dashed",
                        "color": "#E5E5EA"
                    }
                }
            },
            "color": colors,
            "series": []
        }).style("height: 450px; margin: 20px 0;").classes("border rounded-lg shadow-sm")


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
