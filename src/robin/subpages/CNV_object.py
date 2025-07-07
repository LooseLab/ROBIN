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
from robin.subpages.base_analysis import BaseVis, BaseAnalysis
from robin.utilities.break_point_detector import CNVChangeDetectorTracker
from robin.utilities.bed_file import MasterBedTree
import natsort
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
from typing import Optional, Tuple, List, BinaryIO
import re

from collections import Counter
from scipy.ndimage import uniform_filter1d

import math
import time

from robin.core.state import state, ProcessState
import shutil

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


class CNVVis(BaseVis):
    """
    A class for visualizing CNV data.

    """

    def __init__(
        self,
        *args,
        target_panel: Optional[str] = None,
        reference_file: Optional[str] = None,
        bed_file: Optional[str] = None,
        readfish_toml: Optional[Path] = None,
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
        self.cnv_dict = {"bin_width": 0, "variance": 0}
        self.update_cnv_dict = {}
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
        # self.DATA_ARRAY = None
        self.data_array_size = 0
        self.data_array_path = None
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

        self.CNVchangedetector = CNVChangeDetectorTracker(base_proportion=0.02)
        # Add target_table as instance variable
        self.target_table = None
        self.target_table_placeholder = None
        self.last_bed_check = 0  # Track when we last checked for new BED files

        # Add y_axis_log attribute
        self.y_axis_log = False  # Default to linear scale
        # Add show_breakpoints attribute
        self.show_breakpoints = True  # Default to showing breakpoints

        # Add DATA_ARRAY for breakpoint visualization
        self.DATA_ARRAY = None
        self.data_array_path = None

    async def setup_ui(self) -> None:
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
                    self.update_plots(preserve_zoom=True)
                    if e.value == "All"
                    else self.update_plots(gene_target=e.value, preserve_zoom=True)
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
            # Add a toggle switch for color mode
            with ui.row().classes("items-center gap-2"):
                ui.label("Color by:").classes("text-sm")

                # Define callback function before creating the toggle
                def toggle_color_mode(e):
                    self.color_mode = e.value
                    self.update_plots(preserve_zoom=True)

                # Create toggle with options as a dict and pass on_change in constructor
                ui.toggle(
                    options={"chromosome": "Chromosome", "value": "Up/Down"},
                    value="chromosome",
                    on_change=toggle_color_mode,
                ).classes("mt-1")

                # Add y-axis scale toggle
                ui.label("Y-axis scale:").classes("text-sm ml-4")

                def toggle_y_axis_scale(e):
                    self.y_axis_log = e.value == "log"
                    # Update only the first plot (scatter_echart)
                    self.scatter_echart.options["yAxis"][0]["logBase"] = (
                        10 if self.y_axis_log else None
                    )
                    self.scatter_echart.options["yAxis"][0]["type"] = (
                        "log" if self.y_axis_log else "value"
                    )
                    ui.update(self.scatter_echart)

                # Create toggle for y-axis scale
                ui.toggle(
                    options={"linear": "Linear", "log": "Log"},
                    value="linear",  # Default to linear scale
                    on_change=toggle_y_axis_scale,
                ).classes("mt-1")

                # Add breakpoint visualization toggle
                ui.label("Breakpoints:").classes("text-sm ml-4")

                def toggle_breakpoint_visualization(e):
                    self.show_breakpoints = e.value == "show"
                    self.update_plots(preserve_zoom=True)

                ui.toggle(
                    options={"hide": "Hide", "show": "Show"},
                    value="show",  # Default to showing breakpoints
                    on_change=toggle_breakpoint_visualization,
                ).classes("mt-1")

        self.scatter_echart = self.generate_chart(
            title="CNV Scatter Plot", initmin=0, initmax=8
        )

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

        # await ui.context.client.connected()

        if self.browse:
            ui.timer(0.1, lambda: self.show_previous_data(), once=True)
        else:
            ui.timer(15, lambda: self.show_previous_data())

    def update_plots(
        self, gene_target: Optional[str] = None, preserve_zoom: bool = False
    ) -> None:
        """
        Update CNV plots with new data and annotations.

        Parameters
        ----------
        gene_target : Optional[str], optional
            Target gene to focus the plot on. If None, uses current chromosome selection.
        preserve_zoom : bool, default False
            Whether to preserve the current zoom settings. If True, user zoom settings
            are maintained. If False, zoom is reset to show the full data range.
            Set to True for UI control changes (color mode, breakpoints, etc.) and
            False for chromosome/gene selection changes.
        """
        if not gene_target:
            # Reset zoom when switching to "All" chromosomes view
            if self.chrom_select.value == "All":
                self.scatter_echart.options["dataZoom"][0].update(
                    {"startValue": None, "endValue": None}
                )
                self.difference_scatter_echart.options["dataZoom"][0].update(
                    {"startValue": None, "endValue": None}
                )
                # Disable breakpoint visualization in 'All' chromosomes view
                self.show_breakpoints = False
                # If the toggle UI exists, update it to reflect the state
                try:
                    for child in self.chrom_select.parent.children:
                        if (
                            hasattr(child, "label")
                            and getattr(child, "label", None) == "Breakpoints:"
                        ):
                            # The next child should be the toggle
                            idx = self.chrom_select.parent.children.index(child)
                            toggle = self.chrom_select.parent.children[idx + 1]
                            if hasattr(toggle, "value"):
                                toggle.value = "hide"
                                ui.update(toggle)
                except Exception:
                    pass
                # When switching to "All" view, don't preserve zoom
                preserve_zoom = False
            else:
                # In single chromosome view, allow breakpoints to be shown (default to previous state or 'show')
                if not hasattr(self, "_breakpoint_toggle_initialized"):
                    self.show_breakpoints = True
                    self._breakpoint_toggle_initialized = True
            self._update_cnv_plot(
                plot_to_update=self.scatter_echart,
                result=self.result,
                title="CNV",
                preserve_zoom=preserve_zoom,
            )
            self._update_cnv_plot(
                plot_to_update=self.reference_scatter_echart,
                result=self.result2,
                title="Reference CNV",
                preserve_zoom=preserve_zoom,
            )
            self._update_cnv_plot(
                plot_to_update=self.difference_scatter_echart,
                result=self.result3,
                title="Difference CNV",
                min_value="dataMin",
                preserve_zoom=preserve_zoom,
            )
        else:
            self._update_cnv_plot(
                plot_to_update=self.scatter_echart,
                result=self.result,
                gene_target=gene_target,
                title="CNV",
                preserve_zoom=preserve_zoom,
            )
            self._update_cnv_plot(
                plot_to_update=self.reference_scatter_echart,
                result=self.result2,
                gene_target=gene_target,
                title="Reference CNV",
                preserve_zoom=preserve_zoom,
            )
            self._update_cnv_plot(
                plot_to_update=self.difference_scatter_echart,
                result=self.result3,
                gene_target=gene_target,
                title="Difference CNV",
                min_value="dataMin",
                preserve_zoom=preserve_zoom,
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
        # Determine if this is the first plot (scatter_echart)
        is_first_plot = title == "CNV"

        # Set default y-axis limits
        if is_first_plot and self.y_axis_log:
            y_min = 0.1  # Minimum value for log scale (since log(0) is undefined)
            y_max = None  # Will be set dynamically when data is available
        else:
            y_min = initmin if initmin is not None else 0
            y_max = initmax if initmax is not None else None

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
                            "type": (
                                "log"
                                if (is_first_plot and self.y_axis_log)
                                else "value"
                            ),
                            "logBase": (
                                10 if (is_first_plot and self.y_axis_log) else None
                            ),
                            "name": "Ploidy",
                            "position": "left",
                            "min": y_min,
                            "max": y_max,  # Will be set dynamically
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
                            "yAxisIndex": [0, 1],  # Control both y-axes
                            "xAxisIndex": "none",
                            "filterMode": "none",
                            "width": 20,
                            "right": 50,
                            # "start":0,
                            # "end":8,
                            "startValue": y_min,
                            "endValue": y_max,  # Will be set dynamically
                            "minValueSpan": (
                                0.1 if (is_first_plot and self.y_axis_log) else None
                            ),
                            "maxValueSpan": None,  # Will be set dynamically
                            "borderColor": "#E5E5EA",
                            "backgroundColor": "#F5F5F7",
                            "fillerColor": "rgba(0, 122, 255, 0.2)",
                            "handleStyle": {"color": "#007AFF"},
                            "zoomLock": False,
                            "moveHandleSize": 7,
                            "zoomOnMouseWheel": True,
                            "moveOnMouseMove": True,
                            "preventDefaultMouseMove": True,
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
        min_value: Optional[float] = 0,
        ui_mode: bool = True,
        preserve_zoom: bool = False,
    ) -> None:
        """
        Update CNV plots with new data and annotations.

        Parameters
        ----------
        plot_to_update : ui.echart
            The chart component to update.
        result : Optional[Result], optional
            CNV result data. If None, uses self.result.
        gene_target : Optional[str], optional
            Target gene to focus the plot on.
        title : Optional[str], optional
            Title for the plot.
        min_value : Optional[float], default 0
            Minimum value for the x-axis.
        ui_mode : bool, default True
            Whether to update the UI immediately.
        preserve_zoom : bool, default False
            Whether to preserve the current zoom settings. If True, user zoom settings
            are maintained. If False, zoom is reset to show the full data range.
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
            min_value = min_value  # This can be modified later based on gene targeting
            max_value = (
                "dataMax"  # Placeholder to automatically determine the maximum value
            )

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
                min_value = start_pos - 10 * self.cnv_dict["bin_width"]
                max_value = end_pos + 10 * self.cnv_dict["bin_width"]

                # Further adjust the axis limits if the gene region is large
                if start_pos - min_value > 2_000_000:
                    min_value = start_pos - 2_000_000
                if max_value - end_pos > 2_000_000:
                    max_value = end_pos + 2_000_000

                if min_value < 0:
                    min_value = 0  # Ensure the minimum x-axis value is not negative

            # If all chromosomes are selected, prepare the plot to display all chromosomes
            if self.chrom_filter == "All":
                counter = 0  # Reset counter for chromosome loop

                # Update the plot title to reflect that all chromosomes are being shown
                plot_to_update.options["title"]["text"] = f"{title} - All Chromosomes"
                plot_to_update.options["series"] = []

                # Set legend display based on color mode - always hide it
                if "legend" not in plot_to_update.options:
                    plot_to_update.options["legend"] = {
                        "show": False,
                        "right": "10%",
                        "top": "10%",
                    }
                else:
                    plot_to_update.options["legend"]["show"] = False

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
                    plot_to_update.options["xAxis"]["max"] = max_value
                    plot_to_update.options["xAxis"]["min"] = min_value

                    # Append the current chromosome data as a scatter plot series
                    if self.color_mode == "chromosome":
                        # Original coloring method (by chromosome)
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
                                            "xAxis": (
                                                (total) * self.cnv_dict["bin_width"]
                                            ),
                                        },
                                    ],
                                },
                            }
                        )
                    else:
                        # Value-based coloring (red/blue)
                        # Determine expected ploidy value based on chromosome
                        expected_ploidy = 2  # Default for autosomes
                        if (
                            contig in ["chrX", "chrY"]
                            and self.sex_estimate
                            and self.sex_estimate in ["Male", "XY"]
                        ):
                            expected_ploidy = 1

                        # For relative difference plot
                        is_difference_plot = (
                            "Difference" in plot_to_update.options["title"]["text"]
                        )

                        # Calculate global statistics across all autosomes
                        if is_difference_plot:
                            # For difference plots, use 0 as the reference point
                            mean_val = 0
                            # Calculate global std across all autosomes
                            all_autosome_values = []
                            for chrom in result.cnv:
                                if (
                                    chrom.startswith("chr") and chrom[3:].isdigit()
                                ):  # Only autosomes
                                    all_autosome_values.extend(result.cnv[chrom])
                            std_val = (
                                np.std(all_autosome_values)
                                if all_autosome_values
                                else 1.0
                            )
                        else:
                            # For regular CNV plots, calculate global mean and std across autosomes
                            all_autosome_values = []
                            for chrom in result.cnv:
                                if (
                                    chrom.startswith("chr") and chrom[3:].isdigit()
                                ):  # Only autosomes
                                    all_autosome_values.extend(result.cnv[chrom])
                            mean_val = (
                                np.mean(all_autosome_values)
                                if all_autosome_values
                                else expected_ploidy
                            )
                            std_val = (
                                np.std(all_autosome_values)
                                if all_autosome_values
                                else 1.0
                            )

                        # Create three data arrays - one for values above threshold, one for below, and one for normal
                        data_above = []
                        data_below = []
                        data_normal = []

                        for pos, val in data:
                            # Calculate z-score using global statistics
                            z_score = (val - mean_val) / std_val if std_val > 0 else 0

                            if (
                                abs(z_score) > 0.5
                            ):  # Changed from 1.0 to 0.5 standard deviations
                                if z_score > 0:
                                    data_above.append([pos, val])
                                else:
                                    data_below.append([pos, val])
                            else:
                                data_normal.append([pos, val])

                        # Add blue points (values > 1 std above mean)
                        if data_above:
                            plot_to_update.options["series"].append(
                                {
                                    "type": "scatter",
                                    "name": f"Significant {('Gain' if is_difference_plot else 'High')} ({contig})",
                                    "data": data_above,
                                    "symbolSize": 5,
                                    "itemStyle": {"color": "#007AFF"},  # Blue color
                                }
                            )

                        # Add red points (values > 1 std below mean)
                        if data_below:
                            plot_to_update.options["series"].append(
                                {
                                    "type": "scatter",
                                    "name": f"Significant {('Loss' if is_difference_plot else 'Low')} ({contig})",
                                    "data": data_below,
                                    "symbolSize": 5,
                                    "itemStyle": {"color": "#FF3B30"},  # Red color
                                }
                            )

                        # Add gray points for normal values (within 1 std of mean)
                        if data_normal:
                            plot_to_update.options["series"].append(
                                {
                                    "type": "scatter",
                                    "name": f"Normal ({contig})",
                                    "data": data_normal,
                                    "symbolSize": 3,  # Slightly smaller for normal points
                                    "itemStyle": {
                                        "color": "#8E8E93"
                                    },  # iOS system gray
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
                    min_value = min_value
                    max_value = "dataMax"
                else:
                    # If the gene target is "All", adjust axis limits to show the entire chromosome
                    if gene_target == "All":
                        min_value = min_value
                        max_value = "dataMax"
                    else:
                        # Set axis limits to focus on the gene, with a margin of 10 times the bin width
                        min_value = start_pos - 10 * self.cnv_dict["bin_width"]
                        max_value = end_pos + 10 * self.cnv_dict["bin_width"]

                        # Further adjust the axis limits if the gene region is large
                        if start_pos - min_value > 2_000_000:
                            min_value = start_pos - 2_000_000
                        if max_value - end_pos > 2_000_000:
                            max_value = end_pos + 2_000_000

                        if min_value < 0:
                            min_value = (
                                0  # Ensure the minimum x-axis value is not negative
                            )

                # Update the plot title to reflect the selected chromosome
                if "Difference" in title:
                    plot_to_update.options["title"][
                        "text"
                    ] = f"Copy Number Variation (Relative Difference) - {contig}"
                else:
                    plot_to_update.options["title"][
                        "text"
                    ] = f"Copy Number Variation (Absolute) - {contig}"

                plot_to_update.options["xAxis"]["max"] = max_value
                plot_to_update.options["xAxis"]["min"] = min_value
                plot_to_update.options["dataZoom"][1]["endValue"] = ymax

                # Set legend display based on color mode - always hide it
                if "legend" not in plot_to_update.options:
                    plot_to_update.options["legend"] = {
                        "show": False,
                        "right": "10%",
                        "top": "10%",
                    }
                else:
                    plot_to_update.options["legend"]["show"] = False

                # Now initialize the series after we have contig defined
                plot_to_update.options["series"] = []

                # Prepare the data for plotting based on color mode
                if self.color_mode == "chromosome":
                    # Original coloring method
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
                else:
                    # Value-based coloring (red/blue)
                    # Determine expected ploidy value based on chromosome
                    expected_ploidy = 2  # Default for autosomes
                    if (
                        contig in ["chrX", "chrY"]
                        and self.sex_estimate
                        and self.sex_estimate in ["Male", "XY"]
                    ):
                        expected_ploidy = 1

                    # For relative difference plot
                    is_difference_plot = (
                        "Difference" in plot_to_update.options["title"]["text"]
                    )

                    # Calculate global statistics across all autosomes
                    if is_difference_plot:
                        # For difference plots, use 0 as the reference point
                        mean_val = 0
                        # Calculate global std across all autosomes
                        all_autosome_values = []
                        for chrom in result.cnv:
                            if (
                                chrom.startswith("chr") and chrom[3:].isdigit()
                            ):  # Only autosomes
                                all_autosome_values.extend(result.cnv[chrom])
                        std_val = (
                            np.std(all_autosome_values) if all_autosome_values else 1.0
                        )
                    else:
                        # For regular CNV plots, calculate global mean and std across autosomes
                        all_autosome_values = []
                        for chrom in result.cnv:
                            if (
                                chrom.startswith("chr") and chrom[3:].isdigit()
                            ):  # Only autosomes
                                all_autosome_values.extend(result.cnv[chrom])
                        mean_val = (
                            np.mean(all_autosome_values)
                            if all_autosome_values
                            else expected_ploidy
                        )
                        std_val = (
                            np.std(all_autosome_values) if all_autosome_values else 1.0
                        )

                    # Create three data arrays - one for values above threshold, one for below, and one for normal
                    data_above = []
                    data_below = []
                    data_normal = []

                    for pos, val in data:
                        # Calculate z-score using global statistics
                        z_score = (val - mean_val) / std_val if std_val > 0 else 0

                        if (
                            abs(z_score) > 0.25
                        ):  # Changed from 1.0 to 0.25 standard deviations
                            if z_score > 0:
                                data_above.append([pos, val])
                            else:
                                data_below.append([pos, val])
                        else:
                            data_normal.append([pos, val])

                    # Add blue points (values > 1 std above mean)
                    if data_above:
                        plot_to_update.options["series"].append(
                            {
                                "type": "scatter",
                                "name": f"Significant {('Gain' if is_difference_plot else 'High')} ({contig})",
                                "data": data_above,
                                "symbolSize": 5,
                                "itemStyle": {"color": "#007AFF"},  # Blue color
                            }
                        )

                    # Add red points (values > 1 std below mean)
                    if data_below:
                        plot_to_update.options["series"].append(
                            {
                                "type": "scatter",
                                "name": f"Significant {('Loss' if is_difference_plot else 'Low')} ({contig})",
                                "data": data_below,
                                "symbolSize": 5,
                                "itemStyle": {"color": "#FF3B30"},  # Red color
                            }
                        )

                    # Add gray points for normal values (within 1 std of mean)
                    if data_normal:
                        plot_to_update.options["series"].append(
                            {
                                "type": "scatter",
                                "name": f"Normal ({contig})",
                                "data": data_normal,
                                "symbolSize": 3,  # Slightly smaller for normal points
                                "itemStyle": {"color": "#8E8E93"},  # iOS system gray
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
                    try:
                        # Make sure we have a valid series with markArea property
                        if (
                            len(plot_to_update.options["series"]) > 0
                            and "markArea" in plot_to_update.options["series"][0]
                        ):
                            plot_to_update.options["series"][0]["markArea"][
                                "data"
                            ].append(
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
                    except (KeyError, IndexError) as e:
                        logger.warning(f"Error highlighting gene region: {e}")

                # Highlight the centromeres with shaded area:
                for _, row in self.centromere_bed[
                    self.centromere_bed["chrom"] == contig
                ].iterrows():
                    try:
                        # Make sure we have a valid series with markArea property
                        if (
                            len(plot_to_update.options["series"]) > 1
                            and "markArea" in plot_to_update.options["series"][1]
                        ):
                            plot_to_update.options["series"][1]["markArea"][
                                "data"
                            ].append(
                                [
                                    {
                                        "name": row["name"],
                                        "xAxis": row["start_pos"],
                                        "itemStyle": {
                                            "color": "rgba(255, 173, 177, 0.4)"
                                        },
                                        "label": {
                                            "position": "insideBottom",
                                            "distance": 10,
                                        },
                                    },
                                    {
                                        "xAxis": row["end_pos"],
                                    },
                                ]
                            )
                    except (KeyError, IndexError) as e:
                        logger.warning(f"Error highlighting centromere: {e}")

                # Add cytoband highlighting with CNV state colors
                # Check if cytobands data is available
                if self.cytobands_bed is None or self.cytobands_bed.empty:
                    logger.warning("Cytobands data not available for highlighting")
                else:
                    for _, row in self.cytobands_bed[
                        self.cytobands_bed["chrom"] == contig
                    ].iterrows():
                        try:
                            # Initialize cnv_state with a default value
                            cnv_state = "NORMAL"

                            # Verify that the row has the expected columns
                            if not all(
                                col in row.index
                                for col in ["start_pos", "end_pos", "name"]
                            ):
                                logger.warning(
                                    f"Cytoband row missing expected columns: {row.index.tolist()}"
                                )
                                continue

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

                            # Make sure we have a valid series with markArea property
                            if (
                                len(plot_to_update.options["series"]) > 2
                                and "markArea" in plot_to_update.options["series"][2]
                            ):
                                # Add colored regions to the plot
                                plot_to_update.options["series"][2]["markArea"][
                                    "data"
                                ].append(
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
                                                "fontWeight": cnv_state
                                                in ["GAIN", "LOSS"]
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
                        except (KeyError, IndexError) as e:
                            logger.warning(f"Error highlighting cytoband: {e}")
                            # Continue with next cytoband instead of breaking
                            continue

                    # Add horizontal reference lines for CNV thresholds only on the difference plot
                    if "Difference" in plot_to_update.options["title"]["text"]:
                        try:
                            # Make sure we have a valid series with markLine property
                            if (
                                len(plot_to_update.options["series"]) > 2
                                and "markLine" in plot_to_update.options["series"][2]
                            ):
                                plot_to_update.options["series"][2]["markLine"][
                                    "data"
                                ].extend(
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
                        except (KeyError, IndexError) as e:
                            logger.warning(f"Error adding reference lines: {e}")

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

            # Add breakpoint visualization data to the plot
            try:
                # Only show breakpoints if the toggle is enabled
                if not self.show_breakpoints:
                    breakpoint_series = []
                else:
                    # Determine which chromosome to show breakpoints for
                    target_chromosome = None
                    if self.chrom_filter != "All":
                        target_chromosome = valueslist[int(self.chrom_filter)]

                    # Get breakpoint visualization data
                    breakpoint_series = self.get_breakpoint_visualization_data(
                        plot_to_update, target_chromosome
                    )

                # Add breakpoint series to the plot
                for series in breakpoint_series:
                    plot_to_update.options["series"].append(series)

                # Update right y-axis limits for breakpoint data
                if breakpoint_series:
                    # Calculate max value for right y-axis
                    max_breakpoint_value = 0
                    for series in breakpoint_series:
                        if series["data"]:
                            series_max = max([point[1] for point in series["data"]])
                            max_breakpoint_value = max(max_breakpoint_value, series_max)

                    if max_breakpoint_value > 0:
                        # Update right y-axis limits
                        plot_to_update.options["yAxis"][1].update(
                            {
                                "min": 0,
                                "max": max_breakpoint_value * 1.2,  # Add 20% padding
                            }
                        )

            except Exception as e:
                logger.warning(f"Error adding breakpoint visualization: {e}")

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

        # Calculate y-axis limits based on data
        if result and result.cnv:
            all_values = []
            for chrom, cnv in result.cnv.items():
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    all_values.extend(cnv)

            if all_values:
                # Determine if this plot should use log scale
                # Only apply log scale to the absolute CNV plot (scatter_echart), not the difference plot
                is_absolute_plot = plot_to_update == self.scatter_echart
                should_use_log = is_absolute_plot and self.y_axis_log

                # Calculate y-axis limits
                data_min = (
                    max(0.1, min(all_values)) if should_use_log else min(all_values)
                )
                data_max = max(all_values)

                # Add some padding to the max value and round up to nearest integer
                if should_use_log:
                    # For log scale, multiply by a factor and round up
                    y_max = math.ceil(data_max * 1.2)
                    # For log scale, we keep data_min as is since it's a small decimal
                else:
                    # For linear scale, add a percentage and round up
                    y_max = math.ceil(data_max * 1.1)
                    # For linear scale, round down to nearest integer
                    data_min = math.floor(data_min)

                # Update y-axis limits
                if not preserve_zoom:
                    plot_to_update.options["yAxis"][0].update(
                        {"min": data_min, "max": y_max}
                    )

                # Update zoom settings only if not preserving zoom
                if not preserve_zoom:
                    plot_to_update.options["dataZoom"][1].update(
                        {
                            "startValue": data_min,
                            "endValue": y_max,
                            "maxValueSpan": y_max - data_min,
                        }
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
            Estimated genetic sex (XX/Female, XY/Male, or Unknown)
        bin_width : int
            Current bin width used in the analysis
        variance : float
            Current variance value from the analysis
        """
        # Convert old format to new if necessary
        display_sex = xy_estimate
        if xy_estimate == "XX":
            display_sex = "Female (XX)"
        elif xy_estimate == "XY":
            display_sex = "Male (XY)"
        else:
            display_sex = xy_estimate

        # Determine if male or female for icon and styling
        is_male = xy_estimate.split(" ")[0] in ["XY", "Male"]

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
                                    "text-blue-600" if is_male else "text-pink-600"
                                )
                                # status_bg = "bg-blue-100" if is_male else "bg-pink-100"
                                if is_male:
                                    ui.icon("man").classes("text-4xl text-blue-500")
                                else:
                                    ui.icon("woman").classes("text-4xl text-pink-500")
                                ui.label(f"Genetic Sex: {display_sex}").classes(
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
            self.result, self.result2, self.result3, self.cnv_dict = await run.io_bound(
                get_data,
                output,
            )

            self.CNVResults = {}

            async def load_ruptures():
                if self.check_file_time(os.path.join(output, "ruptures.npy")):
                    self.CNVResults = await run.io_bound(
                        np.load,
                        os.path.join(output, "ruptures.npy"),
                        allow_pickle="TRUE",
                    ).item()

                # Load breakpoint data for visualization
                data_array_path = os.path.join(output, "cnv_data_array.npy")
                if os.path.exists(data_array_path):
                    try:
                        # Check file size first
                        file_size = os.path.getsize(data_array_path)
                        logger.info(f"Breakpoint file size: {file_size} bytes")

                        if file_size == 0:
                            logger.warning("Breakpoint file is empty")
                            self.DATA_ARRAY = None
                        else:
                            # Try loading as standard numpy array
                            try:
                                logger.info(
                                    "Attempting to load as standard numpy array..."
                                )
                                self.DATA_ARRAY = np.load(
                                    data_array_path, allow_pickle=True
                                )
                                logger.info(
                                    "Successfully loaded as standard numpy array"
                                )
                            except Exception as array_error:
                                logger.warning(
                                    f"Failed to load breakpoint data: {array_error}"
                                )
                                self.DATA_ARRAY = None

                        # Filter out empty entries (where name is empty string)
                        if self.DATA_ARRAY is not None and len(self.DATA_ARRAY) > 0:
                            # Check if it has the expected structure
                            if (
                                hasattr(self.DATA_ARRAY, "dtype")
                                and "name" in self.DATA_ARRAY.dtype.names
                            ):
                                # Remove entries with empty names
                                valid_mask = self.DATA_ARRAY["name"] != ""
                                self.DATA_ARRAY = self.DATA_ARRAY[valid_mask]
                                logger.info(
                                    f"Loaded breakpoint data: {len(self.DATA_ARRAY)} breakpoints"
                                )
                            else:
                                logger.warning(
                                    "Breakpoint data doesn't have expected structure"
                                )
                                self.DATA_ARRAY = None
                        else:
                            logger.info("No breakpoint data found")
                            self.DATA_ARRAY = None

                    except Exception as e:
                        logger.warning(f"Failed to load breakpoint data: {e}")
                        self.DATA_ARRAY = None
                else:
                    logger.info("Breakpoint data file does not exist")
                    self.DATA_ARRAY = None
                self.update_plots(preserve_zoom=True)

            background_tasks.create(load_ruptures())

            async def load_bedranges():
                if self.check_file_time(os.path.join(output, "bedranges.csv")):
                    self.proportions_df_store = await run.io_bound(
                        pd.read_csv,
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

            await background_tasks.create(load_bedranges())

            if self.summary:
                with self.summary:
                    self.summary.clear()
                    with open(os.path.join(output, "XYestimate.pkl"), "rb") as file:
                        xy_estimate = pickle.load(file)
                        self.sex_estimate = (
                            xy_estimate  # Store the loaded value in sex_estimate
                        )
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

        # Pre-allocate merged_cytobands with a reasonable size
        max_expected_cytobands = len(chromosome_cytobands)
        merged_cytobands = [None] * max_expected_cytobands
        merged_idx = 0
        whole_chr_event = False
        whole_chr_state = "NORMAL"

        # First, analyze the whole chromosome for potential aneuploidy
        if chromosome in cnv_data:
            logger.debug(
                f"\nAnalyzing chromosome {chromosome} for whole chromosome events:"
            )

            # Use boolean mask instead of concatenation
            mask = np.ones(len(cnv_data[chromosome]), dtype=bool)
            centromere = self.centromere_bed[self.centromere_bed["chrom"] == chromosome]
            if not centromere.empty:
                cent_start_bin = int(centromere["start_pos"].iloc[0] / bin_width)
                cent_end_bin = int(centromere["end_pos"].iloc[0] / bin_width)
                mask[cent_start_bin:cent_end_bin] = False
                logger.debug(
                    f"Excluded centromere region: {cent_start_bin}-{cent_end_bin}"
                )
            chr_cnv = cnv_data[chromosome][mask]

            # Calculate chromosome-wide statistics
            chr_mean = np.mean(chr_cnv)
            chr_std = np.std(chr_cnv)
            logger.debug(f"Chromosome-wide mean: {chr_mean:.3f}, std: {chr_std:.3f}")

            # Calculate SD of chromosome means for whole chromosome event detection
            chromosome_means = []
            for chrom in cnv_data:
                if chrom.startswith("chr") and chrom[3:].isdigit():  # Only autosomes
                    # Use mask for centromere exclusion
                    mask = np.ones(len(cnv_data[chrom]), dtype=bool)
                    cent = self.centromere_bed[self.centromere_bed["chrom"] == chrom]
                    if not cent.empty:
                        cent_start = int(cent["start_pos"].iloc[0] / bin_width)
                        cent_end = int(cent["end_pos"].iloc[0] / bin_width)
                        mask[cent_start:cent_end] = False
                    chrom_data = cnv_data[chrom][mask]
                    if len(chrom_data) > 0:
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
                if self.sex_estimate in ["Male", "XY"]:  # Male
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
                if self.sex_estimate in ["Male", "XY"]:  # Male
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

                merged_cytobands[merged_idx] = {
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
                merged_idx += 1

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
                        # For whole chromosome events, only report significant deviations in opposite direction
                        if whole_chr_state == "GAIN":
                            if mean_cnv < chr_mean - (
                                2.0 * chr_std
                            ):  # More stringent threshold for opposite changes
                                state = "LOSS"  # Only report significant losses within gained chromosomes
                            else:
                                state = "NORMAL"  # Don't duplicate gains
                        elif whole_chr_state == "LOSS":
                            if mean_cnv > chr_mean + (
                                2.0 * chr_std
                            ):  # More stringent threshold for opposite changes
                                state = "GAIN"  # Only report significant gains within lost chromosomes
                            else:
                                state = "NORMAL"  # Don't duplicate losses
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

                    # Only add significant changes relative to whole chromosome state
                    if (not whole_chr_event) or (
                        current_group["cnv_state"] != whole_chr_state
                    ):
                        merged_cytobands[merged_idx] = current_group
                        merged_idx += 1

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

                # Only add significant changes relative to whole chromosome state
                if (not whole_chr_event) or (
                    current_group["cnv_state"] != whole_chr_state
                ):
                    merged_cytobands[merged_idx] = current_group
                    merged_idx += 1

        # Trim the pre-allocated list to actual size
        merged_cytobands = merged_cytobands[:merged_idx]

        # Convert to DataFrame and sort
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

    def toggle_y_axis_scale(self, e):
        """Update y-axis scale and zoom settings when toggle changes."""
        self.y_axis_log = e.value == "log"

        # Get current data to calculate new limits
        if hasattr(self, "result") and self.result and self.result.cnv:
            all_values = []
            for chrom, cnv in self.result.cnv.items():
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    all_values.extend(cnv)

            if all_values:
                # Calculate new y-axis limits for absolute plot only
                data_min = (
                    max(0.1, min(all_values)) if self.y_axis_log else min(all_values)
                )
                data_max = max(all_values)

                # Add padding to max value and round up to nearest integer
                if self.y_axis_log:
                    y_max = math.ceil(data_max * 1.2)
                    # For log scale, we keep data_min as is since it's a small decimal
                else:
                    y_max = math.ceil(data_max * 1.1)
                    # For linear scale, round down to nearest integer
                    data_min = math.floor(data_min)

                # Update y-axis settings for absolute plot only
                self.scatter_echart.options["yAxis"][0].update(
                    {
                        "type": "log" if self.y_axis_log else "value",
                        "logBase": 10 if self.y_axis_log else None,
                        "min": data_min,
                        "max": y_max,
                    }
                )

                # Update zoom settings for absolute plot only
                self.scatter_echart.options["dataZoom"][1].update(
                    {
                        "startValue": data_min,
                        "endValue": y_max,
                        "minValueSpan": 0.1 if self.y_axis_log else None,
                        "maxValueSpan": y_max - data_min,
                    }
                )

                ui.update(self.scatter_echart)

        # Trigger a full update of all plots to ensure difference plot stays in linear scale
        self.update_plots(preserve_zoom=True)

    def calculate_breakpoint_metrics(self, chromosome: str = None) -> dict:
        """
        Calculate breakpoint metrics for visualization on the right y-axis.

        Args:
            chromosome (str): Specific chromosome to analyze, or None for all

        Returns:
            dict: Dictionary containing breakpoint metrics for visualization
        """
        if not hasattr(self, "DATA_ARRAY") or self.DATA_ARRAY is None:
            logger.debug("No DATA_ARRAY available for breakpoint metrics")
            return {}

        logger.debug(
            f"DATA_ARRAY type: {type(self.DATA_ARRAY)}, length: {len(self.DATA_ARRAY) if self.DATA_ARRAY is not None else 0}"
        )

        metrics = {
            "density": [],  # Breakpoint density per bin
            "confidence": [],  # Confidence scores
            "positions": [],  # Breakpoint positions
            "timeline": [],  # Detection timeline
        }

        # Get breakpoints for the specified chromosome or all chromosomes
        if chromosome:
            breakpoints = self.DATA_ARRAY[self.DATA_ARRAY["name"] == chromosome]
            logger.debug(
                f"Filtering for chromosome {chromosome}: {len(breakpoints)} breakpoints"
            )
        else:
            breakpoints = self.DATA_ARRAY
            logger.debug(f"Using all breakpoints: {len(breakpoints)} breakpoints")

        if len(breakpoints) == 0:
            logger.debug("No breakpoints found for visualization")
            return metrics

        # Calculate breakpoint density across the genome
        bin_width = self.cnv_dict.get("bin_width", 1_000_000)
        logger.debug(f"Using bin width: {bin_width}")

        # Group breakpoints by chromosome
        unique_chroms = np.unique(breakpoints["name"])
        logger.debug(f"Unique chromosomes in breakpoints: {unique_chroms}")

        for chrom in unique_chroms:
            chrom_breakpoints = breakpoints[breakpoints["name"] == chrom]
            logger.debug(
                f"Processing chromosome {chrom}: {len(chrom_breakpoints)} breakpoints"
            )

            if len(chrom_breakpoints) == 0:
                continue

            # Calculate density per bin
            max_pos = np.max(chrom_breakpoints["end"])
            num_bins = int(max_pos / bin_width) + 1
            logger.debug(f"Chromosome {chrom}: max_pos={max_pos}, num_bins={num_bins}")

            # Initialize density array
            density = np.zeros(num_bins)

            # Count breakpoints in each bin
            for bp in chrom_breakpoints:
                start_bin = int(bp["start"] / bin_width)
                end_bin = int(bp["end"] / bin_width)

                # Ensure bins are within bounds
                start_bin = max(0, min(start_bin, num_bins - 1))
                end_bin = max(0, min(end_bin, num_bins - 1))

                # Add to density (weighted by breakpoint size)
                for bin_idx in range(start_bin, end_bin + 1):
                    density[bin_idx] += 1

            # Convert to positions and values for plotting
            for i, count in enumerate(density):
                if count > 0:
                    pos = i * bin_width
                    metrics["density"].append([pos, count])
                    metrics["positions"].append(pos)

        logger.debug(
            f"Generated metrics - density: {len(metrics['density'])}, positions: {len(metrics['positions'])}"
        )

        # Calculate confidence scores (simplified - based on breakpoint frequency)
        if metrics["density"]:
            max_density = max([d[1] for d in metrics["density"]])
            for pos, density in metrics["density"]:
                confidence = min(1.0, density / max_density) if max_density > 0 else 0
                metrics["confidence"].append([pos, confidence])

        return metrics

    def get_breakpoint_visualization_data(
        self, plot_to_update: ui.echart, chromosome: str = None
    ) -> list:
        """
        Generate breakpoint visualization data for the right y-axis.

        Args:
            plot_to_update: The chart to update
            chromosome (str): Specific chromosome to visualize

        Returns:
            list: List of series data for breakpoint visualization
        """
        logger.debug(
            f"Generating breakpoint visualization for chromosome: {chromosome}"
        )
        metrics = self.calculate_breakpoint_metrics(chromosome)
        series_data = []

        logger.debug(f"Breakpoint metrics: {metrics}")

        if not metrics or not any(metrics.values()):
            logger.debug("No breakpoint metrics available for visualization")
            return series_data

        # Show only points at bins with breakpoints (scatter plot)
        if metrics["density"]:
            unique_breakpoints = sorted(set(metrics["positions"]))
            series_data.append(
                {
                    "type": "scatter",
                    "name": "Breakpoint Density",
                    "yAxisIndex": 1,  # Use right y-axis
                    "data": metrics["density"],
                    "symbolSize": 10,
                    "itemStyle": {
                        "color": "rgba(255, 69, 0, 0.7)",
                        "borderColor": "#FF4500",
                        "borderWidth": 2,
                    },
                    "emphasis": {
                        "itemStyle": {
                            "color": "rgba(255, 69, 0, 0.9)",
                            "borderWidth": 2,
                        }
                    },
                    "markLine": {
                        "symbol": "none",
                        "data": [
                            {
                                "xAxis": pos,
                                "lineStyle": {
                                    "color": "#FF0000",
                                    "type": "dashed",
                                    "width": 2,
                                },
                            }
                            for pos in unique_breakpoints
                        ],
                    },
                }
            )

        logger.debug(f"Generated {len(series_data)} breakpoint series")
        return series_data


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
        self.cnv_dict = {"bin_width": 0, "variance": 0}
        self.update_cnv_dict = {}
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
        self.DATA_ARRAY = None
        self.data_array_size = 0
        self.data_array_path = None
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

        # Initialize CNV detector with proper temp directory
        self.CNVchangedetector = CNVChangeDetectorTracker(
            base_proportion=0.02,
            temp_dir=self.output,  # Will use system temp dir by default
        )

        # Add target_table as instance variable
        self.target_table = None
        self.target_table_placeholder = None
        self.last_bed_check = 0  # Track when we last checked for new BED files

        # Add y_axis_log attribute
        self.y_axis_log = False  # Default to linear scale
        # Add show_breakpoints attribute
        self.show_breakpoints = True  # Default to showing breakpoints

        # Add DATA_ARRAY for breakpoint visualization
        self.DATA_ARRAY = None
        self.data_array_path = None

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
        """Clean up temporary files and state."""
        if self.CNVchangedetector:
            self.CNVchangedetector.storage.clear()
        if self.state_dir and self.state_dir.exists():
            try:
                shutil.rmtree(self.state_dir)
            except Exception as e:
                logger.warning(f"Failed to clean up state directory: {e}")

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
        """
        Detect whole chromosome copy number events.

        Args:
            z_score_threshold (float): Z-score threshold for detecting events.

        Returns:
            List[Dict]: List of detected chromosome events.
        """
        # Initialize list to store events
        events = []

        # Get CNV data
        cnv_data = self.result3.cnv

        # Calculate means for all chromosomes
        chromosome_means = {}
        for chrom, data in cnv_data.items():
            if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                # Exclude centromere region if present
                mask = np.ones_like(data, dtype=bool)
                cent = self.centromere_bed[self.centromere_bed["chrom"] == chrom]
                if not cent.empty:
                    cent_start = int(
                        cent["start_pos"].iloc[0] / self.cnv_dict["bin_width"]
                    )
                    cent_end = int(cent["end_pos"].iloc[0] / self.cnv_dict["bin_width"])
                    if cent_start < len(mask) and cent_end <= len(mask):
                        mask[cent_start:cent_end] = False
                filtered_data = data[mask]

                # Handle empty arrays
                if len(filtered_data) > 0:
                    chromosome_means[chrom] = np.mean(filtered_data)

        # Calculate mean and standard deviation of autosomal chromosomes
        autosomal_means = [
            v
            for k, v in chromosome_means.items()
            if k.startswith("chr") and k[3:].isdigit()
        ]
        if autosomal_means:
            mean_of_means = np.mean(autosomal_means)
            std_of_means = np.std(autosomal_means)
        else:
            mean_of_means = 0
            std_of_means = 0

        # Detect events for each chromosome
        for chrom, mean_cnv in chromosome_means.items():
            # Handle sex chromosomes specially
            if chrom == "chrX":
                if self.sex_estimate in ["Female", "XX"]:  # Female
                    # For females, X should be treated like autosomes
                    z_score = (mean_cnv - mean_of_means) / std_of_means
                    if z_score > z_score_threshold:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "GAIN",
                                "description": f"Chromosome X gain (z-score: {z_score:.2f})",
                            }
                        )
                    elif z_score < -z_score_threshold:
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "LOSS",
                                "description": f"Chromosome X loss (z-score: {z_score:.2f})",
                            }
                        )
                elif self.sex_estimate in ["Male", "XY"]:  # Male
                    # For males, X should have one copy
                    if mean_cnv > 0.5:  # Threshold for X gain in males
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "GAIN",
                                "description": "Chromosome X gain (expected single copy in males)",
                            }
                        )
                    elif mean_cnv < -0.5:  # Threshold for X loss in males
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "LOSS",
                                "description": "Chromosome X loss (expected single copy in males)",
                            }
                        )
            elif chrom == "chrY":
                if self.sex_estimate in ["Male", "XY"]:  # Male
                    # For males, Y should have one copy
                    if mean_cnv > 0.5:  # Threshold for Y gain in males
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "GAIN",
                                "description": "Chromosome Y gain (expected single copy in males)",
                            }
                        )
                    elif mean_cnv < -0.5:  # Threshold for Y loss in males
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "LOSS",
                                "description": "Chromosome Y loss (expected single copy in males)",
                            }
                        )
                elif self.sex_estimate in ["Female", "XX"]:  # Female
                    # For females, Y should not be present
                    if mean_cnv > 0.3:  # Threshold for Y presence in females
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "GAIN",
                                "description": "Unexpected Y chromosome presence (expected absent in females)",
                            }
                        )
                    elif mean_cnv > -0.2 and self.sex_estimate in ["Female", "XX"]:
                        # Not quite a Y gain but stronger signal than expected for females
                        events.append(
                            {
                                "chromosome": chrom,
                                "mean_cnv": mean_cnv,
                                "type": "NORMAL",
                                "description": "Slight Y chromosome signal (check sample quality)",
                            }
                        )
            else:  # Autosomal chromosomes
                # Calculate z-score relative to other chromosomes
                z_score = (mean_cnv - mean_of_means) / std_of_means
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
            if self.CNVchangedetector is None:
                self.CNVchangedetector = CNVChangeDetectorTracker(
                    base_proportion=0.02, temp_dir=self.output
                )

            # Load previous state if it exists
            if self.load_data and (state_dir / "tracker_metadata.pkl").exists():
                try:
                    self.CNVchangedetector.load_state(str(state_dir))
                    logger.info("Successfully loaded CNV detector state")
                except Exception as e:
                    logger.warning(f"Failed to load CNV detector state: {e}")
                    # Initialize fresh state if loading fails
                    if self.CNVchangedetector is None:
                        self.CNVchangedetector = CNVChangeDetectorTracker(
                            base_proportion=0.02, temp_dir=self.output
                        )

            # Initialize data structures
            if self.sampleID not in self.update_cnv_dict.keys():
                self.update_cnv_dict[self.sampleID] = {}

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
            self.data_array_path = os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "cnv_data_array.npy",
            )

            # Initialize or load data array using standard file operations
            if not os.path.exists(self.data_array_path):
                # Initialize with empty array
                self.DATA_ARRAY = np.array([], dtype=self.dtype)
                self.data_array_size = 0
            else:
                # Load existing data
                try:
                    self.DATA_ARRAY = np.load(self.data_array_path, allow_pickle=True)
                    self.data_array_size = len(self.DATA_ARRAY)
                except Exception as e:
                    logger.warning(f"Failed to load existing data array: {e}")
                    self.DATA_ARRAY = np.array([], dtype=self.dtype)
                    self.data_array_size = 0

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
                self.update_cnv_dict[self.sampleID],
                genome_length,
                r2_cnv,
            ) = await run.cpu_bound(
                iterate_bam,
                bamfile,
                self.threads,
                60,
                self.update_cnv_dict[self.sampleID],
                int(logging.ERROR),
                self.ref_cnv_dict,
            )

            self.cnv_dict["bin_width"] = r_bin
            self.cnv_dict["variance"] = r_var

            # Process breakpoints and update CNV detector
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

                    if len(r_cnv[key]) > 3 and self.map_tracker[key] > 2000:
                        paired_changepoints = await run.cpu_bound(
                            run_ruptures,
                            r_cnv[key],
                            5,  # penalty_value
                            self.cnv_dict["bin_width"],
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
                                self.local_data_array = np.append(
                                    self.local_data_array, item
                                )
                                # Append to persistent data array (all breakpoints for this sample)
                                self.DATA_ARRAY = np.append(self.DATA_ARRAY, item)

                                # Save the updated data array to file immediately
                                try:
                                    np.save(self.data_array_path, self.DATA_ARRAY)
                                except Exception as e:
                                    logger.warning(f"Failed to save data array: {e}")

                            # Reset the mapping tracker for this chromosome to prevent reprocessing
                            self.map_tracker[key] = 0
                            # Flag that we have CNV updates to save
                            cnvupdate = True

                            # Extract all breakpoints for this chromosome from persistent storage
                            breakpoints = self.DATA_ARRAY[
                                self.DATA_ARRAY["name"] == key
                            ]
                            # Extract breakpoints from current processing batch
                            local_breakpoints = self.local_data_array[
                                self.local_data_array["name"] == key
                            ]

                            # Only proceed if we have valid breakpoints to process
                            if len(breakpoints) > 0:
                                try:
                                    # Add the detected breakpoints to the CNV change detector
                                    # This updates the internal tracking and generates BED targets
                                    self.CNVchangedetector.add_breakpoints(
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

                    self.CNVchangedetector.save_state(str(state_dir))
                    logger.info("Successfully saved CNV detector state")

                    # Generate and save BED content
                    bedcontent = ""
                    bedcontent2 = ""
                    for chrom in r_cnv.keys():
                        tempbedcontent = self.CNVchangedetector.get_bed_targets(chrom)
                        if len(tempbedcontent) > 0:
                            bedcontent += tempbedcontent + "\n"

                        tempbedcontent2 = (
                            self.CNVchangedetector.get_bed_targets_breakpoints(chrom)
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
                        self.cnv_dict,
                    )

                except Exception as e:
                    logger.error(f"Error saving CNV detector state: {e}")
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

        # Pre-allocate merged_cytobands with a reasonable size
        max_expected_cytobands = len(chromosome_cytobands)
        merged_cytobands = [None] * max_expected_cytobands
        merged_idx = 0
        whole_chr_event = False
        whole_chr_state = "NORMAL"

        # First, analyze the whole chromosome for potential aneuploidy
        if chromosome in cnv_data:
            logger.debug(
                f"\nAnalyzing chromosome {chromosome} for whole chromosome events:"
            )

            # Use boolean mask instead of concatenation
            mask = np.ones(len(cnv_data[chromosome]), dtype=bool)
            centromere = self.centromere_bed[self.centromere_bed["chrom"] == chromosome]
            if not centromere.empty:
                cent_start_bin = int(centromere["start_pos"].iloc[0] / bin_width)
                cent_end_bin = int(centromere["end_pos"].iloc[0] / bin_width)
                mask[cent_start_bin:cent_end_bin] = False
                logger.debug(
                    f"Excluded centromere region: {cent_start_bin}-{cent_end_bin}"
                )
            chr_cnv = cnv_data[chromosome][mask]

            # Calculate chromosome-wide statistics
            chr_mean = np.mean(chr_cnv)
            chr_std = np.std(chr_cnv)
            logger.debug(f"Chromosome-wide mean: {chr_mean:.3f}, std: {chr_std:.3f}")

            # Calculate SD of chromosome means for whole chromosome event detection
            chromosome_means = []
            for chrom in cnv_data:
                if chrom.startswith("chr") and chrom[3:].isdigit():  # Only autosomes
                    # Use mask for centromere exclusion
                    mask = np.ones(len(cnv_data[chrom]), dtype=bool)
                    cent = self.centromere_bed[self.centromere_bed["chrom"] == chrom]
                    if not cent.empty:
                        cent_start = int(cent["start_pos"].iloc[0] / bin_width)
                        cent_end = int(cent["end_pos"].iloc[0] / bin_width)
                        mask[cent_start:cent_end] = False
                    chrom_data = cnv_data[chrom][mask]
                    if len(chrom_data) > 0:
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
                if self.sex_estimate in ["Male", "XY"]:  # Male
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
                if self.sex_estimate in ["Male", "XY"]:  # Male
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

                merged_cytobands[merged_idx] = {
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
                merged_idx += 1

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
                        # For whole chromosome events, only report significant deviations in opposite direction
                        if whole_chr_state == "GAIN":
                            if mean_cnv < chr_mean - (
                                2.0 * chr_std
                            ):  # More stringent threshold for opposite changes
                                state = "LOSS"  # Only report significant losses within gained chromosomes
                            else:
                                state = "NORMAL"  # Don't duplicate gains
                        elif whole_chr_state == "LOSS":
                            if mean_cnv > chr_mean + (
                                2.0 * chr_std
                            ):  # More stringent threshold for opposite changes
                                state = "GAIN"  # Only report significant gains within lost chromosomes
                            else:
                                state = "NORMAL"  # Don't duplicate losses
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

                    # Only add significant changes relative to whole chromosome state
                    if (not whole_chr_event) or (
                        current_group["cnv_state"] != whole_chr_state
                    ):
                        merged_cytobands[merged_idx] = current_group
                        merged_idx += 1

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

                # Only add significant changes relative to whole chromosome state
                if (not whole_chr_event) or (
                    current_group["cnv_state"] != whole_chr_state
                ):
                    merged_cytobands[merged_idx] = current_group
                    merged_idx += 1

        # Trim the pre-allocated list to actual size
        merged_cytobands = merged_cytobands[:merged_idx]

        # Convert to DataFrame and sort
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
