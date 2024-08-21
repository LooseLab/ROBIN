"""
Copy Number Variation (CNV) Analysis Module

This module provides functionality for analyzing copy number variations (CNVs) from BAM files generated during sequencing runs. It includes classes and functions for iterating over BAM files, processing CNV data, and visualizing results using NiceGUI for web-based user interfaces.

Key Components:

1. **Result Class**:
   - A simple class to store CNV results.

2. **Helper Functions**:
   - `iterate_bam`: Iterates over a BAM file and returns CNV data and associated metrics.
   - `iterate_bam_bin`: Iterates over a BAM file with specified bin width and returns CNV data and associated metrics.
   - `reduce_list`: Reduces the length of a list to a specified maximum length by subsampling.
   - `moving_average`: Calculates the moving average of a given data array.
   - `pad_arrays`: Pads two arrays to the same length with a specified value.
   - `ruptures_plotting`: Applies the Kernel Change Point Detection (CPD) algorithm to identify change points in data.

3. **CNVAnalysis Class**:
   - Inherits from `BaseAnalysis` and provides specific methods for CNV analysis.
   - Initializes with specific target panel information and loads reference data.
   - Implements methods to estimate genetic sex (XX or XY) based on CNV data.
   - Processes BAM files to extract CNV data and updates visualizations.
   - Provides a user interface for visualizing CNV data using NiceGUI.

4. **User Interface and Visualization**:
   - `setup_ui`: Sets up the user interface for CNV analysis, including selection controls and plots.
   - `generate_chart`: Generates ECharts objects for displaying CNV scatter plots and difference plots.
   - `_update_cnv_plot`: Updates CNV plots with new data and annotations based on selected chromosomes and genes.
   - `show_previous_data`: Loads and displays previously computed CNV data.

5. **Command-Line Interface**:
   - Uses Click for defining command-line options and arguments to run the CNV analysis application.
   - `test_me`: Helper function to initialize and run the CNV analysis application.
   - `main`: Entry point for the command-line interface to start the CNV analysis application.

Dependencies:
- `cnv_from_bam.iterate_bam_file`
- `robin.subpages.base_analysis.BaseAnalysis`
- `natsort`
- `pandas`
- `logging`
- `numpy`
- `os`
- `sys`
- `nicegui` (ui, app, run)
- `click`
- `pathlib.Path`
- `pickle`
- `ruptures`

Example usage::
    @click.command()
    @click.option("--port", default=12345, help="Port for GUI")
    @click.option("--threads", default=4, help="Number of threads available.")
    @click.argument("watchfolder", type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path), required=False)
    @click.argument("output", type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path), required=False)
    @click.option("--browse", is_flag=True, show_default=True, default=False, help="Browse Historic Data.")
    @click.option("--target_panel", "-t", default="rCNS2", help="Select analysis gene panel from one of these options. Default is rCNS2", type=click.Choice(["rCNS2", "AML"], case_sensitive=True))
    def main(port, threads, watchfolder, output, browse, target_panel):
        # Run the CNV analysis application
        pass

if __name__ in {"__main__", "__mp_main__"}:
    main()
"""

from cnv_from_bam import iterate_bam_file
from robin.subpages.base_analysis import BaseAnalysis
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
from collections import defaultdict

os.environ["CI"] = "1"
# Use the main logger configured in the main application
logger = logging.getLogger(__name__)



class Result:
    """
    A class to store CNV results.
    """

    def __init__(self, cnv_dict: dict) -> None:
        self.cnv = cnv_dict


def run_ruptures(r_cnv, penalty_value, bin_width):
    algo_c = ruptures_plotting(np.array(r_cnv))

    ruptures_result = algo_c.predict(pen=penalty_value)

    scaled_changepoints = [cp * bin_width for cp in ruptures_result]

    return [(cp - bin_width, cp + bin_width) for cp in scaled_changepoints]


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

    return result.cnv, result.bin_width, result.variance, copy_numbers, result.genome_length


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
    len_diff = abs(len(arr1) - len(arr2))
    if len(arr1) < len(arr2):
        arr1 = np.pad(arr1, (0, len_diff), mode="constant", constant_values=pad_value)
    elif len(arr1) > len(arr2):
        arr2 = np.pad(arr2, (0, len_diff), mode="constant", constant_values=pad_value)
    return arr1, arr2


def ruptures_plotting(data: np.ndarray) -> rpt.KernelCPD:
    """
    Apply the Kernel Change Point Detection (CPD) algorithm to identify change points in data.

    Args:
        data (np.ndarray): The data array.

    Returns:
        rpt.KernelCPD: The change point detection algorithm object.
    """
    #x_coords = range(0, (data.size + 1))
    #signal = np.array(list(zip(x_coords, data)))
    #try:
    #    algo_c = rpt.KernelCPD(kernel="linear", min_size=10).fit(signal)
    #except Exception as e:
    #    print(f"e")

    algo_c = rpt.KernelCPD(kernel="rbf").fit(data) #, min_size=10)

    return algo_c


class CNVAnalysis(BaseAnalysis):
    """
    Class for analyzing copy number variations (CNVs) from BAM files.

    Inherits from `BaseAnalysis` and provides specific methods for CNV analysis.
    """

    def __init__(self, *args, target_panel: Optional[str] = None, **kwargs) -> None:
        # self.file_list = []
        self.cnv_dict = {"bin_width": 0, "variance": 0}
        self.update_cnv_dict = {}
        self.result = None
        self.result2 = None
        self.result3 = CNV_Difference()
        self.XYestimate = "Unknown"
        self.dtype = [('name', 'U10'), ('start', 'i8'), ('end', 'i8')]
        self.DATA_ARRAY = np.empty(0, dtype=self.dtype)

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
            usecols=[0,1,2,3],
            names=["chrom", "start_pos", "end_pos", "name"],
            header=None,
            sep="\s+",
        )
        super().__init__(*args, **kwargs)

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
        #try:
        await self.do_cnv_work(bamfile)
        #except Exception as e:
        #    logger.error(e)
        #    logger.error("line 313")

    async def do_cnv_work(self, bamfile: BinaryIO) -> None:
        """
        Perform CNV analysis on a BAM file.

        Args:
            bamfile (BinaryIO): The BAM file to process.
        """
        if self.sampleID not in self.update_cnv_dict.keys():
            self.update_cnv_dict[self.sampleID] = {}
        # print("Running CNV analysis on", self.sampleID, bamfile)
        r_cnv, r_bin, r_var, self.update_cnv_dict[self.sampleID], genome_length = await run.cpu_bound(
            iterate_bam,
            bamfile,
            self.threads,
            60,
            self.update_cnv_dict[self.sampleID],
            int(logging.ERROR),
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


        for key in r_cnv.keys():
            if key != "chrM" and re.match(r'^chr(\d+|X|Y)$', key):
                moving_avg_data1 = await run.cpu_bound(moving_average, r_cnv[key])
                moving_avg_data2 = await run.cpu_bound(moving_average, r2_cnv[key])
                moving_avg_data1, moving_avg_data2 = await run.cpu_bound(
                    pad_arrays, moving_avg_data1, moving_avg_data2
                )
                self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2

                if len(r_cnv[key]) > 3:
                    penalty_value = 5
                    paired_changepoints = await run.cpu_bound(run_ruptures, r_cnv[key], penalty_value, self.cnv_dict["bin_width"])
                    approx_chrom_length=len(r_cnv[key]) * r_bin
                    if len(paired_changepoints) > 0:
                        for (start, end) in paired_changepoints:
                           if start < approx_chrom_length < end:
                               pass
                           elif start < 0 :
                               pass
                           elif start < self.centromere_bed[self.centromere_bed['chrom'].eq(key)]["end_pos"].max() and end > self.centromere_bed[self.centromere_bed['chrom'].eq(key)]["start_pos"].min():
                               pass
                           else:
                               item = np.array([(key, start, end)], dtype=self.dtype)
                               self.DATA_ARRAY = np.append(self.DATA_ARRAY, item)


        self.estimate_XY()

        np.save(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID), "ruptures.npy"
            ),
            self.DATA_ARRAY,
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
        # if not self.browse:
        #    for item in app.storage.general[self.mainuuid]:
        #        if item == 'sample_ids':
        #            for run in app.storage.general[self.mainuuid][item]:
        #                self.sampleID = run
        self.display_row = ui.row().style("width: 100")
        if self.summary:
            with self.summary:
                ui.label("No CNV data available.")
        with self.display_row:
            ui.label("Copy Number Variation").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
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
        if self.browse:
            ui.timer(0.1, lambda: self.show_previous_data(), once=True)
        else:
            ui.timer(15, lambda: self.show_previous_data())

    def generate_chart(
        self,
        title: Optional[str] = None,
        initmax: int = 8,
        initmin: int = 0,
        type: str = "value",
    ) -> ui.echart:
        """
        Generate an ECharts object for displaying CNV scatter plots.

        Args:
            title (Optional[str]): Title of the chart.
            initmax (int): Initial maximum value for the y-axis.
            initmin (int): Initial minimum value for the y-axis.
            type (str): Type of the x-axis.

        Returns:
            ui.echart: The ECharts object configured for the scatter plot.
        """
        return (
            ui.echart(
                {
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": f"{title}"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "tooltip": {"trigger": "axis"},
                    "xAxis": {
                        "type": f"{type}",
                        "max": "dataMax",
                        "splitLine": {"show": False},
                    },
                    "yAxis": [{
                        "type": "value",
                        "logBase": 2,
                        "name": "Ploidy",
                        "position": "left",
                    },
                    {
                        "type": "value",  # Secondary y-axis for the line plot
                        "name": "Candidate Break Points",
                        "position": "right",
                        "axisLine": {"lineStyle": {"color": "blue"}},  # Color the axis line to distinguish it
                        "splitLine": {"show": False},  # Optionally hide the grid lines for the secondary axis
                    },
                    ],
                    "dataZoom": [
                        {"type": "slider", "filterMode": "none"},
                        {
                            "type": "slider",
                            "yAxisIndex": 0,
                            "filterMode": "none",
                            "startValue": initmin,
                            "endValue": initmax,
                        },
                    ],
                    "series": [
                        {
                            "type": "scatter",
                            "symbolSize": 1,
                            "data": [],
                        }
                    ],
                }
            )
            .style("height: 450px")
            .classes("border-double")
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
            valueslist = {"All": "All"}  # Dictionary to store chromosome dropdown options
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
                    valueslist[counter] = contig  # Add chromosome to the dropdown options
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
                plot_to_update.options["series"] = []  # Clear any existing series in the plot

                # Loop through the sorted CNV data and plot each chromosome
                for contig, cnv in natsort.natsorted(result.cnv.items()):
                    # Skip non-standard chromosomes (e.g., mitochondrial DNA or unrecognized chromosomes)
                    if contig == "chrM" or not re.match(r'^chr(\d+|X|Y)$', contig):
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

                    data = reduce_list(data)  # Reduce the data points for more efficient plotting
                    total += len(cnv)  # Update total length with current chromosome data length

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
                plot_to_update.options["series"] = []  # Clear existing series in the plot

                # Update dropdown options for chromosomes
                for counter, contig in enumerate(
                        natsort.natsorted(result.cnv), start=1
                ):
                    valueslist[counter] = contig

                # Define the main chromosomes to be included in the plot
                main_chromosomes = [
                    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
                    "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                    "chr20", "chr21", "chr22", "chrX", "chrY"
                ]

                # Filter CNV data to include only the main chromosomes
                filtered_data = {k: v for k, v in result.cnv.items() if k in main_chromosomes}

                # Select the data for the filtered chromosome based on the user's selection
                contig, cnv = natsort.natsorted(filtered_data.items())[
                    int(self.chrom_filter) - 1
                    ]

                # Prepare the data for plotting
                data = list(
                    zip((np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"], cnv)
                )

                # If no specific gene target is selected, keep the default axis limits
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
                plot_to_update.options["xAxis"]["max"] = max  # Set x-axis maximum
                plot_to_update.options["xAxis"]["min"] = min  # Set x-axis minimum

                # Add the data for the selected chromosome as a scatter plot series
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
                            "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},  # Light blue color
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
                #print(self.DATA_ARRAY[contig])
                # Prepare the data for plotting
                breakpoints = self.DATA_ARRAY[self.DATA_ARRAY['name'] == contig]

                # Use defaultdict to accumulate increments and decrements
                increments = defaultdict(int)
                for entry in breakpoints:
                    increments[entry['start']] += 1
                    increments[entry['end'] + 1] -= 1

                # Generate the coordinates directly
                coordinates = []
                current_value = 0
                sorted_positions = sorted(increments.keys())
                for position in sorted_positions:
                    current_value += increments[position]
                    coordinates.append((position, current_value))

                if len(coordinates)>0:

                    print(f"Max value: {data[-1]}")
                    # Extract y-values
                    x_values, y_values = zip(*coordinates)

                    # Define the threshold
                    #threshold_value = np.percentile(y_values, 99)

                    threshold_value = 5

                    # Find where the y_values cross the threshold
                    above_threshold = np.array(y_values) > threshold_value

                    # Extract the x-values where y-values are above the threshold
                    x_above_threshold = np.array(x_values)[above_threshold]

                    # Determine the x-range where y crosses above the threshold
                    if x_above_threshold.size > 0:
                        x_range_min = np.min(x_above_threshold)
                        x_range_max = np.max(x_above_threshold)
                        x_range = x_range_max - x_range_min

                        print(f"Threshold Value: {threshold_value}")
                        print(f"X-Range where Y crosses above the threshold: {x_range}")
                        print(f"X-Range Min: {x_range_min}")
                        print(f"X-Range Max: {x_range_max}")
                    else:
                        print("No x-values are above the threshold.")

                plot_to_update.options["series"].append(
                    {
                        "type": "line",
                        "name": "Line Plot",
                        "data": coordinates,  # Your line plot data (e.g., [(x1, y1), (x2, y2), ...])
                        "yAxisIndex": 1,  # Specify to use the secondary y-axis (index 1)
                        "lineStyle": {"color": "blue"},  # Color the line to match the y-axis
                        "symbol": "none",  # Optionally remove symbols from the line plot
                    }
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
                return plot_to_update  # If not in UI mode, return the updated plot object


    async def show_previous_data(self) -> None:
        """
        Load and display previously computed CNV data.

        Args:
            output (str): The directory containing previous CNV data.
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
                if key != "chrM" and re.match(r'^chr(\d+|X|Y)$', key):
                    moving_avg_data1 = moving_average(self.result.cnv[key])
                    moving_avg_data2 = moving_average(r2_cnv[key])
                    moving_avg_data1, moving_avg_data2 = pad_arrays(
                        moving_avg_data1, moving_avg_data2
                    )
                    self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2

            if self.check_file_time(os.path.join(output, "ruptures.npy")):
                self.DATA_ARRAY = np.load(
                    os.path.join(output, "ruptures.npy"), allow_pickle="TRUE"
                )



            self.update_plots()
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    with ui.row():
                        with open(os.path.join(output, "XYestimate.pkl"), "rb") as file:
                            XYestimate = pickle.load(file)
                            if XYestimate != "Unknown":
                                if XYestimate == "XY":
                                    ui.icon("man").classes("text-4xl")
                                else:
                                    ui.icon("woman").classes("text-4xl")
                                ui.label(f"Estimated Genetic Sex: {XYestimate}")
                            ui.label(f"Current Bin Width: {self.cnv_dict['bin_width']}")
                            ui.label(
                                f"Current Variance: {round(self.cnv_dict['variance'], 3)}"
                            )


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
