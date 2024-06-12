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

- `cnsmeth.subpages.base_analysis.BaseAnalysis`

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

if name in {"main", "mp_main"}:
    main()

"""

from cnv_from_bam import iterate_bam_file
from cnsmeth.subpages.base_analysis import BaseAnalysis
import natsort
from cnsmeth import theme, resources
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

os.environ["CI"] = "1"


class Result:
    def __init__(self, cnv_dict):
        self.cnv = cnv_dict


def iterate_bam(bamfile, _threads, mapq_filter, copy_numbers, log_level):
    result = iterate_bam_file(
        bamfile,
        _threads=_threads,
        mapq_filter=mapq_filter,
        copy_numbers=copy_numbers,
        log_level=log_level,
    )
    return result.cnv, result.bin_width, result.variance, copy_numbers


def iterate_bam_bin(bamfile, _threads, mapq_filter, copy_numbers, log_level, bin_width):
    result = iterate_bam_file(
        bamfile,
        _threads=_threads,
        mapq_filter=mapq_filter,
        copy_numbers=copy_numbers,
        log_level=log_level,
        bin_width=bin_width,
    )
    return result.cnv, result.bin_width, result.variance, copy_numbers


def reduce_list(lst, max_length=1000):
    while len(lst) > max_length:
        lst = lst[::2]
    return lst


class CNV_Difference:
    def __init__(self, *args, **kwargs):
        self.cnv = {}


def moving_average(data, n=3):
    return np.convolve(data, np.ones(n) / n, mode="same")


def pad_arrays(arr1, arr2, pad_value=0):
    len_diff = abs(len(arr1) - len(arr2))
    if len(arr1) < len(arr2):
        arr1 = np.pad(arr1, (0, len_diff), mode="constant", constant_values=pad_value)
    elif len(arr1) > len(arr2):
        arr2 = np.pad(arr2, (0, len_diff), mode="constant", constant_values=pad_value)
    return arr1, arr2


def ruptures_plotting(data):
    x_coords = range(0, (data.size + 1))  # * 10, bin_slice * 10, )
    signal = np.array(list(zip(x_coords, data)))

    algo_c = rpt.KernelCPD(kernel="linear", min_size=10).fit(signal)  # [:, 1]

    return algo_c


class CNVAnalysis(BaseAnalysis):
    def __init__(self, *args, target_panel=None, **kwargs):
        self.file_list = []
        self.cnv_dict = {}
        self.cnv_dict["bin_width"] = 0
        self.cnv_dict["variance"] = 0
        self.update_cnv_dict = {}
        self.result = None
        self.result2 = None
        self.result3 = CNV_Difference()

        self.XYestimate = "Unknown"

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
            os.path.join(self.gene_bed_file),
            names=["chrom", "start_pos", "end_pos", "gene"],
            header=None,
            sep="\s+",
        )
        super().__init__(*args, **kwargs)

    def estimate_XY(self):
        # We remove zero points as they are likely centromeric.
        X = round(np.average([i for i in self.result3.cnv["chrX"] if i != 0]), 2)
        Y = round(np.average([i for i in self.result3.cnv["chrY"] if i != 0]), 2)
        # print (f"X={X},Y={Y}")
        if X >= 0.1 and Y <= 0.1:
            self.XYestimate = "XX"
        elif X <= 0.1 and Y >= -0.2:
            self.XYestimate = "XY"
        else:
            self.XYestimate = "Unknown"
        file = open(os.path.join(self.output, "XYestimate.pkl"), "wb")
        pickle.dump(self.XYestimate, file)

    async def process_bam(self, bamfile, timestamp):
        self.file_list.append(bamfile)
        await self.do_cnv_work(bamfile)

    async def do_cnv_work(self, bamfile):
        r_cnv, r_bin, r_var, self.update_cnv_dict = await run.cpu_bound(
            iterate_bam,
            bamfile,
            self.threads,
            60,
            self.update_cnv_dict,
            int(logging.ERROR),
        )

        self.cnv_dict["bin_width"] = r_bin  # self.result.bin_width
        self.cnv_dict["variance"] = r_var  # self.result.variance
        r2_cnv, r2_bin, r2_var, self.ref_cnv_dict = await run.cpu_bound(
            iterate_bam_bin,
            bamfile,
            self.threads,
            60,
            self.ref_cnv_dict,
            int(logging.ERROR),
            bin_width=self.cnv_dict["bin_width"],
        )
        # self.result2 = iterate_bam_file(
        #    bam_file_path=None,
        #    _threads=self.threads,
        #    mapq_filter=60,
        #    copy_numbers=self.ref_cnv_dict,
        #    log_level=int(logging.ERROR),
        #    bin_width=self.cnv_dict["bin_width"],
        # )
        # await asyncio.sleep(0)

        for key in r_cnv.keys():
            if key != "chrM":
                moving_avg_data1 = await run.cpu_bound(moving_average, r_cnv[key])
                moving_avg_data2 = await run.cpu_bound(moving_average, r2_cnv[key])
                moving_avg_data1, moving_avg_data2 = await run.cpu_bound(
                    pad_arrays, moving_avg_data1, moving_avg_data2
                )
                self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2
                # print(key, np.mean(self.result3.cnv[key]), np.mean([i for i in self.result3.cnv[key] if i !=0]))
                # if len(self.result3.cnv[key]) > 20:
                #    algo_c = ruptures_plotting(self.result3.cnv[key])
                #    penalty_value = 5  # beta
                #    result = algo_c.predict(pen=penalty_value)
                # print(key, result)

        # await asyncio.sleep(0)

        self.estimate_XY()

        np.save(os.path.join(self.output, "CNV.npy"), r_cnv)
        np.save(os.path.join(self.output, "CNV_dict.npy"), self.cnv_dict)

        self.running = False

    def update_plots(self, gene_target=None):
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
            # ui.update(self.scatter_echart)

    def setup_ui(self):
        self.display_row = ui.row()
        if self.summary:
            with self.summary:
                ui.label("No CNV data available.")
        with self.display_row:
            # self.progrock.visible = False
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
                self.cnv_dict, "variance", backward=lambda n: f"Variance: {round(n,3)}"
            )

        self.scatter_echart = self.generate_chart(title="CNV Scatter Plot")

        self.difference_scatter_echart = self.generate_chart(
            title="Difference Plot", initmin=-2, initmax=2
        )  # , type="log")

        with ui.expansion("See Reference DataSet", icon="loupe").classes("w-full"):
            self.reference_scatter_echart = self.generate_chart(
                title="Reference CNV Scatter Plot"
            )
        if self.browse:
            ui.timer(0.1, lambda: self.show_previous_data(self.output), once=True)
        else:
            ui.timer(15, lambda: self.show_previous_data(self.output))

    def generate_chart(self, title=None, initmax=8, initmin=0, type="value"):
        return (
            ui.echart(
                {
                    "textStyle": {
                        "fontFamily": "Fira Sans, Fira Mono"
                    },
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": f"{title}"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {
                        "type": f"{type}",
                        "max": "dataMax",
                        "splitLine": {"show": False},
                    },
                    "yAxis": {
                        "type": "value",
                        "logBase": 2,
                    },
                    "dataZoom": [
                        {"type": "slider", "filterMode": "none"},
                        {
                            "type": "slider",
                            "yAxisIndex": 0,
                            "filterMode": "none",
                            "startValue": initmin,
                            "endValue": initmax,
                        },
                        # {"type": "inside", "xAxisIndex": 0, "filterMode": "none"},
                        # {"type": "inside", "yAxisIndex": 0, "filterMode": "none"},
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
        plot_to_update=None,
        result=None,
        gene_target=None,
        title=None,
        min=0,
        ui_mode=True,
    ):
        if result or self.result:
            total = 0
            valueslist = {"All": "All"}
            genevalueslist = {"All": "All"}
            try:
                self.chrom_filter = self.chrom_select.value
            except AttributeError:
                self.chrom_filter = "All"

            min = min
            max = "dataMax"

            if gene_target:
                start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                end_pos = self.gene_bed.iloc[int(gene_target)].end_pos
                chrom = self.gene_bed.iloc[int(gene_target)].chrom
                for counter, contig in enumerate(
                    natsort.natsorted(result.cnv), start=1
                ):
                    valueslist[counter] = contig
                    if contig == chrom:
                        break
                self.chrom_filter = counter
                min = start_pos - 10 * self.cnv_dict["bin_width"]
                max = end_pos + 10 * self.cnv_dict["bin_width"]
            if self.chrom_filter == "All":
                counter = 0

                plot_to_update.options["title"]["text"] = f"{title} - All Chromosomes"
                plot_to_update.options["series"] = []
                for contig, cnv in natsort.natsorted(result.cnv.items()):
                    if contig == "chrM":
                        continue
                    counter += 1
                    valueslist[counter] = contig

                    data = list(
                        zip(
                            (np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"],
                            cnv,
                        )
                    )

                    data = reduce_list(data)

                    total += len(cnv)

                    plot_to_update.options["xAxis"]["max"] = max
                    plot_to_update.options["xAxis"]["min"] = min
                    plot_to_update.options["series"].append(
                        {
                            "type": "scatter",
                            "name": contig,
                            "data": data,
                            "symbolSize": 5,
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
                    for index, gene in self.gene_bed[
                        self.gene_bed["chrom"] == contig
                    ].iterrows():
                        genevalueslist[index] = f"{gene.chrom} - {gene.gene}"
            else:
                plot_to_update.options["series"] = []

                for counter, contig in enumerate(
                    natsort.natsorted(result.cnv), start=1
                ):
                    valueslist[counter] = contig

                contig, cnv = natsort.natsorted(result.cnv.items())[
                    int(self.chrom_filter) - 1
                ]

                data = list(
                    zip((np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"], cnv)
                )

                if not gene_target:
                    min = min
                    max = "dataMax"

                else:
                    if gene_target == "All":
                        min = min
                        max = "dataMax"
                    else:
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
                        # self.update_plot(gene_target=[start_pos, end_pos])

                        min = start_pos - 10 * self.cnv_dict["bin_width"]
                        max = end_pos + 10 * self.cnv_dict["bin_width"]

                        if start_pos - min > 2_000_000:
                            min = start_pos - 2_000_000
                        if max - end_pos > 2_000_000:
                            max = end_pos + 2_000_000

                        if min < 0:
                            min = 0

                plot_to_update.options["title"][
                    "text"
                ] = f"Copy Number Variation - {contig}"
                plot_to_update.options["xAxis"]["max"] = max
                plot_to_update.options["xAxis"]["min"] = min
                plot_to_update.options["series"].append(
                    {
                        "type": "scatter",
                        "name": contig,
                        "data": data,
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(255, 173, 177, 0.4)"},
                            "data": [],
                        },
                    }
                )
                for index, gene in self.gene_bed[
                    self.gene_bed["chrom"] == contig
                ].iterrows():
                    genevalueslist[index] = f"{gene.chrom} - {gene.gene}"

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
            try:
                self.chrom_select.set_options(valueslist)
            except AttributeError:
                pass
            try:
                self.gene_select.set_options(genevalueslist)
            except AttributeError:
                pass
            if ui_mode:
                ui.update(plot_to_update)
            else:
                return plot_to_update

    async def show_previous_data(self, output):
        if self.check_file_time(os.path.join(output, "CNV.npy")):
            result = np.load(os.path.join(output, "CNV.npy"), allow_pickle="TRUE").item()
            self.result = Result(result)
            cnv_dict = np.load(
                os.path.join(output, "CNV_dict.npy"), allow_pickle=True
            ).item()
            self.cnv_dict["bin_width"] = cnv_dict["bin_width"]
            self.cnv_dict["variance"] = cnv_dict["variance"]
            # self._update_cnv_plot()
            # self._update_cnv_plot(
            #    plot_to_update=self.scatter_echart, result=self.result, title="CNV"
            # )

            # ToDo: We are recalculating the CNV data here. We should load the data from the output folder.
            # If it does not exist already, we should calculate it and save it.
            # This will reduce the burden when multiple people view the site.

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

            # self.result2 = iterate_bam_file(
            #    bam_file_path=None,
            #    _threads=1,
            #    mapq_filter=60,
            #    copy_numbers=self.ref_cnv_dict,
            #    log_level=int(logging.ERROR),
            #    bin_width=self.cnv_dict["bin_width"],
            # )

            for key in self.result.cnv.keys():
                if key != "chrM":
                    # print(key, np.mean(self.result.cnv[key]))#[i for i in self.result.cnv[key] if i !=0]))
                    moving_avg_data1 = moving_average(self.result.cnv[key])
                    moving_avg_data2 = moving_average(r2_cnv[key])
                    self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2
                    # print(key, np.mean(self.result3.cnv[key]), np.mean([i for i in self.result3.cnv[key] if i !=0]))
                    if len(self.result3.cnv[key]) > 20:
                        algo_c = ruptures_plotting(self.result3.cnv[key])
                        penalty_value = 5  # beta

                        result = algo_c.predict(pen=penalty_value)
                        # print(key, result)

            # self.estimate_XY()
            self.update_plots()
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    with ui.row():
                        file = open(os.path.join(output, "XYestimate.pkl"), "rb")
                        XYestimate = pickle.load(file)
                        if XYestimate != "Unknown":
                            if XYestimate == "XY":
                                ui.icon("man").classes("text-4xl")
                            else:
                                ui.icon("woman").classes("text-4xl")
                            ui.label(f"Estimated Genetic Sex: {XYestimate}")
                        ui.label(f"Current Bin Width: {self.cnv_dict['bin_width']}")
                        ui.label(f"Current Variance: {round(self.cnv_dict['variance'], 3)}")


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
    target_panel: str = "rCNS2",
):
    # my_connection = None
    app.add_static_files("/fonts", str(Path(__file__).parent / "../fonts"))
    with theme.frame("Copy Number Variation Testing."):
        TestObject = CNVAnalysis(
            threads,
            output,
            progress=True,
            # bamqueue=self.bamforcnv,
            # summary=cnvsummary,
            target_panel=target_panel,
        )
        # TestObject = CNVAnalysis(threads, output, progress=True)
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
        # print("Browse mode not implemented.")
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
def main(port, threads, watchfolder, output, browse, target_panel):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
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
            target_panel=target_panel,
            # exclude=exclude,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()
