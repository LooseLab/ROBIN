"""
NanoDX Analysis Module

This module provides functionality for analyzing NanoDX methylation data. The primary
components include running external tools to extract and process methylation sites from BAM files,
visualizing the results using a GUI, and managing the process asynchronously.

Dependencies:
    - pandas: Data manipulation and analysis.
    - os: Interaction with the operating system.
    - sys: System-specific parameters and functions.
    - time: Time-related functions.
    - nicegui: GUI creation.
    - pysam: BAM file processing.
    - shutil: File operations.
    - tempfile: Temporary file creation.
    - click: Command-line interface creation.
    - pathlib: File system paths.
    - typing: Type hinting.
    - cnsmeth: Custom modules for the specific workflow.

Modules:
    - subpages.base_analysis: BaseAnalysis class from cnsmeth.subpages.base_analysis.
    - theme: cnsmeth theme module.
    - resources: cnsmeth resources module.
    - NN_model: NN_classifier class from cnsmeth.submodules.nanoDX.workflow.scripts.NN_model.
    - merge_bedmethyl: Functions for merging, saving, and collapsing bed methylation data.

Environment Variables:
    - CI: Set to "1".

Functions:
    - run_modkit: Executes modkit to extract methylation data from a BAM file.
    - run_samtools_sort: Sorts BAM files using Samtools.
    - classification: Runs classification on the extracted data using a neural network model.

Classes:
    - NanoDX_object: Manages the NanoDX analysis process, including setting up the GUI and handling BAM file processing.

Command-line Interface:
    - run_main: CLI entry point for running the app, using Click for argument parsing.

Usage:
    The module can be run as a script to start the GUI for NanoDX analysis, specifying
    options like the port, number of threads, watch folder, and output directory.
"""

from __future__ import annotations
from cnsmeth.subpages.base_analysis import BaseAnalysis
from nicegui import ui, app, run
import time
import os
import sys
import click
from pathlib import Path
import pysam
import pandas as pd
import shutil
import tempfile
from cnsmeth import models, theme, resources
from cnsmeth.submodules.nanoDX.workflow.scripts.NN_model import NN_classifier
from cnsmeth.utilities.merge_bedmethyl import (
    merge_bedmethyl,
    save_bedmethyl,
    collapse_bedmethyl,
)
from typing import List, Tuple, Optional, Dict, Any
from icecream import ic




def run_modkit(cpgs: str, sortfile: str, temp: str, threads: int) -> None:
    """
    Executes modkit on a bam file and extracts the methylation data.

    Args:
        cpgs (str): Path to the CpG BED file.
        sortfile (str): Path to the sorted BAM file.
        temp (str): Path to the temporary output file.
        threads (int): Number of threads to use.
    """
    try:
        #print (f"modkit pileup --include-bed {cpgs} --filter-threshold 0.73 --combine-mods --only-tabs -t {threads} {sortfile} {temp}")
        os.system(
                f"modkit pileup --include-bed {cpgs} --filter-threshold 0.73 --combine-mods --only-tabs -t {threads} {sortfile} {temp} --suppress-progress >/dev/null 2>&1"
        )
        #shutil.copy(f"{sortfile}","modkit.bam")
        #shutil.copy(f"{temp}", "modkitoutput.bed")
    except Exception as e:
        print(e)


def run_samtools_sort(file: str, tomerge: List[str], sortfile: str, threads: int) -> None:
    """
    Sorts BAM files using Samtools.

    Args:
        file (str): Path to the output BAM file.
        tomerge (List[str]): List of BAM files to merge.
        sortfile (str): Path to the sorted BAM file.
        threads (int): Number of threads to use.
    """
    pysam.cat("-o", file, *tomerge)
    pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)

modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "Capper_et_al_NN.pkl"
        )

NN = NN_classifier(modelfile)

def classification(modelfile: str, test_df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Runs classification on the extracted data using a neural network model.

    Args:
        modelfile (str): Path to the neural network model file.
        test_df (pd.DataFrame): DataFrame containing the test data.

    Returns:
        Tuple[np.ndarray, np.ndarray, int]: Predictions, class labels, and number of features.
    """
    #NN = NN_classifier(modelfile)
    try:
        predictions, class_labels, n_features = NN.predict(test_df)
    except Exception as e:
        ic(e)
        test_df.to_csv("errordf.csv", sep=",", index=False, encoding="utf-8")
        # sys.exit(1)
    return predictions, class_labels, n_features


class NanoDX_object(BaseAnalysis):
    """
    NanoDX_object handles the NanoDX analysis process, including setting up the GUI
    and processing BAM files asynchronously.

    Attributes:
        cpgs_file (str): Path to the CpG BED file.
        cpgs (pd.DataFrame): DataFrame containing CpG data.
        model (str): Name of the neural network model file.
        threshold (float): Threshold for classification confidence.
        nanodx_bam_count (int): Counter for the number of BAM files processed.
        not_first_run (bool): Flag indicating whether it's the first run.
        modelfile (str): Path to the neural network model file.
        nanodx_df_store (pd.DataFrame): DataFrame for storing NanoDX results.
        nanodxfile (tempfile.NamedTemporaryFile): Temporary file for NanoDX results.
    """

    def __init__(self, *args, model: str = "Capper_et_al_NN.pkl", **kwargs) -> None:
        self.cpgs_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "hglft_genome_260e9_91a970_clean.bed",
        )
        self.cpgs = pd.read_csv(
            self.cpgs_file,
            sep="\t",
            header=None,
        )
        self.model = model
        self.threshold = 0.05
        self.nanodx_bam_count = 0
        self.not_first_run = False
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), self.model
        )
        self.nanodx_df_store = pd.DataFrame()
        super().__init__(*args, **kwargs)
        self.nanodxfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".nanodx")

    def setup_ui(self) -> None:
        """
        Sets up the user interface for the NanoDX analysis.
        """
        with ui.card().style("width: 100%"):
            with ui.grid(columns=8).classes("w-full h-auto"):
                with ui.card().classes("col-span-3"):
                    self.create_nanodx_chart("NanoDX")
                with ui.card().classes("col-span-5"):
                    self.create_nanodx_time_chart("NanoDX Time Series")
        if self.summary:
            with self.summary:
                ui.label("NanoDX classification: Unknown")
        if self.browse:
            self.show_previous_data(self.output)
        else:
            ui.timer(5, lambda: self.show_previous_data(self.output))

    def show_previous_data(self, output: str) -> None:
        """
        Displays previously analyzed data from the specified output folder.

        Args:
            output (str): Path to the folder containing previous analysis results.
        """
        if self.check_file_time(os.path.join(output, "nanoDX_scores.csv")):
            self.nanodx_df_store = pd.read_csv(
                os.path.join(os.path.join(self.output, "nanoDX_scores.csv")),
                index_col=0,
            )
            columns_greater_than_threshold = (
                self.nanodx_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = self.nanodx_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()
            self.update_nanodx_time_chart(self.nanodx_df_store.drop(columns=result))
            lastrow = self.nanodx_df_store.iloc[-1].drop("number_probes")
            lastrow_plot = lastrow.sort_values(ascending=False).head(10)
            lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    ui.label(
                        f"NanoDX classification: {lastrow_plot_top.index[0]} - {lastrow_plot_top.values[0]:.2f}"
                    )
            self.update_nanodx_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values),
                "All",
                self.nanodx_df_store.iloc[-1]["number_probes"],
            )

    async def process_bam(self, bamfile: List[Tuple[str, float]]) -> None:
        """
        Processes the BAM files and performs the NanoDX analysis.

        Args:
            bamfile (List[Tuple[str, float]]): List of BAM files with their timestamps.
        """
        tomerge: List[str] = []
        while len(bamfile) > 0:
            self.running = True
            (file, filetime) = bamfile.pop()
            self.nanodx_bam_count += 1
            tomerge.append(file)

            if len(tomerge) > 5:
                break
        app.storage.general[self.mainuuid][self.name]["counters"][
            "bams_in_processing"
        ] += len(tomerge)

        if len(tomerge) > 0:
            tempbam = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")
            sorttempbam = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")
            file = tempbam.name

            temp = tempfile.NamedTemporaryFile(dir=self.output)

            sortfile = sorttempbam.name

            await run.cpu_bound(
                run_samtools_sort, file, tomerge, sortfile, self.threads
            )

            await run.cpu_bound(
                run_modkit, self.cpgs_file, sortfile, temp.name, self.threads
            )

            try:
                os.remove(f"{sortfile}.csi")
            except FileNotFoundError:
                pass

            if self.not_first_run:
                bed_a = pd.read_table(
                    f"{temp.name}",
                    names=[
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
                    ],
                    dtype={
                        "chrom": "category",
                        "start_pos": "int32",
                        "end_pos": "int32",
                        "mod": "category",
                        "score": "int16",
                        "strand": "category",
                        "start_pos2": "int32",
                        "end_pos2": "int32",
                        "colour": "category",
                        "Nvalid": "int16",
                        "fraction": "float16",
                        "Nmod": "int16",
                        "Ncanon": "int16",
                        "Nother": "int16",
                        "Ndel": "int16",
                        "Nfail": "int16",
                        "Ndiff": "int16",
                        "Nnocall": "int16",
                    },
                    header=None,
                    sep="\s+",
                )
                self.merged_bed_file = await run.cpu_bound(
                    merge_bedmethyl, bed_a, self.merged_bed_file
                )
                save_bedmethyl(self.merged_bed_file, self.nanodxfile.name)
            else:
                shutil.copy(f"{temp.name}", self.nanodxfile.name)
                self.merged_bed_file = pd.read_table(
                    self.nanodxfile.name,
                    names=[
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
                    ],
                    dtype={
                        "chrom": "category",
                        "start_pos": "int32",
                        "end_pos": "int32",
                        "mod": "category",
                        "score": "int16",
                        "strand": "category",
                        "start_pos2": "int32",
                        "end_pos2": "int32",
                        "colour": "category",
                        "Nvalid": "int16",
                        "fraction": "float16",
                        "Nmod": "int16",
                        "Ncanon": "int16",
                        "Nother": "int16",
                        "Ndel": "int16",
                        "Nfail": "int16",
                        "Ndiff": "int16",
                        "Nnocall": "int16",
                    },
                    header=None,
                    sep="\s+",
                )
                self.not_first_run = True

            self.merged_bed_file = await run.cpu_bound(
                collapse_bedmethyl, self.merged_bed_file
            )
            test_df = pd.merge(
                self.merged_bed_file,
                self.cpgs,
                left_on=["chrom", "start_pos"],
                right_on=[0, 1],
            )
            test_df.rename(
                columns={3: "probe_id", "fraction": "methylation_call"},
                inplace=True,
            )
            test_df.loc[test_df["methylation_call"] < 60, "methylation_call"] = -1
            test_df.loc[test_df["methylation_call"] >= 60, "methylation_call"] = 1
            predictions, class_labels, n_features = await run.cpu_bound(
                classification, self.modelfile, test_df
            )

            nanoDX_df = pd.DataFrame({"class": class_labels, "score": predictions})
            nanoDX_save = nanoDX_df.set_index("class").T
            nanoDX_save["number_probes"] = n_features
            nanoDX_save["timestamp"] = time.time() * 1000

            self.nanodx_df_store = pd.concat(
                [self.nanodx_df_store, nanoDX_save.set_index("timestamp")]
            )

            self.nanodx_df_store.to_csv(os.path.join(self.output, "nanoDX_scores.csv"))

            app.storage.general[self.mainuuid][self.name]["counters"][
                "bam_processed"
            ] += len(tomerge)
            app.storage.general[self.mainuuid][self.name]["counters"][
                "bams_in_processing"
            ] -= len(tomerge)
        self.running = False

    def create_nanodx_chart(self, title: str) -> None:
        """
        Creates the NanoDX chart.

        Args:
            title (str): Title of the chart.
        """
        self.nanodxchart = self.create_chart(title)

    def update_nanodx_plot(self, x: List[str], y: List[float], count: int, n_features: int) -> None:
        """
        Updates the NanoDX plot with new data.

        Args:
            x (List[str]): List of tumor types.
            y (List[float]): Confidence scores for each tumor type.
            count (int): Number of BAM files used to generate the plot.
            n_features (int): Number of features detected during data analysis.
        """
        self.nanodxchart.options["title"][
            "text"
        ] = f"NanoDX: processed {count} bams and found {int(n_features)} features"
        self.nanodxchart.options["yAxis"]["data"] = x
        self.nanodxchart.options["series"] = [
            {"type": "bar", "name": "NanoDX", "data": y}
        ]
        self.nanodxchart.update()

    def create_nanodx_time_chart(self, title: str) -> None:
        """
        Creates the NanoDX time series chart.

        Args:
            title (str): Title of the chart.
        """
        self.nanodx_time_chart = self.create_time_chart(title)

    def update_nanodx_time_chart(self, datadf: pd.DataFrame) -> None:
        """
        Updates the NanoDX time series chart with new data.

        Args:
            datadf (pd.DataFrame): DataFrame containing the data to plot.
        """
        self.nanodx_time_chart.options["series"] = []
        for series, data in datadf.to_dict().items():
            data_list = [[key, value] for key, value in data.items()]
            if series != "number_probes":
                self.nanodx_time_chart.options["series"].append(
                    {
                        "animation": False,
                        "type": "line",
                        "smooth": True,
                        "name": series,
                        "emphasis": {"focus": "series"},
                        "endLabel": {
                            "show": True,
                            "formatter": "{a}",
                            "distance": 20,
                        },
                        "lineStyle": {
                            "width": 2,
                        },
                        "data": data_list,
                    }
                )
        self.nanodx_time_chart.update()


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
) -> None:
    """
    Sets up and runs the NanoDX analysis application.

    Args:
        port (int): Port number for the server.
        threads (int): Number of threads to use for processing.
        watchfolder (str): Path to the folder to watch for new BAM files.
        output (str): Path to the output directory.
        reload (bool): Flag to reload the application on changes.
        browse (bool): Flag to enable browsing historic data.
    """
    my_connection: Optional[Dict[str, Any]] = None
    with theme.frame("Target Coverage Data", my_connection):
        TestObject = NanoDX_object(threads, output, progress=True, batch=True)
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
def run_main(port: int, threads: int, watchfolder: str, output: str, browse: bool) -> None:
    """
    CLI entry point for running the NanoDX analysis app.

    Args:
        port (int): The port to serve the app on.
        threads (int): Number of threads available for processing.
        watchfolder (str): Directory to watch for new BAM files.
        output (str): Directory to save output files.
        browse (bool): Enable browsing historic data.
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
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    run_main()