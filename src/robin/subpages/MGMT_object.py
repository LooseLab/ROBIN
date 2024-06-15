"""
MGMT Analysis Module

This module provides functionality for analyzing MGMT methylation data. The primary
components include running external tools to extract and process MGMT sites from BAM files,
visualizing the results using a GUI, and managing the process asynchronously.

Dependencies:
    - pandas: Data manipulation and analysis.
    - os: Interaction with the operating system.
    - sys: System-specific parameters and functions.
    - asyncio: Asynchronous I/O.
    - nicegui: GUI creation.
    - pysam: BAM file processing.
    - shutil: File operations.
    - click: Command-line interface creation.
    - pathlib: File system paths.
    - natsort: Natural sorting.
    - tempfile: Temporary file creation.
    - logging: Logging for debugging and monitoring.

Modules:
    - subpages.base_analysis: BaseAnalysis class from robin.subpages.base_analysis.
    - theme: robin theme module.
    - submodules: robin submodules.

Environment Variables:
    - CI: Set to "1".

Constants:
    - HVPATH: Path to the 'hv_rapidCNS2' directory in the submodules.

Functions:
    - run_methylartist(tempmgmtdir: str, plot_out: str) -> None: Executes the methylartist tool to generate plots.
    - run_bedtools(bamfile: str, MGMT_BED: str, tempbamfile: str) -> None: Extracts MGMT sites from BAM files using bedtools.
    - run_modkit(tempmgmtdir: str, MGMTbamfile: str, threads: int) -> None: Processes BAM files with modkit and runs an R script for MGMT prediction.

Classes:
    - MGMT_Object(BaseAnalysis): Manages the MGMT analysis process, including setting up the GUI and handling BAM file processing.

Command-line Interface:
    - main(port: int, threads: int, watchfolder: str, output: str, browse: bool) -> None: CLI entry point for running the app, using Click for argument parsing.

Usage:
    The module can be run as a script to start the GUI for MGMT analysis, specifying
    options like the port, number of threads, watch folder, and output directory.
"""

from robin.subpages.base_analysis import BaseAnalysis
from robin import theme
from robin import submodules
import pandas as pd
import os
import sys
import asyncio
from nicegui import ui, run
import pysam
import shutil
import click
from pathlib import Path
import natsort
import tempfile
import logging
from typing import Optional, Tuple, List

# Configure logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
#logger = logging.getLogger(__name__)

os.environ["CI"] = "1"

HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)


def run_methylartist(tempmgmtdir: str, plot_out: str) -> None:
    """
    Executes the methylartist tool to generate plots for the given BAM file.

    Args:
        tempmgmtdir (str): Temporary directory containing the BAM file.
        plot_out (str): Output path for the plot.

    Returns:
        None
    """
    #logger.debug(
    #    f"Running methylartist with tempmgmtdir={tempmgmtdir}, plot_out={plot_out}"
    #)
    try:
        os.system(
            f"methylartist locus -i chr10:129466536-129467536 -b {os.path.join(tempmgmtdir, 'mgmt.bam')} -o {plot_out} --motif CG --mods m > /dev/null 2>&1"
        )
    except Exception as e:
        #logger.error(f"Error running methylartist: {e}")
        raise


def run_bedtools(bamfile: str, MGMT_BED: str, tempbamfile: str) -> None:
    """
    Extracts the MGMT sites from the BAM file using bedtools.

    Args:
        bamfile (str): Path to the input BAM file.
        MGMT_BED (str): Path to the MGMT BED file.
        tempbamfile (str): Path to the output temporary BAM file.

    Returns:
        None
    """
    #logger.debug(
    #    f"Running bedtools with bamfile={bamfile}, MGMT_BED={MGMT_BED}, tempbamfile={tempbamfile}"
    #)
    try:
        os.system(f"bedtools intersect -a {bamfile} -b {MGMT_BED} > {tempbamfile}")
        pysam.index(tempbamfile, f"{tempbamfile}.bai")
    except Exception as e:
        #logger.error(f"Error running bedtools: {e}")
        raise


def run_modkit(tempmgmtdir: str, MGMTbamfile: str, threads: int) -> None:
    """
    Processes the BAM file with modkit and runs an R script for MGMT prediction.

    Args:
        tempmgmtdir (str): Temporary directory for processing.
        MGMTbamfile (str): Path to the MGMT BAM file.
        threads (int): Number of threads to use.

    Returns:
        None
    """
    #logger.debug(
    #    f"Running modkit with tempmgmtdir={tempmgmtdir}, MGMTbamfile={MGMTbamfile}, threads={threads}"
    #)
    try:
        pysam.sort("-o", os.path.join(tempmgmtdir, "mgmt.bam"), MGMTbamfile)
        pysam.index(
            os.path.join(tempmgmtdir, "mgmt.bam"), f"{tempmgmtdir}/mgmt.bam.bai"
        )
        cmd = f"modkit pileup -t {threads} --filter-threshold 0.73 --combine-mods {os.path.join(tempmgmtdir, 'mgmt.bam')} {os.path.join(tempmgmtdir, 'mgmt.bed')} --suppress-progress >/dev/null 2>&1"
        os.system(cmd)

        cmd = f"Rscript {HVPATH}/bin/mgmt_pred_v0.3.R --input={os.path.join(tempmgmtdir, 'mgmt.bed')} --out_dir={tempmgmtdir} --probes={HVPATH}/bin/mgmt_probes.Rdata --model={HVPATH}/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
        print(cmd)
        os.system(cmd)
    except Exception as e:
        #logger.error(f"Error running modkit: {e}")
        raise


class MGMT_Object(BaseAnalysis):
    """
    MGMT_Object handles the MGMT analysis process, including setting up the GUI
    and processing BAM files asynchronously.

    Attributes:
        MGMTbamfile (Optional[str]): Path to the MGMT BAM file.
        counter (int): Counter for the number of analyses.
        last_seen (int): Last seen timestamp for analysis.
    """

    def __init__(self, *args, **kwargs):
        self.MGMTbamfile: Optional[str] = None
        self.counter: int = 0
        self.last_seen: int = 0
        #logger.debug("Initializing MGMT_Object")
        super().__init__(*args, **kwargs)

    def setup_ui(self) -> None:
        """
        Sets up the user interface for the MGMT analysis.

        Returns:
            None
        """
        #logger.debug("Setting up UI")
        with ui.card().style("width: 100%"):
            ui.label("MGMT Methylation").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            self.mgmtable = ui.row().classes("w-full")
            with self.mgmtable:
                ui.label("Table not yet available.")
            self.mgmtplot = ui.row().style("width: 100%")
            with self.mgmtplot:
                ui.label("Plot not yet available.")
        if self.summary:
            with self.summary:
                ui.label("Current MGMT status: Unknown")
        if self.browse:
            self.show_previous_data(self.output)
        else:
            ui.timer(30, lambda: self.show_previous_data(self.output))

    async def process_bam(self, bamfile: str, timestamp: str) -> None:
        """
        Processes the BAM file to extract and analyze MGMT sites.

        Args:
            bamfile (str): Path to the input BAM file.
            timestamp (str): Timestamp for the analysis.

        Returns:
            None
        """
        #logger.debug(f"Processing BAM file: {bamfile} at {timestamp}")
        MGMT_BED: str = f"{HVPATH}/bin/mgmt_hg38.bed"
        tempbamfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")

        try:
            await run.cpu_bound(run_bedtools, bamfile, MGMT_BED, tempbamfile.name)
        except Exception as e:
            #logger.error(f"Error in process_bam: {e}")
            return

        try:
            if pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True) > 0:
                if not self.MGMTbamfile:
                    self.MGMTbamfile = os.path.join(self.output, "mgmt.bam")
                    shutil.copy2(tempbamfile.name, self.MGMTbamfile)
                    os.remove(f"{tempbamfile.name}.bai")
                else:
                    tempbamholder = tempfile.NamedTemporaryFile(
                        dir=self.output, suffix=".bam"
                    )
                    pysam.cat(
                        "-o", tempbamholder.name, self.MGMTbamfile, tempbamfile.name
                    )
                    shutil.copy2(tempbamholder.name, self.MGMTbamfile)
                    try:
                        os.remove(f"{tempbamholder.name}.bai")
                        os.remove(f"{tempbamfile.name}.bai")
                    except FileNotFoundError:
                        pass
                tempmgmtdir = tempfile.TemporaryDirectory(dir=self.output)

                await run.cpu_bound(
                    run_modkit, tempmgmtdir.name, self.MGMTbamfile, self.threads
                )

                try:
                    results = pd.read_csv(
                        os.path.join(tempmgmtdir.name, "live_analysis_mgmt_status.csv")
                    )
                    self.counter += 1
                    plot_out = os.path.join(self.output, f"{self.counter}_mgmt.png")

                    await run.cpu_bound(run_methylartist, tempmgmtdir.name, plot_out)
                    results.to_csv(
                        os.path.join(self.output, f"{self.counter}_mgmt.csv"),
                        index=False,
                    )
                except Exception as e:
                    #logger.error(f"Error processing results: {e}")
                    raise
            else:
                os.remove(f"{tempbamfile.name}.bai")
        except Exception as e:
            #logger.error(f"Error in BAM file processing: {e}")
            raise
        finally:
            await asyncio.sleep(0.1)
            self.running = False

    def tabulate(self, results: pd.DataFrame) -> None:
        """
        Displays the results in a tabular format.

        Args:
            results (pd.DataFrame): DataFrame containing the analysis results.

        Returns:
            None
        """
        #logger.debug("Tabulating results")
        ui.aggrid.from_pandas(
            results,
            theme="material",
            options={
                "defaultColDef": {
                    "sortable": True,
                    "resizable": True,
                },
                "columnDefs": [
                    {
                        "headerName": "Average",
                        "field": "average",
                        "filter": "agTextColumnFilter",
                        "floatingFilter": False,
                    },
                    {
                        "headerName": "Score",
                        "field": "pred",
                        "filter": "agNumberColumnFilter",
                        "floatingFilter": False,
                    },
                    {
                        "headerName": "Status",
                        "field": "status",
                        "filter": "agNumberColumnFilter",
                        "floatingFilter": False,
                    },
                ],
                "pagination": True,
            },
        ).classes("w-full").style("height: 200px")

    def get_report(self, watchfolder: str) -> Tuple[pd.DataFrame, str, str]:
        """
        Generates a report from the analysis results.

        Args:
            watchfolder (str): Path to the folder containing previous analysis results.

        Returns:
            Tuple[pd.DataFrame, str, str]: DataFrame with results, path to plot, and summary string.
        """
        #logger.debug(f"Generating report from {watchfolder}")
        results = pd.DataFrame()
        plot_out = ""
        summary = ""
        for file in natsort.natsorted(os.listdir(watchfolder)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split('_')[0])
                if count > self.last_seen:
                    results = pd.read_csv(os.path.join(watchfolder, file))
                    plot_out = os.path.join(watchfolder, file.replace(".csv", ".png"))
                    summary = f"Current MGMT status: {results['status'].values[0]}"
        return results, plot_out, summary

    def show_previous_data(self, watchfolder: str) -> None:
        """
        Displays previously analyzed data from the specified watch folder.

        Args:
            watchfolder (str): Path to the folder containing previous analysis results.

        Returns:
            None
        """
        #logger.debug(f"Showing previous data from {watchfolder}")
        if not self.last_seen:
            for file in natsort.natsorted(os.listdir(watchfolder)):
                if file.endswith("_mgmt.csv"):
                    count = int(file.split('_')[0])

                    if count > self.last_seen:
                        results = pd.read_csv(os.path.join(watchfolder, file))
                        plot_out = os.path.join(
                            watchfolder, file.replace(".csv", ".png")
                        )
                        self.mgmtable.clear()
                        with self.mgmtable:
                            self.tabulate(results)
                        if os.path.exists(plot_out):
                            self.mgmtplot.clear()
                            with self.mgmtplot.classes("w-full"):
                                ui.image(plot_out).props("fit=scale-up")
                        if self.summary:
                            with self.summary:
                                self.summary.clear()
                                ui.label(
                                    f"Current MGMT status: {results['status'].values[0]}"
                                )
                        self.last_seen = count


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
) -> None:
    """
    Sets up and runs the MGMT analysis application.

    Args:
        port (int): Port number for the server.
        threads (int): Number of threads to use for processing.
        watchfolder (str): Path to the folder to watch for new BAM files.
        output (str): Path to the output directory.
        reload (bool): Flag to reload the application on changes.
        browse (bool): Flag to enable browsing historic data.

    Returns:
        None
    """
    #logger.debug(
    #    f"Starting MGMT analysis application on port {port} with {threads} threads"
    #)
    my_connection = None
    with theme.frame("MGMT Data", my_connection):
        TestObject = MGMT_Object(threads, output, progress=True)
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
def main(port: int, threads: int, watchfolder: Optional[str], output: Optional[str], browse: bool) -> None:
    """
    CLI entry point for running the MGMT analysis app.

    Args:
        port (int): The port to serve the app on.
        threads (int): Number of threads available for processing.
        watchfolder (Optional[str]): Directory to watch for new BAM files.
        output (Optional[str]): Directory to save output files.
        browse (bool): Enable browsing historic data.

    Returns:
        None
    """
    #logger.debug(
    #    f"Running main function with port={port}, threads={threads}, watchfolder={watchfolder}, output={output}, browse={browse}"
    #)
    if browse:
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=None,
            output=output,
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
    main()
