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

from robin.subpages.base_analysis import BaseAnalysis, BaseVis
from robin import theme
from robin import submodules
import pandas as pd
import os
import sys
from nicegui import ui, run, app, background_tasks
import pysam
import shutil
import click
from pathlib import Path
import natsort
import tempfile
import logging
from typing import Optional, Tuple
from robin.core.state import state, ProcessState

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)

os.environ["CI"] = "1"

HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)


def has_reads(bam_file, chrom, start, end):
    """
    Quickly checks if any reads span a specific genomic locus.

    Args:
        bam_file (str): Path to your BAM file.
        chrom (str): Chromosome name (e.g., 'chr1').
        start (int): Start coordinate (0-based, inclusive).
        end (int): End coordinate (0-based, exclusive).

    Returns:
        bool: True if at least one read spans the region, False otherwise.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for _ in bam.fetch(chrom, start, end):
            return True  # Found at least one read
    return False  # No reads found


def run_methylartist(tempmgmtdir: str, plot_out: str) -> None:
    """
    Executes the methylartist tool to generate plots for the given BAM file.

    Args:
        tempmgmtdir (str): Temporary directory containing the BAM file.
        plot_out (str): Output path for the plot.

    Returns:
        None
    """
    # logger.debug(
    #    f"Running methylartist with tempmgmtdir={tempmgmtdir}, plot_out={plot_out}"
    # )
    try:
        os.system(
            f"methylartist locus -i chr10:129466536-129467536 -b {os.path.join(tempmgmtdir, 'mgmt.bam')} -o {plot_out} --motif CG --mods m > /dev/null 2>&1"
        )
    except Exception as e:
        logger.error(f"Error running methylartist: {e}")
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
    # logger.debug(
    #    f"Running bedtools with bamfile={bamfile}, MGMT_BED={MGMT_BED}, tempbamfile={tempbamfile}"
    # )
    try:
        os.system(f"bedtools intersect -a {bamfile} -b {MGMT_BED} > {tempbamfile}")
        pysam.index(tempbamfile, f"{tempbamfile}.bai")
    except Exception:
        # logger.error(f"Error running bedtools: {e}")
        raise


def run_modkit(
    tempmgmtdir: str, MGMTbamfile: str, threads: int, output_mgmt_bed: str
) -> None:
    """
    Processes the BAM file with modkit and runs an R script for MGMT prediction.

    Args:
        tempmgmtdir (str): Temporary directory for processing.
        MGMTbamfile (str): Path to the MGMT BAM file.
        threads (int): Number of threads to use.

    Returns:
        None
    """
    try:
        pysam.sort("-o", os.path.join(tempmgmtdir, "mgmt.bam"), MGMTbamfile)
        pysam.index(
            os.path.join(tempmgmtdir, "mgmt.bam"), f"{tempmgmtdir}/mgmt.bam.bai"
        )
        # cmd = f"modkit pileup -t {threads} --filter-threshold 0.73 --combine-mods --mixed-delim {os.path.join(tempmgmtdir, 'mgmt.bam')} {os.path.join(tempmgmtdir, 'mgmt.bed')} --suppress-progress >/dev/null 2>&1"
        # hard code to a single thread.
        cmd = f"modkit pileup -t 1 --filter-threshold 0.73 --combine-mods --mixed-delim --interval-size 20000000 --chunk-size 1 {os.path.join(tempmgmtdir, 'mgmt.bam')} {output_mgmt_bed} --suppress-progress >/dev/null 2>&1"
        os.system(cmd)
        if os.path.exists(output_mgmt_bed):
            cmd = f"Rscript {HVPATH}/bin/mgmt_pred_v0.3.R --input={output_mgmt_bed} --out_dir={tempmgmtdir} --probes={HVPATH}/bin/mgmt_probes.Rdata --model={HVPATH}/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
            os.system(cmd)
    except Exception:
        raise


class MGMTVis(BaseVis):
    """
    MGMT_Object handles the MGMT analysis process, including setting up the GUI
    and processing BAM files asynchronously.

    Attributes:
        MGMTbamfile (Optional[str]): Path to the MGMT BAM file.
        counter (int): Counter for the number of analyses.
        last_seen (int): Last seen timestamp for analysis.
    """

    def __init__(self, *args, **kwargs):
        self.MGMTbamfile = {}
        self.counter: int = 0
        self.last_seen: int = 0
        # logger.debug("Initializing MGMT_Object")
        super().__init__(*args, **kwargs)
        # Remove state tracking for MGMT Analysis
        # state.start_process("MGMT Analysis", ProcessType.BATCH)
        # state.set_process_state("MGMT Analysis", ProcessState.WAITING_FOR_DATA)

    def setup_ui(self) -> None:
        """
        Sets up the user interface for the MGMT analysis.

        Returns:
            None
        """
        # logger.debug("Setting up UI")
        with ui.card().style("width: 100%"):
            ui.label("MGMT Methylation").classes("text-sky-600 dark:text-white").style(
                "font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            self.mgmtable = ui.row().classes("w-full")
            with self.mgmtable:
                ui.label("Table not yet available.")
            self.mgmtplot = ui.row().style("width: 100%")
            with self.mgmtplot:
                ui.label("Plot not yet available.")
        if self.summary:
            with self.summary:
                with ui.card().classes("w-full p-4 mb-4"):
                    with ui.row().classes("w-full items-center justify-between"):
                        # Left side - Methylation Status
                        with ui.column().classes("gap-2"):
                            ui.label("MGMT Methylation Analysis").classes(
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
                            ui.label("Methylation Score: --").classes("text-gray-600")
                            ui.label("Average Methylation: --").classes("text-gray-600")

                    # Bottom row - Information
                    with ui.row().classes(
                        "w-full mt-4 text-sm text-gray-500 justify-center"
                    ):
                        ui.label("Methylation status based on MGMT promoter analysis")
        if self.browse:
            self.show_previous_data()
        else:
            ui.timer(30, lambda: self.show_previous_data())

    def tabulate(self, results: pd.DataFrame) -> None:
        """
        Displays the results in a tabular format.

        Args:
            results (pd.DataFrame): DataFrame containing the analysis results.

        Returns:
            None
        """
        # logger.debug("Tabulating results")
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

    def display_bed_data(self, bed_file: str) -> None:
        """
        Displays the modkit pileup bed file data in a table format.

        Args:
            bed_file (str): Path to the bed file.

        Returns:
            None
        """
        try:
            # Read bed file with pandas - adjust for actual file format
            bed_data = pd.read_csv(bed_file, sep="\t", header=None)

            # Check if the file has the expected number of columns
            if bed_data.shape[1] >= 10:
                # Assign column names based on the observed format
                column_names = [
                    "Chromosome",
                    "Start",
                    "End",
                    "Name",
                    "Score",
                    "Strand",
                    "Start2",
                    "End2",
                    "RGB",
                    "Coverage",
                    "Modified_Fraction",
                ]

                # Ensure we only use as many column names as there are columns
                bed_data.columns = column_names[: bed_data.shape[1]]

                # Add any missing columns with default values
                if "Coverage" not in bed_data.columns:
                    bed_data["Coverage"] = 1
                if "Modified_Fraction" not in bed_data.columns:
                    bed_data["Modified_Fraction"] = 0.0
            else:
                # If the format is completely different, log an error and create a minimal dataframe
                logger.error(
                    f"Unexpected bed file format with {bed_data.shape[1]} columns"
                )
                ui.label(
                    f"Unexpected bed file format with {bed_data.shape[1]} columns"
                ).classes("text-red-500")
                return

            # Calculate percentage methylation
            bed_data["Methylation_Percentage"] = (
                bed_data["Modified_Fraction"] * 100
            ).round(2)

            # Display the data in an AG Grid table
            ui.aggrid.from_pandas(
                bed_data,
                theme="material",
                options={
                    "defaultColDef": {
                        "sortable": True,
                        "resizable": True,
                        "filter": True,
                    },
                    "columnDefs": [
                        {"headerName": "Chr", "field": "Chromosome", "width": 80},
                        {"headerName": "Start", "field": "Start", "width": 100},
                        {"headerName": "End", "field": "End", "width": 100},
                        {"headerName": "Coverage", "field": "Coverage", "width": 100},
                        {
                            "headerName": "Modified %",
                            "field": "Methylation_Percentage",
                            "width": 120,
                        },
                    ],
                    "pagination": True,
                    "paginationPageSize": 10,
                },
            ).classes("w-full").style("height: 400px")
        except Exception as e:
            logger.error(f"Error displaying bed file data: {e}")
            ui.label(f"Error reading bed file: {str(e)}").classes("text-red-500")

    def display_specific_sites(self, bed_file: str) -> None:
        """
        Displays the specific methylation sites of interest in a table format.

        Args:
            bed_file (str): Path to the bed file.

        Returns:
            None
        """
        try:
            specific_sites = self.extract_specific_sites(bed_file)

            if specific_sites.empty:
                ui.label("No specific CpG pairs found in the data").classes(
                    "text-amber-500 mt-2"
                )
                ui.label(
                    "Looking for CpG pairs: 129467255/6, 129467258/9, 129467262/3, 129467272/3"
                ).classes("text-gray-600 text-sm mt-1")
                return

            # Save the specific sites data to CSV
            csv_file = self.save_specific_sites_data(specific_sites, bed_file)
            if csv_file:
                ui.label(
                    f"Specific sites data saved to: {os.path.basename(csv_file)}"
                ).classes("text-green-600 text-sm mt-1")

            # Display the specific sites in an AG Grid table
            ui.label("Key CpG Sites for MGMT Analysis").classes(
                "text-lg font-medium mt-4 mb-2"
            )

            ui.aggrid.from_pandas(
                specific_sites,
                theme="material",
                options={
                    "defaultColDef": {
                        "sortable": True,
                        "resizable": True,
                        "filter": True,
                    },
                    "columnDefs": [
                        {"headerName": "Site", "field": "Site_Label", "width": 200},
                        {"headerName": "Chr", "field": "Chromosome", "width": 80},
                        {
                            "headerName": "CpG Position",
                            "field": "Position",
                            "width": 120,
                        },
                        {
                            "headerName": "Forward Coverage",
                            "field": "Coverage_Forward",
                            "width": 140,
                        },
                        {
                            "headerName": "Reverse Coverage",
                            "field": "Coverage_Reverse",
                            "width": 140,
                        },
                        {
                            "headerName": "Total Coverage",
                            "field": "Total_Coverage",
                            "width": 130,
                        },
                        {
                            "headerName": "% Methylation",
                            "field": "Methylation_Percentage",
                            "width": 120,
                        },
                        {
                            "headerName": "Forward % Meth",
                            "field": "Forward_Methylation",
                            "width": 130,
                        },
                        {
                            "headerName": "Reverse % Meth",
                            "field": "Reverse_Methylation",
                            "width": 130,
                        },
                        {"headerName": "Notes", "field": "Notes", "width": 250},
                    ],
                    "pagination": False,
                },
            ).classes(
                "w-full"
            )  # .style("height: 200px")

            # Calculate and display both simple and weighted averages
            if (
                "Methylation_Percentage" in specific_sites.columns
                and not specific_sites["Methylation_Percentage"].isna().all()
            ):
                # Simple arithmetic mean
                simple_avg = specific_sites["Methylation_Percentage"].mean()

                # Coverage-weighted mean
                total_coverage = specific_sites["Total_Coverage"].sum()
                weighted_avg = (
                    (
                        specific_sites["Methylation_Percentage"]
                        * specific_sites["Total_Coverage"]
                    ).sum()
                    / total_coverage
                    if total_coverage > 0
                    else 0
                )

                ui.label("Methylation Summary:").classes(
                    "text-blue-600 font-medium mt-2"
                )
                ui.label(
                    f"• Simple average across CpG pairs: {simple_avg:.2f}%"
                ).classes("text-gray-600 ml-2")
                ui.label(f"• Coverage-weighted average: {weighted_avg:.2f}%").classes(
                    "text-gray-600 ml-2"
                )
                ui.label(f"• Total read coverage: {total_coverage}").classes(
                    "text-gray-600 ml-2"
                )
            else:
                ui.label("Unable to calculate average methylation (no data)").classes(
                    "text-amber-500 mt-2"
                )

        except Exception as e:
            logger.error(f"Error displaying specific sites: {e}")
            ui.label(f"Error analyzing specific sites: {str(e)}").classes(
                "text-red-500"
            )

    def extract_specific_sites(self, bed_file: str) -> pd.DataFrame:
        """
        Extracts specific methylation sites of interest from the bed file.

        The sites of interest are CpG pairs on chromosome 10:
        - 129467255/6 (C+ at 255, C- at 256)
        - 129467258/9 (C+ at 258, C- at 259)
        - 129467262/3 (C+ at 262, C- at 263)
        - 129467272/3 (C+ at 272, C- at 273)

        Note: BED coordinates are 0-based. We look for the exact positions.

        Args:
            bed_file (str): Path to the bed file.

        Returns:
            pd.DataFrame: DataFrame containing the collapsed CpG pair measurements.
        """
        try:
            # Read bed file with pandas - include strand information
            bed_data = pd.read_csv(bed_file, sep="\t", header=None)

            # Check if the file has the expected number of columns
            if bed_data.shape[1] >= 10:
                # Assign column names based on the observed format
                column_names = [
                    "Chromosome",
                    "Start",
                    "End",
                    "Name",
                    "Score",
                    "Strand",
                    "Start2",
                    "End2",
                    "RGB",
                    "Coverage_Info",
                ]

                # Ensure we only use as many column names as there are columns
                bed_data = bed_data.iloc[:, : len(column_names)]
                bed_data.columns = column_names

                # Parse coverage and methylation from Coverage_Info
                bed_data["Coverage"] = (
                    bed_data["Coverage_Info"].str.split().str[0].astype(int)
                )
                bed_data["Modified_Fraction"] = (
                    bed_data["Coverage_Info"].str.split().str[1].astype(float)
                )
            else:
                logger.error(
                    f"Unexpected bed file format with {bed_data.shape[1]} columns"
                )
                return pd.DataFrame()

            # Define the CpG pairs we're interested in (exact positions)
            cpg_pairs = [
                (129467255, 129467256),  # First CpG pair
                (129467258, 129467259),  # Second CpG pair
                (129467262, 129467263),  # Third CpG pair
                (129467272, 129467273),  # Fourth CpG pair
            ]

            # Initialize list to store processed CpG pairs
            processed_cpgs = []

            # Process each CpG pair
            for cpg_pos1, cpg_pos2 in cpg_pairs:
                # Get forward strand C at first position
                fwd_c = bed_data[
                    (bed_data["Chromosome"] == "chr10")
                    & (bed_data["Start"] == cpg_pos1 - 1)  # Convert to 0-based
                    & (bed_data["Strand"] == "+")
                ]

                # Get reverse strand C at second position
                rev_c = bed_data[
                    (bed_data["Chromosome"] == "chr10")
                    & (bed_data["Start"] == cpg_pos2 - 1)  # Convert to 0-based
                    & (bed_data["Strand"] == "-")
                ]

                if not fwd_c.empty and not rev_c.empty:
                    # Calculate combined metrics
                    total_coverage = (
                        fwd_c["Coverage"].iloc[0] + rev_c["Coverage"].iloc[0]
                    )
                    weighted_methylation = (
                        (
                            (
                                fwd_c["Coverage"].iloc[0]
                                * fwd_c["Modified_Fraction"].iloc[0]
                            )
                            + (
                                rev_c["Coverage"].iloc[0]
                                * rev_c["Modified_Fraction"].iloc[0]
                            )
                        )
                        / total_coverage
                        if total_coverage > 0
                        else 0
                    )

                    processed_cpgs.append(
                        {
                            "Chromosome": "chr10",
                            "Position": f"{cpg_pos1}/{cpg_pos2}",
                            "Coverage_Forward": fwd_c["Coverage"].iloc[0],
                            "Coverage_Reverse": rev_c["Coverage"].iloc[0],
                            "Total_Coverage": total_coverage,
                            "Methylation_Percentage": weighted_methylation,  # Already in percentage form
                            "Forward_Methylation": fwd_c["Modified_Fraction"].iloc[
                                0
                            ],  # Already in percentage form
                            "Reverse_Methylation": rev_c["Modified_Fraction"].iloc[
                                0
                            ],  # Already in percentage form
                        }
                    )

            # Create DataFrame from processed CpGs
            result_df = pd.DataFrame(processed_cpgs)

            if not result_df.empty:
                # Add descriptive labels
                def get_site_label(row):
                    pos = row["Position"]
                    site_map = {
                        "129467255/129467256": "Site 1",
                        "129467258/129467259": "Site 2",
                        "129467262/129467263": "Site 3",
                        "129467272/129467273": "Site 4",
                    }
                    return f"{site_map.get(pos, 'Unknown')} (CpG {pos})"

                result_df["Site_Label"] = result_df.apply(get_site_label, axis=1)
                result_df["Notes"] = (
                    "Combined methylation from both strands of CpG pair"
                )

            return result_df

        except Exception as e:
            logger.error(f"Error extracting specific sites: {e}")
            return pd.DataFrame()

    def save_specific_sites_data(
        self, specific_sites: pd.DataFrame, bed_file: str
    ) -> str:
        """
        Saves the specific methylation sites data to a CSV file.

        Args:
            specific_sites (pd.DataFrame): DataFrame containing the specific sites data.
            bed_file (str): Path to the original bed file.

        Returns:
            str: Path to the saved CSV file, or empty string if save failed.
        """
        try:
            if specific_sites.empty:
                logger.info("No specific sites data to save (empty DataFrame)")
                return ""

            # Create a filename based on the original bed file
            output_dir = os.path.dirname(bed_file)
            base_name = os.path.basename(bed_file)

            # Extract count from filename (handle potential format differences)
            try:
                count = base_name.split("_")[0]
                # Verify count is a number
                int(count)
            except (IndexError, ValueError):
                # If we can't extract a valid count, use a timestamp
                count = pd.Timestamp.now().strftime("%Y%m%d%H%M%S")
                logger.warning(
                    f"Could not extract count from filename {base_name}, using timestamp {count}"
                )

            csv_file = os.path.join(output_dir, f"{count}_specific_sites.csv")

            # Add timestamp
            specific_sites["Timestamp"] = pd.Timestamp.now().strftime(
                "%Y-%m-%d %H:%M:%S"
            )

            # Save to CSV
            specific_sites.to_csv(csv_file, index=False)
            logger.info(f"Saved specific sites data to {csv_file}")

            # Also append to a cumulative file for tracking over time
            cumulative_file = os.path.join(output_dir, "cumulative_specific_sites.csv")

            if os.path.exists(cumulative_file):
                # Append to existing file
                specific_sites.to_csv(
                    cumulative_file, mode="a", header=False, index=False
                )
                logger.info(
                    f"Appended data to existing cumulative file {cumulative_file}"
                )
            else:
                # Create new file
                specific_sites.to_csv(cumulative_file, index=False)
                logger.info(f"Created new cumulative file {cumulative_file}")

            return csv_file
        except Exception as e:
            logger.error(f"Error saving specific sites data: {e}")
            return ""

    def show_previous_data(self) -> None:
        """
        Displays previously analyzed data from the specified watch folder.

        Returns:
            None
        """
        # logger.debug(f"Showing previous data from {watchfolder}")
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        logger.info(f"Looking for previous data in {output}")

        for file in natsort.natsorted(os.listdir(output)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > self.last_seen:
                    results = pd.read_csv(os.path.join(output, file))
                    plot_out = os.path.join(output, file.replace(".csv", ".png"))
                    bed_file = os.path.join(output, f"{count}_mgmt.bed")

                    logger.info(
                        f"Processing file {file}, looking for bed file at {bed_file}"
                    )

                    if os.path.exists(plot_out):
                        self.mgmtplot.clear()
                    with self.mgmtplot.classes("w-full"):
                        ui.image(plot_out).props("fit=scale-up")

                    self.mgmtable.clear()
                    with self.mgmtable:
                        self.tabulate(results)

                        # Add bed file visualization if it exists
                        if os.path.exists(bed_file):
                            logger.info(f"Found bed file: {bed_file}")
                            ui.label("MGMT CpG Site Methylation Data").classes(
                                "text-lg font-medium mt-4 mb-2"
                            )
                            # self.display_bed_data(bed_file)

                            # Add specific sites analysis
                            self.display_specific_sites(bed_file)
                        else:
                            logger.warning(f"Bed file not found: {bed_file}")
                            ui.label(
                                f"Bed file not found: {os.path.basename(bed_file)}"
                            ).classes("text-amber-500 mt-2")

                            # Try alternative naming patterns
                            alt_bed_file = os.path.join(
                                output, file.replace(".csv", "_mgmt.bed")
                            )
                            if os.path.exists(alt_bed_file):
                                logger.info(
                                    f"Found alternative bed file: {alt_bed_file}"
                                )
                                ui.label("MGMT CpG Site Methylation Data").classes(
                                    "text-lg font-medium mt-4 mb-2"
                                )
                                # Add specific sites analysis
                                self.display_specific_sites(alt_bed_file)
                                # self.display_bed_data(alt_bed_file)

                    if self.summary:
                        with self.summary:
                            self.summary.clear()
                            with ui.card().classes("w-full p-4 mb-4"):
                                with ui.row().classes(
                                    "w-full items-center justify-between"
                                ):
                                    # Left side - MGMT Status
                                    with ui.column().classes("gap-2"):
                                        if "status" in results.columns:
                                            status = results["status"].values[0]
                                            average = float(
                                                results["average"].values[0]
                                            )
                                            pred = float(results["pred"].values[0])

                                            # Determine status styling
                                            if status.lower() == "methylated":
                                                status_color = "text-blue-600"
                                                status_bg = "bg-blue-100"
                                            else:
                                                status_color = "text-amber-600"
                                                status_bg = "bg-amber-100"

                                            ui.label(
                                                "MGMT Methylation Analysis"
                                            ).classes("text-lg font-medium")
                                            with ui.row().classes("items-center gap-2"):
                                                ui.label(f"Status: {status}").classes(
                                                    f"{status_color} font-medium"
                                                )
                                                ui.label(f"{pred:.1f}%").classes(
                                                    f"px-2 py-1 rounded {status_bg} {status_color}"
                                                )

                                    # Right side - Additional metrics
                                    with ui.column().classes("gap-2 text-right"):
                                        ui.label("Analysis Details").classes(
                                            "font-medium"
                                        )
                                        ui.label(
                                            f"Average Methylation: {average:.1f}%"
                                        ).classes("text-gray-600")
                                        ui.label(
                                            f"Prediction Score: {pred:.1f}%"
                                        ).classes("text-gray-600")

                                # Bottom row - Information
                                with ui.row().classes(
                                    "w-full mt-4 text-sm text-gray-500 justify-center"
                                ):
                                    ui.label(
                                        "MGMT status determined from methylation analysis of 137 CpG sites"
                                    )
                    self.last_seen = count


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
        self.MGMTbamfile = {}
        self.counter: int = 0
        self.last_seen: int = 0
        # logger.debug("Initializing MGMT_Object")
        super().__init__(*args, **kwargs)
        # Remove state tracking for MGMT Analysis
        # state.start_process("MGMT Analysis", ProcessType.BATCH)
        state.set_process_state("MGMT Analysis", ProcessState.WAITING_FOR_DATA)

    async def process_bam(self, bamfile: str, timestamp: str) -> None:
        """
        Processes the BAM file to extract and analyze MGMT sites.

        Args:
            bamfile (str): Path to the input BAM file.
            timestamp (str): Timestamp for the analysis.

        Returns:
            None
        """
        state.set_process_state("MGMT Analysis", ProcessState.RUNNING)
        try:
            # result = await run.cpu_bound(
            result = has_reads(bamfile, "chr10", 129467242, 129467244)
            # )

            if result:
                MGMT_BED: str = f"{HVPATH}/bin/mgmt_hg38.bed"
                tempbamfile = tempfile.NamedTemporaryFile(
                    dir=self.check_and_create_folder(self.output, self.sampleID),
                    suffix=".bam",
                )

                async def bedtools_background_work(bamfile, MGMT_BED, tempbamfile):
                    try:
                        await run.cpu_bound(
                            run_bedtools, bamfile, MGMT_BED, tempbamfile.name
                        )
                    except Exception as e:
                        logger.error(f"Error in bedtools_background_work: {e}")
                        # return

                await background_tasks.create(
                    bedtools_background_work(bamfile, MGMT_BED, tempbamfile)
                )

                try:
                    if (
                        pysam.AlignmentFile(tempbamfile.name, "rb").count(
                            until_eof=True
                        )
                        > 0
                    ):

                        async def mgmt_bam_background_work(
                            sampleID, MGMTbamfile, tempbamfile
                        ):
                            if sampleID not in MGMTbamfile.keys():
                                # if not self.MGMTbamfile:
                                MGMTbamfile[sampleID] = os.path.join(
                                    self.check_and_create_folder(self.output, sampleID),
                                    "mgmt.bam",
                                )
                                shutil.copy2(tempbamfile.name, MGMTbamfile[sampleID])
                                os.remove(f"{tempbamfile.name}.bai")
                            else:
                                tempbamholder = tempfile.NamedTemporaryFile(
                                    dir=self.check_and_create_folder(
                                        self.output, sampleID
                                    ),
                                    suffix=".bam",
                                )
                                pysam.cat(
                                    "-o",
                                    tempbamholder.name,
                                    MGMTbamfile[sampleID],
                                    tempbamfile.name,
                                )
                                shutil.copy2(tempbamholder.name, MGMTbamfile[sampleID])
                                try:
                                    os.remove(f"{tempbamholder.name}.bai")
                                    os.remove(f"{tempbamfile.name}.bai")
                                except FileNotFoundError:
                                    pass

                        await background_tasks.create(
                            mgmt_bam_background_work(
                                self.sampleID, self.MGMTbamfile, tempbamfile
                            )
                        )
                        tempmgmtdir = tempfile.TemporaryDirectory(
                            dir=self.check_and_create_folder(self.output, self.sampleID)
                        )
                        self.counter += 1

                        output_mgmt_bed = os.path.join(
                            self.check_and_create_folder(self.output, self.sampleID),
                            f"{self.counter}_mgmt.bed",
                        )

                        async def modkit_background_work(
                            tempmgmtdir, MGMTbamfile, threads, output_mgmt_bed
                        ):
                            await run.cpu_bound(
                                run_modkit,
                                tempmgmtdir.name,
                                MGMTbamfile,
                                threads,
                                output_mgmt_bed,
                            )

                        await background_tasks.create(
                            modkit_background_work(
                                tempmgmtdir,
                                self.MGMTbamfile[self.sampleID],
                                self.threads,
                                output_mgmt_bed,
                            )
                        )

                        try:
                            if os.path.exists(
                                os.path.join(
                                    tempmgmtdir.name, "live_analysis_mgmt_status.csv"
                                )
                            ):
                                results = pd.read_csv(
                                    os.path.join(
                                        tempmgmtdir.name,
                                        "live_analysis_mgmt_status.csv",
                                    )
                                )
                                # self.counter += 1
                                plot_out = os.path.join(
                                    self.check_and_create_folder(
                                        self.output, self.sampleID
                                    ),
                                    f"{self.counter}_mgmt.png",
                                )

                                async def methylartist_background_work(
                                    tempmgmtdir, plot_out
                                ):
                                    await run.cpu_bound(
                                        run_methylartist, tempmgmtdir.name, plot_out
                                    )

                                await background_tasks.create(
                                    methylartist_background_work(tempmgmtdir, plot_out)
                                )

                                results.to_csv(
                                    os.path.join(
                                        self.check_and_create_folder(
                                            self.output, self.sampleID
                                        ),
                                        f"{self.counter}_mgmt.csv",
                                    ),
                                    index=False,
                                )
                        except Exception:
                            # logger.error(f"Error processing results: {e}")
                            raise
                    else:
                        os.remove(f"{tempbamfile.name}.bai")
                except Exception:
                    # logger.error(f"Error in BAM file processing: {e}")
                    raise
            else:
                self.running = False
        finally:
            state.set_process_state("MGMT Analysis", ProcessState.WAITING_FOR_DATA)

    def get_report(self, watchfolder: str) -> Tuple[pd.DataFrame, str, str]:
        """
        Generates a report from the analysis results.

        Args:
            watchfolder (str): Path to the folder containing previous analysis results.

        Returns:
            Tuple[pd.DataFrame, str, str]: DataFrame with results, path to plot, and summary string.
        """
        # logger.debug(f"Generating report from {watchfolder}")
        results = pd.DataFrame()
        plot_out = ""
        summary = ""
        for file in natsort.natsorted(os.listdir(watchfolder)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > self.last_seen:
                    results = pd.read_csv(os.path.join(watchfolder, file))
                    plot_out = os.path.join(watchfolder, file.replace(".csv", ".png"))
                    summary = f"Current MGMT status: {results['status'].values[0]}"
        return results, plot_out, summary


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
    # logger.debug(
    #    f"Starting MGMT analysis application on port {port} with {threads} threads"
    # )
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
def main(
    port: int,
    threads: int,
    watchfolder: Optional[str],
    output: Optional[str],
    browse: bool,
) -> None:
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
    # logger.debug(
    #    f"Running main function with port={port}, threads={threads}, watchfolder={watchfolder}, output={output}, browse={browse}"
    # )
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
