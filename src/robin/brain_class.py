"""
Main Module for Initializing and Configuring the Brain Tumor Classification Application

This module sets up and configures the Brain Tumor Classification application using the NiceGUI framework (https://nicegui.io/). The application ensures that all long-running processes operate in the background, preventing the main thread from being blocked.

The app creates an instance of the `BrainMeth` class for each specific run, taking arguments from either the command line or a configuration file. Shared arguments are stored in `nicegui.app.storage.general` for accessibility across the application.

Key Components:

1. **BrainMeth Class**:

   - Sets up the application.
   - Manages BAM file processing and classification tasks.
   - Handles the configuration of subpages for various analysis types (e.g., MGMT, Sturgeon, NanoDX).

2. **BAM File Handling**:

   - `check_bam(bamfile)`: Checks a BAM file and returns its attributes.
   - Uses `watchdog.observers.Observer` to monitor a directory for new BAM files.

3. **Subpage Objects**:

   - `MGMT_Object`
   - `Sturgeon_object`
   - `NanoDX_object`
   - `RandomForest_object`
   - `CNVAnalysis`
   - `TargetCoverage`
   - `FusionObject`

4. **Utility Functions**:

   - `check_bam(bamfile)`: Verifies BAM file attributes.
   - `LocalFilePicker`: Allows for local file selection.

5. **Queues for BAM Processing**:

   - `bam_tracking`
   - `bamforcns`
   - `bamforsturgeon`
   - `bamfornanodx`
   - `bamforcnv`
   - `bamfortargetcoverage`
   - `bamformgmt`
   - `bamforfusions`

6. **Configuration**:

   - Uses `click` for command-line interface options.
   - Loads configuration from `config.ini` or specified file.

Dependencies:

- `nicegui` (ui, app, run)
- `robin.utilities.bam_handler.BamEventHandler`
- `robin.subpages` (MGMT_Object, Sturgeon_object, NanoDX_object, RandomForest_object, CNVAnalysis, TargetCoverage, FusionObject)
- `robin.utilities.local_file_picker.LocalFilePicker`
- `robin.utilities.ReadBam.ReadBam`
- `watchdog.observers.Observer`
- `pathlib.Path`
- `queue.Queue`
- `pandas`
- `asyncio`
- `collections.Counter`
- `datetime`
- `os`
"""

import logging
from pathlib import Path
from queue import Queue
import pandas as pd
import asyncio
import pysam
from collections import Counter
from datetime import datetime
from dateutil import parser
import pytz
import os
import shutil
import tempfile
from alive_progress import alive_bar
from nicegui import ui, app, run
import subprocess

from robin import resources


from robin.utilities.bam_handler import BamEventHandler
from robin.subpages.MGMT_object import MGMT_Object
from robin.subpages.Sturgeon_object import (
    Sturgeon_object,
    run_sturgeon_inputtobed,
    pysam_cat,
    run_sturgeon_merge_probes,
    run_probes_methyl_calls,
)
from robin.subpages.NanoDX_object import NanoDX_object
from robin.subpages.RandomForest_object import RandomForest_object
from robin.subpages.CNV_object import CNVAnalysis
from robin.subpages.TargetCoverage_object import TargetCoverage
from robin.subpages.Fusion_object import FusionObject
from robin.utilities.local_file_picker import LocalFilePicker
from robin.utilities.ReadBam import ReadBam
from robin.utilities.bed_file import BedTree
from robin.reporting.report import create_pdf
from robin.reporting.sections.disclaimer_text import EXTENDED_DISCLAIMER_TEXT

from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)

from watchdog.observers import Observer
from typing import List, Tuple, Optional, Dict, Any


import pyranges as pr  # For fast interval-based filtering


def merge_modkit_files(
    new_files: List[str], existing_file: str, output_file: str, filter_bed_file: str
):
    """
    Efficiently merges new modkit pileup files with an existing dataset using incremental updates.
    Filters out any positions that do not overlap with the given BED file using pyranges for speed.

    Parameters:
        new_files (List[str]): List of new input files (modkit format).
        existing_file (str): Path to the previously merged dataset (Parquet format).
        output_file (str): Path to save the updated merged dataset.
        filter_bed_file (str): Path to BED file to filter positions (can be .gz).

    Returns:
        None (Saves the merged DataFrame to `output_file`).
    """
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

    numeric_cols = [
        "chromStart",
        "chromEnd",
        "score_bed",
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

    # Determine if the filter BED file is gzipped
    compression_type = "gzip" if filter_bed_file.endswith(".gz") else None

    # Load filter BED file and convert to pyranges
    logging.debug("🔍 Loading and indexing filter BED file for fast lookup...")
    filter_df = pd.read_csv(
        filter_bed_file,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End", "cg_label"],
        compression=compression_type,
        dtype={"Chromosome": str},
    )

    # Convert to pyranges for fast interval-based filtering
    filter_ranges = pr.PyRanges(filter_df[["Chromosome", "Start", "End"]])

    # Load existing dataset if it exists
    if os.path.exists(existing_file):
        logging.debug(f"🔹 Loading existing dataset: {existing_file}")
        existing_df = pd.read_parquet(existing_file)
    else:
        logging.debug("🔹 No existing dataset found. Creating a new dataset.")
        existing_df = pd.DataFrame(columns=column_names)

    # Read and merge only new files
    new_data = []
    for file in new_files:
        logging.debug(f"📂 Processing file: {file}")

        # Read file efficiently
        modkit_df = pd.read_csv(
            file, sep="\s+", header=None, names=column_names, dtype=str
        )

        # Convert numeric columns properly
        modkit_df[numeric_cols] = modkit_df[numeric_cols].apply(
            pd.to_numeric, errors="coerce"
        )

        # Convert to pyranges format for fast overlap filtering
        modkit_ranges = pr.PyRanges(
            modkit_df.rename(
                columns={
                    "chrom": "Chromosome",
                    "chromStart": "Start",
                    "chromEnd": "End",
                }
            )
        )

        # Perform fast overlap filtering
        filtered_ranges = modkit_ranges.intersect(filter_ranges)

        # Convert back to DataFrame
        modkit_filtered = filtered_ranges.df.rename(
            columns={"Chromosome": "chrom", "Start": "chromStart", "End": "chromEnd"}
        )

        # Append to list (avoids multiple DataFrame copies)
        new_data.append(modkit_filtered)

    if not new_data:
        logging.debug("⚠️ No new data to merge after filtering. Exiting.")
        return

    # Combine new data into a single DataFrame
    new_df = pd.concat(new_data, ignore_index=True)

    # If existing dataset is empty, save new data and exit
    if existing_df.empty:
        logging.debug("✅ No previous data found. Saving only new filtered data.")
        new_df.to_parquet(output_file, index=False)
        return new_df

    # Identify unique positions for faster merging
    existing_positions = set(
        zip(existing_df["chrom"], existing_df["chromStart"], existing_df["chromEnd"])
    )
    new_positions = set(zip(new_df["chrom"], new_df["chromStart"], new_df["chromEnd"]))

    # Find new data that isn't already in the existing dataset
    unique_positions = new_positions - existing_positions
    if not unique_positions:
        logging.debug("⚠️ All new data already exists. No updates needed.")
        return existing_df

    # Filter only truly new data
    new_df = new_df[
        new_df.set_index(["chrom", "chromStart", "chromEnd"]).index.isin(
            unique_positions
        )
    ]

    # Debugging Step: Check for missing columns before aggregation
    logging.debug("\n🔍 Checking DataFrame before aggregation:")
    # logging.debug(new_df.info())

    # **Fix Aggregation Issue: Handle missing categorical values**
    for col in ["mod_code", "strand", "thickStart", "thickEnd", "color"]:
        if col in new_df.columns:
            new_df[col] = new_df[col].ffill()  # Forward fill missing values

    # Merge new data with existing data
    merged_df = pd.concat([existing_df, new_df], ignore_index=True)

    # Aggregate count-based columns only when necessary
    sum_columns = [
        "n_mod",
        "n_canonical",
        "n_othermod",
        "n_delete",
        "n_fail",
        "n_diff",
        "n_nocall",
    ]

    # **⚡ FIX: Use `observed=True` to match behavior and avoid warnings**
    merged_df = (
        merged_df.groupby(
            ["chrom", "chromStart", "chromEnd"], as_index=False, observed=True
        )
        .agg(
            {
                "mod_code": "first",
                "score_bed": "mean",
                "strand": "first",
                "thickStart": "first",
                "thickEnd": "first",
                "color": "first",
                "valid_cov": "sum",
                "percent_modified": "mean",
                **{col: "sum" for col in sum_columns},  # Sum count-based columns
            }
        )
        .reset_index(drop=True)
    )  # Ensure index is reset correctly

    # Save updated dataset
    merged_df.to_parquet(output_file, index=False)
    logging.debug(f"✅ Incremental update saved to: {output_file}")
    return merged_df


def run_modkit(sortfile: str, temp: str, threads: int) -> None:
    """
    Executes modkit on a bam file and extracts the methylation data.

    Args:
        sortfile (str): Path to the sorted BAM file.
        temp (str): Path to the temporary output file.
        threads (int): Number of threads to use.
    """
    cmd = [
        "modkit",
        "pileup",
        "-t",
        str(threads),  # Ensure threads is a string
        "--filter-threshold",
        "0.73",
        "--combine-mods",
        sortfile,
        temp,
        # "--suppress-progress"
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)
    # print("STDOUT:", result.stdout)
    # print("STDERR:", result.stderr)


def run_samtools_sort(
    file: str, tomerge: List[str], sortfile: str, threads: int
) -> None:
    """
    Sorts BAM files using Samtools.

    Args:
        file (str): Path to the output BAM file.
        tomerge (List[str]): List of BAM files to merge.
        sortfile (str): Path to the sorted BAM file.
        threads (int): Number of threads to use.
    """
    # logger.debug(
    #    f"Running samtools sort with the following parameters: file={file}, tomerge={tomerge}, sortfile={sortfile}, threads={threads}"
    # )
    pysam.cat("-o", file, *tomerge)
    pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)
    # logger.debug("samtools sort command executed successfully.")


def check_bam(bamfile):
    """
    Check a BAM file and return its attributes.

    :param bamfile: Path to the BAM file.
    :return: Tuple containing baminfo and bamdata.
    """
    logging.info(f"Checking BAM file: {bamfile}")
    pysam.index(bamfile)
    bam = ReadBam(bamfile)
    baminfo = bam.process_reads()
    bamdata = bam.summary()
    logging.info(f"BAM file processed: {bamfile}")
    return baminfo, bamdata


def sort_bams(files_and_timestamps, watchfolder, file_endings, simtime):
    """
    Sort BAM files by timestamp.

    :param files_and_timestamps: List to store sorted files and timestamps
    :param watchfolder: Folder to watch for BAM files
    :param file_endings: Set of file endings to look for
    :param simtime: Whether to simulate time
    :return: Sorted list of files and timestamps
    """
    import bisect

    def insert_sorted(file_timestamp_tuple):
        file, timestamp, elapsed = file_timestamp_tuple
        datetime_obj = datetime.fromisoformat(timestamp)
        bisect.insort(files_and_timestamps, (datetime_obj, file, elapsed))

    for path, dirs, files in os.walk(watchfolder):
        with alive_bar(len(files)) as bar:
            for f in files:
                if "".join(Path(f).suffix) in file_endings:
                    logging.info(f"Reading BAM file: {os.path.join(path, f)}")
                    if simtime:
                        bam = ReadBam(os.path.join(path, f))
                        baminfo = bam.process_reads()
                        insert_sorted(
                            (
                                os.path.join(path, f),
                                baminfo["last_start"],
                                baminfo["elapsed_time"],
                            )
                        )
                    else:
                        filetime = datetime.fromtimestamp(
                            os.path.getctime(os.path.join(path, f))
                        )
                        elapsedtime = datetime.now() - filetime
                        insert_sorted(
                            (os.path.join(path, f), filetime.isoformat(), elapsedtime)
                        )
                bar()
    return files_and_timestamps


class BrainMeth:
    def __init__(
        self,
        threads=4,
        force_sampleid=None,
        kit=None,
        centreID=None,
        simtime=False,
        watchfolder=None,
        output=None,
        sequencing_summary=None,
        target_panel=None,
        showerrors=False,
        browse=False,
        exclude=[],
        minknow_connection=None,
        reference=None,
        bed_file=None,
        mainuuid=None,
        readfish_toml=None,
    ):
        """
        Initialize the BrainMeth class.

        :param threads: Number of threads to use.
        :param force_sampleid: Force the use of a specific sample ID.
        :param simtime: Simulate time for testing.
        :param watchfolder: Folder to watch for BAM files.
        :param output: Output directory.
        :param sequencing_summary: Path to the sequencing summary file.
        :param target_panel: Target panel for analysis.
        :param showerrors: Show errors during processing.
        :param browse: Enable browsing mode.
        :param exclude: List of analysis types to exclude.
        :param minknow_connection: Connection to MinKNOW.
        :param reference: Reference genome.
        :param mainuuid: Main UUID for the application instance.
        """
        self.mainuuid = mainuuid
        self.force_sampleid = force_sampleid
        self.kit = kit
        self.centreID = centreID
        self.threads = threads
        self.simtime = simtime
        self.watchfolder = watchfolder
        self.output = output
        self.sequencing_summary = sequencing_summary
        self.target_panel = target_panel
        self.showerrors = showerrors
        self.browse = browse
        self.exclude = exclude
        self.minknow_connection = minknow_connection
        self.reference = reference
        self.bed_file = bed_file
        self.observer = None
        if self.browse:
            self.runsfolder = self.output
        self.sampleID = None
        self.readfish_toml = readfish_toml
        self.watchdogbamqueue = Queue()
        self.dataDir = {}
        self.bedDir = {}
        self.first_run = {}
        self.cpgs_master_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "sturg_nanodx_cpgs_0125.bed.gz",
        )
        if self.reference:
            self.NewBed = BedTree(
                preserve_original_tree=True,
                reference_file=f"{self.reference}.fai",
                readfish_toml=self.readfish_toml,
            )
            if self.bed_file:
                self.NewBed.load_from_file(self.bed_file)
        else:
            self.NewBed = None

        logging.info(f"BrainMeth initialized with UUID: {self.mainuuid}")

    def configure_storage(self, sample_id):
        """
        Configure storage for the application.

        :param sample_id: Sample ID to configure storage for.
        """
        logging.info(f"Configuring storage for sample ID: {sample_id}")
        app.storage.general[self.mainuuid]["samples"][sample_id]["file_counters"] = {
            "bam_passed": 0,
            "bam_failed": 0,
            "mapped_count": 0,
            "pass_mapped_count": 0,
            "fail_mapped_count": 0,
            "unmapped_count": 0,
            "pass_unmapped_count": 0,
            "fail_unmapped_count": 0,
            "pass_bases_count": 0,
            "fail_bases_count": 0,
            "bases_count": 0,
            # New counters for read numbers
            "mapped_reads_num": 0,
            "unmapped_reads_num": 0,
            "pass_mapped_reads_num": 0,
            "fail_mapped_reads_num": 0,
            "pass_unmapped_reads_num": 0,
            "fail_unmapped_reads_num": 0,
            # New counters for bases in each category
            "mapped_bases": 0,
            "unmapped_bases": 0,
            "pass_mapped_bases": 0,
            "fail_mapped_bases": 0,
            "pass_unmapped_bases": 0,
            "fail_unmapped_bases": 0,
        }
        app.storage.general[self.mainuuid]["samples"][sample_id]["devices"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["basecall_models"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["run_time"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["flowcell_ids"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["sample_ids"] = []

    def update_counter(self, sample_id, counter_name, value):
        """
        Safely update a counter value.

        :param sample_id: Sample ID to update counter for
        :param counter_name: Name of the counter to update
        :param value: Value to add to the counter
        """
        current = app.storage.general[self.mainuuid]["samples"][sample_id][
            "file_counters"
        ][counter_name]
        new_value = max(0, current + value)  # Ensure we don't go below 0
        app.storage.general[self.mainuuid]["samples"][sample_id]["file_counters"][
            counter_name
        ] = new_value
        logging.debug(
            f"Counter {counter_name} updated for sample {sample_id}: {current} -> {new_value}"
        )

    async def start_background(self):
        """
        Start background tasks for BAM processing and analysis.
        """
        logging.info(f"Starting background tasks for UUID: {self.mainuuid}")
        app.storage.general[self.mainuuid]["bam_count"] = {}
        app.storage.general[self.mainuuid]["bam_count"] = Counter(counter=0)
        app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
        app.storage.general[self.mainuuid]["bam_count"][
            "total_files"
        ] = 0  # Add total files counter
        self.bam_tracking = Queue()
        self.bamforcns = Queue()
        self.bamforsturgeon = Queue()
        self.bamfornanodx = Queue()
        self.bamforpannanodx = Queue()
        self.bamforcnv = Queue()
        self.bamfortargetcoverage = Queue()
        self.bamformgmt = Queue()
        self.bamforfusions = Queue()
        self.bamforbigbadmerge = Queue()
        self.mergecounter = 0  # This is the queue for the big bad merge
        if self.watchfolder:
            logging.info(f"Adding a watchfolder: {self.watchfolder}")
            await self.add_watchfolder(self.watchfolder)

        common_args = {
            "threads": self.threads,
            "output": self.output,
            "progress": True,
            "browse": self.browse,
            "uuid": self.mainuuid,
            "force_sampleid": self.force_sampleid,
        }

        if "sturgeon" not in self.exclude:
            self.Sturgeon = Sturgeon_object(
                analysis_name="STURGEON",
                batch=True,
                bamqueue=self.bamforsturgeon,
                **common_args,
            )
            self.Sturgeon.process_data()

        if "nanodx" not in self.exclude:
            self.NanoDX = NanoDX_object(
                analysis_name="NANODX",
                batch=True,
                bamqueue=self.bamfornanodx,
                **common_args,
            )
            self.NanoDX.process_data()
        if "pannanodx" not in self.exclude:
            self.panNanoDX = NanoDX_object(
                analysis_name="PANNANODX",
                batch=True,
                bamqueue=self.bamforpannanodx,
                model="pancan_devel_v5i_NN.pkl",
                **common_args,
            )
            self.panNanoDX.process_data()
        if "forest" not in self.exclude:
            self.RandomForest = RandomForest_object(
                analysis_name="FOREST",
                batch=True,
                showerrors=self.showerrors,
                bamqueue=self.bamforcns,
                **common_args,
            )
            self.RandomForest.process_data()

        if "cnv" not in self.exclude:
            self.CNV = CNVAnalysis(
                analysis_name="CNV",
                bamqueue=self.bamforcnv,
                target_panel=self.target_panel,
                reference_file=self.reference,
                bed_file=self.bed_file,
                readfish_toml=self.readfish_toml,
                NewBed=self.NewBed,
                **common_args,
            )
            self.CNV.process_data()

        if "coverage" not in self.exclude:
            self.Target_Coverage = TargetCoverage(
                analysis_name="COVERAGE",
                showerrors=self.showerrors,
                bamqueue=self.bamfortargetcoverage,
                target_panel=self.target_panel,
                reference=self.reference,
                **common_args,
            )
            self.Target_Coverage.process_data()

        if "mgmt" not in self.exclude:
            self.MGMT_panel = MGMT_Object(
                analysis_name="MGMT", bamqueue=self.bamformgmt, **common_args
            )
            self.MGMT_panel.process_data()

        if "fusion" not in self.exclude:
            self.Fusion_panel = FusionObject(
                analysis_name="FUSION",
                bamqueue=self.bamforfusions,
                target_panel=self.target_panel,
                NewBed=self.NewBed,
                **common_args,
            )
            self.Fusion_panel.process_data()

    async def render_ui(self, sample_id=None):
        """
        Render the user interface.

        :param sample_id: Sample ID to render the UI for.
        """
        if sample_id:
            self.output = self.check_and_create_folder(self.output, sample_id)

        if not self.browse:
            await self.information_panel(sample_id=sample_id)
        else:
            ui.label("Browse mode enabled. Please choose a folder to see data from.")
            ui.button("Choose file", on_click=self.pick_file, icon="folder")
            self.content = ui.column().classes("w-full")

    async def add_watchfolder(self, watchfolder):
        """
        Add a watch folder to monitor for BAM files.

        :param watchfolder: Path to the folder to watch.
        """
        if not self.observer:
            self.watchfolder = watchfolder
            if "file" not in app.storage.general[self.mainuuid]["bam_count"].keys():
                app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
            # ToDo: Replace this with a queue.
            self.event_handler = BamEventHandler(
                self.watchdogbamqueue
                # app.storage.general[self.mainuuid]["bam_count"]
            )
            self.observer = Observer()
            self.observer.schedule(self.event_handler, self.watchfolder, recursive=True)
            ui.timer(1, callback=self.check_existing_bams, once=True)
            self.observer.start()

            ui.timer(1, callback=self.background_process_bams, once=True)

            logging.info("Watchfolder setup and added")

    async def pick_file(self) -> None:
        """
        Handle file picking for monitoring.

        :return: None
        """
        logging.info("Starting file picker in browse mode")
        ui.notify("Select a folder to monitor")
        result = await LocalFilePicker(f"{self.runsfolder}", multiple=True)
        if result:
            logging.info(f"Selected folder: {result}")
            ui.notify(f"You selected {result}")
            self.content.clear()
            with self.content:
                self.output, self.sampleID = os.path.split(result[0])
                logging.info(
                    f"Split path - output: {self.output}, sampleID: {self.sampleID}"
                )

                # Initialize storage for browse mode
                if self.sampleID not in app.storage.general[self.mainuuid]["samples"]:
                    logging.info(f"Initializing storage for sample {self.sampleID}")
                    app.storage.general[self.mainuuid]["samples"][self.sampleID] = {}
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "file_counters"
                    ] = {
                        "bam_passed": 0,
                        "bam_failed": 0,
                        "mapped_count": 0,
                        "unmapped_count": 0,
                        "pass_mapped_count": 0,
                        "fail_mapped_count": 0,
                        "pass_bases_count": 0,
                        "fail_bases_count": 0,
                        "bases_count": 0,
                    }
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "devices"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "basecall_models"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "run_time"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "flowcell_ids"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "sample_ids"
                    ] = [self.sampleID]
                else:
                    logging.info(f"Storage already exists for sample {self.sampleID}")

                if (
                    self.sampleID
                    not in app.storage.general[self.mainuuid]["sample_list"]
                ):
                    logging.info(f"Adding {self.sampleID} to sample list")
                    app.storage.general[self.mainuuid]["sample_list"].append(
                        self.sampleID
                    )

                # Try to load existing data from master.csv if it exists
                master_csv = os.path.join(result[0], "master.csv")
                logging.info(f"Looking for master.csv at: {master_csv}")
                if os.path.exists(master_csv):
                    try:
                        logging.info("Found master.csv, loading data")
                        df = pd.read_csv(master_csv)
                        # Update storage with data from master.csv
                        for key in df.keys():
                            if (
                                key
                                in app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]
                            ):
                                logging.info(f"Updating {key} from master.csv")
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ][key] = df[key].tolist()
                    except Exception as e:
                        logging.error(
                            f"Error loading master.csv: {str(e)}", exc_info=True
                        )
                else:
                    logging.warning(f"No master.csv found at {master_csv}")

                logging.info("Calling information_panel")
                await self.information_panel(sample_id=self.sampleID)
        else:
            logging.warning("No folder selected in file picker")

    def replay(self):
        """
        Replay the data processing.

        :return: None
        """
        ui.notify("Replaying data")
        self.replaycontrol.visible = False
        self.rcns2_worker.replay_prior_data()
        self.sturgeon_worker.replay_prior_data()

    @property
    def min_start_time(self):
        """
        Get the minimum start time for runs.

        :return: Minimum start time formatted as string.
        """
        if len(self.run_time) > 0:
            dt = datetime.fromisoformat(min(self.run_time))
            formatted_string = dt.strftime("%Y-%m-%d %H:%M")
            return formatted_string
        else:
            return 0

    @ui.refreshable
    def show_list(self):
        """
        Show the list of samples.

        :return: None
        """
        if len(app.storage.general[self.mainuuid]["sample_list"]) == 0:
            with ui.card().classes(
                "max-w-sm rounded overflow-hidden shadow-lg bg-white"
            ).style("box-shadow: 0 0 10px 2px rgba(0, 0, 255, 0.5);"):
                with ui.row():
                    ui.label("No samples found waiting.")
                    ui.spinner("dots", size="lg", color="#6E93D6")
        elif len(app.storage.general[self.mainuuid]["sample_list"]) == 1:
            if self.sampleID != app.storage.general[self.mainuuid]["sample_list"][0]:
                ui.navigate.to(
                    f"/live/{app.storage.general[self.mainuuid]['sample_list'][0]}"
                )
        else:
            myrow = ui.row().classes("w-full")
            with myrow:
                for item in app.storage.general[self.mainuuid]["sample_list"]:
                    if item == self.sampleID:
                        card = (
                            ui.card()
                            .classes("max-w-sm rounded shadow-lg")
                            .style("box-shadow: 0 0 10px 2px rgba(0, 0, 255, 0.5);")
                        )
                    else:
                        card = ui.card().classes("max-w-sm rounded shadow-lg")
                    with card:
                        with ui.expansion(f"{item}", icon="topic").style(
                            "font-size: 120%; font-weight: 300"
                        ).classes("w-full"):
                            for element in [
                                "devices",
                                "basecall_models",
                                "run_time",
                                "flowcell_ids",
                            ]:
                                if (
                                    element
                                    in app.storage.general[self.mainuuid]["samples"][
                                        item
                                    ]
                                ):
                                    ui.html().bind_content_from(
                                        app.storage.general[self.mainuuid]["samples"][
                                            item
                                        ],
                                        element,
                                        backward=lambda n, e=element: (
                                            f"<strong>{e}</strong>:{n[0]}"
                                            if n
                                            else "Default Value"
                                        ),
                                    ).style("font-size: 75%; font-weight: 100").classes(
                                        "font-normal text-sm w-full"
                                    )

                            ui.html().bind_content_from(
                                app.storage.general[self.mainuuid]["samples"][item][
                                    "file_counters"
                                ],
                                "bam_passed",
                                backward=lambda n: f"<span class='inline-block bg-gray-200 rounded-full px-3 py-1 text-sm font-semibold text-gray-700 mr-2 mb-2'>#bam passed: {n}</span>",
                            )
                            ui.html().bind_content_from(
                                app.storage.general[self.mainuuid]["samples"][item][
                                    "file_counters"
                                ],
                                "bam_failed",
                                backward=lambda n: f"<span class='inline-block bg-gray-200 rounded-full px-3 py-1 text-sm font-semibold text-gray-700 mr-2 mb-2'>#bam failed: {n}</span>",
                            )

                            ui.button(
                                "View Data",
                                on_click=lambda i=item: ui.navigate.to(f"/live/{i}"),
                            )

    async def information_panel(self, sample_id=None):
        """
        Display the main information panel with analysis results.
        """
        try:
            logging.info(f"Starting information_panel with sample_id: {sample_id}")
            if sample_id:
                self.sampleID = sample_id
                logging.info(f"Set self.sampleID to {self.sampleID}")
            else:
                logging.warning("No sample_id provided")
                if self.browse:
                    ui.label("Please select a sample folder to view").classes(
                        "text-xl text-gray-600 my-4"
                    )
                    return

            if not self.browse:
                logging.info("Live mode - showing sample list")
                self.show_list()
                app.storage.general[self.mainuuid]["sample_list"].on_change(
                    self.show_list.refresh
                )

            if (
                self.sampleID
                and not self.browse
                and self.sampleID not in app.storage.general[self.mainuuid]["samples"]
            ):
                logging.error(f"Sample {self.sampleID} not found in storage")
                ui.notify(f"Sample {self.sampleID} not found")
                ui.navigate.to("/live")
                return

            if not self.sampleID:
                logging.warning("No sample ID set, cannot render panel")
                return

            logging.info("Starting UI rendering")
            if self.sampleID in app.storage.general[self.mainuuid]["samples"]:
                logging.info(
                    f"Storage state for sample {self.sampleID}: {app.storage.general[self.mainuuid]['samples'][self.sampleID]}"
                )
            else:
                logging.warning(f"No storage found for sample {self.sampleID}")
                return

            with ui.column().classes("w-full px-6"):
                logging.info("Rendering title section")
                # Title and description
                with ui.column().classes("space-y-2 mb-4"):
                    ui.label("CNS Tumor Methylation Classification").classes(
                        "text-2xl font-medium text-gray-900"
                    )
                    ui.label(
                        "This tool enables classification of brain tumors in real time from Oxford Nanopore Data."
                    ).classes("text-gray-600")
                logging.info("Title section rendered")

                # Run Information Card
                with ui.card().classes("w-full p-4 bg-white rounded-lg shadow-sm"):
                    ui.label("Run Information").classes(
                        "text-lg font-medium text-gray-900 mb-4"
                    )
                    with ui.row().classes("items-center gap-6"):
                        # Run Time
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "run_time"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("schedule").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "run_time",
                                    backward=lambda n: (
                                        f"Run: {parser.parse(n[0]).strftime('%Y-%m-%d %H:%M')}"
                                        if n
                                        else None
                                    ),
                                )

                        # Model
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "basecall_models"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("model_training").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "basecall_models",
                                    backward=lambda n: f"Model: {n[0]}" if n else None,
                                )

                        # Device
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "devices"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("memory").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "devices",
                                    backward=lambda n: f"Device: {n[0]}" if n else None,
                                )

                        # Flow Cell
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "flowcell_ids"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("grid_4x4").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "flowcell_ids",
                                    backward=lambda n: (
                                        f"Flow Cell: {n[0]}" if n else None
                                    ),
                                )

                        # Sample
                        if self.sampleID:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("label").classes("text-gray-400")
                                ui.label(f"Sample: {self.sampleID}").classes(
                                    "text-gray-600"
                                )

                # Results Summary Section
                with ui.column().classes("space-y-4"):
                    # Classification Results - Single row with equal width columns
                    with ui.card().classes("w-full p-3 sm:p-4 bg-gray-50 rounded-lg"):
                        ui.label("Classification Summary").classes(
                            "font-medium text-gray-900 mb-3"
                        )
                        with ui.grid().classes(
                            "grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4"
                        ):
                            if "sturgeon" not in self.exclude:
                                sturgeonsummary = ui.column().classes("space-y-2")
                            if "nanodx" not in self.exclude:
                                nanodxsummary = ui.column().classes("space-y-2")
                            if "pannanodx" not in self.exclude:
                                pannanodxsummary = ui.column().classes("space-y-2")
                            if "forest" not in self.exclude:
                                forestsummary = ui.column().classes("space-y-2")

                    # Analysis Results - Two columns grid
                    with ui.grid().classes(
                        "grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4"
                    ):
                        if "coverage" not in self.exclude:
                            coverage = ui.column().classes("space-y-1")
                        if "cnv" not in self.exclude:
                            cnvsummary = ui.column().classes("space-y-1")
                        if "mgmt" not in self.exclude:
                            mgmt = ui.column().classes("space-y-1")
                        if "fusion" not in self.exclude:
                            fusions = ui.column().classes("space-y-1")

                    # Add monitoring information panel - Now positioned after diagnosis cards
                    # Only show in live mode
                    if not self.browse:
                        with ui.card().classes("w-full bg-white rounded-lg shadow-sm"):
                            with ui.expansion(
                                "Monitoring Information", icon="folder"
                            ).classes("w-full") as monitoring_expansion:
                                # Add summary header content
                                with monitoring_expansion.add_slot("header"):
                                    with ui.row().classes(
                                        "items-center justify-between w-full"
                                    ):
                                        # Left side - File counts
                                        with ui.row().classes("items-center gap-4"):
                                            with ui.row().classes("items-center gap-1"):
                                                ui.icon("description").classes(
                                                    "text-blue-600"
                                                )
                                                total_files = ui.label().classes(
                                                    "text-sm text-gray-600"
                                                )

                                                def update_total_files():
                                                    total = (
                                                        app.storage.general[
                                                            self.mainuuid
                                                        ]["bam_count"]["total_files"]
                                                        or 0
                                                    )
                                                    total_files.text = (
                                                        f"{total:,} files"
                                                    )

                                                app.storage.general[self.mainuuid][
                                                    "bam_count"
                                                ].on_change(update_total_files)
                                                update_total_files()

                                            with ui.row().classes("items-center gap-1"):
                                                ui.icon("check_circle").classes(
                                                    "text-green-600"
                                                )
                                                pass_count = ui.label().classes(
                                                    "text-sm text-gray-600"
                                                )

                                                def update_pass_count():
                                                    count = app.storage.general[
                                                        self.mainuuid
                                                    ]["samples"][self.sampleID][
                                                        "file_counters"
                                                    ][
                                                        "bam_passed"
                                                    ]
                                                    pass_count.text = (
                                                        f"{count:,} passed"
                                                        if count is not None
                                                        else "--"
                                                    )

                                                app.storage.general[self.mainuuid][
                                                    "samples"
                                                ][self.sampleID][
                                                    "file_counters"
                                                ].on_change(
                                                    update_pass_count
                                                )
                                                update_pass_count()

                                            with ui.row().classes("items-center gap-1"):
                                                ui.icon("error").classes("text-red-600")
                                                fail_count = ui.label().classes(
                                                    "text-sm text-gray-600"
                                                )

                                                def update_fail_count():
                                                    count = app.storage.general[
                                                        self.mainuuid
                                                    ]["samples"][self.sampleID][
                                                        "file_counters"
                                                    ][
                                                        "bam_failed"
                                                    ]
                                                    fail_count.text = (
                                                        f"{count:,} failed"
                                                        if count is not None
                                                        else "--"
                                                    )

                                                app.storage.general[self.mainuuid][
                                                    "samples"
                                                ][self.sampleID][
                                                    "file_counters"
                                                ].on_change(
                                                    update_fail_count
                                                )
                                                update_fail_count()

                                        # Right side - Title and Sample Name
                                        with ui.row().classes("items-center gap-4"):
                                            ui.label("Monitoring Information:").classes(
                                                "text-sm font-medium text-gray-700"
                                            )
                                            ui.label(self.sampleID).classes(
                                                "text-sm text-gray-600"
                                            )

                                # Existing content
                                with ui.column().classes("space-y-4 p-4"):
                                    # File Paths Section
                                    ui.label("File Paths").classes(
                                        "text-xl font-medium text-gray-900"
                                    )
                                    with ui.column().classes("space-y-4 ml-1"):
                                        # Monitoring Path
                                        with ui.row().classes("items-start gap-2"):
                                            ui.icon("folder_open").classes(
                                                "text-blue-600 shrink-0 mt-1"
                                            )
                                            with ui.column().classes(
                                                "flex-grow min-w-0"
                                            ):
                                                ui.label("Monitoring:").classes(
                                                    "text-gray-600"
                                                )
                                                ui.label(f"{self.watchfolder}").classes(
                                                    "break-all text-gray-900 font-mono"
                                                )

                                        # Output Path
                                        with ui.row().classes("items-start gap-2"):
                                            ui.icon("output").classes(
                                                "text-blue-600 shrink-0 mt-1"
                                            )
                                            with ui.column().classes(
                                                "flex-grow min-w-0"
                                            ):
                                                ui.label("Output:").classes(
                                                    "text-gray-600"
                                                )
                                                ui.label(f"{self.output}").classes(
                                                    "break-all text-gray-900 font-mono"
                                                )

                                    # Summary Statistics Section - Grid layout for BAM File Summary and Sequencing Statistics
                                    with ui.grid().classes("grid-cols-4 gap-6 mt-6"):
                                        # BAM File Summary Card - Takes 1 column (1/4)
                                        with ui.card().classes(
                                            "col-span-1 p-4 bg-gray-50 rounded-lg"
                                        ):
                                            ui.label("BAM File Summary").classes(
                                                "text-lg font-medium text-gray-900 mb-4"
                                            )
                                            with ui.column().classes("space-y-3"):
                                                # Total Files
                                                with ui.row().classes(
                                                    "items-center justify-between w-full"
                                                ):
                                                    with ui.row().classes(
                                                        "items-center gap-2"
                                                    ):
                                                        ui.icon("description").classes(
                                                            "text-blue-600 shrink-0"
                                                        )
                                                        ui.label(
                                                            "Total Files:"
                                                        ).classes("text-gray-600")
                                                    total_files = ui.label().classes(
                                                        "font-medium text-gray-900"
                                                    )

                                                    def update_total_files():
                                                        total = app.storage.general[
                                                            self.mainuuid
                                                        ]["bam_count"]["total_files"]
                                                        total_files.text = (
                                                            f"{total:,}"
                                                            if total > 0
                                                            else "--"
                                                        )

                                                    app.storage.general[self.mainuuid][
                                                        "bam_count"
                                                    ].on_change(update_total_files)
                                                    update_total_files()

                                                # Passed Files
                                                with ui.row().classes(
                                                    "items-center justify-between w-full"
                                                ):
                                                    with ui.row().classes(
                                                        "items-center gap-2"
                                                    ):
                                                        ui.icon("check_circle").classes(
                                                            "text-green-600 shrink-0"
                                                        )
                                                        ui.label("Passed:").classes(
                                                            "text-gray-600"
                                                        )
                                                    passed = ui.label().classes(
                                                        "font-medium text-gray-900"
                                                    )

                                                    def update_passed():
                                                        passed.text = (
                                                            f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_passed']:,}"
                                                            if app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ][
                                                                "bam_passed"
                                                            ]
                                                            is not None
                                                            else "--"
                                                        )

                                                    app.storage.general[self.mainuuid][
                                                        "samples"
                                                    ][self.sampleID][
                                                        "file_counters"
                                                    ].on_change(
                                                        update_passed
                                                    )
                                                    update_passed()

                                                # Failed Files
                                                with ui.row().classes(
                                                    "items-center justify-between w-full"
                                                ):
                                                    with ui.row().classes(
                                                        "items-center gap-2"
                                                    ):
                                                        ui.icon("error").classes(
                                                            "text-red-600 shrink-0"
                                                        )
                                                        ui.label("Failed:").classes(
                                                            "text-gray-600"
                                                        )
                                                    failed = ui.label().classes(
                                                        "font-medium text-gray-900"
                                                    )

                                                    def update_failed():
                                                        failed.text = (
                                                            f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_failed']:,}"
                                                            if app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ][
                                                                "bam_failed"
                                                            ]
                                                            is not None
                                                            else "--"
                                                        )

                                                    app.storage.general[self.mainuuid][
                                                        "samples"
                                                    ][self.sampleID][
                                                        "file_counters"
                                                    ].on_change(
                                                        update_failed
                                                    )
                                                    update_failed()

                                        # Sequencing Statistics Card - Takes 3 columns (3/4)
                                        with ui.card().classes(
                                            "col-span-3 p-4 bg-gray-50 rounded-lg"
                                        ):
                                            ui.label("Sequencing Statistics").classes(
                                                "text-lg font-medium text-gray-900 mb-4"
                                            )
                                            with ui.grid().classes("grid-cols-4 gap-4"):
                                                # Read Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label(
                                                        "Base Count by Read Type"
                                                    ).classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Total Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Total bases in reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            total_reads = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_total_reads():
                                                                mapped = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "mapped_count"
                                                                    ]
                                                                )
                                                                unmapped = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "unmapped_count"
                                                                    ]
                                                                )
                                                                total = (
                                                                    mapped + unmapped
                                                                    if mapped
                                                                    is not None
                                                                    and unmapped
                                                                    is not None
                                                                    else 0
                                                                )
                                                                total_reads.text = (
                                                                    f"{total:,} bp"
                                                                    if total > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_total_reads
                                                            )
                                                            update_total_reads()

                                                        # Mapped Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Bases in mapped reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mapped():
                                                                count = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "mapped_count"
                                                                    ]
                                                                )
                                                                mapped.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mapped
                                                            )
                                                            update_mapped()

                                                        # Unmapped Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Bases in unmapped reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_unmapped():
                                                                count = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "unmapped_count"
                                                                    ]
                                                                )
                                                                unmapped.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_unmapped
                                                            )
                                                            update_unmapped()

                                                # Read Quality Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label(
                                                        "Read Quality Statistics"
                                                    ).classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Total Mapped
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Total mapped reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            total_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_total_mapped():
                                                                pass_mapped = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_count"
                                                                ]
                                                                fail_mapped = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_count"
                                                                ]
                                                                total = (
                                                                    pass_mapped
                                                                    + fail_mapped
                                                                    if pass_mapped
                                                                    is not None
                                                                    and fail_mapped
                                                                    is not None
                                                                    else 0
                                                                )
                                                                total_mapped.text = (
                                                                    f"{total:,} reads"
                                                                    if total > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_total_mapped
                                                            )
                                                            update_total_mapped()

                                                        # High Quality Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "High quality reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            pass_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_pass_mapped():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_count"
                                                                ]
                                                                pass_mapped.text = (
                                                                    f"{count:,} reads"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_pass_mapped
                                                            )
                                                            update_pass_mapped()

                                                        # Low Quality Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Low quality reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            fail_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_fail_mapped():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_count"
                                                                ]
                                                                fail_mapped.text = (
                                                                    f"{count:,} reads"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_fail_mapped
                                                            )
                                                            update_fail_mapped()

                                                # Base Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label("Base Statistics").classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Total Bases
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Total bases:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            total_bases = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_total_bases():
                                                                pass_bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_bases_count"
                                                                ]
                                                                fail_bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_bases_count"
                                                                ]
                                                                total = (
                                                                    pass_bases
                                                                    + fail_bases
                                                                    if pass_bases
                                                                    is not None
                                                                    and fail_bases
                                                                    is not None
                                                                    else 0
                                                                )
                                                                total_bases.text = (
                                                                    f"{total:,} bp"
                                                                    if total > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_total_bases
                                                            )
                                                            update_total_bases()

                                                        # High Quality Bases
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "High quality bases:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            pass_bases = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_pass_bases():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_bases_count"
                                                                ]
                                                                pass_bases.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_pass_bases
                                                            )
                                                            update_pass_bases()

                                                        # Low Quality Bases
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Low quality bases:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            fail_bases = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_fail_bases():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_bases_count"
                                                                ]
                                                                fail_bases.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_fail_bases
                                                            )
                                                            update_fail_bases()

                                                # Read Length Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label(
                                                        "Read Length Statistics"
                                                    ).classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Overall Mean Lengths
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean mapped length:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_mapped():
                                                                bases = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "mapped_bases"
                                                                    ]
                                                                )
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "mapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_mapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_mapped
                                                            )
                                                            update_mean_mapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean unmapped length:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_unmapped():
                                                                bases = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "unmapped_bases"
                                                                    ]
                                                                )
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "unmapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_unmapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_unmapped
                                                            )
                                                            update_mean_unmapped()

                                                        # Pass/Fail Mean Lengths
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean pass mapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_pass_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_pass_mapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_pass_mapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_pass_mapped
                                                            )
                                                            update_mean_pass_mapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean fail mapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_fail_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_fail_mapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_fail_mapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_fail_mapped
                                                            )
                                                            update_mean_fail_mapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean pass unmapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_pass_unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_pass_unmapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_unmapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_unmapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_pass_unmapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_pass_unmapped
                                                            )
                                                            update_mean_pass_unmapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean fail unmapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_fail_unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_fail_unmapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_unmapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_unmapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_fail_unmapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_fail_unmapped
                                                            )
                                                            update_mean_fail_unmapped()

                # Detailed Analysis Tabs
                if sample_id:
                    selectedtab = None
                    with ui.tabs().classes("w-full") as tabs:
                        if not (
                            set(["sturgeon", "pannanodx", "nanodx", "forest"]).issubset(
                                set(self.exclude)
                            )
                        ):
                            methylationtab = ui.tab("Methylation Classification")
                            if not selectedtab:
                                selectedtab = methylationtab
                        if "cnv" not in self.exclude:
                            copy_numbertab = ui.tab("Copy Number Variation")
                            if not selectedtab:
                                selectedtab = copy_numbertab
                        if "coverage" not in self.exclude:
                            coveragetab = ui.tab("Target Coverage")
                            if not selectedtab:
                                selectedtab = coveragetab
                        if "mgmt" not in self.exclude:
                            mgmttab = ui.tab("MGMT")
                            if not selectedtab:
                                selectedtab = mgmttab
                        if "fusion" not in self.exclude:
                            fusionstab = ui.tab("Fusions")
                            if not selectedtab:
                                selectedtab = fusionstab

                    with ui.tab_panels(tabs, value=selectedtab).classes("w-full"):
                        display_args = {
                            "threads": self.threads,
                            "output": self.output,
                            "progress": True,
                            "browse": self.browse,
                            "bamqueue": None,
                            "uuid": self.mainuuid,
                            "force_sampleid": self.force_sampleid,
                            "sample_id": self.sampleID,
                        }
                        if not (
                            set(["sturgeon", "nanodx", "pannanodx", "forest"]).issubset(
                                set(self.exclude)
                            )
                        ):
                            with ui.tab_panel(methylationtab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    ui.label("Methylation Classifications").classes(
                                        "text-sky-600 dark:text-white"
                                    ).style(
                                        "font-size: 150%; font-weight: 300"
                                    ).tailwind(
                                        "drop-shadow", "font-bold"
                                    )
                                    if "sturgeon" not in self.exclude:
                                        self.Sturgeon = Sturgeon_object(
                                            analysis_name="STURGEON",
                                            batch=True,
                                            summary=sturgeonsummary,
                                            **display_args,
                                        )
                                        await self.Sturgeon.render_ui(
                                            sample_id=self.sampleID
                                        )
                                    if "nanodx" not in self.exclude:
                                        self.NanoDX = NanoDX_object(
                                            analysis_name="NANODX",
                                            batch=True,
                                            summary=nanodxsummary,
                                            **display_args,
                                        )
                                        await self.NanoDX.render_ui(
                                            sample_id=self.sampleID
                                        )
                                    if "pannanodx" not in self.exclude:
                                        self.PanNanoDX = NanoDX_object(
                                            analysis_name="PANNANODX",
                                            batch=True,
                                            summary=pannanodxsummary,
                                            model="pancan_devel_v5i_NN.pkl",
                                            **display_args,
                                        )
                                        await self.PanNanoDX.render_ui(
                                            sample_id=self.sampleID
                                        )
                                    if "forest" not in self.exclude:
                                        self.RandomForest = RandomForest_object(
                                            analysis_name="FOREST",
                                            batch=True,
                                            summary=forestsummary,
                                            showerrors=self.showerrors,
                                            **display_args,
                                        )
                                        await self.RandomForest.render_ui(
                                            sample_id=self.sampleID
                                        )

                        if "cnv" not in self.exclude:
                            with ui.tab_panel(copy_numbertab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.CNV = CNVAnalysis(
                                        analysis_name="CNV",
                                        summary=cnvsummary,
                                        target_panel=self.target_panel,
                                        reference_file=self.reference,
                                        bed_file=self.bed_file,
                                        NewBed=self.NewBed,
                                        **display_args,
                                    )
                                    await self.CNV.render_ui(sample_id=self.sampleID)

                        if "coverage" not in self.exclude:
                            with ui.tab_panel(coveragetab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.Target_Coverage = TargetCoverage(
                                        analysis_name="COVERAGE",
                                        summary=coverage,
                                        target_panel=self.target_panel,
                                        reference=self.reference,
                                        **display_args,
                                    )
                                    await self.Target_Coverage.render_ui(
                                        sample_id=self.sampleID
                                    )

                        if "mgmt" not in self.exclude:
                            with ui.tab_panel(mgmttab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.MGMT_panel = MGMT_Object(
                                        analysis_name="MGMT",
                                        summary=mgmt,
                                        **display_args,
                                    )
                                    await self.MGMT_panel.render_ui(
                                        sample_id=self.sampleID
                                    )

                        if "fusion" not in self.exclude:
                            with ui.tab_panel(fusionstab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.Fusion_panel = FusionObject(
                                        analysis_name="FUSION",
                                        summary=fusions,
                                        target_panel=self.target_panel,
                                        NewBed=self.NewBed,
                                        **display_args,
                                    )
                                    await self.Fusion_panel.render_ui(
                                        sample_id=self.sampleID
                                    )

                    async def confirm_report_generation():
                        """
                        Show a confirmation dialog before generating the report.
                        """
                        report_types = {
                            "summary": "Summary Only",
                            "detailed": "Detailed",
                        }
                        state = {"type": "detailed"}  # Store state in a dictionary

                        with ui.dialog() as dialog, ui.card().classes("w-96 p-4"):
                            ui.label("Generate Report").classes(
                                "text-h6 font-bold mb-4"
                            )

                            # Report type selector
                            with ui.column().classes("mb-4"):
                                ui.label("Report Type").classes("font-bold mb-2")
                                ui.toggle(
                                    report_types,
                                    value="detailed",
                                    on_change=lambda e: state.update({"type": e.value}),
                                )

                            # Disclaimer section
                            with ui.column().classes("mb-4"):
                                ui.label("Disclaimer").classes("font-bold mb-2")
                                # Split into paragraphs and join with double breaks for better spacing
                                formatted_text = EXTENDED_DISCLAIMER_TEXT.replace(
                                    "\n\n", "<br><br>"
                                ).replace("\n", " ")
                                ui.label(formatted_text).classes(
                                    "text-sm text-gray-600 mb-4"
                                )

                            ui.label(
                                "Are you sure you want to generate a report?"
                            ).classes("mb-4")

                            # Buttons in a row at the bottom
                            with ui.row().classes("justify-end gap-2"):
                                ui.button(
                                    "No", on_click=lambda: dialog.submit(("No", None))
                                ).props("flat")
                                ui.button(
                                    "Yes",
                                    on_click=lambda: dialog.submit(
                                        ("Yes", state["type"])
                                    ),
                                ).props("color=primary")

                        result, report_type = await dialog
                        if result == "Yes":
                            await download_report(report_type)

                    async def download_report(report_type="detailed"):
                        """
                        Generate and download the report.

                        :param report_type: Type of report to generate ('summary' or 'detailed')
                        :return: None
                        """
                        ui.notify("Generating Report")
                        if not self.browse:
                            for item in app.storage.general[self.mainuuid]:
                                if item == "sample_ids":
                                    for sample in app.storage.general[self.mainuuid][
                                        item
                                    ]:
                                        self.sampleID = sample
                        if self.browse:
                            myfile = await run.io_bound(
                                create_pdf,
                                f"{self.sampleID}_run_report.pdf",
                                self.check_and_create_folder(
                                    self.output, self.sampleID
                                ),
                                report_type,
                            )
                        else:
                            myfile = await run.io_bound(
                                create_pdf,
                                f"{self.sampleID}_run_report.pdf",
                                self.output,
                                report_type,
                            )
                        ui.download(myfile)
                        ui.notify("Report Downloaded")

                    # Replace the direct download_report button with the confirmation dialog
                    ui.button(
                        "Generate Report",
                        on_click=confirm_report_generation,
                        icon="download",
                    )

        except Exception as e:
            logging.error(f"Error rendering analysis UI: {str(e)}", exc_info=True)
            ui.notify(f"Error rendering analysis UI: {str(e)}", type="error")

    async def background_process_bams(self):
        """
        Background task to process BAM files.

        :return: None
        """
        logging.info("Starting background BAM processing")
        self.process_bams_tracker = app.timer(10, self.process_bams)
        self.process_bigbadmerge_tracker = app.timer(10, self.process_bigbadmerge)

    def check_and_create_folder(self, path, folder_name=None):
        """
        Check if a folder exists and create it if necessary.

        :param path: Base path.
        :param folder_name: Folder name to create within the base path.
        :return: Full path to the folder.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"The specified path does not exist: {path}")

        if self.force_sampleid:
            folder_name = self.force_sampleid

        if folder_name:
            # Check if the path already ends with the folder_name
            if path.endswith(folder_name):
                return path
            full_path = os.path.join(path, folder_name)
            if not os.path.exists(full_path):
                os.makedirs(full_path)
                logging.info(f"Folder created: {full_path}")
            return full_path
        else:
            return path

    async def process_bigbadmerge(self):
        # Dictionary to store files by sample ID
        files_by_sample = {}
        latest_files = {}  # Track latest file time per sample

        # Process queue and organize files by sample ID
        while self.bamforbigbadmerge.qsize() > 0:
            file, filetime, sampleID = self.bamforbigbadmerge.get()

            # Initialize containers for new sample IDs
            if sampleID not in files_by_sample:
                files_by_sample[sampleID] = []
                latest_files[sampleID] = 0

                # Create temporary directories if needed
                if sampleID not in self.dataDir:
                    self.dataDir[sampleID] = tempfile.TemporaryDirectory(
                        dir=self.check_and_create_folder(self.output, sampleID)
                    )
                    self.bedDir[sampleID] = tempfile.TemporaryDirectory(
                        dir=self.check_and_create_folder(self.output, sampleID)
                    )

                # Initialize storage for this sample if it doesn't exist
                if sampleID not in app.storage.general[self.mainuuid]:
                    app.storage.general[self.mainuuid][sampleID] = {}

                # Initialize counters for this sample if they don't exist
                analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
                for analysis in analyses:
                    if analysis.lower() not in self.exclude:
                        if analysis not in app.storage.general[self.mainuuid][sampleID]:
                            app.storage.general[self.mainuuid][sampleID][analysis] = {
                                "counters": {
                                    "bams_in_processing": 0,
                                    "bams_processed": 0,
                                    "bams_failed": 0,
                                }
                            }
                        elif (
                            "counters"
                            not in app.storage.general[self.mainuuid][sampleID][
                                analysis
                            ]
                        ):
                            app.storage.general[self.mainuuid][sampleID][analysis][
                                "counters"
                            ] = {
                                "bams_in_processing": 0,
                                "bams_processed": 0,
                                "bams_failed": 0,
                            }

            # Update latest file time for this sample
            if filetime > latest_files[sampleID]:
                latest_files[sampleID] = filetime

            files_by_sample[sampleID].append(file)

            # Process if we have enough files for any sample
            for sample_id in list(files_by_sample.keys()):
                if len(files_by_sample[sample_id]) >= 50:
                    files_to_process = len(files_by_sample[sample_id])
                    logging.info(
                        f"Processing batch of {files_to_process} files for sample {sample_id}"
                    )

                    # Update counters for each non-excluded analysis
                    analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
                    for analysis in analyses:
                        if analysis.lower() not in self.exclude:
                            try:
                                app.storage.general[self.mainuuid][sample_id][analysis][
                                    "counters"
                                ]["bams_in_processing"] += files_to_process
                                logging.debug(
                                    f"Updated {analysis} counter for {sample_id} by {files_to_process}"
                                )
                            except KeyError as e:
                                logging.error(
                                    f"Failed to update counter for {analysis}: {str(e)}"
                                )

                    await self.process_sample_files(
                        sample_id, files_by_sample[sample_id], latest_files[sample_id]
                    )

                    # Clear processed files
                    files_by_sample[sample_id] = []
                    latest_files[sample_id] = 0

        # Process remaining files for each sample
        for sample_id, files in files_by_sample.items():
            if files:  # Only process if there are files
                files_to_process = len(files)
                if files_to_process > 0:
                    logging.info(
                        f"Processing remaining {files_to_process} files for sample {sample_id}"
                    )

                    # Update counters for each non-excluded analysis
                    analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
                    for analysis in analyses:
                        if analysis.lower() not in self.exclude:
                            try:
                                app.storage.general[self.mainuuid][sample_id][analysis][
                                    "counters"
                                ]["bams_in_processing"] += files_to_process
                                logging.debug(
                                    f"Updated {analysis} counter for {sample_id} by {files_to_process}"
                                )
                            except KeyError as e:
                                logging.error(
                                    f"Failed to update counter for {analysis}: {str(e)}"
                                )

                    await self.process_sample_files(
                        sample_id, files, latest_files[sample_id]
                    )

    async def process_sample_files(self, sampleID, tomerge, latest_file):
        """
        Process a batch of BAM files for a specific sample.

        Args:
            sampleID (str): The sample ID being processed
            tomerge (list): List of BAM files to merge
            latest_file (float): Timestamp of the latest file
        """
        if not tomerge:  # Skip if no files to process
            return

        try:
            # Track the number of BAM files seen and merged
            num_bam_files_seen = len(tomerge)
            logging.info(
                f"Processing {num_bam_files_seen} BAM files for sample ID: {sampleID}"
            )

            # Ensure counters exist for this sample
            analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
            for analysis in analyses:
                if analysis.lower() not in self.exclude:
                    if analysis not in app.storage.general[self.mainuuid][sampleID]:
                        app.storage.general[self.mainuuid][sampleID][analysis] = {
                            "counters": {
                                "bams_in_processing": 0,
                                "bams_processed": 0,
                                "bams_failed": 0,
                            }
                        }
                    elif (
                        "counters"
                        not in app.storage.general[self.mainuuid][sampleID][analysis]
                    ):
                        app.storage.general[self.mainuuid][sampleID][analysis][
                            "counters"
                        ] = {
                            "bams_in_processing": 0,
                            "bams_processed": 0,
                            "bams_failed": 0,
                        }

            # Write the length of the tomerge list to the output file FIRST
            tomerge_length_file = os.path.join(
                self.check_and_create_folder(self.output, sampleID),
                "tomerge_length.txt",
            )

            # Initialize the count
            if os.path.exists(tomerge_length_file):
                with open(tomerge_length_file, "r") as f:
                    current_count = int(
                        f.readline().strip().split(": ")[1]
                    )  # Read the current count
            else:
                current_count = 0  # If the file doesn't exist, start from 0

            # Update the count
            new_count = current_count + len(tomerge)

            # Write the updated length of the tomerge list to the output file
            with open(tomerge_length_file, "w") as f:
                f.write(f"Length of tomerge list: {new_count}\n")

            tempbam = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID),
                suffix=".bam",
            )
            sorttempbam = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID),
                suffix=".bam",
            )
            file = tempbam.name
            temp = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
            sortfile = sorttempbam.name

            # Sort and merge BAM files
            await run.cpu_bound(
                run_samtools_sort, file, tomerge, sortfile, self.threads
            )

            with tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            ) as temp2:
                await run.cpu_bound(run_modkit, sortfile, temp.name, self.threads)

                # Create output path specific to this sample
                parquet_path = os.path.join(
                    self.check_and_create_folder(self.output, sampleID),
                    f"{sampleID}.parquet",  # Use sampleID for the output filename
                )

                # Merge modkit files for this sample
                merged_df = await run.cpu_bound(
                    merge_modkit_files,
                    [temp.name],
                    parquet_path,  # Use the parquet_path for the existing file
                    parquet_path,
                    self.cpgs_master_file,  # Use the parquet_path for the output file
                )

                # Log the number of BAM files processed
                logging.info(
                    f"Merged {num_bam_files_seen} BAM files into {parquet_path} for sample ID: {sampleID}"
                )
                self.mergecounter += len(tomerge)

                # Update the processed counters for each analysis type
                analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
                for analysis in analyses:
                    if analysis.lower() not in self.exclude:
                        try:
                            counters = app.storage.general[self.mainuuid][sampleID][
                                analysis
                            ]["counters"]
                            if "bams_in_processing" not in counters:
                                counters["bams_in_processing"] = 0
                            if "bams_processed" not in counters:
                                counters["bams_processed"] = 0

                            # Decrease the in_processing counter and increase the processed counter
                            counters["bams_in_processing"] = max(
                                0, counters["bams_in_processing"] - num_bam_files_seen
                            )
                            counters["bams_processed"] += num_bam_files_seen

                            logging.debug(
                                f"Updated {analysis} processed counter for {sampleID} by {num_bam_files_seen}"
                            )
                        except KeyError as e:
                            logging.error(
                                f"Failed to update processed counter for {analysis}: {str(e)}"
                            )
                            # Initialize counters if they don't exist
                            app.storage.general[self.mainuuid][sampleID][analysis] = {
                                "counters": {
                                    "bams_in_processing": 0,
                                    "bams_processed": num_bam_files_seen,
                                    "bams_failed": 0,
                                }
                            }
                        except Exception as e:
                            logging.error(
                                f"Error updating counters for {analysis}: {str(e)}"
                            )

        except Exception as e:
            logging.error(f"Error in process_sample_files: {str(e)}", exc_info=True)
            # Update the failed counters for each analysis type
            analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
            for analysis in analyses:
                if analysis.lower() not in self.exclude:
                    try:
                        app.storage.general[self.mainuuid][sampleID][analysis][
                            "counters"
                        ]["bams_in_processing"] -= num_bam_files_seen
                        app.storage.general[self.mainuuid][sampleID][analysis][
                            "counters"
                        ]["bams_failed"] += num_bam_files_seen
                        logging.debug(
                            f"Updated {analysis} failed counter for {sampleID} by {num_bam_files_seen}"
                        )
                    except Exception as nested_e:
                        logging.error(
                            f"Error updating failed counters for {analysis}: {str(nested_e)}"
                        )
            raise

    async def process_bams(self) -> None:
        """
        Process BAM files and update the application's storage.

        :return: None
        """
        logging.info("Process Bam Starting")
        self.process_bams_tracker.active = False
        counter = 0
        if "file" in app.storage.general[self.mainuuid]["bam_count"]:
            while self.watchdogbamqueue.qsize() > 0:
                filename, timestamp = self.watchdogbamqueue.get()
                app.storage.general[self.mainuuid]["bam_count"]["counter"] += 1
                app.storage.general[self.mainuuid]["bam_count"][
                    "total_files"
                ] += 1  # Increment total files
                app.storage.general[self.mainuuid]["bam_count"]["file"][
                    filename
                ] = timestamp

            while len(app.storage.general[self.mainuuid]["bam_count"]["file"]) > 0:
                self.nofiles = False
                file = (
                    k := next(
                        iter(app.storage.general[self.mainuuid]["bam_count"]["file"])
                    ),
                    app.storage.general[self.mainuuid]["bam_count"]["file"].pop(k),
                )
                baminfo, bamdata = await run.cpu_bound(check_bam, file[0])
                sample_id = baminfo["sample_id"]
                if sample_id not in app.storage.general[self.mainuuid]["samples"]:
                    app.storage.general[self.mainuuid]["samples"][sample_id] = {}
                    self.configure_storage(sample_id)
                    app.storage.general[self.mainuuid]["sample_list"].append(sample_id)
                if baminfo["state"] == "pass":
                    self.update_counter(sample_id, "bam_passed", 1)
                    self.update_counter(
                        sample_id, "pass_mapped_count", bamdata["mapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "pass_unmapped_count", bamdata["unmapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "pass_bases_count", bamdata["yield_tracking"]
                    )
                    # Add pass read numbers and bases
                    self.update_counter(
                        sample_id, "pass_mapped_reads_num", bamdata["mapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id,
                        "pass_unmapped_reads_num",
                        bamdata["unmapped_reads_num"],
                    )
                    self.update_counter(
                        sample_id, "pass_mapped_bases", bamdata["mapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "pass_unmapped_bases", bamdata["unmapped_bases"]
                    )
                else:
                    self.update_counter(sample_id, "bam_failed", 1)
                    self.update_counter(
                        sample_id, "fail_mapped_count", bamdata["mapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "fail_unmapped_count", bamdata["unmapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "fail_bases_count", bamdata["yield_tracking"]
                    )
                    # Add fail read numbers and bases
                    self.update_counter(
                        sample_id, "fail_mapped_reads_num", bamdata["mapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id,
                        "fail_unmapped_reads_num",
                        bamdata["unmapped_reads_num"],
                    )
                    self.update_counter(
                        sample_id, "fail_mapped_bases", bamdata["mapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "fail_unmapped_bases", bamdata["unmapped_bases"]
                    )

                # Update total counts
                self.update_counter(
                    sample_id,
                    "mapped_count",
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_mapped_count"]
                    + app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_mapped_count"],
                )
                self.update_counter(
                    sample_id,
                    "unmapped_count",
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_unmapped_count"]
                    + app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_unmapped_count"],
                )
                self.update_counter(
                    sample_id,
                    "bases_count",
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_bases_count"]
                    + app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_bases_count"],
                )

                # Update total read numbers and bases
                self.update_counter(
                    sample_id,
                    "mapped_reads_num",
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_mapped_reads_num"]
                    + app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_mapped_reads_num"],
                )
                self.update_counter(
                    sample_id,
                    "unmapped_reads_num",
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_unmapped_reads_num"]
                    + app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_unmapped_reads_num"],
                )
                self.update_counter(
                    sample_id,
                    "mapped_bases",
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_mapped_bases"]
                    + app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_mapped_bases"],
                )
                self.update_counter(
                    sample_id,
                    "unmapped_bases",
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_unmapped_bases"]
                    + app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_unmapped_bases"],
                )

                if (
                    baminfo["device_position"]
                    not in app.storage.general[self.mainuuid]["samples"][sample_id][
                        "devices"
                    ]
                ):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "devices"
                    ].append(baminfo["device_position"])
                if (
                    baminfo["basecall_model"]
                    not in app.storage.general[self.mainuuid]["samples"][sample_id][
                        "basecall_models"
                    ]
                ):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "basecall_models"
                    ].append(baminfo["basecall_model"])
                if not self.force_sampleid:
                    if (
                        baminfo["sample_id"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "sample_ids"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "sample_ids"
                        ].append(baminfo["sample_id"])
                else:
                    if (
                        self.force_sampleid
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "sample_ids"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "sample_ids"
                        ].append(self.force_sampleid)
                if (
                    baminfo["flow_cell_id"]
                    not in app.storage.general[self.mainuuid]["samples"][sample_id][
                        "flowcell_ids"
                    ]
                ):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "flowcell_ids"
                    ].append(baminfo["flow_cell_id"])
                if (
                    baminfo["time_of_run"]
                    not in app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_time"
                    ]
                ):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_time"
                    ].append(baminfo["time_of_run"])

                mydf = pd.DataFrame.from_dict(app.storage.general)
                if not self.force_sampleid:
                    sample_id = baminfo["sample_id"]
                else:
                    sample_id = self.force_sampleid
                mydf.to_csv(
                    os.path.join(
                        self.check_and_create_folder(self.output, sample_id),
                        "master.csv",
                    )
                )

                counter += 1
                analyses = ["forest", "sturgeon", "nanodx", "pannanodx"]
                if any(analysis.lower() not in self.exclude for analysis in analyses):
                    self.bamforbigbadmerge.put([file[0], file[1], sample_id])

                if "forest" not in self.exclude:
                    self.bamforcns.put([file[0], file[1], sample_id])
                if "sturgeon" not in self.exclude:
                    self.bamforsturgeon.put([file[0], file[1], sample_id])
                if "nanodx" not in self.exclude:
                    self.bamfornanodx.put([file[0], file[1], sample_id])
                if "pannanodx" not in self.exclude:
                    self.bamforpannanodx.put([file[0], file[1], sample_id])
                if "cnv" not in self.exclude:
                    self.bamforcnv.put([file[0], file[1], sample_id])
                if "coverage" not in self.exclude:
                    self.bamfortargetcoverage.put([file[0], file[1], sample_id])
                if "mgmt" not in self.exclude:
                    self.bamformgmt.put([file[0], file[1], sample_id])
                if "fusion" not in self.exclude:
                    self.bamforfusions.put([file[0], file[1], sample_id])

            self.nofiles = True
        logging.info("Process Bam Finishing")
        self.process_bams_tracker.active = True

    async def check_existing_bams(self, sequencing_summary=None):
        """
        Check and process existing BAM files based on a sequencing summary.

        :param sequencing_summary: Path to the sequencing summary file.
        :return: None
        """
        file_endings = {".bam"}
        """
        if sequencing_summary:
            logging.info("Using sequencing summary for existing BAM files")
            df = pd.read_csv(
                sequencing_summary,
                delimiter="\t",
                usecols=["filename_bam", "template_start", "template_duration"],
            )
            df["template_end"] = df["template_start"] + df["template_duration"]

            df.drop(columns=["template_start", "template_duration"], inplace=True)
            latest_timestamps = (
                df.groupby("filename_bam")["template_end"]
                .max()
                .reset_index()
                .sort_values(by="template_end")
                .reset_index(drop=True)
            )

            latest_timestamps["full_path"] = ""
            latest_timestamps["file_produced"] = latest_timestamps["template_end"]
            for path, dirs, files in os.walk(self.watchfolder):
                for f in files:
                    if "".join(Path(f).suffixes) in file_endings:
                        latest_timestamps.loc[
                            latest_timestamps["filename_bam"] == f, "full_path"
                        ] = os.path.join(path, f)

            step_size = 20

            if "sturgeon" not in self.exclude:
                self.Sturgeon.playback(latest_timestamps, step_size=step_size)
            if "nanodx" not in self.exclude:
                self.NanoDX.playback(latest_timestamps, step_size=step_size)
            if "forest" not in self.exclude:
                self.RandomForest.playback(latest_timestamps, step_size=step_size)
            if "cnv" not in self.exclude:
                self.CNV.playback(latest_timestamps, step_size=step_size)
            if "coverage" not in self.exclude:
                self.Target_Coverage.playback(
                    latest_timestamps, step_size=step_size, start_time=self.run_time
                )
            if "fusion" not in self.exclude:
                self.Fusion_panel.playback(latest_timestamps, step_size=step_size)
            if "mgmt" not in self.exclude:
                self.MGMT_panel.playback(latest_timestamps, step_size=step_size)

            for index, row in latest_timestamps.iterrows():
                if len(row["full_path"]) > 0:
                    self.check_bam(row["full_path"])
            self.runfinished = True
        else:
        """
        files_and_timestamps = []
        files_and_timestamps = await run.cpu_bound(
            sort_bams,
            files_and_timestamps,
            self.watchfolder,
            file_endings,
            self.simtime,
        )

        # Iterate through the sorted list
        for timestamp, f, elapsed_time in files_and_timestamps:
            logging.debug(f"Processing existing BAM file: {f} at {timestamp}")
            app.storage.general[self.mainuuid]["bam_count"]["counter"] += 1
            app.storage.general[self.mainuuid]["bam_count"][
                "total_files"
            ] += 1  # Increment total files counter
            if "file" not in app.storage.general[self.mainuuid]["bam_count"]:
                app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
            app.storage.general[self.mainuuid]["bam_count"]["file"][
                f
            ] = timestamp.timestamp()
            if self.simtime:
                await asyncio.sleep(1)
            else:
                await asyncio.sleep(0)
