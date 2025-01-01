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
from alive_progress import alive_bar
from nicegui import ui, app, run


from robin.utilities.bam_handler import BamEventHandler
from robin.subpages.MGMT_object import MGMT_Object
from robin.subpages.Sturgeon_object import Sturgeon_object
from robin.subpages.NanoDX_object import NanoDX_object
from robin.subpages.RandomForest_object import RandomForest_object
from robin.subpages.CNV_object import CNVAnalysis
from robin.subpages.TargetCoverage_object import TargetCoverage
from robin.subpages.Fusion_object import FusionObject
from robin.utilities.local_file_picker import LocalFilePicker
from robin.utilities.ReadBam import ReadBam
from robin.reporting.report import create_pdf

from watchdog.observers import Observer


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
    import bisect

    # Function to insert into sorted list
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
                            (
                                os.path.join(path, f),
                                filetime.isoformat(),
                                elapsedtime,
                            )
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

        logging.info(f"BrainMeth initialized with UUID: {self.mainuuid}")

    def configure_storage(self, sample_id):
        """
        Configure storage for the application.

        :param sample_id: Sample ID to configure storage for.
        """
        logging.info(f"Configuring storage for sample ID: {sample_id}")
        app.storage.general[self.mainuuid]["samples"][sample_id]["file_counters"] = (
            Counter(
                bam_passed=0,
                bam_failed=0,
                mapped_count=0,
                pass_mapped_count=0,
                fail_mapped_count=0,
                unmapped_count=0,
                pass_unmapped_count=0,
                fail_unmapped_count=0,
                pass_bases_count=0,
                fail_bases_count=0,
                bases_count=0,
            )
        )
        app.storage.general[self.mainuuid]["samples"][sample_id]["devices"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["basecall_models"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["run_time"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["flowcell_ids"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["sample_ids"] = []

    async def start_background(self):
        """
        Start background tasks for BAM processing and analysis.
        """
        logging.info(f"Starting background tasks for UUID: {self.mainuuid}")
        app.storage.general[self.mainuuid]["bam_count"] = {}
        app.storage.general[self.mainuuid]["bam_count"] = Counter(counter=0)
        app.storage.general[self.mainuuid]["bam_count"]["files"] = {}
        self.bam_tracking = Queue()
        self.bamforcns = Queue()
        self.bamforsturgeon = Queue()
        self.bamfornanodx = Queue()
        self.bamforpannanodx = Queue()
        self.bamforcnv = Queue()
        self.bamfortargetcoverage = Queue()
        self.bamformgmt = Queue()
        self.bamforfusions = Queue()

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
                logging.info(f"Split path - output: {self.output}, sampleID: {self.sampleID}")
                
                # Initialize storage for browse mode
                if self.sampleID not in app.storage.general[self.mainuuid]["samples"]:
                    logging.info(f"Initializing storage for sample {self.sampleID}")
                    app.storage.general[self.mainuuid]["samples"][self.sampleID] = {}
                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"] = {
                        "bam_passed": 0,
                        "bam_failed": 0,
                        "mapped_count": 0,
                        "unmapped_count": 0,
                        "pass_mapped_count": 0,
                        "fail_mapped_count": 0,
                        "pass_bases_count": 0,
                        "fail_bases_count": 0,
                        "bases_count": 0
                    }
                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["devices"] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["basecall_models"] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["run_time"] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["flowcell_ids"] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["sample_ids"] = [self.sampleID]
                else:
                    logging.info(f"Storage already exists for sample {self.sampleID}")
                
                if self.sampleID not in app.storage.general[self.mainuuid]["sample_list"]:
                    logging.info(f"Adding {self.sampleID} to sample list")
                    app.storage.general[self.mainuuid]["sample_list"].append(self.sampleID)

                # Try to load existing data from master.csv if it exists
                master_csv = os.path.join(result[0], "master.csv")
                logging.info(f"Looking for master.csv at: {master_csv}")
                if os.path.exists(master_csv):
                    try:
                        logging.info("Found master.csv, loading data")
                        df = pd.read_csv(master_csv)
                        # Update storage with data from master.csv
                        for key in df.keys():
                            if key in app.storage.general[self.mainuuid]["samples"][self.sampleID]:
                                logging.info(f"Updating {key} from master.csv")
                                app.storage.general[self.mainuuid]["samples"][self.sampleID][key] = df[key].tolist()
                    except Exception as e:
                        logging.error(f"Error loading master.csv: {str(e)}", exc_info=True)
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
                                    ).style(
                                        "font-size: 75%; font-weight: 100"
                                    ).classes(
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
                    ui.label("Please select a sample folder to view").classes('text-xl text-gray-600 my-4')
                    return
                
            if not self.browse:
                logging.info("Live mode - showing sample list")
                self.show_list()
                app.storage.general[self.mainuuid]["sample_list"].on_change(
                    self.show_list.refresh
                )

            if self.sampleID and not self.browse and self.sampleID not in app.storage.general[self.mainuuid]["samples"]:
                logging.error(f"Sample {self.sampleID} not found in storage")
                ui.notify(f"Sample {self.sampleID} not found")
                ui.navigate.to("/live")
                return

            if not self.sampleID:
                logging.warning("No sample ID set, cannot render panel")
                return

            logging.info("Starting UI rendering")
            if self.sampleID in app.storage.general[self.mainuuid]["samples"]:
                logging.info(f"Storage state for sample {self.sampleID}: {app.storage.general[self.mainuuid]['samples'][self.sampleID]}")
            else:
                logging.warning(f"No storage found for sample {self.sampleID}")
                return

            with ui.column().classes('w-full px-6'):
                logging.info("Rendering title section")
                # Title and description
                with ui.column().classes('space-y-2 mb-4'):
                    ui.label("CNS Tumor Methylation Classification").classes('text-2xl font-medium text-gray-900')
                    ui.label("This tool enables classification of brain tumors in real time from Oxford Nanopore Data.").classes('text-gray-600')
                logging.info("Title section rendered")

                # Run Information - Compact display at the top
                if self.sampleID and self.sampleID in app.storage.general[self.mainuuid]["samples"]:
                    logging.info("Rendering run information display")
                    with ui.row().classes('w-full py-3 text-sm text-gray-600 items-center justify-between bg-gray-50 rounded-lg px-4 mb-4'):
                        if self.browse:
                            # In browse mode, only show sample ID
                            with ui.row().classes('items-center gap-2'):
                                ui.icon('label').classes('text-gray-400')
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][self.sampleID],
                                    "sample_ids",
                                    backward=lambda n: f"Sample: {str(n[0])}" if n else "--"
                                )
                        else:
                            # In live mode, show all information
                            # Left column - Device and Run Info
                            with ui.column().classes('space-y-1'):
                                with ui.row().classes('items-center gap-2'):
                                    ui.icon('devices').classes('text-gray-400')
                                    ui.label().bind_text_from(
                                        app.storage.general[self.mainuuid]["samples"][self.sampleID],
                                        "devices",
                                        backward=lambda n: f"Device: {str(n[0])}" if n else "--"
                                    )
                                with ui.row().classes('items-center gap-2'):
                                    ui.icon('memory').classes('text-gray-400')
                                    ui.label().bind_text_from(
                                        app.storage.general[self.mainuuid]["samples"][self.sampleID],
                                        "basecall_models",
                                        backward=lambda n: f"Model: {str(n[0])}" if n else "--"
                                    )
                            
                            # Middle column - Flow Cell and Sample
                            with ui.column().classes('space-y-1'):
                                with ui.row().classes('items-center gap-2'):
                                    ui.icon('view_module').classes('text-gray-400')
                                    ui.label().bind_text_from(
                                        app.storage.general[self.mainuuid]["samples"][self.sampleID],
                                        "flowcell_ids",
                                        backward=lambda n: f"Flow Cell: {str(n[0])}" if n else "--"
                                    )
                                with ui.row().classes('items-center gap-2'):
                                    ui.icon('label').classes('text-gray-400')
                                    ui.label().bind_text_from(
                                        app.storage.general[self.mainuuid]["samples"][self.sampleID],
                                        "sample_ids",
                                        backward=lambda n: f"Sample: {str(n[0])}" if n else "--"
                                    )
                            
                            # Right column - Time and Stats
                            with ui.column().classes('space-y-1'):
                                with ui.row().classes('items-center gap-2'):
                                    ui.icon('schedule').classes('text-gray-400')
                                    ui.label().bind_text_from(
                                        app.storage.general[self.mainuuid]["samples"][self.sampleID],
                                        "run_time",
                                        backward=lambda n: f"Run: {parser.parse(n[0]).strftime('%Y-%m-%d %H:%M')}" if n else "--"
                                    )
                    logging.info("Run information display rendered")

                # Results Summary Section
                with ui.column().classes('space-y-4'):
                    # Classification Results - Single row with equal width columns
                    with ui.grid().classes('grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4'):
                        if "sturgeon" not in self.exclude:
                            sturgeonsummary = ui.column().classes('space-y-1')
                        if "nanodx" not in self.exclude:
                            nanodxsummary = ui.column().classes('space-y-1')
                        if "pannanodx" not in self.exclude:
                            pannanodxsummary = ui.column().classes('space-y-1')
                        if "forest" not in self.exclude:
                            forestsummary = ui.column().classes('space-y-1')

                    # Analysis Results - Two columns grid
                    with ui.grid().classes('grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4'):
                        if "coverage" not in self.exclude:
                            coverage = ui.column().classes('space-y-1')
                        if "cnv" not in self.exclude:
                            cnvsummary = ui.column().classes('space-y-1')
                        if "mgmt" not in self.exclude:
                            mgmt = ui.column().classes('space-y-1')
                        if "fusion" not in self.exclude:
                            fusions = ui.column().classes('space-y-1')

                    # Add monitoring information panel - Now positioned after diagnosis cards
                    # Only show in live mode
                    if not self.browse:
                        with ui.expansion("Monitoring Information", icon="folder"): #.classes('mt-4 sm:mt-5 md:mt-6 mb-2 sm:mb-3 md:mb-4'):
                            with ui.column():#.classes('p-2 sm:p-3 md:p-4 space-y-2 sm:space-y-3 md:space-y-4 bg-white rounded-lg shadow-sm'):
                                # Paths Section
                                with ui.column(): #.classes('space-y-3 pb-4 border-b border-gray-200'):
                                    ui.label("File Paths")#.classes('break-all font-medium text-gray-900')
                                    # Monitoring Path
                                    with ui.column():#.classes('space-y-2'):
                                        with ui.row():#.classes('items-start gap-2 flex-wrap'):
                                            ui.icon('folder_open')#.classes('text-blue-600 shrink-0 mt-1')
                                            with ui.column():#classes('flex-grow min-w-0'):
                                                ui.label("Monitoring:").classes('text-gray-600')
                                                ui.label(f"{self.watchfolder}").classes('break-all text-gray-900 font-mono')
                                        
                                        # Output Path
                                        with ui.row().classes('items-start gap-2 flex-wrap'):
                                            ui.icon('output').classes('text-blue-600 shrink-0 mt-1')
                                            with ui.column():#.classes('flex-grow min-w-0'):
                                                ui.label("Output:").classes('text-gray-600')
                                                ui.label(f"{self.output}").classes('break-all text-gray-900 font-mono')
                                
                                
                                # BAM File Statistics
                                with ui.column().classes('space-y-3 pb-4 border-b border-gray-200 w-full'):
                                    ui.label("BAM File Summary").classes('font-medium text-gray-900')
                                    with ui.card().classes('w-full p-3 sm:p-4 bg-gray-50 rounded-lg'):
                                        with ui.column().classes('space-y-3'):
                                            # Total Files
                                            with ui.row().classes('items-center justify-between w-full'):
                                                with ui.row().classes('items-center gap-2'):
                                                    ui.icon('description').classes('text-blue-600 shrink-0')
                                                    ui.label("Total Files:").classes('text-gray-600')
                                                total_files = ui.label().classes('font-medium text-gray-900')
                                                def update_total_files():
                                                    total = len(app.storage.general[self.mainuuid]["bam_count"]["file"])
                                                    total_files.text = f"{total:,}" if total > 0 else "--"
                                                app.storage.general[self.mainuuid]["bam_count"].on_change(update_total_files)
                                                update_total_files()
                                            
                                            # Passed Files
                                            with ui.row().classes('items-center justify-between w-full'):
                                                with ui.row().classes('items-center gap-2'):
                                                    ui.icon('check_circle').classes('text-green-600 shrink-0')
                                                    ui.label("Passed:").classes('text-gray-600')
                                                passed = ui.label().classes('font-medium text-gray-900')
                                                def update_passed():
                                                    passed.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_passed']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_passed'] is not None else "--"
                                                app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_passed)
                                                update_passed()
                                            
                                            # Failed Files
                                            with ui.row().classes('items-center justify-between w-full'):
                                                with ui.row().classes('items-center gap-2'):
                                                    ui.icon('error').classes('text-red-600 shrink-0')
                                                    ui.label("Failed:").classes('text-gray-600')
                                                failed = ui.label().classes('font-medium text-gray-900')
                                                def update_failed():
                                                    failed.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_failed']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_failed'] is not None else "--"
                                                app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_failed)
                                                update_failed()

                                # Detailed Statistics Grid
                                with ui.column().classes('space-y-3 w-full'):
                                    ui.label("Sequencing Statistics").classes('font-medium text-gray-900')
                                    with ui.grid().classes('grid-cols-1 sm:grid-cols-3 gap-2 sm:gap-4 w-full'):
                                        # Read Statistics
                                        with ui.card().classes('w-full p-2 sm:p-3 bg-gray-50 rounded-lg'):
                                            ui.label("Read Statistics").classes('text-sm font-medium text-gray-700 mb-2')
                                            with ui.column().classes('space-y-1.5'):
                                                # Total Reads
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Total:").classes('text-gray-600 text-sm')
                                                    total_reads = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_total_reads():
                                                        mapped = app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"]["mapped_count"]
                                                        unmapped = app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"]["unmapped_count"]
                                                        total = mapped + unmapped if mapped is not None and unmapped is not None else 0
                                                        total_reads.text = f"{total:,}" if total > 0 else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_total_reads)
                                                    update_total_reads()

                                                # Mapped Reads
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Mapped:").classes('text-gray-600 text-sm')
                                                    mapped = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_mapped():
                                                        mapped.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['mapped_count']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['mapped_count'] is not None else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_mapped)
                                                    update_mapped()

                                                # Unmapped Reads
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Unmapped:").classes('text-gray-600 text-sm')
                                                    unmapped = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_unmapped():
                                                        unmapped.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['unmapped_count']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['unmapped_count'] is not None else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_unmapped)
                                                    update_unmapped()

                                        # Mapping Quality
                                        with ui.card().classes('p-2 sm:p-3 bg-gray-50 rounded-lg'):
                                            ui.label("Mapping Quality").classes('text-sm font-medium text-gray-700 mb-2')
                                            with ui.column().classes('space-y-1.5'):
                                                # Total Mapped
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Total:").classes('text-gray-600 text-sm')
                                                    total_mapped = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_total_mapped():
                                                        pass_mapped = app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"]["pass_mapped_count"]
                                                        fail_mapped = app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"]["fail_mapped_count"]
                                                        total = pass_mapped + fail_mapped if pass_mapped is not None and fail_mapped is not None else 0
                                                        total_mapped.text = f"{total:,}" if total > 0 else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_total_mapped)
                                                    update_total_mapped()

                                                # Pass Mapped
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Pass:").classes('text-gray-600 text-sm')
                                                    pass_mapped = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_pass_mapped():
                                                        pass_mapped.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['pass_mapped_count']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['pass_mapped_count'] is not None else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_pass_mapped)
                                                    update_pass_mapped()

                                                # Fail Mapped
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Fail:").classes('text-gray-600 text-sm')
                                                    fail_mapped = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_fail_mapped():
                                                        fail_mapped.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['fail_mapped_count']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['fail_mapped_count'] is not None else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_fail_mapped)
                                                    update_fail_mapped()

                                        # Base Statistics
                                        with ui.card().classes('p-2 sm:p-3 bg-gray-50 rounded-lg'):
                                            ui.label("Base Statistics").classes('text-sm font-medium text-gray-700 mb-2')
                                            with ui.column().classes('space-y-1.5'):
                                                # Total Bases
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Total:").classes('text-gray-600 text-sm')
                                                    total_bases = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_total_bases():
                                                        total = app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"]["bases_count"]
                                                        total_bases.text = f"{total:,}" if total > 0 else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_total_bases)
                                                    update_total_bases()

                                                # Pass Bases
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Pass:").classes('text-gray-600 text-sm')
                                                    pass_bases = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_pass_bases():
                                                        pass_bases.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['pass_bases_count']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['pass_bases_count'] is not None else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_pass_bases)
                                                    update_pass_bases()

                                                # Fail Bases
                                                with ui.row().classes('justify-between items-center'):
                                                    ui.label("Fail:").classes('text-gray-600 text-sm')
                                                    fail_bases = ui.label().classes('font-medium text-gray-900 text-sm')
                                                    def update_fail_bases():
                                                        fail_bases.text = f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['fail_bases_count']:,}" if app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['fail_bases_count'] is not None else "--"
                                                    app.storage.general[self.mainuuid]["samples"][self.sampleID]["file_counters"].on_change(update_fail_bases)
                                                    update_fail_bases()
                                    
                # Detailed Analysis Tabs
                if sample_id:
                    selectedtab = None
                    with ui.tabs().classes('w-full') as tabs:
                        if not (set(["sturgeon", "pannanodx", "nanodx", "forest"]).issubset(set(self.exclude))):
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
                
                    with ui.tab_panels(tabs, value=selectedtab).classes('w-full'):
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
                            set(["sturgeon", "nanodx", "forest"]).issubset(set(self.exclude))
                        ):
                            with ui.tab_panel(methylationtab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    ui.label("Methylation Classifications").classes('text-sky-600 dark:text-white').style(
                                        "font-size: 150%; font-weight: 300"
                                    ).tailwind("drop-shadow", "font-bold")
                                    if "sturgeon" not in self.exclude:
                                        self.Sturgeon = Sturgeon_object(
                                            analysis_name="STURGEON",
                                            batch=True,
                                            summary=sturgeonsummary,
                                            **display_args,
                                        )
                                        await self.Sturgeon.render_ui(sample_id=self.sampleID)
                                    if "nanodx" not in self.exclude:
                                        self.NanoDX = NanoDX_object(
                                            analysis_name="NANODX",
                                            batch=True,
                                            summary=nanodxsummary,
                                            **display_args,
                                        )
                                        await self.NanoDX.render_ui(sample_id=self.sampleID)
                                    if "pannanodx" not in self.exclude:
                                        self.PanNanoDX = NanoDX_object(
                                            analysis_name="PANNANODX",
                                            batch=True,
                                            summary=pannanodxsummary,
                                            model="pancan_devel_v5i_NN.pkl",
                                            **display_args,
                                        )
                                        await self.PanNanoDX.render_ui(sample_id=self.sampleID)
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
                                        bed_file = self.bed_file,
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
                                        analysis_name="MGMT", summary=mgmt, **display_args
                                    )
                                    await self.MGMT_panel.render_ui(sample_id=self.sampleID)

                        if "fusion" not in self.exclude:
                            with ui.tab_panel(fusionstab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.Fusion_panel = FusionObject(
                                        analysis_name="FUSION",
                                        summary=fusions,
                                        target_panel=self.target_panel,
                                        **display_args,
                                    )
                                    await self.Fusion_panel.render_ui(sample_id=self.sampleID)

                    async def download_report():
                        '''
                        Generate and download the report.

                        :return: None
                        '''
                        ui.notify("Generating Report")
                        if not self.browse:
                            for item in app.storage.general[self.mainuuid]:
                                if item == "sample_ids":
                                    for sample in app.storage.general[self.mainuuid][item]:
                                        self.sampleID = sample
                        if self.browse:
                            myfile = await run.io_bound(
                                create_pdf,
                                f"{self.sampleID}_run_report.pdf",
                                self.check_and_create_folder(self.output, self.sampleID),
                            )
                        else:
                            myfile = await run.io_bound(
                                create_pdf, f"{self.sampleID}_run_report.pdf", self.output
                            )
                        ui.download(myfile)
                        ui.notify("Report Downloaded")
                        
                

                ui.button("Generate Report", on_click=download_report, icon="download")

        except Exception as e:
            logging.error(f"Error rendering analysis UI: {str(e)}", exc_info=True)
            ui.notify(f"Error rendering analysis UI: {str(e)}", type="error")

    async def background_process_bams(self):
        """
        Background task to process BAM files.

        :return: None
        """
        logging.info("Starting background BAM processing")
        self.process_bams_tracker = ui.timer(10, self.process_bams)

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
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["bam_passed"] += 1
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_mapped_count"] += bamdata["mapped_reads"]
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_unmapped_count"] += bamdata["unmapped_reads"]
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["pass_bases_count"] += bamdata["yield_tracking"]
                else:
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["bam_failed"] += 1
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_mapped_count"] += bamdata["mapped_reads"]
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_unmapped_count"] += bamdata["unmapped_reads"]
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "file_counters"
                    ]["fail_bases_count"] += bamdata["yield_tracking"]
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
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "file_counters"
                ]["mapped_count"] += bamdata["mapped_reads"]
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "file_counters"
                ]["unmapped_count"] += bamdata["unmapped_reads"]
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "file_counters"
                ]["bases_count"] += bamdata["yield_tracking"]

                mydf = pd.DataFrame.from_dict(app.storage.general)
                if not self.force_sampleid:
                    sample_id = baminfo["sample_id"]
                else:
                    sample_id = self.force_sampleid
                mydf.to_csv(
                    os.path.join(
                        self.check_and_create_folder(self.output, sample_id),
                        "master.csv"
                    )
                )

                counter += 1
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
            if "file" not in app.storage.general[self.mainuuid]["bam_count"]:
                app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
            app.storage.general[self.mainuuid]["bam_count"]["file"][
                f
            ] = timestamp.timestamp()  # time.time()
            if self.simtime:
                await asyncio.sleep(1)
            else:
                await asyncio.sleep(0)
