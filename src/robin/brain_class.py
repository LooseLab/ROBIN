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
        ui.notify("Select a folder to monitor")
        result = await LocalFilePicker(f"{self.runsfolder}", multiple=True)
        if result:
            ui.notify(f"You selected {result}")
            self.content.clear()
            with self.content:
                ui.label(f"Monitoring the path: {result}").tailwind(
                    "drop-shadow", "font-bold"
                )
                self.output, self.sampleID = os.path.split(result[0])
                await self.information_panel(sample_id=self.sampleID)

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
        Render the information panel for a sample.

        :param sample_id: Sample ID to render the panel for.
        :return: None
        """
        if sample_id:
            self.sampleID = sample_id
        if not self.browse:
            self.show_list()
            app.storage.general[self.mainuuid]["sample_list"].on_change(
                self.show_list.refresh
            )

        if sample_id:
            self.frontpage = ui.card().classes("rounded w-full")
            with (self.frontpage):
                if not self.browse:
                    if (
                        self.sampleID
                        not in app.storage.general[self.mainuuid]["samples"]
                    ):
                        ui.notify(f"Sample {self.sampleID} not found")
                        ui.navigate.to("/live")
                        return

                ui.label(
                    f"CNS Tumor Methylation Classification - {self.sampleID}"
                ).classes('text-sky-600 dark:text-white').style(
                    "font-size: 150%; font-weight: 300"
                ).tailwind("drop-shadow", "font-bold")

                ui.label(
                    "This tool enables classification of brain tumors in real time from Oxford Nanopore Data."
                ).style("font-size: 100%; font-weight: 300").classes(
                    "drop-shadow font-bold"
                )

                with ui.row():
                    ui.label("Results Summary").classes('text-sky-600 dark:text-white').style(
                    "font-size: 150%; font-weight: 300"
                ).tailwind("drop-shadow", "font-bold")

                with ui.row().style(
                    "font-size: 120%; font-weight: 300"
                ):
                    if "sturgeon" not in self.exclude:
                        sturgeonsummary = ui.column()
                    if "nanodx" not in self.exclude:
                        nanodxsummary = ui.column()
                    if "pannanodx" not in self.exclude:
                        pannanodxsummary = ui.column()
                    if "forest" not in self.exclude:
                        forestsummary = ui.column()

                with ui.row():
                    if not self.browse:
                        ui.label().bind_text_from(
                            app.storage.general[self.mainuuid]["samples"][
                                self.sampleID
                            ],
                            "devices",
                            backward=lambda n: [f"Devices: {str(item)}" for item in n],
                        ).style("font-size: 100%; font-weight: 300")
                        ui.label().bind_text_from(
                            app.storage.general[self.mainuuid]["samples"][
                                self.sampleID
                            ],
                            "basecall_models",
                            backward=lambda n: [
                                f"Basecall Models: {str(item)}" for item in n
                            ],
                        ).style("font-size: 100%; font-weight: 300")
                        ui.label().bind_text_from(
                            app.storage.general[self.mainuuid]["samples"][
                                self.sampleID
                            ],
                            "flowcell_ids",
                            backward=lambda n: [
                                f"Flowcell IDs: {str(item)}" for item in n
                            ],
                        ).style("font-size: 100%; font-weight: 300")
                        ui.label().bind_text_from(
                            app.storage.general[self.mainuuid]["samples"][
                                self.sampleID
                            ],
                            "run_time",
                            backward=lambda n: [
                                f"Run Start Time: {parser.parse(date).astimezone(pytz.UTC).strftime('%Y-%m-%d %H:%M:%S %Z')}"
                                for date in n
                            ],
                        ).style("font-size: 100%; font-weight: 300")
                        ui.label().bind_text_from(
                            app.storage.general[self.mainuuid]["samples"][
                                self.sampleID
                            ],
                            "sample_ids",
                            backward=lambda n: [
                                f"Sample ID: {str(item)}" for item in n
                            ],
                        ).style("font-size: 100%; font-weight: 300")

                if "mgmt" not in self.exclude:
                    with ui.row().style(
                        "font-size: 120%; font-weight: 300"
                    ):
                        mgmt = ui.column()
                if "cnv" not in self.exclude:
                    with ui.row().style(
                        "font-size: 120%; font-weight: 300"
                    ):
                        cnvsummary = ui.column()
                if "fusion" not in self.exclude:
                    with ui.row().style(
                        "font-size: 120%; font-weight: 300"
                    ):
                        fusions = ui.column()
                if "coverage" not in self.exclude:
                    with ui.row().style(
                        "font-size: 120%; font-weight: 300"
                    ):
                        coverage = ui.column()

                if not self.browse:
                    with ui.expansion(icon="work").bind_text_from(
                        self,
                        "watchfolder",
                        backward=lambda n: f"Monitoring the path: {n}",
                    ):
                        with ui.row():
                            ui.label(f"Outputting to: {self.output}").style(
                                "font-size: 100%; font-weight: 300"
                            )
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["bam_count"],
                                "counter",
                                backward=lambda n: f"BAM files seen: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "bam_passed",
                                backward=lambda n: f"BAM pass: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "bam_failed",
                                backward=lambda n: f"BAM fail: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")

                        with ui.row():
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "mapped_count",
                                backward=lambda n: f"Mapped Read Count: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "unmapped_count",
                                backward=lambda n: f"Unmapped Read Count: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "pass_mapped_count",
                                backward=lambda n: f"Pass Mapped Read Count: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "fail_mapped_count",
                                backward=lambda n: f"Fail Mapped Read Count: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")

                        with ui.row():
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "bases_count",
                                backward=lambda n: f"Total Mapped Bases: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "pass_bases_count",
                                backward=lambda n: f"Total Mapped Pass Bases: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]["file_counters"],
                                "fail_bases_count",
                                backward=lambda n: f"Total Mapped Fail Bases: {n:,}",
                            ).style("font-size: 100%; font-weight: 300")

                        with ui.row():
                            if "forest" not in self.exclude:
                                ui.label().bind_text_from(
                                    self,
                                    "bamforcns",
                                    backward=lambda n: f"BAM files for CNS: {n.qsize()}",
                                ).style(
                                    "font-size: 100%; font-weight: 300"
                                )
                            if "sturgeon" not in self.exclude:
                                ui.label().bind_text_from(
                                    self,
                                    "bamforsturgeon",
                                    backward=lambda n: f"BAM files for Sturgeon: {n.qsize()}",
                                ).style(
                                    "font-size: 100%; font-weight: 300"
                                )
                            if "nanodx" not in self.exclude:
                                ui.label().bind_text_from(
                                    self,
                                    "bamfornanodx",
                                    backward=lambda n: f"BAM files for NanoDX: {n.qsize()}",
                                ).style(
                                    "font-size: 100%; font-weight: 300"
                                )
                            if "pannanodx" not in self.exclude:
                                ui.label().bind_text_from(
                                    self,
                                    "bamforpannanodx",
                                    backward=lambda n: f"BAM files for Pan NanoDX: {n.qsize()}",
                                ).style(
                                    "font-size: 100%; font-weight: 300"
                                )

        if sample_id:
            selectedtab = None
            with ui.tabs().classes("w-full") as tabs:
                if not (
                    set(["sturgeon", "nanodx", "forest"]).issubset(set(self.exclude))
                ):
                    methylationtab = ui.tab("Methylation Classification").classes('text-sky-600 dark:text-white').style(
                        "font-size: 150%; font-weight: 600"
                    )
                    if not selectedtab:
                        selectedtab = methylationtab
                if "cnv" not in self.exclude:
                    copy_numbertab = ui.tab("Copy Number Variation").classes('text-sky-600 dark:text-white').style(
                        "font-size: 150%; font-weight: 600"
                    )
                    if not selectedtab:
                        selectedtab = copy_numbertab
                if "coverage" not in self.exclude:
                    coveragetab = ui.tab("Target Coverage").classes('text-sky-600 dark:text-white').style(
                        "font-size: 150%; font-weight: 600"
                    )
                    if not selectedtab:
                        selectedtab = coveragetab
                if "mgmt" not in self.exclude:
                    mgmttab = ui.tab("MGMT").classes('text-sky-600 dark:text-white').style(
                        "font-size: 150%; font-weight: 600"
                    )
                    if not selectedtab:
                        selectedtab = mgmttab
                if "fusion" not in self.exclude:
                    fusionstab = ui.tab("Fusions").classes('text-sky-600 dark:text-white').style(
                        "font-size: 150%; font-weight: 600"
                    )
                    if not selectedtab:
                        selectedtab = fusionstab
            with ui.tab_panels(tabs, value=selectedtab).classes("w-screen"):
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
                """
                Generate and download the report.

                :return: None
                """
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
                        "master.csv",
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
