"""
Main Module for Initializing and Configuring the Brain Tumour Classification Application

This module sets up and configures the Brain Tumour Classification application using the NiceGUI framework (https://nicegui.io/). The application ensures that all long-running processes operate in the background, preventing the main thread from being blocked.

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

- `nicegui` (ui, app)

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
from robin.reporting.generate_report import create_pdf

from watchdog.observers import Observer
from pathlib import Path
from queue import Queue
import pandas as pd
import asyncio
import pysam
import logging
from collections import Counter



import time
from datetime import datetime
import pytz
import os

# Configure logging
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
#logger = logging.getLogger(__name__)


def check_bam(bamfile):
    """
    This function checks a bam file and returns a dictionary of its attributes.
    :param bamfile:
    :return:
    """
    pysam.index(bamfile)
    bam = ReadBam(bamfile)
    baminfo = bam.process_reads()
    bamdata = bam.summary()
    return baminfo, bamdata


class BrainMeth:
    def __init__(
        self,
        threads=4,
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
        mainuuid=None,
    ):
        self.mainuuid = mainuuid
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
        self.observer = None
        if self.browse:
            self.runsfolder = self.output
        self.sampleID=None

    async def start_background(self):
        #logger.debug("Init Brain Class")
        app.storage.general[self.mainuuid]["bam_count"] = {}
        app.storage.general[self.mainuuid]["bam_count"] = Counter(counter=0)
        app.storage.general[self.mainuuid]["bam_count"]["files"] = {}
        app.storage.general[self.mainuuid]["file_counters"] = Counter(
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
        app.storage.general[self.mainuuid]["devices"] = []  # set
        app.storage.general[self.mainuuid]["basecall_models"] = []  # set()
        app.storage.general[self.mainuuid]["run_time"] = []  # set()
        app.storage.general[self.mainuuid]["flowcell_ids"] = []  # set()
        app.storage.general[self.mainuuid]["sample_ids"] = []  # set()

        self.bam_tracking = Queue()
        self.bamforcns = Queue()
        self.bamforsturgeon = Queue()
        self.bamfornanodx = Queue()
        self.bamforcnv = Queue()
        self.bamfortargetcoverage = Queue()
        self.bamformgmt = Queue()
        self.bamforfusions = Queue()

        if self.watchfolder:
            print(f"Adding a watchfolder {self.watchfolder}")
            await self.add_watchfolder(self.watchfolder)

        if "sturgeon" not in self.exclude:
            self.Sturgeon = Sturgeon_object(
                self.threads,
                self.output,
                "STURGEON",
                progress=True,
                batch=True,
                bamqueue=self.bamforsturgeon,
                browse=self.browse,
                uuid=self.mainuuid,
            )
            self.Sturgeon.process_data()

        if "nanodx" not in self.exclude:
            self.NanoDX = NanoDX_object(
                self.threads,
                self.output,
                "NANODX",
                progress=True,
                batch=True,
                bamqueue=self.bamfornanodx,
                browse=self.browse,
                uuid=self.mainuuid,
            )
            self.NanoDX.process_data()

        if "forest" not in self.exclude:
            self.RandomForest = RandomForest_object(
                self.threads,
                self.output,
                "FOREST",
                progress=True,
                batch=True,
                showerrors=self.showerrors,
                bamqueue=self.bamforcns,
                browse=self.browse,
                uuid=self.mainuuid,
            )
            self.RandomForest.process_data()

        if "cnv" not in self.exclude:
            self.CNV = CNVAnalysis(
                self.threads,
                self.output,
                "CNV",
                progress=True,
                bamqueue=self.bamforcnv,
                # summary=cnvsummary,
                target_panel=self.target_panel,
                browse=self.browse,
                uuid=self.mainuuid,
            )
            self.CNV.process_data()

        if "coverage" not in self.exclude:
            self.Target_Coverage = TargetCoverage(
                self.threads,
                self.output,
                "COVERAGE",
                progress=True,
                bamqueue=self.bamfortargetcoverage,
                #    summary=coverage,
                target_panel=self.target_panel,
                reference=self.reference,
                browse=self.browse,
                uuid=self.mainuuid,
            )
            self.Target_Coverage.process_data()

        if "mgmt" not in self.exclude:
            self.MGMT_panel = MGMT_Object(
                self.threads,
                self.output,
                "MGMT",
                progress=True,
                bamqueue=self.bamformgmt,
                browse=self.browse,
                uuid=self.mainuuid,
            )
            self.MGMT_panel.process_data()

        if "fusion" not in self.exclude:
            self.Fusion_panel = FusionObject(
                self.threads,
                self.output,
                "FUSION",
                progress=True,
                bamqueue=self.bamforfusions,
                browse=self.browse,
                target_panel=self.target_panel,
                uuid=self.mainuuid,
            )
            self.Fusion_panel.process_data()

        # if self.minknow_connection:
        #    while not self.minknow_connection.connected:
        #        await asyncio.sleep(1)
        #    while not self.minknow_connection.selected_position:
        #        await asyncio.sleep(1)
        #    await self.minknow_connection.access_device.clicked()

    async def render_ui(self):
        if not self.browse:
            await self.information_panel()

            # if self.watchfolder:
            #    self.add_watchfolder(self.watchfolder)

        else:
            ui.label("Browse mode enabled. Please choose a folder to see data from.")
            ui.button("Choose file", on_click=self.pick_file, icon="folder")

            self.content = ui.column().classes("w-full")

    async def add_watchfolder(self, watchfolder):
        if not self.observer:
            self.watchfolder = watchfolder
            if "file" not in app.storage.general[self.mainuuid]["bam_count"].keys():
                app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
            self.event_handler = BamEventHandler(
                app.storage.general[self.mainuuid]["bam_count"]
            )
            self.observer = Observer()
            self.observer.schedule(self.event_handler, self.watchfolder, recursive=True)
            self.observer.start()

            ui.timer(1, callback=self.background_process_bams, once=True)

            ui.timer(1, callback=self.check_existing_bams, once=True)

            print("watchfolder setup and added")

    async def waitforclick(self):
        self.minknow_connection.check_connection()
        await self.minknow_connection.access_device.clicked()

    async def pick_file(self) -> None:
        result = await LocalFilePicker(f"{self.runsfolder}", multiple=True)
        if result:
            ui.notify(f"You selected {result}")
            self.content.clear()
            with self.content:
                ui.label(f"Monitoring the path:{result}").tailwind(
                    "drop-shadow", "font-bold"
                )
                self.output = result[0]
                await self.information_panel()

    def replay(self):
        ui.notify("Replaying data")
        self.replaycontrol.visible = False
        self.rcns2_worker.replay_prior_data()
        self.sturgeon_worker.replay_prior_data()

    def replay_cnv(self):
        ui.notify("Replaying CNV data")

    @property
    def min_start_time(self):
        if len(self.run_time) > 0:
            dt = datetime.fromisoformat(min(self.run_time))
            formatted_string = dt.strftime("%Y-%m-%d %H:%M")
            return formatted_string
        else:
            return 0

    async def information_panel(self):
        self.frontpage = ui.card().classes("w-screen")
        with self.frontpage:
            ui.label("CNS Tumour Methylation Classification").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            ui.label(
                "This tool enables classification of brain tumours in real time from Oxford Nanopore Data."
            ).style("color: #000000; font-size: 100%; font-weight: 300").tailwind(
                "drop-shadow", "font-bold"
            )
            with ui.row():
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid],
                    "devices",
                    backward=lambda n: [f"Devices: {str(item)}" for item in n],
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid],
                    "basecall_models",
                    backward=lambda n: [
                        f"Basecall Models: {str(item)}" for item in n
                    ],
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid],
                    "flowcell_ids",
                    backward=lambda n: [f"Flowcell IDs: {str(item)}" for item in n],
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid],
                    "run_time",
                    backward=lambda n: [
                        f"Run Start Time: {datetime.strptime(date.replace('Z', ''), '%Y-%m-%dT%H:%M:%S').replace(tzinfo=pytz.UTC).strftime('%Y-%m-%d %H:%M:%S %Z')}"
                        for date in n
                    ],
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid],
                    "sample_ids",
                    backward=lambda n: [f"Sample ID: {str(item)}" for item in n],
                ).style("color: #000000; font-size: 100%; font-weight: 300")

            with ui.row():
                ui.label("Results Summary").style(
                    "color: #6E93D6; font-size: 125%; font-weight: 300"
                ).tailwind("drop-shadow", "font-bold")

            with ui.row().style("color: #000000; font-size: 100%; font-weight: 300"):
                if "sturgeon" not in self.exclude:
                    sturgeonsummary = ui.column()
                if "nanodx" not in self.exclude:
                    nanodxsummary = ui.column()
                if "forest" not in self.exclude:
                    forestsummary = ui.column()
            if "mgmt" not in self.exclude:
                with ui.row().style(
                    "color: #000000; font-size: 100%; font-weight: 300"
                ):
                    mgmt = ui.column()
            if "cnv" not in self.exclude:
                with ui.row().style(
                    "color: #000000; font-size: 100%; font-weight: 300"
                ):
                    if "cnv" not in self.exclude:
                        cnvsummary = ui.column()
            if "fusion" not in self.exclude:
                with ui.row().style(
                    "color: #000000; font-size: 100%; font-weight: 300"
                ):
                    if "fusion" not in self.exclude:
                        fusions = ui.column()
            if "coverage" not in self.exclude:
                with ui.row().style(
                    "color: #000000; font-size: 100%; font-weight: 300"
                ):
                    if "coverage" not in self.exclude:
                        coverage = ui.column()

            with ui.expansion(icon="work").bind_text_from(
                self, "watchfolder", backward=lambda n: f"Monitoring the path: {n}"
            ).classes(""):
                with ui.row():
                    ui.label().bind_text_from(
                        self,
                        "watchfolder",
                        backward=lambda n: f"Monitoring the path: {n}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label(f"Outputting to:{self.output}").style(
                        "color: #000000; font-size: 100%; font-weight: 300"
                    )
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["bam_count"],
                        "counter",
                        backward=lambda n: f"BAM files seen: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "bam_passed",
                        backward=lambda n: f"BAM pass: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "bam_failed",
                        backward=lambda n: f"BAM fail: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                with ui.row():
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "mapped_count",
                        backward=lambda n: f"Mapped Read Count: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "unmapped_count",
                        backward=lambda n: f"Unmapped Read Count: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "pass_mapped_count",
                        backward=lambda n: f"Pass Mapped Read Count: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "fail_mapped_count",
                        backward=lambda n: f"Fail Mapped Read Count: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                with ui.row():
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "bases_count",
                        backward=lambda n: f"Total Mapped Bases: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "pass_bases_count",
                        backward=lambda n: f"Total Mapped Pass Bases: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid]["file_counters"],
                        "fail_bases_count",
                        backward=lambda n: f"Total Mapped Fail Bases: {n:,}",
                    ).style("color: #000000; font-size: 100%; font-weight: 300")
                with ui.row():
                    if "forest" not in self.exclude:
                        ui.label().bind_text_from(
                            self,
                            "bamforcns",
                            backward=lambda n: f"BAM files for CNS: {n.qsize()}",
                        ).style("color: #000000; font-size: 100%; font-weight: 300")
                    if "sturgeon" not in self.exclude:
                        ui.label().bind_text_from(
                            self,
                            "bamforsturgeon",
                            backward=lambda n: f"BAM files for Sturgeon: {n.qsize()}",
                        ).style("color: #000000; font-size: 100%; font-weight: 300")
                    if "nanodx" not in self.exclude:
                        ui.label().bind_text_from(
                            self,
                            "bamfornanodx",
                            backward=lambda n: f"BAM files for NanoDX: {n.qsize()}",
                        ).style("color: #000000; font-size: 100%; font-weight: 300")

        selectedtab = None
        await ui.context.client.connected()
        with ui.tabs().classes("w-full") as tabs:
            if not (set(["sturgeon", "nanodx", "forest"]).issubset(set(self.exclude))):
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
        with ui.tab_panels(tabs, value=selectedtab).classes("w-screen"):
            if not (set(["sturgeon", "nanodx", "forest"]).issubset(set(self.exclude))):
                with ui.tab_panel(methylationtab).classes("w-full"):
                    with ui.card().style("width: 100%"):
                        ui.label("Methylation Classifications").style(
                            "color: #6E93D6; font-size: 150%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        if "sturgeon" not in self.exclude:
                            self.Sturgeon = Sturgeon_object(
                                self.threads,
                                self.output,
                                "STURGEON",
                                progress=True,
                                batch=True,
                                bamqueue=None,
                                summary=sturgeonsummary,
                                browse=self.browse,
                                uuid=self.mainuuid,
                            )
                            await self.Sturgeon.render_ui()
                        if "nanodx" not in self.exclude:
                            self.NanoDX = NanoDX_object(
                                self.threads,
                                self.output,
                                "NANODX",
                                progress=True,
                                batch=True,
                                bamqueue=None,
                                summary=nanodxsummary,
                                browse=self.browse,
                                uuid=self.mainuuid,
                            )
                            await self.NanoDX.render_ui()
                        if "forest" not in self.exclude:
                            self.RandomForest = RandomForest_object(
                                self.threads,
                                self.output,
                                "FOREST",
                                progress=True,
                                batch=True,
                                bamqueue=None,
                                summary=forestsummary,
                                showerrors=self.showerrors,
                                browse=self.browse,
                                uuid=self.mainuuid,
                            )
                            await self.RandomForest.render_ui()
            if "cnv" not in self.exclude:
                with ui.tab_panel(copy_numbertab).classes("w-full"):
                    with ui.card().classes("w-full"):
                        self.CNV = CNVAnalysis(
                            self.threads,
                            self.output,
                            "CNV",
                            progress=True,
                            bamqueue=None,  # self.bamforcnv,
                            summary=cnvsummary,
                            target_panel=self.target_panel,
                            browse=self.browse,
                            uuid=self.mainuuid,
                        )
                        await self.CNV.render_ui()
            if "coverage" not in self.exclude:
                with ui.tab_panel(coveragetab).classes("w-full"):
                    with ui.card().classes("w-full"):
                        self.Target_Coverage = TargetCoverage(
                            self.threads,
                            self.output,
                            "COVERAGE",
                            progress=True,
                            bamqueue=None,  # self.bamfortargetcoverage,
                            summary=coverage,
                            target_panel=self.target_panel,
                            reference=self.reference,
                            browse=self.browse,
                            uuid=self.mainuuid,
                        )
                        await self.Target_Coverage.render_ui()
            if "mgmt" not in self.exclude:
                with ui.tab_panel(mgmttab).classes("w-full"):
                    with ui.card().classes("w-full"):
                        self.MGMT_panel = MGMT_Object(
                            self.threads,
                            self.output,
                            "MGMT",
                            progress=True,
                            bamqueue=None,
                            summary=mgmt,
                            browse=self.browse,
                            uuid=self.mainuuid,
                        )
                        await self.MGMT_panel.render_ui()
            if "fusion" not in self.exclude:
                with ui.tab_panel(fusionstab).classes("w-full"):
                    with ui.card().classes("w-full"):
                        self.Fusion_panel = FusionObject(
                            self.threads,
                            self.output,
                            "FUSION",
                            progress=True,
                            bamqueue=None,
                            summary=fusions,
                            target_panel=self.target_panel,
                            browse=self.browse,
                            uuid=self.mainuuid,
                        )
                        await self.Fusion_panel.render_ui()

        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == 'sample_ids':
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
        ui.button("Generate Report", on_click=lambda: create_pdf(f"{self.sampleID}_run_report.pdf", self.check_and_create_folder(self.output,self.sampleID)))

    async def background_process_bams(self):
        await asyncio.sleep(5)
        self.process_bams_tracker = ui.timer(10, self.process_bams)

    def check_and_create_folder(self, path, folder_name=None):
        # Check if the path exists
        if not os.path.exists(path):
            raise FileNotFoundError(f"The specified path does not exist: {path}")

        # If folder_name is provided
        if folder_name:
            full_path = os.path.join(path, folder_name)
            # Create the folder if it doesn't exist
            if not os.path.exists(full_path):
                os.makedirs(full_path)
            return full_path
        else:
            return path

    async def process_bams(self) -> None:
        """
        This function processes the bam files and adds them to the bamforcns and bamforsturgeon lists.
        These lists are then processed by the rapid_cns2 and sturgeon functions.
        #ToDo: switch to using threadsafe queues.
        :param self:
        :return:
        """
        self.process_bams_tracker.active = False
        counter = 0
        # while True:
        # print(f"We're processing a bam man - {self.mainuuid}")
        if "file" in app.storage.general[self.mainuuid]["bam_count"]:
            # print (app.storage.general[self.mainuuid]['bam_count']["file"])
            while len(app.storage.general[self.mainuuid]["bam_count"]["file"]) > 0:
                self.nofiles = False
                file = app.storage.general[self.mainuuid]["bam_count"]["file"].popitem()
                # ToDo: Check if the file is still being written to.
                #if file[1] > time.time() - 5:
                #    break
                # ToDo: This function should be moved to a background task.
                # self.check_bam(file[0])
                #self.bam_tracking.put(file[0])

                baminfo, bamdata = await run.cpu_bound(check_bam, file[0])
                if baminfo["state"] == "pass":
                    app.storage.general[self.mainuuid]["file_counters"][
                        "bam_passed"
                    ] += 1
                    app.storage.general[self.mainuuid]["file_counters"][
                        "pass_mapped_count"
                    ] += bamdata['mapped_reads']
                    app.storage.general[self.mainuuid]["file_counters"][
                        "pass_unmapped_count"
                    ] += bamdata['unmapped_reads']
                    app.storage.general[self.mainuuid]["file_counters"][
                        "pass_bases_count"
                    ] += bamdata['yield_tracking']
                else:
                    app.storage.general[self.mainuuid]["file_counters"][
                        "bam_failed"
                    ] += 1
                    app.storage.general[self.mainuuid]["file_counters"][
                        "fail_mapped_count"
                    ] += bamdata['mapped_reads']
                    app.storage.general[self.mainuuid]["file_counters"][
                        "fail_unmapped_count"
                    ] += bamdata['unmapped_reads']
                    app.storage.general[self.mainuuid]["file_counters"][
                        "fail_bases_count"
                    ] += bamdata['yield_tracking']
                    # self.basecall_models.add(baminfo["basecall_model"])
                if (
                        baminfo["device_position"]
                        not in app.storage.general[self.mainuuid]["devices"]
                ):
                    app.storage.general[self.mainuuid]["devices"].append(
                        baminfo["device_position"]
                    )
                if (
                        baminfo["basecall_model"]
                        not in app.storage.general[self.mainuuid]["basecall_models"]
                ):
                    app.storage.general[self.mainuuid]["basecall_models"].append(
                        baminfo["basecall_model"]
                    )
                if (
                        baminfo["sample_id"]
                        not in app.storage.general[self.mainuuid]["sample_ids"]
                ):
                    app.storage.general[self.mainuuid]["sample_ids"].append(
                        baminfo["sample_id"]
                    )
                if (
                        baminfo["flow_cell_id"]
                        not in app.storage.general[self.mainuuid]["flowcell_ids"]
                ):
                    app.storage.general[self.mainuuid]["flowcell_ids"].append(
                        baminfo["flow_cell_id"]
                    )
                if (
                        baminfo["time_of_run"]
                        not in app.storage.general[self.mainuuid]["run_time"]
                ):
                    app.storage.general[self.mainuuid]["run_time"].append(
                        baminfo["time_of_run"]
                    )
                # self.sample_ids.add(baminfo["sample_id"])
                # self.flowcell_ids.add(baminfo["flow_cell_id"])
                # self.run_time.add(baminfo["time_of_run"])
                app.storage.general[self.mainuuid]["file_counters"][
                    "mapped_count"
                ] += bamdata['mapped_reads']
                app.storage.general[self.mainuuid]["file_counters"][
                    "unmapped_count"
                ] += bamdata['unmapped_reads']
                app.storage.general[self.mainuuid]["file_counters"][
                    "bases_count"
                ] += bamdata['yield_tracking']

                mydf = pd.DataFrame.from_dict(app.storage.general)

                mydf.to_csv(os.path.join(self.check_and_create_folder(self.output,baminfo["sample_id"]), "master.csv"))



                counter += 1
                if "forest" not in self.exclude:
                    self.bamforcns.put([file[0], file[1],baminfo["sample_id"]])
                if "sturgeon" not in self.exclude:
                    self.bamforsturgeon.put([file[0], file[1],baminfo["sample_id"]])
                if "nanodx" not in self.exclude:
                    self.bamfornanodx.put([file[0], file[1],baminfo["sample_id"]])
                if "cnv" not in self.exclude:
                    self.bamforcnv.put([file[0], file[1],baminfo["sample_id"]])
                if "coverage" not in self.exclude:
                    self.bamfortargetcoverage.put([file[0], file[1],baminfo["sample_id"]])
                if "mgmt" not in self.exclude:
                    self.bamformgmt.put([file[0], file[1],baminfo["sample_id"]])
                if "fusion" not in self.exclude:
                    self.bamforfusions.put([file[0], file[1],baminfo["sample_id"]])

            self.nofiles = True
        self.process_bams_tracker.active = True

    async def check_existing_bams(self, sequencing_summary=None):
        file_endings = {".bam"}
        if sequencing_summary:
            print("Using sequencing summary")
            # with self.frontpage:
            #    ui.notify(
            #        "Checking for existing bams against sequencing summary",
            #        type="info",
            #        position="top-right",
            #    )
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
            latest_timestamps["file_produced"] = latest_timestamps[
                "template_end"
            ]  # + now
            for path, dirs, files in os.walk(self.watchfolder):
                for f in files:
                    if "".join(Path(f).suffixes) in file_endings:
                        latest_timestamps.loc[
                            latest_timestamps["filename_bam"] == f, "full_path"
                        ] = os.path.join(path, f)

            step_size = 20

            # with self.frontpage:
            #    ui.notify("Target Playback Started", type="info", position="top-right")

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
                    # ToDo: This needs to be called in the background.
                    self.check_bam(row["full_path"])

            self.runfinished = True
        else:
            for path, dirs, files in os.walk(self.watchfolder):
                for f in files:
                    if "".join(Path(f).suffixes) in file_endings:
                        app.storage.general[self.mainuuid]["bam_count"]["counter"] += 1
                        if (
                            "file"
                            not in app.storage.general[self.mainuuid]["bam_count"]
                        ):
                            app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
                        app.storage.general[self.mainuuid]["bam_count"]["file"][
                            os.path.join(path, f)
                        ] = time.time()
                        if self.simtime:
                            await asyncio.sleep(1)
                        else:
                            await asyncio.sleep(0)
