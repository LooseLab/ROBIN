from nicegui import ui

from cnsmeth.utilities.bam_handler import BamEventHandler

# from cnsmeth.Sturgeon_worker import Sturgeon_worker
# from cnsmeth.RCNS2_worker import RCNS2_worker

# from cnsmeth.CrossNN_worker import CrossNN_worker
from cnsmeth.subpages.MGMT_object import MGMT_Object
from cnsmeth.subpages.Sturgeon_object import Sturgeon_object
from cnsmeth.subpages.NanoDX_object import NanoDX_object
from cnsmeth.subpages.RandomForest_object import RandomForest_object

# from cnsmeth.subpages.fusion_panel import Fusion_Panel
from cnsmeth.subpages.CNV_object import CNVAnalysis
from cnsmeth.subpages.TargetCoverage_object import TargetCoverage
from cnsmeth.subpages.Fusion_object import Fusion_object
from cnsmeth.utilities.local_file_picker import local_file_picker
from cnsmeth.utilities.ReadBam import ReadBam

from watchdog.observers import Observer
from pathlib import Path
from queue import Queue
import pandas as pd

import threading
import time
from datetime import datetime
import os


class BrainMeth:
    def __init__(
        self,
        threads=4,
        simtime=False,
        watchfolder=None,
        output=None,
        sequencing_summary=None,
        showerrors=False,
        browse=False,
        exclude=[],
    ):
        self.threads = threads
        self.simtime = simtime
        self.watchfolder = watchfolder
        self.output = output
        self.sequencing_summary = sequencing_summary
        self.showerrors = showerrors
        self.browse = browse
        self.exclude = exclude

        self.bam_count = {"counter": 0}
        self.bam_passed = {"counter": 0}
        self.bam_failed = {"counter": 0}
        self.devices = set()
        self.basecall_models = set()
        self.run_time = set()
        self.flowcell_ids = set()
        self.sample_ids = set()
        self.bamforcns = Queue()
        self.bamforsturgeon = Queue()
        self.bamfornanodx = Queue()
        self.bamforcnv = Queue()
        self.bamfortargetcoverage = Queue()
        self.bamformgmt = Queue()
        self.bamforfusions = Queue()

        if not self.browse:
            self.information_panel()
            self.event_handler = BamEventHandler(self.bam_count)
            self.observer = Observer()
            self.observer.schedule(self.event_handler, self.watchfolder, recursive=True)
            self.observer.start()

            self.bam_processing = threading.Thread(target=self.process_bams, args=())
            self.bam_processing.daemon = True
            self.bam_processing.start()

            self.check_for_existing_bams = threading.Thread(
                target=self.check_existing_bams,
                args=(),
                kwargs={"sequencing_summary": self.sequencing_summary},
            )
            self.check_for_existing_bams.daemon = True
            self.check_for_existing_bams.start()

        else:
            ui.label("Browse mode enabled. Please choose a folder to see data from.")
            ui.button("Choose file", on_click=self.pick_file, icon="folder")

            self.content = ui.column().classes("w-full")

    async def pick_file(self) -> None:
        result = await local_file_picker("/", multiple=True)
        print(result)
        if result:
            ui.notify(f"You selected {result}")
            self.content.clear()
            with self.content:
                ui.label(f"Monitoring the path:{result}").tailwind(
                    "drop-shadow", "font-bold"
                )
                ui.label("Not yet implemented.")
                """
                with ui.tabs().classes("w-full") as tabs:
                    methylation = ui.tab("Methylation Classification")
                    copy_numer = ui.tab("Copy Number Variation")
                    coverage = ui.tab("Target Coverage")
                    mgmt = ui.tab("MGMT")
                    fusions = ui.tab("Fusions")
                with ui.tab_panels(tabs, value=methylation).classes("w-full"):
                    with ui.tab_panel(methylation).classes("w-full"):
                        with ui.card().classes("w-full"):
                            with ui.column().classes("w-full"):
                                self.rcns2_worker = RCNS2_worker(
                                    self.bamforcns,
                                    self.cnv,
                                    self.target_coverage,
                                    self.mgmt_panel,
                                    self.fusion_panel,
                                    threads=self.threads,
                                    output_folder=result[0],
                                    showerrors=self.showerrors,
                                    browse=self.browse,
                                )
                                self.sturgeon_worker = Sturgeon_worker(
                                    self.bamforsturgeon,
                                    threads=self.threads,
                                    output_folder=result[0],
                                    showerrors=self.showerrors,
                                    browse=self.browse,
                                )

                                with ui.card().classes("w-full h-auto"):
                                    with ui.grid(columns=3).classes("w-full h-auto"):
                                        with ui.column():
                                            with ui.card().style("width: 100%"):
                                                self.rcns2_worker.status_rcns2()
                                                self.rcns2_worker.create_rcns2_chart(
                                                    "RapidCNS2 (Random Forest)"
                                                )

                                        with ui.column():
                                            with ui.card().style("width: 100%"):
                                                self.sturgeon_worker.status_sturgeon()
                                                self.sturgeon_worker.create_sturgeon_chart(
                                                    "Sturgeon"
                                                )
                                        with ui.column():
                                            with ui.card().style("width: 100%"):
                                                # self.sturgeon_worker.status_sturgeon()
                                                self.sturgeon_worker.create_nanodx_chart(
                                                    "NanoDX"
                                                )

                                with ui.card().style("width: 100%"):
                                    self.rcns2_worker.create_rcns2_time_chart()

                                with ui.card().style("width: 100%"):
                                    self.sturgeon_worker.create_sturgeon_time_chart()

                    with ui.tab_panel(copy_numer).classes("w-full"):
                        with ui.card().style("width: 100%"):
                            self.cnv.create_cnv_scatter("CNV Scatter")

                    with ui.tab_panel(coverage).classes("w-full"):
                        pass

                    with ui.tab_panel(mgmt).classes("w-full"):
                        self.mgmt_panel.setup_ui(mgmt)

                    with ui.tab_panel(fusions).classes("w-full"):
                        self.fusion_panel.setup_ui()

                    self.rcns2_worker.load_prior_data()
                    self.sturgeon_worker.load_prior_data()

                self.replaycontrol = ui.button(
                    "replay data", on_click=self.replay, icon="replay"
                )
                """
            # self.rcns2_worker.replay_prior_data()

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

    def information_panel(self):
        self.frontpage = ui.card().style("width: 100%")
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
                ui.label(f"Monitoring the path:{self.watchfolder}").style(
                    "color: #000000; font-size: 100%; font-weight: 300"
                )
                ui.label(f"Outputting to:{self.output}").style(
                    "color: #000000; font-size: 100%; font-weight: 300"
                )
                ui.label().bind_text_from(
                    self.bam_count, "counter", backward=lambda n: f"BAM files seen: {n}"
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    self.bam_passed, "counter", backward=lambda n: f"BAM pass: {n}"
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    self.bam_failed, "counter", backward=lambda n: f"BAM fail: {n}"
                ).style("color: #000000; font-size: 100%; font-weight: 300")

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
            with ui.row():
                ui.label().bind_text_from(
                    self,
                    "devices",
                    backward=lambda n: f"Devices: {str(n)}",
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    self,
                    "basecall_models",
                    backward=lambda n: f"Basecall Models: {str(n)}",
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    self,
                    "flowcell_ids",
                    backward=lambda n: f"Flowcell IDs: {str(n)}",
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    self,
                    "min_start_time",
                    backward=lambda n: f"Run Start Time: {n}",
                ).style("color: #000000; font-size: 100%; font-weight: 300")
                ui.label().bind_text_from(
                    self,
                    "sample_ids",
                    backward=lambda n: f"Sample ID: {str(n)}",
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
                with ui.row().style("color: #000000; font-size: 100%; font-weight: 300"):
                    if "mgmt" not in self.exclude:
                        mgmt = ui.column()
            if "cnv" not in self.exclude:
                with ui.row().style("color: #000000; font-size: 100%; font-weight: 300"):
                    if "cnv" not in self.exclude:
                        cnvsummary = ui.column()
            if "fusion" not in self.exclude:
                with ui.row().style("color: #000000; font-size: 100%; font-weight: 300"):
                    if "fusion" not in self.exclude:
                        fusions = ui.column()
            if "coverage" not in self.exclude:
                with ui.row().style("color: #000000; font-size: 100%; font-weight: 300"):
                    if "coverage" not in self.exclude:
                        coverage = ui.column()

        with ui.card().style("width: 100%"):
            ui.label("Methylation Classifications").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            if "sturgeon" not in self.exclude:
                self.Sturgeon = Sturgeon_object(
                    self.threads,
                    self.output,
                    progress=True,
                    batch=True,
                    bamqueue=self.bamforsturgeon,
                    summary=sturgeonsummary,
                )
            if "nanodx" not in self.exclude:
                self.NanoDX = NanoDX_object(
                    self.threads,
                    self.output,
                    progress=True,
                    batch=True,
                    bamqueue=self.bamfornanodx,
                    summary=nanodxsummary,
                )
            if "forest" not in self.exclude:
                self.RandomForest = RandomForest_object(
                    self.threads,
                    self.output,
                    progress=True,
                    batch=True,
                    bamqueue=self.bamforcns,
                    summary=forestsummary,
                    showerrors=self.showerrors,
                )
            pass
        #    with ui.tab_panel(copy_numer).classes("w-full"):
        if "cnv" not in self.exclude:
            with ui.card().style("width: 100%"):
                self.CNV = CNVAnalysis(
                    self.threads,
                    self.output,
                    progress=True,
                    bamqueue=self.bamforcnv,
                    summary=cnvsummary,
                )
            pass
        #    with ui.tab_panel(coverage).classes("w-full"):
        if "coverage" not in self.exclude:
            with ui.card().style("width: 100%"):
                self.Target_Coverage = TargetCoverage(
                    self.threads,
                    self.output,
                    progress=True,
                    bamqueue=self.bamfortargetcoverage,
                    summary=coverage,
                )
            pass
        #    with ui.tab_panel(mgmt).classes("w-full"):
        if "mgmt" not in self.exclude:
            with ui.card().style("width: 100%"):
                self.MGMT_panel = MGMT_Object(
                    self.threads,
                    self.output,
                    progress=True,
                    bamqueue=self.bamformgmt,
                    summary=mgmt,
                )
            pass
        #    with ui.tab_panel(fusions).classes("w-full"):
        if "fusion" not in self.exclude:
            with ui.card().style("width: 100%"):
                self.Fusion_panel = Fusion_object(
                    self.threads,
                    self.output,
                    progress=True,
                    bamqueue=self.bamforfusions,
                    summary=fusions,
                )
            pass

    def process_bams(self) -> None:
        """
        This function processes the bam files and adds them to the bamforcns and bamforsturgeon lists.
        These lists are then processed by the rapid_cns2 and sturgeon functions.
        #ToDo: switch to using threadsafe queues.
        :param self:
        :return:
        """
        while True:
            if "file" in self.bam_count:
                while len(self.bam_count["file"]) > 0:
                    self.nofiles = False
                    file = self.bam_count["file"].popitem()
                    if file[1] > time.time() - 5:
                        time.sleep(5)
                    self.check_bam(file[0])
                    if "forest" not in self.exclude:
                        self.bamforcns.put([file[0], file[1]])
                    if "sturgeon" not in self.exclude:
                        self.bamforsturgeon.put([file[0], file[1]])
                    if "nanodx" not in self.exclude:
                        self.bamfornanodx.put([file[0], file[1]])
                    if "cnv" not in self.exclude:
                        self.bamforcnv.put([file[0], file[1]])
                    if "coverage" not in self.exclude:
                        self.bamfortargetcoverage.put([file[0], file[1]])
                    if "mgmt" not in self.exclude:
                        self.bamformgmt.put([file[0], file[1]])
                    if "fusion" not in self.exclude:
                        self.bamforfusions.put([file[0], file[1]])

                time.sleep(1)
                self.nofiles = True
            time.sleep(1)

    def check_bam(self, bamfile):
        """
        This function checks a bam file and returns a dictionary of its attributes.
        :param bamfile:
        :return:
        """
        baminfo = ReadBam(bamfile).process_reads()
        if baminfo["state"] == "pass":
            self.bam_passed["counter"] += 1
        else:
            self.bam_failed["counter"] += 1
        self.basecall_models.add(baminfo["basecall_model"])
        self.devices.add(baminfo["device_position"])
        self.sample_ids.add(baminfo["sample_id"])
        self.flowcell_ids.add(baminfo["flow_cell_id"])
        self.run_time.add(baminfo["time_of_run"])

    def check_existing_bams(self, sequencing_summary=None):
        file_endings = {".bam"}
        if sequencing_summary:

            with self.frontpage:
                ui.notify("Checking for existing bams against sequencing summary")
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
                        self.check_bam(os.path.join(path, f))
            with self.frontpage:
                ui.notify("Bams Checked")

            step_size = 20

            with self.frontpage:
                ui.notify("Target Playback Started")

            if "sturgeon" not in self.exclude:
                self.Sturgeon.playback(latest_timestamps, step_size=step_size)
            if "nanodx" not in self.exclude:
                self.NanoDX.playback(latest_timestamps, step_size=step_size)
            if "forest" not in self.exclude:
                self.RandomForest.playback(latest_timestamps, step_size=step_size)
            if "cnv" not in self.exclude:
                self.CNV.playback(latest_timestamps, step_size=step_size)
            if "coverage" not in self.exclude:
                self.Target_Coverage.playback(latest_timestamps, step_size=step_size)
            if "fusion" not in self.exclude:
                self.Fusion_panel.playback(latest_timestamps, step_size=step_size)
            if "mgmt" not in self.exclude:
                self.MGMT_panel.playback(latest_timestamps, step_size=step_size)

            self.runfinished = True
        else:
            for path, dirs, files in os.walk(self.watchfolder):
                for f in files:
                    if "".join(Path(f).suffixes) in file_endings:
                        self.bam_count["counter"] += 1
                        if "file" not in self.bam_count:
                            self.bam_count["file"] = {}
                        self.bam_count["file"][os.path.join(path, f)] = time.time()
                        if self.simtime:
                            time.sleep(30)
