from nicegui import ui

from cnsmeth.utilities.bam_handler import BamEventHandler

# from cnsmeth.Sturgeon_worker import Sturgeon_worker
from cnsmeth.RCNS2_worker import RCNS2_worker

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

from watchdog.observers import Observer
from pathlib import Path
from queue import Queue
import pandas as pd

import threading
import time
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
    ):
        self.threads = threads
        self.simtime = simtime
        self.watchfolder = watchfolder
        self.output = output
        self.sequencing_summary = sequencing_summary
        self.showerrors = showerrors
        self.browse = browse

        self.bam_count = {"counter": 0}
        self.bamforcns = Queue()
        self.bamforsturgeon = Queue()
        self.bamfornanodx = Queue()
        self.bamforcnv = Queue()
        self.bamfortargetcoverage = Queue()
        self.bamformgmt = Queue()
        self.bamforfusions = Queue()

        # self.cnv = CNV_Plot()
        # self.mgmt_panel = MGMT_Panel()
        # self.fusion_panel = Fusion_Panel()

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
            # self.sturgeon_worker = Sturgeon_worker(
            #    self.bamforsturgeon, threads=self.threads, output_folder=self.output
            # )

            # self.nanodx_worker = CrossNN_worker(
            #    self.bamfornanodx, threads=self.threads, output_folder=self.output
            # )

            # self.rcns2_worker = RCNS2_worker(
            #    self.bamforcns,
            #    self.cnv,
            #    #self.mgmt_panel,
            #    #self.fusion_panel,
            #    threads=self.threads,
            #    output_folder=self.output,
            #    showerrors=self.showerrors,
            #    browse=self.browse,
            # )

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

            # self.rcns2_worker.replay_prior_data()

    def replay(self):
        ui.notify("Replaying data")
        self.replaycontrol.visible = False
        self.rcns2_worker.replay_prior_data()
        self.sturgeon_worker.replay_prior_data()

    def replay_cnv(self):
        ui.notify("Replaying CNV data")

    def information_panel(self):
        self.frontpage = ui.card().style("width: 100%")
        with self.frontpage:
            ui.label(
                "This tool enables classification of brain tumours in real time from Oxford Nanopore Data."
            ).tailwind("drop-shadow", "font-bold")
            with ui.row():
                ui.label(f"Monitoring the path:{self.watchfolder}").tailwind(
                    "drop-shadow"
                )
                ui.label(f"Outputting to:{self.output}").tailwind("drop-shadow")
                ui.label().bind_text_from(
                    self.bam_count, "counter", backward=lambda n: f"BAM files seen: {n}"
                ).tailwind("drop-shadow")
                ui.label().bind_text_from(
                    self,
                    "bamforcns",
                    backward=lambda n: f"BAM files for CNS: {n.qsize()}",
                ).tailwind("drop-shadow")
                ui.label().bind_text_from(
                    self,
                    "bamforsturgeon",
                    backward=lambda n: f"BAM files for Sturgeon: {n.qsize()}",
                ).tailwind("drop-shadow")
                ui.label().bind_text_from(
                    self,
                    "bamfornanodx",
                    backward=lambda n: f"BAM files for NanoDX: {n.qsize()}",
                ).tailwind("drop-shadow")
        # with ui.tabs().classes("w-full") as tabs:
        #    methylation = ui.tab("Methylation Classification")
        #    copy_numer = ui.tab("Copy Number Variation")
        #    coverage = ui.tab("Target Coverage")
        #    mgmt = ui.tab("MGMT")
        #    fusions = ui.tab("Fusions")
        # with ui.tab_panels(tabs, value=methylation).classes("w-full"):
        #    with ui.tab_panel(methylation).classes("w-full"):
        with ui.card().style("width: 100%"):
            self.Sturgeon = Sturgeon_object(self.threads, progress=True, batch=True, bamqueue=self.bamforsturgeon)
            self.NanoDX = NanoDX_object(self.threads, progress=True, batch=True, bamqueue=self.bamfornanodx)
            self.RandomForest = RandomForest_object(self.threads, progress=True, batch=True, bamqueue=self.bamforcns)
            pass
        #    with ui.tab_panel(copy_numer).classes("w-full"):
        with ui.card().style("width: 100%"):
            self.CNV = CNVAnalysis(self.threads, progress=True, bamqueue=self.bamforcnv)
            pass

        #    with ui.tab_panel(coverage).classes("w-full"):
        with ui.card().style("width: 100%"):
            self.Target_Coverage = TargetCoverage(self.threads, progress=True,bamqueue=self.bamfortargetcoverage)
            pass

        #    with ui.tab_panel(mgmt).classes("w-full"):
        with ui.card().style("width: 100%"):
            self.MGMT_panel = MGMT_Object(self.threads, progress=True, bamqueue=self.bamformgmt)
            pass

        #    with ui.tab_panel(fusions).classes("w-full"):
        with ui.card().style("width: 100%"):
            self.Fusion_panel = Fusion_object(self.threads, progress=True, bamqueue=self.bamforfusions)
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
                    self.bamforcns.put([file[0],file[1]])
                    self.bamforsturgeon.put([file[0],file[1]])
                    self.bamfornanodx.put([file[0],file[1]])
                    self.bamforcnv.put([file[0],file[1]])
                    self.bamfortargetcoverage.put([file[0],file[1]])
                    self.bamformgmt.put([file[0],file[1]])
                    self.bamforfusions.put([file[0],file[1]])

                time.sleep(1)
                self.nofiles = True
            time.sleep(1)

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
            with self.frontpage:
                ui.notify("Bams Checked")

            step_size = 20

            with self.frontpage:
                ui.notify("Target Playback Started")

            self.Sturgeon.playback(latest_timestamps,step_size=step_size)
            self.NanoDX.playback(latest_timestamps, step_size=step_size)
            self.RandomForest.playback(latest_timestamps, step_size=step_size)
            self.CNV.playback(latest_timestamps, step_size=step_size)
            self.Target_Coverage.playback(latest_timestamps, step_size=step_size)
            self.Fusion_panel.playback(latest_timestamps, step_size=step_size)
            self.MGMT_panel.playback(latest_timestamps,step_size=step_size)

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
