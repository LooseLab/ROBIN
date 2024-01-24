from nicegui import Tailwind, ui, app
from nicegui.events import ValueChangeEventArguments
from io import StringIO

from methnicegui.bam_handler import BamEventHandler
from methnicegui.Sturgeon_worker import Sturgeon_worker
from methnicegui.RCNS2_worker import RCNS2_worker
from methnicegui.copy_number_component import CNV_Plot
from methnicegui.target_coverage import TargetCoverage
from methnicegui.mgmt_panel import MGMT_Panel
from methnicegui.local_file_picker import local_file_picker

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

        print("BrainMeth Loaded")
        self.cnv = CNV_Plot()
        self.target_coverage = TargetCoverage()
        self.mgmt_panel = MGMT_Panel()

        if not self.browse:
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
            self.sturgeon_worker = Sturgeon_worker(
                self.bamforsturgeon, threads=self.threads, output_folder=self.output
            )

            self.rcns2_worker = RCNS2_worker(
                self.bamforcns,
                self.cnv,
                self.target_coverage,
                self.mgmt_panel,
                threads=self.threads,
                output_folder=self.output,
                showerrors=self.showerrors,
                browse=self.browse,
            )

            self.information_panel()
        else:
            ui.label("Browse mode enabled. Please choose a folder to see data from.")
            ui.button("Choose file", on_click=self.pick_file, icon="folder")

            self.content = ui.column().classes("w-full")

    async def pick_file(self) -> None:
        result = await local_file_picker("~", multiple=True)
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
            with ui.tab_panels(tabs, value=methylation).classes("w-full"):
                with ui.tab_panel(methylation).classes("w-full"):
                    with ui.card().classes("w-full"):
                        with ui.column().classes("w-full"):
                            self.rcns2_worker = RCNS2_worker(
                                self.bamforcns,
                                self.cnv,
                                self.target_coverage,
                                self.mgmt_panel,
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
                                with ui.grid(columns=2).classes("w-full h-auto"):
                                    with ui.column():
                                        with ui.card().style("width: 100%"):
                                            self.rcns2_worker.status_rcns2()
                                            self.rcns2_worker.create_rcns2_chart(
                                                "RapidCNS2"
                                            )

                                    with ui.column():
                                        with ui.card().style("width: 100%"):
                                            self.sturgeon_worker.status_sturgeon()
                                            self.sturgeon_worker.create_sturgeon_chart(
                                                "Sturgeon"
                                            )

                            with ui.card().style("width: 100%"):
                                self.rcns2_worker.create_rcns2_time_chart()

                            with ui.card().style("width: 100%"):
                                self.sturgeon_worker.create_sturgeon_time_chart()
                            ui.button(
                                "replay data", on_click=self.replay, icon="replay"
                            )

                with ui.tab_panel(copy_numer).classes("w-full"):
                    with ui.card().style("width: 100%"):
                        self.cnv.create_cnv_scatter("CNV Scatter")
                        pass
                    ui.button(
                        "replay cnv data", on_click=self.replay_cnv, icon="replay"
                    )

                with ui.tab_panel(coverage).classes("w-full"):
                    self.target_coverage.setup_ui()

                with ui.tab_panel(mgmt).classes("w-full"):
                    self.mgmt_panel.setup_ui(mgmt)

                self.rcns2_worker.load_prior_data()
                self.sturgeon_worker.load_prior_data()

                # self.rcns2_worker.replay_prior_data()

    def replay(self):
        self.rcns2_worker.replay_prior_data()
        self.sturgeon_worker.replay_prior_data()

    def replay_cnv(self):
        ui.notify("Replaying CNV data")

    def information_panel(self):
        with ui.card().style("width: 100%"):
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
        with ui.tabs().classes("w-full") as tabs:
            methylation = ui.tab("Methylation Classification")
            copy_numer = ui.tab("Copy Number Variation")
            coverage = ui.tab("Target Coverage")
            mgmt = ui.tab("MGMT")
        with ui.tab_panels(tabs, value=methylation).classes("w-full"):
            with ui.tab_panel(methylation).classes("w-full"):
                with ui.card().style("width: 100%"):
                    with ui.grid(columns=2).classes("w-full h-auto"):
                        with ui.column():
                            with ui.card().style("width: 100%"):
                                self.rcns2_worker.status_rcns2()
                                self.rcns2_worker.create_rcns2_chart("RapidCNS2")

                        with ui.column():
                            with ui.card().style("width: 100%"):
                                self.sturgeon_worker.status_sturgeon()
                                self.sturgeon_worker.create_sturgeon_chart(
                                    "Sturgeon Prediction"
                                )

                with ui.card().style("width: 100%"):
                    self.sturgeon_worker.create_sturgeon_time_chart()

                with ui.card().style("width: 100%"):
                    self.rcns2_worker.create_rcns2_time_chart()

            with ui.tab_panel(copy_numer).classes("w-full"):
                with ui.card().style("width: 100%"):
                    self.cnv.create_cnv_scatter("CNV Scatter")

            with ui.tab_panel(coverage).classes("w-full"):
                self.target_coverage.setup_ui()

            with ui.tab_panel(mgmt).classes("w-full"):
                self.mgmt_panel.setup_ui(mgmt)

    def process_bams(self) -> None:
        """
        This function processes the bam files and adds them to the bamforcns and bamforsturgeon lists.
        These lists are then processed by the rapid_cns2 and sturgeon functions.
        #ToDo: switch to using threadsafe queues.
        :param self:
        :return:
        """

        # First we look for existing files:
        # self.check_existing_bams()

        while True:
            # self.log("Processing bam files")
            # self.log(self.bam_count)
            # self.bam_count = self.bam_count
            if "file" in self.bam_count:
                while len(self.bam_count["file"]) > 0:
                    self.nofiles = False
                    file = self.bam_count["file"].popitem()
                    if file[1] > time.time() - 5:
                        time.sleep(5)
                    self.bamforcns.put(file[0])
                    self.bamforsturgeon.put(file[0])
                time.sleep(1)
                self.nofiles = True
            time.sleep(1)

    def check_existing_bams(self, sequencing_summary=None):
        file_endings = {".bam"}
        if sequencing_summary:
            # get the current time
            now = time.time()
            print("Checking for existing bams against sequencing summary")
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
            latest_timestamps["file_produced"] = latest_timestamps["template_end"] + now
            print(latest_timestamps)
            for path, dirs, files in os.walk(self.watchfolder):
                for f in files:
                    if "".join(Path(f).suffixes) in file_endings:
                        latest_timestamps.loc[
                            latest_timestamps["filename_bam"] == f, "full_path"
                        ] = os.path.join(path, f)
                        # latest_timestamps[latest_timestamps['filename_bam'] == f]['full_path'] = os.path.join(path, f)
            print(latest_timestamps)
            for index, row in latest_timestamps.iterrows():
                current_time = time.time()
                time_diff = row["file_produced"] - current_time
                print(f"Sleeping for {time_diff} seconds.")
                if time_diff > 0:
                    time.sleep(time_diff)
                self.bam_count["counter"] += 1
                if "file" not in self.bam_count:
                    self.bam_count["file"] = {}
                self.bam_count["file"][row["full_path"]] = time.time()
                # print (index,row)
            self.runfinished = True
            # os.kill(os.getpid(), signal.SIGINT)
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
