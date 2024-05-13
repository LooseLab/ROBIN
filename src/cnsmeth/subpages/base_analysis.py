"""
This file provides a base class for analysis of bam files output during a sequencing run.
The base class provides a queue to receive bam files and a background thread to process the data.
"""

import queue
from nicegui import ui
from typing import BinaryIO
import pandas as pd
import time
import asyncio
import threading


class BaseAnalysis:
    def __init__(
        self,
        threads,
        outputfolder,
        summary=None,
        bamqueue=None,
        progress=False,
        batch=False,
        start_time=None,
        browse=False,
        *args,
        **kwargs,
    ):
        if bamqueue:
            self.bamqueue = bamqueue
        else:
            self.bamqueue = queue.Queue()
        self.start_time = start_time
        self.batch = batch
        self.output = outputfolder
        self.summary = summary
        self.browse = browse
        self.setup_ui()
        if progress and not self.browse:
            self.progress()
        self.bam_count = 0
        self.bam_processed = 0
        self.bams_in_processing = 0
        self.running = False
        if threads > 1:
            self.threads = int(threads / 2)
        else:
            self.threads = threads
        if batch:
            self.bams = []
            self.batch_timer_run()
        else:
            self.timer_run()

    def timer_run(self):
        self.timer = ui.timer(0.1, self._worker)

    async def _worker(self):
        """
        This function takes reads from the queue and adds them to the background thread for processing.
        """
        self.timer.active = False
        if not self.bamqueue.empty() and not self.running:
            self.running = True
            bamfile, timestamp = self.bamqueue.get()
            self.bam_count += 1
            if not timestamp:
                timestamp = time.time()
            # print (bamfile,timestamp)
            await self.process_bam(bamfile, timestamp)
            self.bam_processed += 1
        else:
            await asyncio.sleep(1)
        self.timer.active = True

    def add_bam(self, bamfile: BinaryIO, timestamp=None):
        """
        Adds a bam file to the queue for processing.
        :param bamfile: The path to the bam file to process.
        :return:
        """
        self.bamqueue.put([bamfile, timestamp])
        # self.bam_count += 1

    def batch_timer_run(self):
        self.timer = ui.timer(1, self._batch_worker)

    async def _batch_worker(self):
        """
        This function takes bam files from a queue in batches and adds them to a background thread for processing.
        """
        self.timer.active = False
        while self.bamqueue.qsize() > 0:
            self.bams.append((self.bamqueue.get()))
            self.bam_count += 1
            # self.bams_in_processing += 1
        if not self.running and len(self.bams) > 0:
            self.running = True
            await self.process_bam(self.bams)
        else:
            await asyncio.sleep(1)
        self.timer.active = True

    @property
    def _progress(self):
        """
        This property generates a progress bar indicating the number of files that have been successfully processed.
        """
        if self.bam_count == 0:
            return 0
        return (self.bam_processed) / self.bam_count

    @property
    def _progress2(self):
        if self.bam_count == 0:
            return 0
        return (self.bams_in_processing) / self.bam_count

    @property
    def _not_analysed(self):
        if self.bam_count == 0:
            return 0
        return (
            self.bam_count - self.bams_in_processing - self.bam_processed
        ) / self.bam_count

    def progress(self):
        """
        Show a progress bar for the number of bam files processed.
        :return:
        """
        self.progress_trackers = ui.card().classes("w-full")
        with self.progress_trackers:
            with ui.row():
                ui.label("File Tracker").tailwind("drop-shadow", "font-bold")
                ui.label().bind_text_from(
                    self, "bam_count", backward=lambda n: f"Bam files seen: {n}"
                )
                if self.batch:
                    ui.label().bind_text_from(
                        self,
                        "bams_in_processing",
                        backward=lambda n: f"Bam files being processed: {n}",
                    )
                ui.label().bind_text_from(
                    self,
                    "bam_processed",
                    backward=lambda n: f"Bam files processed: {n}",
                )

            progressbar3 = (
                ui.linear_progress(
                    size="10px",
                    show_value=False,
                    value=0,
                    color="red",
                )
                .tooltip("Indicates read files not yet processed.")
                .props("instant-feedback")
            )
            ui.timer(1, callback=lambda: progressbar3.set_value(self._not_analysed))
            if self.batch:
                progressbar2 = (
                    ui.linear_progress(
                        size="10px",
                        show_value=False,
                        value=0,
                        color="amber",
                    )
                    .tooltip("Indicates read files being processed.")
                    .props("instant-feedback")
                )
                ui.timer(1, callback=lambda: progressbar2.set_value(self._progress2))

            progressbar = (
                ui.linear_progress(
                    size="10px",
                    show_value=False,
                    value=0,
                    color="green",
                )
                .tooltip("Indicates read files processed.")
                .props("instant-feedback")
            )
            ui.timer(1, callback=lambda: progressbar.set_value(self._progress))

    def playback(self, data: pd.DataFrame, step_size=2, start_time=None):
        self.data = data
        playback = threading.Thread(
            target=self.playback_bams, args=([step_size, start_time])
        )
        playback.daemon = True
        playback.start()

    def playback_bams(self, step_size=2, start_time=None):
        """
        This function plays back the processing of bam files from a pandas dataframe.
        To simulate the behaviour of a run it adds files to the queue in the order in which they were produced.
        The run time for each individual processing loop is monitored and the reads are added to the queue
        accounting for this delay.
        The delay is calculated by the offset parameter which starts as 0.

        """
        latest_timestamps = self.data
        self.offset = 0
        # print (start_time)
        playback_start_time = time.time()
        for index, row in latest_timestamps.iterrows():
            elapsed_time = (time.time() - playback_start_time) + self.offset
            time_diff = row["file_produced"]  # - elapsed_time
            # print (elapsed_time, time_diff, row["file_produced"])
            if self.offset == 0:
                self.offset = time_diff
                time_diff = 0
            while elapsed_time < time_diff:
                if self.running:
                    time.sleep(1)
                    elapsed_time = (time.time() - playback_start_time) + self.offset
                else:
                    time.sleep(1)
                    self.offset += step_size
                    elapsed_time += self.offset
            # print("out the loop")
            # print (f"elapsed time now {elapsed_time}")
            if len(row["full_path"]) > 0:
                # Here we check that the path to the bam file absolutely exists.
                self.add_bam(
                    row["full_path"], playback_start_time + row["file_produced"]
                )

    def create_time_chart(self):
        """
        This function creates a time chart for any module which needs it.
        :return: an echarts object
        """
        return (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": "NanoDX Over Time"},
                    "toolbox": {"show": True, "feature": {"dataZoom": {"yAxisIndex": "none"},"restore": {},"saveAsImage": {}}},
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    "series": [],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def create_chart(self, title):
        """
        This function creates a time chart for any module which needs it.
        :return: an echarts object
        """
        return (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "value", "max": 1},
                    "yAxis": {"type": "category", "data": [], "inverse": True},
                    #'legend': {},
                    "series": [],
                }
            )
            .style("color: #6E93D6; font-size: 150%; font-weight: 300; height: 350px")
            .classes("border-double")
        )

    def process_bam(self, bamfile, timestamp):
        """
        Process a bam file.
        :param bamfile: The path to the bam file to process.
        :param timestamp: A time stamp indicating when the file was generated.
        :return:
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def setup_ui(self):
        """
        Set up the user interface for the analysis.
        This function should be overridden by subclasses to create the user interface for the analysis.
        It should create the necessary widgets and add them to the display.
        The function should do no work of its own. That should all be done in the process_bam method.
        The function may call other methods to trigger updates to the display.
        :return:
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def cleanup(self):
        """
        This function is called when the app is closed.
        It should be overridden by subclasses to clean up any resources used by the analysis.
        :return:
        """
        pass

    def check_resources(self):
        """
        This function is called to check the resources required for the analysis are present.
        It should be overridden by subclasses to check the resources required for the analysis are present.
        :return:
        """
        pass
