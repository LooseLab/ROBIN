"""
This file provides a base class for analysis of bam files output during a sequencing run.
The base class provides a queue to receive bam files and a background thread to process the data.
"""
import queue
from nicegui import ui, run
from typing import BinaryIO
import pandas as pd
import time
import asyncio
import threading

class BaseAnalysis:
    def __init__(self, progress=False, batch=False, *args, **kwargs):
        self.queue = queue.Queue()
        self.batch = batch
        self.setup_ui()
        if progress:
            self.progress()
        self.bam_count = 0
        self.bam_processed = 0
        self.bams_in_processing = 0
        self.running = False

        if batch:
            self.bams = []
            self.batch_timer_run()
        else:
            self.timer_run()


    def timer_run(self):
        print("Running timer")
        timer = ui.timer(0.1,self._worker)

    def batch_timer_run(self):
        print("Running batch timer")
        print (type(self).__name__)
        ui.timer(1,self._batch_worker)

    def add_bam(self, bamfile: BinaryIO, timestamp=None):
        """
        Adds a bam file to the queue for processing.
        :param bamfile: The path to the bam file to process.
        :return:
        """
        self.queue.put([bamfile, timestamp])
        self.bam_count += 1


    async def _batch_worker(self):
        """
        This function takes bam files from a queue in batches and adds them to a backround thread for processing.
        """
        #print("Running batch worker")
        while self.queue.qsize() > 0:
            self.bams.append((self.queue.get()))
            self.bams_in_processing+=1
        if not self.running and len(self.bams) > 0:
            await self.process_bam(self.bams)
        else:
            await asyncio.sleep(1)


    async def _worker(self):
        """
        This function takes reads from the queue and adds them to the background thread for processing.
        """
        # Indicate that the thread is running

        if not self.queue.empty() and not self.running:
            self.running = True
            bamfile, timestamp = self.queue.get()
            self.process_bam(bamfile, timestamp)
            self.bam_processed += 1
        else:
            await asyncio.sleep(1)
        #self.running = False


    @property
    def _progress(self):
        if self.bam_count == 0:
            return 0
        else:
            return ((self.bam_count - self.queue.qsize()-self.bams_in_processing) / self.bam_count)

    @property
    def _progress2(self):
        if self.bam_count == 0:
            return 0
        else:
            return ((self.bams_in_processing) / self.bam_count)



    def progress(self):
        """
        Show a progress bar for the number of bam files processed.
        :return:
        """
        with ui.card().classes("w-full"):
            with ui.row():
                ui.label("File Tracker").tailwind("drop-shadow", "font-bold")
                ui.label().bind_text_from(self, 'bam_count', backward=lambda n: f'Bam files seen: {n}')
                if self.batch:
                    ui.label().bind_text_from(self, 'bams_in_processing', backward=lambda n: f'Bam files being processed: {n}')
                ui.label().bind_text_from(self, 'bam_processed', backward=lambda n: f'Bam files processed: {n}')
            if self.batch:
                ui.timer(0.1, callback=lambda: progressbar2.set_value(self._progress2))
                progressbar2 = ui.linear_progress(size="10px", show_value=False, value=0, color="orange").props('instant-feedback')
            ui.timer(0.1, callback=lambda: progressbar.set_value(self._progress))
            progressbar = ui.linear_progress(size="10px", show_value=False, value=0).props('instant-feedback')

    def playback(self, data: pd.DataFrame, step_size=2):
        self.data = data
        playback = threading.Thread(target=self.playback_bams, args=([step_size]))
        playback.daemon = True
        playback.start()

    def playback_bams(self, step_size=2):
        latest_timestamps = self.data
        self.offset = 0
        for index, row in latest_timestamps.iterrows():
            current_time = time.time()
            time_diff = row["file_produced"] - current_time - self.offset
            if self.offset == 0:
                self.offset = time_diff
                time_diff = 0
            while row["file_produced"] - current_time - self.offset > 0:
                time.sleep(0.1)
                if not self.running:
                    self.offset += step_size
            self.add_bam(row["full_path"], row["file_produced"])

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



