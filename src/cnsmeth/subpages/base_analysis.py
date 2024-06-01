# """
# This file provides a base class for analysis of bam files output during a sequencing run.
# The base class provides a queue to receive bam files and a background thread to process the data.
# The class also provides a render pipeline. These are separate. The background thread should only be called in
# one class instance. Rendering can be called in as many classes as make connections.
# """
#
# import queue
# from nicegui import ui, app
# from typing import BinaryIO
# import pandas as pd
# import time
# import asyncio
# import threading
# from collections import Counter
#
#
# class BaseAnalysis:
#     def __init__(
#         self,
#         threads,
#         outputfolder,
#         analysis_name,
#         summary=None,
#         bamqueue=None,
#         progress=False,
#         batch=False,
#         start_time=None,
#         browse=False,
#         uuid=None,
#         *args,
#         **kwargs,
#     ):
#         # Create a unique ID that we can use to register this specific item.
#         self.mainuuid = uuid
#         if bamqueue:
#             self.bamqueue = bamqueue
#         else:
#             self.bamqueue = queue.Queue()
#         self.name = analysis_name
#         self.start_time = start_time
#         self.batch = batch
#         self.output = outputfolder
#         self.summary = summary
#         self.browse = browse
#         self.progress = progress
#         if self.name not in app.storage.general[self.mainuuid]:
#             app.storage.general[self.mainuuid][self.name] = {}
#             app.storage.general[self.mainuuid][self.name]["counters"] = Counter(
#                 bam_count=0, bam_processed=0, bams_in_processing=0
#             )
#         # self.bam_count = 0
#         # self.bam_processed = 0
#         # self.bams_in_processing = 0
#         self.running = False
#         if threads > 1:
#             self.threads = int(threads / 2)
#         else:
#             self.threads = threads
#
#     def render_ui(self):
#         self.setup_ui()
#         if self.progress and not self.browse:
#             self.progress_bars()
#
#     def process_data(self):
#         if self.batch:
#             self.bams = []
#             self.batch_timer_run()
#         else:
#             self.timer_run()
#
#     def timer_run(self):
#         self.timer = ui.timer(1, self._worker)
#
#     async def _worker(self):
#         """
#         This function takes reads from the queue and adds them to the background thread for processing.
#         """
#         self.timer.active = False
#         if not self.bamqueue.empty() and not self.running:
#             self.running = True
#             bamfile, timestamp = self.bamqueue.get()
#             app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"] += 1
#
#             if not timestamp:
#                 timestamp = time.time()
#             await self.process_bam(bamfile, timestamp)
#             app.storage.general[self.mainuuid][self.name]["counters"][
#                 "bam_processed"
#             ] += 1
#         #else:
#         #    await asyncio.sleep(1)
#         self.timer.active = True
#
#     def add_bam(self, bamfile: BinaryIO, timestamp=None):
#         """
#         Adds a bam file to the queue for processing.
#         :param bamfile: The path to the bam file to process.
#         :return:
#         """
#         self.bamqueue.put([bamfile, timestamp])
#         # self.bam_count += 1
#
#     def batch_timer_run(self):
#         self.timer = ui.timer(5, self._batch_worker)
#
#     async def _batch_worker(self):
#         """
#         This function takes bam files from a queue in batches and adds them to a background thread for processing.
#         """
#         self.timer.active = False
#         count = 0
#         while self.bamqueue.qsize() > 0:
#             self.bams.append((self.bamqueue.get()))
#             count += 1
#             if count >= 200:
#                 break
#         app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"] += count
#
#         if not self.running and len(self.bams) > 0:
#             self.running = True
#             await self.process_bam(self.bams)
#         # else:
#         # await asyncio.sleep(1)
#         self.timer.active = True
#
#     @property
#     def _progress(self):
#         """
#         This property generates a progress bar indicating the number of files that have been successfully processed.
#         """
#         if app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"] == 0:
#             return 0
#         return (
#             app.storage.general[self.mainuuid][self.name]["counters"]["bam_processed"]
#         ) / app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"]
#
#     @property
#     def _progress2(self):
#         if app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"] == 0:
#             return 0
#         return (
#             app.storage.general[self.mainuuid][self.name]["counters"][
#                 "bams_in_processing"
#             ]
#         ) / app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"]
#
#     @property
#     def _not_analysed(self):
#         if app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"] == 0:
#             return 0
#         return (
#             app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"]
#             - app.storage.general[self.mainuuid][self.name]["counters"][
#                 "bams_in_processing"
#             ]
#             - app.storage.general[self.mainuuid][self.name]["counters"]["bam_processed"]
#         ) / app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"]
#
#     def progress_bars(self):
#         """
#         Show a progress bar for the number of bam files processed.
#         :return:
#         """
#         self.progress_trackers = ui.card().classes("w-full")
#         with self.progress_trackers:
#             with ui.row():
#                 ui.label("File Tracker").tailwind("drop-shadow", "font-bold")
#                 ui.label().bind_text_from(
#                     app.storage.general[self.mainuuid][self.name]["counters"],
#                     "bam_count",
#                     backward=lambda n: f"Bam files seen: {n}",
#                 )
#                 if self.batch:
#                     ui.label().bind_text_from(
#                         app.storage.general[self.mainuuid][self.name]["counters"],
#                         "bams_in_processing",
#                         backward=lambda n: f"Bam files being processed: {n}",
#                     )
#                 ui.label().bind_text_from(
#                     app.storage.general[self.mainuuid][self.name]["counters"],
#                     "bam_processed",
#                     backward=lambda n: f"Bam files processed: {n}",
#                 )
#
#             progressbar3 = (
#                 ui.linear_progress(
#                     size="10px",
#                     show_value=False,
#                     value=0,
#                     color="red",
#                 )
#                 .tooltip("Indicates read files not yet processed.")
#                 .props("instant-feedback")
#             )
#             ui.timer(1, callback=lambda: progressbar3.set_value(self._not_analysed))
#             if self.batch:
#                 progressbar2 = (
#                     ui.linear_progress(
#                         size="10px",
#                         show_value=False,
#                         value=0,
#                         color="amber",
#                     )
#                     .tooltip("Indicates read files being processed.")
#                     .props("instant-feedback")
#                 )
#                 ui.timer(1, callback=lambda: progressbar2.set_value(self._progress2))
#
#             progressbar = (
#                 ui.linear_progress(
#                     size="10px",
#                     show_value=False,
#                     value=0,
#                     color="green",
#                 )
#                 .tooltip("Indicates read files processed.")
#                 .props("instant-feedback")
#             )
#             ui.timer(1, callback=lambda: progressbar.set_value(self._progress))
#
#     def playback(self, data: pd.DataFrame, step_size=2, start_time=None):
#         self.data = data
#         playback = threading.Thread(
#             target=self.playback_bams, args=([step_size, start_time])
#         )
#         playback.daemon = True
#         playback.start()
#
#     def playback_bams(self, step_size=2, start_time=None):
#         """
#         This function plays back the processing of bam files from a pandas dataframe.
#         To simulate the behaviour of a run it adds files to the queue in the order in which they were produced.
#         The run time for each individual processing loop is monitored and the reads are added to the queue
#         accounting for this delay.
#         The delay is calculated by the offset parameter which starts as 0.
#
#         """
#         latest_timestamps = self.data
#         self.offset = 0
#         playback_start_time = time.time()
#         for index, row in latest_timestamps.iterrows():
#             elapsed_time = (time.time() - playback_start_time) + self.offset
#             time_diff = row["file_produced"]  # - elapsed_time
#             if self.offset == 0:
#                 self.offset = time_diff
#                 time_diff = 0
#             while elapsed_time < time_diff:
#                 if self.running:
#                     time.sleep(1)
#                     elapsed_time = (time.time() - playback_start_time) + self.offset
#                 else:
#                     time.sleep(1)
#                     self.offset += step_size
#                     elapsed_time += self.offset
#             if len(row["full_path"]) > 0:
#                 # Here we check that the path to the bam file absolutely exists.
#                 self.add_bam(
#                     row["full_path"], playback_start_time + row["file_produced"]
#                 )
#
#     def create_time_chart(self, title):
#         """
#         This function creates a time chart for any module which needs it.
#         :return: an echarts object
#         """
#         return (
#             ui.echart(
#                 {
#                     "animation": False,
#                     "grid": {"containLabel": True},
#                     "title": {"text": title},
#                     "toolbox": {
#                         "show": True,
#                         "feature": {
#                             "dataZoom": {"yAxisIndex": "none"},
#                             "restore": {},
#                             "saveAsImage": {},
#                         },
#                     },
#                     "xAxis": {"type": "time"},
#                     "yAxis": {"type": "value", "data": [], "inverse": False},
#                     "series": [],
#                 }
#             )
#             .style("height: 350px")
#             .classes("border-double")
#         )
#
#     def create_chart(self, title):
#         """
#         This function creates a time chart for any module which needs it.
#         :return: an echarts object
#         """
#         return (
#             ui.echart(
#                 {
#                     "animation": False,
#                     "grid": {"containLabel": True},
#                     "title": {"text": title},
#                     "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
#                     "xAxis": {"type": "value", "max": 1},
#                     "yAxis": {"type": "category", "data": [], "inverse": True},
#                     #'legend': {},
#                     "series": [],
#                 }
#             )
#             .style("color: #6E93D6; font-size: 150%; font-weight: 300; height: 350px")
#             .classes("border-double")
#         )
#
#     def process_bam(self, bamfile, timestamp):
#         """
#         Process a bam file.
#         :param bamfile: The path to the bam file to process.
#         :param timestamp: A time stamp indicating when the file was generated.
#         :return:
#         """
#         raise NotImplementedError("Subclasses must implement this method.")
#
#     def setup_ui(self):
#         """
#         Set up the user interface for the analysis.
#         This function should be overridden by subclasses to create the user interface for the analysis.
#         It should create the necessary widgets and add them to the display.
#         The function should do no work of its own. That should all be done in the process_bam method.
#         The function may call other methods to trigger updates to the display.
#         :return:
#         """
#         raise NotImplementedError("Subclasses must implement this method.")
#
#     def cleanup(self):
#         """
#         This function is called when the app is closed.
#         It should be overridden by subclasses to clean up any resources used by the analysis.
#         :return:
#         """
#         pass
#
#     def check_resources(self):
#         """
#         This function is called to check the resources required for the analysis are present.
#         It should be overridden by subclasses to check the resources required for the analysis are present.
#         :return:
#         """
#         pass


"""
This file provides a base class for analysis of bam files output during a sequencing run.
The base class provides a queue to receive bam files and a background thread to process the data.
The class also provides a render pipeline. These are separate. The background thread should only be called in
one class instance. Rendering can be called in as many classes as make connections.
"""

import queue
from nicegui import ui, app
from typing import BinaryIO, Optional, List, Tuple
import pandas as pd
import time
import threading
from collections import Counter
import logging

# Configure logging
logger = logging.getLogger(__name__)


class BaseAnalysis:
    """
    A base class for analysis of BAM files output during a sequencing run.

    This class provides a queue to receive BAM files and a background thread to process the data.
    It also provides a render pipeline for displaying progress and results.

    Args:
        threads (int): Number of threads to use for processing.
        outputfolder (str): Directory to store output files.
        analysis_name (str): Name of the analysis.
        summary (Optional[str]): Summary of the analysis.
        bamqueue (Optional[queue.Queue]): Queue to hold BAM files.
        progress (bool): Flag to display progress bars.
        batch (bool): Flag to enable batch processing.
        start_time (Optional[float]): Start time of the analysis.
        browse (bool): Flag to enable browsing of results.
        uuid (Optional[str]): Unique identifier for the analysis instance.
    """

    def __init__(
        self,
        threads: int,
        outputfolder: str,
        analysis_name: str,
        summary: Optional[str] = None,
        bamqueue: Optional[queue.Queue] = None,
        progress: bool = False,
        batch: bool = False,
        start_time: Optional[float] = None,
        browse: bool = False,
        uuid: Optional[str] = None,
        *args,
        **kwargs,
    ) -> None:
        self.mainuuid = uuid
        self.bamqueue = bamqueue if bamqueue else queue.Queue()
        self.name = analysis_name
        self.start_time = start_time
        self.batch = batch
        self.output = outputfolder
        self.summary = summary
        self.browse = browse
        self.progress = progress
        if self.name not in app.storage.general.get(self.mainuuid, {}):
            app.storage.general.setdefault(self.mainuuid, {})[self.name] = {
                "counters": Counter(bam_count=0, bam_processed=0, bams_in_processing=0)
            }
        self.running = False
        self.threads = max(1, threads // 2)

    def render_ui(self) -> None:
        """
        Set up the user interface for the analysis and display progress bars if enabled.
        """
        self.setup_ui()
        if self.progress and not self.browse:
            self.progress_bars()

    def process_data(self) -> None:
        """
        Start processing BAM files either in batch mode or in a continuous timer mode.
        """
        if self.batch:
            self.bams: List[Tuple[BinaryIO, Optional[float]]] = []
            self.batch_timer_run()
        else:
            self.timer_run()

    def timer_run(self) -> None:
        """
        Set up a timer to periodically run the _worker method for processing BAM files.
        """
        self.timer = ui.timer(1, self._worker)

    async def _worker(self) -> None:
        """
        Process BAM files from the queue in the background.

        This function takes reads from the queue and adds them to the background thread for processing.
        """
        self.timer.active = False
        if not self.bamqueue.empty() and not self.running:
            self.running = True
            bamfile, timestamp = self.bamqueue.get()
            app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"] += 1

            if not timestamp:
                timestamp = time.time()
            await self.process_bam(bamfile, timestamp)
            app.storage.general[self.mainuuid][self.name]["counters"][
                "bam_processed"
            ] += 1

        self.timer.active = True

    def add_bam(self, bamfile: BinaryIO, timestamp: Optional[float] = None) -> None:
        """
        Adds a BAM file to the queue for processing.

        Args:
            bamfile (BinaryIO): The BAM file to process.
            timestamp (Optional[float]): The timestamp indicating when the file was generated.
        """
        self.bamqueue.put((bamfile, timestamp))

    def batch_timer_run(self) -> None:
        """
        Set up a timer to periodically run the _batch_worker method for batch processing of BAM files.
        """
        self.timer = ui.timer(5, self._batch_worker)

    async def _batch_worker(self) -> None:
        """
        Process BAM files from the queue in batches in the background.

        This function takes BAM files from a queue in batches and adds them to a background thread for processing.
        """
        self.timer.active = False
        count = 0
        while self.bamqueue.qsize() > 0:
            self.bams.append(self.bamqueue.get())
            count += 1
            if count >= 200:
                break
        app.storage.general[self.mainuuid][self.name]["counters"]["bam_count"] += count

        if not self.running and self.bams:
            self.running = True
            await self.process_bam(self.bams)

        self.timer.active = True

    @property
    def _progress(self) -> float:
        """
        Calculate the progress of BAM file processing.

        Returns:
            float: The progress as a fraction of processed files over total files.
        """
        counters = app.storage.general[self.mainuuid][self.name]["counters"]
        if counters["bam_count"] == 0:
            return 0.0
        return counters["bam_processed"] / counters["bam_count"]

    @property
    def _progress2(self) -> float:
        """
        Calculate the progress of BAM files currently being processed.

        Returns:
            float: The progress as a fraction of files being processed over total files.
        """
        counters = app.storage.general[self.mainuuid][self.name]["counters"]
        if counters["bam_count"] == 0:
            return 0.0
        return counters["bams_in_processing"] / counters["bam_count"]

    @property
    def _not_analysed(self) -> float:
        """
        Calculate the fraction of BAM files not yet analysed.

        Returns:
            float: The fraction of files not yet processed over total files.
        """
        counters = app.storage.general[self.mainuuid][self.name]["counters"]
        if counters["bam_count"] == 0:
            return 0.0
        return (
            counters["bam_count"]
            - counters["bams_in_processing"]
            - counters["bam_processed"]
        ) / counters["bam_count"]

    def progress_bars(self) -> None:
        """
        Show a progress bar for the number of BAM files processed.
        """
        self.progress_trackers = ui.card().classes("w-full")
        with self.progress_trackers:
            with ui.row():
                ui.label("File Tracker").tailwind("drop-shadow", "font-bold")
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid][self.name]["counters"],
                    "bam_count",
                    backward=lambda n: f"Bam files seen: {n}",
                )
                if self.batch:
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid][self.name]["counters"],
                        "bams_in_processing",
                        backward=lambda n: f"Bam files being processed: {n}",
                    )
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid][self.name]["counters"],
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

    def playback(
        self, data: pd.DataFrame, step_size: int = 2, start_time: Optional[float] = None
    ) -> None:
        """
        Simulate the processing of BAM files from a dataframe in playback mode.

        Args:
            data (pd.DataFrame): The dataframe containing BAM file information.
            step_size (int): The step size for the playback.
            start_time (Optional[float]): The start time of the playback.
        """
        self.data = data
        playback_thread = threading.Thread(
            target=self.playback_bams, args=(step_size, start_time)
        )
        playback_thread.daemon = True
        playback_thread.start()

    def playback_bams(
        self, step_size: int = 2, start_time: Optional[float] = None
    ) -> None:
        """
        Simulate the processing of BAM files from a pandas dataframe.

        This function plays back the processing of BAM files from a pandas dataframe.
        To simulate the behavior of a run, it adds files to the queue in the order in which they were produced.
        The run time for each individual processing loop is monitored, and the reads are added to the queue
        accounting for this delay. The delay is calculated by the offset parameter which starts as 0.

        Args:
            step_size (int): The step size for the playback.
            start_time (Optional[float]): The start time of the playback.
        """
        latest_timestamps = self.data
        self.offset = 0
        playback_start_time = time.time()
        for index, row in latest_timestamps.iterrows():
            elapsed_time = (time.time() - playback_start_time) + self.offset
            time_diff = row["file_produced"]
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
            if len(row["full_path"]) > 0:
                # Here we check that the path to the bam file absolutely exists.
                self.add_bam(
                    row["full_path"], playback_start_time + row["file_produced"]
                )

    def create_time_chart(self, title: str) -> ui.echart:
        """
        Create a time chart for visualizing data.

        Args:
            title (str): The title of the chart.

        Returns:
            ui.echart: The Echarts object configured for time series data.
        """
        return (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {
                        "show": True,
                        "feature": {
                            "dataZoom": {"yAxisIndex": "none"},
                            "restore": {},
                            "saveAsImage": {},
                        },
                    },
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    "series": [],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def create_chart(self, title: str) -> ui.echart:
        """
        Create a categorical chart for visualizing data.

        Args:
            title (str): The title of the chart.

        Returns:
            ui.echart: The Echarts object configured for categorical data.
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
                    "series": [],
                }
            )
            .style("color: #6E93D6; font-size: 150%; font-weight: 300; height: 350px")
            .classes("border-double")
        )

    def process_bam(self, bamfile: BinaryIO, timestamp: float) -> None:
        """
        Process a BAM file.

        Args:
            bamfile (BinaryIO): The BAM file to process.
            timestamp (float): A timestamp indicating when the file was generated.
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def setup_ui(self) -> None:
        """
        Set up the user interface for the analysis.

        This function should be overridden by subclasses to create the user interface for the analysis.
        It should create the necessary widgets and add them to the display.
        The function should do no work of its own. That should all be done in the process_bam method.
        The function may call other methods to trigger updates to the display.
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def cleanup(self) -> None:
        """
        Clean up any resources used by the analysis.

        This function is called when the app is closed.
        It should be overridden by subclasses to clean up any resources used by the analysis.
        """
        pass

    def check_resources(self) -> None:
        """
        Check the resources required for the analysis.

        This function is called to check the resources required for the analysis are present.
        It should be overridden by subclasses to check the resources required for the analysis are present.
        """
        pass
