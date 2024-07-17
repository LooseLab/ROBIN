"""
Base Class for Analysis of BAM Files Output During a Sequencing Run

This module provides a base class (`BaseAnalysis`) for the analysis of BAM files generated during a sequencing run. The class includes a queue to receive BAM files and a background thread for data processing. It also features a render pipeline for displaying progress and results.

Key Components:

1. **BaseAnalysis Class**:
   - Initializes analysis parameters and sets up queues for BAM file processing.
   - Provides methods to process BAM files either in batch mode or continuous mode.
   - Includes user interface elements for displaying progress and results.

2. **Queue Handling**:
   - `bamqueue`: Queue for holding BAM files to be processed.
   - `_worker`: Background worker function for processing BAM files.
   - `_batch_worker`: Batch processing worker function for handling multiple BAM files.

3. **Progress Tracking**:
   - Properties (`_progress`, `_progress2`, `_not_analysed`) for tracking the progress of BAM file processing.
   - `progress_bars`: Method to render progress bars in the user interface.

4. **Playback Mode**:
   - `playback`: Simulates the processing of BAM files from a pandas DataFrame.
   - `playback_bams`: Handles the playback of BAM files to simulate real-time processing.

5. **User Interface**:
   - `render_ui`: Sets up the user interface and displays progress bars.
   - `create_time_chart`: Creates a time-based chart for visualization.
   - `create_chart`: Creates a categorical chart for visualization.

6. **Abstract Methods**:
   - `process_bam`: Must be implemented by subclasses to define how BAM files are processed.
   - `setup_ui`: Must be implemented by subclasses to set up the user interface for the analysis.
   - `cleanup`: Optional method for cleaning up resources when the app is closed.
   - `check_resources`: Optional method for checking required resources for the analysis.

Dependencies:
- `queue`
- `nicegui` (ui, app)
- `typing` (BinaryIO, Optional, List, Tuple)
- `pandas`
- `time`
- `asyncio`
- `threading`
- `collections.Counter`
- `logging`

Example usage::

    from base_analysis import BaseAnalysis

    class CustomAnalysis(BaseAnalysis):
        def process_bam(self, bamfile, timestamp):
            # Implement custom BAM file processing logic
            pass

        def setup_ui(self):
            # Implement custom UI setup logic
            pass

    analysis = CustomAnalysis(threads=4, outputfolder='/path/to/output', analysis_name='ExampleAnalysis')
    analysis.process_data()

"""

import queue
from nicegui import ui, app
from typing import BinaryIO, Optional, List, Tuple, Dict
import pandas as pd
import time
import threading
from collections import Counter
import os
import logging
from datetime import timedelta

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
        threads: int = 1,
        output: str = None,
        analysis_name: str = None,
        summary: Optional[str] = None,
        bamqueue: Optional[queue.Queue] = None,
        progress: bool = False,
        batch: bool = False,
        start_time: Optional[float] = None,
        browse: bool = False,
        uuid: Optional[str] = None,
        force_sampleid: Optional[str] = None,
        sampleID: Optional[str] = None,
        *args,
        **kwargs,
    ) -> None:
        self.mainuuid = uuid
        self.bamqueue = bamqueue if bamqueue else queue.Queue()
        self.name = analysis_name
        self.start_time = start_time
        self.batch = batch
        self.output = output
        self.summary = summary
        self.browse = browse
        self.progress = progress
        self.file_mod_times = {}
        self.MENU_BREAKPOINT = 520
        self.running = False
        self.force_sampleid = force_sampleid
        self.threads = max(1, threads // 2)
        if sampleID:
            self.sampleID = sampleID
            # print(f"SampleID: {self.sampleID}")
        else:
            self.sampleID = None
        self.module_start_time = time.time()
        self.track_elapsed_time = 0
        self.five_minutes = 0
            # print("No SampleID provided")

    def check_file_time(self, file_path: str) -> bool:
        """
        Check if the file exists and whether it has been modified since last seen.

        Args:
            file_path (str): Path to the file.

        Returns:
            bool: True if the file exists and has been modified, False otherwise.
        """
        if not os.path.exists(file_path):
            return False

        current_mod_time = os.path.getmtime(file_path)

        if file_path not in self.file_mod_times:
            self.file_mod_times[file_path] = current_mod_time
            return True

        if self.file_mod_times[file_path] == current_mod_time:
            return False

        self.file_mod_times[file_path] = current_mod_time
        return True

    def parse_timedelta(self,time_str):
        # This handles strings like "1 day, 2:03:04" or "2:03:04"
        if ',' in time_str:
            days, time = time_str.split(',')
            days = int(days.split()[0])
        else:
            days = 0
            time = time_str
        hours, minutes, seconds = map(int, time.strip().split(':'))
        return timedelta(days=days, hours=hours, minutes=minutes, seconds=seconds)

    def check_and_create_folder(self, path, folder_name=None):
        # Check if the path exists
        if not os.path.exists(path):
            raise FileNotFoundError(f"The specified path does not exist: {path}")

        if self.force_sampleid:
            folder_name = self.force_sampleid

        # If folder_name is provided
        if folder_name:
            full_path = os.path.join(path, folder_name)
            # Create the folder if it doesn't exist
            if not os.path.exists(full_path):
                os.makedirs(full_path)
            return full_path
        else:
            return path

    async def render_ui(self, sample_id=None) -> None:
        """
        Set up the user interface for the analysis and display progress bars if enabled.
        """
        if sample_id:
            self.sampleID = sample_id
        # print(f"Rendering UI for {self.name} and {self.sampleID}")
        self.setup_ui()
        if self.progress and not self.browse:
            self.progress_bars()

    def process_data(self) -> None:
        """
        Start processing BAM files either in batch mode or in a continuous timer mode.
        """
        if self.batch:
            # self.bams: List[Tuple[BinaryIO, Optional[float]]] = []
            self.bams: Dict[str, List[Tuple[BinaryIO, Optional[float]]]] = {}
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
            bamfile, timestamp, sampleID = self.bamqueue.get()
            self.sampleID = sampleID
            if self.sampleID not in app.storage.general[self.mainuuid]:
                app.storage.general[self.mainuuid][self.sampleID] = {}

            if self.name not in app.storage.general[self.mainuuid].get(
                self.sampleID, {}
            ):
                app.storage.general[self.mainuuid].setdefault(self.sampleID, {})[
                    self.name
                ] = {
                    "counters": Counter(
                        bam_count=0, bam_processed=0, bams_in_processing=0
                    )
                }
            app.storage.general[self.mainuuid][self.sampleID][self.name]["counters"][
                "bam_count"
            ] += 1

            if not timestamp:
                timestamp = time.time()
            try:
                await self.process_bam(bamfile, timestamp)
            except Exception as e:
                print(f"Error processing BAM files: {e}")
            app.storage.general[self.mainuuid][self.sampleID][self.name]["counters"][
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

        # Different bam files may come from different runs (sampleIDs). Therefore we must process reads according to the run they have come from.
        while self.bamqueue.qsize() > 0:
            bamfile, timestamp, sampleID = self.bamqueue.get()
            if sampleID not in self.bams:
                self.bams[sampleID] = []
            self.bams[sampleID].append((bamfile, timestamp))
            count += 1
            if sampleID not in app.storage.general[self.mainuuid]:
                app.storage.general[self.mainuuid][sampleID] = {}

            if self.name not in app.storage.general[self.mainuuid].get(sampleID, {}):
                app.storage.general[self.mainuuid].setdefault(sampleID, {})[
                    self.name
                ] = {
                    "counters": Counter(
                        bam_count=0, bam_processed=0, bams_in_processing=0
                    )
                }
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bam_count"
            ] += 1

            if count >= 100:
                break

        for sample_id, data_list in self.bams.items():
            if not self.running and len(data_list) > 0:
                self.running = True
                self.sampleID = sample_id
                if self.sampleID not in app.storage.general[self.mainuuid]:
                    app.storage.general[self.mainuuid][self.sampleID] = {}

                if self.name not in app.storage.general[self.mainuuid].get(
                    self.sampleID, {}
                ):
                    app.storage.general[self.mainuuid].setdefault(self.sampleID, {})[
                        self.name
                    ] = {
                        "counters": Counter(
                            bam_count=0, bam_processed=0, bams_in_processing=0
                        )
                    }
                #app.storage.general[self.mainuuid][self.sampleID][self.name][
                #    "counters"
                #]["bams_in_processing"] += len(data_list)

                try:
                    await self.process_bam(data_list)
                except Exception as e:
                    print(f"Error processing BAM files: {e}")

                #app.storage.general[self.mainuuid][self.sampleID][self.name][
                #        "counters"
                #    ]["bams_in_processing"] -= len(data_list)
                #app.storage.general[self.mainuuid][self.sampleID][self.name][
                #        "counters"
                #    ]["bam_processed"] += len(data_list)
                    # self.running = False

        self.timer.active = True

    @property
    def _progress(self) -> float:
        """
        Calculate the progress of BAM file processing.

        Returns:
            float: The progress as a fraction of processed files over total files.
        """
        counters = app.storage.general[self.mainuuid][self.sampleID][self.name][
            "counters"
        ]
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
        counters = app.storage.general[self.mainuuid][self.sampleID][self.name][
            "counters"
        ]
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
        counters = app.storage.general[self.mainuuid][self.sampleID][self.name][
            "counters"
        ]
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
        # print(f"Progress bars for {self.name} and {self.sampleID}")
        if self.sampleID not in app.storage.general[self.mainuuid]:
            app.storage.general[self.mainuuid][self.sampleID] = {}

        if self.name not in app.storage.general[self.mainuuid].get(self.sampleID, {}):
            app.storage.general[self.mainuuid].setdefault(self.sampleID, {})[
                self.name
            ] = {
                "counters": Counter(bam_count=0, bam_processed=0, bams_in_processing=0)
            }
        with self.progress_trackers:
            with ui.row():
                ui.label("File Tracker").tailwind("drop-shadow", "font-bold")
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid][self.sampleID][self.name][
                        "counters"
                    ],
                    "bam_count",
                    backward=lambda n: f"Bam files seen: {n}",
                )
                if self.batch:
                    ui.label().bind_text_from(
                        app.storage.general[self.mainuuid][self.sampleID][self.name][
                            "counters"
                        ],
                        "bams_in_processing",
                        backward=lambda n: f"Bam files being processed: {n}",
                    )
                ui.label().bind_text_from(
                    app.storage.general[self.mainuuid][self.sampleID][self.name][
                        "counters"
                    ],
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
                    # time.sleep(1)
                    elapsed_time = (time.time() - playback_start_time) + self.offset
                else:
                    # time.sleep(1)
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
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    'tooltip': {
                        'order': 'valueDesc',
                        'trigger': 'axis'
                    },
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
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
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
