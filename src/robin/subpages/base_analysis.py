"""
Base Class for Analysis of BAM Files Output During a Sequencing Run

This module provides a base class (`BaseAnalysisOrig`) for the analysis of BAM files generated during a sequencing run. The class includes a queue to receive BAM files and a background thread for data processing. It also features a render pipeline for displaying progress and results.

Key Components:

1. **BaseAnalysisOrig Class**:
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
- `nicegui` (ui, app, run)
- `typing` (BinaryIO, Optional, List, Tuple, Dict)
- `pandas`
- `time`
- `asyncio`
- `threading`
- `collections.Counter`
- `logging`

Example usage::

    from base_analysis import BaseAnalysisOrig

    class CustomAnalysis(BaseAnalysisOrig):
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
import asyncio

from robin.core.state import state, ProcessType, ProcessState

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


class BaseVis:
    """
    Base class for visualization components of BAM file analysis.
    Handles all UI and visualization-related functionality.
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
        high_confidence_threshold: float = 0.9,
        medium_confidence_threshold: float = 0.7,
        *args,
        **kwargs,
    ) -> None:
        self.name = analysis_name
        self.mainuuid = uuid
        self.file_mod_times = {}
        self.high_confidence_threshold = high_confidence_threshold
        self.medium_confidence_threshold = medium_confidence_threshold
        self.summary = summary
        self.MENU_BREAKPOINT = 1200
        self.progress = progress
        self.browse = browse
        self.force_sampleid = force_sampleid
        self.output = output
        self.bamqueue = bamqueue
        self.batch = batch
        self.start_time = start_time
        self.threads = threads
        self.sampleID = sampleID

    async def render_ui(self, sample_id=None) -> None:
        """Set up the visualization components of the UI"""
        if sample_id:
            self.sampleID = sample_id
        await self.setup_ui()
        if self.progress and not self.browse:
            self.progress_bars()

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
                logger.info(f"Folder created: {full_path}")
            return full_path
        else:
            return path

    def setup_ui(self) -> None:
        """Set up visualization-specific UI components"""
        raise NotImplementedError("Subclasses must implement this method.")

    def create_time_chart(self, title: str) -> ui.echart:
        """Create a time-based chart for visualizing data."""
        return (
            ui.echart(
                {
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "tooltip": {"order": "valueDesc", "trigger": "axis"},
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
            .classes("border-double text-sky-600 dark:text-white dark:bg-black")
            .style("height: 350px")
        )

    def create_chart(self, title: str) -> ui.echart:
        """Create a categorical chart for visualizing data."""
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
            .classes("border-double text-sky-600 dark:text-white dark:bg-black")
            .style("font-size: 150%; font-weight: 300; height: 350px")
        )

    def create_summary_card(
        self,
        classification_text,
        confidence_value,
        model_name=None,
        features_found=None,
    ):
        """Create a standardized summary card for displaying classification results."""
        with ui.card().classes("w-full shadow-sm rounded-lg").style(
            "background-color: #F5F5F7; padding: 16px; margin-bottom: 16px;"
        ):
            with ui.row().classes("w-full items-center justify-between gap-4"):
                with ui.column().classes("gap-2 flex-grow"):
                    ui.label(classification_text.split(" - ")[0]).classes(
                        "text-lg font-medium text-gray-900"
                    )
                    with ui.row().classes("gap-4 items-center"):
                        if model_name:
                            ui.label(f"Model: {model_name}").classes(
                                "text-sm text-gray-500"
                            )
                        if features_found:
                            ui.label(f"Features: {features_found:,}").classes(
                                "text-sm text-gray-500"
                            )

                confidence_percent = confidence_value * 100
                with ui.column().classes("items-end gap-1 min-w-fit"):
                    ui.label(f"{confidence_percent:.1f}%").classes(
                        "text-xl font-semibold"
                    ).style(f"color: {self._get_confidence_color(confidence_value)}")
                    ui.label(self._get_confidence_text(confidence_value)).classes(
                        "text-sm text-gray-500"
                    )

    def _get_confidence_color(self, confidence):
        """Get color based on confidence level."""
        if confidence >= self.high_confidence_threshold:
            return "#34C759"  # Green for high confidence
        elif confidence >= self.medium_confidence_threshold:
            return "#007AFF"  # Blue for medium confidence
        else:
            return "#FF9500"  # Orange for low confidence

    def _get_confidence_text(self, confidence):
        """Get descriptive text based on confidence level."""
        if confidence >= self.high_confidence_threshold:
            return "High confidence"
        elif confidence >= self.medium_confidence_threshold:
            return "Medium confidence"
        else:
            return "Low confidence"

    def progress_bars(self) -> None:
        """Show progress bars for file processing."""
        try:
            self.progress_trackers = ui.card().classes("w-full p-2")
            logger.debug(f"Progress bars for {self.name} and {self.sampleID}")

            self._initialize_counters(self.sampleID)

            with self.progress_trackers:
                with ui.expansion().classes("w-full") as expansion:
                    with expansion.add_slot("header"):
                        with ui.row().classes(
                            "w-full items-center justify-between gap-1"
                        ):
                            with ui.column().classes("flex-grow gap-0"):
                                ui.label("File Processing Status").classes(
                                    "text-sm font-medium mb-0 pb-0"
                                )
                                with ui.row().classes("w-full items-center gap-1"):
                                    progress_summary = (
                                        ui.linear_progress(
                                            size="4px",
                                            show_value=False,
                                            value=0,
                                            color="primary",
                                        )
                                        .classes("flex-grow my-1")
                                        .props("instant-feedback")
                                    )
                                    ui.label().bind_text_from(
                                        app.storage.general[self.mainuuid][
                                            self.sampleID
                                        ][self.name]["counters"],
                                        "bam_processed",
                                        backward=lambda n: f"{n}/{self._get_counter_value('bam_count', 0)}",
                                    ).classes("text-xs min-w-fit")

                    with ui.column().classes("w-full gap-0 mt-1"):
                        with ui.row().classes("w-full items-center gap-1 py-0"):
                            ui.label("Not processed").classes("text-xs min-w-[80px]")
                            progressbar3 = (
                                ui.linear_progress(
                                    size="4px",
                                    show_value=False,
                                    value=0,
                                    color="red",
                                )
                                .tooltip("Indicates read files not yet processed.")
                                .props("instant-feedback")
                                .classes("flex-grow my-0")
                            )
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid][self.sampleID][
                                    self.name
                                ]["counters"],
                                "bam_count",
                                backward=lambda n: str(
                                    max(
                                        0,
                                        n
                                        # - self._get_counter_value("bam_processed", 0)
                                        - self._get_counter_value(
                                            "bams_in_processing", 0
                                        )
                                        - self._get_counter_value("bam_processed", 0),
                                    )
                                ),
                            ).classes("text-xs min-w-[30px] text-right")

                        with ui.row().classes("w-full items-center gap-1 py-0"):
                            ui.label("Completed").classes("text-xs min-w-[80px]")
                            progressbar = (
                                ui.linear_progress(
                                    size="4px",
                                    show_value=False,
                                    value=0,
                                    color="green",
                                )
                                .tooltip("Indicates read files processed.")
                                .props("instant-feedback")
                                .classes("flex-grow my-0")
                            )
                            ui.label().bind_text_from(
                                app.storage.general[self.mainuuid][self.sampleID][
                                    self.name
                                ]["counters"],
                                "bam_processed",
                                backward=lambda n: str(n or 0),
                            ).classes("text-xs min-w-[30px] text-right")

            if not state.shutdown_event:
                ui.timer(
                    1,
                    callback=lambda: progress_summary.set_value(
                        self._progress_processed
                    ),
                )  # total
                ui.timer(
                    1, callback=lambda: progressbar3.set_value(self._not_analysed)
                )  # not processed
                ui.timer(
                    1, callback=lambda: progressbar.set_value(self._progress_processed)
                )  # processed

        except Exception as e:
            logging.error(f"Error creating progress bars for {self.name}: {str(e)}")
            with ui.card().classes("w-full p-2 text-red-500"):
                ui.label(f"Error displaying progress: {str(e)}")

    def _initialize_counters(self, sample_id: str) -> None:
        """Initialize counters for a specific sample if they don't exist."""
        try:
            if not sample_id:
                logging.warning(f"No sample ID provided for {self.name}")
                return

            if not self.mainuuid:
                logging.warning(f"No UUID provided for {self.name}")
                return

            if self.mainuuid not in app.storage.general:
                app.storage.general[self.mainuuid] = {}
                logging.debug(f"Initialized storage for UUID {self.mainuuid}")

            if sample_id not in app.storage.general[self.mainuuid]:
                app.storage.general[self.mainuuid][sample_id] = {}
                logging.debug(f"Initialized storage for sample {sample_id}")

            if self.name not in app.storage.general[self.mainuuid][sample_id]:
                app.storage.general[self.mainuuid][sample_id][self.name] = {
                    "counters": Counter(
                        bam_count=0, bam_processed=0, bams_in_processing=0
                    )
                }
                logging.debug(
                    f"Initialized counters for {self.name} in sample {sample_id}"
                )
            elif (
                "counters"
                not in app.storage.general[self.mainuuid][sample_id][self.name]
            ):
                app.storage.general[self.mainuuid][sample_id][self.name]["counters"] = (
                    Counter(bam_count=0, bam_processed=0, bams_in_processing=0)
                )
                logging.debug(
                    f"Initialized missing counters for {self.name} in sample {sample_id}"
                )

        except Exception as e:
            logging.error(
                f"Error initializing counters for {self.name} in sample {sample_id}: {str(e)}"
            )
            try:
                app.storage.general[self.mainuuid] = {
                    sample_id: {
                        self.name: {
                            "counters": Counter(
                                bam_count=0, bam_processed=0, bams_in_processing=0
                            )
                        }
                    }
                }
            except Exception as nested_e:
                logging.error(f"Failed to create fallback counters: {str(nested_e)}")

    def _get_counter_value(self, counter_name: str, default: int = 0) -> int:
        """Safely get a counter value from storage."""
        try:
            return app.storage.general[self.mainuuid][self.sampleID][self.name][
                "counters"
            ].get(counter_name, default)
        except (KeyError, AttributeError):
            return default

    @property
    def _progress(self) -> float:
        """Calculate the progress of BAM file processing."""
        try:
            self._initialize_counters(self.sampleID)
            counters = app.storage.general[self.mainuuid][self.sampleID][self.name][
                "counters"
            ]
            if counters.get("bam_processed", 0) == 0:
                return 0.0
            return counters.get("bam_processed", 0) / counters.get("bam_count", 1)
        except (KeyError, TypeError):
            return 0.0

    @property
    def _progress_processed(self) -> float:
        """Calculate the progress of BAM file processing."""
        try:
            self._initialize_counters(self.sampleID)
            counters = app.storage.general[self.mainuuid][self.sampleID][self.name][
                "counters"
            ]
            if counters.get("bam_processed", 0) == 0:
                return 0.0
            if counters.get("bam_count") == 0:
                return 0.0
            return counters.get("bam_processed", 0) / counters.get("bam_count", 1)
        except (KeyError, TypeError):
            return 0.0

    @property
    def _not_analysed(self) -> float:
        """Calculate the fraction of BAM files not yet analysed."""
        try:
            self._initialize_counters(self.sampleID)
            counters = app.storage.general[self.mainuuid][self.sampleID][self.name][
                "counters"
            ]
            if counters.get("bam_count", 0) == 0:
                return 0.0
            not_analysed = max(
                0,
                (
                    counters.get("bam_count", 0)
                    - counters.get("bams_in_processing", 0)
                    - counters.get("bam_processed", 0)
                ),
            )
            if not_analysed == 0 or counters.get("bam_count") == 0:
                return 0.0
            else:
                return not_analysed / counters.get("bam_count", 1)
        except (KeyError, TypeError):
            return 0.0


def sort_bams(data_list):
    sorted_data = sorted(
        data_list,
        key=lambda x: x[1] if x[1] is not None else float("inf"),
    )
    bamfiles = [item[0] for item in sorted_data]
    latest_timestamp = sorted_data[-1][1] if sorted_data else None
    return bamfiles, latest_timestamp


class BaseAnalysis:
    """
    Base class for analysis of BAM files output during a sequencing run.
    Handles core analysis functionality and file processing.
    """

    def __init__(
        self,
        threads: int = 1,
        output: str = None,
        analysis_name: str = None,
        summary: Optional[str] = None,
        bamqueue: Optional[queue.Queue] = None,
        parquetqueue: Optional[queue.Queue] = None,
        progress: bool = False,
        batch: bool = False,
        start_time: Optional[float] = None,
        browse: bool = False,
        uuid: Optional[str] = None,
        force_sampleid: Optional[str] = None,
        sampleID: Optional[str] = None,
        high_confidence_threshold: float = 0.9,
        medium_confidence_threshold: float = 0.7,
        *args,
        **kwargs,
    ) -> None:
        # Initialize core analysis properties
        self.bamqueue = bamqueue if bamqueue else queue.Queue()
        self.parquetqueue = parquetqueue if parquetqueue else queue.Queue()
        self.name = analysis_name
        self.mainuuid = uuid
        self.start_time = start_time
        self.batch = batch
        self.output = output
        self.summary = summary
        self.browse = browse
        self.progress = progress
        self.file_mod_times = {}
        self.running = False
        self.terminate = False
        self.force_sampleid = force_sampleid
        self.threads = max(1, threads // 2)
        self.sampleID = sampleID

        # Initialize visualization component
        # self.vis = BaseVis(
        #    analysis_name=analysis_name,
        #    uuid=uuid,
        #    high_confidence_threshold=high_confidence_threshold,
        #    medium_confidence_threshold=medium_confidence_threshold
        # )
        # self.vis.progress = progress
        # self.vis.browse = browse

    async def process_data(self) -> None:
        """Start processing BAM files either in batch mode or in a continuous timer mode."""
        state.start_process(
            self.name, ProcessType.BATCH if self.batch else ProcessType.PER_FILE
        )
        if self.batch:
            self.bams: Dict[str, List[Tuple[BinaryIO, Optional[float]]]] = {}
            self.batch_timer_run()
        else:
            self.timer_run()

    def timer_run(self) -> None:
        """Set up a timer to periodically run the _worker method for processing BAM files."""
        self.timer = app.timer(0.1, self._worker, once=True)

    async def _worker(self) -> None:
        """Process BAM files from the queue in the background."""
        while not self.terminate:
            # self.timer.active = False
            if state.shutdown_event:
                print(f"shutdown_event is True, stopping _worker in {self.name}")
                await self.stop_analysis()
                return

            if not self.bamqueue.empty() and not self.running:
                self.running = True
                state.set_process_state(self.name, ProcessState.RUNNING)
                bamfile, timestamp, sampleID = self.bamqueue.get()
                logger.debug(
                    f"Processing BAM file: {bamfile} with timestamp: {timestamp} and sampleID: {sampleID}"
                )

                # Validate sampleID - if None, use a fallback
                if sampleID is None:
                    if hasattr(bamfile, "name"):
                        sampleID = f"unknown_sample_{os.path.basename(bamfile.name)}"
                    else:
                        sampleID = f"unknown_sample_{id(bamfile)}"
                    logger.warning(
                        f"BAM file has no sampleID in queue, using fallback: {sampleID}"
                    )

                self.sampleID = sampleID
                # Initialize counters before accessing them
                self._initialize_counters(sampleID)
                app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                    "bam_count"
                ] += 1
                if not timestamp:
                    timestamp = time.time()
                try:
                    logger.debug(
                        f"Calling process_bam with bamfile: {bamfile}, timestamp: {timestamp}"
                    )
                    await self.process_bam(bamfile, timestamp)
                except ValueError as e:
                    if "invalid literal for int() with base 10" in str(e):
                        logger.error(
                            f"Error processing BAM file: {e}. Skipping this file and continuing."
                        )
                        logger.debug(f"BAM file that caused error: {bamfile}")
                        logger.debug(f"Timestamp that caused error: {timestamp}")
                        logger.debug(f"SampleID that caused error: {sampleID}")
                    else:
                        logger.error(f"Error processing BAM files: {e}")
                except Exception as e:
                    logger.error(f"Error processing BAM files: {e}")
                    logger.debug(
                        f"Unexpected error occurred while processing BAM file: {bamfile}"
                    )
                    logger.debug(f"Error details: {str(e)}")
                    logger.debug(f"Error type: {type(e)}")
                finally:
                    if self.sampleID:
                        self._initialize_counters(self.sampleID)
                        self._update_bam_processed_counter(1)
                    else:
                        logger.warning(
                            f"No sample ID available for {self.name} in finally block - mainuuid: {getattr(self, 'mainuuid', 'None')}, sampleID: {getattr(self, 'sampleID', 'None')}"
                        )
                self.running = False
            state.set_process_state(self.name, ProcessState.WAITING_FOR_DATA)
            # if not self.terminate:
            # self.timer.active = True
            await asyncio.sleep(1)

    def add_bam(
        self,
        bamfile: BinaryIO,
        timestamp: Optional[float] = None,
        sampleID: Optional[str] = None,
    ) -> None:
        """Adds a BAM file to the queue for processing."""
        if sampleID is None:
            sampleID = self.sampleID

        # Validate sampleID - if still None, use a fallback
        if sampleID is None:
            if hasattr(bamfile, "name"):
                sampleID = f"unknown_sample_{os.path.basename(bamfile.name)}"
            else:
                sampleID = f"unknown_sample_{id(bamfile)}"
            logger.warning(f"BAM file has no sampleID, using fallback: {sampleID}")

        # Initialize counters and increment bam_count for single processing mode
        if not self.batch:
            self._initialize_counters(sampleID)
            try:
                app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                    "bam_count"
                ] += 1
            except Exception as e:
                logger.warning(f"Could not increment bam_count counter: {e}")
                # Try to re-initialize and increment again
                try:
                    self._initialize_counters(sampleID)
                    app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                        "bam_count"
                    ] += 1
                except Exception as retry_e:
                    logger.error(
                        f"Failed to increment bam_count counter after retry: {retry_e}"
                    )

        self.bamqueue.put((bamfile, timestamp, sampleID))

    def batch_timer_run(self) -> None:
        """Set up a timer to periodically run the _batch_worker method for batch processing."""
        self.timer = app.timer(5, self._batch_worker)

    async def _batch_worker(self) -> None:
        """Process BAM files from the queue in batches in the background."""
        if state.shutdown_event:
            print(f"shutdown_event is True, stopping _batch_worker in {self.name}")
            await self.stop_analysis()
            return
        self.timer.active = False
        latest_timestamp = time.time()
        await self.process_bam("SasuagesNone", latest_timestamp)

        """
        if self.bamqueue.qsize() > 0:
            count = 0
            while self.bamqueue.qsize() > 0:
                state.set_process_state(self.name, ProcessState.RUNNING)

                bamfile, timestamp, sampleID = self.bamqueue.get()
                print(f"Batch worker: {sampleID} {bamfile} {timestamp}")

                if sampleID not in self.bams:
                    self.bams[sampleID] = []
                    self._initialize_counters(sampleID)

                self.bams[sampleID].append((bamfile, timestamp))
                count += 1

                app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                    "bam_count"
                ] += 1

                #if count >= 50:
                #    break
            
            for sample_id, data_list in self.bams.items():
                if not self.running and len(data_list) > 0:
                    state.set_process_state(self.name, ProcessState.RUNNING)
                    self.running = True
                    self.sampleID = sample_id
                    batch_size = len(data_list)
                    print(f"Batch worker: {sample_id} {batch_size}")
                    try:
                        bamfiles, latest_timestamp = await run.cpu_bound(
                            sort_bams,
                            data_list,
                        )

                        await self.process_bam(bamfiles, latest_timestamp)

                    except Exception as e:
                        logger.error(f"Error processing BAM files: {e}")
                        print(f"Error processing BAM files in {self.name}: {e}")
                    finally:
                        # Always update bam_processed counter for each file in the batch
                        # regardless of success or failure
                        #self._update_bam_processed_counter(batch_size)
                        # Clear the processed batch
                        self.bams[sampleID] = []
                        self.running = False
        """
        state.set_process_state(self.name, ProcessState.WAITING_FOR_DATA)
        self.timer.active = True

    def check_file_time(self, file_path: str) -> bool:
        """Check if the file exists and whether it has been modified since last seen."""
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

    def check_and_create_folder(self, path, folder_name=None):
        """Check if the path exists and create the folder if it doesn't."""
        if not os.path.exists(path):
            raise FileNotFoundError(f"The specified path does not exist: {path}")

        if self.force_sampleid:
            folder_name = self.force_sampleid

        if folder_name:
            full_path = os.path.join(path, folder_name)
            if not os.path.exists(full_path):
                os.makedirs(full_path)
                logger.info(f"Folder created: {full_path}")
            return full_path
        else:
            return path

    def _update_bam_processed_counter(self, increment: int = 1) -> None:
        """Safely update the bam_processed counter."""
        if (
            hasattr(self, "mainuuid")
            and hasattr(self, "sampleID")
            and hasattr(self, "name")
        ):
            try:
                # Ensure sampleID is valid
                if not self.sampleID:
                    logger.warning(
                        f"No sample ID available for {self.name} - mainuuid: {getattr(self, 'mainuuid', 'None')}, sampleID: {getattr(self, 'sampleID', 'None')}"
                    )
                    return

                # Ensure counters are initialized before updating
                self._initialize_counters(self.sampleID)
                app.storage.general[self.mainuuid][self.sampleID][self.name][
                    "counters"
                ]["bam_processed"] += increment
            except Exception as e:
                logger.warning(f"Could not update bam_processed counter: {e}")
                # Try to re-initialize and update again
                try:
                    if self.sampleID:
                        self._initialize_counters(self.sampleID)
                        app.storage.general[self.mainuuid][self.sampleID][self.name][
                            "counters"
                        ]["bam_processed"] += increment
                except Exception as retry_e:
                    logger.error(
                        f"Failed to update bam_processed counter after retry: {retry_e}"
                    )

    def _initialize_counters(self, sample_id: str) -> None:
        """
        Initialize counters for a specific sample if they don't exist.

        Args:
            sample_id (str): The ID of the sample to initialize counters for.
        """
        try:
            if not sample_id:
                logging.warning(f"No sample ID provided for {self.name}")
                return

            if not self.mainuuid:
                logging.warning(f"No UUID provided for {self.name}")
                return

            # Ensure sample_id is a string
            sample_id = str(sample_id)

            if self.mainuuid not in app.storage.general:
                app.storage.general[self.mainuuid] = {}
                logging.debug(f"Initialized storage for UUID {self.mainuuid}")

            if sample_id not in app.storage.general[self.mainuuid]:
                app.storage.general[self.mainuuid][sample_id] = {}
                logging.debug(f"Initialized storage for sample {sample_id}")

            if self.name not in app.storage.general[self.mainuuid][sample_id]:
                app.storage.general[self.mainuuid][sample_id][self.name] = {
                    "counters": Counter(
                        bam_count=0, bam_processed=0, bams_in_processing=0
                    )
                }
                logging.debug(
                    f"Initialized counters for {self.name} in sample {sample_id}"
                )
            elif (
                "counters"
                not in app.storage.general[self.mainuuid][sample_id][self.name]
            ):
                app.storage.general[self.mainuuid][sample_id][self.name]["counters"] = (
                    Counter(bam_count=0, bam_processed=0, bams_in_processing=0)
                )
                logging.debug(
                    f"Initialized missing counters for {self.name} in sample {sample_id}"
                )

        except Exception as e:
            logging.error(
                f"Error initializing counters for {self.name} in sample {sample_id}: {str(e)}"
            )
            # Create empty counter to prevent further errors
            try:
                app.storage.general[self.mainuuid] = {
                    sample_id: {
                        self.name: {
                            "counters": Counter(
                                bam_count=0, bam_processed=0, bams_in_processing=0
                            )
                        }
                    }
                }
            except Exception as nested_e:
                logging.error(f"Failed to create fallback counters: {str(nested_e)}")

    async def stop_analysis(self):
        """Stop the analysis and clean up resources."""
        state.set_process_state(self.name, ProcessState.STOPPING)
        self.running = False
        self.batch_running = False
        print(
            f"Analysis stopped for {self.name} and {self.sampleID}. Active processes: {list(state.process_states.keys())}"
        )
        if self.name in state.process_states:
            state.stop_process(self.name)

    def process_bam(self, bamfile: BinaryIO, timestamp: float) -> None:
        """Process a BAM file."""
        raise NotImplementedError("Subclasses must implement this method.")

    def setup_ui(self) -> None:
        """Set up the user interface for the analysis."""
        raise NotImplementedError("Subclasses must implement this method.")

    def cleanup(self) -> None:
        """Clean up any resources used by the analysis."""
        pass

    def check_resources(self) -> None:
        """Check the resources required for the analysis."""
        pass

    def playback(
        self, data: pd.DataFrame, step_size: int = 2, start_time: Optional[float] = None
    ) -> None:
        """Simulate the processing of BAM files from a dataframe in playback mode."""
        self.data = data
        playback_thread = threading.Thread(
            target=self.playback_bams, args=(step_size, start_time)
        )
        playback_thread.daemon = True
        playback_thread.start()

    def playback_bams(
        self, step_size: int = 2, start_time: Optional[float] = None
    ) -> None:
        """Simulate the processing of BAM files from a pandas dataframe."""
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
                    elapsed_time = (time.time() - playback_start_time) + self.offset
                else:
                    self.offset += step_size
                    elapsed_time += self.offset
            if len(row["full_path"]) > 0:
                self.add_bam(
                    row["full_path"], playback_start_time + row["file_produced"]
                )
