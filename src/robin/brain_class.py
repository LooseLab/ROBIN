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
   - `MNPFlex_Object`

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
- `robin.subpages` (MGMT_Object, Sturgeon_object, NanoDX_object, RandomForest_object, CNVAnalysis, TargetCoverage, FusionObject, MNPFlex_Object)
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
from queue import Queue
import pandas as pd
import asyncio
from collections import Counter
from datetime import datetime
from dateutil import parser
import os
import psutil
import tempfile
from nicegui import ui, app, run, background_tasks
import concurrent.futures
import time
import weakref
import atexit

from robin import resources


from robin.utilities.bam_handler import BamEventHandler
from robin.subpages.MGMT_object import MGMTVis, MGMT_Object
from robin.subpages.Sturgeon_object import (
    Sturgeon_object,
    SturgeonVis,
)
from robin.subpages.NanoDX_object import NanoDX_object, NanoDXVis
from robin.subpages.RandomForest_object import (
    RandomForest_object,
    RandomForestVis,
)
from robin.subpages.CNVObjectClass import CNVAnalysis
from robin.subpages.CNV_object import CNVVis
from robin.subpages.TargetCoverage_object import TargetCoverage, TargetCoverageVis
from robin.subpages.Fusion_object import FusionVis
from robin.subpages.FusionObjectClass import FusionObject
from robin.subpages.MNPFlex_object import MNPFlex_Object
from robin.utilities.local_file_picker import LocalFilePicker
from robin.utilities.bed_file import MasterBedTree
from robin.reporting.report import create_pdf
from robin.reporting.sections.disclaimer_text import EXTENDED_DISCLAIMER_TEXT


from watchdog.observers import Observer


from .core.state import state, ProcessType, ProcessState

import json

from robin.utils import (
    merge_modkit_files,
    run_matkit,
    run_samtools_sort,
    check_bam,
    sort_bams,
)

# Configure logging
logger = logging.getLogger(__name__)


class FileHandleManager:
    """
    Manages file handles to prevent leaks and ensure proper cleanup.
    """
    
    def __init__(self):
        self._file_handles = weakref.WeakSet()
        self._temp_files = weakref.WeakSet()
        self._temp_directories = weakref.WeakSet()
        atexit.register(self.cleanup_all)
    
    def register_file_handle(self, file_handle):
        """Register a file handle for tracking and cleanup."""
        self._file_handles.add(file_handle)
        return file_handle
    
    def register_temp_file(self, temp_file):
        """Register a temporary file for tracking and cleanup."""
        self._temp_files.add(temp_file)
        return temp_file
    
    def register_temp_directory(self, temp_directory):
        """Register a temporary directory for tracking and cleanup."""
        self._temp_directories.add(temp_directory)
        return temp_directory
    
    def cleanup_file_handles(self):
        """Clean up all registered file handles."""
        for file_handle in list(self._file_handles):
            try:
                if not file_handle.closed:
                    file_handle.close()
                    logging.debug(f"Closed file handle: {file_handle}")
            except Exception as e:
                logging.error(f"Error closing file handle: {str(e)}", exc_info=True)
        self._file_handles.clear()
    
    def cleanup_temp_files(self):
        """Clean up all registered temporary files."""
        for temp_file in list(self._temp_files):
            try:
                if hasattr(temp_file, 'cleanup'):
                    temp_file.cleanup()
                elif hasattr(temp_file, 'close'):
                    temp_file.close()
                logging.debug(f"Cleaned up temporary file: {temp_file}")
            except Exception as e:
                logging.error(f"Error cleaning up temporary file: {str(e)}", exc_info=True)
        self._temp_files.clear()
    
    def cleanup_temp_directories(self):
        """Clean up all registered temporary directories."""
        for temp_dir in list(self._temp_directories):
            try:
                if hasattr(temp_dir, 'cleanup'):
                    temp_dir.cleanup()
                logging.debug(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                logging.error(f"Error cleaning up temporary directory: {str(e)}", exc_info=True)
        self._temp_directories.clear()
    
    def cleanup_all(self):
        """Clean up all registered resources."""
        logging.info("Cleaning up all file handles and temporary resources")
        self.cleanup_file_handles()
        self.cleanup_temp_files()
        self.cleanup_temp_directories()
    
    def get_stats(self):
        """Get statistics about managed resources."""
        return {
            'file_handles': len(self._file_handles),
            'temp_files': len(self._temp_files),
            'temp_directories': len(self._temp_directories)
        }


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
        mnpflex_config=None,
        enable_snp_calling=False,
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
        :param mnpflex_config: Configuration for MNP-FLEX integration.
        :param enable_snp_calling: Whether to enable SNP calling functionality.
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
        self.mnpflex_config = mnpflex_config  # Store mnpflex_config
        self.enable_snp_calling = enable_snp_calling  # Store enable_snp_calling
        if self.browse:
            self.runsfolder = self.output
        self.sampleID = None
        self.readfish_toml = readfish_toml
        self.watchdogbamqueue = Queue()
        self.dataDir = {}
        self.bedDir = {}
        self.first_run = {}
        self.terminate = False
        self.finished = False
        self.file_handle_manager = FileHandleManager()
        self.cpgs_master_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "sturg_nanodx_cpgs_0125.bed.gz",
        )

        if self.reference:
            # ToDo: We need to pass through an instance of the MasterBedTree class here.
            self.master_bed_tree = MasterBedTree(
                default_preserve_original_tree=True,
                default_reference_file=f"{self.reference}.fai",
                default_readfish_toml=self.readfish_toml,
            )

            # self.NewBed = BedTree(
            #    preserve_original_tree=True,
            #    reference_file=f"{self.reference}.fai",
            #     readfish_toml=self.readfish_toml,
            # )
            # if self.bed_file:
            #    self.NewBed.load_from_file(self.bed_file)
        else:
            self.master_bed_tree = MasterBedTree(
                default_preserve_original_tree=True,
            )

        logging.info(f"BrainMeth initialized with UUID: {self.mainuuid}")

    def configure_storage(self, sample_id):
        """
        Configure storage for the application.

        :param sample_id: Sample ID to configure storage for.
        """
        logging.info(f"Configuring storage for sample ID: {sample_id}")
        app.storage.general[self.mainuuid]["samples"][sample_id]["file_counters"] = {
            "bam_passed": 0,
            "bam_failed": 0,
            "mapped_count": 0,
            "pass_mapped_count": 0,
            "fail_mapped_count": 0,
            "unmapped_count": 0,
            "pass_unmapped_count": 0,
            "fail_unmapped_count": 0,
            "pass_bases_count": 0,
            "fail_bases_count": 0,
            "bases_count": 0,
            # New counters for read numbers
            "mapped_reads_num": 0,
            "unmapped_reads_num": 0,
            "pass_mapped_reads_num": 0,
            "fail_mapped_reads_num": 0,
            "pass_unmapped_reads_num": 0,
            "fail_unmapped_reads_num": 0,
            # New counters for bases in each category
            "mapped_bases": 0,
            "unmapped_bases": 0,
            "pass_mapped_bases": 0,
            "fail_mapped_bases": 0,
            "pass_unmapped_bases": 0,
            "fail_unmapped_bases": 0,
        }
        app.storage.general[self.mainuuid]["samples"][sample_id]["devices"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["basecall_models"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["run_time"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["flowcell_ids"] = []
        app.storage.general[self.mainuuid]["samples"][sample_id]["sample_ids"] = []

    def update_counter(self, sample_id, counter_name, value):
        """
        Safely update a counter value.

        :param sample_id: Sample ID to update counter for
        :param counter_name: Name of the counter to update
        :param value: Value to add to the counter
        """
        try:
            # Ensure sample exists in storage
            if sample_id not in app.storage.general[self.mainuuid]["samples"]:
                app.storage.general[self.mainuuid]["samples"][sample_id] = {}
                self.configure_storage(sample_id)
                app.storage.general[self.mainuuid]["sample_list"].append(sample_id)
                logging.info(f"Initialized storage for new sample: {sample_id}")

            # Ensure file_counters exists
            if (
                "file_counters"
                not in app.storage.general[self.mainuuid]["samples"][sample_id]
            ):
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "file_counters"
                ] = {
                    "bam_passed": 0,
                    "bam_failed": 0,
                    "mapped_count": 0,
                    "pass_mapped_count": 0,
                    "fail_mapped_count": 0,
                    "unmapped_count": 0,
                    "pass_unmapped_count": 0,
                    "fail_unmapped_count": 0,
                    "pass_bases_count": 0,
                    "fail_bases_count": 0,
                    "bases_count": 0,
                    "mapped_reads_num": 0,
                    "unmapped_reads_num": 0,
                    "pass_mapped_reads_num": 0,
                    "fail_mapped_reads_num": 0,
                    "pass_unmapped_reads_num": 0,
                    "fail_unmapped_reads_num": 0,
                    "mapped_bases": 0,
                    "unmapped_bases": 0,
                    "pass_mapped_bases": 0,
                    "fail_mapped_bases": 0,
                    "pass_unmapped_bases": 0,
                    "fail_unmapped_bases": 0,
                }
                logging.info(f"Initialized file_counters for sample: {sample_id}")

            # Ensure counter exists
            if (
                counter_name
                not in app.storage.general[self.mainuuid]["samples"][sample_id][
                    "file_counters"
                ]
            ):
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "file_counters"
                ][counter_name] = 0
                logging.info(
                    f"Initialized counter {counter_name} for sample: {sample_id}"
                )

            current = app.storage.general[self.mainuuid]["samples"][sample_id][
                "file_counters"
            ][counter_name]
            new_value = max(0, current + value)  # Ensure we don't go below 0
            app.storage.general[self.mainuuid]["samples"][sample_id]["file_counters"][
                counter_name
            ] = new_value
            logging.debug(
                f"Counter {counter_name} updated for sample {sample_id}: {current} -> {new_value}"
            )

        except Exception as e:
            logging.error(
                f"Error updating counter {counter_name} for sample {sample_id}: {str(e)}",
                exc_info=True,
            )
            # Re-raise the exception to ensure we don't silently fail
            raise

    def shutdown_background(self):
        """
        Shutdown the application gracefully.
        """
        print("Shutting down ROBIN... from brainclas?")
        print(
            "Here we need to do some very graceful shutdown to make sure we don't leave any threads running and we don't leave any files open."
        )
        self.terminate = True

        # First stop all timers
        print("Stop observer")
        if hasattr(self, "observer"):
            state.set_process_state("File Observer", ProcessState.STOPPING)
            state.stop_process("File Observer")
            self.observer.stop()
            self.observer.join()
        print("Observer stopped")

        print("Stop check_existing_bams_timer")
        if hasattr(self, "check_existing_bams_timer"):
            self.check_existing_bams_timer.cancel()
        print("check_existing_bams_timer stopped")

        print("Stop background_process_bams_timer")
        if hasattr(self, "background_process_bams_timer"):
            self.background_process_bams_timer.cancel()
        print("background_process_bams_timer stopped")

        # Stop additional timers that may have been created
        print("Stop app_state_timer")
        if hasattr(self, "app_state_timer"):
            self.app_state_timer.cancel()
        print("app_state_timer stopped")

        print("Stop process_bams_tracker")
        if hasattr(self, "process_bams_tracker"):
            self.process_bams_tracker.cancel()
        print("process_bams_tracker stopped")

        print("Stop process_bigbadmerge_tracker")
        if hasattr(self, "process_bigbadmerge_tracker"):
            self.process_bigbadmerge_tracker.cancel()
        print("process_bigbadmerge_tracker stopped")

        print("Stop check_and_log_memory")
        if hasattr(self, "check_and_log_memory"):
            self.check_and_log_memory.cancel()
        print("check_and_log_memory stopped")

        # Clean up all temporary directories and file handles
        print("Cleaning up temporary directories and file handles")
        self.cleanup_all_temp_directories()
        self.file_handle_manager.cleanup_all()
        print("Temporary directories and file handles cleaned up")

        # Wait for any pending tasks to complete
        while not self.finished:
            print("Waiting for finished")
            print(f"Value of terminate: {self.terminate}")
            print(f"Value of finished: {self.finished}")
            time.sleep(0.1)

        if self.mainuuid in app.storage.general:
            app.storage.general.pop(self.mainuuid)
        print("Shutdown complete")

        state.shutdown_event = False

    async def start_background(self):
        """
        Start background tasks for BAM processing and analysis.
        """
        logging.info(f"Starting background tasks for UUID: {self.mainuuid}")

        async def start_background_tasks():
            # Ensure storage directory exists
            storage_dir = os.path.join(os.path.expanduser("~"), ".nicegui")
            os.makedirs(storage_dir, exist_ok=True)

            # Initialize storage if it doesn't exist
            if not hasattr(app.storage, "general"):
                app.storage.general = {}

            if self.mainuuid not in app.storage.general:
                app.storage.general[self.mainuuid] = {}

            app.storage.general[self.mainuuid]["bam_count"] = Counter(counter=0)
            app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
            app.storage.general[self.mainuuid]["bam_count"]["total_files"] = 0
            app.storage.general[self.mainuuid]["samples"] = {}
            app.storage.general[self.mainuuid]["sample_list"] = []

            # Initialize queues
            logging.info("Initializing processing queues")
            self.bam_tracking = Queue()
            self.bamforcns = Queue()

            self.bamforcnv = Queue()
            self.bamfortargetcoverage = Queue()
            self.bamformgmt = Queue()
            self.bamforfusions = Queue()
            self.bamforbigbadmerge = Queue()
            #self.mergecounter = 0

            if self.watchfolder:
                logging.info(f"Adding watchfolder: {self.watchfolder}")
                await self.add_watchfolder(self.watchfolder)

            common_args = {
                "threads": self.threads,
                "output": self.output,
                "progress": True,
                "browse": self.browse,
                "uuid": self.mainuuid,
                "force_sampleid": self.force_sampleid,
            }

            # Initialize analysis objects with batch processing
            if "sturgeon" not in self.exclude:
                logging.info("Initializing Sturgeon analysis with batch processing")
                self.parquetqueuesturgeon = Queue()
                self.bamforsturgeon = Queue()
                self.Sturgeon = Sturgeon_object(
                    analysis_name="STURGEON",
                    batch=True,
                    bamqueue=self.bamforsturgeon,
                    parquetqueue=self.parquetqueuesturgeon,
                    **common_args,
                )

                await self.Sturgeon.process_data()

            if "nanodx" not in self.exclude:
                logging.info("Initializing NanoDX analysis with batch processing")
                self.bamfornanodx = Queue()
                self.parquetqueuenanodx = Queue()
                self.NanoDX = NanoDX_object(
                    analysis_name="NANODX",
                    batch=True,
                    bamqueue=self.bamfornanodx,
                    parquetqueue=self.parquetqueuenanodx,
                    **common_args,
                )
                await self.NanoDX.process_data()

            if "pannanodx" not in self.exclude:
                logging.info("Initializing PanNanoDX analysis with batch processing")
                self.bamforpannanodx = Queue()
                self.parquetqueuepannanodx = Queue()
                self.panNanoDX = NanoDX_object(
                    analysis_name="PANNANODX",
                    batch=True,
                    bamqueue=self.bamforpannanodx,
                    model="pancan_devel_v5i_NN.pkl",
                    parquetqueue=self.parquetqueuepannanodx,
                    **common_args,
                )
                await self.panNanoDX.process_data()

            if "forest" not in self.exclude:
                logging.info("Initializing RandomForest analysis with batch processing")
                self.bamforcns = Queue()
                self.parquetqueuecns = Queue()
                self.RandomForest = RandomForest_object(
                    analysis_name="FOREST",
                    batch=True,
                    showerrors=self.showerrors,
                    bamqueue=self.bamforcns,
                    parquetqueue=self.parquetqueuecns,
                    **common_args,
                )
                await self.RandomForest.process_data()

            if "cnv" not in self.exclude:
                logging.info("Initializing CNV analysis")
                self.CNV = CNVAnalysis(
                    analysis_name="CNV",
                    bamqueue=self.bamforcnv,
                    target_panel=self.target_panel,
                    reference_file=self.reference,
                    bed_file=self.bed_file,
                    readfish_toml=self.readfish_toml,  # ToDo: This assumes a single sample per CNV analysis.
                    # NewBed=self.NewBed, #ToDo: This assumes a single sample per CNV analysis.
                    master_bed_tree=self.master_bed_tree,
                    **common_args,
                )
                await self.CNV.process_data()

            if "coverage" not in self.exclude:
                logging.info("Initializing Coverage analysis")
                self.Target_Coverage = TargetCoverage(
                    analysis_name="COVERAGE",
                    showerrors=self.showerrors,
                    bamqueue=self.bamfortargetcoverage,
                    target_panel=self.target_panel,
                    reference=self.reference,
                    enable_snp_calling=self.enable_snp_calling
                    and self.reference is not None,
                    **common_args,
                )
                await self.Target_Coverage.process_data()

            if "mgmt" not in self.exclude:
                logging.info("Initializing MGMT analysis")
                self.MGMT_panel = MGMT_Object(
                    analysis_name="MGMT",
                    bamqueue=self.bamformgmt,
                    **common_args,
                )
                await self.MGMT_panel.process_data()

            if "fusion" not in self.exclude:
                logging.info("Initializing Fusion analysis")
                self.Fusion_panel = FusionObject(
                    analysis_name="FUSION",
                    bamqueue=self.bamforfusions,
                    target_panel=self.target_panel,
                    reference_file=self.reference,
                    bed_file=self.bed_file,
                    readfish_toml=self.readfish_toml,  # ToDo: This assumes a single sample per CNV analysis.
                    # NewBed=self.NewBed,
                    master_bed_tree=self.master_bed_tree,
                    **common_args,
                )
                await self.Fusion_panel.process_data()

        background_tasks.create(start_background_tasks())

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
            self.check_existing_bams_timer = app.timer(
                1, callback=self.check_existing_bams, once=True
            )

            # Add status tracking for the observer
            state.start_process("File Observer", ProcessType.SYSTEM)
            state.set_process_state("File Observer", ProcessState.STARTING)
            self.observer.start()
            state.set_process_state("File Observer", ProcessState.RUNNING)

            self.background_process_bams_timer = app.timer(
                1, callback=self.background_process_bams, once=True
            )

            logging.info("Watchfolder setup and added")

    async def pick_file(self) -> None:
        """
        Handle file picking for monitoring.

        :return: None
        """
        logging.info("Starting file picker in browse mode")
        ui.notify("Select a folder to monitor")
        result = await LocalFilePicker(f"{self.runsfolder}", multiple=True)
        if result:
            logging.info(f"Selected folder: {result}")
            ui.notify(f"You selected {result}")
            self.content.clear()
            with self.content:
                self.output, self.sampleID = os.path.split(result[0])
                logging.info(
                    f"Split path - output: {self.output}, sampleID: {self.sampleID}"
                )

                # Initialize storage for browse mode
                if self.sampleID not in app.storage.general[self.mainuuid]["samples"]:
                    logging.info(f"Initializing storage for sample {self.sampleID}")
                    app.storage.general[self.mainuuid]["samples"][self.sampleID] = {}
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "file_counters"
                    ] = {
                        "bam_passed": 0,
                        "bam_failed": 0,
                        "mapped_count": 0,
                        "unmapped_count": 0,
                        "pass_mapped_count": 0,
                        "fail_mapped_count": 0,
                        "pass_bases_count": 0,
                        "fail_bases_count": 0,
                        "bases_count": 0,
                    }
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "devices"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "basecall_models"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "run_time"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "flowcell_ids"
                    ] = []
                    app.storage.general[self.mainuuid]["samples"][self.sampleID][
                        "sample_ids"
                    ] = [self.sampleID]
                else:
                    logging.info(f"Storage already exists for sample {self.sampleID}")

                if (
                    self.sampleID
                    not in app.storage.general[self.mainuuid]["sample_list"]
                ):
                    logging.info(f"Adding {self.sampleID} to sample list")
                    app.storage.general[self.mainuuid]["sample_list"].append(
                        self.sampleID
                    )

                # Try to load existing data from master.csv if it exists
                master_csv = os.path.join(result[0], "master.csv")
                logging.info(f"Looking for master.csv at: {master_csv}")
                if os.path.exists(master_csv):
                    try:
                        logging.info("Found master.csv, loading data")
                        df = pd.read_csv(master_csv)
                        # Update storage with data from master.csv
                        for key in df.keys():
                            if (
                                key
                                in app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ]
                            ):
                                logging.info(f"Updating {key} from master.csv")
                                app.storage.general[self.mainuuid]["samples"][
                                    self.sampleID
                                ][key] = df[key].tolist()
                    except Exception as e:
                        logging.error(
                            f"Error loading master.csv: {str(e)}", exc_info=True
                        )
                else:
                    logging.warning(f"No master.csv found at {master_csv}")

                logging.info("Calling information_panel")
                await self.information_panel(sample_id=self.sampleID)
        else:
            logging.warning("No folder selected in file picker")

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
                                    ).style("font-size: 75%; font-weight: 100").classes(
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

    def create_dynamic_classification_chart(self, key, y_max=1):
        """Create a chart for a specific classification key (e.g., superfamily, family) with dynamic y-axis max."""
        return ui.echart(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": f"{key.title()} Classification Over Time",
                    "left": "center",
                    "top": 10,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "normal",
                        "color": "#000000",
                    },
                },
                "tooltip": {
                    "trigger": "axis",
                    "axisPointer": {"type": "line"},
                    "textStyle": {"fontSize": 14},
                },
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "5%",
                    "top": "25%",
                    "containLabel": True,
                },
                "legend": {
                    "type": "scroll",
                    "orient": "horizontal",
                    "top": 45,
                    "width": "90%",
                    "left": "center",
                    "textStyle": {"fontSize": 12, "color": "#666666"},
                    "pageButtonPosition": "end",
                    "pageButtonGap": 5,
                    "pageButtonItemGap": 5,
                    "pageIconColor": "#666666",
                    "pageIconInactiveColor": "#aaa",
                    "pageIconSize": 12,
                    "pageTextStyle": {"color": "#666666"},
                    "itemGap": 25,
                    "itemWidth": 14,
                    "itemHeight": 14,
                    "selectedMode": True,
                },
                "xAxis": {
                    "type": "time",
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{yyyy}-{MM}-{dd} {HH}:{mm}",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "yAxis": {
                    "type": "value",
                    "min": 0,
                    "max": y_max,
                    "interval": y_max / 5 if y_max > 0 else 1,
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{value}",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "series": [],
            }
        ).classes("w-full h-64")

    def create_mgmt_chart(self):
        """Create a chart for displaying MGMT methylation trends."""
        return ui.echart(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": "MGMT Methylation Over Time",
                    "left": "center",
                    "top": 10,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "normal",
                        "color": "#000000",
                    },
                },
                "tooltip": {
                    "trigger": "axis",
                    "axisPointer": {"type": "line"},
                    "textStyle": {"fontSize": 14},
                },
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "5%",
                    "top": "25%",
                    "containLabel": True,
                },
                "xAxis": {
                    "type": "time",
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{yyyy}-{MM}-{dd} {HH}:{mm}",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "yAxis": {
                    "type": "value",
                    "min": 0,
                    "max": 100,
                    "interval": 20,
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{value}%",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "series": [
                    {
                        "name": "Methylation",
                        "type": "line",
                        "smooth": True,
                        "animation": False,
                        "symbolSize": 6,
                        "emphasis": {
                            "focus": "series",
                            "itemStyle": {"borderWidth": 2},
                        },
                        "lineStyle": {"width": 2, "color": "#007AFF"},
                        "itemStyle": {"color": "#007AFF"},
                        "data": [],
                    }
                ],
            }
        ).classes("w-full h-64")

    def update_extracted_data_charts(self):
        """Update dynamic charts with data from extracted_data.json."""
        try:
            logging.info("Starting chart update")
            if not hasattr(self, "classification_charts") or not hasattr(
                self, "mgmt_chart"
            ):
                logging.warning("Charts not initialized")
                return

            json_path = os.path.join(
                self.output, "extracted_data", "extracted_data.json"
            )
            logging.info(f"Looking for data file at: {json_path}")

            if not os.path.exists(json_path):
                logging.warning(f"No data file found at {json_path}")
                return

            with self.safe_open(json_path, "r") as f:
                data = json.load(f)

            if not data:
                logging.warning("No data found in JSON file")
                return

            logging.info(f"Processing {len(data)} data points")

            # Build a dict: key -> label -> list of (timestamp, score)
            key_label_series = {}
            mgmt_data = []
            colors = [
                "#007AFF",
                "#34C759",
                "#FF9500",
                "#FF2D55",
                "#5856D6",
                "#FF3B30",
                "#5AC8FA",
                "#4CD964",
            ]
            for entry in data:
                timestamp = entry.get("timestamp")
                # Robustly handle missing classification
                for key, info in entry.get("classification", {}).items():
                    label = info.get("label")
                    score = info.get("score")
                    if label is None or score is None:
                        continue
                    if key not in key_label_series:
                        key_label_series[key] = {}
                    if label not in key_label_series[key]:
                        key_label_series[key][label] = []
                    key_label_series[key][label].append([timestamp, score])
                # Robustly handle missing mgmt_info
                mgmt = entry.get("mgmt_info", {})
                avg = mgmt.get("average_methylation")
                if avg is not None and timestamp is not None:
                    mgmt_data.append([timestamp, avg])
            # Sort all series by timestamp
            for key in key_label_series:
                for label in key_label_series[key]:
                    key_label_series[key][label].sort()
            mgmt_data.sort()
            logging.info(f"Found classification keys: {list(key_label_series.keys())}")
            for key in key_label_series:
                for label in key_label_series[key]:
                    logging.info(
                        f"Key '{key}', label '{label}': {len(key_label_series[key][label])} points"
                    )
            logging.info(f"Found {len(mgmt_data)} MGMT data points")

            # Dynamically update or create charts for each key
            for key, label_dict in key_label_series.items():
                # Find the max score for this key
                all_scores = [
                    score for points in label_dict.values() for _, score in points
                ]
                y_max = max(all_scores) if all_scores else 1
                y_max = max(1, int(y_max + 0.5))  # Round up for clarity
                if key not in self.classification_charts:
                    self.classification_charts[key] = (
                        self.create_dynamic_classification_chart(key, y_max=y_max)
                    )
                chart = self.classification_charts[key]
                # Update y-axis max and interval dynamically
                chart.options["yAxis"]["max"] = y_max
                chart.options["yAxis"]["interval"] = y_max / 5 if y_max > 0 else 1
                chart_series = []
                for idx, (label, data_points) in enumerate(label_dict.items()):
                    chart_series.append(
                        {
                            "name": label,
                            "type": "line",
                            "smooth": True,
                            "animation": False,
                            "symbolSize": 6,
                            "emphasis": {
                                "focus": "series",
                                "itemStyle": {"borderWidth": 2},
                            },
                            "lineStyle": {
                                "width": 2,
                                "color": colors[idx % len(colors)],
                            },
                            "itemStyle": {"color": colors[idx % len(colors)]},
                            "data": data_points,
                        }
                    )
                chart.options["series"] = chart_series
                chart.update()

            # Update MGMT chart
            if mgmt_data:
                logging.info("Updating MGMT chart with %d points", len(mgmt_data))
                self.mgmt_chart.options["series"][0]["data"] = mgmt_data
                self.mgmt_chart.update()
        except Exception as e:
            logging.error(
                f"Error updating extracted data charts: {str(e)}", exc_info=True
            )

    async def information_panel(self, sample_id=None):
        """Display the main information panel with analysis results."""
        try:
            logging.info(f"Starting information_panel with sample_id: {sample_id}")
            if sample_id:
                self.sampleID = sample_id
                logging.info(f"Set self.sampleID to {self.sampleID}")
            else:
                logging.warning("No sample_id provided")
                if self.browse:
                    ui.label("Please select a sample folder to view").classes(
                        "text-xl text-gray-600 my-4"
                    )
                    return

            if not self.browse:
                logging.info("Live mode - showing sample list")
                self.show_list()
                app.storage.general[self.mainuuid]["sample_list"].on_change(
                    self.show_list.refresh
                )

            if (
                self.sampleID
                and not self.browse
                and self.sampleID not in app.storage.general[self.mainuuid]["samples"]
            ):
                logging.error(f"Sample {self.sampleID} not found in storage")
                ui.notify(f"Sample {self.sampleID} not found")
                ui.navigate.to("/live")
                return

            if not self.sampleID:
                logging.warning("No sample ID set, cannot render panel")
                return

            logging.info("Starting UI rendering")
            if self.sampleID in app.storage.general[self.mainuuid]["samples"]:
                logging.info(
                    f"Storage state for sample {self.sampleID}: {app.storage.general[self.mainuuid]['samples'][self.sampleID]}"
                )
            else:
                logging.warning(f"No storage found for sample {self.sampleID}")
                return

            with ui.column().classes("w-full px-6"):
                logging.info("Rendering title section")
                # Title and description
                with ui.column().classes("space-y-2 mb-4"):
                    ui.label("CNS Tumor Methylation Classification").classes(
                        "text-2xl font-medium text-gray-900"
                    )
                    ui.label(
                        "This tool enables classification of brain tumors in real time from Oxford Nanopore Data."
                    ).classes("text-gray-600")
                logging.info("Title section rendered")

                # Run Information Card
                with ui.card().classes("w-full p-4 bg-white rounded-lg shadow-sm"):
                    ui.label("Run Information").classes(
                        "text-lg font-medium text-gray-900 mb-4"
                    )
                    with ui.row().classes("items-center gap-6"):
                        # Run Time
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "run_time"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("schedule").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "run_time",
                                    backward=lambda n: (
                                        f"Run: {parser.parse(n[0]).strftime('%Y-%m-%d %H:%M')}"
                                        if n
                                        else None
                                    ),
                                )

                        # Model
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "basecall_models"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("model_training").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "basecall_models",
                                    backward=lambda n: f"Model: {n[0]}" if n else None,
                                )

                        # Device
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "devices"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("memory").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "devices",
                                    backward=lambda n: f"Device: {n[0]}" if n else None,
                                )

                        # Flow Cell
                        if app.storage.general[self.mainuuid]["samples"][self.sampleID][
                            "flowcell_ids"
                        ]:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("grid_4x4").classes("text-gray-400")
                                ui.label().bind_text_from(
                                    app.storage.general[self.mainuuid]["samples"][
                                        self.sampleID
                                    ],
                                    "flowcell_ids",
                                    backward=lambda n: (
                                        f"Flow Cell: {n[0]}" if n else None
                                    ),
                                )

                        # Sample
                        if self.sampleID:
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("label").classes("text-gray-400")
                                ui.label(f"Sample: {self.sampleID}").classes(
                                    "text-gray-600"
                                )

                # Results Summary Section
                with ui.column().classes("w-full space-y-4"):
                    # Classification Results - Single row with equal width columns
                    # Only show if at least one classifier is not excluded
                    if not all(
                        classifier in self.exclude
                        for classifier in ["sturgeon", "nanodx", "pannanodx", "forest"]
                    ):
                        with ui.card().classes(
                            "w-full p-3 sm:p-4 bg-gray-50 rounded-lg"
                        ):
                            ui.label("Classification Summary").classes(
                                "font-medium text-gray-900 mb-3"
                            )
                            with ui.grid().classes(
                                "grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4"
                            ):
                                if "sturgeon" not in self.exclude:
                                    sturgeonsummary = ui.column().classes("space-y-2")
                                if "nanodx" not in self.exclude:
                                    nanodxsummary = ui.column().classes("space-y-2")
                                if "pannanodx" not in self.exclude:
                                    pannanodxsummary = ui.column().classes("space-y-2")
                                if "forest" not in self.exclude:
                                    forestsummary = ui.column().classes("space-y-2")

                with ui.column().classes("w-full space-y-4"):
                    # Analysis Results - Two columns grid
                    with ui.grid().classes(
                        "grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4"
                    ):
                        if "coverage" not in self.exclude:
                            coverage = ui.column().classes("space-y-1")
                        if "cnv" not in self.exclude:
                            cnvsummary = ui.column().classes("space-y-1")
                        if "mgmt" not in self.exclude:
                            mgmt = ui.column().classes("space-y-1")
                        if "fusion" not in self.exclude:
                            fusions = ui.column().classes("space-y-1")

                        # Add summary containers for each analysis type above the tab panel
                        if self.mnpflex_config and all(self.mnpflex_config.values()):
                            mnpflexsummary = ui.column().classes("space-y-1")

                    # Add monitoring information panel - Now positioned after diagnosis cards
                    # Only show in live mode
                    if not self.browse:
                        with ui.card().classes("w-full bg-white rounded-lg shadow-sm"):
                            with ui.expansion(
                                "Monitoring Information", icon="folder"
                            ).classes("w-full") as monitoring_expansion:
                                # Add summary header content
                                with monitoring_expansion.add_slot("header"):
                                    with ui.row().classes(
                                        "items-center justify-between w-full"
                                    ):
                                        # Left side - File counts
                                        with ui.row().classes("items-center gap-4"):
                                            with ui.row().classes("items-center gap-1"):
                                                ui.icon("description").classes(
                                                    "text-blue-600"
                                                )
                                                ui.label().bind_text_from(
                                                    app.storage.general[self.mainuuid][
                                                        "bam_count"
                                                    ],
                                                    "total_files",
                                                    backward=lambda n: (
                                                        f"{n} total" if n else 0
                                                    ),
                                                ).classes("text-sm text-gray-600")

                                            with ui.row().classes("items-center gap-1"):
                                                ui.icon("check_circle").classes(
                                                    "text-green-600"
                                                )
                                                ui.label().bind_text_from(
                                                    app.storage.general[self.mainuuid][
                                                        "samples"
                                                    ][self.sampleID]["file_counters"],
                                                    "bam_passed",
                                                    backward=lambda n: (
                                                        f"{n} passed" if n else 0
                                                    ),
                                                ).classes("text-sm text-gray-600")

                                            with ui.row().classes("items-center gap-1"):
                                                ui.icon("error").classes("text-red-600")
                                                fail_count = ui.label().classes(
                                                    "text-sm text-gray-600"
                                                )

                                                def update_fail_count():
                                                    count = app.storage.general[
                                                        self.mainuuid
                                                    ]["samples"][self.sampleID][
                                                        "file_counters"
                                                    ][
                                                        "bam_failed"
                                                    ]
                                                    fail_count.text = (
                                                        f"{count:,} failed"
                                                        if count is not None
                                                        else "--"
                                                    )

                                                app.storage.general[self.mainuuid][
                                                    "samples"
                                                ][self.sampleID][
                                                    "file_counters"
                                                ].on_change(
                                                    update_fail_count
                                                )
                                                update_fail_count()

                                        # Right side - Title and Sample Name
                                        with ui.row().classes("items-center gap-4"):
                                            ui.label("Monitoring Information:").classes(
                                                "text-sm font-medium text-gray-700"
                                            )
                                            ui.label(self.sampleID).classes(
                                                "text-sm text-gray-600"
                                            )

                                # Existing content
                                with ui.column().classes("space-y-4 p-4"):
                                    # File Paths Section
                                    ui.label("File Paths").classes(
                                        "text-xl font-medium text-gray-900"
                                    )
                                    with ui.column().classes("space-y-4 ml-1"):
                                        # Monitoring Path
                                        with ui.row().classes("items-start gap-2"):
                                            ui.icon("folder_open").classes(
                                                "text-blue-600 shrink-0 mt-1"
                                            )
                                            with ui.column().classes(
                                                "flex-grow min-w-0"
                                            ):
                                                ui.label("Monitoring:").classes(
                                                    "text-gray-600"
                                                )
                                                ui.label(f"{self.watchfolder}").classes(
                                                    "break-all text-gray-900 font-mono"
                                                )

                                        # Output Path
                                        with ui.row().classes("items-start gap-2"):
                                            ui.icon("output").classes(
                                                "text-blue-600 shrink-0 mt-1"
                                            )
                                            with ui.column().classes(
                                                "flex-grow min-w-0"
                                            ):
                                                ui.label("Output:").classes(
                                                    "text-gray-600"
                                                )
                                                ui.label(f"{self.output}").classes(
                                                    "break-all text-gray-900 font-mono"
                                                )

                                    # Summary Statistics Section - Grid layout for BAM File Summary and Sequencing Statistics
                                    with ui.grid().classes("grid-cols-6 gap-6 mt-6"):
                                        # BAM File Summary Card - Takes 1 column (1/4)
                                        with ui.card().classes(
                                            "col-span-1 p-4 bg-gray-50 rounded-lg"
                                        ):
                                            ui.label("BAM File Summary").classes(
                                                "text-lg font-medium text-gray-900 mb-4"
                                            )
                                            with ui.column().classes("space-y-3"):
                                                # Total Files
                                                with ui.row().classes(
                                                    "items-center justify-between w-full"
                                                ):
                                                    with ui.row().classes(
                                                        "items-center gap-2"
                                                    ):
                                                        ui.icon("description").classes(
                                                            "text-blue-600 shrink-0"
                                                        )
                                                        ui.label(
                                                            "Total Files:"
                                                        ).classes("text-gray-600")
                                                    total_files = ui.label().classes(
                                                        "font-medium text-gray-900"
                                                    )

                                                    def update_total_files():
                                                        total = app.storage.general[
                                                            self.mainuuid
                                                        ]["bam_count"]["total_files"]
                                                        total_files.text = (
                                                            f"{total:,}"
                                                            if total > 0
                                                            else "--"
                                                        )

                                                    app.storage.general[self.mainuuid][
                                                        "bam_count"
                                                    ].on_change(update_total_files)
                                                    update_total_files()

                                                # Passed Files
                                                with ui.row().classes(
                                                    "items-center justify-between w-full"
                                                ):
                                                    with ui.row().classes(
                                                        "items-center gap-2"
                                                    ):
                                                        ui.icon("check_circle").classes(
                                                            "text-green-600 shrink-0"
                                                        )
                                                        ui.label("Passed:").classes(
                                                            "text-gray-600"
                                                        )
                                                    passed = ui.label().classes(
                                                        "font-medium text-gray-900"
                                                    )

                                                    def update_passed():
                                                        passed.text = (
                                                            f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_passed']:,}"
                                                            if app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ][
                                                                "bam_passed"
                                                            ]
                                                            is not None
                                                            else "--"
                                                        )

                                                    app.storage.general[self.mainuuid][
                                                        "samples"
                                                    ][self.sampleID][
                                                        "file_counters"
                                                    ].on_change(
                                                        update_passed
                                                    )
                                                    update_passed()

                                                # Failed Files
                                                with ui.row().classes(
                                                    "items-center justify-between w-full"
                                                ):
                                                    with ui.row().classes(
                                                        "items-center gap-2"
                                                    ):
                                                        ui.icon("error").classes(
                                                            "text-red-600 shrink-0"
                                                        )
                                                        ui.label("Failed:").classes(
                                                            "text-gray-600"
                                                        )
                                                    failed = ui.label().classes(
                                                        "font-medium text-gray-900"
                                                    )

                                                    def update_failed():
                                                        failed.text = (
                                                            f"{app.storage.general[self.mainuuid]['samples'][self.sampleID]['file_counters']['bam_failed']:,}"
                                                            if app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ][
                                                                "bam_failed"
                                                            ]
                                                            is not None
                                                            else "--"
                                                        )

                                                    app.storage.general[self.mainuuid][
                                                        "samples"
                                                    ][self.sampleID][
                                                        "file_counters"
                                                    ].on_change(
                                                        update_failed
                                                    )
                                                    update_failed()

                                        # Sequencing Statistics Card - Takes 3 columns (3/4)
                                        with ui.card().classes(
                                            "col-span-3 p-4 bg-gray-50 rounded-lg"
                                        ):
                                            ui.label("Sequencing Statistics").classes(
                                                "text-lg font-medium text-gray-900 mb-4"
                                            )
                                            with ui.grid().classes("grid-cols-4 gap-4"):
                                                # Read Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label(
                                                        "Base Count by Read Type"
                                                    ).classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Total Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Total bases in reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            total_reads = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_total_reads():
                                                                mapped = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "mapped_bases"
                                                                    ]
                                                                )
                                                                unmapped = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "unmapped_bases"
                                                                    ]
                                                                )
                                                                total = (
                                                                    mapped + unmapped
                                                                    if mapped
                                                                    is not None
                                                                    and unmapped
                                                                    is not None
                                                                    else 0
                                                                )
                                                                total_reads.text = (
                                                                    f"{total:,} bp"
                                                                    if total > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_total_reads
                                                            )
                                                            update_total_reads()

                                                        # Mapped Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Bases in mapped reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mapped():
                                                                count = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "mapped_bases"
                                                                    ]
                                                                )
                                                                mapped.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mapped
                                                            )
                                                            update_mapped()

                                                        # Unmapped Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Bases in unmapped reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_unmapped():
                                                                count = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "unmapped_bases"
                                                                    ]
                                                                )
                                                                unmapped.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_unmapped
                                                            )
                                                            update_unmapped()

                                                # Read Quality Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label(
                                                        "Read Quality Statistics"
                                                    ).classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Total Mapped
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Total mapped reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            total_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_total_mapped():
                                                                pass_mapped = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_count"
                                                                ]
                                                                fail_mapped = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_count"
                                                                ]
                                                                total = (
                                                                    pass_mapped
                                                                    + fail_mapped
                                                                    if pass_mapped
                                                                    is not None
                                                                    and fail_mapped
                                                                    is not None
                                                                    else 0
                                                                )
                                                                total_mapped.text = (
                                                                    f"{total:,} reads"
                                                                    if total > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_total_mapped
                                                            )
                                                            update_total_mapped()

                                                        # High Quality Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "High quality reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            pass_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_pass_mapped():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_count"
                                                                ]
                                                                pass_mapped.text = (
                                                                    f"{count:,} reads"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_pass_mapped
                                                            )
                                                            update_pass_mapped()

                                                        # Low Quality Reads
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Low quality reads:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            fail_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_fail_mapped():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_count"
                                                                ]
                                                                fail_mapped.text = (
                                                                    f"{count:,} reads"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_fail_mapped
                                                            )
                                                            update_fail_mapped()

                                                # Base Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label("Base Statistics").classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Total Bases
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Total bases:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            total_bases = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_total_bases():
                                                                pass_bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_bases_count"
                                                                ]
                                                                fail_bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_bases_count"
                                                                ]
                                                                total = (
                                                                    pass_bases
                                                                    + fail_bases
                                                                    if pass_bases
                                                                    is not None
                                                                    and fail_bases
                                                                    is not None
                                                                    else 0
                                                                )
                                                                total_bases.text = (
                                                                    f"{total:,} bp"
                                                                    if total > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_total_bases
                                                            )
                                                            update_total_bases()

                                                        # High Quality Bases
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "High quality bases:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            pass_bases = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_pass_bases():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_bases_count"
                                                                ]
                                                                pass_bases.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_pass_bases
                                                            )
                                                            update_pass_bases()

                                                        # Low Quality Bases
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Low quality bases:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            fail_bases = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_fail_bases():
                                                                count = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_bases_count"
                                                                ]
                                                                fail_bases.text = (
                                                                    f"{count:,} bp"
                                                                    if count is not None
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_fail_bases
                                                            )
                                                            update_fail_bases()

                                                # Read Length Statistics
                                                with ui.card().classes(
                                                    "w-full p-2 sm:p-3 bg-white rounded-lg"
                                                ):
                                                    ui.label(
                                                        "Read Length Statistics"
                                                    ).classes(
                                                        "text-sm font-medium text-gray-700 mb-2"
                                                    )
                                                    with ui.column().classes(
                                                        "space-y-1.5"
                                                    ):
                                                        # Overall Mean Lengths
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean mapped length:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_mapped():
                                                                bases = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "mapped_bases"
                                                                    ]
                                                                )
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "mapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_mapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_mapped
                                                            )
                                                            update_mean_mapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean unmapped length:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_unmapped():
                                                                bases = (
                                                                    app.storage.general[
                                                                        self.mainuuid
                                                                    ]["samples"][
                                                                        self.sampleID
                                                                    ][
                                                                        "file_counters"
                                                                    ][
                                                                        "unmapped_bases"
                                                                    ]
                                                                )
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "unmapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_unmapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_unmapped
                                                            )
                                                            update_mean_unmapped()

                                                        # Pass/Fail Mean Lengths
                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean pass mapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_pass_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_pass_mapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_mapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_pass_mapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_pass_mapped
                                                            )
                                                            update_mean_pass_mapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean fail mapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_fail_mapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_fail_mapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_mapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_fail_mapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_fail_mapped
                                                            )
                                                            update_mean_fail_mapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean pass unmapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_pass_unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_pass_unmapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_unmapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "pass_unmapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_pass_unmapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_pass_unmapped
                                                            )
                                                            update_mean_pass_unmapped()

                                                        with ui.row().classes(
                                                            "justify-between items-center"
                                                        ):
                                                            ui.label(
                                                                "Mean fail unmapped:"
                                                            ).classes(
                                                                "text-gray-600 text-sm"
                                                            )
                                                            mean_fail_unmapped = ui.label().classes(
                                                                "font-medium text-gray-900 text-sm"
                                                            )

                                                            def update_mean_fail_unmapped():
                                                                bases = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_unmapped_bases"
                                                                ]
                                                                reads = app.storage.general[
                                                                    self.mainuuid
                                                                ][
                                                                    "samples"
                                                                ][
                                                                    self.sampleID
                                                                ][
                                                                    "file_counters"
                                                                ][
                                                                    "fail_unmapped_reads_num"
                                                                ]
                                                                mean = (
                                                                    bases / reads
                                                                    if reads > 0
                                                                    else 0
                                                                )
                                                                mean_fail_unmapped.text = (
                                                                    f"{mean:,.0f} bp"
                                                                    if mean > 0
                                                                    else "--"
                                                                )

                                                            app.storage.general[
                                                                self.mainuuid
                                                            ]["samples"][self.sampleID][
                                                                "file_counters"
                                                            ].on_change(
                                                                update_mean_fail_unmapped
                                                            )
                                                            update_mean_fail_unmapped()

                # Detailed Analysis Tabs
                if sample_id:
                    selectedtab = None
                    with ui.tabs().classes(
                        "w-full"
                    ) as tabs:  # .on('click', lambda: ui.notify('You clicked a tab.'))
                        if not (
                            set(["sturgeon", "pannanodx", "nanodx", "forest"]).issubset(
                                set(self.exclude)
                            )
                        ):
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
                        if self.mnpflex_config and all(self.mnpflex_config.values()):
                            mnpflextab = ui.tab("MNPFlex")
                            if not selectedtab:
                                selectedtab = mnpflextab

                    with ui.tab_panels(tabs, value=selectedtab).classes(
                        "w-full"
                    ):  # , on_change=lambda e: print(e))
                        display_args = {
                            "threads": self.threads,
                            "output": self.output,
                            "progress": True,
                            "browse": self.browse,
                            "uuid": self.mainuuid,
                            "force_sampleid": self.force_sampleid,
                        }
                        if not (
                            set(["sturgeon", "nanodx", "pannanodx", "forest"]).issubset(
                                set(self.exclude)
                            )
                        ):
                            with ui.tab_panel(methylationtab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    ui.label("Methylation Classifications").classes(
                                        "text-sky-600 dark:text-white"
                                    ).style(
                                        "font-size: 150%; font-weight: 300"
                                    ).tailwind(
                                        "drop-shadow", "font-bold"
                                    )
                                    if "sturgeon" not in self.exclude:
                                        self.Sturgeon = SturgeonVis(
                                            analysis_name="STURGEON",
                                            batch=True,
                                            summary=sturgeonsummary,
                                            **display_args,
                                        )
                                        await self.Sturgeon.render_ui(
                                            sample_id=self.sampleID
                                        )
                                    if "nanodx" not in self.exclude:
                                        self.NanoDX = NanoDXVis(
                                            analysis_name="NANODX",
                                            batch=True,
                                            summary=nanodxsummary,
                                            **display_args,
                                        )
                                        await self.NanoDX.render_ui(
                                            sample_id=self.sampleID
                                        )
                                    if "pannanodx" not in self.exclude:
                                        self.PanNanoDX = NanoDXVis(
                                            analysis_name="PANNANODX",
                                            batch=True,
                                            summary=pannanodxsummary,
                                            model="pancan_devel_v5i_NN.pkl",
                                            **display_args,
                                        )
                                        await self.PanNanoDX.render_ui(
                                            sample_id=self.sampleID
                                        )
                                    if "forest" not in self.exclude:
                                        self.RandomForest = RandomForestVis(
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
                                    self.CNV = CNVVis(
                                        analysis_name="CNV",
                                        summary=cnvsummary,
                                        target_panel=self.target_panel,
                                        reference_file=self.reference,
                                        bed_file=self.bed_file,
                                        readfish_toml=self.readfish_toml,  # ToDo: This assumes a single sample per CNV analysis.
                                        # NewBed=self.NewBed, #ToDo: This assumes a single sample per CNV analysis.
                                        master_bed_tree=self.master_bed_tree,
                                        **display_args,
                                    )
                                    await self.CNV.render_ui(sample_id=self.sampleID)

                        if "coverage" not in self.exclude:
                            with ui.tab_panel(coveragetab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.Target_Coverage = TargetCoverageVis(
                                        analysis_name="COVERAGE",
                                        summary=coverage,
                                        target_panel=self.target_panel,
                                        reference=self.reference,
                                        enable_snp_calling=self.enable_snp_calling
                                        and self.reference is not None,
                                        **display_args,
                                    )
                                    await self.Target_Coverage.render_ui(
                                        sample_id=self.sampleID
                                    )

                        if "mgmt" not in self.exclude:
                            with ui.tab_panel(mgmttab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.MGMT_panel = MGMTVis(
                                        analysis_name="MGMT",
                                        summary=mgmt,
                                        **display_args,
                                    )
                                    await self.MGMT_panel.render_ui(
                                        sample_id=self.sampleID
                                    )

                        if "fusion" not in self.exclude:
                            with ui.tab_panel(fusionstab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.Fusion_panel = FusionVis(
                                        analysis_name="FUSION",
                                        summary=fusions,
                                        target_panel=self.target_panel,
                                        reference_file=self.reference,
                                        bed_file=self.bed_file,
                                        readfish_toml=self.readfish_toml,  # ToDo: This assumes a single sample per CNV analysis.
                                        # NewBed=self.NewBed,
                                        master_bed_tree=self.master_bed_tree,
                                        **display_args,
                                    )
                                    await self.Fusion_panel.render_ui(
                                        sample_id=self.sampleID
                                    )

                        if self.mnpflex_config and all(self.mnpflex_config.values()):
                            with ui.tab_panel(mnpflextab).classes("w-full"):
                                with ui.card().classes("rounded w-full"):
                                    self.MNPFlex = MNPFlex_Object(
                                        analysis_name="MNPFLEX",
                                        mnpflex_config=self.mnpflex_config,
                                        summary=mnpflexsummary,
                                        **display_args,
                                    )
                                    await self.MNPFlex.render_ui(
                                        sample_id=self.sampleID
                                    )

                    async def download_mnpflex_report(self, report_name):
                        """Download the MNPFlex report."""
                        report_path = os.path.join(
                            self.output, self.sampleID, report_name
                        )
                        if os.path.exists(report_path):
                            ui.download(report_path)
                        else:
                            ui.notify(
                                f"Report {report_name} not found", type="negative"
                            )

                    async def confirm_report_generation():
                        """
                        Show a confirmation dialog before generating the report.
                        """
                        report_types = {
                            "summary": "Summary Only",
                            "detailed": "Detailed",
                        }
                        state = {"type": "detailed"}  # Store state in a dictionary

                        with ui.dialog() as dialog, ui.card().classes("w-96 p-4"):
                            ui.label("Generate Report").classes(
                                "text-h6 font-bold mb-4"
                            )

                            # Report type selector
                            with ui.column().classes("mb-4"):
                                ui.label("Report Type").classes("font-bold mb-2")
                                ui.toggle(
                                    report_types,
                                    value="detailed",
                                    on_change=lambda e: state.update({"type": e.value}),
                                )

                            # Disclaimer section
                            with ui.column().classes("mb-4"):
                                ui.label("Disclaimer").classes("font-bold mb-2")
                                # Split into paragraphs and join with double breaks for better spacing
                                formatted_text = EXTENDED_DISCLAIMER_TEXT.replace(
                                    "\n\n", "<br><br>"
                                ).replace("\n", " ")
                                ui.label(formatted_text).classes(
                                    "text-sm text-gray-600 mb-4"
                                )

                            ui.label(
                                "Are you sure you want to generate a report?"
                            ).classes("mb-4")

                            # Buttons in a row at the bottom
                            with ui.row().classes("justify-end gap-2"):
                                ui.button(
                                    "No", on_click=lambda: dialog.submit(("No", None))
                                ).props("flat")
                                ui.button(
                                    "Yes",
                                    on_click=lambda: dialog.submit(
                                        ("Yes", state["type"])
                                    ),
                                ).props("color=primary")

                        dialog_result = await dialog
                        if dialog_result is None:
                            return  # Dialog was closed without submitting

                        result, report_type = dialog_result
                        if result == "Yes":
                            await download_report(report_type)

                    async def download_report(report_type="detailed"):
                        """
                        Generate and download the report.

                        :param report_type: Type of report to generate ('summary' or 'detailed')
                        :return: None
                        """
                        ui.notify("Generating Report")
                        if not self.browse:
                            for item in app.storage.general[self.mainuuid]:
                                if item == "sample_ids":
                                    for sample in app.storage.general[self.mainuuid][
                                        item
                                    ]:
                                        self.sampleID = sample
                        if (
                            self.browse
                        ):  # This is correct use of io_bound as we are in the main thread.
                            myfile = await run.io_bound(
                                create_pdf,
                                f"{self.sampleID}_run_report.pdf",
                                self.check_and_create_folder(
                                    self.output, self.sampleID
                                ),
                                report_type,
                            )
                        else:
                            myfile = await run.io_bound(
                                create_pdf,
                                f"{self.sampleID}_run_report.pdf",
                                self.output,
                                report_type,
                            )
                        ui.download(myfile)
                        ui.notify("Report Downloaded")

                    # Replace the direct download_report button with the confirmation dialog
                    ui.button(
                        "Generate Report",
                        on_click=confirm_report_generation,
                        icon="download",
                    )

                    # Add report generation button
        except Exception as e:
            logging.error(f"Error rendering analysis UI: {str(e)}", exc_info=True)
            ui.notify(f"Error rendering analysis UI: {str(e)}", type="error")

    async def background_process_bams(self):
        """
        Background task to process BAM files.

        :return: None
        """
        logging.info("Starting background BAM processing")

        # Set up a timer to check application state every 10 seconds
        # This timer monitors for shutdown events and manages process states
        self.app_state_timer = app.timer(10, self.check_app_state)

        # Set up a timer to process BAM files every 10 seconds
        # This handles individual BAM file processing in a serial manner
        self.process_bams_tracker = app.timer(10, self.process_bams)

        # Register the serial BAM analysis process in the state management system
        # This allows the UI to track and display its status
        state.start_process("Serial Bam Analysis", ProcessType.BACKGROUND)

        # Set the initial state of serial BAM analysis to waiting for data
        # This indicates the process is ready but not actively processing yet
        state.set_process_state("Serial Bam Analysis", ProcessState.WAITING_FOR_DATA)

        # Set up a timer to handle merging of BAM files every 10 seconds
        # This manages the merging of multiple BAM files into consolidated files
        self.process_bigbadmerge_tracker = app.timer(10, self.process_bigbadmerge)

        # Register the merge BAM analysis process in the state management system
        # This allows the UI to track and display its status
        state.start_process("Merge Bam Analysis", ProcessType.BACKGROUND)

        # Set the initial state of merge BAM analysis to waiting for data
        # This indicates the process is ready but not actively merging yet
        state.set_process_state("Merge Bam Analysis", ProcessState.WAITING_FOR_DATA)

        # Define a list of external processes to monitor for memory usage
        # These are important system processes that run alongside the BAM processing
        otherprocs = ["dorado", "minknow", "docker"]

        # Set up a timer to check and log memory usage every 15 seconds
        # This monitors the memory consumption of both the application and external processes
        self.check_and_log_memory = app.timer(
            15, lambda: self.get_memory_usage(process_names_to_monitor=otherprocs)
        )

    async def check_app_state(self):
        if state.shutdown_event:
            print("shutdown_event is True, stopping background processes")
            if hasattr(self, "observer"):
                state.set_process_state("File Observer", ProcessState.STOPPING)
                state.stop_process("File Observer")
                self.observer.stop()
                self.observer.join()
            print("Observer stopped")
            while state.get_running_process_count() > 0:
                await asyncio.sleep(1)
                print(
                    f"Waiting for {state.get_running_process_count()} processes to finish"
                )
                print(
                    f"Running processes: {[p for p, s in state.process_states.items() if s == ProcessState.RUNNING]}"
                )
            self.finished = True
            self.shutdown_background()
        else:
            # Periodic cleanup of old temporary directories (every 10 checks = 100 seconds)
            if not hasattr(self, '_cleanup_counter'):
                self._cleanup_counter = 0
            self._cleanup_counter += 1
            
            if self._cleanup_counter >= 10:
                logging.info("Performing periodic cleanup of old temporary directories")
                self.cleanup_old_temp_directories()
                self.log_file_handle_stats()  # Log file handle statistics
                self._cleanup_counter = 0

    def check_and_create_folder(self, path, folder_name=None):
        """Check if a folder exists and create it if it doesn't."""
        if folder_name:
            full_path = os.path.join(path, folder_name)
        else:
            full_path = path
        os.makedirs(full_path, exist_ok=True)
        return full_path

    async def download_mnpflex_report(self, report_name):
        """Download the MNPFlex report."""
        report_path = os.path.join(self.output, self.sampleID, report_name)
        if os.path.exists(report_path):
            ui.download(report_path)
        else:
            ui.notify(f"Report {report_name} not found", type="negative")

    def get_memory_usage(self, process_names_to_monitor=None):
        """
        Get memory usage statistics for the current process and optionally other specified processes.

        Args:
            process_names_to_monitor (list): Optional list of process names to monitor (e.g., ['modkit', 'samtools'])

        Returns:
            dict: Dictionary containing memory usage statistics
        """
        try:
            # Get the current process
            # print('a. links:', len(binding.active_links))
            # print(binding.active_links)

            current_process = psutil.Process(os.getpid())
            memory = psutil.virtual_memory()
            available_memory_gb = memory.available / (1024**3)
            total_memory_gb = memory.total / (1024**3)
            mem_info = current_process.memory_info()
            rss = mem_info.rss / (1024**3)
            vms = mem_info.vms / (1024**3)

            # Create memory usage data dictionary for current process
            memory_data = {
                "system": {
                    "total": total_memory_gb,
                    "avail": available_memory_gb,
                },
                "main_process": {
                    "name": current_process.name(),
                    "pid": current_process.pid,
                    "rss": rss,  # Resident Set Size (physical memory)
                    "vms": vms,  # Virtual Memory Size
                    "cpu_percent": current_process.cpu_percent(),
                },
                "monitored_processes": {},
            }

            # Monitor other processes if specified
            if process_names_to_monitor:
                for proc in psutil.process_iter(
                    ["pid", "name", "memory_info", "cpu_percent"]
                ):
                    try:
                        for monitor_name in process_names_to_monitor:
                            if monitor_name.lower() in proc.info["name"].lower():
                                if proc.info["memory_info"] is not None:
                                    proc_mem = proc.info["memory_info"]
                                    proc_rss = proc_mem.rss / (1024**3)
                                    proc_vms = proc_mem.vms / (1024**3)

                                    # Add or append to process list in memory data
                                    if (
                                        monitor_name
                                        not in memory_data["monitored_processes"]
                                    ):
                                        memory_data["monitored_processes"][
                                            monitor_name
                                        ] = []

                                    memory_data["monitored_processes"][
                                        monitor_name
                                    ].append(
                                        {
                                            "pid": proc.info["pid"],
                                            "rss": proc_rss,
                                            "vms": proc_vms,
                                            "cpu_percent": proc.info["cpu_percent"],
                                        }
                                    )
                    except (
                        psutil.NoSuchProcess,
                        psutil.AccessDenied,
                        psutil.ZombieProcess,
                    ):
                        continue

            # Create memory log file path
            memory_log_file = os.path.join(
                self.output, f"memory_usage_{self.mainuuid}.csv"
            )

            # Get current timestamp
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Create header if file doesn't exist
            if not os.path.exists(memory_log_file):
                header = [
                    "timestamp",
                    "total_memory_gb",
                    "available_memory_gb",
                    "main_process_name",
                    "main_process_pid",
                    "main_process_rss_gb",
                    "main_process_vms_gb",
                    "main_process_cpu_percent",
                ]

                if process_names_to_monitor:
                    for proc_name in process_names_to_monitor:
                        header.extend(
                            [
                                f"{proc_name}_pid",
                                f"{proc_name}_rss_gb",
                                f"{proc_name}_vms_gb",
                                f"{proc_name}_cpu_percent",
                            ]
                        )

                with self.safe_open(memory_log_file, "w") as f:
                    f.write(",".join(header) + "\n")

            # Prepare log entry
            log_data = [
                timestamp,
                f"{total_memory_gb:.2f}",
                f"{available_memory_gb:.2f}",
                current_process.name(),
                str(current_process.pid),
                f"{rss:.2f}",
                f"{vms:.2f}",
                f"{current_process.cpu_percent():.1f}",
            ]

            # Add monitored processes data
            if process_names_to_monitor:
                for proc_name in process_names_to_monitor:
                    if (
                        proc_name in memory_data["monitored_processes"]
                        and memory_data["monitored_processes"][proc_name]
                    ):
                        # Use the first instance if multiple processes with same name
                        proc_data = memory_data["monitored_processes"][proc_name][0]
                        log_data.extend(
                            [
                                str(proc_data["pid"]),
                                f"{proc_data['rss']:.2f}",
                                f"{proc_data['vms']:.2f}",
                                f"{proc_data['cpu_percent']:.1f}",
                            ]
                        )
                    else:
                        # If process not found, add placeholder values
                        log_data.extend(["NA", "0.00", "0.00", "0.0"])

            # Append memory data with timestamp
            with self.safe_open(memory_log_file, "a") as f:
                f.write(",".join(log_data) + "\n")

            return memory_data

        except psutil.NoSuchProcess:
            logging.error("Error: Process not found.")
            return None
        except Exception as e:
            logging.error(f"Error getting memory usage: {str(e)}")
            return None

    async def process_bigbadmerge(self):
        self.process_bigbadmerge_tracker.active = False
        # Dictionary to store files by sample ID
        files_by_sample = {}
        latest_files = {}  # Track latest file time per sample
        state.start_process("Merge Bam Analysis", ProcessType.BACKGROUND)
        state.set_process_state("Merge Bam Analysis", ProcessState.WAITING_FOR_DATA)
        
        # Stopping Analysis to check for memory
        
        analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
        # Process queue and organize files by sample ID
        while self.bamforbigbadmerge.qsize() > 0:
            state.set_process_state("Merge Bam Analysis", ProcessState.RUNNING)
            file, filetime, sampleID = self.bamforbigbadmerge.get()

            # Validate sampleID - if None, use a fallback
            if sampleID is None:
                sampleID = f"unknown_sample_{os.path.basename(file)}"
                logging.warning(
                    f"BAM file {file} has no sampleID in queue, using fallback: {sampleID}"
                )

            # Initialize containers for new sample IDs
            if sampleID not in files_by_sample:
                files_by_sample[sampleID] = []
                latest_files[sampleID] = 0
                # Create temporary directories if needed
                if sampleID not in self.dataDir:
                    logging.info(f"Creating temporary directories for sample {sampleID}")
                    data_temp_dir = tempfile.TemporaryDirectory(
                        dir=self.check_and_create_folder(self.output, sampleID)
                    )
                    bed_temp_dir = tempfile.TemporaryDirectory(
                        dir=self.check_and_create_folder(self.output, sampleID)
                    )
                    
                    # Register with file handle manager
                    self.file_handle_manager.register_temp_directory(data_temp_dir)
                    self.file_handle_manager.register_temp_directory(bed_temp_dir)
                    
                    self.dataDir[sampleID] = data_temp_dir
                    self.bedDir[sampleID] = bed_temp_dir
                    logging.info(f"Created temporary directories for sample {sampleID}: data={self.dataDir[sampleID].name}, bed={self.bedDir[sampleID].name}")

                # Initialize storage for this sample if it doesn't exist
                if sampleID not in app.storage.general[self.mainuuid]:
                    app.storage.general[self.mainuuid][sampleID] = {}

                # Initialize counters for this sample if they don't exist

                for analysis in analyses:
                    if analysis.lower() not in self.exclude:
                        if analysis not in app.storage.general[self.mainuuid][sampleID]:
                            app.storage.general[self.mainuuid][sampleID][analysis] = {}
                            app.storage.general[self.mainuuid][sampleID][analysis][
                                "counters"
                            ] = Counter(
                                bam_count=0, bam_processed=0, bams_in_processing=0
                            )
                        elif (
                            "counters"
                            not in app.storage.general[self.mainuuid][sampleID][
                                analysis
                            ]
                        ):
                            app.storage.general[self.mainuuid][sampleID][analysis][
                                "counters"
                            ] = Counter(
                                bam_count=0, bam_processed=0, bams_in_processing=0
                            )

            for analysis in analyses:
                if analysis.lower() not in self.exclude:
                    app.storage.general[self.mainuuid][sampleID][analysis]["counters"][
                        "bam_count"
                    ] += 1
            
            # Add file to the list for this sample
            files_by_sample[sampleID].append(file)
            latest_files[sampleID] = max(latest_files[sampleID], filetime or 0)

            # Process if we have enough files for any sample
            for sample_id in list(files_by_sample.keys()):
                # The length of bams we allow to be processed at once will influence memory usage.
                if (
                    len(files_by_sample[sample_id]) >= 1
                ):  # This can be changed at any time
                    files_to_process = len(files_by_sample[sample_id])
                    logging.info(
                        f"Processing batch of {files_to_process} files for sample {sample_id}"
                    )

                    for analysis in analyses:
                        if analysis.lower() not in self.exclude:
                            try:
                                # Ensure the counter structure exists
                                app.storage.general[self.mainuuid][sample_id][analysis][
                                    "counters"
                                ]["bams_in_processing"] += files_to_process
                                logging.debug(
                                    f"Updated {analysis} counter for {sample_id} by {files_to_process}"
                                )
                            except Exception as e:
                                logging.error(
                                    f"Failed to update counter for {analysis}: {str(e)}"
                                )
            
                    await self.process_sample_files(
                        sample_id, files_by_sample[sample_id], latest_files[sample_id]
                    )

                    # Clear processed files
                    files_by_sample[sample_id] = []
                    latest_files[sample_id] = 0
            
        """
        # Process remaining files for each sample
        for sample_id, files in files_by_sample.items():
            if files:  # Only process if there are files
                files_to_process = len(files)
                if files_to_process > 0:
                    logging.info(
                        f"Processing remaining {files_to_process} files for sample {sample_id}"
                    )

                    # Update counters for each non-excluded analysis
                    analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
                    for analysis in analyses:
                        if analysis.lower() not in self.exclude:
                            try:
                                # Ensure the counter structure exists
                                if (
                                    analysis
                                    not in app.storage.general[self.mainuuid][sample_id]
                                ):
                                    app.storage.general[self.mainuuid][sample_id][
                                        analysis
                                    ] = {
                                        "counters": {
                                            "bams_in_processing": 0,
                                            "bam_processed": 0,
                                            "bams_failed": 0,
                                        }
                                    }
                                elif (
                                    "counters"
                                    not in app.storage.general[self.mainuuid][
                                        sample_id
                                    ][analysis]
                                ):
                                    app.storage.general[self.mainuuid][sample_id][
                                        analysis
                                    ]["counters"] = {
                                        "bams_in_processing": 0,
                                        "bam_processed": 0,
                                        "bams_failed": 0,
                                    }

                                app.storage.general[self.mainuuid][sample_id][analysis][
                                    "counters"
                                ]["bams_in_processing"] += files_to_process
                                logging.debug(
                                    f"Updated {analysis} counter for {sample_id} by {files_to_process}"
                                )
                            except Exception as e:
                                logging.error(
                                    f"Failed to update counter for {analysis}: {str(e)}"
                                )

                    await self.process_sample_files(
                        sample_id, files, latest_files[sample_id]
                    )
        """
        if not state.shutdown_event:
            self.process_bigbadmerge_tracker.active = True
            state.set_process_state("Merge Bam Analysis", ProcessState.WAITING_FOR_DATA)
        else:
            state.set_process_state("Merge Bam Analysis", ProcessState.STOPPING)
            state.stop_process("Merge Bam Analysis")
            self.finished = True

    async def process_sample_files(self, sampleID, tomerge, latest_file):
        """Process sample files using modkit."""
        if state.shutdown_event:
            logging.info("Shutdown event detected, skipping file processing")
            state.stop_process("Merge Bam Analysis")
            return

        try:
            # Track the number of BAM files seen and merged
            num_bam_files_seen = len(tomerge)
            logging.info(
                f"Processing {num_bam_files_seen} BAM files for sample ID: {sampleID}"
            )

            # Ensure counters exist for this sample
            analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
            for analysis in analyses:
                if analysis.lower() not in self.exclude:
                    if analysis not in app.storage.general[self.mainuuid][sampleID]:
                        app.storage.general[self.mainuuid][sampleID][analysis] = {
                            "counters": {
                                "bams_in_processing": 0,
                                "bam_processed": 0,
                                "bams_failed": 0,
                            }
                        }
                    elif (
                        "counters"
                        not in app.storage.general[self.mainuuid][sampleID][analysis]
                    ):
                        app.storage.general[self.mainuuid][sampleID][analysis][
                            "counters"
                        ] = {
                            "bams_in_processing": 0,
                            "bam_processed": 0,
                            "bams_failed": 0,
                        }

            # Write the length of the tomerge list to the output file FIRST
            tomerge_length_file = os.path.join(
                self.check_and_create_folder(self.output, sampleID),
                "tomerge_length.txt",
            )

            # Initialize the count
            if os.path.exists(tomerge_length_file):
                with self.safe_open(tomerge_length_file, "r") as f:
                    current_count = int(
                        f.readline().strip().split(": ")[1]
                    )  # Read the current count
            else:
                current_count = 0  # If the file doesn't exist, start from 0

            # Update the count
            new_count = current_count + len(tomerge)

            # Write the updated length of the tomerge list to the output file
            with self.safe_open(tomerge_length_file, "w") as f:
                f.write(f"Length of tomerge list: {new_count}\n")
            
            if len(tomerge) > 1:
                # Create temporary files with proper registration
                tempbam = tempfile.NamedTemporaryFile(
                    dir=self.check_and_create_folder(self.output, sampleID),
                    suffix=".bam",
                    delete=False
                )
                sorttempbam = tempfile.NamedTemporaryFile(
                    dir=self.check_and_create_folder(self.output, sampleID),
                    suffix=".bam",
                    delete=False
                )
                temp = tempfile.NamedTemporaryFile(
                    dir=self.check_and_create_folder(self.output, sampleID),
                    delete=False
                )
                
                # Register with file handle manager
                self.file_handle_manager.register_temp_file(tempbam)
                self.file_handle_manager.register_temp_file(sorttempbam)
                self.file_handle_manager.register_temp_file(temp)
                
                file = tempbam.name
                sortfile = sorttempbam.name
                temp_name = temp.name

                # Set process state to running
                state.set_process_state("Merge Bam Analysis", ProcessState.RUNNING)

                # Sort and merge BAM files
                try:  # Here we shouldn't need to use cpu_bound as we should be in the background.
                    await run.cpu_bound(
                        run_samtools_sort, file, tomerge, sortfile, self.threads
                    )
                except concurrent.futures.process.BrokenProcessPool:
                    logger.warning(
                        "Process pool was terminated. This is normal during shutdown."
                    )
                    state.set_process_state("Merge Bam Analysis", ProcessState.STOPPED)
                    return
            else:
                sortfile = tomerge[0]
                temp = tempfile.NamedTemporaryFile(
                    dir=self.check_and_create_folder(self.output, sampleID),
                    delete=False
                )
                
                # Register with file handle manager
                self.file_handle_manager.register_temp_file(temp)
                
                temp_name = temp.name
            

            with tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            ):
                try:  # Here we shouldn't need to use cpu_bound as we should be in the background.
                    print(f"Running modkit on {sortfile}")
                    # await run.cpu_bound(run_modkit, sortfile, temp_name, self.threads)
                    await run.cpu_bound(run_matkit, sortfile, temp_name)
                    # print("Modkit run complete")
                    # run_modkit(sortfile, temp_name, self.threads)
                except concurrent.futures.process.BrokenProcessPool:
                    logger.warning(
                        "Process pool was terminated. This is normal during shutdown."
                    )
                    state.set_process_state("Merge Bam Analysis", ProcessState.STOPPED)
                    return
                
                # Create output path specific to this sample
                parquet_path = os.path.join(
                    self.check_and_create_folder(self.output, sampleID),
                    f"{sampleID}.parquet",  # Use sampleID for the output filename
                )
                
                # Merge modkit files for this sample
                try:  # Here we shouldn't need to use cpu_bound as we should be in the background.
                    await run.cpu_bound(
                        merge_modkit_files,
                        [temp_name],
                        parquet_path,  # Use the parquet_path for the existing file
                        parquet_path,
                        self.cpgs_master_file,  # Use the parquet_path for the output file
                        sampleID,
                        self.output,
                        dict(self.mnpflex_config),
                        num_bam_files_seen,  # Pass the number of BAM files that contributed
                    )
                    
                    for analysis in analyses:
                        if analysis.lower() not in self.exclude:
                            if analysis == "STURGEON":
                                self.parquetqueuesturgeon.put(
                                    (parquet_path, sampleID, num_bam_files_seen)
                                )
                            elif analysis == "NANODX":
                                self.parquetqueuenanodx.put(
                                    (parquet_path, sampleID, num_bam_files_seen)
                                )
                            elif analysis == "PANNANODX":
                                self.parquetqueuepannanodx.put(
                                    (parquet_path, sampleID, num_bam_files_seen)
                                )
                            elif analysis == "FOREST":
                                self.parquetqueuecns.put(
                                    (parquet_path, sampleID, num_bam_files_seen)
                                )
                except concurrent.futures.process.BrokenProcessPool:
                    logger.warning(
                        "Process pool was terminated. This is normal during shutdown."
                    )
                    state.set_process_state("Merge Bam Analysis", ProcessState.STOPPED)
                    return
                
                # Log the number of BAM files processed
                logging.info(
                    f"Merged {num_bam_files_seen} BAM files into {parquet_path} for sample ID: {sampleID}"
                )
                #self.mergecounter += len(tomerge)
                
                # Update the processed counters for each analysis type
                analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
                for analysis in analyses:
                    if analysis.lower() not in self.exclude:
                        try:
                            if self.mainuuid not in app.storage.general:
                                app.storage.general[self.mainuuid] = {}
                            if sampleID not in app.storage.general[self.mainuuid]:
                                app.storage.general[self.mainuuid][sampleID] = {}
                            if (
                                analysis
                                not in app.storage.general[self.mainuuid][sampleID]
                            ):
                                app.storage.general[self.mainuuid][sampleID][
                                    analysis
                                ] = {
                                    "counters": {
                                        "bams_in_processing": 0,
                                        "bam_processed": 0,
                                        "bams_failed": 0,
                                    }
                                }

                            counters = app.storage.general[self.mainuuid][sampleID][
                                analysis
                            ]["counters"]

                            # Decrease the in_processing counter and increase the processed counter
                            counters["bams_in_processing"] = max(
                                0, counters["bams_in_processing"] - num_bam_files_seen
                            )

                            logging.debug(
                                f"Updated {analysis} processed counter for {sampleID} by {num_bam_files_seen}"
                            )
                        except Exception as e:
                            logging.error(
                                f"Error updating counters for {analysis}: {str(e)}"
                            )
                        if state.shutdown_event:
                            state.set_process_state(
                                "Merge Bam Analysis", ProcessState.STOPPED
                            )
                            return
                
                # Clean up temporary directories for this sample after processing is complete
                logging.info(f"Cleaning up temporary directories for sample {sampleID}")
                self.cleanup_temp_directories_for_sample(sampleID)
                
                # Clean up temporary files using file handle manager
                logging.info(f"Cleaning up temporary files for sample {sampleID}")
                self.cleanup_temp_files_for_sample(sampleID)
                
                # Also clean up the file paths manually for safety
                try:
                    if len(tomerge) > 1:
                        # Clean up temporary BAM files
                        if os.path.exists(file):
                            os.unlink(file)
                            logging.debug(f"Cleaned up temporary BAM file: {file}")
                        if os.path.exists(sortfile):
                            os.unlink(sortfile)
                            logging.debug(f"Cleaned up temporary sorted BAM file: {sortfile}")
                    if os.path.exists(temp_name):
                        os.unlink(temp_name)
                        logging.debug(f"Cleaned up temporary modkit file: {temp_name}")
                except Exception as e:
                    logging.error(f"Error cleaning up temporary files for sample {sampleID}: {str(e)}", exc_info=True)
            

        except concurrent.futures.process.BrokenProcessPool:
            logger.warning(
                "Process pool was terminated. This is normal during shutdown."
            )
            state.set_process_state("Merge Bam Analysis", ProcessState.STOPPED)
        except Exception as e:
            logging.error(f"Error in process_sample_files: {str(e)}", exc_info=True)
            state.set_process_state("Merge Bam Analysis", ProcessState.STOPPED)
            # Update the failed counters for each analysis type
            analyses = ["STURGEON", "NANODX", "PANNANODX", "FOREST"]
            for analysis in analyses:
                if analysis.lower() not in self.exclude:
                    try:
                        if (
                            self.mainuuid in app.storage.general
                            and sampleID in app.storage.general[self.mainuuid]
                            and analysis in app.storage.general[self.mainuuid][sampleID]
                        ):
                            app.storage.general[self.mainuuid][sampleID][analysis][
                                "counters"
                            ]["bams_in_processing"] -= num_bam_files_seen
                            app.storage.general[self.mainuuid][sampleID][analysis][
                                "counters"
                            ]["bams_failed"] += num_bam_files_seen
                            logging.debug(
                                f"Updated {analysis} failed counter for {sampleID} by {num_bam_files_seen}"
                            )
                    except Exception as nested_e:
                        logging.error(
                            f"Error updating failed counters for {analysis}: {str(nested_e)}"
                        )
            raise

    async def process_bams(self) -> None:
        """
        Process BAM files in the watch folder.
        """
        logging.info("Process Bam Starting")
        state.start_process("Serial Bam Analysis", ProcessType.BACKGROUND)
        state.set_process_state("Serial Bam Analysis", ProcessState.WAITING_FOR_DATA)
        self.process_bams_tracker.active = False
        if state.shutdown_event:
            print("shutdown_event is True, stopping process_bams")
            state.set_process_state("Serial Bam Analysis", ProcessState.STOPPING)
            self.process_bams_tracker.active = False
            state.stop_process("Serial Bam Analysis")
            return

        if "file" in app.storage.general[self.mainuuid]["bam_count"]:
            state.set_process_state("Serial Bam Analysis", ProcessState.RUNNING)
            while self.watchdogbamqueue.qsize() > 0:
                filename, timestamp = self.watchdogbamqueue.get()
                logging.info(f"Processing new BAM file from watchdog queue: {filename}")
                # First check the BAM file to get its sample ID
                baminfo, bamdata = await run.cpu_bound(check_bam, filename)
                sample_id = (
                    baminfo["sample_id"]
                    if not self.force_sampleid
                    else self.force_sampleid
                )

                # Validate sample_id - if None, use a fallback or skip
                if sample_id is None:
                    if self.force_sampleid:
                        sample_id = self.force_sampleid
                        logging.warning(
                            f"BAM file {filename} has no sample_id, using force_sampleid: {sample_id}"
                        )
                    else:
                        # Generate a fallback sample ID based on filename
                        sample_id = f"unknown_sample_{os.path.basename(filename)}"
                        logging.warning(
                            f"BAM file {filename} has no sample_id, using fallback: {sample_id}"
                        )

                logging.info(f"BAM file {filename} belongs to sample: {sample_id}")

                # Initialize sample-specific counters if they don't exist
                if sample_id not in app.storage.general[self.mainuuid]["samples"]:
                    logging.info(f"Initializing new sample storage for: {sample_id}")
                    app.storage.general[self.mainuuid]["samples"][sample_id] = {}
                    self.configure_storage(sample_id)
                    app.storage.general[self.mainuuid]["sample_list"].append(sample_id)
                    # Initialize sample-specific BAM tracking
                    if (
                        "bam_tracking"
                        not in app.storage.general[self.mainuuid]["samples"][sample_id]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "bam_tracking"
                        ] = {"counter": 0, "total_files": 0, "files": {}}

                # Update run information from BAM file
                if (
                    "run_info"
                    not in app.storage.general[self.mainuuid]["samples"][sample_id]
                ):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ] = {}

                # Update both the run_info dictionary and the arrays
                if baminfo.get("time_of_run"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["run_time"] = baminfo["time_of_run"]
                    if (
                        baminfo["time_of_run"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "run_time"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "run_time"
                        ].append(baminfo["time_of_run"])

                if baminfo.get("device_position"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["device"] = baminfo["device_position"]
                    if (
                        baminfo["device_position"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "devices"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "devices"
                        ].append(baminfo["device_position"])

                if baminfo.get("basecall_model"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["model"] = baminfo["basecall_model"]
                    if (
                        baminfo["basecall_model"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "basecall_models"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "basecall_models"
                        ].append(baminfo["basecall_model"])

                if baminfo.get("flow_cell_id"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["flow_cell"] = baminfo["flow_cell_id"]
                    if (
                        baminfo["flow_cell_id"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "flowcell_ids"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "flowcell_ids"
                        ].append(baminfo["flow_cell_id"])

                # Update sample-specific counters
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "bam_tracking"
                ]["counter"] += 1
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "bam_tracking"
                ]["total_files"] += 1
                app.storage.general[self.mainuuid]["samples"][sample_id][
                    "bam_tracking"
                ]["files"][filename] = timestamp

                # Also update global counters for overall tracking
                app.storage.general[self.mainuuid]["bam_count"]["counter"] += 1
                app.storage.general[self.mainuuid]["bam_count"]["total_files"] += 1
                app.storage.general[self.mainuuid]["bam_count"]["file"][
                    filename
                ] = timestamp

                # Inside the process_bams method, after updating pass_* or fail_* counters:
                # Update the combined counters
                self.update_counter(sample_id, "mapped_bases", bamdata["mapped_bases"])
                self.update_counter(
                    sample_id, "unmapped_bases", bamdata["unmapped_bases"]
                )
                self.update_counter(
                    sample_id, "mapped_reads_num", bamdata["mapped_reads_num"]
                )
                self.update_counter(
                    sample_id, "unmapped_reads_num", bamdata["unmapped_reads_num"]
                )
                if state.shutdown_event:
                    print("shutdown_event is True, stopping process_bams")
                    self.process_bams_tracker.active = False
                    state.stop_process("Serial Bam Analysis")
                    return

            # Process the files
            while len(app.storage.general[self.mainuuid]["bam_count"]["file"]) > 0:
                self.nofiles = False
                file = (
                    k := next(
                        iter(app.storage.general[self.mainuuid]["bam_count"]["file"])
                    ),
                    app.storage.general[self.mainuuid]["bam_count"]["file"].pop(k),
                )

                logging.info(f"Processing BAM file from global queue: {file[0]}")
                baminfo, bamdata = await run.cpu_bound(check_bam, file[0])
                sample_id = (
                    baminfo["sample_id"]
                    if not self.force_sampleid
                    else self.force_sampleid
                )

                # Validate sample_id - if None, use a fallback or skip
                if sample_id is None:
                    if self.force_sampleid:
                        sample_id = self.force_sampleid
                        logging.warning(
                            f"BAM file {file[0]} has no sample_id, using force_sampleid: {sample_id}"
                        )
                    else:
                        # Generate a fallback sample ID based on filename
                        sample_id = f"unknown_sample_{os.path.basename(file[0])}"
                        logging.warning(
                            f"BAM file {file[0]} has no sample_id, using fallback: {sample_id}"
                        )

                logging.info(f"BAM file {file[0]} belongs to sample: {sample_id}")

                # Remove from sample-specific tracking as well
                if sample_id in app.storage.general[self.mainuuid]["samples"]:
                    if (
                        "bam_tracking"
                        in app.storage.general[self.mainuuid]["samples"][sample_id]
                    ):
                        if (
                            file[0]
                            in app.storage.general[self.mainuuid]["samples"][sample_id][
                                "bam_tracking"
                            ]["files"]
                        ):
                            app.storage.general[self.mainuuid]["samples"][sample_id][
                                "bam_tracking"
                            ]["files"].pop(file[0])

                # Initialize sample storage if it doesn't exist
                if sample_id not in app.storage.general[self.mainuuid]["samples"]:
                    logging.info(f"Initializing new sample storage for: {sample_id}")
                    app.storage.general[self.mainuuid]["samples"][sample_id] = {}
                    self.configure_storage(sample_id)
                    app.storage.general[self.mainuuid]["sample_list"].append(sample_id)
                    # Initialize sample-specific BAM tracking
                    if (
                        "bam_tracking"
                        not in app.storage.general[self.mainuuid]["samples"][sample_id]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "bam_tracking"
                        ] = {"counter": 0, "total_files": 0, "files": {}}

                # Update run information from BAM file
                if (
                    "run_info"
                    not in app.storage.general[self.mainuuid]["samples"][sample_id]
                ):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ] = {}

                # Update both the run_info dictionary and the arrays
                if baminfo.get("time_of_run"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["run_time"] = baminfo["time_of_run"]
                    if (
                        baminfo["time_of_run"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "run_time"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "run_time"
                        ].append(baminfo["time_of_run"])

                if baminfo.get("device_position"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["device"] = baminfo["device_position"]
                    if (
                        baminfo["device_position"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "devices"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "devices"
                        ].append(baminfo["device_position"])

                if baminfo.get("basecall_model"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["model"] = baminfo["basecall_model"]
                    if (
                        baminfo["basecall_model"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "basecall_models"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "basecall_models"
                        ].append(baminfo["basecall_model"])

                if baminfo.get("flow_cell_id"):
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "run_info"
                    ]["flow_cell"] = baminfo["flow_cell_id"]
                    if (
                        baminfo["flow_cell_id"]
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "flowcell_ids"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "flowcell_ids"
                        ].append(baminfo["flow_cell_id"])

                # Process state and update counters
                if baminfo["state"] == "pass":
                    # logging.info(f"BAM file {filename} passed quality checks")
                    self.update_counter(sample_id, "bam_passed", 1)
                    self.update_counter(
                        sample_id, "pass_mapped_count", bamdata["mapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "pass_unmapped_count", bamdata["unmapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "pass_bases_count", bamdata["yield_tracking"]
                    )
                    self.update_counter(
                        sample_id, "pass_mapped_reads_num", bamdata["mapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id,
                        "pass_unmapped_reads_num",
                        bamdata["unmapped_reads_num"],
                    )
                    self.update_counter(
                        sample_id, "pass_mapped_bases", bamdata["mapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "pass_unmapped_bases", bamdata["unmapped_bases"]
                    )

                    # Update the combined counters for pass files
                    self.update_counter(
                        sample_id, "mapped_bases", bamdata["mapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "unmapped_bases", bamdata["unmapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "mapped_reads_num", bamdata["mapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id, "unmapped_reads_num", bamdata["unmapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id, "mapped_count", bamdata["mapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "unmapped_count", bamdata["unmapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "bases_count", bamdata["yield_tracking"]
                    )
                else:
                    # logging.info(f"BAM file {filename} failed quality checks")
                    self.update_counter(sample_id, "bam_failed", 1)
                    self.update_counter(
                        sample_id, "fail_mapped_count", bamdata["mapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "fail_unmapped_count", bamdata["unmapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "fail_bases_count", bamdata["yield_tracking"]
                    )
                    self.update_counter(
                        sample_id, "fail_mapped_reads_num", bamdata["mapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id,
                        "fail_unmapped_reads_num",
                        bamdata["unmapped_reads_num"],
                    )
                    self.update_counter(
                        sample_id, "fail_mapped_bases", bamdata["mapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "fail_unmapped_bases", bamdata["unmapped_bases"]
                    )

                    # Update the combined counters for fail files
                    self.update_counter(
                        sample_id, "mapped_bases", bamdata["mapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "unmapped_bases", bamdata["unmapped_bases"]
                    )
                    self.update_counter(
                        sample_id, "mapped_reads_num", bamdata["mapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id, "unmapped_reads_num", bamdata["unmapped_reads_num"]
                    )
                    self.update_counter(
                        sample_id, "mapped_count", bamdata["mapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "unmapped_count", bamdata["unmapped_reads"]
                    )
                    self.update_counter(
                        sample_id, "bases_count", bamdata["yield_tracking"]
                    )

                # Save master data to CSV
                try:
                    # Create a flattened dictionary for the DataFrame
                    sample_data = app.storage.general[self.mainuuid]["samples"][
                        sample_id
                    ]
                    flattened_data = {}

                    # Handle counters
                    if "file_counters" in sample_data:
                        for key, value in sample_data["file_counters"].items():
                            flattened_data[f"counter_{key}"] = value

                    # Handle arrays - convert to comma-separated strings
                    array_fields = [
                        "devices",
                        "basecall_models",
                        "run_time",
                        "flowcell_ids",
                    ]
                    for field in array_fields:
                        if field in sample_data:
                            flattened_data[field] = ",".join(
                                map(str, sample_data[field])
                            )

                    # Handle run_info dictionary
                    if "run_info" in sample_data:
                        for key, value in sample_data["run_info"].items():
                            flattened_data[f"run_info_{key}"] = value

                    # Handle bam_tracking
                    if "bam_tracking" in sample_data:
                        for key, value in sample_data["bam_tracking"].items():
                            if key == "files":
                                # Convert files dictionary to comma-separated string
                                flattened_data["bam_files"] = ",".join(
                                    sample_data["bam_tracking"]["files"].keys()
                                )
                            else:
                                flattened_data[f"bam_tracking_{key}"] = value

                    # Create DataFrame from flattened data
                    mydf = pd.DataFrame([flattened_data])

                    # Save to CSV
                    master_csv_path = os.path.join(
                        self.check_and_create_folder(self.output, sample_id),
                        "master.csv",
                    )
                    mydf.to_csv(master_csv_path, index=False)
                    logging.info(f"Updated master.csv for sample {sample_id}")
                except Exception as e:
                    logging.error(
                        f"Error writing master.csv for sample {sample_id}: {str(e)}"
                    )

                # Route to analysis queues
                analyses = ["forest", "sturgeon", "nanodx", "pannanodx"]
                if any(analysis.lower() not in self.exclude for analysis in analyses):
                    logging.info(
                        f"Routing {file[0]} to bigbadmerge queue for batch processing"
                    )
                    self.bamforbigbadmerge.put([file[0], file[1], sample_id])
                """
                # Individual analysis queue routing
                if "forest" not in self.exclude:
                    logging.info(f"Routing {file[0]} to forest analysis queue")
                    self.bamforcns.put([file[0], file[1], sample_id])
                if "sturgeon" not in self.exclude:
                    logging.info(f"Routing {file[0]} to sturgeon analysis queue")
                    self.bamforsturgeon.put([file[0], file[1], sample_id])
                if "nanodx" not in self.exclude:
                    logging.info(f"Routing {file[0]} to nanodx analysis queue")
                    self.bamfornanodx.put([file[0], file[1], sample_id])
                if "pannanodx" not in self.exclude:
                    logging.info(f"Routing {file[0]} to pannanodx analysis queue")
                    self.bamforpannanodx.put([file[0], file[1], sample_id])
                """
                if "cnv" not in self.exclude:
                    logging.info(f"Routing {file[0]} to cnv analysis queue")
                    self.bamforcnv.put([file[0], file[1], sample_id])
                if "coverage" not in self.exclude:
                    logging.info(f"Routing {file[0]} to coverage analysis queue")
                    self.bamfortargetcoverage.put([file[0], file[1], sample_id])
                if "mgmt" not in self.exclude:
                    logging.info(f"Routing {file[0]} to mgmt analysis queue")
                    self.bamformgmt.put([file[0], file[1], sample_id])
                if "fusion" not in self.exclude:
                    logging.info(f"Routing {file[0]} to fusion analysis queue")
                    self.bamforfusions.put([file[0], file[1], sample_id])

                self.nofiles = True

                if state.shutdown_event:
                    print("shutdown_event is True, stopping process_bams")
                    self.process_bams_tracker.active = False
                    state.stop_process("Serial Bam Analysis")
                    return

            logging.info("Process Bam Finishing")
            if not state.shutdown_event:
                self.process_bams_tracker.active = True
                state.set_process_state(
                    "Serial Bam Analysis", ProcessState.WAITING_FOR_DATA
                )
            else:
                self.finished = True
                state.set_process_state("Serial Bam Analysis", ProcessState.STOPPING)
                state.stop_process("Serial Bam Analysis")

    async def check_existing_bams(self, sequencing_summary=None):
        # ToDo: THis needs to be a background process.
        """
        Check and process existing BAM files based on a sequencing summary.

        :param sequencing_summary: Path to the sequencing summary file.
        :return: None
        """
        file_endings = {".bam"}

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
            if state.shutdown_event:
                self.terminate = True
                return
            app.storage.general[self.mainuuid]["bam_count"]["counter"] += 1
            app.storage.general[self.mainuuid]["bam_count"][
                "total_files"
            ] += 1  # Increment total files counter
            if "file" not in app.storage.general[self.mainuuid]["bam_count"]:
                app.storage.general[self.mainuuid]["bam_count"]["file"] = {}
            app.storage.general[self.mainuuid]["bam_count"]["file"][
                f
            ] = timestamp.timestamp()
            if self.simtime:
                await asyncio.sleep(0.1)
            else:
                await asyncio.sleep(0.1)

    def cleanup_temp_directories_for_sample(self, sample_id):
        """
        Clean up temporary directories for a specific sample.
        
        :param sample_id: The sample ID to clean up temporary directories for.
        """
        try:
            if sample_id in self.dataDir:
                logging.info(f"Cleaning up data temporary directory for sample {sample_id}")
                temp_dir = self.dataDir[sample_id]
                temp_dir.cleanup()
                del self.dataDir[sample_id]
                logging.info(f"Cleaned up data temporary directory for sample {sample_id}")
            
            if sample_id in self.bedDir:
                logging.info(f"Cleaning up bed temporary directory for sample {sample_id}")
                temp_dir = self.bedDir[sample_id]
                temp_dir.cleanup()
                del self.bedDir[sample_id]
                logging.info(f"Cleaned up bed temporary directory for sample {sample_id}")
                
        except Exception as e:
            logging.error(f"Error cleaning up temporary directories for sample {sample_id}: {str(e)}", exc_info=True)
    
    def cleanup_temp_files_for_sample(self, sample_id):
        """
        Clean up temporary files for a specific sample.
        
        :param sample_id: The sample ID to clean up temporary files for.
        """
        try:
            # This will be handled by the file handle manager
            logging.info(f"Cleaning up temporary files for sample {sample_id}")
            # The file handle manager will handle the cleanup automatically
        except Exception as e:
            logging.error(f"Error cleaning up temporary files for sample {sample_id}: {str(e)}", exc_info=True)

    def cleanup_all_temp_directories(self):
        """
        Clean up all temporary directories for all samples.
        """
        try:
            logging.info("Cleaning up all temporary directories")
            
            # Clean up data directories
            for sample_id in list(self.dataDir.keys()):
                try:
                    logging.info(f"Cleaning up data temporary directory for sample {sample_id}")
                    temp_dir = self.dataDir[sample_id]
                    temp_dir.cleanup()
                    del self.dataDir[sample_id]
                except Exception as e:
                    logging.error(f"Error cleaning up data temporary directory for sample {sample_id}: {str(e)}", exc_info=True)
            
            # Clean up bed directories
            for sample_id in list(self.bedDir.keys()):
                try:
                    logging.info(f"Cleaning up bed temporary directory for sample {sample_id}")
                    temp_dir = self.bedDir[sample_id]
                    temp_dir.cleanup()
                    del self.bedDir[sample_id]
                except Exception as e:
                    logging.error(f"Error cleaning up bed temporary directory for sample {sample_id}: {str(e)}", exc_info=True)
            
            logging.info("Completed cleanup of all temporary directories")
            
        except Exception as e:
            logging.error(f"Error in cleanup_all_temp_directories: {str(e)}", exc_info=True)

    def cleanup_old_temp_directories(self, active_samples=None):
        """
        Clean up temporary directories for samples that are no longer active.
        
        :param active_samples: List of sample IDs that are currently active. If None, assumes all samples in storage are active.
        """
        try:
            if active_samples is None:
                # If no active samples provided, use all samples in storage
                if self.mainuuid in app.storage.general:
                    active_samples = list(app.storage.general[self.mainuuid]["samples"].keys())
                else:
                    active_samples = []
            
            logging.info(f"Cleaning up old temporary directories. Active samples: {active_samples}")
            
            # Clean up data directories for inactive samples
            for sample_id in list(self.dataDir.keys()):
                if sample_id not in active_samples:
                    try:
                        logging.info(f"Cleaning up old data temporary directory for inactive sample {sample_id}")
                        temp_dir = self.dataDir[sample_id]
                        temp_dir.cleanup()
                        del self.dataDir[sample_id]
                    except Exception as e:
                        logging.error(f"Error cleaning up old data temporary directory for sample {sample_id}: {str(e)}", exc_info=True)
            
            # Clean up bed directories for inactive samples
            for sample_id in list(self.bedDir.keys()):
                if sample_id not in active_samples:
                    try:
                        logging.info(f"Cleaning up old bed temporary directory for inactive sample {sample_id}")
                        temp_dir = self.bedDir[sample_id]
                        temp_dir.cleanup()
                        del self.bedDir[sample_id]
                    except Exception as e:
                        logging.error(f"Error cleaning up old bed temporary directory for sample {sample_id}: {str(e)}", exc_info=True)
            
            logging.info("Completed cleanup of old temporary directories")
            
        except Exception as e:
            logging.error(f"Error in cleanup_old_temp_directories: {str(e)}", exc_info=True)

    def get_temp_directory_info(self):
        """
        Get information about current temporary directory usage for debugging.
        
        :return: Dictionary with information about temporary directories
        """
        info = {
            'data_directories': {},
            'bed_directories': {},
            'total_data_dirs': len(self.dataDir),
            'total_bed_dirs': len(self.bedDir)
        }
        
        for sample_id, temp_dir in self.dataDir.items():
            info['data_directories'][sample_id] = {
                'path': temp_dir.name,
                'exists': os.path.exists(temp_dir.name)
            }
        
        for sample_id, temp_dir in self.bedDir.items():
            info['bed_directories'][sample_id] = {
                'path': temp_dir.name,
                'exists': os.path.exists(temp_dir.name)
            }
        
        return info
    
    def get_file_handle_stats(self):
        """
        Get statistics about file handle usage for debugging.
        
        :return: Dictionary with file handle statistics
        """
        return self.file_handle_manager.get_stats()
    
    def log_file_handle_stats(self):
        """
        Log current file handle statistics for debugging.
        """
        stats = self.get_file_handle_stats()
        logging.info(f"File handle statistics: {stats}")
        
        temp_dir_info = self.get_temp_directory_info()
        logging.info(f"Temporary directory statistics: {temp_dir_info}")
    
    def safe_open(self, file_path, mode='r', **kwargs):
        """
        Safely open a file with automatic registration in the file handle manager.
        
        :param file_path: Path to the file to open
        :param mode: File open mode
        :param kwargs: Additional arguments to pass to open()
        :return: File handle registered with the file handle manager
        """
        file_handle = open(file_path, mode, **kwargs)
        return self.file_handle_manager.register_file_handle(file_handle)
    
    def safe_temp_file(self, **kwargs):
        """
        Safely create a temporary file with automatic registration.
        
        :param kwargs: Arguments to pass to tempfile.NamedTemporaryFile
        :return: Temporary file registered with the file handle manager
        """
        temp_file = tempfile.NamedTemporaryFile(**kwargs)
        return self.file_handle_manager.register_temp_file(temp_file)
    
    def safe_temp_directory(self, **kwargs):
        """
        Safely create a temporary directory with automatic registration.
        
        :param kwargs: Arguments to pass to tempfile.TemporaryDirectory
        :return: Temporary directory registered with the file handle manager
        """
        temp_dir = tempfile.TemporaryDirectory(**kwargs)
        return self.file_handle_manager.register_temp_directory(temp_dir)
    
    def force_cleanup_all_resources(self):
        """
        Force cleanup of all file handles and temporary resources.
        This should be called when memory issues are detected.
        """
        logging.warning("Forcing cleanup of all file handles and temporary resources")
        try:
            # Clean up all temporary directories
            self.cleanup_all_temp_directories()
            
            # Clean up all file handles and temporary files
            self.file_handle_manager.cleanup_all()
            
            # Log statistics after cleanup
            self.log_file_handle_stats()
            
            logging.info("Forced cleanup completed")
        except Exception as e:
            logging.error(f"Error during forced cleanup: {str(e)}", exc_info=True)
