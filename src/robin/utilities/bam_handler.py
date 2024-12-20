"""
Provides a class that handles events from the file system watcher.
"""

from watchdog.events import FileSystemEventHandler, FileSystemEvent
import time
from collections import defaultdict
from typing import Dict, DefaultDict, Any
from queue import Queue
import logging

# Create a logger for this module
logger = logging.getLogger(__name__)


def configure_logging(level=logging.INFO):
    """
    Configure the logging for this module.

    Args:
        level (int): The logging level to set. Defaults to logging.INFO.
    """
    logger.setLevel(level)

    # Create a console handler if no handlers are set
    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)


class BamEventHandler(FileSystemEventHandler):
    """
    A class to handle file system events specifically for ".bam" files.

    This class extends the FileSystemEventHandler from the watchdog library
    to handle creation and move events for ".bam" files within a watch folder.

    Attributes:
        bam_count (Queue): A queue to store BAM file paths and their detection timestamps.
    """

    def __init__(self, bam_count: Queue):
        """
        Initializes the BamEventHandler with a bam_count queue.

        Args:
            bam_count (Queue): A queue to store BAM file paths and their timestamps.
        """
        super().__init__()
        self.bam_count = bam_count

    def _handle_bam_file(self, file_path: str) -> None:
        """
        Helper method to process BAM files consistently across different events.

        This method checks if the file has a .bam extension and if so, adds it to
        the processing queue with its current timestamp.

        Args:
            file_path (str): Path to the file to process
        """
        if file_path.endswith(".bam"):
            self.bam_count.put((file_path, time.time()))
            logger.info(f"New .bam file detected: {file_path}")

    def on_created(self, event: FileSystemEvent) -> None:
        """
        Handles the file creation event.

        This method is called when a new file is created in the watch folder.
        If the created file has a ".bam" extension, it processes the file.

        Args:
            event (FileSystemEvent): The event object representing the file system event.
        """
        if not event.is_directory:
            self._handle_bam_file(event.src_path)

    def on_moved(self, event: FileSystemEvent) -> None:
        """
        Handles file move/rename events.

        This method is called when a file is moved or renamed in the watch folder,
        which catches rsync's final move operation when copying files into the
        watched directory.

        Args:
            event (FileSystemEvent): The event object representing the move/rename event,
                                   containing both source and destination paths.
        """
        if not event.is_directory:
            self._handle_bam_file(event.dest_path)


def create_bam_count() -> DefaultDict[str, Dict[str, Any]]:
    """
    Creates and returns a defaultdict for tracking BAM file counts and details.

    Returns:
        DefaultDict[str, Dict[str, Any]]: A dictionary with 'counter' initialized to 0 and an empty 'file' dict.
    """
    return defaultdict(lambda: {"counter": 0, "file": {}})


if __name__ == "__main__":
    # Configure logging
    configure_logging(level=logging.DEBUG)

    # Example of how to initialize and use BamEventHandler
    bam_count = create_bam_count()
    event_handler = BamEventHandler(bam_count=bam_count)

    logger.debug("BamEventHandler initialized and ready to process events.")
