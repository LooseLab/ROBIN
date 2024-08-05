"""
Provides a class that handles events from the file system watcher.
"""

from watchdog.events import FileSystemEventHandler, FileSystemEvent
import time
from collections import defaultdict
from typing import Dict, DefaultDict, Any, Tuple
from queue import Queue
import logging

# Create a logger for this module
logger = logging.getLogger(__name__)


def configure_logging(level=logging.INFO, log_file=None):
    """
    Configure the logging for this module.

    Args:
        level (int): The logging level to set. Defaults to logging.INFO.
        log_file (str): The path to the log file. If None, no file handler will be added.
    """
    logger.setLevel(level)

    # Create a console handler if no handlers are set
    if not logger.handlers:
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)


class BamEventHandler(FileSystemEventHandler):
    """
    A class to handle file system events specifically for ".bam" files.

    This class extends the FileSystemEventHandler from the watchdog library
    to handle creation events for ".bam" files within a watch folder.

    Attributes:
        bam_queue (Queue[Tuple[str, float]]): A queue to store the file paths and creation timestamps of ".bam" files.
    """

    def __init__(self, bam_queue: Queue[Tuple[str, float]]):
        """
        Initializes the BamEventHandler with a bam_queue.

        Args:
            bam_queue (Queue[Tuple[str, float]]): A queue to store the file paths and creation timestamps of ".bam" files.
        """
        super().__init__()
        self.bam_queue = bam_queue

    def on_created(self, event: FileSystemEvent) -> None:
        """
        Handles the file creation event.

        This method is called when a new file is created in the watch folder.
        If the created file has a ".bam" extension, it puts the file path and creation timestamp into the bam_queue
        and logs an info message.

        Args:
            event (FileSystemEvent): The event object representing the file system event.
        """
        if not event.is_directory and event.src_path.endswith(".bam"):
            try:
                self.bam_queue.put((event.src_path, time.time()))
                logger.info(f"New .bam file detected: {event.src_path}")
            except Exception as e:
                logger.error(f"Error processing new .bam file: {e}")


if __name__ == "__main__":
    # Configure logging
    configure_logging(level=logging.DEBUG)

    # Example of how to initialize and use BamEventHandler
    bam_queue = Queue()
    event_handler = BamEventHandler(bam_queue=bam_queue)

    logger.debug("BamEventHandler initialized and ready to process events.")