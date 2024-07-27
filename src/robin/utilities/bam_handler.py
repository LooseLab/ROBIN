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
    to handle creation events for ".bam" files within a watch folder.

    Attributes:
        bam_count (DefaultDict[str, Any]): A dictionary to keep track of the number of ".bam" files
                                           and their creation timestamps.
    """

    def __init__(self, bam_count: Queue):
        """
        Initializes the BamEventHandler with a bam_count dictionary.

        Args:
            bam_count (DefaultDict[str, Dict[str, Any]]): A dictionary to store the count and details of ".bam" files.
                                                          Expected to have a "counter" key initialized to zero.
        """
        super().__init__()
        self.bam_count = bam_count

    def on_created(self, event: FileSystemEvent) -> None:
        """
        Handles the file creation event.

        This method is called when a new file is created in the watch folder.
        If the created file has a ".bam" extension, it increments the counter in bam_count
        and records the file's creation timestamp.

        Args:
            event (FileSystemEvent): The event object representing the file system event.
        """
        if event.is_directory:
            return

        if event.src_path.endswith(".bam"):
            self.bam_count.put((event.src_path,time.time()))
            logger.info(
                #f"New .bam file detected: {event.src_path}, total count: {self.bam_count['counter']}"
                 f"New .bam file detected: {event.src_path}" #, total count: {self.bam_count['counter']}"
            )
            #except Exception as e:
            #    print("Bam error in handler.")
            #    print(e)


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
