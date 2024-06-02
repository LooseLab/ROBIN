"""
Provides a class that handles events from the file system watcher.
"""

from watchdog.events import FileSystemEventHandler, FileSystemEvent
import time
from collections import defaultdict
from typing import Dict, DefaultDict
import logging

# Configure logging
logger = logging.getLogger(__name__)


class BamEventHandler(FileSystemEventHandler):
    """
    A class to handle file system events specifically for ".bam" files.

    This class extends the FileSystemEventHandler from the watchdog library
    to handle creation events for ".bam" files within a watch folder.

    Attributes:
        bam_count (DefaultDict[str, Any]): A dictionary to keep track of the number of ".bam" files
                                           and their creation timestamps.
    """

    def __init__(self, bam_count: DefaultDict[str, Dict[str, float]]):
        """
        Initializes the BamEventHandler with a bam_count dictionary.

        Args:
            bam_count (DefaultDict[str, Dict[str, float]]): A dictionary to store the count and details of ".bam" files.
                                                            Expected to have a "counter" key initialized to zero.
        """
        self.bam_count = bam_count

    def on_created(self, event: FileSystemEvent):
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
            self.bam_count["counter"] += 1
            self.bam_count["file"][event.src_path] = time.time()
            logger.debug(
                f"New .bam file detected: {event.src_path}, total count: {self.bam_count['counter']}"
            )


if __name__ == "__main__":
    # Example of how to initialize and use BamEventHandler
    bam_count = defaultdict(lambda: {"counter": 0, "file": {}})
    event_handler = BamEventHandler(bam_count=bam_count)
