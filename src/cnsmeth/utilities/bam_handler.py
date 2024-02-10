"""
Provides a class that handles events from the file system watcher.
"""

from watchdog.events import FileSystemEventHandler
import time


class BamEventHandler(FileSystemEventHandler):
    def __init__(self, bam_count):
        """This class handles events from the file system watcher."""
        self.bam_count = bam_count

    def on_created(self, event):
        """This will add a file which is added to the watchfolder to the creates and the info file."""
        if event.src_path.endswith(".bam"):
            self.bam_count["counter"] += 1
            if "file" not in self.bam_count:
                self.bam_count["file"] = {}
            self.bam_count["file"][event.src_path] = time.time()
