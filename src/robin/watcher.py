"""File watching functionality using Watchdog."""

import os
import subprocess
import time
from pathlib import Path
from typing import List, Optional

from tqdm import tqdm
from watchdog.events import FileSystemEvent, FileSystemEventHandler
from watchdog.observers import Observer


class FileChangeHandler(FileSystemEventHandler):
    """Handler for file system events."""

    def __init__(
        self,
        patterns: List[str],
        ignore_patterns: List[str],
        command: Optional[str],
        verbose: bool = False,
    ) -> None:
        """Initialize the file change handler.

        Args:
            patterns: List of file patterns to watch
            ignore_patterns: List of file patterns to ignore
            command: Command to run when files change
            verbose: Enable verbose output
        """
        self.patterns = patterns
        self.ignore_patterns = ignore_patterns
        self.command = command
        self.verbose = verbose
        self.last_modified = 0
        self.debounce_time = 0.5  # seconds

    def _should_process_event(self, event: FileSystemEvent) -> bool:
        """Check if the event should be processed.

        Args:
            event: The file system event

        Returns:
            True if the event should be processed, False otherwise
        """
        if event.is_directory:
            return False

        # Check if file matches any ignore patterns
        for pattern in self.ignore_patterns:
            if Path(event.src_path).match(pattern):
                return False

        # Check if file matches any watch patterns
        for pattern in self.patterns:
            if Path(event.src_path).match(pattern):
                return True

        return False

    def _debounce(self) -> bool:
        """Implement debouncing to prevent multiple rapid events.

        Returns:
            True if the event should be processed, False if it should be debounced
        """
        current_time = time.time()
        if current_time - self.last_modified < self.debounce_time:
            return False
        self.last_modified = current_time
        return True

    def on_created(self, event: FileSystemEvent) -> None:
        """Handle file creation events."""
        if self._should_process_event(event) and self._debounce():
            self._handle_event("created", event.src_path)

    def on_modified(self, event: FileSystemEvent) -> None:
        """Handle file modification events."""
        if self._should_process_event(event) and self._debounce():
            self._handle_event("modified", event.src_path)

    def on_deleted(self, event: FileSystemEvent) -> None:
        """Handle file deletion events."""
        if self._should_process_event(event) and self._debounce():
            self._handle_event("deleted", event.src_path)

    def on_moved(self, event: FileSystemEvent) -> None:
        """Handle file move/rename events."""
        if self._should_process_event(event) and self._debounce():
            self._handle_event("moved", event.src_path)

    def _handle_event(self, event_type: str, file_path: str) -> None:
        """Handle a file system event.

        Args:
            event_type: Type of event (created, modified, deleted, moved)
            file_path: Path to the affected file
        """
        if self.verbose:
            print(f"[{event_type.upper()}] {file_path}")

        if self.command:
            self._run_command(file_path)

    def _run_command(self, file_path: str) -> None:
        """Run the specified command when a file changes.

        Args:
            file_path: Path to the changed file
        """
        try:
            # Replace {file} placeholder with the actual file path
            cmd = self.command.replace("{file}", file_path)

            if self.verbose:
                print(f"Running command: {cmd}")

            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=30,
            )

            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(f"Error: {result.stderr}", file=os.sys.stderr)

        except subprocess.TimeoutExpired:
            print(f"Command timed out: {self.command}", file=os.sys.stderr)
        except Exception as e:
            print(f"Error running command: {e}", file=os.sys.stderr)


class FileWatcher:
    """File watcher using Watchdog."""

    def __init__(
        self,
        path: Path,
        recursive: bool = False,
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        command: Optional[str] = None,
        verbose: bool = False,
        show_progress: bool = True,
    ) -> None:
        """Initialize the file watcher.

        Args:
            path: Path to watch
            recursive: Whether to watch directories recursively
            patterns: List of file patterns to watch
            ignore_patterns: List of file patterns to ignore
            command: Command to run when files change
            verbose: Enable verbose output
            show_progress: Show progress bars for file processing
        """
        self.path = path
        self.recursive = recursive
        self.patterns = patterns or ["*"]
        self.ignore_patterns = ignore_patterns or []
        self.command = command
        self.verbose = verbose
        self.show_progress = show_progress

        self.observer = Observer()
        self.handler = FileChangeHandler(
            patterns=self.patterns,
            ignore_patterns=self.ignore_patterns,
            command=self.command,
            verbose=self.verbose,
        )

    def process_existing_files(self) -> None:
        """Process existing files that match the patterns."""
        if not self.command:
            return

        if self.verbose:
            print(f"Processing existing files in: {self.path}")
            print(f"Patterns: {self.patterns}")
            print(f"Ignore patterns: {self.ignore_patterns}")

        existing_files = []

        # Find all existing files that match patterns
        for pattern in self.patterns:
            if pattern == "*":
                # Default pattern - all files
                if self.recursive:
                    pattern_files = list(self.path.rglob("*"))
                else:
                    pattern_files = list(self.path.glob("*"))
            elif "**" in pattern:
                # Recursive pattern like **/*.bam
                pattern_files = list(self.path.rglob(pattern.replace("**/", "")))
            else:
                # Simple pattern like *.bam
                if self.recursive:
                    pattern_files = list(self.path.rglob(pattern))
                else:
                    pattern_files = list(self.path.glob(pattern))

            # Filter to only files (not directories) and apply ignore patterns
            for file_path in pattern_files:
                if file_path.is_file():
                    # Check if file should be ignored
                    should_ignore = False
                    for ignore_pattern in self.ignore_patterns:
                        if file_path.match(ignore_pattern):
                            should_ignore = True
                            break

                    if not should_ignore:
                        existing_files.append(file_path)

        # Remove duplicates and sort
        existing_files = sorted(set(existing_files))

        if existing_files:
            if self.verbose:
                print(f"Found {len(existing_files)} existing file(s) to process:")
                for file_path in existing_files:
                    print(f"  - {file_path}")
                print()

            # Process each existing file with progress bar
            with tqdm(
                total=len(existing_files),
                desc="Processing existing files",
                unit="file",
                disable=not self.show_progress,
            ) as pbar:
                for file_path in existing_files:
                    if self.show_progress:
                        pbar.set_postfix_str(f"Processing: {file_path.name}")
                    self.handler._run_command(str(file_path))
                    pbar.update(1)
        else:
            if self.verbose:
                print("No existing files found matching the patterns.")

    def start(self, process_existing: bool = True) -> None:
        """Start watching for file changes.

        Args:
            process_existing: Whether to process existing files before starting to watch
        """
        if self.verbose:
            print(f"Starting file watcher for: {self.path}")
            print(f"Patterns: {self.patterns}")
            print(f"Ignore patterns: {self.ignore_patterns}")
            if self.command:
                print(f"Command: {self.command}")

        # Process existing files first if requested
        if process_existing and self.command:
            self.process_existing_files()
            if self.verbose:
                print("Now watching for new changes...\n")

        self.observer.schedule(
            self.handler,
            str(self.path),
            recursive=self.recursive,
        )
        self.observer.start()

        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            self.stop()

    def stop(self) -> None:
        """Stop watching for file changes."""
        self.observer.stop()
        self.observer.join()
        if self.verbose:
            print("File watcher stopped.")
