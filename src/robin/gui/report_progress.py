"""
report_progress.py

Queue-based progress tracking for report generation using NiceGUI patterns.
"""

import logging
import queue
import threading
from datetime import datetime
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


def normalize_report_progress(progress: Optional[float]) -> Optional[float]:
    """Clamp report progress to 0–1 and round to whole percent for UI display."""
    if progress is None:
        return None
    percent = max(0, min(100, int(round(float(progress) * 100))))
    return percent / 100.0


class ReportProgressManager:
    """Manages report generation progress using a queue-based system."""

    def __init__(self):
        """Initialize the progress manager."""
        self.progress_queue = queue.Queue()
        self.active_reports: Dict[str, Dict[str, Any]] = {}
        self._notification_ids: Dict[str, str] = {}

    def start_report(self, sample_id: str) -> str:
        """Start tracking a new report generation.

        Args:
            sample_id: ID of the sample being processed

        Returns:
            Notification ID for tracking
        """
        try:
            # Store tracking info
            self.active_reports[sample_id] = {
                "notification_id": None,  # Will be set when notification is shown
                "stage": "initializing",
                "start_time": datetime.now(),
                "sample_id": sample_id,
                "progress": 0.0,
            }

            # Queue the initial notification
            self.progress_queue.put({"type": "start", "sample_id": sample_id})

            logger.info(f"Started tracking report generation for {sample_id}")
            return sample_id  # Return sample_id as identifier

        except Exception as e:
            logger.error(f"Error starting report tracking: {e}")
            return None

    def update_progress(
        self, sample_id: str, stage: str, message: str, progress: float = None
    ):
        """Update progress for a specific report.

        Args:
            sample_id: ID of the sample
            stage: Current stage (initializing, processing_sections, building_pdf, etc.)
            message: Human-readable message
            progress: Progress value (0.0 to 1.0)
        """
        if sample_id not in self.active_reports:
            logger.warning(f"No active report found for {sample_id}")
            return

        # Put update in queue for UI thread processing
        try:
            self.progress_queue.put(
                {
                    "type": "update",
                    "sample_id": sample_id,
                    "stage": stage,
                    "message": message,
                    "progress": normalize_report_progress(progress),
                }
            )
        except Exception as e:
            logger.error(f"Error queuing progress update: {e}")

    def complete_report(self, sample_id: str, filename: str = None):
        """Mark a report as completed.

        Args:
            sample_id: ID of the sample
            filename: Name of the generated file
        """
        if sample_id not in self.active_reports:
            logger.warning(f"No active report found for {sample_id}")
            return

        # Put completion in queue for UI thread processing
        try:
            self.progress_queue.put(
                {"type": "complete", "sample_id": sample_id, "filename": filename}
            )
        except Exception as e:
            logger.error(f"Error queuing completion: {e}")

    def error_report(self, sample_id: str, error_message: str):
        """Mark a report as failed.

        Args:
            sample_id: ID of the sample
            error_message: Error description
        """
        if sample_id not in self.active_reports:
            logger.warning(f"No active report found for {sample_id}")
            return

        # Put error in queue for UI thread processing
        try:
            self.progress_queue.put(
                {
                    "type": "error",
                    "sample_id": sample_id,
                    "error_message": error_message,
                }
            )
        except Exception as e:
            logger.error(f"Error queuing error: {e}")

    def process_queue(self):
        """Process queued progress updates on the UI thread."""
        # This method is now handled by the GUI launcher
        # to ensure proper UI context
        pass

    def _remove_report(self, sample_id: str):
        """Remove a report from active tracking."""
        if sample_id in self.active_reports:
            del self.active_reports[sample_id]
        if sample_id in self._notification_ids:
            del self._notification_ids[sample_id]
        logger.debug(f"Removed {sample_id} from active reports")

    def get_active_reports(self) -> Dict[str, Dict[str, Any]]:
        """Get all currently active reports."""
        return self.active_reports.copy()


# Global instance
progress_manager = ReportProgressManager()


def create_progress_callback(sample_id: str):
    """Create a progress callback function for a specific sample.

    Args:
        sample_id: ID of the sample being processed

    Returns:
        Callback function that can be passed to report generation
    """
    # Start tracking this report
    progress_manager.start_report(sample_id)

    def progress_callback(progress_data: Dict[str, Any]):
        """Progress callback function."""
        stage = progress_data.get("stage", "unknown")
        message = progress_data.get("message", "")
        progress = normalize_report_progress(progress_data.get("progress"))

        progress_manager.update_progress(
            sample_id=sample_id, stage=stage, message=message, progress=progress
        )

    return progress_callback
