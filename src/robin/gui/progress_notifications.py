"""
progress_notifications.py

This module handles progress notifications for report generation in the NiceGUI interface.
"""

import logging
from datetime import datetime
from typing import Any, Callable, Dict, Optional

logger = logging.getLogger(__name__)


class ReportProgressNotifier:
    """Handles progress notifications for report generation."""

    def __init__(self):
        """Initialize the progress notifier."""
        self.active_reports: Dict[str, Dict[str, Any]] = {}
        self.notification_ids: Dict[str, str] = {}

    def handle_progress_update(self, progress_data: Dict[str, Any]):
        """Handle a progress update event.

        Args:
            progress_data: Dictionary containing progress information
        """
        sample_id = progress_data.get("sample_id")
        stage = progress_data.get("stage")
        message = progress_data.get("message")
        progress = progress_data.get("progress", 0.0)

        if not sample_id:
            logger.warning("Progress update missing sample_id")
            return

        # Update active reports tracking
        self.active_reports[sample_id] = {
            "stage": stage,
            "message": message,
            "progress": progress,
            "timestamp": datetime.now(),
            "completed_sections": progress_data.get("completed_sections", 0),
            "total_sections": progress_data.get("total_sections", 0),
            "current_section": progress_data.get("current_section", ""),
        }

        # Send update to GUI using existing system
        self._send_gui_update(sample_id, progress_data)

        # Show notification
        self._update_notification(sample_id, progress_data)

        # Clean up completed/error reports
        if stage in ["completed", "error"]:
            if sample_id in self.active_reports:
                del self.active_reports[sample_id]

    def _send_gui_update(self, sample_id: str, progress_data: Dict[str, Any]):
        """Send progress update to GUI using existing update system."""
        try:
            from robin.gui_launcher import UpdateType, send_gui_update

            # Send progress update to GUI
            send_gui_update(
                UpdateType.PROGRESS_UPDATE,
                {
                    "sample_id": sample_id,
                    "report_progress": progress_data,
                    "active_reports": self.active_reports.copy(),
                },
                priority=1,  # High priority for progress updates
            )
        except Exception as e:
            logger.error(f"Error sending GUI update: {e}")

    def _update_notification(self, sample_id: str, progress_data: Dict[str, Any]):
        """Update or create a progress notification."""
        stage = progress_data.get("stage")
        message = progress_data.get("message")
        progress = progress_data.get("progress", 0.0)

        # Determine notification type and content
        if stage == "completed":
            self._show_completion_notification(sample_id, progress_data)
        elif stage == "error":
            self._show_error_notification(sample_id, progress_data)
        else:
            self._show_progress_notification(sample_id, progress_data)

    def _show_progress_notification(
        self, sample_id: str, progress_data: Dict[str, Any]
    ):
        """Show a progress notification."""
        stage = progress_data.get("stage")
        message = progress_data.get("message")
        progress = progress_data.get("progress", 0.0)
        completed_sections = progress_data.get("completed_sections", 0)
        total_sections = progress_data.get("total_sections", 0)
        current_section = progress_data.get("current_section", "")

        # Create detailed message
        if total_sections > 0:
            section_info = f" ({completed_sections}/{total_sections} sections)"
            if current_section:
                section_info += f" - {current_section}"
        else:
            section_info = ""

        detailed_message = f"{message}{section_info}"

        # Calculate progress percentage
        progress_percent = int(progress * 100)

        # Create notification with progress
        notification_type = "info"
        if stage in ["initializing", "loading_data"]:
            notification_type = "info"
        elif stage == "processing_sections":
            notification_type = "ongoing"
        elif stage in ["building_pdf", "exporting_csv", "creating_zip"]:
            notification_type = "info"

        # Show notification with progress information
        try:
            from nicegui import ui

            ui.notify(
                f"[{sample_id}] {detailed_message} ({progress_percent}%)",
                type=notification_type,
                timeout=(
                    0 if stage != "completed" else 5000
                ),  # Persistent until completion
                position="top-right",
            )
        except Exception as e:
            logger.error(f"Error showing notification: {e}")

    def _show_completion_notification(
        self, sample_id: str, progress_data: Dict[str, Any]
    ):
        """Show a completion notification."""
        filename = progress_data.get("filename", "report.pdf")

        try:
            from nicegui import ui

            ui.notify(
                f"[{sample_id}] Report generation completed: {filename}",
                type="positive",
                timeout=5000,
                position="top-right",
            )
        except Exception as e:
            logger.error(f"Error showing completion notification: {e}")

    def _show_error_notification(self, sample_id: str, progress_data: Dict[str, Any]):
        """Show an error notification."""
        error_message = progress_data.get("error_message", "Unknown error")
        error_details = progress_data.get("error_details", "")

        full_message = f"[{sample_id}] Report generation failed: {error_message}"
        if error_details:
            full_message += f" ({error_details})"

        try:
            from nicegui import ui

            ui.notify(
                full_message, type="negative", timeout=10000, position="top-right"
            )
        except Exception as e:
            logger.error(f"Error showing error notification: {e}")

    def get_active_reports(self) -> Dict[str, Dict[str, Any]]:
        """Get information about currently active report generations."""
        return self.active_reports.copy()

    def is_report_active(self, sample_id: str) -> bool:
        """Check if a report is currently being generated."""
        return sample_id in self.active_reports


# Global instance
progress_notifier = ReportProgressNotifier()


def create_progress_callback(sample_id: str):
    """Create a progress callback function for a specific sample.

    Args:
        sample_id: ID of the sample being processed

    Returns:
        Callback function that can be passed to report generation
    """

    def progress_callback(progress_data: Dict[str, Any]):
        """Progress callback function."""
        # Ensure sample_id is set
        progress_data["sample_id"] = sample_id
        progress_notifier.handle_progress_update(progress_data)

    return progress_callback
