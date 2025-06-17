from typing import Optional
from pathlib import Path
import logging
from nicegui import ui
import os

from .config import BrainMethConfig
from .metrics import MetricsManager
from .storage import StorageManager
from ..services.background import BackgroundProcessor
from ..ui.components import (
    MetricsDisplay,
    SampleList,
    InformationPanel,
    ReportGenerator,
)
from ..utils.file_utils import check_and_create_folder
from robin.utilities.performance_metrics import PerformanceMetrics


class BrainMeth:
    """Main class coordinating all components of the BrainMeth application."""

    def __init__(self, config: BrainMethConfig):
        """Initialize the BrainMeth application.

        Args:
            config: Application configuration
        """
        self.config = config
        self.logger = logging.getLogger(__name__)

        # Initialize managers
        self.metrics_manager = MetricsManager(config.mainuuid)
        self.storage_manager = StorageManager(config.mainuuid)

        # Initialize background processor
        self.background_processor = BackgroundProcessor(
            config=config,
            metrics_manager=self.metrics_manager,
            storage_manager=self.storage_manager,
        )

        # Initialize UI components
        self.metrics_display = MetricsDisplay(self.metrics_manager)
        self.sample_list = SampleList(self.storage_manager)
        self.information_panel = InformationPanel(
            storage_manager=self.storage_manager, metrics_manager=self.metrics_manager
        )
        self.report_generator = ReportGenerator(
            storage_manager=self.storage_manager, metrics_manager=self.metrics_manager
        )

        # Initialize performance metrics
        self.performance_metrics = PerformanceMetrics(
            metrics_dir=os.path.join(
                self.config.get_app_data_dir(), "performance_metrics"
            )
        )

        # Start RAM tracking with 30-second intervals
        self.performance_metrics.start_ram_tracking(interval_seconds=30)

        # Set up UI update timer for performance metrics
        ui.timer(1.0, self._update_performance_metrics_ui)

    async def start(self):
        """Start the application."""
        try:
            # Start background processing
            await self.background_processor.start()

            # Render UI
            await self.render_ui()
        except Exception as e:
            self.logger.error(f"Error starting application: {e}")
            raise

    async def stop(self):
        """Stop the application."""
        try:
            await self.background_processor.stop()
        except Exception as e:
            self.logger.error(f"Error stopping application: {e}")

    async def render_ui(self, sample_id: Optional[str] = None):
        """Render the main UI.

        Args:
            sample_id: Optional sample ID to display
        """
        try:
            with ui.row():
                with ui.column():
                    self.metrics_display.render()
                    self.sample_list.render()

                with ui.column():
                    await self.information_panel.render(sample_id)
                    self.report_generator.render()
        except Exception as e:
            self.logger.error(f"Error rendering UI: {e}")
            raise

    async def add_watchfolder(self, watchfolder: str):
        """Add a watch folder for monitoring.

        Args:
            watchfolder: Path to the folder to watch
        """
        try:
            watchfolder_path = Path(watchfolder)
            if not watchfolder_path.exists():
                check_and_create_folder(watchfolder_path)

            self.config.watchfolder = watchfolder_path
            self.logger.info(f"Added watch folder: {watchfolder}")
        except Exception as e:
            self.logger.error(f"Error adding watch folder: {e}")
            raise

    async def pick_file(self) -> None:
        """Open file picker dialog."""
        try:
            # Implementation of file picking
            pass
        except Exception as e:
            self.logger.error(f"Error picking file: {e}")
            raise

    def replay(self):
        """Replay the current session."""
        try:
            # Implementation of replay functionality
            pass
        except Exception as e:
            self.logger.error(f"Error replaying session: {e}")
            raise

    @property
    def min_start_time(self):
        """Get the minimum start time of the current session."""
        try:
            # Implementation of min start time calculation
            return 0
        except Exception as e:
            self.logger.error(f"Error getting min start time: {e}")
            return 0

    def _update_performance_metrics_ui(self):
        """Update the performance metrics UI in the Process Status panel."""
        if (
            hasattr(self.config, "process_status")
            and self.config.process_status.performance_metrics_container
        ):
            # Clear existing content
            self.config.process_status.performance_metrics_container.clear()
            # Create new performance metrics UI in the correct container
            self.performance_metrics.create_ui(
                parent=self.config.process_status.performance_metrics_container
            )

    def shutdown(self):
        """Clean up resources before shutdown."""
        # Stop RAM tracking
        self.performance_metrics.stop_ram_tracking()

        # Save final metrics
        self.performance_metrics.save_metrics()

        # Existing shutdown code
        self.background_processor.shutdown()
