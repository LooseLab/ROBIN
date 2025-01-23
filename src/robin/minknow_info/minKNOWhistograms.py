"""
MinKNOW Histograms Module

This module provides real-time visualization of MinKNOW sequencing data through
interactive histograms. It handles the collection, processing, and display of
read length distributions and other sequencing metrics.

Key Features:
- Real-time histogram updates
- Read length distribution visualization
- Adaptive sampling monitoring
- Interactive data exploration
- Multiple data series display

The module includes classes for:
- StatusChip: Visual status indicator
- MinknowHistograms: Main histogram management
- StackedBarPlot: Interactive plot generation

Dependencies:
- nicegui: Web interface components
- minknow_api: MinKNOW data access
- numpy: Numerical computations
- logging: Application logging
"""

# Python imports.
from __future__ import annotations
from nicegui import ui

import threading

import time


from robin import theme

from minknow_api.manager import Manager
import minknow_api.manager_pb2 as manager_pb2
from minknow_api.statistics_pb2 import (
    ReadLengthHistogramSplit,
    ReadEndReason,
    DataSelection,
)
import logging
import numpy as np
import math

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


class StatusChip(ui.chip):
    """
    A custom UI chip component for displaying status information.

    This component extends the basic UI chip with color-coding based on
    the status text (Enabled/Disabled).

    Attributes:
        Inherits all attributes from ui.chip
    """

    def _handle_text_change(self, text: str) -> None:
        """
        Handle changes to the chip's text and update its appearance.

        Args:
            text (str): The new text to display

        Note:
            Colors the chip green for "Enabled" and red for other states
        """
        super()._handle_text_change(text)
        if "Enabled" in text:
            self.props("color=green-500")
        else:
            self.props("color=red-500")


class MinknowHistograms:
    """
    Class for managing and displaying MinKNOW sequencing histograms.

    This class handles the collection and visualization of read length
    distributions and other sequencing metrics in real-time.

    Attributes:
        Position: MinKNOW position object
        Live_Run (bool): Whether a run is currently active
        Run_ID: Current run identifier
        color (str): UI color theme
        _live_data (list): Collected data points
        padding (int): Data padding value
        basecalling_enabled (bool): Basecalling status
        adaptive_sampling (bool): Adaptive sampling status
        adaptive_sampling_ignore_chunk (int): Chunk size to ignore
    """

    # basecalling_enabled: reactive[bool] = reactive(False)
    # purpose: reactive[str] = reactive("")
    # events_to_base_ratio: reactive[float] = reactive(0.0)
    # sample_rate: reactive[int] = reactive(0)
    # run_data: reactive[dict] = reactive(dict())
    # filter: reactive[str] = reactive("All")
    # adaptive_sampling: reactive[bool] = reactive(False)
    # adaptive_sampling_ignore_chunk: reactive[int] = reactive(0)
    # padding: reactive[int] = reactive(0)

    def __init__(self, position):
        """
        Initialize a MinknowHistograms instance.

        Args:
            position: MinKNOW position object to monitor
        """
        super().__init__()
        self.Position = position
        self.Live_Run = True
        self.Run_ID = None
        self.color = "text-blue"
        self.connect_me()
        self._live_data = []
        self.padding = 0
        self.basecalling_enabled = False
        self.adaptive_sampling = False
        self.adaptive_sampling_ignore_chunk = 0
        self._stream_histogram_info = None
        self.stream_current_run_thread = threading.Thread(
            target=self.stream_current_run, args=()
        )
        self.stream_current_run_thread.daemon = True
        self.stream_current_run_thread.start()
        self.stream_histogram_info_thread = threading.Thread(
            target=self.stream_histogram_info, args=()
        )
        self.stream_histogram_info_thread.daemon = True
        self.stream_histogram_info_thread.start()
        self.renderme()

    def renderme(self):
        # Main container with full width
        with ui.card().classes("w-full").bind_visibility_from(self, "Live_Run"):
            with ui.card_section():
                with ui.row().classes("items-center"):
                    ui.icon("bar_chart", size="lg").classes("text-primary q-mr-md")
                    with ui.column():
                        ui.label("Read Length Distribution").classes("text-h6")
                        ui.label(
                            "Information on read length distributions taken from MinKNOW"
                        ).classes(f"text-caption {self.color}")

            # Main content area with full width
            with ui.card_section().classes("w-full").bind_visibility_from(
                self, "Live_Run"
            ):
                # Plot container with full width and fixed height
                with ui.card().classes("w-full").style("min-height: 500px"):
                    self.myplot = StackedBarPlot()

            # Information grid
            with ui.grid(columns=4):
                # Run Status
                with ui.card().classes("q-pa-sm"):
                    ui.label("Run Status").classes("text-caption text-secondary")
                    ui.label().bind_text_from(
                        self, "Run_ID", backward=lambda n: f"Run ID: {n}"
                    ).classes("text-body2")
                    ui.label().bind_text_from(
                        self, "Live_Run", backward=lambda n: f"Live Run: {n}"
                    ).classes("text-body2")

                # Technical Details
                with ui.card():  # .classes('q-pa-sm'):
                    ui.label("Technical Details").classes("text-caption text-secondary")
                    ui.label().bind_text_from(
                        self, "sample_rate", backward=lambda n: f"Sample Rate: {n}"
                    ).classes("text-body2")
                    ui.label().bind_text_from(
                        self,
                        "events_to_base_ratio",
                        backward=lambda n: f"Events to Base Ratio: {n}",
                    ).classes("text-body2")

                # Configuration
                with ui.card():  # .classes('q-pa-sm'):
                    ui.label("Configuration").classes("text-caption text-secondary")
                    StatusChip(text="").bind_text_from(
                        self,
                        "basecalling_enabled",
                        backward=lambda v: (
                            "Basecalling Enabled" if v else "Basecalling Disabled"
                        ),
                    )
                    ui.label().bind_text_from(
                        self, "purpose", backward=lambda n: f"Purpose: {n}"
                    ).classes("text-body2")

                # Adaptive Sampling
                with ui.card():  # .classes('q-pa-sm'):
                    ui.label("Adaptive Sampling").classes("text-caption text-secondary")
                    StatusChip(text="").bind_text_from(
                        self,
                        "adaptive_sampling",
                        backward=lambda v: (
                            "Adaptive Sampling Enabled"
                            if v
                            else "Adaptive Sampling Disabled"
                        ),
                    )
                    ui.label().bind_text_from(
                        self,
                        "adaptive_sampling_ignore_chunk",
                        backward=lambda n: f"Chunk: {n}",
                    ).classes("text-body2")

                # Padding info at the bottom
                ui.label().bind_text_from(
                    self, "padding", backward=lambda n: f"Suggested padding: {n} bases"
                ).classes("text-body1 q-mt-md")

            # Empty state message
            with ui.column().classes("q-pa-lg text-center").bind_visibility_from(
                self, "Live_Run", backward=lambda v: not v
            ):
                ui.label(
                    "Set up a run using MinKNOW to see histogram information."
                ).classes("text-secondary")

    def connect_me(self) -> None:
        """
        Establish connection to the MinKNOW position.
        """
        self.connection = self.Position.connect()

    def stream_current_run(self) -> None:
        """
        Monitor the current sequencing run.

        This method continuously monitors the run status and updates
        the Run_ID and Live_Run status accordingly.
        """
        while True:
            self.stream = self.connection.acquisition.watch_current_acquisition_run()
            for info in self.stream:
                try:
                    if info.state == manager_pb2.FlowCellPosition.State.STATE_RUNNING:
                        self.Run_ID = info.run_id
                        if not self.Live_Run:
                            self.Live_Run = True
                    else:
                        self.Run_ID = None
                        self.Live_Run = False
                except KeyboardInterrupt:
                    self.stream.cancel()
                    logger.info("KeyboardInterrupt")
                    break

    def stream_histogram_info(self, step_size=None, end=None) -> None:
        """
        Stream and process histogram information from MinKNOW.

        This method collects read length distribution data and updates
        the visualization in real-time.

        Args:
            step_size: Time step between data points
            end: End time for data collection
        """
        while True:
            while self.Live_Run:
                if self.Run_ID:
                    self.purpose = (
                        self.connection.protocol.get_protocol_purpose().purpose
                    )
                    if self.purpose == "sequencing_run":
                        self.basecalling_enabled = (
                            self.connection.acquisition.get_current_acquisition_run().config_summary.basecalling_enabled
                        )
                        rltype = 2 if self.basecalling_enabled else 1
                        self.events_to_base_ratio = round(
                            self.connection.acquisition.get_current_acquisition_run().config_summary.events_to_base_ratio,
                            2,
                        )
                        self.sample_rate = int(
                            self.connection.acquisition.get_current_acquisition_run().config_summary.sample_rate
                        )
                        self._stream_histogram_info = (
                            self.connection.statistics.stream_read_length_histogram(
                                acquisition_run_id=self.Run_ID,
                                bucket_value_type=1,
                                read_length_type=rltype,
                                poll_time_seconds=5,
                                data_selection=DataSelection(step=step_size, end=end),
                                split=ReadLengthHistogramSplit(read_end_reason=True),
                                discard_outlier_percent=0.05,
                            )
                        )

                        for info in self._stream_histogram_info:
                            _rundata = dict()
                            _rundata["read_end_reason"] = list()
                            _rundata["bucket_values"] = list()
                            _rundata["n50"] = list()
                            _rundata["histdata"] = list()

                            for hist in info.histogram_data:
                                if sum(hist.bucket_values) > 0:
                                    _rundata["read_end_reason"].append(
                                        f"{ReadEndReason.Name(hist.filtering[0].read_end_reason)}"
                                    )
                                    _rundata["n50"].append(hist.n50)
                                    _rundata["histdata"].append(hist.bucket_values)

                            for bucket in info.bucket_ranges:
                                _rundata["bucket_values"].append(bucket.end)

                            self.run_data = _rundata
                            self.myplot.add_series(self.run_data)
                            self.calculate_sequenced_read_lengths()
                        time.sleep(1)
                    time.sleep(2)
            time.sleep(5)

    def trim_trailing_zeros(self, lst: list) -> None:
        while lst and lst[-1] == 0:
            lst.pop()
        return lst

    def calculate_sequenced_read_lengths(self):
        """
        Using the histogram plots, calculate the approximate read lengths for each read end reason.
        In addition calculate a corrected read length for each read end reason taking into account the
        events to base ratio if basecalling not enabled.
        Finally calculate the mean read lengths accounting for adaptive sampling.
        :return:
        """
        if not self.adaptive_sampling:
            for i, item in enumerate(self.run_data["read_end_reason"]):
                if item.startswith("DataServiceUnblockMuxChange"):
                    self.adaptive_sampling = True
                    self.adaptive_sampling_index = (
                        i  # This is the index of the DataServiceUnblockMuxChange
                    )
        # When calculating the mean read length, we should ignore MuxChange reads as well as UnblockMuxChange reads
        # as these are not true reads and will skew the mean read length.
        # We should also ignore reads that are in the DataServiceMuxUnblock range.
        # This is because we will see proportionally more of these reads in a run with adaptive sampling enabled
        # which does not reflect the true underlying mean.
        if self.adaptive_sampling:
            # Get the adaptive sampling bin max value.
            # print(
            #    len(
            #        self.trim_trailing_zeros(
            #            self.run_data["histdata"][self.adaptive_sampling_index]
            #        )
            #    )
            # )
            self.adaptive_sampling_ignore_chunk = len(
                self.trim_trailing_zeros(
                    self.run_data["histdata"][self.adaptive_sampling_index]
                )
            )
        read_count = 0
        read_length = 0
        for i, end_reason in enumerate(self.run_data["read_end_reason"]):
            if end_reason.startswith("Signal"):
                read_count += np.sum(
                    np.array(
                        self.run_data["histdata"][i][
                            self.adaptive_sampling_ignore_chunk :
                        ]
                    )
                    / np.array(
                        self.run_data["bucket_values"][
                            self.adaptive_sampling_ignore_chunk :
                        ]
                    )
                )
                read_length += np.sum(
                    np.array(
                        self.run_data["histdata"][i][
                            self.adaptive_sampling_ignore_chunk :
                        ]
                    )
                )
        if read_count > 0:
            if not self.basecalling_enabled:
                self.padding = (
                    math.ceil(
                        read_length / read_count * self.events_to_base_ratio / 1000
                    )
                    * 1000
                )
            else:
                self.padding = math.ceil(read_length / read_count / 1000) * 1000
            # print(f"Suggested padding: {self.padding}")


class StackedBarPlot:
    """
    Class for creating and managing stacked bar plots.

    This class handles the creation and updating of interactive
    stacked bar plots using the echarts library.

    Attributes:
        echart: The echarts plot object
    """

    def __init__(self):
        """Initialize a stacked bar plot with proper sizing and layout"""
        # Keep track of selected state
        self.selected_state = {}
        self.series_names = set()

        self.echart = (
            ui.echart(
                {
                    "grid": {
                        "left": "10%",  # Increased left padding for y-axis labels
                        "right": "10%",  # Increased right padding for second y-axis
                        "bottom": "15%",  # Increased bottom padding for x-axis labels
                        "top": "25%",  # Increased top padding for title and legend
                        "containLabel": True,
                    },
                    "toolbox": {
                        "show": True,
                        "feature": {"saveAsImage": {"title": "Save Image"}},
                    },
                    "tooltip": {"trigger": "axis", "axisPointer": {"type": "cross"}},
                    "yAxis": {
                        "type": "value",
                        "name": "Read Count",
                        "nameLocation": "middle",
                        "nameGap": 50,
                        "axisLabel": {"formatter": "{value}", "margin": 10},
                    },
                    "xAxis": {
                        "type": "category",
                        "name": "Read Length",
                        "nameLocation": "middle",
                        "nameGap": 40,
                        "data": [],
                        "axisLabel": {"rotate": 45, "margin": 15},
                    },
                    "legend": {
                        "type": "scroll",
                        "top": "50px",  # Move legend below title
                        "padding": [10, 10],  # Add padding around legend
                        "selected": self.selected_state,  # Link to our state tracker
                        "selectedMode": "multiple",  # Allow multiple selection
                    },
                    "series": [],
                    "animation": False,
                }
            )
            .classes("w-full")
            .style("min-height: 500px")
        )

    def add_series(self, rundata: dict) -> None:
        """Add or update data series in the plot with proper filtering"""
        if not rundata or "bucket_values" not in rundata:
            return

        # Update x-axis data
        self.echart.options["xAxis"]["data"] = rundata["bucket_values"]

        # Get current series by name for quick lookup
        current_series = {
            series["name"]: i for i, series in enumerate(self.echart.options["series"])
        }

        # Color mapping for consistency
        color_map = {
            "SignalPositive": "#91CC75",  # Green
            "SignalNegative": "#EE6666",  # Red
            "DataServiceUnblockMuxChange": "#5470C6",  # Blue
            "UnblockMuxChange": "#FAC858",  # Yellow
            "Other": "#73C0DE",  # Light blue
        }

        # Update or add series
        for i, value in enumerate(rundata["read_end_reason"]):
            # Skip series with all zero values
            if all(v == 0 for v in rundata["histdata"][i]):
                continue

            if value in current_series:
                # Just update the data of the existing series
                series_index = current_series[value]
                self.echart.options["series"][series_index]["data"] = list(
                    rundata["histdata"][i]
                )
            else:
                # Only add new series if it doesn't exist
                color = next((v for k, v in color_map.items() if k in value), "#73C0DE")
                self.echart.options["series"].append(
                    {
                        "type": "bar",
                        "name": value,
                        "stack": "total",
                        "data": list(rundata["histdata"][i]),
                        "itemStyle": {"color": color},
                        "emphasis": {"focus": "series"},
                        "animation": False,
                    }
                )

        # Update the chart
        self.echart.update()


def index_page() -> None:
    initial_ip = "127.0.0.1"
    my_connection = Manager(host=initial_ip)
    with theme.frame("Test Histograms"):
        # my_connection.connect_to_minknow()
        positions = list(my_connection.flow_cell_positions())
        ui.label(f"{positions[0]}")
        MinknowHistograms(positions[0])


def run_class(port: int, reload: bool):
    """
    Helper function to run the histogram visualization application.

    Args:
        port (int): Port number to serve the application
        reload (bool): Whether to enable auto-reload for development

    Returns:
        None
    """
    index_page()
    ui.run(port=port, reload=reload, title="Readfish NiceGUI")


def main():  # , threads, simtime, watchfolder, output, sequencing_summary):
    """
    Entrypoint for when GUI is launched directly.
    :return: None
    """
    run_class(port=12398, reload=False)


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    if __name__ == "__mp_main__":
        print("GUI launched by auto-reload")
    run_class(port=12398, reload=True)
