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

import numpy as np
import math


class MinknowHistograms:
    """Class that displays a histogram from MinKNOW."""

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
        with ui.card().classes("flat border-[2px] no-shadow"):
            with ui.card().tight().classes("flat border-[2px] no-shadow"):
                ui.label("Histograms").classes("text-h6")
                with ui.row().props("align-middle"):
                    ui.avatar("bar_chart", color="grey", square=False)
                    ui.label(
                        "Information on read length distributions taken from MinNKOW."
                    ).classes(f"text-h6 text-overline {self.color}")
            with ui.column().classes("w-full").bind_visibility_from(self, "Live_Run"):
                self.myplot = StackedBarPlot()
                with ui.row():
                    ui.label().bind_text_from(
                        self,
                        "padding",
                        backward=lambda n: f"Suggested padding: {n} bases",
                    )
                    self.checkbox = ui.checkbox("Basecalling Enabled").bind_value_from(
                        self, "basecalling_enabled"
                    )
                    ui.label().bind_text_from(
                        self,
                        "events_to_base_ratio",
                        backward=lambda n: f"Events to Base Ratio: {n}",
                    )
                    ui.label().bind_text_from(
                        self, "sample_rate", backward=lambda n: f"Sample Rate: {n}"
                    )
                    ui.label().bind_text_from(
                        self, "purpose", backward=lambda n: f"Purpose: {n}"
                    )
                    ui.label().bind_text_from(
                        self, "Run_ID", backward=lambda n: f"Run ID: {n}"
                    )
                    ui.label().bind_text_from(
                        self, "Live_Run", backward=lambda n: f"Live Run: {n}"
                    )
                    ui.label().bind_text_from(
                        self,
                        "adaptive_sampling",
                        backward=lambda n: f"Adaptive Sampling: {n}",
                    )
                    ui.label().bind_text_from(
                        self,
                        "adaptive_sampling_ignore_chunk",
                        backward=lambda n: f"Adaptive Sampling Ignore Chunk: {n}",
                    )
            with ui.column().bind_visibility_from(
                self, "Live_Run", backward=lambda v: not v
            ):
                ui.label("Set up a run using MinKNOW to see histogram information.")

    def connect_me(self) -> None:
        self.connection = self.Position.connect()

    def stream_current_run(self) -> None:
        # worker = get_current_worker()
        while True:
            self.stream = self.connection.acquisition.watch_current_acquisition_run()
            for info in self.stream:
                try:
                    if info.state == manager_pb2.FlowCellPosition.State.STATE_RUNNING:
                        self.Run_ID = info.run_id
                        # print("Run ID: " + str(self.Run_ID))
                        if not self.Live_Run:
                            self.Live_Run = True
                    else:
                        self.Run_ID = None
                        self.Live_Run = False
                except KeyboardInterrupt:
                    self.stream.cancel()
                    print("KeyboardInterrupt")
                    break

    def stream_histogram_info(self, step_size=None, end=None) -> None:
        """Stream histogram info."""
        # worker = get_current_worker()
        while True:
            while self.Live_Run:
                # print("Trying to Streaming Histogram Info")
                if self.Run_ID:
                    self.purpose = (
                        self.connection.protocol.get_protocol_purpose().purpose
                    )
                    # print(f"Purpose: {self.purpose}")
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

                        # print(
                        #    f"Streaming info with step_size={step_size}, end={end}"
                        # )

                        for info in self._stream_histogram_info:
                            # print(str(info))
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
                            # print("Sending Histogram Info")
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
    def __init__(self):

        self.echart = ui.echart(
            {
                "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                #'title': {'text': 'Read Length Histogram'},
                "yAxis": {"type": "value", "name": "Yield"},
                "xAxis": {"type": "category", "name": "Bin Size", "data": []},
                "legend": {"textStyle": {"color": "gray"}},
                "series": [
                    # {'type': 'bar', 'name': 'Alpha', 'stack': 'x', 'data': [0.1, 0.2]},
                    # {'type': 'bar', 'name': 'Beta', 'stack': 'x',  'data': [0.3, 0.4]},
                ],
            }
        ).classes("w-full")

    def add_series(self, rundata: dict) -> None:
        # TODO: Maintain data visibility state between updates.
        # print(dir(self.echart.props))
        self.echart.options["xAxis"]["data"] = rundata["bucket_values"]
        # print(self.echart.options['series'])
        for i, value in enumerate(rundata["read_end_reason"]):
            dataexists = False
            for seri in self.echart.options["series"]:
                if seri["name"] == value:
                    dataexists = True
                    # print (f"Found {value}")
                    # print (seri)
                    seri["data"] = list(rundata["histdata"][i])
            if not dataexists:
                # print(f"Not Found {value} in {[x['name'] for x in self.echart.options['series']]}")
                self.echart.options["series"].append(
                    {
                        "type": "bar",
                        "name": f"{value}",
                        "stack": "x",
                        "data": list(rundata["histdata"][i]),
                    }
                )
        # print(self.echart.options['legend'])
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
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    index_page()
    ui.run(
        port=port, reload=reload, title="Readfish NiceGUI"
    )  # , native=True, fullscreen=False, window_size=(1200, 1000))


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
