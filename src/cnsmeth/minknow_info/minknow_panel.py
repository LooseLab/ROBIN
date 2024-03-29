import queue
from nicegui import ui, run
from typing import BinaryIO
import pandas as pd
import time
import asyncio
import threading

from datetime import datetime

from contextlib import contextmanager
from nicegui.events import ValueChangeEventArguments
from io import StringIO
import os
import ipaddress
import grpc
# MinKNOW API Imports
from minknow_api.manager import Manager
import minknow_api.manager_pb2 as manager_pb2
from minknow_api.protocol_pb2 import ProtocolPhase
from readfish._utils import get_device

@contextmanager
def disable(button: ui.button):
    button.disable()
    try:
        yield
    finally:
        button.enable()


class Position():
    def __init__(self, position, ip):
        super().__init__()
        self.ip = ip
        self.position = get_device(position.name, host=self.ip)
        self.connection = self.position.connect()
        self.Run_ID = None
        self.setup()
        self.check_flowcell = threading.Thread(target=self.watch_flowcell, args=())
        self.check_flowcell.daemon = True
        self.check_flowcell.start()
        self.check_run = threading.Thread(target=self.stream_current_run, args=())
        self.check_run.daemon = True
        self.check_run.start()
        self.check_instance = threading.Thread(target=self.stream_instance_activity, args=())
        self.check_instance.daemon = True
        self.check_instance.start()


    def setup(self):
        with ui.splitter().classes("w-full") as splitter:
            with splitter.before:
                with ui.card().classes("w-full p-8"):
                    ui.label('Readfish Control.').classes('mr-2')
                    #self.readfishconfig = ReadfishArgs(self.position)
            with splitter.after:
                #minknow_info(self.position)
                #with ui.row():#.classes("p-4 border-dashed border-2 w-full"):
                    #self.my_object = MinknowHistograms(self.position)
                with ui.card().classes("w-full p-8"):
                    with ui.card().classes("drop-shadow"):
                        ui.label('MinKNOW Monitoring.').classes('ml-2')
                        with ui.row():
                            with ui.card().classes("drop-shadow"):
                                ui.label().bind_text_from(self, "Run_ID", backward=lambda n: f'Run ID: {n}')
                                ui.label().bind_text_from(self.position, "name", backward=lambda n: f'Name: {n}')
                                ui.label().bind_text_from(self.position, "device_type", backward=lambda n: f'Device Type: {n}')
                                ui.label().bind_text_from(self.position, "host", backward=lambda n: f'Host: {n}')
                                ui.label().bind_text_from(self.position, "running", backward=lambda n: f'Running: {n}')
                                ui.label().bind_text_from(self.position, "state", backward=lambda n: f'State: {n}')
                            self.protocol_run_info = ui.column().classes("p-4 border-dashed border-2")
                        self.flowcell_info = ui.row()

                        self.yield_summary_info = ui.row().classes("p-4 border-dashed border-2")

                        self.acquisition_run_info = ui.column()





    def watch_flowcell(self):
        while True:
            if self.connection:
                #with self.flowcell_info:
                for stuff in self.connection.device.stream_flow_cell_info():
                    self.flowcell_info.clear()
                    with self.flowcell_info:
                        with ui.card().classes("drop-shadow"):
                            for field in stuff.DESCRIPTOR.fields:
                                ui.label(f"{field.name}: {getattr(stuff, field.name)}")


                pass



            else:
                time.sleep(1)

    def stream_current_run(self) -> None:
        """
        This function will stream the current run info from the minknow api.
        It runs as a worker to not block the main thread.
        """
        first_run = True
        while True:
            self._stream_current_run = (
                self.connection.acquisition.watch_current_acquisition_run()
            )
            for info in self._stream_current_run:
                if info.state == manager_pb2.FlowCellPosition.State.STATE_RUNNING:
                    #print(info)
                    self.Run_ID = info.run_id
                    self.Live_Run = True
                    #print(f"Run {self.Run_ID} seen on {self.position}")
                    #self.not_running.styles.display = "none"
                    #self.running.styles.display = "block"
                    msg = f"Run {self.Run_ID} seen on {self.position}"
                    title = "MinKNOW Info"
                    severity = "information"
                    timeout = 5
                else:
                    self.Run_ID = None
                    self.Live_Run = False
                    #self.not_running.styles.display = "block"
                    #self.running.styles.display = "none"
                    if first_run:
                        msg = f"Connected to {self.position}."
                        title = "MinKNOW Info"
                        severity = "information"
                        timeout = 5
                    else:
                        msg = f"Run {self.Run_ID} finished on {self.position}."
                        title = "MinKNOW Info"
                        severity = "warning"
                        timeout = 10
                if first_run:
                    first_run = False
                #self.notify(msg, title=title, severity=severity, timeout=timeout)

    def stream_instance_activity(self) -> None:
        """
        This function will stream instance activity from the minknow api.
        We configure a connection and give it a handler.
        This allows us to call cancel on the handler on exit to escape the loop.
        """

        self._stream_instance_activity = (
            self.connection.instance.stream_instance_activity()
        )
        for info in self._stream_instance_activity:
            #self.log("Seen info")
            #self.query_one("Static#ElapsedTime").update(
            #    f"Elapsed Time: {datetime.now() - self.Start_Time}"
            #)
            if info.HasField("device_info"):
                self.Channel_Count = info.device_info.device_info.max_channel_count
            if info.HasField("acquisition_run_info"):
                self.acquisition_run_info.clear()
                with self.acquisition_run_info:
                    for field in info.acquisition_run_info.DESCRIPTOR.fields:
                        ui.label(f"{field.name}: {getattr(info.acquisition_run_info, field.name)}")
                ### This isn't giving us the start time for some reason?!
                self.start_time = datetime.fromtimestamp(
                    info.acquisition_run_info.start_time.seconds
                    + info.acquisition_run_info.start_time.nanos / 1e9
                )
            if info.HasField("basecall_speed"):
                #self.log(f"{info.basecall_speed.mean_basecall_speed}")
                #self.log(type(info.basecall_speed.mean_basecall_speed))
                self.Mean_Basecall_Speed = info.basecall_speed.mean_basecall_speed
            # if info.HasField("flow_cell_health"):
            # if info.HasField("flow_cell_info"):
            if info.HasField("n50"):
                self.N50 = info.n50.n50
                self.Estimated_N50 = info.n50.estimated_n50
            if info.HasField("protocol_run_info"):
                self.protocol_run_info.clear()
                with self.protocol_run_info:
                #    for field in info.protocol_run_info.DESCRIPTOR.fields:
                #        ui.label(f"{field.name}: {getattr(info.protocol_run_info, field.name)}")
                    self.Experiment_Group = info.protocol_run_info.user_info.protocol_group_id.value
                    ui.label(f"Experiment Group: {self.Experiment_Group}")
                    self.Sample_ID = info.protocol_run_info.user_info.sample_id.value
                    ui.label(f"Sample ID: {self.Sample_ID}")
                    self.Current_Output_Directory = info.protocol_run_info.output_path
                    ui.label(f"Current Output Directory: {self.Current_Output_Directory}")
                    self.Kit = info.protocol_run_info.meta_info.tags["kit"].string_value
                    ui.label(f"Kit: {self.Kit}")
                    self.Phase = ProtocolPhase.Name(info.protocol_run_info.phase)
                    ui.label(f"Phase: {self.Phase}")
                    self.Start_Time = datetime.fromtimestamp(
                        info.protocol_run_info.start_time.seconds
                        + info.protocol_run_info.start_time.nanos / 1e9
                    )
                    ui.label(f"Start Time: {self.Start_Time}")
            if info.HasField("yield_summary"):
                self.yield_summary_info.clear()
                with self.yield_summary_info:
                #    for field in info.yield_summary.DESCRIPTOR.fields:
                #        ui.label(f"{field.name}: {getattr(info.yield_summary, field.name)}")
                #self.log("yield_summary")
                    self.Read_Count = info.yield_summary.read_count
                    ui.label(f"Read Count: {self.Read_Count}")
                    self.Percent_Basecalled = info.yield_summary.fraction_basecalled
                    ui.label(f"Percent Basecalled: {self.Percent_Basecalled}")
                    self.Pass_Read_Count = info.yield_summary.basecalled_pass_read_count
                    ui.label(f"Pass Read Count: {self.Pass_Read_Count}")
                    self.Fail_Read_Count = info.yield_summary.basecalled_fail_read_count
                    ui.label(f"Fail Read Count: {self.Fail_Read_Count}")
                    self.Pass_Bases = info.yield_summary.basecalled_pass_bases
                    ui.label(f"Pass Bases: {self.Pass_Bases}")
                    self.Fail_Bases = info.yield_summary.basecalled_fail_bases
                    ui.label(f"Fail Bases: {self.Fail_Bases}")

class MinKNOWFish():
    """
    A way of handling a minknow connection.
    """

    def __init__(self, **kwargs):
        """
        This is called when the app is initialized.
        """
        super().__init__(**kwargs)
        self.tabs = None
        self.manager = None
        self.positions = None
        self.selected_position = None
        self.connection_ip: str | None = None
        self.connected=False
        #self.pos_monitor = threading.Thread(target=self.watch_positions)
        #self.pos_monitor.daemon = True
        #self.pos_monitor.start()



    def _check_ip(self, ip: str) -> None:
        """
        Called by check_connection to validate the IP address whilst it is being typed.
        :param ip:
        :return: bool
        """
        try:
            ipaddress.ip_address(ip)
            return True
        except ValueError:
            return False

    def _action_ip(self, e):
        if self.ipinput.validate():
            self.connect_now.visible = True
            self.result.set_text('Connect to: ' + e.value)
            self.connection_ip = e.value
        else:
            self.result.set_text('Invalid IP address: ' + e.value)
            self.connect_now.visible = False

    def check_connection(self):
        """
        This is called when the app is initialized and checks
        to see if we have an active connection to a minKNOW instance.
        To do this it requests a URL for the manager to connect to.
        If the URL is already set it will not launch the dialog.
        :return:
        """
        if not self.connection_ip:
            with ui.dialog().props('persistent') as self.connectiondialog, ui.card():
                ui.label("No connection to minKNOW instance found.")
                self.connect_local = ui.button("Connect to Localhost", on_click=self.connect_to_localhost).props("outline").classes("shadow-lg")
                ui.label("Please select an IP address to connect to.")
                self.ipinput = ui.input(label='IP Address', placeholder='start typing',
                            on_change=lambda e: self._action_ip(e),validation={'Invalid IP': lambda value: self._check_ip(value)}).props('clearable')
                self.result = ui.label()
                self.connect_now = ui.button("Connect", on_click=self.add_positions).props("outline").classes("shadow-lg")
                self.connect_now.visible = False
            with ui.dialog().props('persistent') as self.position_choose, ui.card():
                ui.label("Connected to minKNOW.")
                ui.label("Please select a position to use.")
                self.choices=ui.row()
                self.access_device = ui.button("Choose Device", on_click=self.choose_device).props("outline").classes("shadow-lg")
                self.access_device.visible = False
            self.connectiondialog.open()
        pass

    def choose_device(self):
        self.position_choose.close()
        with self.holder:
            Position(self.selected_position, self.connection_ip)




    def enable_access(self, e):
        self.access_device.visible = True
        #self.position_choose.open()

    async def connect_to_localhost(self):
        self.connection_ip = "127.0.0.1"
        await self.add_positions()

    async def add_positions(self) -> None:
        """
        This is called when the connectionip variable is changed.
        It adds a tab for each position on the flowcell.
        It does this using a thread so that the app doesn't freeze
        while it waits for minKNOW to respond.
        """
        with disable(self.connect_now):
            ui.notify("Trying to connect to MinKNOW")
            await run.io_bound(self._connect_positions)
            if self.connected:
                ui.notify('Connection Successful.', type="positive")
                with self.connect_now:
                    ui.notify("Connected to MinKNOW - getting positions.")
                self.positions = list(self.manager.flow_cell_positions())
                self.connectiondialog.close()
                with self.choices:
                    ui.radio({item: str(item) for index, item in enumerate(self.positions)}, on_change=self.enable_access).bind_value(self, 'selected_position').props('inline')
                self.position_choose.open()
            else:
                self.connection_ip = None
                with self.connect_now:
                    ui.notify(
                        "Unable to connect.\nPlease check MinKNOW is running on this IP. ",
                        #title="Uh Oh!",
                        type="warning",
                        #timeout=10,
                    )
            ui.notify("Finished trying to connect to MinKNOW")

    def _connect_positions(self) -> None:
        if not self.manager and self.connection_ip:
            try:
                self.manager = Manager(host=self.connection_ip)
                self.connected = True
            except Exception:
                self.connected = False

    def watch_positions(self):
        while True:
            if self.manager:
                print(self.positions)
                for item in self.manager.rpc.watch_flow_cell_positions():
                    print("hello")
                    if self.tabs:
                        with self.tabs:
                            for position in item.additions:
                                ui.notify(f'{position.name} connected.', type="positive")
                                ui.tab(f"{position.name}")
                        with self.all_tab_panels:
                            for position in item.additions:
                                with ui.tab_panel(position.name):
                                    Position(position,self.connection_ip)
                        self.tabs.set_value(f"{position.name}")
            else:
                time.sleep(1)

    def setup_ui(self):
        with ui.card().classes("w-full"):
            self.label = ui.label("MinKNOW Information")
            self.holder = ui.row().classes("w-full")

