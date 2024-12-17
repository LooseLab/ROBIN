from nicegui import binding, ui, run, app
import time
import threading


# MinKNOW API Imports
import minknow_api.manager_pb2 as manager_pb2
from minknow_api.protocol_pb2 import ProtocolPhase, ProtocolRunInfo, ProtocolState
from readfish._utils import get_device
from robin.minknow_info.minKNOWhistograms import MinknowHistograms
from robin.minknow_info.minknow_info import Minknow_Info

from robin import theme

from contextlib import contextmanager
import ipaddress

# MinKNOW API Imports
from minknow_api.manager import Manager

from typing import Sequence
from pathlib import Path
import uuid

from datetime import datetime


# We need `find_protocol` to search for the required protocol given a kit + product code.
from minknow_api.tools import protocols


UNIQUE_ID: str = str(uuid.uuid4())


class ExperimentSpec(object):
    def __init__(self, position):
        self.position = position
        self.protocol_id = ""


ExperimentSpecs = Sequence[ExperimentSpec]


class ErrorChecker:
    def __init__(self, *elements) -> None:
        self.elements = elements

    @property
    def no_errors(self) -> bool:
        return all(
            validation(element.value)
            for element in self.elements
            for validation in element.validation.values()
        )


# Determine which protocol to run for each experiment, and add its ID to experiment_specs
def add_protocol_ids(experiment_specs, kit, basecall_config, expected_flowcell_id):
    for spec in experiment_specs:
        # Connect to the sequencing position:
        position_connection = spec.position.connect()
        flow_cell_info = position_connection.device.get_flow_cell_info()
        if flow_cell_info.flow_cell_id != expected_flowcell_id:
            ui.notify(
                f"Flowcell {expected_flowcell_id} is not found in position {spec.position}. Please check.",
                type="negative",
            )
            return
        if not flow_cell_info.has_flow_cell:
            ui.notify(
                "No flow cell present in position {}".format(spec.position),
                type="negative",
            )
            return

        product_code = flow_cell_info.user_specified_product_code
        if not product_code:
            product_code = flow_cell_info.product_code
        """
        device_connection: Connection,
        product_code: str,
        kit: str,
        basecalling: bool = False,
        basecall_config: Optional[str] = None,
        barcoding: bool = False,
        barcoding_kits: Optional[List[str]] = None,
        force_reload: bool = False,
        experiment_type: str = "sequencing",
        ) -> Optional[str]:
        """
        # Find the protocol identifier for the required protocol:
        protocol_info = protocols.find_protocol(
            position_connection,
            product_code=product_code,
            kit=kit,
            basecalling=True,
            basecall_config=basecall_config,
            barcoding=False,
            barcoding_kits=[],
            force_reload=True,
            experiment_type="sequencing",
        )

        if not protocol_info:
            ui.notify(
                "Failed to find protocol for position %s" % (spec.position),
                type="negative",
            )

            # print("Requested protocol:")
            # print("  product-code: %s" % product_code)
            # print("  kit: %s" % kit)
            # ui.notify("Failed to find protocol for position %s" % (spec.position))
            # ui.notify("Requested protocol:")
            # ui.notify("  product-code: %s" % product_code)
            # ui.notify("  kit: %s" % kit)
            return

        # Store the identifier for later:
        spec.protocol_id = protocol_info.identifier

    return True


@contextmanager
def disable(button: ui.button):
    button.disable()
    try:
        yield
    finally:
        button.enable()


class MinKNOWFish:
    """
    A way of handling a minknow connection.
    """

    def __init__(
        self,
        kit=None,
        centreID=None,
        experiment_duration=24,
        bed_file=None,
        basecall_config=None,
        reference=None,
        **kwargs,
    ):
        """
        This is called when the app is initialized.
        """
        # super().__init__(**kwargs)
        self.tabs = None
        self.manager = None
        self.positions = None
        self.selected_position = None
        self.connection_ip: str | None = None
        self.connected = False
        self.watchfolder = None
        self.devices = set()
        self.kit = kit
        self.basecall_config = basecall_config
        self.reference = reference
        self.bed_file = bed_file
        self.experiment_duration = experiment_duration
        self.centreID = centreID

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
            self.result.set_text("Connect to: " + e.value)
            self.connection_ip = e.value
        else:
            self.result.set_text("Invalid IP address: " + e.value)
            self.connect_now.visible = False

    async def auto_connect(self):
        self.connection_ip = "127.0.0.1"
        # await self.connect_to_localhost()
        await run.io_bound(self._connect_positions)
        self.positions = list(self.manager.flow_cell_positions())
        tabdict = {}
        self.focuspos = None

        def alert_click(position):
            if tabdict[position]["Counter"] != self.focuspos:
                ui.notify(f"Switching to {position}.")
                if self.focuspos is not None:
                    if self.temppos:
                        self.temppos.stop()
                    tabpanels.clear()
                self.focuspos = tabdict[position]["Counter"]
                with tabpanels:
                    with ui.tab_panel(position):
                        self.temppos = Position(
                            self.positions[tabdict[position]["Counter"]],
                            self.connection_ip,
                            self.centreID,
                            self.kit,
                            self.reference,
                            self.basecall_config,
                            self.bed_file,
                            self.experiment_duration,
                        )  # self.tabdict[position]["Position"].setup()
                        self.temppos.setup()
            else:
                ui.notify(f"{position} already selected.")

        tabs = ui.tabs().classes("w-full")

        tabpanels = ui.tab_panels(tabs).classes("w-full")

        for counter, position in enumerate(self.positions):
            tabdict[position.name] = {}
            with tabs:
                tabdict[position.name]["Tab"] = ui.tab(position.name).on(
                    "click", lambda: alert_click(tabs.value)
                )
                tabdict[position.name]["Counter"] = counter

        tabs.set_value(self.positions[0].name)
        alert_click(self.positions[0].name)

    def check_connection(self):
        """
        This is called when the app is initialized and checks
        to see if we have an active connection to a minKNOW instance.
        To do this it requests a URL for the manager to connect to.
        If the URL is already set it will not launch the dialog.
        :return:
        """
        if not self.connection_ip:
            with ui.dialog().props("persistent") as self.connectiondialog, ui.card():
                ui.label("No connection to minKNOW instance found.")
                self.connect_local = (
                    ui.button(
                        "Connect to Localhost", on_click=self.connect_to_localhost
                    )
                    .props("outline")
                    .classes("shadow-lg")
                )
                # ui.label("Please select an IP address to connect to.")
                # self.ipinput = ui.input(label='IP Address', placeholder='start typing',
                #            on_change=lambda e: self._action_ip(e),validation={'Invalid IP': lambda value: self._check_ip(value)}).props('clearable')
                # self.result = ui.label()
                self.connect_now = (
                    ui.button("Connect", on_click=self.add_positions)
                    .props("outline")
                    .classes("shadow-lg")
                )
                self.connect_now.visible = False
            with ui.dialog().props("persistent") as self.position_choose, ui.card():
                ui.label("Connected to minKNOW.")
                ui.label("Please select a position to use.")
                self.choices = ui.row()
                self.access_device = (
                    ui.button("Choose Device", on_click=self.choose_device)
                    .props("outline")
                    .classes("shadow-lg")
                )
                self.access_device.visible = False
            self.connectiondialog.open()
            # self.watchfolder = "/Users/mattloose/datasets/ds1305_Intraop0034_b"

        pass

    def choose_device(self):
        self.position_choose.close()
        with self.holder:
            self.minKNOW_display = Position(self.selected_position, self.connection_ip)
            self.minKNOW_display.setup()

    def enable_access(self, e):
        self.access_device.visible = True
        # self.position_choose.open()

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
                ui.notify("Connection Successful.", type="positive")
                with self.connect_now:
                    ui.notify("Connected to MinKNOW - getting positions.")
                self.positions = list(self.manager.flow_cell_positions())
                self.connectiondialog.close()
                with self.choices:
                    ui.radio(
                        {item: str(item) for index, item in enumerate(self.positions)},
                        on_change=self.enable_access,
                    ).bind_value(self, "selected_position").props("inline")
                self.position_choose.open()
            else:
                self.connection_ip = None
                with self.connect_now:
                    ui.notify(
                        "Unable to connect.\nPlease check MinKNOW is running on this IP. ",
                        # title="Uh Oh!",
                        type="warning",
                        # timeout=10,
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
                for item in self.manager.rpc.watch_flow_cell_positions():
                    if self.tabs:
                        with self.tabs:
                            for position in item.additions:
                                ui.notify(
                                    f"{position.name} connected.", type="positive"
                                )
                                ui.tab(f"{position.name}")
                        with self.all_tab_panels:
                            for position in item.additions:
                                with ui.tab_panel(position.name):
                                    Position(position, self.connection_ip)
                        self.tabs.set_value(f"{position.name}")
            else:
                time.sleep(1)

    def setup_ui(self):
        self.holder = ui.row().classes("w-full")


class Position(MinKNOWFish):
    def __init__(
        self,
        position,
        ip,
        centreID,
        kit,
        reference,
        basecall_config,
        bed_file,
        experiment_duration,
    ):
        # super().__init__()
        print("starting position")
        self.ip = ip
        self.position = get_device(position.name, host=self.ip)
        self.connection = self.position.connect()
        self.Run_ID = None
        self.yield_summary_info = None
        self.acquisition_run_info = None
        self.protocol_run_info = None
        self.flowcell_info = None
        self.minknow_info_pane = None
        self.watchfolder = None
        self.centreID = centreID
        self.kit = kit
        self.reference = reference
        self.basecall_config = basecall_config
        self.bed_file = bed_file
        self.experiment_duration = experiment_duration
        self.shutdown_event = threading.Event()

        self.check_flowcell = threading.Thread(target=self.watch_flowcell, args=())
        self.check_flowcell.daemon = True
        self.check_flowcell.start()

        self.check_run = threading.Thread(target=self.stream_current_run, args=())
        self.check_run.daemon = True
        self.check_run.start()

        self.check_instance = threading.Thread(
            target=self.stream_instance_activity, args=()
        )
        self.check_instance.daemon = True
        self.check_instance.start()

        self.timer = ui.timer(1, self._worker)

        # Convert yield_data to use bindable properties
        self.yield_data = {
            'read_count': binding.BindableProperty(0),
            'selected_raw_samples': binding.BindableProperty(0),
            'selected_events': binding.BindableProperty(0),
            'estimated_selected_bases': binding.BindableProperty(0),
            'basecalled_pass_read_count': binding.BindableProperty(0),
            'basecalled_fail_read_count': binding.BindableProperty(0),
            'basecalled_pass_bases': binding.BindableProperty(0),
            'basecalled_fail_bases': binding.BindableProperty(0),
            'fraction_basecalled': binding.BindableProperty(0)
        }
        
        # Store metrics UI elements
        self.metric_elements = {}

    def stop(self):
        self.shutdown_event.set()
        # print(self.shutdown_event)

    def stopped(self):
        return self.shutdown_event.isSet()

    def _worker(self):
        if self.minknow_info_pane:
            if hasattr(self.minknow_info_pane, "output_folder"):
                if self.minknow_info_pane.show:
                    if self.minknow_info_pane.output_folder != self.watchfolder:
                        self.watchfolder = self.minknow_info_pane.output_folder

    def setup(self):
        with ui.card().classes("w-full"):
            with ui.card().classes("w-full drop-shadow"):
                """
                self,
                position,
                centreID,
                kit,
                reference,
                basecall_config,
                bed_file,
                experiment_duration,
                dev=False,
                """
                self.minknow_info_pane = Minknow_Info(
                    self.position,
                    self.centreID,
                    self.kit,
                    self.reference,
                    self.basecall_config,
                    self.bed_file,
                    self.experiment_duration,
                )

                with ui.column().bind_visibility_from(self.minknow_info_pane, "show"):
                    with ui.row():
                        with ui.card().classes("w-full p-4 drop-shadow-lg"):
                            with ui.column().classes("gap-2"):
                                ui.label("Device Information") \
                                    .classes("text-lg font-medium mb-2")
                                
                                with ui.grid().classes("grid-cols-2 gap-4"):
                                    ui.label().bind_text_from(
                                        self, "Run_ID",
                                        backward=lambda n: f"Run ID: {n}"
                                    ).classes("text-sm text-gray-600")
                                    
                                    ui.label().bind_text_from(
                                        self.position, "name",
                                        backward=lambda n: f"Name: {n}"
                                    ).classes("text-sm text-gray-600")
                                    
                                    ui.label().bind_text_from(
                                        self.position, "device_type",
                                        backward=lambda n: f"Device Type: {n}"
                                    ).classes("text-sm text-gray-600")
                                    
                                    ui.label().bind_text_from(
                                        self.position, "host",
                                        backward=lambda n: f"Host: {n}" 
                                    ).classes("text-sm text-gray-600")
                                    
                                    ui.label().bind_text_from(
                                        self.position, "running",
                                        backward=lambda n: f"Running: {n}"
                                    ).classes("text-sm font-medium").props('color=green' if self.position.running else 'color=red')
                                    
                                    ui.label().bind_text_from(
                                        self.position, "state",
                                        backward=lambda n: f"State: {n}"
                                    ).classes("text-sm font-medium")

                        self.protocol_run_info = ui.column().classes(
                            "p-4 border border-gray-200 rounded-lg bg-gray-50"
                        )

                    self.flowcell_info = ui.row()

                    self.yield_summary_info = ui.row().classes(
                        "p-4 border-dashed border-2"
                    )

                    self.acquisition_run_info = ui.column()
                    
                    self.minknowhistogram = MinknowHistograms(self.position)


    def watch_flowcell(self):
        while not self.stopped():
            if self.connection:
                if self.flowcell_info:
                    for stuff in self.connection.device.stream_flow_cell_info():
                        self.flowcell_info.clear()
                        with self.flowcell_info:
                            with ui.card().classes("drop-shadow"):
                                for field in stuff.DESCRIPTOR.fields:
                                    ui.label(
                                        f"{field.name}: {getattr(stuff, field.name)}"
                                    ).classes("text-sm font-medium")
            else:
                time.sleep(1)
        print("watch flowcell saw a shutdown")

    def stream_current_run(self) -> None:
        first_run = True
        while not self.stopped():
            self._stream_current_run = (
                self.connection.acquisition.watch_current_acquisition_run()
            )
            for info in self._stream_current_run:
                if info.state == manager_pb2.FlowCellPosition.State.STATE_RUNNING:
                    self.Run_ID = info.run_id
                    self.Live_Run = True
                    msg = f"Run {self.Run_ID} seen on {self.position}"
                    title = "MinKNOW Info"
                    severity = "information"
                    timeout = 5
                else:
                    self.Run_ID = None
                    self.Live_Run = False
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

    def stream_instance_activity(self) -> None:
        while not self.stopped():
            self._stream_instance_activity = (
                self.connection.instance.stream_instance_activity()
            )
            for info in self._stream_instance_activity:
                if self.shutdown_event.is_set():
                    break

                # Update yield data when we receive yield summary
                if info.HasField("yield_summary"):
                    yield_info = info.yield_summary
                    print(yield_info)
                    # Update bindable properties
                    for field in self.yield_data:
                        print(field, getattr(yield_info, field))
                        self.yield_data[field].value = getattr(yield_info, field)
                    print(self.yield_data)

                    # Update history
                    if hasattr(self, '_yield_history'):
                        self._yield_history.append({
                            'timestamp': datetime.now().strftime('%H:%M:%S'),
                            'reads': self.yield_data['read_count'].value,
                            'bases': self.yield_data['estimated_selected_bases'].value
                        })
                    else:
                        self._yield_history = [{
                            'timestamp': datetime.now().strftime('%H:%M:%S'),
                            'reads': self.yield_data['read_count'].value,
                            'bases': self.yield_data['estimated_selected_bases'].value
                        }]

                # Update acquisition run info display when we receive that update
                if self.acquisition_run_info:
                    if info.HasField("acquisition_run_info"):
                        self.acquisition_run_info.clear()
                        with self.acquisition_run_info:
                            # Main stats card
                            with ui.card().classes("w-full p-4 drop-shadow-lg"):
                                with ui.row().classes("justify-between"):
                                    # Left column - Run Info
                                    with ui.column().classes("gap-2 flex-1"):
                                        ui.label("Run Information").classes("text-lg font-medium")
                                        
                                        run_info = info.acquisition_run_info
                                        print(run_info)
                                        ui.label(f"Run ID: {run_info.run_id}").classes("text-sm")
                                        ui.label(f"State: {run_info.state}").classes("text-sm font-medium")
                                        
                                        start_time = datetime.fromtimestamp(
                                            run_info.start_time.seconds + run_info.start_time.nanos / 1e9
                                        )
                                        ui.label(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}").classes("text-sm")
                                        
                                        if hasattr(run_info, 'data_read_start_time'):
                                            read_start = datetime.fromtimestamp(
                                                run_info.data_read_start_time.seconds + 
                                                run_info.data_read_start_time.nanos / 1e9
                                            )
                                            ui.label(f"Reading Since: {read_start.strftime('%Y-%m-%d %H:%M:%S')}").classes("text-sm")

                                    # Right column - Yield Summary
                                    if hasattr(run_info, 'yield_summary'):
                                        #print(run_info)
                                        with ui.column().classes("gap-2 flex-1"):
                                            ui.label("Yield Summary").classes("text-lg font-medium")
                                            yield_info = run_info.yield_summary
                                            print(yield_info)
                                            
                                            metrics = [
                                                ("Reads", 'read_count'),
                                                ("Raw Samples", 'selected_raw_samples', lambda x: f"{x:,}"),
                                                ("Events", 'selected_events', lambda x: f"{x:,}"),
                                                ("Est. Bases", 'estimated_selected_bases', lambda x: f"{x:,}"),
                                                ("Pass Reads", 'basecalled_pass_read_count', lambda x: f"{x:,}"),
                                                ("Fail Reads", 'basecalled_fail_read_count', lambda x: f"{x:,}"),
                                                ("Pass Bases", 'basecalled_pass_bases', lambda x: f"{x:,}"),
                                                ("Fail Bases", 'basecalled_fail_bases', lambda x: f"{x:,}"),
                                                ("Percent Basecalled", 'fraction_basecalled', lambda x: f"{x:.1%}")
                                            ]
                                            
                                            for label, bind_key, *formatter in metrics:
                                                format_func = formatter[0] if formatter else str
                                                with ui.row().classes("justify-between w-full"):
                                                    ui.label(f"{label}:").classes("text-sm")
                                                    ui.label().bind_text_from(
                                                        self.yield_data[bind_key], 
                                                        'value',
                                                        backward=format_func
                                                    ).classes("text-sm font-medium")

                            # Add charts if we have historical data
                            if hasattr(self, '_yield_history'):
                                self._yield_history.append({
                                    'timestamp': datetime.now().strftime('%H:%M:%S'),
                                    'reads': yield_info.read_count,
                                    'bases': yield_info.estimated_selected_bases
                                })
                            else:
                                self._yield_history = [{
                                    'timestamp': datetime.now().strftime('%H:%M:%S'),
                                    'reads': yield_info.read_count,
                                    'bases': yield_info.estimated_selected_bases
                                }]

                            # Plot yield over time
                            with ui.card().classes("w-full p-4 mt-4 drop-shadow-lg"):
                                ui.label("Sequencing Progress").classes("text-lg font-medium mb-4")
                                
                                ui.echart({
                                    'tooltip': {'trigger': 'axis'},
                                    'legend': {'data': ['Reads', 'Bases (M)']},
                                    'xAxis': {
                                        'type': 'category',
                                        'data': [d['timestamp'] for d in self._yield_history]
                                    },
                                    'yAxis': [
                                        {'type': 'value', 'name': 'Reads'},
                                        {'type': 'value', 'name': 'Bases (M)', 'alignTicks': True}
                                    ],
                                    'series': [
                                        {
                                            'name': 'Reads',
                                            'type': 'line',
                                            'data': [d['reads'] for d in self._yield_history],
                                            'smooth': True
                                        },
                                        {
                                            'name': 'Bases (M)',
                                            'type': 'line',
                                            'yAxisIndex': 1,
                                            'data': [d['bases']/1_000_000 for d in self._yield_history],
                                            'smooth': True
                                        }
                                    ]
                                }).classes('w-full h-64')

                            # Configuration Summary if available
                            if hasattr(run_info, 'config_summary'):
                                with ui.card().classes("w-full p-4 mt-4 drop-shadow-lg"):
                                    ui.label("Configuration").classes("text-lg font-medium")
                                    config = run_info.config_summary
                                    ui.label(f"Purpose: {config.purpose}").classes("text-sm")
                                    ui.label(f"Output Directory: {config.reads_directory}").classes("text-sm")
                                    ui.label(f"Sample Rate: {config.sample_rate} Hz").classes("text-sm")
                                    ui.label(f"Channel Count: {config.channel_count}").classes("text-sm")

                        self.start_time = start_time
                        
                if info.HasField("basecall_speed"):
                    self.Mean_Basecall_Speed = info.basecall_speed.mean_basecall_speed
                else:
                    self.Mean_Basecall_Speed = None
                if info.HasField("n50"):
                    self.N50 = info.n50.n50
                    self.Estimated_N50 = info.n50.estimated_n50
                if self.protocol_run_info:
                    if info.HasField("protocol_run_info"):
                        self.protocol_run_info.clear()
                        with self.protocol_run_info:
                            with ui.card().classes("w-full p-4 drop-shadow-lg"):
                                with ui.column().classes("gap-2"):
                                    # Header section
                                    ui.label("Protocol Information").classes("text-lg font-medium mb-2")
                                    
                                    # Main info grid
                                    with ui.grid().classes("grid-cols-2 gap-4"):
                                        if info.protocol_run_info.state == ProtocolState.PROTOCOL_COMPLETED:
                                            ui.label("Status: Protocol Completed").classes("text-sm font-medium text-green-600")
                                        else:
                                            # Basic Protocol Info
                                            ui.label(f"Protocol ID: {info.protocol_run_info.protocol_id}").classes("text-sm text-gray-600")
                                            ui.label(f"Run ID: {info.protocol_run_info.run_id}").classes("text-sm text-gray-600")
                                            
                                            # User Info
                                            if hasattr(info.protocol_run_info, 'user_info'):
                                                ui.label(f"Group ID: {info.protocol_run_info.user_info.protocol_group_id.value}").classes("text-sm text-gray-600")
                                                ui.label(f"Sample ID: {info.protocol_run_info.user_info.sample_id.value}").classes("text-sm text-gray-600")
                                                ui.label(f"Flow Cell ID: {info.protocol_run_info.user_info.user_specified_flow_cell_id.value}").classes("text-sm text-gray-600")
                                            
                                            # Meta Info
                                            if hasattr(info.protocol_run_info, 'meta_info'):
                                                meta = info.protocol_run_info.meta_info
                                                if hasattr(meta, 'tags'):
                                                    for tag in ['kit', 'flow cell', 'speed']:
                                                        if tag in meta.tags:
                                                            value = meta.tags[tag].string_value or meta.tags[tag].int_value
                                                            ui.label(f"{tag.title()}: {value}").classes("text-sm text-gray-600")
                                            
                                            # Phase and Output Info
                                            ui.label(f"Phase: {info.protocol_run_info.phase}").classes("text-sm font-medium")
                                            ui.label(f"Output Path: {info.protocol_run_info.output_path}").classes("text-sm text-gray-600")
                                            
                                            # Software Versions (detailed collapsible section)
                                            if hasattr(info.protocol_run_info, 'software_versions'):
                                                with ui.expansion("Software Versions", icon="developer_board").classes("w-full mt-4"):
                                                    with ui.grid().classes("grid-cols-2 gap-2 mt-2"):
                                                        sw = info.protocol_run_info.software_versions
                                                        
                                                        # MinKNOW version
                                                        with ui.card().classes("p-2"):
                                                            ui.label("MinKNOW").classes("text-sm font-medium")
                                                            ui.label(f"Version: {sw.minknow.full}").classes("text-xs text-gray-600")
                                                            ui.label(f"Major: {sw.minknow.major}").classes("text-xs text-gray-600")
                                                            ui.label(f"Minor: {sw.minknow.minor}").classes("text-xs text-gray-600")
                                                            ui.label(f"Patch: {sw.minknow.patch}").classes("text-xs text-gray-600")
                                                        
                                                        # Bream version
                                                        with ui.card().classes("p-2"):
                                                            ui.label("Bream").classes("text-sm font-medium")
                                                            ui.label(f"Version: {sw.bream}").classes("text-xs text-gray-600")
                                                        
                                                        # Distribution info
                                                        with ui.card().classes("p-2"):
                                                            ui.label("Distribution").classes("text-sm font-medium")
                                                            ui.label(f"Version: {sw.distribution_version}").classes("text-xs text-gray-600")
                                                            ui.label(f"Status: {sw.distribution_status}").classes("text-xs text-gray-600")
                                                            ui.label(f"Type: {sw.installation_type}").classes("text-xs text-gray-600")
                                                        
                                                        # Basecaller info
                                                        with ui.card().classes("p-2"):
                                                            ui.label("Basecaller").classes("text-sm font-medium")
                                                            ui.label(f"Connected: {sw.basecaller_connected_version}").classes("text-xs text-gray-600")
                                                            ui.label(f"Build: {sw.basecaller_build_version}").classes("text-xs text-gray-600")
                                                        
                                                        # Protocol config
                                                        with ui.card().classes("p-2"):
                                                            ui.label("Protocol").classes("text-sm font-medium")
                                                            ui.label(f"Configuration: {sw.protocol_configuration}").classes("text-xs text-gray-600")
                                        
                                        # Arguments section (collapsible)
                                        if hasattr(info.protocol_run_info, 'args') and info.protocol_run_info.args:
                                            with ui.expansion("Protocol Arguments", icon="settings").classes("w-full mt-4"):
                                                with ui.column().classes("gap-1"):
                                                    for arg in info.protocol_run_info.args:
                                                        ui.label(arg).classes("text-xs text-gray-500")
                

    # def shutdown(self):
    #    self.shutdown_event.set()
    #    self.check_flowcell.join()
    #    self.check_run.join()
    #    self.check_instance.join()
    #    self.timer.stop()


@ui.page("/", response_timeout=30)
async def index_page() -> None:
    # initial_ip = "127.0.0.1"
    # my_connection = Manager(host=initial_ip)
    with theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ):
        # my_connection.connect_to_minknow()
        # positions = list(my_connection.flow_cell_positions())
        # ui.label(f"{positions[0]}")
        display_object = MinKNOWFish()
        display_object.setup_ui()
        await display_object.auto_connect()
        # ui.label().bind_text_from(
        #    self.minknow_connection,
        ##    "connection_ip",
        #    backward=lambda n: f"Connected to: {n}",
        # )


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    # Add some custom CSS because - why not!
    ui.add_css(
        """
        .shadows-into light-regular {
            font-family: "Shadows Into Light", cursive;
            font-weight: 800;
            font-style: normal;
        }
    """
    )
    # Register some fonts that we might need later on.
    app.add_static_files("/fonts", str(Path(__file__).parent.parent / "fonts"))

    ui.run(
        port=port,
        reload=reload,
        title="Readfish NiceGUI",
        storage_secret="waynesworld",
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
