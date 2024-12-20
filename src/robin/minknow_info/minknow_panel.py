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
from minknow_api.tools import protocols
import grpc

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
    ):
        self.connection_ip = None
        self.manager = None
        self.connected = False
        self.positions = []
        self.position_list = []
        self.selected_position = None
        
        # Configuration parameters
        self.centreID = centreID
        self.kit = kit
        self.reference = reference
        self.basecall_config = basecall_config
        self.bed_file = bed_file
        self.experiment_duration = experiment_duration
        
        # UI elements that will be initialized later
        self.position_choose = None
        self.holder = None
        self.connect_now = None
        self.choices = None
        self.connectiondialog = None
        self.access_device = None
        self.tabs = None
        self.all_tab_panels = None

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
        print(f"[choose_device] Setting up device for position: {self.selected_position}")
        self.position_choose.close()
        
        try:
            print("[choose_device] Creating Position instance")
            with self.holder:
                self.minKNOW_display = Position(
                    self.selected_position,
                    self.connection_ip,
                    self.centreID,
                    self.kit,
                    self.reference,
                    self.basecall_config,
                    self.bed_file,
                    self.experiment_duration,
                )
                print("[choose_device] Position instance created, calling setup")
                self.minKNOW_display.setup(self.holder)
                print("[choose_device] Setup completed")
        except Exception as e:
            print(f"[choose_device] Error: {str(e)}")
            print(f"[choose_device] Error type: {type(e)}")
            import traceback
            traceback.print_exc()
            ui.notify(f"Error setting up device: {str(e)}", type="negative")

    def enable_access(self, e):
        self.access_device.visible = True
        # self.position_choose.open()

    async def connect_to_localhost(self):
        print("[connect_to_localhost] Attempting to connect to localhost...")
        try:
            self.connection_ip = "127.0.0.1"
            print("[connect_to_localhost] Set connection_ip to 127.0.0.1")
            await self.add_positions()
            print("[connect_to_localhost] Successfully connected to localhost")
            return True
        except Exception as e:
            print(f"[connect_to_localhost] Error: {str(e)}")
            print(f"[connect_to_localhost] Error type: {type(e)}")
            import traceback
            traceback.print_exc()
            return False

    async def add_positions(self) -> None:
        print("[add_positions] Starting...")
        with disable(self.connect_now):
            ui.notify("Trying to connect to MinKNOW")
            try:
                print("[add_positions] Calling _connect_positions...")
                await run.io_bound(self._connect_positions)
                print("[add_positions] _connect_positions completed")
                
                if self.connected:
                    print("[add_positions] Connection successful")
                    ui.notify("Connection Successful.", type="positive")
                    self.connectiondialog.close()
                    
                    print("[add_positions] Creating radio buttons for positions")
                    with self.choices:
                        ui.radio(
                            {item: str(item) for item in self.positions},
                            on_change=self.enable_access,
                        ).bind_value(self, "selected_position").props("inline")
                    print("[add_positions] Opening position chooser")
                    self.position_choose.open()
                else:
                    print("[add_positions] Connection failed")
                    self.connection_ip = None
                    ui.notify(
                        "Unable to connect.\nPlease check MinKNOW is running on this IP.",
                        type="warning",
                    )
            except Exception as e:
                print(f"[add_positions] Error: {str(e)}")
                print(f"[add_positions] Error type: {type(e)}")
                import traceback
                traceback.print_exc()
                ui.notify(f"Connection error: {str(e)}", type="negative")
            finally:
                print("[add_positions] Finished")
                ui.notify("Finished trying to connect to MinKNOW")

    def _connect_positions(self):
        print("Starting MinKNOWFish connection...")
        try:
            print(f"Attempting to connect to MinKNOW at {self.connection_ip}")
            self.manager = Manager(host=self.connection_ip)
            print("Manager connection established")
            
            try:
                print("Getting version info...")
                version_info = self.manager.version
                print(f"MinKNOW version info: {version_info}")
            except Exception as ve:
                print(f"Warning: Could not get version info: {str(ve)}")
                print("Continuing anyway...")
            
            print("Getting flow cell positions...")
            self.positions = list(self.manager.flow_cell_positions())
            print(f"Found positions: {self.positions}")
            self.connected = True
            
        except Exception as e:
            print(f"MinKNOWFish connection error: {str(e)}")
            print(f"Error type: {type(e)}")
            import traceback
            traceback.print_exc()
            self.connected = False
            raise

    def watch_positions(self):
        """Watch for new flow cell positions and update UI accordingly"""
        print("[watch_positions] Starting position watcher...")
        try:
            while True:
                if self.manager:
                    print("[watch_positions] Watching for position changes...")
                    for item in self.manager.rpc.watch_flow_cell_positions():
                        if self.tabs:
                            print(f"[watch_positions] Processing position changes: {item}")
                            with self.tabs:
                                for position in item.additions:
                                    ui.notify(
                                        f"{position.name} connected.", type="positive"
                                    )
                                    ui.tab(f"{position.name}")
                            with self.all_tab_panels:
                                for position in item.additions:
                                    with ui.tab_panel(position.name):
                                        # Instead of creating a new Position instance,
                                        # just update the UI to show the position is available
                                        ui.label(f"Position {position.name} is available")
                                        ui.button(
                                            "Connect to this position",
                                            on_click=lambda p=position: self.connect_to_position(p)
                                        )
                            self.tabs.set_value(f"{position.name}")
                else:
                    print("[watch_positions] No manager available, waiting...")
                    time.sleep(1)
        except Exception as e:
            print(f"[watch_positions] Error: {str(e)}")
            print(f"[watch_positions] Error type: {type(e)}")
            import traceback
            traceback.print_exc()

    def connect_to_position(self, position):
        """Connect to a specific position when requested"""
        print(f"[connect_to_position] Connecting to position: {position.name}")
        try:
            self.selected_position = position
            self.choose_device()
        except Exception as e:
            print(f"[connect_to_position] Error: {str(e)}")
            print(f"[connect_to_position] Error type: {type(e)}")
            import traceback
            traceback.print_exc()
            ui.notify(f"Error connecting to position: {str(e)}", type="negative")

    def setup_ui(self):
        self.holder = ui.row().classes("w-full")


class Position(MinKNOWFish):
    def __init__(
        self,
        position,
        ip,
        centreID=None,
        kit=None,
        reference=None,
        basecall_config=None,
        bed_file=None,
        experiment_duration=24,
    ):
        print("[Position.__init__] Starting position initialization...")
        # Call parent class init with configuration parameters
        print("[Position.__init__] Calling parent class init")
        super().__init__(
            kit=kit,
            centreID=centreID,
            experiment_duration=experiment_duration,
            bed_file=bed_file,
            basecall_config=basecall_config,
            reference=reference,
        )
        print("[Position.__init__] Parent class init completed")
        
        self.ip = ip
        self.yield_data = []  # Initialize yield_data
        self.threads = []  # Keep track of background threads
        self._initialized = False  # Track initialization state
        
        try:
            print(f"[Position.__init__] Using position object: {position}")
            print(f"[Position.__init__] Position type: {type(position)}")
            print(f"[Position.__init__] Position attributes: {dir(position)}")
            
            # Use the position object directly
            self.position = position
            print("[Position.__init__] Attempting to connect to position...")
            self.connection = self.position.connect()
            print("[Position.__init__] Position connection established")
            
            # Initialize other attributes
            print("[Position.__init__] Initializing attributes")
            self.Run_ID = None
            self.yield_summary_info = None
            self.acquisition_run_info = None
            self.protocol_run_info = None
            self.flowcell_info = None
            self.minknow_info_pane = None
            self.watchfolder = None
            self.shutdown_event = threading.Event()
            self._initialized = True
            print("[Position.__init__] Position initialization completed successfully")
        except Exception as e:
            print(f"[Position.__init__] Error: {str(e)}")
            print(f"[Position.__init__] Error type: {type(e)}")
            import traceback
            traceback.print_exc()
            raise

    def setup(self, container=None):
        """Setup UI elements in the correct context"""
        print("[Position.setup] Starting setup...")
        if not self._initialized:
            print("[Position.setup] Error: Position not properly initialized")
            raise RuntimeError("Position not properly initialized")
            
        try:
            print("[Position.setup] Creating UI elements")
            # Create timer in the global UI context
            self.timer = ui.timer(1, self._worker)
            
            # Create UI elements in the provided container or global context
            target = container if container is not None else ui
            print(f"[Position.setup] Using container: {target}")
            
            print("[Position.setup] Creating UI cards")
            with target.card().classes("w-full"):
                # Status and basic info section
                with ui.card().classes("w-full p-4 mb-4"):
                    ui.label("Status Information").classes("text-h6 mb-2")
                    self.status_label = ui.label("Status: Initializing...").classes("text-body1")
                    self.flowcell_label = ui.label("Flowcell: -").classes("text-body1")
                    self.run_id_label = ui.label("Run ID: -").classes("text-body1")
                
                # Yield statistics section
                with ui.card().classes("w-full p-4 mb-4"):
                    ui.label("Yield Statistics").classes("text-h6 mb-2")
                    with ui.row().classes("justify-between"):
                        with ui.column().classes("w-1/2"):
                            self.read_count_label = ui.label("Reads: 0").classes("text-body1")
                            self.bases_label = ui.label("Bases: 0").classes("text-body1")
                        with ui.column().classes("w-1/2"):
                            self.pass_reads_label = ui.label("Pass Reads: 0").classes("text-body1")
                            self.fail_reads_label = ui.label("Fail Reads: 0").classes("text-body1")
                
                # MinKNOW info section
                with ui.card().classes("w-full drop-shadow"):
                    print("[Position.setup] Creating Minknow_Info instance")
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
                        
                        print("[Position.setup] Creating MinknowHistograms instance")
                        self.minknowhistogram = MinknowHistograms(self.position)
                
            print("[Position.setup] Starting background tasks")
            self._start_background_tasks()
            print("[Position.setup] Setup completed successfully")
            
        except Exception as e:
            print(f"[Position.setup] Error: {str(e)}")
            print(f"[Position.setup] Error type: {type(e)}")
            import traceback
            traceback.print_exc()
            raise

    def _start_background_tasks(self):
        """Start background tasks with proper error handling"""
        print("[Position._start_background_tasks] Starting...")
        if not self._initialized:
            print("[Position._start_background_tasks] Error: Position not properly initialized")
            raise RuntimeError("Position not properly initialized")
            
        tasks = [
            (self.watch_flowcell, "Flowcell Monitor"),
            (self.stream_current_run, "Run Monitor"),
            (self.stream_instance_activity, "Activity Monitor")
        ]
        
        for task_func, task_name in tasks:
            try:
                print(f"[Position._start_background_tasks] Starting {task_name}")
                thread = threading.Thread(
                    target=self._run_task_with_error_handling,
                    args=(task_func, task_name),
                    daemon=True
                )
                thread.start()
                self.threads.append(thread)
                print(f"[Position._start_background_tasks] Started {task_name} thread")
            except Exception as e:
                print(f"[Position._start_background_tasks] Error starting {task_name}: {str(e)}")
                print(f"[Position._start_background_tasks] Error type: {type(e)}")
                import traceback
                traceback.print_exc()
                
    def _run_task_with_error_handling(self, task_func, task_name):
        """Wrapper to run a task with error handling"""
        print(f"[Position._run_task_with_error_handling] Starting {task_name}")
        if not self._initialized:
            print(f"[Position._run_task_with_error_handling] Warning: Attempting to run {task_name} before initialization")
            return
            
        try:
            print(f"[Position._run_task_with_error_handling] Executing {task_name}")
            task_func()
            print(f"[Position._run_task_with_error_handling] {task_name} completed")
        except Exception as e:
            print(f"[Position._run_task_with_error_handling] Error in {task_name}: {str(e)}")
            print(f"[Position._run_task_with_error_handling] Error type: {type(e)}")
            import traceback
            traceback.print_exc()

    def stop(self):
        """Stop all background tasks"""
        print("Stopping position tasks...")
        self.shutdown_event.set()
        for thread in self.threads:
            thread.join(timeout=1.0)
        print("Position tasks stopped")

    def _worker(self):
        """Timer callback to update UI elements"""
        print("[Position._worker] Starting worker update")
        try:
            if not self._initialized:
                print("[Position._worker] Warning: Position not initialized")
                return
                
            if not hasattr(self, 'status_label'):
                print("[Position._worker] Warning: UI elements not initialized")
                return
                
            # Update status
            try:
                state = self.position.state
                self.status_label.text = f"Status: {state}"
            except Exception as e:
                print(f"[Position._worker] Error updating status: {str(e)}")
                
            # Update flowcell info
            try:
                if hasattr(self, 'flowcell_info') and self.flowcell_info:
                    self.flowcell_label.text = f"Flowcell: {self.flowcell_info}"
            except Exception as e:
                print(f"[Position._worker] Error updating flowcell info: {str(e)}")
                
            # Update run ID
            try:
                if hasattr(self, 'Run_ID') and self.Run_ID:
                    self.run_id_label.text = f"Run ID: {self.Run_ID}"
            except Exception as e:
                print(f"[Position._worker] Error updating run ID: {str(e)}")
                
            # Update yield statistics
            try:
                if hasattr(self, 'yield_summary_info') and self.yield_summary_info:
                    # Convert yield summary info to a dictionary if it's not already
                    yield_info = {}
                    if hasattr(self.yield_summary_info, '__dict__'):
                        yield_info = self.yield_summary_info.__dict__
                    elif hasattr(self.yield_summary_info, 'DESCRIPTOR'):
                        # Handle protobuf message
                        for field in self.yield_summary_info.DESCRIPTOR.fields:
                            yield_info[field.name] = getattr(self.yield_summary_info, field.name)
                    else:
                        # Try to use it as is
                        yield_info = self.yield_summary_info

                    # Update labels with safe access
                    def safe_get(d, key, default=0):
                        try:
                            return d.get(key, default) if hasattr(d, 'get') else getattr(d, key, default)
                        except:
                            return default

                    self.read_count_label.text = f"Reads: {safe_get(yield_info, 'read_count'):,}"
                    self.bases_label.text = f"Bases: {safe_get(yield_info, 'selected_raw_samples'):,}"
                    self.pass_reads_label.text = f"Pass Reads: {safe_get(yield_info, 'basecalled_pass_read_count'):,}"
                    self.fail_reads_label.text = f"Fail Reads: {safe_get(yield_info, 'basecalled_fail_read_count'):,}"
            except Exception as e:
                print(f"[Position._worker] Error updating yield stats: {str(e)}")
                print(f"[Position._worker] Yield info type: {type(self.yield_summary_info)}")
                print(f"[Position._worker] Yield info: {self.yield_summary_info}")
                import traceback
                traceback.print_exc()
                
        except Exception as e:
            print(f"[Position._worker] Error in worker: {str(e)}")
            print(f"[Position._worker] Error type: {type(e)}")
            import traceback
            traceback.print_exc()

    def stopped(self):
        """Check if the position has been stopped"""
        return self.shutdown_event.is_set()

    def watch_flowcell(self):
        """Monitor flowcell information"""
        print("[Position.watch_flowcell] Starting flowcell monitoring")
        try:
            while not self.stopped():
                if self.connection:
                    if hasattr(self, 'flowcell_info') and self.flowcell_info is not None:
                        for stuff in self.connection.device.stream_flow_cell_info():
                            if self.stopped():
                                break
                            self.flowcell_info.clear()
                            with self.flowcell_info:
                                with ui.card().classes("drop-shadow"):
                                    for field in stuff.DESCRIPTOR.fields:
                                        ui.label(
                                            f"{field.name}: {getattr(stuff, field.name)}"
                                        ).classes("text-sm font-medium")
                else:
                    time.sleep(1)
            print("[Position.watch_flowcell] Flowcell monitoring stopped")
        except Exception as e:
            print(f"[Position.watch_flowcell] Error: {str(e)}")
            print(f"[Position.watch_flowcell] Error type: {type(e)}")
            import traceback
            traceback.print_exc()

    def stream_current_run(self) -> None:
        """Monitor current run information"""
        print("[Position.stream_current_run] Starting run monitoring")
        try:
            first_run = True
            while not self.stopped():
                try:
                    self._stream_current_run = (
                        self.connection.acquisition.watch_current_acquisition_run()
                    )
                    for info in self._stream_current_run:
                        if self.stopped():
                            break
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
                except Exception as e:
                    print(f"[Position.stream_current_run] Inner loop error: {str(e)}")
                    time.sleep(1)
            print("[Position.stream_current_run] Run monitoring stopped")
        except Exception as e:
            print(f"[Position.stream_current_run] Error: {str(e)}")
            print(f"[Position.stream_current_run] Error type: {type(e)}")
            import traceback
            traceback.print_exc()

    def stream_instance_activity(self) -> None:
        """Monitor instance activity"""
        print("[Position.stream_instance_activity] Starting instance monitoring")
        try:
            while not self.stopped():
                try:
                    self._stream_instance_activity = (
                        self.connection.instance.stream_instance_activity()
                    )
                    for info in self._stream_instance_activity:
                        if self.stopped():
                            break
                        if info.HasField("device_info"):
                            self.Channel_Count = info.device_info.device_info.max_channel_count
                        if info.HasField("acquisition_run_info"):
                            self.basecalling_config_filename = (
                                info.acquisition_run_info.config_summary.basecalling_config_filename
                            )
                            self.start_time = datetime.fromtimestamp(
                                info.acquisition_run_info.start_time.seconds
                                + info.acquisition_run_info.start_time.nanos / 1e9
                            )
                            self.output_folder = (
                                info.acquisition_run_info.config_summary.reads_directory
                            )
                        if info.HasField("basecall_speed"):
                            self.Mean_Basecall_Speed = info.basecall_speed.mean_basecall_speed
                        if info.HasField("n50"):
                            self.N50 = info.n50.n50
                            self.Estimated_N50 = info.n50.estimated_n50
                        if info.HasField("protocol_run_info"):
                            self.running_kit = info.protocol_run_info.meta_info.tags[
                                "kit"
                            ].string_value
                            self.Flowcell_Type = info.protocol_run_info.meta_info.tags[
                                "flow cell"
                            ].string_value
                            if info.protocol_run_info.phase != 0:
                                self.show = True
                            else:
                                self.show = False
                except Exception as e:
                    print(f"[Position.stream_instance_activity] Inner loop error: {str(e)}")
                    time.sleep(1)
            print("[Position.stream_instance_activity] Instance monitoring stopped")
        except Exception as e:
            print(f"[Position.stream_instance_activity] Error: {str(e)}")
            print(f"[Position.stream_instance_activity] Error type: {type(e)}")
            import traceback
            traceback.print_exc()


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
