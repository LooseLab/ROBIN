from nicegui import binding, ui, run, app
import time
import threading
import logging

# Configure logger
logger = logging.getLogger(__name__)

# Create console handler with formatting
console_handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

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

from datetime import datetime, timedelta


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
        logger.debug(f"Setting up device for position: {self.selected_position}")
        self.position_choose.close()
        
        try:
            logger.debug("Creating Position instance")
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
                logger.debug("Position instance created, calling setup")
                self.minKNOW_display.setup(self.holder)
                logger.debug("Setup completed")
        except Exception as e:
            logger.error(f"Error setting up device: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            logger.error(traceback.format_exc())
            ui.notify(f"Error setting up device: {str(e)}", type="negative")

    def enable_access(self, e):
        self.access_device.visible = True
        # self.position_choose.open()

    async def connect_to_localhost(self):
        logger.debug("Attempting to connect to localhost...")
        try:
            self.connection_ip = "127.0.0.1"
            logger.debug("Set connection_ip to 127.0.0.1")
            await self.add_positions()
            logger.debug("Successfully connected to localhost")
            return True
        except Exception as e:
            logger.error(f"Error connecting to localhost: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            logger.error(traceback.format_exc())
            return False

    async def add_positions(self) -> None:
        logger.debug("Starting add_positions...")
        with disable(self.connect_now):
            ui.notify("Trying to connect to MinKNOW")
            try:
                logger.debug("Calling _connect_positions...")
                await run.io_bound(self._connect_positions)
                logger.debug("_connect_positions completed")
                
                if self.connected:
                    logger.info("Connection successful")
                    ui.notify("Connection Successful.", type="positive")
                    self.connectiondialog.close()
                    
                    logger.debug("Creating radio buttons for positions")
                    with self.choices:
                        ui.radio(
                            {item: str(item) for item in self.positions},
                            on_change=self.enable_access,
                        ).bind_value(self, "selected_position").props("inline")
                    logger.debug("Opening position chooser")
                    self.position_choose.open()
                else:
                    logger.warning("Connection failed")
                    self.connection_ip = None
                    ui.notify(
                        "Unable to connect.\nPlease check MinKNOW is running on this IP.",
                        type="warning",
                    )
            except Exception as e:
                logger.error(f"Error in add_positions: {str(e)}")
                logger.error(f"Error type: {type(e)}")
                logger.error(traceback.format_exc())
                ui.notify(f"Connection error: {str(e)}", type="negative")
            finally:
                logger.debug("add_positions finished")
            ui.notify("Finished trying to connect to MinKNOW")

    def _connect_positions(self):
        logger.debug("Starting MinKNOWFish connection...")
        try:
            logger.debug(f"Attempting to connect to MinKNOW at {self.connection_ip}")
            self.manager = Manager(host=self.connection_ip)
            logger.debug("Manager connection established")
            
            try:
                logger.debug("Getting version info...")
                version_info = self.manager.version
                logger.debug(f"MinKNOW version info: {version_info}")
            except Exception as ve:
                logger.warning(f"Could not get version info: {str(ve)}")
                logger.debug("Continuing anyway...")
            
            logger.debug("Getting flow cell positions...")
            self.positions = list(self.manager.flow_cell_positions())
            logger.debug(f"Found positions: {self.positions}")
            self.connected = True
            
        except Exception as e:
            logger.error(f"MinKNOWFish connection error: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            logger.error(traceback.format_exc())
            self.connected = False
            raise

    def watch_positions(self):
        """Watch for new flow cell positions and update UI accordingly"""
        logger.debug("Starting position watcher...")
        try:
            while True:
                if self.manager:
                    logger.debug("Watching for position changes...")
                    for item in self.manager.rpc.watch_flow_cell_positions():
                        if self.tabs:
                            logger.debug(f"Processing position changes: {item}")
                            with self.tabs:
                                for position in item.additions:
                                    ui.notify(f"{position.name} connected.", type="positive")
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
                    logger.debug("No manager available, waiting...")
                    time.sleep(1)
        except Exception as e:
            logger.error(f"Error in watch_positions: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            logger.error(traceback.format_exc())

    def connect_to_position(self, position):
        """Connect to a specific position when requested"""
        logger.debug(f"Connecting to position: {position.name}")
        try:
            self.selected_position = position
            self.choose_device()
        except Exception as e:
            logger.error(f"Error connecting to position: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            logger.error(traceback.format_exc())
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
        logger.debug("Starting position initialization...")
        # Call parent class init with configuration parameters
        logger.debug("Calling parent class init")
        super().__init__(
            kit=kit,
            centreID=centreID,
            experiment_duration=experiment_duration,
            bed_file=bed_file,
            basecall_config=basecall_config,
            reference=reference,
        )
        logger.debug("Parent class init completed")
        
        self.ip = ip
        self.yield_data = []  # Initialize yield_data
        self.threads = []  # Keep track of background threads
        self._initialized = False  # Track initialization state
        
        # Add state variables for UI updates
        self._current_state = None
        self._current_run_info = None
        self._needs_ui_update = False
        
        try:
            logger.debug(f"Using position object: {position}")
            logger.debug(f"Position type: {type(position)}")
            logger.debug(f"Position attributes: {dir(position)}")
            
            # Use the position object directly
            self.position = position
            logger.debug("Attempting to connect to position...")
            self.connection = self.position.connect()
            logger.debug("Position connection established")
            
            # Initialize other attributes
            logger.debug("Initializing attributes")
            self.Run_ID = None
            self.yield_summary_info = None
            self.acquisition_run_info = None
            self.protocol_run_info = None
            self.flowcell_info = None
            self.minknow_info_pane = None
            self.watchfolder = None
            self.shutdown_event = threading.Event()
            self._initialized = True
            logger.debug("Position initialization completed successfully")
        except Exception as e:
            logger.error(f"Error: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            import traceback
            logger.error(traceback.format_exc())
            raise

    def setup(self, container=None):
        """Setup UI elements in the correct context"""
        logger.debug("Starting setup...")
        if not self._initialized:
            logger.error("Error: Position not properly initialized")
            raise RuntimeError("Position not properly initialized")
        
        try:
            logger.debug("Creating UI elements")
            # Create timer in the global UI context
            self.timer = ui.timer(1, self._worker)

            # Create UI elements in the provided container or global context
            target = container if container is not None else ui
            logger.debug(f"Using container: {target}")
            
            logger.debug("Creating UI cards")
            with target.card().classes("w-full"):
                # Device Header
                with ui.card().classes("w-full p-4 mb-4 bg-gray-50"):
                    with ui.row().classes("items-center"):
                        # Device Icon
                        with ui.avatar(size='xl').classes("mr-4"):
                            ui.icon('biotech', size='xl').classes("text-gray-600")
                        # Device Info
                        with ui.column().classes("flex-grow"):
                            with ui.row().classes("items-baseline"):
                                ui.label(f"MinION - {self.position.name}").classes("text-xl font-bold text-gray-900")
                                ui.label("MinKNOW Monitoring").classes("ml-2 text-sm text-gray-600")
                            ui.label(f"Device Type: {self.position.device_type}").classes("text-sm text-gray-700")
                            # Add Sample ID and Experiment ID
                            self.sample_id_label = ui.label("Sample ID: -").classes("text-lg font-semibold text-blue-700 mt-2")
                            self.experiment_id_label = ui.label("Experiment ID: -").classes("text-lg font-semibold text-blue-700")
                            with ui.row().classes("mt-1"):
                                ui.badge("RUNNING", color="positive" if self.position.running else "negative")
                                ui.badge(str(self.position.state), color="gray")

                # Status and Hardware Information
                with ui.card().classes("w-full p-4 mb-4 bg-gray-50"):
                    ui.label("Device Status").classes("text-lg font-bold mb-3 text-gray-900")
                    with ui.grid().classes("grid-cols-3 gap-4"):
                        # Device Info
                        with ui.column().classes("space-y-2"):
                            ui.label("Device Information").classes("text-md font-semibold mb-2 text-gray-800")
                            self.status_label = ui.label("Status: Initializing...").classes("text-sm text-gray-700")
                            ui.label().bind_text_from(
                                self.position, "device_type",
                                backward=lambda n: f"Device Type: {n}"
                            ).classes("text-sm text-gray-700")
                            ui.label().bind_text_from(
                                self.position, "host",
                                backward=lambda n: f"Host: {n}"
                            ).classes("text-sm text-gray-700")
                            self.run_id_label = ui.label("Run ID: -").classes("text-sm text-gray-700")
                            self.run_state_label = ui.label("Run State: -").classes("text-sm text-gray-700")
                            self.run_start_time_label = ui.label("Start Time: -").classes("text-sm text-gray-700")
                            self.acquisition_duration_label = ui.label("Run Duration: -").classes("text-sm text-gray-700")
                            self.expected_duration_label = ui.label("Expected Duration: -").classes("text-sm text-gray-700")
                            self.expected_finish_label = ui.label("Expected Finish: -").classes("text-sm text-gray-700")
                            self.run_purpose_label = ui.label("Purpose: -").classes("text-sm text-gray-700")
                            
                        # Flowcell Hardware
                        with ui.column().classes("space-y-2"):
                            ui.label("Flowcell Hardware").classes("text-md font-semibold mb-2 text-gray-800")
                            self.flowcell_label = ui.label("Flowcell: -").classes("text-sm text-gray-700")
                            self.flowcell_id_label = ui.label("Flow Cell ID: -").classes("text-sm text-gray-700")
                            self.flowcell_type_label = ui.label("Type: -").classes("text-sm text-gray-700")
                            self.asic_label = ui.label("ASIC: -").classes("text-sm text-gray-700")
                            self.channel_label = ui.label("Channels: -").classes("text-sm text-gray-700")
                            self.user_specified_product_code_label = ui.label("User Product Code: -").classes("text-sm text-gray-700")
                            self.user_specified_flow_cell_id_label = ui.label("User Flow Cell ID: -").classes("text-sm text-gray-700")
                            
                        # Flowcell Status
                        with ui.column().classes("space-y-2"):
                            ui.label("Flowcell Status").classes("text-md font-semibold mb-2 text-gray-800")
                            self.adapter_status = ui.label("Adapter: -").classes("text-sm text-gray-700")
                            self.temp_label = ui.label("Temperature: -").classes("text-sm text-gray-700")
                            self.usage_label = ui.label("Usage Count: -").classes("text-sm text-gray-700")
                            self.barcode_label = ui.label("Barcode Kit: -").classes("text-sm text-gray-700")
                            self.has_adapter_label = ui.label("Has Adapter: -").classes("text-sm text-gray-700")
                            self.adapter_id_label = ui.label("Adapter ID: -").classes("text-sm text-gray-700")
                            self.adapter_insertion_time_label = ui.label("Adapter Inserted: -").classes("text-sm text-gray-700")

                # Yield Statistics
                with ui.card().classes("w-full p-4 mb-4 bg-white"):
                    ui.label("Yield Statistics").classes("text-lg font-bold mb-3 text-gray-900")
                    with ui.grid().classes("grid-cols-3 gap-6"):
                        # Read Statistics
                        with ui.card().classes("p-4 bg-gray-50 rounded-lg"):
                            ui.label("Read Statistics").classes("text-md font-semibold mb-2 text-gray-800")
                            with ui.column().classes("space-y-2"):
                                self.read_count_label = ui.label("Reads: 0").classes("text-sm text-gray-700")
                                self.bases_label = ui.label("Bases: 0").classes("text-sm text-gray-700")
                                self.n50_label = ui.label("N50: -").classes("text-sm text-gray-700")
                                self.estimated_n50_label = ui.label("Estimated N50: -").classes("text-sm text-gray-700")
                        
                        # Quality Metrics
                        with ui.card().classes("p-4 bg-gray-50 rounded-lg"):
                            ui.label("Quality Metrics").classes("text-md font-semibold mb-2 text-gray-800")
                            with ui.column().classes("space-y-2"):
                                self.pass_reads_label = ui.label("Pass Reads: 0").classes("text-sm text-gray-700")
                                self.fail_reads_label = ui.label("Fail Reads: 0").classes("text-sm text-gray-700")
                                self.mean_basecall_speed_label = ui.label("Basecall Speed: -").classes("text-sm text-gray-700")
                        
                        # Output Information
                        with ui.card().classes("p-4 bg-gray-50 rounded-lg"):
                            ui.label("Output Settings").classes("text-md font-semibold mb-2 text-gray-800")
                            with ui.column().classes("space-y-2"):
                                self.output_folder_label = ui.label("Output: -").classes("text-sm text-gray-700 break-all")
                                self.basecall_config_label = ui.label("Basecall Config: -").classes("text-sm text-gray-700 break-all")

                # MinKNOW Details
                with ui.card().classes("w-full p-4 bg-white"):
                    logger.debug("Creating Minknow_Info instance")
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
                        # Protocol Information
                        self.protocol_run_info = ui.column().classes(
                            "w-full p-4 bg-gray-50 rounded-lg mt-4"
                        )
                        # Flowcell Information
                        self.flowcell_info = ui.column().classes(
                            "w-full p-4 bg-gray-50 rounded-lg mt-4"
                        )
                        # Yield Summary
                        self.yield_summary_info = ui.column().classes(
                            "w-full p-4 bg-gray-50 rounded-lg mt-4"
                        )
                        # Acquisition Information
                        self.acquisition_run_info = ui.column().classes(
                            "w-full p-4 bg-gray-50 rounded-lg mt-4"
                        )
                        
                        logger.debug("Creating MinknowHistograms instance")
                        self.minknowhistogram = MinknowHistograms(self.position)

            # Start background tasks
            logger.debug("Starting background tasks")
            self._start_background_tasks()
            logger.debug("Setup completed successfully")
                
        except Exception as e:
            logger.error(f"Error in setup: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            import traceback
            logger.error(traceback.format_exc())
            raise

    def _start_background_tasks(self):
        """Start background tasks with proper error handling"""
        logger.debug("Starting background tasks...")
        if not self._initialized:
            logger.error("Error: Position not properly initialized")
            raise RuntimeError("Position not properly initialized")
            
        tasks = [
            (self.watch_flowcell, "Flowcell Monitor"),
            (self.stream_current_run, "Run Monitor"),
            (self.stream_instance_activity, "Activity Monitor")
        ]
        
        for task_func, task_name in tasks:
            try:
                logger.debug(f"Starting {task_name}")
                thread = threading.Thread(
                    target=self._run_task_with_error_handling,
                    args=(task_func, task_name),
                    daemon=True
                )
                thread.start()
                self.threads.append(thread)
                logger.debug(f"Started {task_name} thread")
            except Exception as e:
                logger.error(f"Error starting {task_name}: {str(e)}")
                logger.error(f"Error type: {type(e)}")
                import traceback
                logger.error(traceback.format_exc())
                
    def _run_task_with_error_handling(self, task_func, task_name):
        """Wrapper to run a task with error handling"""
        logger.debug(f"Starting {task_name}")
        if not self._initialized:
            logger.warning(f"Attempting to run {task_name} before initialization")
            return
            
        try:
            logger.debug(f"Executing {task_name}")
            task_func()
            logger.debug(f"{task_name} completed")
        except Exception as e:
            logger.error(f"Error in {task_name}: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            import traceback
            logger.error(traceback.format_exc())

    def stop(self):
        """Stop all background tasks"""
        logger.debug("Stopping position tasks...")
        self.shutdown_event.set()
        for thread in self.threads:
            thread.join(timeout=1.0)
        logger.debug("Position tasks stopped")

    def _worker(self):
        """Timer callback to update UI elements"""
        logger.debug("Starting worker update")
        try:
            if not self._initialized:
                logger.warning("Warning: Position not initialized")
                return
                
            if not hasattr(self, 'status_label'):
                logger.warning("Warning: UI elements not initialized")
                return
            
            # Get current acquisition info
            try:
                # Check if acquisition is running first
                current_status = self.connection.acquisition.current_status()
                logger.debug(f"Current status: {current_status}")
                # MinKNOW status constants from acquisition_pb2
                if current_status.status == 2:  # PROCESSING is 2 in MinKNOW API
                    # Only try to get acquisition info if we're actually processing
                    current_acquisition = self.connection.acquisition.get_current_acquisition_run()
                    if current_acquisition and hasattr(current_acquisition, 'run_info'):
                        run_info = current_acquisition.run_info
                        
                        # Get progress information
                        progress = self.connection.acquisition.get_progress()
                        if progress:
                            logger.debug(f"Progress info: {progress}")
                            
                            # Update read statistics from progress info
                            if hasattr(progress, 'raw_per_channel'):
                                self.read_count_label.text = f"Reads: {progress.raw_per_channel.acquired:,}"
                            
                            if hasattr(progress, 'basecalled_per_channel'):
                                basecalled = progress.basecalled_per_channel
                                if hasattr(basecalled, 'pass_'):
                                    self.pass_reads_label.text = f"Pass Reads: {basecalled.pass_:,}"
                                if hasattr(basecalled, 'fail'):
                                    self.fail_reads_label.text = f"Fail Reads: {basecalled.fail:,}"
                            
                            # Update bases count and rate
                            if hasattr(progress, 'basecalled_bases'):
                                base_text = f"Bases: {progress.basecalled_bases:,}"
                                if hasattr(progress, 'estimated_bases_per_second'):
                                    base_text += f" ({progress.estimated_bases_per_second:.1f} bases/s)"
                                self.bases_label.text = base_text
                # Don't update sample ID and experiment ID here - let stream_instance_activity handle it
            except Exception as e:
                if "No acquisition running" not in str(e):
                    logger.error(f"Error getting acquisition info: {str(e)}")
            
            # Update run information if needed
            if self._needs_ui_update and self._current_run_info:
                info = self._current_run_info
                state_name = manager_pb2.FlowCellPosition.State.Name(info.state)
                self.run_state_label.text = f"Run State: {state_name}"
                
                if info.state == manager_pb2.FlowCellPosition.State.STATE_RUNNING:
                    # Calculate run duration and update start time if start time is available
                    if hasattr(info, 'start_time'):
                        start_time = datetime.fromtimestamp(
                            info.start_time.seconds +
                            info.start_time.nanos / 1e9
                        )
                        # Update start time label with formatted datetime
                        self.run_start_time_label.text = f"Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
                        
                        # Calculate and update duration
                        duration = datetime.now() - start_time
                        hours = duration.total_seconds() / 3600
                        self.acquisition_duration_label.text = f"Run Duration: {hours:.1f} hours"
                        
                        # Get protocol run info for expected duration
                        try:
                            # Get current acquisition info
                            current_acquisition = self.connection.acquisition.get_current_acquisition_run()
                            if current_acquisition:
                                logger.debug(f"Current acquisition: {current_acquisition}")
                                
                                # Access config_summary directly from AcquisitionRunInfo
                                if hasattr(current_acquisition, 'config_summary'):
                                    config = current_acquisition.config_summary
                                    logger.debug(f"Acquisition config summary: {config}")
                                    
                                    if hasattr(config, 'target_run_until'):
                                        target = config.target_run_until
                                        logger.debug(f"Target run until criteria: {target}")
                                        
                                        if hasattr(target, 'run_time'):
                                            # run_time is in seconds
                                            expected_hours = target.run_time / 3600
                                            self.expected_duration_label.text = f"Expected Duration: {expected_hours:.1f} hours"
                                            
                                            # Calculate and display expected finish time
                                            finish_time = start_time + timedelta(hours=expected_hours)
                                            self.expected_finish_label.text = f"Expected Finish: {finish_time.strftime('%Y-%m-%d %H:%M:%S')}"
                                        else:
                                            self.expected_duration_label.text = "Expected Duration: Not Set (No run_time)"
                                            self.expected_finish_label.text = "Expected Finish: Not Available"
                                    else:
                                        self.expected_duration_label.text = "Expected Duration: Not Set (No target)"
                                        self.expected_finish_label.text = "Expected Finish: Not Available"
                                else:
                                    self.expected_duration_label.text = "Expected Duration: Not Available (No config)"
                                    self.expected_finish_label.text = "Expected Finish: Not Available"
                        except Exception as e:
                            logger.error(f"Error getting acquisition run info: {str(e)}")
                            logger.error(f"Error type: {type(e)}")
                            import traceback
                            logger.error(traceback.format_exc())
                            self.expected_duration_label.text = "Expected Duration: Error"
                            self.expected_finish_label.text = "Expected Finish: Error"
                    
                    # Update acquisition info
                    if hasattr(info, 'config_summary'):
                        config = info.config_summary
                        config_info = []
                        
                        if hasattr(config, 'basecalling_enabled'):
                            config_info.append(f"Basecalling: {'Enabled' if config.basecalling_enabled else 'Disabled'}")
                        
                        if hasattr(config, 'reads_directory'):
                            config_info.append(f"Output: {config.reads_directory}")
                            
                        if hasattr(config, 'basecalling_config_filename'):
                            config_info.append(f"Config: {config.basecalling_config_filename}")
                            
                        if config_info:
                            self.run_purpose_label.text = " | ".join(config_info)
                        else:
                            self.run_purpose_label.text = "No config information available"
                else:
                    self.acquisition_duration_label.text = "Run Duration: -"
                    self.run_start_time_label.text = "Start Time: -"
                    self.expected_duration_label.text = "Expected Duration: -"
                    self.expected_finish_label.text = "Expected Finish: -"
                    self.run_purpose_label.text = "No active run"
                
                self._needs_ui_update = False
                
            # Update status
            try:
                state = self.position.state
                self.status_label.text = f"Status: {state}"
            except Exception as e:
                logger.error(f"Error updating status: {str(e)}")
                
            # Update run ID
            try:
                if hasattr(self, 'Run_ID') and self.Run_ID:
                    self.run_id_label.text = f"Run ID: {self.Run_ID}"
            except Exception as e:
                logger.error(f"Error updating run ID: {str(e)}")
                
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
                logger.error(f"Error updating yield stats: {str(e)}")
                logger.error(f"Yield info type: {type(self.yield_summary_info)}")
                logger.error(f"Yield info: {self.yield_summary_info}")
                import traceback
                logger.error(traceback.format_exc())
                
        except Exception as e:
            logger.error(f"Error in worker: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            import traceback
            logger.error(traceback.format_exc())

    def stopped(self):
        """Check if the position has been stopped"""
        return self.shutdown_event.is_set()

    def watch_flowcell(self):
        """Monitor flowcell information"""
        logger.debug("Starting flowcell monitoring")
        try:
            while not self.stopped():
                if self.connection:
                    for info in self.connection.device.stream_flow_cell_info():
                        if self.stopped():
                            break
                            
                        # Update all flowcell information
                        if info.has_flow_cell:
                            # Update flowcell label with connection status
                            self.flowcell_label.text = "Flowcell: Connected"
                            
                            # Update hardware information
                            self.flowcell_id_label.text = f"Flow Cell ID: {info.flow_cell_id if info.flow_cell_id else 'Not Set'}"
                            self.flowcell_type_label.text = f"Type: {info.product_code if info.product_code else 'Standard'}"
                            self.asic_label.text = f"ASIC: {info.asic_version} ({info.asic_id_str[:8]}...)"
                            self.channel_label.text = f"Channels: {info.channel_count:,} ({info.wells_per_channel} wells/ch)"
                            
                            # Update user-specified information
                            user_product = info.user_specified_product_code if info.user_specified_product_code else "Not Set"
                            self.user_specified_product_code_label.text = f"User Product Code: {user_product}"
                            
                            user_flow_cell = info.user_specified_flow_cell_id if info.user_specified_flow_cell_id else "Not Set"
                            self.user_specified_flow_cell_id_label.text = f"User Flow Cell ID: {user_flow_cell}"
                            
                            # Update adapter information
                            self.has_adapter_label.text = f"Has Adapter: {'Yes' if info.has_adapter else 'No'}"
                            if info.has_adapter:
                                self.adapter_status.text = "Adapter: Connected"
                                self.adapter_status.classes("text-sm text-green-700")
                                self.adapter_id_label.text = f"Adapter ID: {info.adapter_id}"
                                if hasattr(info, 'adapter_insertion_time'):
                                    insertion_time = datetime.fromtimestamp(
                                        info.adapter_insertion_time.seconds +
                                        info.adapter_insertion_time.nanos / 1e9
                                    )
                                    self.adapter_insertion_time_label.text = f"Adapter Inserted: {insertion_time.strftime('%Y-%m-%d %H:%M:%S')}"
                            else:
                                self.adapter_status.text = "Adapter: Not Connected"
                                self.adapter_status.classes("text-sm text-red-700")
                                self.adapter_id_label.text = "Adapter ID: -"
                                self.adapter_insertion_time_label.text = "Adapter Inserted: -"
                            
                            # Update temperature and usage information
                            if info.temperature_offset:
                                temp_offset = round(info.temperature_offset, 1)
                                self.temp_label.text = f"Temperature Offset: {temp_offset}°C"
                            else:
                                self.temp_label.text = "Temperature Offset: Not Available"
                            
                            self.usage_label.text = f"Usage Count: {info.use_count}"
                            
                            # Update barcode information
                            if info.barcode_kit:
                                self.barcode_label.text = f"Barcode Kit: {info.barcode_kit}"
                                if info.barcodes:
                                    barcodes_str = ', '.join(info.barcodes)
                                    self.barcode_label.text += f" ({barcodes_str})"
                            else:
                                self.barcode_label.text = "Barcode Kit: None"
                        else:
                            # Clear all labels when no flowcell is present
                            self.flowcell_label.text = "Flowcell: Not Connected"
                            self.flowcell_id_label.text = "Flow Cell ID: -"
                            self.flowcell_type_label.text = "Type: -"
                            self.asic_label.text = "ASIC: -"
                            self.channel_label.text = "Channels: -"
                            self.user_specified_product_code_label.text = "User Product Code: -"
                            self.user_specified_flow_cell_id_label.text = "User Flow Cell ID: -"
                            self.adapter_status.text = "Adapter: -"
                            self.has_adapter_label.text = "Has Adapter: -"
                            self.adapter_id_label.text = "Adapter ID: -"
                            self.adapter_insertion_time_label.text = "Adapter Inserted: -"
                            self.temp_label.text = "Temperature: -"
                            self.usage_label.text = "Usage Count: -"
                            self.barcode_label.text = "Barcode Kit: -"
                            
                            # Reset adapter status color
                            self.adapter_status.classes("text-sm text-gray-700")
                else:
                    time.sleep(1)
            logger.debug("Flowcell monitoring stopped")
        except Exception as e:
            logger.error(f"Error in watch_flowcell: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            import traceback
            logger.error(traceback.format_exc())

    def stream_current_run(self) -> None:
        """Monitor current run information"""
        logger.debug("Starting run monitoring")
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

                        # Update run state and ID
                        self.Run_ID = info.run_id
                        
                        # Store current run info and flag for UI update
                        self._current_run_info = info
                        self._needs_ui_update = True
                        
                        # Update run status
                        if info.state == manager_pb2.FlowCellPosition.State.STATE_RUNNING:
                            self.Live_Run = True
                            if first_run:
                                logger.info(f"Run {self.Run_ID} seen on {self.position}")
                        else:
                            self.Live_Run = False
                            self.Run_ID = None
                            if first_run:
                                logger.info(f"Connected to {self.position}.")
                            else:
                                logger.warning(f"Run {self.Run_ID} finished on {self.position}.")
                        
                        first_run = False
                except Exception as e:
                    logger.error(f"Inner loop error: {str(e)}")
                    time.sleep(1)
            logger.debug("Run monitoring stopped")
        except Exception as e:
            logger.error(f"Error: {str(e)}")
            logger.error(f"Error type: {type(e)}")
            import traceback
            logger.error(traceback.format_exc())

    def stream_instance_activity(self) -> None:
        """Monitor instance activity"""
        logger.debug("Starting instance monitoring")
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
                            # Log the entire protocol run info for debugging
                            logger.debug(f"Protocol run info: {info.protocol_run_info}")
                            
                            # Get sample ID and experiment ID from user_info
                            if hasattr(info.protocol_run_info, 'user_info'):
                                user_info = info.protocol_run_info.user_info
                                logger.debug(f"User info: {user_info}")
                                
                                # Update sample ID
                                if hasattr(user_info, 'sample_id'):
                                    sample_id = user_info.sample_id.value if hasattr(user_info.sample_id, 'value') else 'Not Set'
                                    self.sample_id_label.text = f"Sample ID: {sample_id}"
                                else:
                                    self.sample_id_label.text = "Sample ID: Not Set"
                                
                                # Update protocol group ID (experiment ID)
                                if hasattr(user_info, 'protocol_group_id'):
                                    protocol_group_id = user_info.protocol_group_id.value if hasattr(user_info.protocol_group_id, 'value') else 'Not Set'
                                    self.experiment_id_label.text = f"Experiment ID: {protocol_group_id}"
                                else:
                                    self.experiment_id_label.text = "Experiment ID: Not Set"
                            
                            # Update other existing fields from meta_info
                            meta_info = info.protocol_run_info.meta_info
                            if hasattr(meta_info, 'tags'):
                                if 'kit' in meta_info.tags:
                                    self.running_kit = meta_info.tags['kit'].string_value
                                if 'flow cell' in meta_info.tags:
                                    self.Flowcell_Type = meta_info.tags['flow cell'].string_value
                            
                            if info.protocol_run_info.phase != 0:
                                self.show = True
                            else:
                                self.show = False
                except Exception as e:
                    logger.error(f"[Position.stream_instance_activity] Inner loop error: {str(e)}")
                    time.sleep(1)
            logger.debug("Instance monitoring stopped")
        except Exception as e:
            logger.error(f"[Position.stream_instance_activity] Error: {str(e)}")
            logger.error(f"[Position.stream_instance_activity] Error type: {type(e)}")
            import traceback
            logger.error(traceback.format_exc())


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
        logging.info("GUI launched by auto-reload")

    run_class(port=12398, reload=True)
