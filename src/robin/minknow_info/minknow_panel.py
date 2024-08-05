from nicegui import ui, run, app
import time
import threading


# MinKNOW API Imports
import minknow_api.manager_pb2 as manager_pb2
from minknow_api.protocol_pb2 import ProtocolPhase
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
        kit="SQK-RAD114",
        centreID="NUH",
        experiment_duration=24,
        bed_file="/home/deepseq/panel_adaptive_nogenenames_20122021_hg38.bed",
        basecall_config="dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg",
        reference="mockref.ref",
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

        # function using _stop function

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
                # print (self.reference)
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
                    self.minknowhistogram = c (self.position)

                    with ui.row():
                        with ui.card().classes("drop-shadow"):
                            ui.label().bind_text_from(
                                self, "Run_ID", backward=lambda n: f"Run ID: {n}"
                            )
                            ui.label().bind_text_from(
                                self.position, "name", backward=lambda n: f"Name: {n}"
                            )
                            ui.label().bind_text_from(
                                self.position,
                                "device_type",
                                backward=lambda n: f"Device Type: {n}",
                            )
                            ui.label().bind_text_from(
                                self.position, "host", backward=lambda n: f"Host: {n}"
                            )
                            ui.label().bind_text_from(
                                self.position,
                                "running",
                                backward=lambda n: f"Running: {n}",
                            )
                            ui.label().bind_text_from(
                                self.position, "state", backward=lambda n: f"State: {n}"
                            )
                        self.protocol_run_info = ui.column().classes(
                            "p-4 border-dashed border-2"
                        )

                    self.flowcell_info = ui.row()

                    self.yield_summary_info = ui.row().classes(
                        "p-4 border-dashed border-2"
                    )

                    self.acquisition_run_info = ui.column()

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
                                    )
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
                if info.HasField("device_info"):
                    self.Channel_Count = info.device_info.device_info.max_channel_count
                if self.acquisition_run_info:
                    if info.HasField("acquisition_run_info"):
                        self.acquisition_run_info.clear()
                        with self.acquisition_run_info:
                            for field in info.acquisition_run_info.DESCRIPTOR.fields:
                                ui.label(
                                    f"{field.name}: {getattr(info.acquisition_run_info, field.name)}"
                                )
                        self.start_time = datetime.fromtimestamp(
                            info.acquisition_run_info.start_time.seconds
                            + info.acquisition_run_info.start_time.nanos / 1e9
                        )
                if info.HasField("basecall_speed"):
                    self.Mean_Basecall_Speed = info.basecall_speed.mean_basecall_speed
                if info.HasField("n50"):
                    self.N50 = info.n50.n50
                    self.Estimated_N50 = info.n50.estimated_n50
                if self.protocol_run_info:
                    if info.HasField("protocol_run_info"):
                        self.protocol_run_info.clear()
                        with self.protocol_run_info:
                            self.Experiment_Group = (
                                info.protocol_run_info.user_info.protocol_group_id.value
                            )
                            ui.label(f"Experiment Group: {self.Experiment_Group}")
                            self.Sample_ID = (
                                info.protocol_run_info.user_info.sample_id.value
                            )
                            ui.label(f"Sample ID: {self.Sample_ID}")
                            self.Current_Output_Directory = (
                                info.protocol_run_info.output_path
                            )
                            ui.label(
                                f"Current Output Directory: {self.Current_Output_Directory}"
                            )
                            self.Kit = info.protocol_run_info.meta_info.tags[
                                "kit"
                            ].string_value
                            ui.label(f"Kit: {self.Kit}")
                            self.Phase = ProtocolPhase.Name(
                                info.protocol_run_info.phase
                            )
                            ui.label(f"Phase: {self.Phase}")
                            self.Start_Time = datetime.fromtimestamp(
                                info.protocol_run_info.start_time.seconds
                                + info.protocol_run_info.start_time.nanos / 1e9
                            )
                            ui.label(f"Start Time: {self.Start_Time}")
                if info.HasField("yield_summary"):
                    if self.yield_summary_info:
                        self.yield_summary_info.clear()
                        with self.yield_summary_info:
                            self.Read_Count = info.yield_summary.read_count
                            ui.label(f"Read Count: {self.Read_Count}")
                            self.Percent_Basecalled = (
                                info.yield_summary.fraction_basecalled
                            )
                            ui.label(f"Percent Basecalled: {self.Percent_Basecalled}")
                            self.Pass_Read_Count = (
                                info.yield_summary.basecalled_pass_read_count
                            )
                            ui.label(f"Pass Read Count: {self.Pass_Read_Count}")
                            self.Fail_Read_Count = (
                                info.yield_summary.basecalled_fail_read_count
                            )
                            ui.label(f"Fail Read Count: {self.Fail_Read_Count}")
                            self.Pass_Bases = info.yield_summary.basecalled_pass_bases
                            ui.label(f"Pass Bases: {self.Pass_Bases}")
                            self.Fail_Bases = info.yield_summary.basecalled_fail_bases
                            ui.label(f"Fail Bases: {self.Fail_Bases}")

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
