# Python imports.
from __future__ import annotations
from nicegui import ui, app
from robin import images
from pathlib import Path
import threading
from datetime import datetime
from robin import theme
from minknow_api.manager import Manager
import os
import asyncio
import re
from contextlib import contextmanager

# MinKNOW API Imports
from robin.utilities.camera import Camera
from typing import Sequence
import uuid

import cv2
import zxingcpp

import numpy as np


import base64

# We need `find_protocol` to search for the required protocol given a kit + product code.
from minknow_api.tools import protocols

UNIQUE_ID: str = str(uuid.uuid4())


class ExperimentSpec(object):
    def __init__(self, position):
        self.position = position
        self.protocol_id = ""


ExperimentSpecs = Sequence[ExperimentSpec]


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


# Define the path to the image file used in the header and footer
IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
)


DEVICEDICT = {
    "FLONGLE": "Flongle",
    "GRIDION": "GridION",
    "MINION": "MinION",
    "PROMETHION": "PromethION",
    "P2_INTEGRATED": "P2Integrated",
    "P2_SOLO": "P2Solo",
}


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


def identify_device(device_type, device_name):
    if device_type == "MINION":
        return f"{DEVICEDICT[device_type]}", "MinION-flow-cell-lights-on-whiteBG.jpg"
    elif device_type == "GRIDION":
        return (
            f"{DEVICEDICT[device_type]}",
            "GridION_Front_Square_Elevated_Closed_Flow Cells 2_White.jpg",
        )
    elif device_type == "P2_INTEGRATED":
        return (
            f"{DEVICEDICT[device_type]}",
            "p2_-left-45_screen-up-1_splash_transparent.png",
        )
    elif device_type == "P2_SOLO":
        if device_name.startswith("P2S"):
            return "P2 Solo", "P2_Solo_Left-45_Open_Full.png"
        else:
            return (f"{device_type}", None)
    elif device_type == "PROMETHION":
        return (
            f"{DEVICEDICT[device_type]}",
            "PromethION_24_Right 45 Elevated_Closed_Full_white.jpg",
        )
    else:
        return (f"{device_type}", None)


class Minknow_Info:
    def __init__(
        self,
        position,
        centreID,
        kit,
        reference,
        basecall_config,
        bed_file,
        experiment_duration,
        dev=False,
    ):
        self.dev = dev
        self.color = "text-blue"
        self.position = position
        self.connection = self.position.connect()
        self.device = self.position.device_type
        self.show = False

        self.basecall_config = basecall_config
        self.centreID = centreID
        self.kit = kit
        self.reference = reference
        self.bed_file = bed_file
        self.experiment_duration = experiment_duration

        self.name, self.image = identify_device(
            self.position.device_type, self.position.name
        )
        if self.image:
            self.deviceicon = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)),
                "ONTimages",
                self.image,
            )
        else:
            self.deviceicon = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)),
                "ONTimages",
                "unknown_sequencer.png",
            )
        self.check_instance = threading.Thread(
            target=self.stream_instance_activity, args=()
        )
        self.check_instance.daemon = True
        self.check_instance.start()
        self.render_me()

    def render_me(self):
        # Initialize camera handlers
        sample_camera = Camera(
            icon="photo_camera",
            icon_color="primary",
            background_color="rgba(66, 165, 245, 0.3)",
            for_id="samplefileupload",
            canvas_id="samplecanvas",
            on_change=lambda: generate_sampleID(),
        )

        flowcell_camera = Camera(
            icon="photo_camera",
            icon_color="primary",
            background_color="rgba(66, 165, 245, 0.3)",
            for_id="flowcellfileupload",
            canvas_id="flowcellcanvas",
            on_change=lambda: generate_flowcellID(),
        )

        # Create containers for both panels with full width
        with ui.column().classes('w-full h-full'):
            # Setup Panel - only visible when not running
            setup_container = ui.element('div').classes('w-full h-full').bind_visibility_from(self, 'show', value=False)
            monitor_container = ui.element('div').classes('w-full h-full').bind_visibility_from(self, 'show')

            with setup_container:
                with ui.card().classes('w-full h-full'):
                    with ui.card_section():
                        ui.label("Sample Setup").classes('text-h5 text-primary q-mb-md')
                    
                    with ui.card_section().classes('q-pa-md'):
                        with ui.stepper().props('vertical').classes('w-full') as stepper:
                            # Sample ID Step
                            with ui.step('Sample ID').classes('text-body1 text-weight-medium'):
                                ui.label('Sample ID Guidelines').classes('text-subtitle1 text-weight-medium q-mb-sm')
                                ui.label('IDs should not include human identifiable information').classes('text-caption q-mb-md')
                                
                                with ui.row().classes('items-center w-full q-gutter-md'):
                                    sampleid = ui.input(
                                        placeholder='Enter sample ID',
                                        validation={
                                            "Too short": lambda value: len(value) >= 5,
                                            "No Flowcell IDs": lambda value: re.match(r"^[A-Za-z]{3}\d{5}$", value) is None,
                                            "Alphanumeric only": lambda value: value.isalnum() and not any(char.isspace() for char in value),
                                        }
                                    ).props('outlined dense').classes('w-full')
                                    
                                    sample_camera.show_camera()
                                    
                                with ui.stepper_navigation():
                                    checker = ErrorChecker(sampleid)
                                    ui.button('Next', on_click=stepper.next).bind_enabled_from(checker, 'no_errors')

                            # Flowcell ID Step
                            with ui.step('Flowcell').classes('text-body1 text-weight-medium'):
                                ui.label('Flowcell ID Entry').classes('text-subtitle1 text-weight-medium q-mb-sm')
                                
                                with ui.row().classes('items-center w-full q-gutter-md'):
                                    flowcellid = ui.input(
                                        placeholder='Enter flowcell ID',
                                        validation={
                                            "Not a valid flowcell ID": lambda value: re.match(
                                                r"^[A-Za-z]{3}\d{5}$", value
                                            ) is not None
                                        },
                                    ).props('outlined dense').classes('w-full')
                                    
                                    flowcell_camera.show_camera()
                                
                                with ui.stepper_navigation():
                                    ui.button('Back', on_click=stepper.previous).props('flat')
                                    checkerflowcell = ErrorChecker(flowcellid)
                                    ui.button('Next', on_click=stepper.next).bind_enabled_from(
                                        checkerflowcell, 'no_errors'
                                    )

                            # Device Position Step
                            with ui.step('Device Position').classes('text-body1 text-weight-medium'):
                                ui.label('Position Selection').classes('text-subtitle1 text-weight-medium q-mb-sm')
                                ui.label('Select the correct device position').classes('text-caption q-mb-md')
                                
                                run_button = ui.button(
                                    'Start Run',
                                    on_click=lambda: self.start_run(
                                        position=self.position.name,
                                        reference=self.reference,
                                        sample_id=sampleid.value,
                                        flowcell_id=flowcellid.value,
                                        kit=self.kit,
                                        basecall_config=self.basecall_config,
                                        centreID=self.centreID,
                                        experiment_duration=self.experiment_duration,
                                        bed_file=self.bed_file,
                                    )
                                ).props('color=primary')

                                ui.radio(
                                    [self.position.name], 
                                    value=self.position.name
                                ).on('update:model-value', lambda: run_button.enable()
                                ).classes('q-mt-md')

                                with ui.stepper_navigation():
                                    ui.button('Back', on_click=stepper.previous).props('flat')
                                    run_button.disable()

                        # Settings Summary
                        with ui.expansion('Fixed Settings', icon='settings').classes('w-full q-mt-md'):
                            with ui.card().classes('q-pa-md'):
                                self._render_setting_item('Device', self.device)
                                current_date = datetime.now()
                                self._render_setting_item(
                                    'Experiment Group ID', 
                                    f"{self.centreID}_{current_date.strftime('%B')}_{current_date.year}"
                                )
                                self._render_setting_item('Kit', self.kit)
                                self._render_setting_item('Reference', self.reference)

            # Monitor Panel - only visible when running
            with monitor_container:
                with ui.card().classes('w-full h-full'):
                    # Device Header
                    with ui.row().classes('items-center q-pa-md'):
                        with ui.avatar(size='xl').classes('q-mr-md'):
                            ui.image(self.deviceicon)
                        with ui.column():
                            ui.label(f'{self.name} - {self.position.name}').classes('text-h6')
                            ui.label('MinKNOW Monitoring').classes(f'text-caption {self.color}')
                            ui.label(f'Device Type: {self.position.device_type}').classes('text-body2')

                    # Status Indicators
                    with ui.row().classes('q-pa-md q-gutter-md justify-between'):
                        self._render_status_chip('Running', self, 'show')
                        self._render_status_chip('Basecalling', self, 'Basecall_Speed', 
                                              lambda v: 'positive' if v > 0 else 'negative')

                    # Monitoring Grid
                    with ui.grid(columns=3).classes('q-pa-md q-gutter-md').bind_visibility_from(self, 'show'):
                        # Basic Info
                        self._render_monitoring_tile('Experiment Group', 'Experiment_Group')
                        self._render_monitoring_tile('Sample ID', 'Sample_ID')
                        self._render_monitoring_tile('Flowcell Type', 'Flowcell_Type')
                        
                        # Performance Metrics
                        #self._render_monitoring_tile('Read Count', 'Read_Count', format_func=lambda n: f'{n:,}')
                        #self._render_monitoring_tile('N50', 'N50', format_func=lambda n: f'{n:,} bp')
                        #self._render_monitoring_tile('Basecall Speed', 'Mean_Basecall_Speed', 
                                                    #format_func=lambda n: f'{n:.1f} samples/s')
                        
                        # Quality Metrics
                        #self._render_monitoring_tile('Pass Reads', 'Pass_Read_Count', format_func=lambda n: f'{n:,}')
                        #self._render_monitoring_tile('Pass Bases', 'Pass_Bases', format_func=lambda n: f'{n:,}')
                        #self._render_monitoring_tile('Channel Count', 'Channel_Count')

                    # Time Information
                    with ui.row().classes('q-pa-md q-gutter-md justify-between'):
                        with ui.card().classes('col'):
                            ui.label('Start Time').classes(f'text-caption {self.color}')
                            ui.label().bind_text_from(
                                self, 'start_time', 
                                lambda t: t.strftime('%Y-%m-%d %H:%M:%S') if t else 'Not Started'
                            )

    def _render_setting_item(self, label: str, value: str):
        """Helper method to render consistent setting items"""
        with ui.row().classes('items-center justify-between w-full q-py-sm'):
            ui.label(label).classes('text-body2')
            ui.label(str(value)).classes('text-body2 text-weight-medium')

    def _render_monitoring_tile(self, label: str, bind_property: str, format_func=str):
        """Helper method to render monitoring information tiles"""
        with ui.card().classes('q-pa-sm'):
            ui.label(label).classes(f'text-caption {self.color}')
            ui.label('--').bind_text_from(
                self, bind_property,
                lambda n: format_func(n) if n is not None else '--'
            ).classes('text-body1 text-weight-medium')

    def _render_status_chip(self, label: str, obj, bind_property: str, 
                          color_func=lambda v: 'positive' if v else 'negative'):
        """Helper method to render status indicator chips"""
        status_button = ui.button(
            label,
            icon='circle'
        ).props('flat dense').classes('q-ma-xs')
        
        def update_status():
            if hasattr(obj, bind_property):
                value = getattr(obj, bind_property)
                color = color_func(value)
                status_button.props(f'color={color}')
        
        ui.timer(1, update_status)
        return status_button

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
            # print(info)
            if info.HasField("device_info"):
                self.Channel_Count = info.device_info.device_info.max_channel_count
                # print (info.device_info.device_info)
            if info.HasField("acquisition_run_info"):
                # print(info.acquisition_run_info)
                # self.acquisition_run_info.clear()
                # with self.acquisition_run_info:
                # for field in info.acquisition_run_info.DESCRIPTOR.fields:
                # print (f"{field.name}: {getattr(info.acquisition_run_info, field.name)}")
                # ui.label(f"{field.name}: {getattr(info.acquisition_run_info, field.name)}")
                ### This isn't giving us the start time for some reason?!
                # print (info.acquisition_run_info.config_summary)
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
                # self.log(f"{info.basecall_speed.mean_basecall_speed}")
                # self.log(type(info.basecall_speed.mean_basecall_speed))
                self.Mean_Basecall_Speed = info.basecall_speed.mean_basecall_speed
            # if info.HasField("flow_cell_health"):
            # if info.HasField("flow_cell_info"):
            if info.HasField("n50"):
                self.N50 = info.n50.n50
                self.Estimated_N50 = info.n50.estimated_n50
            if info.HasField("protocol_run_info"):
                # print(info.protocol_run_info.meta_info.tags)
                self.running_kit = info.protocol_run_info.meta_info.tags[
                    "kit"
                ].string_value
                self.Flowcell_Type = info.protocol_run_info.meta_info.tags[
                    "flow cell"
                ].string_value
                # print(info.protocol_run_info.phase)
                if info.protocol_run_info.phase != 0:
                    self.show = True
                else:
                    self.show = False

                # self.protocol_run_info.clear()
                # with self.protocol_run_info:
                #    for field in info.protocol_run_info.DESCRIPTOR.fields:
                #        ui.label(f"{field.name}: {getattr(info.protocol_run_info, field.name)}")
                self.Experiment_Group = (
                    info.protocol_run_info.user_info.protocol_group_id.value
                )
                self.Sample_ID = info.protocol_run_info.user_info.sample_id.value
                # print(f"Experiment Group: {self.Experiment_Group}")
                """
                    #ui.label(f"Experiment Group: {self.Experiment_Group}")
                    s
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
                    # self.log("yield_summary")
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
            """

    async def start_run(
        self,
        position=None,
        reference=None,
        sample_id=None,
        centreID=None,
        flowcell_id=None,
        kit=None,
        basecall_config=None,
        experiment_duration=None,
        bed_file=None,
    ):
        ui.notify(f"Starting Run {sample_id} on {flowcell_id}!", type="positive")
        # ToDo: At every stage we need to confirm that the correct values have been entered.
        # position = "1B"
        kit = kit
        ###Memo to self... basecall config must not include cfg
        basecall_config = basecall_config
        alignment_reference = reference
        bed_file = bed_file
        experiment_duration = experiment_duration
        current_date = datetime.now()
        centreID = centreID
        experiment_group_id = (
            f"{centreID}_{current_date.strftime('%B')}_{current_date.year}"
        )
        # sample_id = "SAMPLE_ID"
        experiment_specs = []
        self.connection_ip = "127.0.0.1"
        self.manager = Manager(host=self.connection_ip)
        # Add all the positions to the list:
        for pos in self.manager.flow_cell_positions():
            if pos.name == position:
                experiment_specs.append(ExperimentSpec(position=pos))
        # Check if the flowcell ID is correct

        if add_protocol_ids(experiment_specs, kit, basecall_config, flowcell_id):

            # Build arguments for starting protocol:
            alignment_args = protocols.AlignmentArgs(
                reference_files=[alignment_reference],
                bed_file=bed_file,
            )

            basecalling_args = protocols.BasecallingArgs(
                config=basecall_config,
                barcoding=None,
                alignment=alignment_args,
            )

            read_until_args = protocols.ReadUntilArgs(
                filter_type="enrich",
                reference_files=[alignment_reference],
                bed_file=bed_file,
                first_channel=None,
                last_channel=None,
            )

            bam_arguments = protocols.OutputArgs(
                reads_per_file=4000,
                batch_duration="1",
            )
            pod5_arguments = protocols.OutputArgs(
                reads_per_file=4000, batch_duration="1"
            )

            # Now start the protocol(s):
            for spec in experiment_specs:
                position_connection = spec.position.connect()

                # Generate stop criteria for use by Run Until
                # The `runtime` is in seconds, while the `experiment_duration` is in hours
                stop_criteria = protocols.CriteriaValues(
                    runtime=int(experiment_duration * 60 * 60)
                )

                run_id = protocols.start_protocol(
                    position_connection,
                    identifier=spec.protocol_id,
                    sample_id=sample_id,
                    experiment_group=experiment_group_id,
                    barcode_info=None,
                    basecalling=basecalling_args,
                    read_until=read_until_args,
                    fastq_arguments=None,
                    fast5_arguments=None,
                    pod5_arguments=pod5_arguments,
                    bam_arguments=bam_arguments,
                    disable_active_channel_selection=False,
                    mux_scan_period=1.5,
                    stop_criteria=stop_criteria,
                    args=[],  # Any extra args passed.
                )

                flow_cell_info = position_connection.device.get_flow_cell_info()

                ui.notify(
                    f"Started protocol:\n    run_id={run_id}\n    position={spec.position.name}\n    flow_cell_id={flow_cell_info.flow_cell_id}\n",
                    multi_line=True,
                    type="positive",
                )
        else:
            ui.notify("Run Start Failed", type="negative")


@ui.page("/", response_timeout=30)
def index_page() -> None:
    initial_ip = "127.0.0.1"
    my_connection = Manager(host=initial_ip)
    with theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ):
        # my_connection.connect_to_minknow()
        positions = list(my_connection.flow_cell_positions())
        ui.label(f"{positions[0]}")
        display_object = Minknow_Info(positions[0])


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    ui.add_css("""
        .monitoring-tile {
            transition: all 0.3s ease;
        }
        .monitoring-tile:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }
        .status-chip {
            transition: background-color 0.3s ease;
        }
        .device-header {
            border-bottom: 1px solid #e0e0e0;
        }
        @media (max-width: 600px) {
            .monitoring-grid {
                grid-template-columns: repeat(2, 1fr) !important;
            }
        }
    """)
    
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
