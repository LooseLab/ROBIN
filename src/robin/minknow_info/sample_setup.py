from nicegui import ui, run, app

import asyncio
import re
from robin import theme

from contextlib import contextmanager
import ipaddress

# MinKNOW API Imports
from minknow_api.manager import Manager
from robin.utilities.camera import Camera

from typing import Sequence
from pathlib import Path
import sys
import uuid

from datetime import datetime


UNIQUE_ID: str = str(uuid.uuid4())



import cv2
import zxingcpp

import numpy as np


import base64

# We need `find_protocol` to search for the required protocol given a kit + product code.
from minknow_api.tools import protocols


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
            ui.notify(f"Flowcell {expected_flowcell_id} is not found in position {spec.position}. Please check.", type="negative")
            return
        if not flow_cell_info.has_flow_cell:
            ui.notify("No flow cell present in position {}".format(spec.position), type="negative")
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
            ui.notify("Failed to find protocol for position %s" % (spec.position),type='negative')

            #print("Requested protocol:")
            #print("  product-code: %s" % product_code)
            #print("  kit: %s" % kit)
            #ui.notify("Failed to find protocol for position %s" % (spec.position))
            #ui.notify("Requested protocol:")
            #ui.notify("  product-code: %s" % product_code)
            #ui.notify("  kit: %s" % kit)
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
        self.connected = False
        self.watchfolder = None
        self.devices = set()

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

    async def connect_to_localhost(self):
        self.connection_ip = "127.0.0.1"
        await self.add_positions()

    async def connect_to_remotehost(self, remotehost):
        self.connection_ip = remotehost
        await self.add_positions()

    async def add_positions(self) -> None:
        """
        This is called when the connectionip variable is changed.
        It adds a tab for each position on the flowcell.
        It does this using a thread so that the app doesn't freeze
        while it waits for minKNOW to respond.
        """
        # with disable(self.connect_now):
        ui.notify("Trying to connect to MinKNOW")
        await run.io_bound(self._connect_positions)
        if self.connected:
            ui.notify("Connection Successful.", type="positive")
            # with self.connect_now:
            ui.notify("Connected to MinKNOW - getting positions.")
            self.positions = list(self.manager.flow_cell_positions())
            for position in self.positions:
                self.devices.add(position.device_type)
            # self.connectiondialog.close()
            # with self.choices:
            #    ui.radio(
            #        {item: str(item) for index, item in enumerate(self.positions)},
            #        on_change=self.enable_access,
            #    ).bind_value(self, "selected_position").props("inline")
            # self.position_choose.open()
        else:
            self.connection_ip = None
            # with self.connect_now:
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

    async def start_run(
        self,
        position=None,
        reference=None,
        sample_id=None,
        experiment_group_id=None,
        flowcell_id=None,
    ):
        ui.notify(f"Starting Run {sample_id} on {flowcell_id}!", type="positive")
        #ToDo: At every stage we need to confirm that the correct values have been entered.
        # position = "1B"
        kit = "SQK-RAD114"
        ###Memo to self... basecall config must not include cfg
        basecall_config = "dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg"
        alignment_reference = reference
        bed_file = "/home/deepseq/panel_adaptive_nogenenames_20122021_hg38.bed"
        experiment_duration = 24
        current_date = datetime.now()
        centreID = "NUH"
        experiment_group_id = f"{centreID}_{current_date.strftime('%B')}_{current_date.year}"
        # sample_id = "SAMPLE_ID"
        experiment_specs = []
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
                reads_per_file=4000,
                batch_duration="1"
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

                ui.notify(f"Started protocol:\n    run_id={run_id}\n    position={spec.position.name}\n    flow_cell_id={flow_cell_info.flow_cell_id}\n",
                    multi_line = True,
                    type="positive",
                )
        else:
            ui.notify("Run Start Failed", type="negative")


class ErrorChecker:
    def __init__(self, *elements) -> None:
        self.elements = elements

    @property
    def no_errors(self) -> bool:
        return all(validation(element.value) for element in self.elements for validation in
                   element.validation.values())

@ui.page("/")
async def content():
    with ((theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ))):
        sample_camera = Camera(
            icon="photo_camera",
            icon_color="blue-5",
            background_color="rgba(66, 165, 245, 0.3)",
            for_id="samplefileupload",
            canvas_id="samplecanvas",
            on_change=lambda: generate_sampleID(),
        )
        async def generate_sampleID():
            await asyncio.sleep(0.1)
            print("calling camera sample")
            image = await sample_camera.get_image()
            base64_data = image.split(",")[1]
            image_data = base64.b64decode(base64_data)
            nparr = np.frombuffer(image_data, np.uint8)
            img = cv2.imdecode(nparr, cv2.IMREAD_COLOR)
            barcodes = zxingcpp.read_barcodes(img)
            if len(barcodes) == 0:
                sampleid.value = "Could not find any barcode."
            elif len(barcodes) == 1:
                sampleid.value = barcodes[0].text
            elif len(barcodes) >= 1:
                sampleid.value = "Too many barcodes - please try again."

        flowcell_camera = Camera(
            icon="photo_camera",
            icon_color="blue-5",
            background_color="rgba(66, 165, 245, 0.3)",
            for_id="flowcellfileupload",
            canvas_id="flowcellcanvas",
            on_change=lambda: generate_flowcellID(),
        )

        async def generate_flowcellID():
            await asyncio.sleep(0.1)
            image = await flowcell_camera.get_image()
            base64_data = image.split(",")[1]
            image_data = base64.b64decode(base64_data)
            nparr = np.frombuffer(image_data, np.uint8)
            img = cv2.imdecode(nparr, cv2.IMREAD_COLOR)
            barcodes = zxingcpp.read_barcodes(img)
            if len(barcodes) == 0:
                flowcellid.value = "Could not find any barcode."
            elif len(barcodes) == 1:
                flowcellid.value = barcodes[0].text
            elif len(barcodes) >= 1:
                flowcellid.value = "Too many barcodes - please try again."

        minknow = MinKNOWFish()
        await minknow.connect_to_remotehost("10.157.252.30")
        #await minknow.connect_to_localhost()
        with ui.dialog() as startup, ui.card().classes("w-full"):
            ui.label("Sample Setup").style(
                "color: #6E93D6; font-size: 200%; font-weight: 300"
            )
            ui.separator()
            with ui.stepper().props("vertical").classes("w-full") as stepper:
                step = ui.step("Sample ID").style('color: #000000; font-size: 100%; font-weight: 600')
                with step:
                    ui.label(
                        "Remember that sample IDs may be shared with others and should not include human identifiable information."
                    ).style('color: #000000; font-size: 100%; font-weight: 600')
                    ui.label("Enter the sample ID:").style('color: #000000; font-size: 80%; font-weight: 300')
                    sampleid = ui.input(
                        placeholder="start typing",
                        validation = {
                            'Too short': lambda value: len(value) >= 5,
                            'Do not use Flowcell IDs as sample IDs': lambda value: re.match(r"^[A-Za-z]{3}\d{5}$", value) is None,
                            'No whitespace allowed': lambda value: not any(char.isspace() for char in value),
                            'Only alphanumeric characters allowed': lambda value: value.isalnum(),
                        }
                    ).style('color: #000000; font-size: 80%; font-weight: 300').props("rounded outlined dense")
                    sample_camera.show_camera()
                    with ui.stepper_navigation():
                        checker = ErrorChecker(sampleid)
                        c = ui.button("Next", on_click=stepper.next).bind_enabled_from(checker, 'no_errors')

                with ui.step("Flowcell").style('color: #000000; font-size: 100%; font-weight: 600'):
                    ui.label("Enter the flowcell ID.").style('color: #000000; font-size: 80%; font-weight: 300')
                    flowcellid = ui.input(
                        placeholder="start typing",
                        validation = {'Not a valid flowcell ID': lambda value: re.match(r"^[A-Za-z]{3}\d{5}$", value) is not None}
                    ).style('color: #000000; font-size: 80%; font-weight: 300').props("rounded outlined dense")
                    flowcell_camera.show_camera()
                    with ui.stepper_navigation():
                        ui.button("Back", on_click=stepper.previous).props("flat")
                        checkerflowcell = ErrorChecker(flowcellid)
                        d=ui.button("Next", on_click=stepper.next).bind_enabled_from(checkerflowcell, 'no_errors')
                with ui.step("Device Position").style('color: #000000; font-size: 100%; font-weight: 600'):
                    ui.label("Select the device position to be used.").style('color: #000000; font-size: 80%; font-weight: 300')
                    run_button = ui.button(
                        "Start Run",
                        on_click=lambda: minknow.start_run(
                            position=position.value.description.name,
                            reference="/home/deepseq/refs/hg38_simple.fa",
                            sample_id=sampleid.value,
                            flowcell_id=flowcellid.value,
                        ),
                    )
                    def dostuff():
                        run_button.enable()

                    position = ui.radio(
                        {
                            item: str(item)
                            for index, item in enumerate(minknow.positions)
                        },
                    ).on(
                        "update:model-value", dostuff
                    ).style('color: #000000; font-size: 80%; font-weight: 300')
                    finalstep = ui.stepper_navigation()
                    with finalstep:
                        ui.button("Back", on_click=stepper.previous).props("flat")
                        run_button.move(finalstep)
                        run_button.disable()

            ui.separator()
            ui.label("Fixed Settings:").style(
                "color: #6E93D6; font-size: 100%; font-weight: 400"
            )
            ui.label(" ".join(f"Device: {device}" for device in minknow.devices))
            current_date = datetime.now()
            centreID = "NUH"
            ui.label(f"experiment_group_id = {centreID}_{current_date.strftime('%B')}_{current_date.year}")
            ui.label("kit = SQK-RAD114")
            ui.label("reference = {reference}")
            ui.separator()
            ui.label("User Settings:").style(
                "color: #6E93D6; font-size: 100%; font-weight: 400"
            )
            ui.label().bind_text_from(
                sampleid,
                "value",
                backward=lambda sid: (
                    f"Sample ID: {sid}"
                    if sid is not None
                    else "Sample ID: Not specified"
                ),
            )
            ui.label().bind_text_from(
                flowcellid,
                "value",
                backward=lambda fid: (
                    f"Flowcell ID: {fid}"
                    if fid is not None
                    else "Flowcell ID: Not specified"
                ),
            )
            ui.label().bind_text_from(
                position,
                "value",
                backward=lambda pid: (
                    f"Position: {pid.description.name}"
                    if pid is not None
                    else "Position: Not selected"
                ),
            )

            ui.button("Close", on_click=startup.close)

        ui.button("Run Setup and Start", on_click=startup.open)


def main():
    # with theme.frame(
    #        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
    #        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    # ):
    #    ui.label("hello")
    app.add_static_files("/fonts", str(Path(__file__).parent / "../fonts"))
    ui.run(
        storage_secret="UNIQUE_ID",
    )


if __name__ in {"__main__", "__mp_main__"}:
    main()