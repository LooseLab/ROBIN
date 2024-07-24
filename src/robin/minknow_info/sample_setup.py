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
def add_protocol_ids(experiment_specs, kit, basecall_config):
    for spec in experiment_specs:
        # Connect to the sequencing position:
        position_connection = spec.position.connect()
        print (position_connection)
        # Check if a flowcell is available for sequencing
        flow_cell_info = position_connection.device.get_flow_cell_info()
        if not flow_cell_info.has_flow_cell:
            print("No flow cell present in position {}".format(spec.position))
            ui.notify("No flow cell present in position {}".format(spec.position))
            return

        print(flow_cell_info)

        product_code = flow_cell_info.user_specified_product_code
        if not product_code:
            product_code = flow_cell_info.product_code

        print(product_code)

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

        print (protocol_info)

        if not protocol_info:
            print("Failed to find protocol for position %s" % (spec.position))
            print("Requested protocol:")
            print("  product-code: %s" % product_code)
            print("  kit: %s" % kit)
            #ui.notify("Failed to find protocol for position %s" % (spec.position))
            #ui.notify("Requested protocol:")
            #ui.notify("  product-code: %s" % product_code)
            #ui.notify("  kit: %s" % kit)
            return


        # Store the identifier for later:
        spec.protocol_id = protocol_info.identifier


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
        # position = "1B"
        kit = "SQK-RAD114"
        ###Memo to self... basecall config must not include cfg
        basecall_config = "dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac"
        alignment_reference = reference
        bed_file = "/home/deepseq/panel_adaptive_nogenenames_20122021_hg38.bed"
        experiment_duration = 24
        experiment_group_id = "experiment_group"
        # sample_id = "SAMPLE_ID"
        experiment_specs = []
        # Add all the positions to the list:
        ui.notify(f'looking for position {position}')
        for pos in self.manager.flow_cell_positions():
            print(pos.name, position)
            if pos.name == position:
                experiment_specs.append(ExperimentSpec(position=pos))
                print(f"Found {pos}")
        print(experiment_specs)
        print(basecall_config)
        add_protocol_ids(experiment_specs, kit, basecall_config)

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

        print(basecalling_args)

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
        print("Starting protocol on %s positions" % len(experiment_specs))
        for spec in experiment_specs:
            position_connection = spec.position.connect()

            # Generate stop criteria for use by Run Until
            # The `runtime` is in seconds, while the `experiment_duration` is in hours
            stop_criteria = protocols.CriteriaValues(
                runtime=int(experiment_duration * 60 * 60)
            )
            print("stop_criteria")
            print(stop_criteria)

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
                args=['--verbose'],  # Any extra args passed.
            )

            flow_cell_info = position_connection.device.get_flow_cell_info()

            print("Started protocol:")
            print("    run_id={}".format(run_id))
            print("    position={}".format(spec.position.name))
            print("    flow_cell_id={}".format(flow_cell_info.flow_cell_id))
            print(
                "    user_specified_flow_cell_id={}".format(
                    flow_cell_info.user_specified_flow_cell_id
                )
            )


@ui.page("/")
async def content():
    with theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ):
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
        #await minknow.connect_to_remotehost("10.157.252.30")
        await minknow.connect_to_localhost()
        with ui.dialog() as startup, ui.card().classes("w-full"):
            ui.label("Sample Setup").style(
                "color: #6E93D6; font-size: 200%; font-weight: 300"
            )
            ui.separator()
            with ui.stepper().props("vertical").classes("w-full") as stepper:
                step = ui.step("Sample ID")
                with step:
                    ui.label("Enter the sample ID:")
                    ui.label(
                        "Remember that sample IDs may be shared with others and should not include human identifiable information."
                    )
                    sampleid = ui.input(
                        placeholder="start typing",
                        validation=lambda value: (
                            "Too short" if len(value) < 5 else None
                        ),
                    ).props("rounded outlined dense")
                    sample_camera.show_camera()
                    with ui.stepper_navigation():
                        ui.button("Next", on_click=stepper.next)
                with ui.step("Flowcell"):
                    ui.label("Enter the flowcell ID.")
                    flowcellid = ui.input(
                        placeholder="start typing",
                        validation=lambda value: (
                            "Not a valid flowcell ID"
                            if not re.match(r"^[A-Za-z]{3}\d{5}$", value)
                            else None
                        ),
                    ).props("rounded outlined dense")
                    flowcell_camera.show_camera()
                    with ui.stepper_navigation():
                        ui.button("Back", on_click=stepper.previous).props("flat")
                        ui.button("Next", on_click=stepper.next)
                with ui.step("Device Position"):
                    ui.label("Select the device position to be used.")
                    position = ui.radio(
                        {
                            item: str(item)
                            for index, item in enumerate(minknow.positions)
                        },
                        # on_change=minknow.enable_access,
                    )  # .bind_value(self, "selected_position").props("inline")
                    with ui.stepper_navigation():
                        ui.button("Back", on_click=stepper.previous).props("flat")
                        ui.button(
                            "Start Run",
                            on_click=lambda: minknow.start_run(
                                position=position.value.description.name,
                                reference="/home/deepseq/refs/hg38_simple.fa",
                                sample_id=sampleid.value,
                                flowcell_id=flowcellid.value,
                            ),
                        )
            ui.separator()
            ui.label("Sample Settings:").style(
                "color: #6E93D6; font-size: 100%; font-weight: 300"
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

"""

import logging
from pathlib import Path
import sys

# minknow_api.manager supplies "Manager" a wrapper around MinKNOW's Manager gRPC API with utilities
# for querying sequencing positions + offline basecalling tools.
from enum import Enum
from typing import Sequence

from minknow_api.examples.load_sample_sheet import (
    ParsedSampleSheetEntry,
    SampleSheetParseError,
    load_sample_sheet_csv,
)
from minknow_api.manager import Manager

# We need `find_protocol` to search for the required protocol given a kit + product code.
from minknow_api.tools import protocols

class ExperimentSpec(object):
    def __init__(self, position):
        self.position = position
        self.protocol_id = ""


ExperimentSpecs = Sequence[ExperimentSpec]


# Determine which protocol to run for each experiment, and add its ID to experiment_specs
def add_protocol_ids(experiment_specs, kit, basecall_config):
    for spec in experiment_specs:
        # Connect to the sequencing position:
        position_connection = spec.position.connect()

        # Check if a flowcell is available for sequencing
        flow_cell_info = position_connection.device.get_flow_cell_info()
        if not flow_cell_info.has_flow_cell:
            print("No flow cell present in position {}".format(spec.position))
            sys.exit(1)

        product_code = flow_cell_info.user_specified_product_code
        if not product_code:
            product_code = flow_cell_info.product_code

        # Find the protocol identifier for the required protocol:
        protocol_info = protocols.find_protocol(
            position_connection,
            product_code=product_code,
            kit=kit,
            basecalling=True,
            basecall_config=basecall_config,
            barcoding=False,
            barcoding_kits=None,
        )

        if not protocol_info:
            print("Failed to find protocol for position %s" % (spec.position))
            print("Requested protocol:")
            print("  product-code: %s" % product_code)
            print("  kit: %s" % kit)
            sys.exit(1)

        # Store the identifier for later:
        spec.protocol_id = protocol_info.identifier


def main():
    position = "1B"
    kit = "SQK-RAD114"
    basecall_config = "dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac.cfg"
    alignment_reference = "/data/lambda.fasta"
    bed_file = "/data/test.bed"
    experiment_duration = 72
    experiment_group_id = "experiment_group"
    sample_id = "SAMPLE_ID"

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    # Construct a manager using the host + port provided:
    manager = Manager(
        host="localhost",
        port=None,
    )

    experiment_specs = []
    # Add all the positions to the list:
    for pos in manager.flow_cell_positions():
        if pos.name == position:
            experiment_specs.append(ExperimentSpec(position=pos))
    add_protocol_ids(experiment_specs, kit, basecall_config)

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
        last_channel=None
    )

    bam_arguments = protocols.OutputArgs(
        reads_per_file=None, # 4000,
        batch_duration="60",
    )

    # Now start the protocol(s):
    print("Starting protocol on %s positions" % len(experiment_specs))
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
            pod5_arguments=None,
            bam_arguments=bam_arguments,
            disable_active_channel_selection=False,
            mux_scan_period=1.5,
            stop_criteria=stop_criteria,
            args=None,  # Any extra args passed.
        )

        flow_cell_info = position_connection.device.get_flow_cell_info()

        print("Started protocol:")
        print("    run_id={}".format(run_id))
        print("    position={}".format(spec.position.name))
        print("    flow_cell_id={}".format(flow_cell_info.flow_cell_id))
        print(
            "    user_specified_flow_cell_id={}".format(
                flow_cell_info.user_specified_flow_cell_id
            )
        )


if __name__ == "__main__":
    main()


"""
