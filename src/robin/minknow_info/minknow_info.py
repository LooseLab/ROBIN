# Python imports.
from __future__ import annotations
from nicegui import ui, app

from robin import images
from robin import __about__

from pathlib import Path

import threading

from datetime import datetime


from robin import theme

from minknow_api.manager import Manager

import os

DEVICEDICT = {
    "FLONGLE": "Flongle",
    "GRIDION": "GridION",
    "MINION": "MinION",
    "PROMETHION": "PromethION",
    "P2_INTEGRATED": "P2Integrated"
}


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
    elif device_type == "PROMETHION":
        if device_name.startswith("PS2"):
            return "P2 Solo", "P2_Solo_Left-45_Open_Full.png"
        else:
            return (
                f"{DEVICEDICT[device_type]}",
                "PromethION_24_Right 45 Elevated_Closed_Full_white.jpg",
            )
    else:
        return (f"{device_type}", None)


class Minknow_Info:
    def __init__(self, position, dev=False):
        self.dev = dev
        self.color = "text-blue"
        self.position = position
        self.connection = self.position.connect()
        self.device = self.position.device_type
        self.show = False
        self.name, self.image = identify_device(
            self.position.device_type, self.position.name
        )
        if self.image:
            self.deviceicon = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)), "ONTimages", self.image
            )
        else:
            self.deviceicon = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)), "ONTimages", "unknown_sequencer.png"
            )
        self.check_instance = threading.Thread(
            target=self.stream_instance_activity, args=()
        )
        self.check_instance.daemon = True
        self.check_instance.start()

        self.render_me()

    def render_me(self):
        if self.dev:
            ui.label(f"{dir(self.position.connect)}")
            ui.label(f"{self.position.credentials}")
            ui.label(f"{self.position.description}")
            ui.label(f"{self.position.device_type}")
            ui.label(f"{self.position.host}")
            ui.label(f"{self.position.name}")
            ui.label(f"{self.position.running}")
            ui.label(f"{self.position.state}")
        with ui.card().classes("w-full"):  # .classes("flat border-[2px] no-shadow"):
            with ui.card().tight().classes("flat border-[2px] no-shadow"):
                ui.label(f"{self.name} - {self.position}").classes("text-h6")
                with ui.row().props("align-middle"):
                    with ui.avatar(square=False, color="white"):
                        ui.image(self.deviceicon).classes("w-full h-full")
                    ui.label("MinKNOW Monitoring.").classes(
                        f"text-overline {self.color}"
                    )
                ui.label("Data from the current experiment.").classes("text-subtitle")
            with ui.card().tight().classes("flat no-shadow").bind_visibility_from(
                self, "show"
            ):
                with ui.grid(columns=2).classes("gap-0 p-0"):
                    with ui.column().classes("gap-0 p-0"):
                        ui.label("Experiment Group:").classes(
                            f"gap-0 p-0 text-overline {self.color}"
                        )
                        ui.label("None").bind_text_from(
                            self, "Experiment_Group", backward=lambda n: f"{n}"
                        ).classes("gap-0 p-0")
                    with ui.column().classes("gap-0 p-0"):
                        ui.label("Sample Name:").classes(
                            f"gap-0 p-0 gap-0 text-overline {self.color}"
                        )
                        ui.label("None").bind_text_from(
                            self, "Sample_ID", backward=lambda n: f"{n}"
                        ).classes("gap-0 p-0")
                    with ui.column().classes("gap-0 p-0"):
                        ui.label("Flowcell:").classes(
                            f"gap-0 p-0 text-overline {self.color}"
                        )
                        ui.label("None").bind_text_from(
                            self, "Flowcell_Type", backward=lambda n: f"{n}"
                        ).classes("gap-0 p-0")
                    with ui.column().classes("gap-0 p-0"):
                        ui.label("Kit ID:").classes(
                            f"gap-0 p-0 text-overline {self.color}"
                        )
                        ui.label("None").bind_text_from(
                            self, "kit", backward=lambda n: f"{n}"
                        ).classes("gap-0 p-0")
                    with ui.column().classes("gap-0 p-0"):
                        ui.label("Output Directory:").classes(
                            f"gap-0 p-0 text-overline {self.color}"
                        )
                        ui.label("None").bind_text_from(
                            self, "output_folder", backward=lambda n: f"{n}"
                        ).classes("gap-0 p-0")
                    with ui.column().classes("gap-0 p-0"):
                        ui.label("Basecall Model:").classes(
                            f"gap-0 p-0 text-overline {self.color}"
                        )
                        ui.label("None").bind_text_from(
                            self,
                            "basecalling_config_filename",
                            backward=lambda n: f"{n}",
                        ).classes("gap-0 p-0")
            with ui.column().bind_visibility_from(
                self, "show", backward=lambda v: not v
            ):
                ui.label("Set up a run using MinKNOW to see more information.")

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
                self.kit = info.protocol_run_info.meta_info.tags["kit"].string_value
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
        port=port, reload=reload, title="Readfish NiceGUI", storage_secret="waynesworld",
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
