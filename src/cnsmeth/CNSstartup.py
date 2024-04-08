import click
from configparser import ConfigParser
from cnsmeth.brain_class import BrainMeth
from pathlib import Path
from nicegui import ui, app
import os
import sys
import asyncio
from cnsmeth import images

from cnsmeth import theme
from cnsmeth.minknow_info.minknow_panel import MinKNOWFish

DEFAULT_CFG = "config.ini"


def configure(ctx, param, filename):
    cfg = ConfigParser()
    cfg.read(filename)
    try:
        options = dict(cfg["options"])
    except KeyError:
        options = {}
    print(options)
    ctx.default_map = options


class Methnice:
    def __init__(
        self,
        threads: int,
        simtime: bool,
        watchfolder: Path,
        output: Path,
        sequencing_summary: Path,
        target_panel: str,
        showerrors: bool,
        browse: bool,
        exclude: list,
        reference: Path,
        # minknow_connection,
    ):
        self.threads = threads
        self.simtime = simtime
        self._watchfolder = watchfolder
        self.output = output
        self.sequencing_summary = sequencing_summary
        self.target_panel = target_panel
        self.showerrors = showerrors
        self.browse = browse
        self.exclude = exclude
        self.reference = reference
        self.minknow_connection = None
        self.timer = ui.timer(1, self._worker)

    @property
    def watchfolder(self):
        return self._watchfolder

    @watchfolder.setter
    def watchfolder(self, value):
        self._watchfolder = value
        self.on_watchfolder_changed()

    def _worker(self):
        if hasattr(self.minknow_connection, "minKNOW_display"):
            if self.minknow_connection.minKNOW_display.watchfolder != self.watchfolder:
                self.watchfolder = self.minknow_connection.minKNOW_display.watchfolder

    def on_watchfolder_changed(self):
        # define your function here
        if self.watchfolder:
            print(f"watchfolder value has been changed! {self.watchfolder}")
            self.robin.add_watchfolder(self.watchfolder)

    @ui.page("/home")
    async def index_page(self) -> None:
        with theme.frame(
            "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>"
        ):
            with ui.column().classes("w-full"):
                if self.watchfolder is None:
                    self.minknow_connection = MinKNOWFish()
                else:
                    self.minknow_connection = None
                if self.minknow_connection:
                    with ui.splitter(value=10).classes("w-full h-full") as splitter:
                        with splitter.before:
                            with ui.tabs().props("vertical").classes("w-full") as tabs:
                                analysis = ui.tab("Analysis", icon="analytics")
                                readfish = ui.tab("ReadFish", icon="phishing")
                                minknow = ui.tab("minKNOW", icon="set_meal")
                        with splitter.after:
                            with ui.tab_panels(tabs, value=minknow).props(
                                "vertical"
                            ).classes("w-full h-full"):
                                self.analysis_tab_pane = ui.tab_panel(analysis)
                                with self.analysis_tab_pane:
                                    with ui.row():
                                        ui.icon("analytics", color="primary").classes(
                                            "text-h4"
                                        )
                                        ui.label("Real Time Analysis").classes(
                                            "text-h4"
                                        )
                                with ui.tab_panel(readfish):
                                    with ui.row():
                                        ui.icon("phishing", color="primary").classes(
                                            "text-h4"
                                        )
                                        ui.label("ReadFish Data").classes("text-h4")
                                    ui.label("To Be Updated.")
                                self.minknow_tab_pane = ui.tab_panel(minknow)
                                with self.minknow_tab_pane:
                                    with ui.row():
                                        ui.icon("set_meal", color="primary").classes(
                                            "text-h4"
                                        )
                                        ui.label("MinKNOW Data").classes("text-h4")

                else:
                    self.analysis_tab_pane = ui.row()

            if self.minknow_connection:
                with self.minknow_tab_pane:
                    self.minknow_connection.setup_ui()
                    self.minknow_connection.check_connection()
                    ui.label().bind_text_from(
                        self.minknow_connection,
                        "connection_ip",
                        backward=lambda n: f"Connected to: {n}",
                    )

            with self.analysis_tab_pane:
                self.robin = BrainMeth(
                    threads=self.threads,
                    simtime=self.simtime,
                    watchfolder=self.watchfolder,
                    output=self.output,
                    sequencing_summary=self.sequencing_summary,
                    target_panel=self.target_panel,
                    showerrors=self.showerrors,
                    browse=self.browse,
                    exclude=self.exclude,
                    minknow_connection=self.minknow_connection,
                    reference=self.reference,
                )
                await self.robin.init()
                # self.watchfolder = self.minknow_connection.watchfolder


def run_class(
    port: int,
    reload: bool,
    threads: int,
    simtime: bool,
    watchfolder: Path,
    output: Path,
    sequencing_summary: Path,
    target_panel: str,
    showerrors: bool,
    browse: bool,
    exclude: list,
    reference: Path,
    # minknow_connection,
):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    iconfile = os.path.join(
        os.path.dirname(os.path.abspath(images.__file__)), "favicon.ico"
    )
    mainpage = Methnice(
        threads=threads,
        simtime=simtime,
        watchfolder=watchfolder,
        output=output,
        sequencing_summary=sequencing_summary,
        target_panel=target_panel,
        showerrors=showerrors,
        browse=browse,
        exclude=exclude,
        reference=reference,
        # minknow_connection=minknow_connection,
    )
    # ui.add_head_html(r'''
    # <style>
    # @import url('https://fonts.googleapis.com/css2?family=Shadows+Into+Light&display=swap')
    # </style>
    #''')

    ui.add_style(
        """
        .shadows-into light-regular {
            font-family: "Shadows Into Light", cursive;
            font-weight: 400;
            font-style: normal;
        }
    """
    )
    app.add_static_files("/fonts", str(Path(__file__).parent / "fonts"))
    app.on_startup(mainpage.index_page)
    # app.on_startup(startup)
    ui.run(
        port=port, reload=reload, title="ROBIN", favicon=iconfile, on_air=False
    )  # , native=True, fullscreen=False, window_size=(1200, 1000))


@click.command()
@click.option(
    "-c",
    "--config",
    type=click.Path(dir_okay=False),
    default=DEFAULT_CFG,
    callback=configure,
    is_eager=True,
    expose_value=False,
    help="Read option defaults from the specified INI file",
    show_default=True,
)
@click.option(
    "--port",
    default=8081,
    help="Port for GUI",
)
@click.option(
    "--threads", default=4, help="Number of threads available.", required=True
)
@click.option(
    "--simtime",
    default=False,
    help="If set, will simulate the addition of existing files to the pipeline based on read data.",
)
# @click.option(
#    "--showerrors",
#    default=False,
#    help="If set, will display all errors in running R.",
# )
@click.option(
    "--showerrors/--noerrors",
    default=False,
    help="If set, will display all errors in running R.",
)
@click.option(
    "--sequencing_summary",
    default=None,
    help="Path to sequencing summary file. If provided, timestamps will be taken from this file.",
)
@click.option(
    "--target_panel",
    "-t",
    default="rCNS2",
    help="Select analysis gene panel from one of these options. Default is rCNS2",
    type=click.Choice(
        ["rCNS2", "AML"],
        case_sensitive=True,
    ),
)
@click.option(
    "--reference",
    "-r",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, resolve_path=True, path_type=Path
    ),
    help="Path to the reference genome and index.",
    required=True,
)
@click.option(
    "--browse",
    is_flag=True,
    show_default=True,
    default=False,
    help="Browse Historic Data.",
)
@click.option(
    "--exclude",
    "-e",
    multiple=True,
    help="Exclude analysis types with one or more of these options.",
    type=click.Choice(
        ["sturgeon", "forest", "nanodx", "cnv", "fusion", "coverage", "mgmt"],
        case_sensitive=False,
    ),
)
@click.option(
    "--watchfolder",
    "-w",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.argument(
    "output",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
def package_run(
    port,
    threads,
    simtime,
    showerrors,
    sequencing_summary,
    target_panel,
    watchfolder,
    output,
    browse,
    exclude,
    reference,
):  # , threads, simtime, watchfolder, output, sequencing_summary):
    """
    Entrypoint for when GUI is launched directly.
    :return: None
    """
    if browse:
        # Handle the case when --browse is set
        click.echo("Browse mode is enabled. Watchfolder and output are not required.")
        run_class(
            port=port,
            reload=False,
            threads=threads,
            simtime=simtime,
            watchfolder=None,
            output=None,
            sequencing_summary=sequencing_summary,
            target_panel=target_panel,
            showerrors=showerrors,
            browse=browse,
            exclude=exclude,
        )
        # Your logic for browse mode
    else:
        # Handle the case when --browse is not set
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        # if watchfolder is None or output is None:
        #    click.echo("Watchfolder and output are required when --browse is not set.")
        #    sys.exit(1)
        if output is None:
            click.echo("Output is required when --browse is not set.")
            sys.exit(1)

        # Your logic for non-browse mode

        run_class(
            port=port,
            reload=False,
            threads=threads,
            simtime=simtime,
            watchfolder=watchfolder,
            output=output,
            sequencing_summary=sequencing_summary,
            target_panel=target_panel,
            showerrors=showerrors,
            browse=browse,
            exclude=exclude,
            reference=reference,
            # minknow_connection=minknow_connection,
        )


def startup():
    loop = asyncio.get_running_loop()
    loop.set_debug(True)
    loop.slow_callback_duration = 0.2


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    run_class(port=12398, reload=True)
