"""
This is the main page which sets up the ROBIN real time application.

The code implements a NiceGUI app (https://nicegui.io/) which requires all long running processes to
operate in the background. The fundamental principle is that nothing should block the 
main thread.

On launch, the app creates an instance of the run class for the specific run, taking arguments either
from the command line or via a config file.

The run class in turn sets up the application. Arguments to be shared across the application
are written to niceguis app.storage.general to allow them to be accessed by all users of this
instance.

The run class then configures the app with a specific on_startup function. This function handles
all slow running background processes. However, the function shares the main event loop with the ui
process. It is important to exploit run.cpu_bound and run.io_bound wherever necessary.

The result of this is that a single instance of the Methnice class is created which will process data 
in the background. Whenever a user connects to the site, the user will receive individual 
instances of these classes which will enable visualisation of the data but will not themselves
run any analysis.

The concept of one class providing both analysis and ui elements is maintained throughout the codebase.

Hopefully this simplifies code maintenance and future improvements.


"""

import click
import os
import sys
import signal
import uuid
import asyncio

from pathlib import Path
from nicegui import ui, app
from cnsmeth import theme
from cnsmeth import images
from configparser import ConfigParser

from cnsmeth.brain_class import BrainMeth
from cnsmeth.minknow_info.minknow_panel import MinKNOWFish

### The location for a default config file.
# ToDo: Confirm if this actually works!!!
DEFAULT_CFG = "config.ini"

### This ID provides an access key for this specific run of NiceGUI.
### Data will be written to the general store using this key.
# ToDo: Confirm if multiple instances of ROBIN can run in the same folder.
UNIQUE_ID = str(uuid.uuid4())


@ui.page("/", response_timeout=10)
async def index():
    """
    This function defines the main index page for the ROBIN site.
    We setup an instance of Methnice which will provide the user with
    the "splash screen" and allow them to choose which pages to browse
    to.
    """
    GUI = Methnice(
        threads=app.storage.general[UNIQUE_ID]["threads"],
        simtime=app.storage.general[UNIQUE_ID]["simtime"],
        watchfolder=app.storage.general[UNIQUE_ID]["watchfolder"],
        output=app.storage.general[UNIQUE_ID]["output"],
        sequencing_summary=app.storage.general[UNIQUE_ID]["sequencing_summary"],
        target_panel=app.storage.general[UNIQUE_ID]["target_panel"],
        showerrors=app.storage.general[UNIQUE_ID]["showerrors"],
        browse=app.storage.general[UNIQUE_ID]["browse"],
        exclude=app.storage.general[UNIQUE_ID]["exclude"],
        reference=app.storage.general[UNIQUE_ID]["reference"],
        unique_id=UNIQUE_ID,
    )
    GUI.setup()
    await GUI.splash_screen()
    ui.link("Home", "/")
    ui.link("Live", "/live")
    ui.link("Browse", "/browse")


@ui.page("/live", response_timeout=10)
async def live():
    """
    This page is served at /live and is the interaction with live data from ROBIN.
    This page is created for every user who visits.
    #ToDo: Many of the variables being passed to Methnice are most likely redundant.
    """
    GUI = Methnice(
        threads=app.storage.general[UNIQUE_ID]["threads"],
        simtime=app.storage.general[UNIQUE_ID]["simtime"],
        watchfolder=app.storage.general[UNIQUE_ID]["watchfolder"],
        output=app.storage.general[UNIQUE_ID]["output"],
        sequencing_summary=app.storage.general[UNIQUE_ID]["sequencing_summary"],
        target_panel=app.storage.general[UNIQUE_ID]["target_panel"],
        showerrors=app.storage.general[UNIQUE_ID]["showerrors"],
        browse=app.storage.general[UNIQUE_ID]["browse"],
        exclude=app.storage.general[UNIQUE_ID]["exclude"],
        reference=app.storage.general[UNIQUE_ID]["reference"],
        unique_id=UNIQUE_ID,
    )
    GUI.setup()
    await GUI.index_page()


@ui.page("/browse", response_timeout=10)
async def test():
    """
    #ToDo: This code is to be implemented.
    The idea is that we will implement the browse functions that allow someone
    to explore previous data even when a live run is in progress.
    """
    ui.link("Home", "/")
    ui.link("Live", "/live")
    ui.link("Browse", "/browse")
    ui.label("Good to be back. This code is still to be implemented.")


def clean_up():
    """
    This function serves to clean up the stored data at the end of a run.
    #ToDo: We should wait to pop the data until after all the subrunning processes have
    been wound up.
    """
    print("App shutdown detected.")
    print("Resetting Storage.")
    # We have seen the app stop. Therefore we want to remove data important for this instance.
    app.storage.general.pop(UNIQUE_ID)
    app.shutdown()


async def startup():
    """
        This function starts data processing in the main application loop.
    The function is only started once. It must not be run more than once!
    """

    #loop = asyncio.get_running_loop()
    #loop.set_debug(True)
    #loop.slow_callback_duration = 0.2
    print(f"Setting up {UNIQUE_ID}.")
    MAINPAGE = Methnice(
        threads=app.storage.general[UNIQUE_ID]["threads"],
        simtime=app.storage.general[UNIQUE_ID]["simtime"],
        watchfolder=app.storage.general[UNIQUE_ID]["watchfolder"],
        output=app.storage.general[UNIQUE_ID]["output"],
        sequencing_summary=app.storage.general[UNIQUE_ID]["sequencing_summary"],
        target_panel=app.storage.general[UNIQUE_ID]["target_panel"],
        showerrors=app.storage.general[UNIQUE_ID]["showerrors"],
        browse=app.storage.general[UNIQUE_ID]["browse"],
        exclude=app.storage.general[UNIQUE_ID]["exclude"],
        reference=app.storage.general[UNIQUE_ID]["reference"],
        unique_id=UNIQUE_ID,
    )
    MAINPAGE.setup()
    await MAINPAGE.start_analysis()

    # handle any occasions when a user stops with a ctrl-c
    # ToDo: surely this can be replaced by the clean_up funtion itself.
    def handler(*args):
        clean_up()

    signal.signal(signal.SIGINT, handler=handler)


class Methnice:
    """
    This class handles configuration of all the pages and page contents.

    #ToDo: See if we can reduce or remove this.
    """

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
        unique_id: str,
        # minknow_connection,
    ):
        self.MAINID = unique_id
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
        # self.timer = ui.timer(1, self._worker)

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

    async def on_watchfolder_changed(self):
        # define your function here
        if self.watchfolder:
            print(f"watchfolder value has been changed! {self.watchfolder}")
            await self.robin.add_watchfolder(self.watchfolder)

    def setup(self):
        self.robin = BrainMeth(
            mainuuid=self.MAINID,
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

    async def start_analysis(self):
        self.timer = ui.timer(1, self._worker)

        await self.robin.init()

    # @ui.page("/browse_page")
    async def browse_page(self):
        with theme.frame(
            "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>"
        ):
            ui.label(f"private page with ID {uuid.uuid4()}")
            self.analysis_tab_pane = ui.row().classes("w-full")

            with self.analysis_tab_pane:
                self.robin_browse = BrainMeth(
                    mainuuid=self.MAINID,
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
                await self.robin_browse.init()

    async def splash_screen(self) -> None:
        with theme.frame(
            "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>"
        ):
            with ui.column().classes("w-full"):
                ui.label("Welcome to Robin.")

    async def index_page(self) -> None:
        with theme.frame(
            "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>"
        ):
            await ui.context.client.connected()
            with ui.column().classes("w-full"):
                if self.watchfolder is None and not self.browse:
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
                    self.analysis_tab_pane = ui.row().classes("w-full")

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
                self.robin.render_ui()


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
    # Make sure no previous data are in the repository.
    # ToDo: This will delete everything for any run of cnsmeth.
    try:
        app.storage.general.clear()
    except Exception as e:
        print(e)
    # Configure the icon for the GUI
    iconfile = os.path.join(
        os.path.dirname(os.path.abspath(images.__file__)), "favicon.ico"
    )

    # Here is where we do anything which we want done once and once only.
    # create a custom storage area for this run
    app.storage.general[UNIQUE_ID] = dict()
    app.storage.general[UNIQUE_ID]["threads"] = threads
    app.storage.general[UNIQUE_ID]["simtime"] = simtime
    app.storage.general[UNIQUE_ID]["watchfolder"] = watchfolder
    app.storage.general[UNIQUE_ID]["output"] = output
    app.storage.general[UNIQUE_ID]["sequencing_summary"] = sequencing_summary
    app.storage.general[UNIQUE_ID]["target_panel"] = target_panel
    app.storage.general[UNIQUE_ID]["showerrors"] = showerrors
    app.storage.general[UNIQUE_ID]["browse"] = browse
    app.storage.general[UNIQUE_ID]["exclude"] = exclude
    app.storage.general[UNIQUE_ID]["reference"] = reference

    # Add some custom CSS because - why not!
    ui.add_css(
        """
        .shadows-into light-regular {
            font-family: "Shadows Into Light", cursive;
            font-weight: 400;
            font-style: normal;
        }
    """
    )
    # Register some fonts that we might need later on.
    app.add_static_files("/fonts", str(Path(__file__).parent / "fonts"))

    # We need to register things to do when the application actually starts up but
    # not when a new user connects to the website.
    app.on_startup(startup)

    # At last we can actually start ROBIN
    ui.run(port=port, reload=reload, title="ROBIN", favicon=iconfile, on_air=False)


def configure(ctx, param, filename):
    cfg = ConfigParser()
    cfg.read(filename)
    try:
        options = dict(cfg["options"])
    except KeyError:
        options = {}
    ctx.default_map = options


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
@click.option(
    "--showerrors/--noerrors",
    default=False,
    help="If set, will display all errors in running R.",
)
@click.option(
    "--sequencing_summary",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, resolve_path=True, path_type=Path
    ),
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
    if sequencing_summary:
        sequencing_summary = click.format_filename(sequencing_summary)
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
            reference=click.format_filename(reference),
        )
    else:
        # Handle the case when --browse is not set
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if output is None:
            click.echo("Output is required when --browse is not set.")
            sys.exit(1)

        run_class(
            port=port,
            reload=False,
            threads=threads,
            simtime=simtime,
            watchfolder=click.format_filename(watchfolder),
            output=click.format_filename(output),
            sequencing_summary=sequencing_summary,
            target_panel=target_panel,
            showerrors=showerrors,
            browse=browse,
            exclude=exclude,
            reference=click.format_filename(reference),
            # minknow_connection=minknow_connection,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    package_run()
