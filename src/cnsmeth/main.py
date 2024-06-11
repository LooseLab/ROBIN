"""
Main Page for Setting Up the ROBIN Real-Time Application

This module initializes and configures the ROBIN real-time application using the NiceGUI framework (https://nicegui.io/). The application ensures that all long-running processes operate in the background, preventing the main thread from being blocked.

The app creates an instance of the `run` class for each specific run, taking arguments from either the command line or a configuration file. Shared arguments are stored in `nicegui.app.storage.general` for accessibility across the application.

Key Components:

1. **Run Class**:

   - Sets up the application.

   - Shares arguments across the app through `app.storage.general`.

2. **Methnice Class**:

   - Processes data in the background.

   - Each user connecting to the site gets individual instances for data visualization without running analysis.

3. **Pages**:

   - `index`: Main page displaying the splash screen and navigation options.

   - `live`: Page for live data interaction.

   - `test`: Placeholder for browsing historic data (to be implemented).

4. **Utility Functions**:

   - `setup_logging(level)`: Configures global logging level.

   - `clean_up()`: Cleans up stored data at the end of a run.

   - `startup()`: Starts data processing in the main application loop.

5. **Configuration**:

   - Uses `click` for command-line interface options.

   - Loads configuration from `config.ini` or specified file.

Constants:

- `DEFAULT_CFG`: Default path to the configuration file.

- `UNIQUE_ID`: Unique identifier for each run, used as an access key in `app.storage.general`.

Dependencies:

- `click`

- `os`

- `sys`

- `signal`

- `uuid`

- `logging`

- `pathlib.Path`

- `nicegui` (ui, app)

- `cnsmeth.theme`

- `cnsmeth.images`

- `configparser.ConfigParser`

Example usage::

    python main.py --config config.ini --port 8081 --threads 4 --log-level INFO

"""

import click
import os
import sys
import signal
import uuid
import logging

from pathlib import Path
from nicegui import ui, app, Client
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


def setup_logging(level):
    """Configure global logging level."""
    numeric_level = getattr(logging, level.upper(), None)
    if numeric_level is None:
        raise ValueError(f"Invalid log level: {level}")

    # Clear any existing handlers to ensure basicConfig can set up the new configuration
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Logging a message after configuration to check the level
    logging.info(f"Logging configured to level: {level}")
    logging.debug("Debug logging enabled.")


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


@ui.page("/live", response_timeout=20)
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
    def clean_up_handler(thingtokill):
        print (f"Killing {thingtokill}")
        del thingtokill
        print(f"Killed")

    ui.context.client.on_disconnect(lambda: clean_up_handler(GUI))
    await GUI.index_page()


@ui.page("/browse", response_timeout=10)
async def test():
    """
    #ToDo: This code is to be implemented.
    The idea is that we will implement the browse functions that allow someone
    to explore previous data even when a live run is in progress.
    """
    GUI_browse = Methnice(
        threads=1,  # app.storage.general[UNIQUE_ID]["threads"],
        simtime=False,  # app.storage.general[UNIQUE_ID]["simtime"],
        watchfolder=None,  # app.storage.general[UNIQUE_ID]["watchfolder"],
        output=None,  # app.storage.general[UNIQUE_ID]["output"],
        sequencing_summary=None,  # app.storage.general[UNIQUE_ID]["sequencing_summary"],
        target_panel=app.storage.general[UNIQUE_ID]["target_panel"],
        showerrors=app.storage.general[UNIQUE_ID]["showerrors"],
        browse=True,  # app.storage.general[UNIQUE_ID]["browse"],
        exclude=app.storage.general[UNIQUE_ID]["exclude"],
        reference=app.storage.general[UNIQUE_ID]["reference"],
        unique_id=UNIQUE_ID,
    )
    GUI_browse.setup()

    def clean_up_handler(thingtokill):
        print (f"Killing {thingtokill}")
        del thingtokill
        print(f"Killed")

    ui.context.client.on_disconnect(lambda: clean_up_handler(GUI_browse))
    await GUI_browse.browse_page()


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
    # loop = asyncio.get_running_loop()
    # loop.set_debug(True)
    # loop.slow_callback_duration = 0.2
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
            self.analysis_tab_pane = ui.row().classes("w-full")
            self.robin_browse = BrainMeth(
                mainuuid=self.MAINID,
                threads=self.threads,
                simtime=self.simtime,
                watchfolder=None,
                output=self.output,
                sequencing_summary=self.sequencing_summary,
                target_panel=self.target_panel,
                showerrors=self.showerrors,
                browse=True,
                exclude=self.exclude,
                minknow_connection=self.minknow_connection,
                reference=self.reference,
            )
            await self.robin_browse.init()

            with self.analysis_tab_pane:
                await self.robin_browse.render_ui()

    async def splash_screen(self) -> None:
        with theme.frame(
            "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>"
        ):
            self.frontpage = ui.card().classes("w-full")
            with self.frontpage:
                ui.label("Welcome to R.O.B.I.N").style(
                    "color: #6E93D6; font-size: 150%; font-weight: 300"
                ).tailwind("drop-shadow", "font-bold")
                ui.label(
                    "This tool enables classification of brain tumours in real time from Oxford Nanopore Data."
                ).style("color: #000000; font-size: 100%; font-weight: 300").tailwind(
                    "drop-shadow", "font-bold"
                )
                with ui.button(on_click=lambda: ui.navigate.to("/live")).props(
                    "color=green"
                ):
                    ui.label("View Live Data")
                    ui.image(
                        os.path.join(
                            os.path.dirname(os.path.abspath(images.__file__)),
                            "favicon.ico",
                        )
                    ).classes("rounded-full w-16 h-16 ml-4")
                with ui.button(on_click=lambda: ui.navigate.to("/browse")).props(
                    "color=green"
                ):
                    ui.label("Browse Historic Data")
                    ui.image(
                        os.path.join(
                            os.path.dirname(os.path.abspath(images.__file__)),
                            "favicon.ico",
                        )
                    ).classes("rounded-full w-16 h-16 ml-4")

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
                await self.robin.render_ui()


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
    #logger = logging.getLogger(__name__)
    #logger.info("Starting ROBIN.")
    #logger.debug("Starting DEBUG ROBIN")
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

    #def clean_up(client: Client):
    #    """
    #    This function will be called when the app is disconnected.
    #    """
    #    print(f"Client disconnect detected. {client}")
    #    print(f"dir{dir(client)}")
    #    print(f"dir{dir(app)}")
        #print("Resetting Storage.")
        # We have seen the app stop. Therefore we want to remove data important for this instance.
        #app.storage.general.pop(UNIQUE_ID)
        #app.shutdown()

    #app.on_disconnect(clean_up)

    # At last we can actually start ROBIN
    ui.run(port=port, reload=reload, title="ROBIN", favicon=iconfile, on_air=False, show=False)


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
    "--log-level",
    default="WARNING",
    type=click.Choice(
        ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        case_sensitive=True,
    ),
    help="Set the logging level (e.g., DEBUG, INFO, WARNING).",
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
    log_level,
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
    #setup_logging(log_level)
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
            watchfolder=(
                click.format_filename(watchfolder) if watchfolder is not None else None
            ),
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
