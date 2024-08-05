import uuid
import click
import os
import sys
import signal
import logging

from typing import Optional, List
from pathlib import Path
from nicegui import ui, app, core, observables
import nicegui.air
from robin import theme
from robin import images
from configparser import ConfigParser

from robin.brain_class import BrainMeth
from robin.minknow_info.minknow_panel import MinKNOWFish

from robin.__about__ import __version__


DEFAULT_CFG: str = "config.ini"
UNIQUE_ID: str = str(uuid.uuid4())


def setup_logging(level: str, log_file: Path) -> None:
    """
    Configure global logging level and log file.

    :param level: Logging level as a string (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    :param log_file: Path to the log file.
    """
    numeric_level = getattr(logging, level.upper(), None)
    if numeric_level is None:
        raise ValueError(f"Invalid log level: {level}")

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    try:
        logging.basicConfig(
            level=numeric_level,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
        )
        logging.info(f"Logging configured to level: {level}")
        logging.info(f"Logging to file: {log_file}")
        logging.debug("Debug logging enabled.")
    except PermissionError:
        print(f"Error: Unable to write to log file {log_file}. Check permissions.")
        sys.exit(1)
    except Exception as e:
        print(f"Error setting up logging: {str(e)}")
        sys.exit(1)


def clean_up_handler(thingtokill):
    del thingtokill


@ui.page("/", response_timeout=30)
async def index() -> None:
    """
    Main index page for the ROBIN site.
    Sets up an instance of Methnice for the splash screen and navigation.
    """
    GUI = Methnice(
        force_sampleid=app.storage.general[UNIQUE_ID]["force_sampleid"],
        kit=app.storage.general[UNIQUE_ID]["kit"],
        centreID=app.storage.general[UNIQUE_ID]["centreID"],
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

    def use_on_air():
        """
        Enable or disable remote access based on the value of the general argument.

        Example:
            >>> args = events.ValueChangeEventArguments(value=True)
            >>> use_on_air(args)
            None
        """
        if "use_on_air" in app.storage.general.keys():
            if app.storage.general["use_on_air"]:
                if any("on-air" in url for url in app.urls):
                    pass
                else:
                    if core.air is None:
                        core.ai = nicegui.air.Air("")
                    nicegui.air.connect()
            else:
                nicegui.air.disconnect()

    ui.timer(10, use_on_air)
    ui.context.client.on_disconnect(lambda: clean_up_handler(GUI))
    await GUI.splash_screen()


@ui.page("/live/{sample_id}", response_timeout=30)
async def live_sample(sample_id: str) -> None:
    """
    Page for live data interaction with a specific sample ID.
    Sets up an instance of Methnice for live data visualization.
    """
    GUI = Methnice(
        force_sampleid=app.storage.general[UNIQUE_ID]["force_sampleid"],
        kit=app.storage.general[UNIQUE_ID]["kit"],
        centreID=app.storage.general[UNIQUE_ID]["centreID"],
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
        sample_id=sample_id,
    )
    GUI.setup()
    ui.context.client.on_disconnect(lambda: clean_up_handler(GUI))
    await GUI.index_page()


@ui.page("/live", response_timeout=30)
async def live() -> None:
    """
    Page for live data interaction.
    Sets up an instance of Methnice for live data visualization.
    """
    GUI = Methnice(
        force_sampleid=app.storage.general[UNIQUE_ID]["force_sampleid"],
        kit=app.storage.general[UNIQUE_ID]["kit"],
        centreID=app.storage.general[UNIQUE_ID]["centreID"],
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
    ui.context.client.on_disconnect(lambda: clean_up_handler(GUI))
    await GUI.index_page()


@ui.page("/browse", response_timeout=60)
async def test() -> None:
    """
    Placeholder for browsing historic data.
    """
    GUI_browse = Methnice(
        threads=1,
        force_sampleid=None,
        kit=None,
        centreID=None,
        simtime=False,
        watchfolder=None,
        output=app.storage.general[UNIQUE_ID]["output"],
        sequencing_summary=None,
        target_panel=app.storage.general[UNIQUE_ID]["target_panel"],
        showerrors=app.storage.general[UNIQUE_ID]["showerrors"],
        browse=True,
        exclude=app.storage.general[UNIQUE_ID]["exclude"],
        reference=app.storage.general[UNIQUE_ID]["reference"],
        unique_id=UNIQUE_ID,
    )
    GUI_browse.setup()
    ui.context.client.on_disconnect(lambda: clean_up_handler(GUI_browse))
    await GUI_browse.browse_page()


def clean_up() -> None:
    """
    Clean up stored data at the end of a run.
    """
    logging.info("App shutdown detected.")
    logging.info("Resetting Storage.")
    app.storage.general.pop(UNIQUE_ID)
    app.shutdown()


async def startup() -> None:
    """
    Start data processing in the main application loop.
    """
    logging.info(f"Setting up {UNIQUE_ID}.")
    MAINPAGE = Methnice(
        force_sampleid=app.storage.general[UNIQUE_ID]["force_sampleid"],
        kit=app.storage.general[UNIQUE_ID]["kit"],
        centreID=app.storage.general[UNIQUE_ID]["centreID"],
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

    def handler(*args):
        clean_up()

    signal.signal(signal.SIGINT, handler=handler)


class Methnice:
    """
    This class handles configuration of all the pages and page contents.
    """

    def __init__(
        self,
        force_sampleid: Optional[str],
        kit: Optional[str],
        centreID: Optional[str],
        threads: int,
        simtime: bool,
        watchfolder: Optional[Path],
        output: Optional[Path],
        sequencing_summary: Optional[Path],
        target_panel: str,
        showerrors: bool,
        browse: bool,
        exclude: List[str],
        reference: Path,
        unique_id: str,
        sample_id: Optional[str] = None,
    ):
        self.force_sampleid = force_sampleid
        self.kit = kit
        self.centreID = centreID
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
        if sample_id:
            self.sample_id = sample_id
        else:
            self.sample_id = None

    @property
    def watchfolder(self) -> Optional[Path]:
        """
        Getter for watchfolder property.
        """
        return self._watchfolder

    @watchfolder.setter
    def watchfolder(self, value: Optional[Path]) -> None:
        """
        Setter for watchfolder property.
        Triggers on_watchfolder_changed when watchfolder is set.
        """
        self._watchfolder = value
        self.on_watchfolder_changed()

    def _worker(self) -> None:
        """
        Worker method to handle background tasks.
        """
        if hasattr(self.minknow_connection, "minKNOW_display"):
            if self.minknow_connection.minKNOW_display.watchfolder != self.watchfolder:
                self.watchfolder = self.minknow_connection.minKNOW_display.watchfolder

    async def on_watchfolder_changed(self) -> None:
        """
        Async method to handle changes to the watchfolder.
        """
        if self.watchfolder:
            logging.info(f"watchfolder value has been changed! {self.watchfolder}")
            try:
                await self.robin.add_watchfolder(self.watchfolder)
            except Exception as e:
                logging.error(f"Error adding watchfolder: {str(e)}")

    def setup(self) -> None:
        """
        Setup method for initializing BrainMeth.
        """
        try:
            self.robin = BrainMeth(
                force_sampleid=self.force_sampleid,
                kit=self.kit,
                centreID=self.centreID,
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
        except Exception as e:
            logging.error(f"Error initializing BrainMeth: {str(e)}")
            raise

    async def start_analysis(self) -> None:
        """
        Async method to start analysis.
        """
        try:
            await self.robin.start_background()
        except Exception as e:
            logging.error(f"Error starting analysis: {str(e)}")
            # Consider how to handle this error (e.g., show an error message to the user)

    async def browse_page(self) -> None:
        """
        Async method for rendering the browse page for historic data.

        This method sets up a simplified version of the main page, focused on
        browsing historic data. It includes:
        1. Creating the main frame with the ROBIN title.
        2. Setting up a BrainMeth instance for browsing.
        3. Rendering the UI for the browse functionality.

        The browse page does not include MinKNOW connection or real-time analysis features.

        Raises:
            Any exceptions raised during UI rendering are not caught here and should be
            handled by the caller.
        """
        with theme.frame(
            "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
            smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
        ):
            self.analysis_tab_pane = ui.row().classes("w-full")
            self.robin_browse = BrainMeth(
                force_sampleid=self.force_sampleid,
                kit=self.kit,
                centreID=self.centreID,
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

            with self.analysis_tab_pane:
                await self.robin_browse.render_ui()

    async def splash_screen(self) -> None:
        """
        Async method for rendering the splash screen.
        """
        with theme.frame(
            "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
            smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
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
                            "ROBIN_logo_small.png",
                        )
                    ).classes("rounded-full w-16 h-16 ml-4")
                with ui.button(on_click=lambda: ui.navigate.to("/browse")).props(
                    "color=green"
                ):
                    ui.label("Browse Historic Data")
                    ui.image(
                        os.path.join(
                            os.path.dirname(os.path.abspath(images.__file__)),
                            "ROBIN_logo_small.png",
                        )
                    ).classes("rounded-full w-16 h-16 ml-4")

    async def index_page(self) -> None:
        """
        Async method for rendering the main index page of the application.

        This method sets up the primary user interface, including:
        1. Creating the main frame with the ROBIN title.
        2. Setting up the MinKNOW connection if necessary.
        3. Creating a splitter layout with tabs for Analysis and MinKNOW data.
        4. Rendering the analysis tab content.
        5. Setting up and connecting to the MinKNOW panel if available.

        The layout and content of the page depend on whether a watchfolder is set
        and if the application is in browse mode.

        Raises:
            Any exceptions raised during UI rendering or MinKNOW connection are not
            caught here and should be handled by the caller.
        """
        try:
            with theme.frame(
                "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
                smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
            ):
                await ui.context.client.connected()
                with ui.column().classes("w-full"):
                    # Set up MinKNOW connection if necessary
                    if self.watchfolder is None and not self.browse:
                        self.minknow_connection = MinKNOWFish(
                            kit=self.kit,
                            reference=self.reference,
                            centreID=self.centreID,
                            basecall_config = self.basecall_config,
                            experiment_duration = self.experiment_duration,
                            bed_file = self.bed_file,
                        )
                    else:
                        self.minknow_connection = None
                    
                    # Create splitter layout if MinKNOW connection is available
                    if self.minknow_connection:
                        with ui.splitter(value=10).classes("w-full h-full") as splitter:
                            with splitter.before:
                                with ui.tabs().props("vertical").classes("w-full") as tabs:
                                    analysis = ui.tab("Analysis", icon="analytics")
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
                                    self.minknow_tab_pane = ui.tab_panel(minknow)
                                    with self.minknow_tab_pane:
                                        with ui.row():
                                            ui.icon("set_meal", color="primary").classes(
                                                "text-h4"
                                            )
                                            ui.label("MinKNOW Data").classes("text-h4")
                    else:
                        self.analysis_tab_pane = ui.row().classes("w-full")

                # Set up MinKNOW panel if connection is available
                if self.minknow_connection:
                    with self.minknow_tab_pane:
                        self.minknow_connection.setup_ui()
                        try:
                            await self.minknow_connection.auto_connect()
                        except Exception as e:
                            logging.error(f"Error connecting to MinKNOW: {str(e)}")
                            ui.notify(f"Failed to connect to MinKNOW: {str(e)}", color="negative")
                        ui.label().bind_text_from(
                            self.minknow_connection,
                            "connection_ip",
                            backward=lambda n: f"Connected to: {n}",
                        )

                # Render the analysis UI
                with self.analysis_tab_pane:
                    try:
                        await self.robin.render_ui(sample_id=self.sample_id)
                    except Exception as e:
                        logging.error(f"Error rendering analysis UI: {str(e)}")
                        ui.notify("Failed to render analysis UI", color="negative")

        except Exception as e:
            logging.error(f"Error rendering index page: {str(e)}")
            ui.notify("An error occurred while loading the page", color="negative")


def run_class(
    port: int,
    force_sampleid: Optional[str],
    kit: str,
    centreID: str,
    reload: bool,
    threads: int,
    simtime: bool,
    watchfolder: Optional[Path],
    output: Optional[Path],
    sequencing_summary: Optional[Path],
    target_panel: str,
    showerrors: bool,
    browse: bool,
    exclude: List[str],
    reference: Path,
) -> None:
    """
    Set up and run the ROBIN application.

    :param port: Port for the GUI.
    :param reload: Boolean indicating if the server should reload on changes.
    :param threads: Number of threads available.
    :param simtime: Boolean indicating if simulation mode is enabled.
    :param watchfolder: Path to the watchfolder.
    :param output: Path to the output directory.
    :param sequencing_summary: Path to the sequencing summary file. If provided, timestamps will be taken from this file.
    :param target_panel: Analysis gene panel.
    :param showerrors: Boolean indicating if errors should be displayed.
    :param browse: Boolean indicating if browse mode is enabled.
    :param exclude: List of analysis types to exclude.
    :param reference: Path to the reference genome and index.
    """
    try:
        app.storage.general.clear()
    except Exception as e:
        logging.error(f"Error clearing app storage: {str(e)}")

    try:
        iconfile = os.path.join(
            os.path.dirname(os.path.abspath(images.__file__)), "favicon.ico"
        )
        if not os.path.exists(iconfile):
            logging.warning(f"Favicon file not found: {iconfile}")
    except Exception as e:
        logging.error(f"Error locating favicon: {str(e)}")
        iconfile = None

    app.storage.general[UNIQUE_ID] = {
        "threads": threads,
        "force_sampleid": force_sampleid,
        "kit": kit,
        "centreID": centreID,
        "simtime": simtime,
        "watchfolder": watchfolder,
        "output": output,
        "sequencing_summary": sequencing_summary,
        "target_panel": target_panel,
        "showerrors": showerrors,
        "browse": browse,
        "exclude": exclude,
        "reference": reference,
    }
    app.storage.general[UNIQUE_ID]["samples"] = {}
    app.storage.general[UNIQUE_ID]["sample_list"] = observables.ObservableList([])

    ui.add_css(
        """
        .shadows-into light-regular {
            font-family: "Shadows Into Light", cursive;
            font-weight: 400;
            font-style: normal;
        }
    """
    )

    try:
        app.add_static_files("/fonts", str(Path(__file__).parent / "fonts"))
    except Exception as e:
        logging.error(f"Error adding static files: {str(e)}")

    app.on_startup(startup)

    try:
        ui.run(
            port=port,
            reload=reload,
            title="ROBIN",
            favicon=iconfile,
            on_air=False,
            show=False,
            storage_secret="UNIQUE_ID",
            reconnect_timeout=60,
        )
    except Exception as e:
        logging.error(f"Error running the application: {str(e)}")
        sys.exit(1)


def configure(ctx: click.Context, param: click.Parameter, filename: str) -> None:
    """
    Configure the application based on the provided INI file.

    :param ctx: Click context.
    :param param: Click parameter.
    :param filename: Path to the configuration file.
    """
    cfg = ConfigParser()
    try:
        cfg.read(filename)
        options = dict(cfg["options"])
    except FileNotFoundError:
        logging.warning(f"Configuration file not found: {filename}")
        options = {}
    except KeyError:
        logging.warning("No 'options' section found in the configuration file")
        options = {}
    except Exception as e:
        logging.error(f"Error reading configuration file: {str(e)}")
        options = {}
    
    ctx.default_map = options


@click.command()
@click.version_option(version=__version__, message="%(version)s")
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
    "--force_sampleid",
    default=None,
    help="Force a specific sampleID.",
    required=False,
    type=str,
)
@click.option(
    "--kit",
    default="SQK-RAD114",
    help="Specify sequencing kit.",
    required=False,
    type=str,
)
@click.option(
    "--centreID",
    help="Provide an identifier to be used in the experiment name field.",
    required=True,
    type=str,
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
    "--log-file",
    type=click.Path(dir_okay=False, writable=True, resolve_path=True, path_type=Path),
    default="ROBIN.log",
    help="Path to the log file.",
    show_default=True,
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
    port: int,
    force_sampleid: Optional[str],
    kit: str,
    centreid: str,
    threads: int,
    log_level: str,
    log_file: Path,
    simtime: bool,
    showerrors: bool,
    sequencing_summary: Optional[Path],
    target_panel: str,
    watchfolder: Optional[Path],
    output: Optional[Path],
    browse: bool,
    exclude: List[str],
    reference: Path,
) -> None:
    """
    Entrypoint for when GUI is launched directly.
    """
    setup_logging(log_level, log_file)

    if sequencing_summary:
        sequencing_summary = click.format_filename(sequencing_summary)
    if browse:
        logging.info("Browse mode is enabled. Watchfolder and output are not required.")
        run_class(
            port=port,
            force_sampleid=force_sampleid,
            kit=kit,
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
        logging.info(f"Watchfolder: {watchfolder}, Output: {output}")
        if output is None:
            logging.error("Output is required when --browse is not set.")
            sys.exit(1)

        run_class(
            port=port,
            force_sampleid=force_sampleid,
            kit=kit,
            centreID=centreid,
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
        )


if __name__ in {"__main__", "__mp_main__"}:
    logging.info("GUI launched by auto-reload function.")
    package_run()
