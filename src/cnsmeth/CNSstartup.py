import click
from cnsmeth.brain_class import BrainMeth
from pathlib import Path
from nicegui import ui, app
import os
import sys
from cnsmeth import images

from cnsmeth import theme


class Methnice:
    def __init__(
        self,
        threads: int,
        simtime: bool,
        watchfolder: Path,
        output: Path,
        sequencing_summary: Path,
        showerrors: bool,
        browse: bool,
    ):
        self.threads = threads
        self.simtime = simtime
        self.watchfolder = watchfolder
        self.output = output
        self.sequencing_summary = sequencing_summary
        self.showerrors = showerrors
        self.browse = browse

    @ui.page('/home')
    def index_page(self) -> None:
        my_connection = None
        with theme.frame("Real Time Brain Tumour Classification", my_connection):
            # my_connection.connect_to_minknow()
            BrainMeth(
                threads=self.threads,
                simtime=self.simtime,
                watchfolder=self.watchfolder,
                output=self.output,
                sequencing_summary=self.sequencing_summary,
                showerrors=self.showerrors,
                browse=self.browse,
            )


def run_class(
    port: int,
    reload: bool,
    threads: int,
    simtime: bool,
    watchfolder: Path,
    output: Path,
    sequencing_summary: Path,
    showerrors: bool,
    browse: bool,
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
        showerrors=showerrors,
        browse=browse,
    )
    app.on_startup(mainpage.index_page)
    ui.run(
        port=port, reload=reload, title="RCBTC", favicon=iconfile
    )  # , native=True, fullscreen=False, window_size=(1200, 1000))


@click.command()
@click.option(
    "--port",
    default=8081,
    help="Port for GUI",
)
@click.option("--threads", default=4, help="Number of threads available.")
@click.option(
    "--simtime",
    default=False,
    help="If set, will simulate the addition of existing files to the pipeline based on read data.",
)
@click.option(
    "--showerrors",
    default=False,
    help="If set, will display all errors in running R.",
)
@click.option(
    "--sequencing_summary",
    default=None,
    help="Path to sequencing summary file. If provided, timestamps will be taken from this file.",
)
@click.option(
    "--browse",
    is_flag=True,
    show_default=True,
    default=False,
    help="Browse Historic Data.",
)
@click.argument(
    "watchfolder",
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
    port, threads, simtime, showerrors, sequencing_summary, watchfolder, output, browse
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
            showerrors=showerrors,
            browse=browse,
        )
        # Your logic for browse mode
    else:
        # Handle the case when --browse is not set
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
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
            showerrors=showerrors,
            browse=browse,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    run_class(port=12398, reload=True)