"""
This module contains the theme for the whole application.
"""

from contextlib import contextmanager

from nicegui import ui, app
from nicegui.events import ValueChangeEventArguments

from cnsmeth import images

import os

IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "MethBrain_small.png"
)


@contextmanager
def frame(navtitle: str, myconnection):
    """Custom page frame to share the same styling and behavior across all pages"""
    with ui.dialog().props('persistent') as quitdialog, ui.card():
        ui.label("Quitting the app will stop running methylation analysis. Are you sure?")
        ui.label("If you want to keep analysis running, click Cancel.")
        ui.label("You can safely close this window and analysis will keep running in the background.")
        ui.button("Cancel", on_click=quitdialog.close).props("outline").classes(
            "shadow-lg"
        )
        ui.button(
            "Really Quit", icon="logout", on_click=cleanup_and_exit
        ).props("outline").classes("shadow-lg")
    with ui.header(fixed=True).classes(replace="row items-center p-2").style(
        "box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1)"
    ):
        with ui.grid(columns=2).style("width: 100%"):
            ui.label(navtitle).tailwind("text-2xl font-bold font-italic drop-shadow")
            with ui.row().classes("ml-auto"):
                ui.switch("Dark Mode", on_change=dark_mode).classes(
                    "ml-4 bg-transparent"
                ).props('color="black"')
                ui.image(IMAGEFILE).style("width: 50px")
    with ui.column().classes("w-full"):
        yield
    with ui.footer():
        with ui.dialog() as dialog, ui.card():
            ui.label("Useful Information.").tailwind(
                "text-2xl font-bold font-italic drop-shadow"
            )
            ui.separator()
            ui.link("Code on GitHub", "https://github.com/looselab/cnsmeth")
            ui.link(
                "Rapid CNS2 Paper",
                "https://link.springer.com/article/10.1007/s00401-022-02415-6",
            )
            ui.link(
                "Sturgeon Classifier",
                "https://www.nature.com/articles/s41586-023-06615-2",
            )
            ui.link(
                "Protocol",
                "https://www.protocols.io/view/intra-operative-nanopore-sequencing-to-classify-br-c65qzg5w",
            )
            ui.link("Oxford Nanopore", "https://nanoporetech.com/")
            ui.link("epi2me labs", "https://labs.epi2me.io/")
            ui.link("Looselab", "https://looselab.github.io/")
            ui.button("Close", on_click=dialog.close)
        ui.image(IMAGEFILE).style("width: 30px")
        ui.button("More Information", on_click=dialog.open)
        ui.button(
            "Quit", icon="logout", on_click=quitdialog.open
        )  # .classes('ml-4')#.props('outline') #.classes('shadow-lg')
        ui.label(
            "Some aspects of this application are ©Looselab - all analyses provided for research use only."
        ).tailwind("text-sm font-italic")


def cleanup_and_exit():
    """
    This function should handle anything that needs to be closed down before the app exits.
    :return: None
    """
    app.shutdown()


def dark_mode(event: ValueChangeEventArguments):
    """
    This function handles toggling dark mode for the interface.
    :param event:
    :return: None
    """
    if event.value:
        ui.dark_mode().enable()
    else:
        ui.dark_mode().disable()
