"""
This module contains the theme for the whole application.
"""

from contextlib import contextmanager

from nicegui import ui, app, events, core
import nicegui.air

from pathlib import Path

from cnsmeth import images

import os

IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "MethBrain_small.png"
)
HEADER_HTML = (Path(__file__).parent / "static" / "header.html").read_text()
STYLE_CSS = (Path(__file__).parent / "static" / "styles.css").read_text()


@contextmanager
def frame(navtitle: str):
    ui.add_head_html(HEADER_HTML + f"<style>{STYLE_CSS}</style>")
    """Custom page frame to share the same styling and behavior across all pages"""
    with ui.dialog().props("persistent") as quitdialog, ui.card():
        ui.label(
            "Quitting the app will stop running methylation analysis. Are you sure?"
        )
        ui.label("If you want to keep analysis running, click Cancel.")
        ui.label(
            "You can safely close this window and analysis will keep running in the background."
        )
        ui.button("Cancel", on_click=quitdialog.close).props("outline").classes(
            "shadow-lg"
        )
        ui.button("Really Quit", icon="logout", on_click=cleanup_and_exit).props(
            "outline"
        ).classes("shadow-lg")
    with ui.header(elevated=True).classes("items-center duration-200 p-0 px-4 no-wrap"):
        # .style('background-color: #4F9153'):

        with ui.grid(columns=2).style("width: 100%"):
            ui.html(navtitle).classes("shadows-into").style(
                "color: #FFFFFF; font-size: 200%; font-weight: 150"
            ).tailwind(
                "drop-shadow", "font-bold"
            )  # .tailwind("text-2xl font-bold font-italic drop-shadow")
            with ui.row().classes("ml-auto align-top"):
                # ui.switch(
                #    "allow remote access", value=False, on_change=use_on_air
                # ).classes("ml-4 bg-transparent object-center").props('color="black"')
                with ui.button(icon="menu"):
                    with ui.menu() as menu:
                        ui.menu_item("Home", lambda: ui.navigate.to("/"))
                        ui.menu_item("Live Data", lambda: ui.navigate.to("/live"))
                        ui.menu_item(
                            "Browse Historic Data", lambda: ui.navigate.to("/browse")
                        )
                        ui.separator()
                        # ToDo: Dark mode doesn't work on every element
                        # ui.switch("Dark Mode", on_change=dark_mode)
                        # ToDo: Air mode needs to be a global level switch - not per page
                        # ui.switch(
                        #    "allow remote access", value=False, on_change=use_on_air
                        # )
                        ui.separator()
                        ui.menu_item("Close", menu.close)
                ui.image(IMAGEFILE).style("width: 50px")

    with ui.footer().style("background-color: #4F9153"):
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
        ui.colors(primary="#555")
        ui.button("More Information", on_click=dialog.open)
        ui.button(
            "Quit", icon="logout", on_click=quitdialog.open
        )  # .classes('ml-4')#.props('outline') #.classes('shadow-lg')

        ui.label().bind_text_from(
            app, "urls", backward=lambda n: f"Available urls: {n}"
        )
        ui.label(
            "Some aspects of this application are ©Looselab - all analyses provided for research use only."
        ).tailwind("text-sm font-italic")
    yield


def cleanup_and_exit():
    """
    This function should handle anything that needs to be closed down before the app exits.
    :return: None
    """
    app.shutdown()


def dark_mode(event: events.ValueChangeEventArguments):
    """
    This function handles toggling dark mode for the interface.
    :param event:
    :return: None
    """
    if event.value:
        ui.dark_mode().enable()
    else:
        ui.dark_mode().disable()


def use_on_air(args: events.ValueChangeEventArguments):
    if args.value:
        if core.air is None:
            core.air = nicegui.air.Air("")
        nicegui.air.connect()
    else:
        nicegui.air.disconnect()
