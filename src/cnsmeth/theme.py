"""
This module contains the theme for the whole application.
"""

from contextlib import contextmanager

from nicegui import ui, app, events, core
import nicegui.air

from pathlib import Path

from cnsmeth import images

import os

# Define the path to the image file used in the header and footer
IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "MethBrain_small.png"
)

# Read the HTML content for the header
HEADER_HTML = (Path(__file__).parent / "static" / "header.html").read_text()

# Read the CSS styles for the application
STYLE_CSS = (Path(__file__).parent / "static" / "styles.css").read_text()

@contextmanager
def frame(navtitle: str):
    """
    Context manager to create a custom page frame with consistent styling and behavior across all pages.
    
    Args:
        navtitle (str): The title to display in the navigation header.
    
    Yields:
        None
    """
    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(HEADER_HTML + f"<style>{STYLE_CSS}</style>")
    
    # Create a persistent dialog for quitting the app
    with ui.dialog().props("persistent") as quitdialog, ui.card():
        ui.label("Quitting the app will stop running methylation analysis. Are you sure?")
        ui.label("If you want to keep analysis running, click Cancel.")
        ui.label("You can safely close this window and analysis will keep running in the background.")
        ui.button("Cancel", on_click=quitdialog.close).props("outline").classes("shadow-lg")
        ui.button("Really Quit", icon="logout", on_click=cleanup_and_exit).props("outline").classes("shadow-lg")
    
    # Create a header with navigation title and menu
    with ui.header(elevated=True).classes("items-center duration-200 p-0 px-4 no-wrap"):
        with ui.grid(columns=2).style("width: 100%"):
            ui.html(navtitle).classes("shadows-into").style("color: #FFFFFF; font-size: 200%; font-weight: 150").tailwind("drop-shadow", "font-bold")
            with ui.row().classes("ml-auto align-top"):
                with ui.button(icon="menu"):
                    with ui.menu() as menu:
                        ui.menu_item("Home", lambda: ui.navigate.to("/"))
                        ui.menu_item("Live Data", lambda: ui.navigate.to("/live"))
                        ui.menu_item("Browse Historic Data", lambda: ui.navigate.to("/browse"))
                        ui.separator()
                        ui.separator()
                        ui.menu_item("Close", menu.close)
                ui.image(IMAGEFILE).style("width: 50px")

    # Create a footer with useful information and quit button
    with ui.footer().style("background-color: #4F9153"):
        with ui.dialog() as dialog, ui.card():
            ui.label("Useful Information.").tailwind("text-2xl font-bold font-italic drop-shadow")
            ui.separator()
            ui.link("Code on GitHub", "https://github.com/looselab/cnsmeth")
            ui.link("Rapid CNS2 Paper", "https://link.springer.com/article/10.1007/s00401-022-02415-6")
            ui.link("Sturgeon Classifier", "https://www.nature.com/articles/s41586-023-06615-2")
            ui.link("Protocol", "https://www.protocols.io/view/intra-operative-nanopore-sequencing-to-classify-br-c65qzg5w")
            ui.link("Oxford Nanopore", "https://nanoporetech.com/")
            ui.link("epi2me labs", "https://labs.epi2me.io/")
            ui.link("Looselab", "https://looselab.github.io/")
            ui.button("Close", on_click=dialog.close)
        ui.image(IMAGEFILE).style("width: 30px")
        ui.colors(primary="#555")
        ui.button("More Information", on_click=dialog.open)
        ui.button("Quit", icon="logout", on_click=quitdialog.open)
        ui.label().bind_text_from(app, "urls", backward=lambda n: f"Available urls: {n}")
        ui.label("Some aspects of this application are ©Looselab - all analyses provided for research use only.").tailwind("text-sm font-italic")
    
    yield

def cleanup_and_exit():
    """
    Handle any necessary cleanup operations before exiting the application and then shut down the application.
    
    Returns:
        None
    
    Example:
        >>> cleanup_and_exit()
        None
    """
    app.shutdown()

def dark_mode(event: events.ValueChangeEventArguments):
    """
    Toggle dark mode for the interface based on the value of the event argument.
    
    Args:
        event (events.ValueChangeEventArguments): The event argument containing the value for dark mode toggle.
    
    Returns:
        None
    
    Example:
        >>> event = events.ValueChangeEventArguments(value=True)
        >>> dark_mode(event)
        None
    """
    if event.value:
        ui.dark_mode().enable()
    else:
        ui.dark_mode().disable()

def use_on_air(args: events.ValueChangeEventArguments):
    """
    Enable or disable remote access based on the value of the event argument.
    
    Args:
        args (events.ValueChangeEventArguments): The event argument containing the value for remote access toggle.
    
    Returns:
        None
    
    Example:
        >>> args = events.ValueChangeEventArguments(value=True)
        >>> use_on_air(args)
        None
    """
    if args.value:
        if core.air is None:
            core.air = nicegui.air.Air("")
        nicegui.air.connect()
    else:
        nicegui.air.disconnect()