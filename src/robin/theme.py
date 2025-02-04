"""
Module: theme

This module defines the theme and layout for the entire application, ensuring a consistent look and feel across all pages.

It includes:

- A context manager `frame` to create a custom page frame with navigation, header, and footer.
- Utility functions to handle dark mode and remote access toggling.
- Paths to external resources like images, HTML, and CSS files for styling.

Functions:

- frame(navtitle: str): Context manager for creating a consistent page layout with header, footer, and navigation.
- cleanup_and_exit(): Handles cleanup operations before shutting down the application.
- dark_mode(event: events.ValueChangeEventArguments): Toggles dark mode based on the event argument.
- use_on_air(args: events.ValueChangeEventArguments): Toggles remote access based on the event argument.

Constants:

- IMAGEFILE: Path to the image file used in the header and footer.
- HEADER_HTML: HTML content for the header.
- STYLE_CSS: CSS styles for the application.

External Dependencies:

- contextlib.contextmanager
- nicegui (ui, app, events, core, air)
- pathlib.Path
- robin.images
- os
- psutil
- platform
"""

from contextlib import contextmanager
from packaging import version
import requests
import asyncio

from nicegui import ui, app, events, core
import nicegui.air

from pathlib import Path

from robin import images
from robin import __about__

import os
import psutil
import platform

async def check_version():
    """
    Check the current version against the remote version on GitHub.
    Shows a notification or dialog to the user about their version status.
    """
    # Check if version has already been checked in this session
    if app.storage.tab.get("version_checked", False):
        return

    try:
        # Get the remote version from GitHub
        response = requests.get('https://raw.githubusercontent.com/LooseLab/ROBIN/main/src/robin/__about__.py')
        response.raise_for_status()
        
        # Extract version from the response text
        remote_version_str = None
        for line in response.text.split('\n'):
            if line.startswith('__version__'):
                remote_version_str = line.split('=')[1].strip().strip('"').strip("'")
                break
        
        if not remote_version_str:
            with ui.dialog() as dialog, ui.card():
                ui.label('Version Check Error').classes('text-h6')
                ui.label('Could not determine remote version. Please check manually.')
                ui.button('OK', on_click=dialog.close)
            dialog.open()
            return

        local_version = version.parse(__about__.__version__)
        remote_version = version.parse(remote_version_str)

        if local_version == remote_version:
            ui.notify('Your ROBIN installation is up to date!', type='positive')
        elif local_version < remote_version:
            with ui.dialog() as dialog, ui.card():
                ui.label('Update Available!').classes('text-h6')
                ui.label(f'Your version: {local_version}')
                ui.label(f'Latest version: {remote_version}')
                ui.label('Would you like to visit the GitHub repository to update?')
                with ui.row():
                    ui.button('Continue with current version', on_click=dialog.close)
                    ui.button('Visit GitHub', on_click=lambda: ui.open('https://github.com/LooseLab/ROBIN')).classes('bg-primary')
            dialog.open()
        else:
            with ui.dialog() as dialog, ui.card():
                ui.label('Development Version').classes('text-h6')
                ui.label(f'You are running a development version ({local_version}).')
                ui.label(f'Latest release: {remote_version}')
                ui.label('This version may be unstable and is only for testing purposes. It is not recommended for production use.')
                ui.label('Please consider using the latest release instead.')
                ui.button('OK', on_click=dialog.close)
            dialog.open()
    
    except requests.RequestException:
        with ui.dialog() as dialog, ui.card():
            ui.label('Connection Error').classes('text-h6')
            ui.label('Could not check for updates.')
            ui.label('Either you are not connected to the internet or you cannot access https://www.github.com/looselab/robin.')
            ui.label('Please manually check for updates.')
            ui.button('OK', on_click=dialog.close)
        dialog.open()
    except Exception as e:
        with ui.dialog() as dialog, ui.card():
            ui.label('Error').classes('text-h6')
            ui.label(f'Error checking version: {str(e)}')
            ui.button('OK', on_click=dialog.close)
        dialog.open()
    
    # Mark version as checked for this session
    app.storage.tab["version_checked"] = True

# Define the path to the image file used in the header and footer
IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
)


MENU_BREAKPOINT = 1200

# Read the HTML content for the header
HEADER_HTML = (Path(__file__).parent / "static" / "header.html").read_text()

# Read the CSS styles for the application
STYLE_CSS = (Path(__file__).parent / "static" / "styles.css").read_text()


@contextmanager
def frame(navtitle: str, batphone=False, smalltitle=None):
    """
    Context manager to create a custom page frame with consistent styling and behavior across all pages.

    Args:
        navtitle (str): The title to display in the navigation header.

    Yields:
        None

    Example:
        >>> with frame("Home"):
        ...     ui.label("Welcome to the Application")
    """
    if batphone:
        navtitle = f"BATMAN & {navtitle}"
    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>'
    )
    ui.add_head_html(HEADER_HTML + f"<style>{STYLE_CSS}</style>")
    ui.add_head_html(
        """
        <script>
        function emitSize() {
            emitEvent('resize', {
                width: document.body.offsetWidth,
                height: document.body.offsetHeight,
            });
        }
        window.onload = emitSize;
        window.onresize = emitSize;
        </script>
    """
    )

    # Create disclaimer dialog that appears on first visit
    async def show_disclaimer():
        await ui.context.client.connected()
        if not app.storage.tab.get("disclaimer_acknowledged", False):
            with ui.dialog().props(
                "persistent"
            ) as disclaimer_dialog, ui.card().classes("w-160"):
                ui.label("DISCLAIMER").classes("text-h5 text-weight-bold q-mb-md")
                ui.label(
                    "This tool and the data generated by it are intended for research use only and should not be used for "
                    "direct diagnostic purposes. The methylation-based classifications and other analyses provided here may "
                    "be considered by neuropathologists as supplementary information in the context of comprehensive "
                    "diagnostic assessment, which should include clinical history, radiological findings, and complete "
                    "histopathological and molecular evaluation. The final interpretation and diagnosis should always be "
                    "made by qualified healthcare professionals based on all available information."
                ).classes("text-body1 q-mb-md")

                if batphone:
                    ui.label("BATMAN Mode").classes("text-h5 text-weight-bold q-mb-md")
                    ui.label(
                        "You are running this tool in BATMAN mode. "
                        "This is a beta version of the tool and may not be fully functional. "
                        "BATMAN means: Breakpoint Adaptive Targeting alongside Methylation Analysis on Nanopore. "
                        "This means that the target regions will be updated in real-time based on detected breakpoints. "
                        "This code only works with ReadFish at this time. "
                    )

                def acknowledge():
                    app.storage.tab["disclaimer_acknowledged"] = True
                    disclaimer_dialog.close()

                ui.button("I agree", on_click=acknowledge).props("color=primary")
            disclaimer_dialog.open()

    ui.timer(0.1, show_disclaimer, once=True)
    
    # Add version check timer
    ui.timer(1.0, check_version, once=True)

    # Create a persistent dialog for quitting the app
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

    # Create a header with navigation title and menu
    header_classes = "items-center duration-200 p-0 px-4 no-wrap"
    if batphone:
        header_classes += " batphone"
    with ui.header(elevated=True).classes(header_classes):
        with ui.grid(columns=2).style("width: 100%"):
            with ui.row().classes(
                f"max-[{MENU_BREAKPOINT}px]:hidden items-center align-left"
            ):  # .classes('items-left m-auto'):
                ui.html(navtitle).classes("shadows-into").style(
                    "font-size: 150%; font-weight: 300"
                ).tailwind("drop-shadow", "font-bold")
            with ui.row().classes(
                f"min-[{MENU_BREAKPOINT+1}px]:hidden items-center align-left"
            ):
                ui.html(smalltitle).style("font-size: 150%; font-weight: 300").tailwind(
                    "drop-shadow", "font-bold"
                )
            with ui.row().classes("ml-auto align-top"):
                with ui.row().classes("items-center m-auto"):
                    ui.label(f"Viewing: {platform.node()}").classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )
                    ui.label("CPU").classes(f"max-[{MENU_BREAKPOINT}px]:hidden")
                    cpu_activity = ui.circular_progress(max=100).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )
                    ui.timer(
                        1.0,
                        lambda: cpu_activity.set_value(
                            f"{psutil.getloadavg()[1] / os.cpu_count() * 100:.1f}"
                        ),
                    )
                    ui.label("RAM").classes(f"max-[{MENU_BREAKPOINT}px]:hidden")
                    ram_utilisation = ui.circular_progress(max=100).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )
                    ui.timer(
                        1.0,
                        lambda: ram_utilisation.set_value(
                            f"{psutil.virtual_memory()[2]:.1f}"
                        ),
                    )
                    with ui.button(icon="menu"):
                        with ui.menu() as menu:
                            ui.menu_item("Home", lambda: ui.navigate.to("/"))
                            ui.menu_item("Live Data", lambda: ui.navigate.to("/live"))
                            ui.menu_item(
                                "Browse Historic Data",
                                lambda: ui.navigate.to("/browse"),
                            )
                            ui.separator()
                            ui.switch("Allow Remote Access").classes(
                                "ml-4 bg-transparent"
                            ).props('color="black"').bind_value(
                                app.storage.general, "use_on_air"
                            )
                            ui.separator()
                            # ui.switch(
                            #    "allow remote access", value=False, on_change=use_on_air
                            # ).classes("ml-4 bg-transparent").props('color="black"')
                            # ui.switch("Dark Mode", on_change=dark_mode).classes(
                            #    "ml-4 bg-transparent"
                            # ).props('color="black"')
                            ui.switch("Dark Mode").classes("ml-4 bg-transparent").props(
                                'color="black"'
                            ).bind_value(app.storage.browser, "dark_mode")
                            ui.dark_mode().bind_value(app.storage.browser, "dark_mode")
                            ui.separator()
                            ui.menu_item("Close", menu.close)
                            ui.button("Quit", icon="logout", on_click=quitdialog.open)
                    ui.image(IMAGEFILE).style("width: 50px")

    # Create a footer with useful information and quit button
    # Create a header with navigation title and menu
    footer_classes = "items-center"
    if batphone:
        footer_classes += " batphone"
    with ui.footer().classes(footer_classes):
        # ui.on('resize', lambda e: print(f'resize: {e.args}'))

        with ui.dialog() as dialog, ui.card():
            ui.label("Links").tailwind("text-2xl font-bold font-italic drop-shadow")
            ui.separator()
            ui.link("Code on GitHub", "https://github.com/looselab/robin")
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
        ui.image(IMAGEFILE).style("width: 40px")
        ui.colors(primary="#555")
        ui.button("Links", on_click=dialog.open)

        with ui.button(icon="info"):
            with ui.menu() as menu:
                ui.label().bind_text_from(
                    app, "urls", backward=lambda n: f"Available urls: {n}"
                )
                ui.label("Version: " + __about__.__version__)
        ui.label(
            "Some aspects of this application are ©Looselab - all analyses provided for research use only."
        ).classes(f"max-[{MENU_BREAKPOINT}px]:hidden").tailwind("text-sm font-italic")
        ui.label("©Looselab").classes(f"min-[{MENU_BREAKPOINT+1}px]:hidden").tailwind(
            "text-sm font-italic"
        )
        ui.label("Not for diagnostic use.").classes(
            f"min-[{MENU_BREAKPOINT+1}px]:hidden"
        ).tailwind("text-sm font-italic")

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


@ui.page("/")
def my_page():
    with frame(
        "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>",
        smalltitle="<strong>R.O.B.I.N</strong>",
    ):
        ui.label("Welcome to the Application")


def main():
    """
    Main function to test the theme by creating a simple page using the frame context manager.

    Example:
        >>> main()
        None
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
    app.add_static_files("/fonts", str(Path(__file__).parent / "fonts"))

    ui.run(storage_secret="robin")


if __name__ in {"__main__", "__mp_main__"}:
    # import doctest
    # doctest.testmod()
    main()
