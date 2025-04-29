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
import logging

from nicegui import ui, app, events, core
import nicegui.air

from pathlib import Path

from robin import images
from robin import __about__
from robin.core.state import state, ProcessState, ProcessType

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
        response = requests.get(
            "https://raw.githubusercontent.com/LooseLab/ROBIN/main/src/robin/__about__.py"
        )
        response.raise_for_status()

        # Extract version from the response text
        remote_version_str = None
        for line in response.text.split("\n"):
            if line.startswith("__version__"):
                remote_version_str = line.split("=")[1].strip().strip('"').strip("'")
                break

        if not remote_version_str:
            with ui.dialog() as dialog, ui.card():
                ui.label("Version Check Error").classes("text-h6")
                ui.label("Could not determine remote version. Please check manually.")
                ui.button("OK", on_click=dialog.close)
            dialog.open()
            return

        local_version = version.parse(__about__.__version__)
        remote_version = version.parse(remote_version_str)

        if local_version == remote_version:
            ui.notify("Your ROBIN installation is up to date!", type="positive")
        elif local_version < remote_version:
            with ui.dialog() as dialog, ui.card():
                ui.label("Update Available!").classes("text-h6")
                ui.label(f"Your version: {local_version}")
                ui.label(f"Latest version: {remote_version}")
                ui.label("Would you like to visit the GitHub repository to update?")
                with ui.row():
                    ui.button("Continue with current version", on_click=dialog.close)
                    ui.button(
                        "Visit GitHub",
                        on_click=lambda: ui.open("https://github.com/LooseLab/ROBIN"),
                    ).classes("bg-primary")
            dialog.open()
        else:
            with ui.dialog() as dialog, ui.card():
                ui.label("Development Version").classes("text-h6")
                ui.label(f"You are running a development version ({local_version}).")
                ui.label(f"Latest release: {remote_version}")
                ui.label(
                    "This version may be unstable and is only for testing purposes. It is not recommended for production use."
                )
                ui.label("Please consider using the latest release instead.")
                ui.button("OK", on_click=dialog.close)
            dialog.open()

    except requests.RequestException:
        with ui.dialog() as dialog, ui.card():
            ui.label("Connection Error").classes("text-h6")
            ui.label("Could not check for updates.")
            ui.label(
                "Either you are not connected to the internet or you cannot access https://www.github.com/looselab/robin."
            )
            ui.label("Please manually check for updates.")
            ui.button("OK", on_click=dialog.close)
        dialog.open()
    except Exception as e:
        with ui.dialog() as dialog, ui.card():
            ui.label("Error").classes("text-h6")
            ui.label(f"Error checking version: {str(e)}")
            ui.button("OK", on_click=dialog.close)
        dialog.open()

    # Mark version as checked for this session
    app.storage.tab["version_checked"] = True


# Define the path to the image file used in the header and footer
IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
)

# Module-level variables
quitdialog = None

MENU_BREAKPOINT = 1200

# Read the HTML content for the header
HEADER_HTML = (Path(__file__).parent / "static" / "header.html").read_text()

# Read the CSS styles for the application
STYLE_CSS = (Path(__file__).parent / "static" / "styles.css").read_text()


def create_activity_monitor():
    """
    Creates an activity monitor component that shows the status of processes using a clean, streamlined layout.
    Uses Quasar's color system for consistent styling.

    Returns:
        ui.card: The activity monitor card component
    """
    with ui.card().classes("w-full shadow-lg rounded-xl") as monitor:
        # Create a data model for process status
        class ProcessStatus:
            def __init__(self):
                self.processes = {}  # Dictionary to store process states
                self.columns = {}  # Dictionary to store column containers
                self.indicators = {}  # Dictionary to store process indicators
                self.is_shutting_down = False

                # Initialize process states
                for process_type in ProcessType:
                    self.processes[process_type] = {}

        status = ProcessStatus()

        # Create the main expansion panel
        with ui.expansion("Process Status").classes(
            "w-full text-sky-600 dark:text-white"
        ).style("font-size: 120%; font-weight: 500") as status_panel:
            # Status key in a compact format
            with ui.row().classes(
                "w-full px-4 pb-4 gap-2 justify-center items-center text-sm"
            ):
                # Using Quasar's color system
                for status_type, color_class, is_spinning in [
                    ("Running", "primary", True),  # Blue
                    ("Waiting", "warning", False),  # Orange
                    ("Starting", "info", True),  # Light Blue
                    ("Error", "negative", False),  # Red
                    ("Finished", "positive", False),  # Green
                ]:
                    with ui.row().classes("items-center gap-1"):
                        progress = ui.circular_progress(size="xs", show_value=False)
                        progress.props(f"color={color_class}")
                        if is_spinning:
                            progress.props("indeterminate")
                        else:
                            progress.value = 100
                        ui.label(status_type).style("font-size: 0.75rem")

            # Create the grid for process types
            grid = ui.grid(columns=4).classes("w-full gap-4")

            # Create columns for each process type
            for process_type in ProcessType:
                with grid:
                    with ui.column().classes("w-full"):
                        ui.label(process_type.name).classes("text-lg font-bold mb-2")
                        status.columns[process_type] = ui.column().classes("gap-1")

        def create_process_row(process_name, process_type):
            """Creates a new process row with bound indicators"""
            with status.columns[process_type]:
                with ui.row().classes("items-center gap-2 min-w-[200px]").props(
                    "no-wrap"
                ):
                    display_name = " ".join(
                        word.capitalize() for word in process_name.split("_")
                    )
                    ui.label(display_name).classes("text-sm")
                    progress = ui.circular_progress(size="xs", show_value=False)
                    status.indicators[process_name] = progress
                    return progress

        def get_process_state_info(process):
            """
            Determines the color and animation state for a process
            Returns: (color, should_spin)
            """
            try:
                if process in state.process_states:
                    process_state = state.get_process_state(process)
                    state_map = {
                        ProcessState.RUNNING: ("primary", True),  # Blue, spinning
                        ProcessState.WAITING_FOR_DATA: (
                            "warning",
                            False,
                        ),  # Orange, static
                        ProcessState.STARTING: ("info", True),  # Light Blue, spinning
                        ProcessState.STOPPING: ("negative", False),  # Red, static
                    }
                    return state_map.get(
                        process_state, ("grey-7", False)
                    )  # Gray for unknown
                elif process in state.finished_processes:
                    return ("positive", False)  # Green for finished
                return ("negative", False)  # Red for error
            except Exception as e:
                logging.error(f"Error getting state for {process}: {str(e)}")
                return ("negative", False)  # Red for error

        async def update_process_status():
            """Updates the process status display using binding"""
            # Check for shutdown
            if state.shutdown_event or status.is_shutting_down:
                if not status.is_shutting_down:
                    status.is_shutting_down = True
                    logging.info("Status monitor detected shutdown, cleaning up...")
                    # Set all indicators to stopped state
                    for indicator in status.indicators.values():
                        try:
                            # Explicitly remove indeterminate state and set color
                            indicator.props("indeterminate=false color=grey-7")
                            indicator.value = 100
                        except Exception:
                            pass  # Ignore errors during shutdown
                return

            try:
                # Get current processes
                current_processes = set(state.process_states.keys())
                finished_processes = set(state.finished_processes)
                all_processes = current_processes | finished_processes

                # Remove indicators for processes that no longer exist
                removed_processes = set(status.indicators.keys()) - all_processes
                for process in removed_processes:
                    if process in status.indicators:
                        try:
                            status.indicators[process].delete()
                            del status.indicators[process]
                        except Exception:
                            pass  # Ignore cleanup errors

                # Update or create indicators for current processes
                for process in all_processes:
                    try:
                        process_type = state.get_process_type(process)
                        if not process_type or process_type not in ProcessType:
                            continue

                        # Create new indicator if needed
                        if process not in status.indicators:
                            status.indicators[process] = create_process_row(
                                process, process_type
                            )

                        # Get process state information
                        color, should_spin = get_process_state_info(process)

                        # Update indicator
                        indicator = status.indicators[process]

                        if should_spin:
                            # Set spinning state with color
                            indicator.props(f"indeterminate color={color}")
                            indicator.value = None
                        else:
                            # Explicitly remove indeterminate state and set color
                            indicator.props(f"indeterminate=false color={color}")
                            indicator.value = 100

                    except Exception as e:
                        if (
                            not status.is_shutting_down
                        ):  # Only log errors if not shutting down
                            logging.error(f"Error updating process {process}: {str(e)}")

            except Exception as e:
                if not status.is_shutting_down:  # Only log errors if not shutting down
                    logging.error(f"Error in update_process_status: {str(e)}")
                    import traceback

                    logging.error(traceback.format_exc())

        # Update status every 0.2 seconds
        ui.timer(0.2, update_process_status)

    return monitor


@contextmanager
def frame(navtitle: str, batphone=False, smalltitle=None):
    """
    Context manager to create a custom page frame with consistent styling and behavior across all pages.

    Args:
        navtitle (str): The title to display in the navigation header.
        batphone (bool): Whether to show the BATMAN mode title.
        smalltitle (str): The title to display on small screens.

    Yields:
        None
    """
    global quitdialog
    if batphone:
        navtitle = f"BATMAN & {navtitle}"

    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@3.2.0/dist/igv.min.js"></script>'
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
    quitdialog = ui.dialog().props("persistent")

    async def quit_app():
        quitdialog.close()
        await cleanup_and_exit()

    with quitdialog, ui.card():
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
        ui.button("Really Quit", icon="logout", on_click=quit_app).props(
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
            ):
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
                    ui.label("RAM").classes(f"max-[{MENU_BREAKPOINT}px]:hidden")
                    ram_utilisation = ui.circular_progress(max=100).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )

                    # Create a data model for system metrics
                    class SystemMetrics:
                        def __init__(self):
                            self.cpu = 0
                            self.ram = 0

                    metrics = SystemMetrics()

                    # Bind the progress indicators to the model
                    cpu_activity.bind_value_from(metrics, "cpu")
                    ram_utilisation.bind_value_from(metrics, "ram")

                    # Single timer to update both metrics
                    def update_metrics():
                        metrics.cpu = round(
                            psutil.getloadavg()[1] / os.cpu_count() * 100, 1
                        )
                        metrics.ram = round(psutil.virtual_memory()[2], 1)

                    ui.timer(1.0, update_metrics)

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
                            ui.switch("Dark Mode").classes("ml-4 bg-transparent").props(
                                'color="black"'
                            ).bind_value(app.storage.browser, "dark_mode")
                            ui.dark_mode().bind_value(app.storage.browser, "dark_mode")
                            ui.separator()
                            ui.menu_item("Close", menu.close)
                            ui.button("Quit", icon="logout", on_click=quitdialog.open)
                    ui.image(IMAGEFILE).style("width: 50px")

    # Add the process status monitor between main content and footer
    with ui.column().classes(
        "w-full p-4 border-t border-gray-200 dark:border-gray-700"
    ):
        create_activity_monitor()

    # Create a footer with useful information and quit button
    footer_classes = "items-center"
    if batphone:
        footer_classes += " batphone"
    with ui.footer().classes(footer_classes):
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


async def cleanup_and_exit():
    """
    Handle any necessary cleanup operations before exiting the application and then shut down the application.

    Returns:
        None

    Example:
        >>> cleanup_and_exit()
        None
    """
    logging.info("User initiated shutdown via UI")

    # Create and show shutdown modal
    with ui.dialog().props("persistent") as shutdown_dialog, ui.card().classes("w-96"):
        ui.label("Shutting Down").classes("text-h5 text-weight-bold q-mb-md")
        ui.label("ROBIN is shutting down. Please wait while we clean up...").classes(
            "text-body1 q-mb-md"
        )
        with ui.row().classes("w-full justify-center"):
            ui.spinner(size="lg", color="primary")

    # Close the quit dialog if it's open
    if quitdialog:
        quitdialog.close()

    # Show the shutdown dialog
    shutdown_dialog.open()
    logging.info("Shutdown dialog opened")

    # Set shutdown event
    state.shutdown_event = True
    logging.info("Shutdown event set")

    # Wait a moment to ensure the dialog is visible
    await asyncio.sleep(1.0)
    logging.info("Initial wait complete")

    # Perform cleanup
    logging.info("Performing cleanup operations...")
    print("Shutting down ROBIN... from theme.py")
    print(
        "Here we need to do some very graceful shutdown to make sure we don't leave any threads running and we don't leave any files open."
    )

    # Wait for cleanup to complete
    while state.shutdown_event and state.get_running_process_count() > 0:
        await asyncio.sleep(3.0)
    logging.info("Cleanup wait complete")

    # Close the shutdown dialog
    shutdown_dialog.close()
    logging.info("Shutdown dialog closed")

    # Shutdown the application
    logging.info("Application shutdown initiated")
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


# ToDo: The shutdown should be handled in the main.py file. So a gui triggered shutdown needs to set a flag in the main.py file.
#      Then the main.py file can trigger the shutdown and handle the cleanup.
#      The issue is that the nicegui shutdown event is immedately called when the app.shutdown() is called and I need to complete some running code first.


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
