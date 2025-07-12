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
import subprocess
import importlib.metadata
from datetime import datetime

from nicegui import ui, app, events, core, run
import nicegui.air

from pathlib import Path

from robin import images
from robin import __about__
from robin.core.state import state, ProcessState, ProcessType

import os
import psutil
import platform


def get_version_from_github():
    response = requests.get(
        "https://raw.githubusercontent.com/LooseLab/ROBIN/main/src/robin/__about__.py"
    )
    response.raise_for_status()
    remote_version_str = None
    for line in response.text.split("\n"):
        if line.startswith("__version__"):
            remote_version_str = line.split("=")[1].strip().strip('"').strip("'")
            break
    return remote_version_str


async def check_version():
    """
    Check the current version against the remote version on GitHub.
    Shows a notification or dialog to the user about their version status.
    """
    # Check if version has already been checked in this session
    if app.storage.tab.get("version_checked", False):
        return

    try:
        remote_version_str = await run.io_bound(get_version_from_github)

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

# Global RAM history for background tracking (now only for ROBIN)
ram_history = {
    "robin": [],
    "peak": [],
    "available": [],  # Available system memory
    "free": [],  # Free system memory
    "swap_used": [],  # Used swap memory
    "swap_total": [],  # Total swap memory
    "timestamps": [],
    "max_points": 2880,  # 24 hours at 30s intervals (24 * 60 * 2)
}

# Track the current peak in a mutable container
current_peak = [0]


def get_robin_ram_usage():
    """
    Get accurate RAM usage for ROBIN process and its children on Linux.
    Uses private memory (USS) to avoid double-counting shared memory.
    Falls back to RSS if USS is not available.
    """
    try:
        # Get the main ROBIN process
        main_process = psutil.Process(os.getpid())
        
        # Method 1: Try to get USS using psutil's memory_full_info() (newer versions)
        try:
            total_uss = 0
            total_rss = 0
            
            # Get actual child processes (not just all processes in the tree)
            children = _get_real_children(main_process)
            processes_to_check = [main_process] + children
            
            for proc in processes_to_check:
                try:
                    # Try memory_full_info() first (available in newer psutil versions)
                    mem_info = proc.memory_full_info()
                    if hasattr(mem_info, 'uss'):
                        total_uss += mem_info.uss
                    else:
                        # Fallback to smaps parsing
                        uss = _get_uss_from_smaps(proc.pid)
                        if uss is not None:
                            total_uss += uss
                        else:
                            # Final fallback to RSS
                            total_uss += proc.memory_info().rss
                    
                    # Also track RSS for comparison
                    total_rss += proc.memory_info().rss
                    
                except (psutil.NoSuchProcess, psutil.AccessDenied, FileNotFoundError):
                    continue
                except AttributeError:
                    # memory_full_info() not available, try smaps
                    uss = _get_uss_from_smaps(proc.pid)
                    if uss is not None:
                        total_uss += uss
                    else:
                        total_uss += proc.memory_info().rss
                    total_rss += proc.memory_info().rss
            
            # Log for debugging
            if total_uss / (1024 * 1024 * 1024) > 5:  # If more than 5GB
                logging.debug(f"Memory usage: USS={total_uss/(1024**3):.2f}GB, RSS={total_rss/(1024**3):.2f}GB, Children={len(children)}")
            
            return round(total_uss / (1024 * 1024 * 1024), 2)  # Convert to GB
            
        except Exception as e:
            logging.debug(f"USS method failed, falling back to RSS: {e}")
            
            # Method 2: Fallback to RSS method with better child filtering
            ram_gb = main_process.memory_info().rss / (1024 * 1024 * 1024)
            children = _get_real_children(main_process)
            for child in children:
                try:
                    ram_gb += child.memory_info().rss / (1024 * 1024 * 1024)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            return round(ram_gb, 2)
            
    except Exception as e:
        logging.warning(f"Error getting RAM usage: {e}")
        return 0.0


def debug_memory_measurement():
    """
    Debug function to compare different memory measurement methods.
    This helps verify the accuracy of our memory reporting.
    """
    try:
        main_process = psutil.Process(os.getpid())
        
        print("=== Memory Measurement Debug ===")
        print(f"Process: {main_process.name()} (PID: {main_process.pid})")
        
        # Method 1: RSS (old method)
        rss_bytes = main_process.memory_info().rss
        rss_gb = rss_bytes / (1024 * 1024 * 1024)
        print(f"RSS: {rss_gb:.2f} GB ({rss_bytes:,} bytes)")
        
        # Method 2: USS via memory_full_info() if available
        try:
            mem_info = main_process.memory_full_info()
            if hasattr(mem_info, 'uss'):
                uss_bytes = mem_info.uss
                uss_gb = uss_bytes / (1024 * 1024 * 1024)
                print(f"USS (psutil): {uss_gb:.2f} GB ({uss_bytes:,} bytes)")
            else:
                print("USS not available in psutil")
        except Exception as e:
            print(f"USS (psutil) failed: {e}")
        
        # Method 3: USS via smaps
        uss_smaps = _get_uss_from_smaps(main_process.pid)
        if uss_smaps is not None:
            uss_smaps_gb = uss_smaps / (1024 * 1024 * 1024)
            print(f"USS (smaps): {uss_smaps_gb:.2f} GB ({uss_smaps:,} bytes)")
        else:
            print("USS (smaps) not available")
        
        # Method 4: Our new function
        new_method = get_robin_ram_usage()
        print(f"New method: {new_method:.2f} GB")
        
        # Show ALL child processes (unfiltered)
        all_children = list(main_process.children(recursive=True))
        print(f"\nALL child processes ({len(all_children)}):")
        total_all_rss = 0
        total_all_uss = 0
        
        for i, child in enumerate(all_children):
            try:
                child_rss = child.memory_info().rss / (1024 * 1024 * 1024)
                total_all_rss += child_rss
                
                # Try to get USS for child
                try:
                    child_mem_info = child.memory_full_info()
                    if hasattr(child_mem_info, 'uss'):
                        child_uss = child_mem_info.uss / (1024 * 1024 * 1024)
                        total_all_uss += child_uss
                        print(f"  {i+1}. {child.name()} (PID: {child.pid}): RSS={child_rss:.2f}GB, USS={child_uss:.2f}GB")
                    else:
                        print(f"  {i+1}. {child.name()} (PID: {child.pid}): RSS={child_rss:.2f}GB, USS=unknown")
                except:
                    print(f"  {i+1}. {child.name()} (PID: {child.pid}): RSS={child_rss:.2f}GB, USS=error")
                    
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                print(f"  {i+1}. {child.name()} (PID: {child.pid}): <access denied>")
        
        print(f"\nALL child process totals: RSS={total_all_rss:.2f}GB, USS={total_all_uss:.2f}GB")
        print(f"Main + ALL Children RSS: {rss_gb + total_all_rss:.2f}GB")
        print(f"Main + ALL Children USS: {uss_smaps_gb + total_all_uss:.2f}GB" if uss_smaps is not None else "Main + ALL Children USS: unknown")
        
        # Show FILTERED child processes (only ROBIN-related)
        filtered_children = _get_real_children(main_process)
        print(f"\nFILTERED child processes ({len(filtered_children)}):")
        total_filtered_rss = 0
        total_filtered_uss = 0
        
        for i, child in enumerate(filtered_children):
            try:
                child_rss = child.memory_info().rss / (1024 * 1024 * 1024)
                total_filtered_rss += child_rss
                
                # Try to get USS for child
                try:
                    child_mem_info = child.memory_full_info()
                    if hasattr(child_mem_info, 'uss'):
                        child_uss = child_mem_info.uss / (1024 * 1024 * 1024)
                        total_filtered_uss += child_uss
                        print(f"  {i+1}. {child.name()} (PID: {child.pid}): RSS={child_rss:.2f}GB, USS={child_uss:.2f}GB")
                    else:
                        print(f"  {i+1}. {child.name()} (PID: {child.pid}): RSS={child_rss:.2f}GB, USS=unknown")
                except:
                    print(f"  {i+1}. {child.name()} (PID: {child.pid}): RSS={child_rss:.2f}GB, USS=error")
                    
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                print(f"  {i+1}. {child.name()} (PID: {child.pid}): <access denied>")
        
        print(f"\nFILTERED child process totals: RSS={total_filtered_rss:.2f}GB, USS={total_filtered_uss:.2f}GB")
        print(f"Main + FILTERED Children RSS: {rss_gb + total_filtered_rss:.2f}GB")
        print(f"Main + FILTERED Children USS: {uss_smaps_gb + total_filtered_uss:.2f}GB" if uss_smaps is not None else "Main + FILTERED Children USS: unknown")
        
        # Compare with smem output
        print(f"\nComparison with smem:")
        print(f"smem RSS: 3023.1MB = {3023.1/1024:.2f}GB")
        print(f"Our filtered total RSS: {rss_gb + total_filtered_rss:.2f}GB")
        print(f"Our filtered total USS: {uss_smaps_gb + total_filtered_uss:.2f}GB" if uss_smaps is not None else "Our filtered total USS: unknown")
        
        print("=== End Debug ===")
        
    except Exception as e:
        print(f"Debug failed: {e}")


def _get_uss_from_smaps(pid):
    """
    Get USS (Unique Set Size) from /proc/[pid]/smaps on Linux.
    This is the most accurate way to measure process memory usage.
    """
    try:
        uss = 0
        smaps_path = f"/proc/{pid}/smaps"
        
        if not os.path.exists(smaps_path):
            return None
            
        with open(smaps_path, 'r') as f:
            for line in f:
                if line.startswith('Private_Clean:') or line.startswith('Private_Dirty:'):
                    # Extract the size in KB and convert to bytes
                    size_kb = int(line.split()[1])
                    uss += size_kb * 1024
                    
        return uss
    except (FileNotFoundError, PermissionError, ValueError):
        return None


def _get_real_children(process):
    """
    Get only the actual child processes, filtering out shared libraries and unrelated processes.
    """
    try:
        children = []
        for child in process.children(recursive=True):
            try:
                # Only count processes that are actually related to ROBIN
                if _is_robin_related_process(child):
                    children.append(child)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        return children
    except Exception:
        return []


def _is_robin_related_process(process):
    """
    Check if a process is actually related to ROBIN.
    """
    try:
        # Check if it's a Python process with ROBIN-related command line
        if process.name() in ['python', 'python3', 'python3.9', 'robin']:
            cmdline = process.cmdline()
            if any('robin' in arg.lower() for arg in cmdline):
                return True
            
        # Check if it's a direct child (not grandchild) of the main process
        if process.ppid() == os.getpid():
            return True
            
        # Check if it's a subprocess we spawned for analysis
        if process.name() in ['modkit', 'matkit', 'sturgeon', 'R', 'Rscript']:
            return True
            
        return False
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return False


def reset_peak_memory():
    """Reset the peak memory tracking to current value."""
    current_peak[0] = get_robin_ram_usage()


def collect_ram_usage():
    robin_ram = get_robin_ram_usage()
    current_time = datetime.now().strftime("%H:%M:%S")

    # Get system memory info
    mem = psutil.virtual_memory()
    available_gb = round(mem.available / (1024 * 1024 * 1024), 2)  # Convert to GB
    free_gb = round(mem.free / (1024 * 1024 * 1024), 2)  # Convert to GB

    # Get swap memory info
    swap = psutil.swap_memory()
    swap_used_gb = round(swap.used / (1024 * 1024 * 1024), 2)  # Convert to GB
    swap_total_gb = round(swap.total / (1024 * 1024 * 1024), 2)  # Convert to GB

    # Debug logging for memory discrepancy investigation
    if robin_ram > 10:  # Log if memory usage is suspiciously high
        logging.warning(f"High memory usage detected: {robin_ram:.2f}GB at {current_time}")
        # Run debug function to get detailed breakdown
        debug_memory_measurement()

    ram_history["robin"].append(robin_ram)
    ram_history["available"].append(available_gb)
    ram_history["free"].append(free_gb)
    ram_history["swap_used"].append(swap_used_gb)
    ram_history["swap_total"].append(swap_total_gb)
    ram_history["timestamps"].append(current_time)

    # Update peak - FIXED: Don't reset to 0, properly track the maximum
    current_peak[0] = max(current_peak[0], robin_ram)
    ram_history["peak"].append(current_peak[0])
    # Remove the buggy line: current_peak[0] = 0  # Reset for next interval

    # Keep only the last N points
    max_points = ram_history["max_points"]
    for key in [
        "robin",
        "peak",
        "available",
        "free",
        "swap_used",
        "swap_total",
        "timestamps",
    ]:
        ram_history[key] = ram_history[key][-max_points:]


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
                self.performance_metrics_container = None
                self.ram_history = {
                    "robin": [],
                    "timestamps": [],
                }  # Store RAM usage history
                self.max_history_points = 60  # Store 30 minutes of data (30s * 60)

                # Initialize process states
                for process_type in ProcessType:
                    self.processes[process_type] = {}

        status = ProcessStatus()

        # Create the main expansion panel
        with ui.expansion("Process Status").classes(
            "w-full text-sky-600 dark:text-white"
        ).style("font-size: 120%; font-weight: 500"):
            # Status key in a compact format
            with ui.row().classes(
                "w-full px-4 pb-4 gap-2 justify-center items-center text-sm"
            ):
                # Using Quasar's color system
                for status_type, color_class, is_spinning in [
                    ("Running", "primary", True),  # Blue
                    ("Waiting", "warning", False),  # Orange
                    ("Starting", "info", True),  # Light Blue, spinning
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

            # Add performance metrics UI
            with ui.card().classes("w-full p-4 mb-4"):
                ui.label("System Performance").classes("text-lg font-medium mb-4")
                # Create a container for performance metrics that will be populated by the Brain class
                status.performance_metrics_container = ui.column().classes(
                    "w-full gap-4"
                )

                # Add RAM usage graph using ECharts
                with ui.card().classes("w-full p-4"):
                    ui.label("Memory Usage").classes("text-lg font-medium mb-4")

                    # Initial ECharts options for two lines (current and peak)
                    ram_echart_options = {
                        "backgroundColor": "transparent",
                        "title": {"text": "Memory Usage (24h)", "left": "center"},
                        "tooltip": {"trigger": "axis"},
                        "legend": {
                            "data": [
                                "ROBIN",
                                "Peak",
                                "Available",
                                "Free",
                                "Swap Used",
                                "Swap Total",
                            ],
                            "top": 30,
                        },
                        "xAxis": {
                            "type": "category",
                            "name": "Time",
                            "axisLabel": {"rotate": 45, "formatter": "{value}"},
                            "data": [],
                        },
                        "yAxis": [
                            {
                                "type": "value",
                                "name": "Memory (GB)",
                                "axisLabel": {"formatter": "{value} GB"},
                                "min": 0,
                                "max": "dataMax",
                                "position": "left",
                            },
                            {
                                "type": "value",
                                "name": "Swap (GB)",
                                "axisLabel": {"formatter": "{value} GB"},
                                "min": 0,
                                "max": "dataMax",
                                "position": "right",
                            },
                        ],
                        "dataZoom": [
                            {
                                "type": "slider",
                                "show": True,
                                "xAxisIndex": [0],
                                "start": 0,
                                "end": 100,
                                "height": 20,
                                "bottom": 0,
                            },
                            {
                                "type": "inside",
                                "xAxisIndex": [0],
                                "start": 0,
                                "end": 100,
                            },
                        ],
                        "series": [
                            {
                                "name": "ROBIN",
                                "type": "line",
                                "smooth": True,
                                "data": [],
                                "showSymbol": False,
                                "sampling": "lttb",  # Use LTTB sampling for better performance
                                "yAxisIndex": 0,
                            },
                            {
                                "name": "Peak",
                                "type": "line",
                                "smooth": True,
                                "data": [],
                                "showSymbol": False,
                                "lineStyle": {"type": "dashed"},
                                "sampling": "lttb",  # Use LTTB sampling for better performance
                                "yAxisIndex": 0,
                            },
                            {
                                "name": "Available",
                                "type": "line",
                                "smooth": True,
                                "data": [],
                                "showSymbol": False,
                                "sampling": "lttb",
                                "lineStyle": {"type": "dotted"},
                                "yAxisIndex": 0,
                            },
                            {
                                "name": "Free",
                                "type": "line",
                                "smooth": True,
                                "data": [],
                                "showSymbol": False,
                                "sampling": "lttb",
                                "lineStyle": {"type": "dotted"},
                                "yAxisIndex": 0,
                            },
                            {
                                "name": "Swap Used",
                                "type": "line",
                                "smooth": True,
                                "data": [],
                                "showSymbol": False,
                                "sampling": "lttb",
                                "lineStyle": {"type": "dashed"},
                                "yAxisIndex": 1,
                                "itemStyle": {"color": "#ff6b6b"},  # Red color for swap
                            },
                            {
                                "name": "Swap Total",
                                "type": "line",
                                "smooth": True,
                                "data": [],
                                "showSymbol": False,
                                "sampling": "lttb",
                                "lineStyle": {"type": "dashed"},
                                "yAxisIndex": 1,
                                "itemStyle": {
                                    "color": "#ff8787"
                                },  # Lighter red for total swap
                            },
                        ],
                    }

                    ram_chart = ui.echart(ram_echart_options).classes("w-full h-64")

                    def update_ram_echart():
                        ram_chart.options["xAxis"]["data"] = ram_history["timestamps"]
                        ram_chart.options["series"][0]["data"] = ram_history["robin"]
                        ram_chart.options["series"][1]["data"] = ram_history["peak"]
                        ram_chart.options["series"][2]["data"] = ram_history[
                            "available"
                        ]
                        ram_chart.options["series"][3]["data"] = ram_history["free"]
                        ram_chart.options["series"][4]["data"] = ram_history[
                            "swap_used"
                        ]
                        ram_chart.options["series"][5]["data"] = ram_history[
                            "swap_total"
                        ]
                        ram_chart.update()

                    # UI timer just for refreshing the plot
                    ui.timer(10.0, update_ram_echart)

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
        #await ui.context.client.connected()
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

    ui.timer(0.5, show_disclaimer, once=True)

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
                            ui.menu_item(
                                "Workflow",
                                lambda: ui.navigate.to("/workflow"),
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

    with ui.column().classes("w-full h-full") as main_content:
        pass

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

    with main_content:
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


@ui.page("/workflow")
def workflow_page():
    """Display the ROBIN workflow diagram."""
    with frame(
        "ROBIN Workflow",
        smalltitle="Workflow",
    ):
        # Get versions for the diagram
        modkit_version = get_modkit_version()
        sturgeon_version = get_sturgeon_version()
        crossnn_version = get_crossnn_version()
        cnv_from_bam_version = get_cnv_from_bam_version()

        ui.mermaid(
            f"""
flowchart TD
    %% Style definitions
    classDef minKNOW fill:#f5fafd80,stroke:#339af0,stroke-width:1px,color:#1c7ed6,font-size:14px,font-weight:500
    classDef robin fill:#f6fcf7,stroke:#495057,stroke-width:1px,color:#388e3c,font-size:14px,font-weight:500
    classDef classifier fill:#fff0f699,stroke:#e64980,stroke-width:1px,color:#c2255c,font-size:14px,font-weight:500
    classDef analysis fill:#fff4e699,stroke:#fd7e14,stroke-width:1px,color:#e8590c,font-size:14px,font-weight:500
    classDef output fill:#ebfbee99,stroke:#40c057,stroke-width:1px,color:#2b8a3e,font-size:14px,font-weight:500
    classDef headerLabel fill:#ffffff00,stroke:#ffffff00,color:#222,font-size:18px,font-weight:700

    %% Custom header nodes
    robinLabel["<b>R.O.B.I.N v{__about__.__version__}</b>"]:::headerLabel
    minKNOWLabel["<b>MinKNOW Pipeline</b>"]:::headerLabel

    %% MinKNOW pipeline
    subgraph MinKNOW[" "]
        sequencing["Sequencing"]
        adaptive["Adaptive Sampling"]
        alignment["Alignment"]
        bam["BAM Files"]
    end
    minKNOWLabel -.-> MinKNOW

    %% R.O.B.I.N. pipeline
    subgraph ROBIN[" "]
        runInfo["Extract Run Information"]
        groupRuns["Group by Run"]
        mergeBam["Merge BAM Files"]
        methylation["Extract Methylation<br>(modkit v{modkit_version})"]
        individualBam["Process BAMs Individually"]
        cnv["CNV Analysis<br>(CNV from BAM v{cnv_from_bam_version})"]
        fusion["Fusion Analysis"]
        mgmt["MGMT Analysis"]
        coverage["Coverage Analysis"]
        snv["Optional SNP/V Analysis<br>(ClairS-To)"]
        crossnnCNS["CrossNN - CNS<br>({crossnn_version})"]
        crossnnPan["CrossNN - PanCancer<br>({crossnn_version})"]
        sturgeon["Sturgeon<br>({sturgeon_version})"]
        randomForest["Random Forest"]
        integrate["Integrate Results"]
        report["Report Generation"]
    end
    robinLabel -.-> ROBIN

    %% Connections
    sequencing --> adaptive
    adaptive --> alignment
    alignment --> bam
    bam --> runInfo
    runInfo --> groupRuns

    %% Split after groupRuns
    groupRuns -->|"For methylation"| mergeBam
    groupRuns -->|"For other analyses"| individualBam

    mergeBam --> methylation
    methylation --> crossnnCNS & crossnnPan & sturgeon & randomForest
    crossnnCNS --> integrate
    crossnnPan --> integrate
    sturgeon --> integrate
    randomForest --> integrate

    individualBam --> cnv & fusion & mgmt & coverage
    coverage --> snv

    integrate --> report
    cnv --> report
    fusion --> report
    mgmt --> report
    snv --> report

    %% Styling
    class sequencing,adaptive,alignment,bam minKNOW
    class runInfo,groupRuns,mergeBam,individualBam,integrate robin
    class crossnnCNS,crossnnPan,sturgeon,randomForest classifier
    class cnv,fusion,mgmt,coverage,snv analysis
    class report output
    %% Subgraph background coloring
    style MinKNOW fill:#f5fafd80,stroke:#339af0,stroke-width:2px
    style ROBIN fill:#f6fcf7,stroke:#388e3c,stroke-width:2px
""",
            config={
                "theme": "redux",
                "look": "neo",
                "flowchart": {"curve": "basis", "defaultRenderer": "elk"},
            },
        ).classes("w-full")


def get_modkit_version():
    """Get the version of modkit installed."""
    try:
        result = subprocess.run(["modkit", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        return "unknown"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return "unknown"


def get_sturgeon_version():
    """Get the version of sturgeon installed."""
    try:
        result = subprocess.run(
            ["sturgeon", "--version"], capture_output=True, text=True
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return "unknown"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        try:
            return importlib.metadata.version("sturgeon")
        except importlib.metadata.PackageNotFoundError:
            return "unknown"


def get_crossnn_version():
    """Get the version of CrossNN installed."""
    try:
        return importlib.metadata.version("nanoDX")
    except importlib.metadata.PackageNotFoundError:
        try:
            result = subprocess.run(
                ["git", "describe", "--tags"],
                cwd=os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "submodules/nanoDX",
                ),
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return "unknown"
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return "unknown"


def get_cnv_from_bam_version():
    """Get the version of CNV from BAM installed."""
    try:
        return importlib.metadata.version("cnv_from_bam")
    except importlib.metadata.PackageNotFoundError:
        try:
            result = subprocess.run(
                ["git", "describe", "--tags"],
                cwd=os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "submodules/cnv_from_bam",
                ),
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return "unknown"
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return "unknown"


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


def get_process_ram_usage():
    """Get RAM usage for Python and R processes."""
    python_ram = 0
    r_ram = 0

    for proc in psutil.process_iter(["pid", "name", "memory_info"]):
        try:
            # Get process name and memory info
            name = proc.info["name"].lower()
            mem_info = proc.info["memory_info"]

            # Skip if memory_info is None
            if mem_info is None:
                continue

            # Calculate RAM usage in GB
            ram_gb = mem_info.rss / (1024 * 1024 * 1024)  # Convert bytes to GB

            # Only count processes that are actually using memory
            if ram_gb > 0:
                if "python" in name:
                    python_ram += ram_gb
                elif "r" in name or "Rscript" in name:
                    r_ram += ram_gb

        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            continue
        except Exception as e:
            logging.debug(
                f"Error getting memory info for process {proc.info.get('name', 'unknown')}: {str(e)}"
            )
            continue

    return round(python_ram, 2), round(
        r_ram, 2
    )  # Round to 2 decimal places for cleaner display


# Start the background timer ONCE at app startup
ui.timer(10.0, collect_ram_usage)


def debug_process_tree():
    """
    Debug function to understand why we're picking up incorrect child processes.
    Shows the full process tree and explains why each process is being counted.
    """
    try:
        main_process = psutil.Process(os.getpid())
        
        print("=== Process Tree Debug ===")
        print(f"Main ROBIN process: {main_process.name()} (PID: {main_process.pid})")
        print(f"Command line: {' '.join(main_process.cmdline())}")
        print(f"Parent PID: {main_process.ppid()}")
        
        # Get all processes in the system
        all_processes = []
        for proc in psutil.process_iter(['pid', 'name', 'ppid', 'cmdline']):
            try:
                all_processes.append(proc.info)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        
        # Find all processes that could be considered "children"
        potential_children = []
        for proc_info in all_processes:
            pid = proc_info['pid']
            ppid = proc_info['ppid']
            
            # Check if this process is in our process tree
            if _is_in_process_tree(pid, main_process.pid, all_processes):
                potential_children.append(proc_info)
        
        print(f"\nFound {len(potential_children)} processes in ROBIN's process tree:")
        
        for i, proc_info in enumerate(potential_children):
            pid = proc_info['pid']
            name = proc_info['name']
            ppid = proc_info['ppid']
            cmdline = ' '.join(proc_info['cmdline']) if proc_info['cmdline'] else 'N/A'
            
            # Determine why this process is being counted
            reason = _explain_process_inclusion(proc_info, main_process.pid)
            
            print(f"\n{i+1}. {name} (PID: {pid}, PPID: {ppid})")
            print(f"   Command: {cmdline}")
            print(f"   Reason: {reason}")
            
            # Show memory usage if available
            try:
                proc = psutil.Process(pid)
                rss = proc.memory_info().rss / (1024 * 1024 * 1024)
                print(f"   Memory: {rss:.2f}GB RSS")
            except:
                print(f"   Memory: <access denied>")
        
        print("\n=== End Process Tree Debug ===")
        
    except Exception as e:
        print(f"Process tree debug failed: {e}")


def _is_in_process_tree(target_pid, root_pid, all_processes):
    """
    Check if a process is in the process tree starting from root_pid.
    """
    if target_pid == root_pid:
        return True
    
    # Find the process
    target_proc = None
    for proc_info in all_processes:
        if proc_info['pid'] == target_pid:
            target_proc = proc_info
            break
    
    if not target_proc:
        return False
    
    # Recursively check if parent is in the tree
    return _is_in_process_tree(target_proc['ppid'], root_pid, all_processes)


def _explain_process_inclusion(proc_info, main_pid):
    """
    Explain why a process is being included in our child process count.
    """
    pid = proc_info['pid']
    ppid = proc_info['ppid']
    name = proc_info['name']
    cmdline = proc_info['cmdline']
    
    reasons = []
    
    # Direct child
    if ppid == main_pid:
        reasons.append("Direct child of main ROBIN process")
    
    # Python process with ROBIN in command line
    if name in ['python', 'python3', 'python3.9', 'robin']:
        if cmdline and any('robin' in arg.lower() for arg in cmdline):
            reasons.append("Python process with ROBIN in command line")
    
    # Analysis tool
    if name in ['modkit', 'matkit', 'sturgeon', 'R', 'Rscript', 'bedtools']:
        reasons.append("Analysis tool spawned by ROBIN")
    
    # Process in tree but not direct child
    if ppid != main_pid and ppid != 1:
        reasons.append("Process in ROBIN's process tree (indirect child)")
    
    # Reparented process
    if ppid == 1:
        reasons.append("Reparented to init/systemd (orphaned process)")
    
    # Shared library or dependency
    if name in ['libc', 'libpython', 'libssl', 'libcrypto']:
        reasons.append("Shared library process")
    
    if not reasons:
        reasons.append("Unknown reason - process in tree but unclear relationship")
    
    return "; ".join(reasons)
