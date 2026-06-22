from __future__ import annotations

from typing import Any, Optional, Dict
from pathlib import Path

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

# We intentionally import the heavy launcher class from the legacy module
# to keep code changes minimal while we migrate features into this package.
from ..gui_launcher import GUILauncher  # type: ignore

_current_gui_launcher: Optional[GUILauncher] = None


def launch_gui(
    host: str = "0.0.0.0",
    port: int = 8081,
    show: bool = False,
    workflow_runner: Any = None,
    workflow_steps: Optional[list] = None,
    monitored_directory: str = "",
    reload: bool = False,
    center: str = None,
    workflow_toml: Optional[str] = None,
) -> GUILauncher:
    """Central entrypoint for launching the robin GUI (refactored).

    This function mirrors the legacy signature and behavior while living under
    the package `robin.gui`. It returns the launcher instance and sets a
    global reference for update dispatching.
    """
    if ui is None:
        raise ImportError(
            "NiceGUI is not available. Please install it with: pip install nicegui"
        )

    launcher = GUILauncher(host, port, reload)

    # Normalize monitored directory to an absolute path to avoid CWD issues
    abs_mon_dir = (
        str(Path(monitored_directory).resolve()) if monitored_directory else ""
    )
    success = launcher.launch_gui(
        workflow_runner,
        workflow_steps,
        abs_mon_dir,
        reload,
        center,
        workflow_toml=workflow_toml,
    )
    if not success:
        raise RuntimeError("Failed to launch GUI")

    global _current_gui_launcher
    _current_gui_launcher = launcher

    if show:
        import webbrowser  # local to avoid unnecessary import at module import time

        webbrowser.open(launcher.get_gui_url())

    return launcher


def get_gui_launcher() -> Optional[GUILauncher]:
    return _current_gui_launcher


def send_gui_update(update_type: Any, data: Dict[str, Any], priority: int = 0) -> None:
    launcher = get_gui_launcher()
    if launcher:
        launcher.send_update(update_type, data, priority)
