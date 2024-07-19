import platform
from pathlib import Path
from typing import Optional, Any, Dict
from collections import defaultdict
import logging

from nicegui import events, ui

from robin.utilities.bam_handler import BamEventHandler

# Configure logging
logger = logging.getLogger(__name__)


def configure_logging(level: int = logging.INFO) -> None:
    """
    Configure the logging for this module.

    Args:
        level (int): The logging level to set. Defaults to logging.INFO.
    """
    logger.setLevel(level)

    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)


class LocalFilePicker(ui.dialog):
    """
    A class to create a local file picker dialog using NiceGUI.

    This class provides a user interface to select files from the local filesystem.

    Attributes:
        path (Path): The current directory path.
        upper_limit (Optional[Path]): The directory limit up to which navigation is allowed.
        show_hidden_files (bool): Whether to show hidden files in the file picker.
        grid (ui.aggrid): The grid UI component displaying the files and directories.
        drives_toggle (Optional[ui.toggle]): The toggle UI component for selecting drives (Windows only).
    """

    def __init__(
        self,
        directory: str,
        *,
        upper_limit: Optional[str] = None,
        multiple: bool = False,
        show_hidden_files: bool = False,
    ) -> None:
        """
        Initializes the LocalFilePicker with specified parameters.

        Args:
            directory (str): The directory to start in.
            upper_limit (Optional[str]): The directory to stop at (None: no limit, default: same as the starting directory).
            multiple (bool): Whether to allow multiple files to be selected.
            show_hidden_files (bool): Whether to show hidden files.
        """
        super().__init__()

        self.path = Path(directory).expanduser()
        self.upper_limit = (
            Path(directory if upper_limit == "..." else upper_limit).expanduser()
            if upper_limit
            else None
        )
        self.show_hidden_files = show_hidden_files
        self.multiple = multiple

        with self, ui.card():
            self.drives_toggle = self.add_drives_toggle()
            self.grid = self.create_grid()
            with ui.row().classes("w-full justify-end"):
                ui.button("Cancel", on_click=self.close).props("outline")
                ui.button("Ok", on_click=self._handle_ok)
        self.update_grid()

    def add_drives_toggle(self) -> Optional[ui.toggle]:
        """
        Adds a toggle for selecting drives if the platform is Windows.

        Returns:
            Optional[ui.toggle]: The created toggle component or None if not on Windows.
        """
        if platform.system() == "Windows":
            import win32api

            drives = win32api.GetLogicalDriveStrings().split("\000")[:-1]
            return ui.toggle(drives, value=drives[0], on_change=self.update_drive)
        return None

    def create_grid(self) -> ui.aggrid:
        """
        Creates and returns the grid component for displaying files and directories.

        Returns:
            ui.aggrid: The created grid component.
        """
        return (
            ui.aggrid(
                {
                    "columnDefs": [{"field": "name", "headerName": "File"}],
                    "rowSelection": "multiple" if self.multiple else "single",
                },
                html_columns=[0],
            )
            .classes("w-96")
            .on("cellDoubleClicked", self.handle_double_click)
        )

    def update_drive(self) -> None:
        """
        Updates the current path based on the selected drive and refreshes the grid.
        """
        if self.drives_toggle:
            self.path = Path(self.drives_toggle.value).expanduser()
            self.update_grid()
            logger.info(f"Drive changed to: {self.path}")

    def update_grid(self) -> None:
        """
        Updates the grid to display the contents of the current directory.
        Filters hidden files if the option is set.
        """
        try:
            paths = list(self.path.glob("*"))
            if not self.show_hidden_files:
                paths = [p for p in paths if not p.name.startswith(".")]
            paths.sort(key=lambda p: (not p.is_dir(), p.name.lower()))

            self.grid.options["rowData"] = [
                {
                    "name": f"📁 <strong>{p.name}</strong>" if p.is_dir() else p.name,
                    "path": str(p),
                }
                for p in paths
            ]

            if self.can_navigate_up():
                self.grid.options["rowData"].insert(
                    0,
                    {
                        "name": "📁 <strong>..</strong>",
                        "path": str(self.path.parent),
                    },
                )
            self.grid.update()
            logger.debug(f"Grid updated for path: {self.path}")
        except Exception as e:
            logger.error(f"Error updating grid: {e}")

    def can_navigate_up(self) -> bool:
        """
        Determines if navigation to the parent directory is allowed.

        Returns:
            bool: True if navigation up is allowed, False otherwise.
        """
        return (self.upper_limit is None and self.path != self.path.parent) or (
            self.upper_limit is not None and self.path != self.upper_limit
        )

    def handle_double_click(self, e: events.GenericEventArguments) -> None:
        """
        Handles double-click events on the grid.

        If a directory is double-clicked, navigates into the directory and updates the grid.
        If a file is double-clicked, submits the selected file.

        Args:
            e (events.GenericEventArguments): The event arguments from the double-click event.
        """
        clicked_path = Path(e.args["data"]["path"])
        if clicked_path.is_dir():
            self.path = clicked_path
            self.update_grid()
            logger.info(f"Navigated to directory: {self.path}")
        else:
            self.submit([str(clicked_path)])
            logger.info(f"File selected: {clicked_path}")

    async def _handle_ok(self) -> None:
        """
        Handles the OK button click event.

        Submits the selected file(s) from the grid.
        """
        try:
            rows = await ui.run_javascript(
                f"getElement({self.grid.id}).gridOptions.api.getSelectedRows()"
            )
            selected_paths = [r["path"] for r in rows]
            self.submit(selected_paths)
            logger.info(f"Selected files: {selected_paths}")
        except Exception as e:
            logger.error(f"Error handling OK button: {e}")


def create_bam_count() -> Dict[str, Dict[str, Any]]:
    """
    Creates and returns a defaultdict for tracking BAM file counts and details.

    Returns:
        Dict[str, Dict[str, Any]]: A dictionary with 'counter' initialized to 0 and an empty 'file' dict.
    """
    return defaultdict(lambda: {"counter": 0, "file": {}})


# Example of how to initialize and use LocalFilePicker
if __name__ == "__main__":
    configure_logging(level=logging.DEBUG)

    bam_count = create_bam_count()
    event_handler = BamEventHandler(bam_count=bam_count)

    async def open_file_picker():
        result = await LocalFilePicker("~/", multiple=True).open()
        if result:
            logger.info(f"Selected files: {result}")

    ui.button("Open File Picker", on_click=open_file_picker)
    ui.run()
