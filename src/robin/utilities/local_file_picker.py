import platform
from pathlib import Path
from typing import Optional
from collections import defaultdict

from nicegui import events, ui

from robin.utilities.bam_handler import BamEventHandler

# Configure logging
# logger = logging.getLogger(__name__)


class LocalFilePicker(ui.dialog):
    """
    A class to create a local file picker dialog using NiceGUI.

    This class provides a user interface to select files from the local filesystem.

    Attributes:
        path (Path): The current directory path.
        upper_limit (Optional[Path]): The directory limit up to which navigation is allowed.
        show_hidden_files (bool): Whether to show hidden files in the file picker.
        grid (Any): The grid UI component displaying the files and directories.
        drives_toggle (Optional[Any]): The toggle UI component for selecting drives (Windows only).
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
        if upper_limit is None:
            self.upper_limit = None
        else:
            self.upper_limit = Path(
                directory if upper_limit == "..." else upper_limit
            ).expanduser()
        self.show_hidden_files = show_hidden_files

        with self, ui.card():
            self.add_drives_toggle()
            self.grid = (
                ui.aggrid(
                    {
                        "columnDefs": [{"field": "name", "headerName": "File"}],
                        "rowSelection": "multiple" if multiple else "single",
                    },
                    html_columns=[0],
                )
                .classes("w-96")
                .on("cellDoubleClicked", self.handle_double_click)
            )
            with ui.row().classes("w-full justify-end"):
                ui.button("Cancel", on_click=self.close).props("outline")
                ui.button("Ok", on_click=self._handle_ok)
        self.update_grid()

    def add_drives_toggle(self) -> None:
        """
        Adds a toggle for selecting drives if the platform is Windows.
        """
        if platform.system() == "Windows":
            import win32api

            drives = win32api.GetLogicalDriveStrings().split("\000")[:-1]
            self.drives_toggle = ui.toggle(
                drives, value=drives[0], on_change=self.update_drive
            )

    def update_drive(self) -> None:
        """
        Updates the current path based on the selected drive and refreshes the grid.
        """
        self.path = Path(self.drives_toggle.value).expanduser()
        self.update_grid()

    def update_grid(self) -> None:
        """
        Updates the grid to display the contents of the current directory.
        Filters hidden files if the option is set.
        """
        paths = list(self.path.glob("*"))
        if not self.show_hidden_files:
            paths = [p for p in paths if not p.name.startswith(".")]
        paths.sort(key=lambda p: p.name.lower())
        paths.sort(key=lambda p: not p.is_dir())

        self.grid.options["rowData"] = [
            {
                "name": f"📁 <strong>{p.name}</strong>" if p.is_dir() else p.name,
                "path": str(p),
            }
            for p in paths
        ]
        if (
            self.upper_limit is None
            and self.path != self.path.parent
            or self.upper_limit is not None
            and self.path != self.upper_limit
        ):
            self.grid.options["rowData"].insert(
                0,
                {
                    "name": "📁 <strong>..</strong>",
                    "path": str(self.path.parent),
                },
            )
        self.grid.update()

    def handle_double_click(self, e: events.GenericEventArguments) -> None:
        """
        Handles double-click events on the grid.

        If a directory is double-clicked, navigates into the directory and updates the grid.
        If a file is double-clicked, submits the selected file.

        Args:
            e (events.GenericEventArguments): The event arguments from the double-click event.
        """
        self.path = Path(e.args["data"]["path"])
        if self.path.is_dir():
            self.update_grid()
        else:
            self.submit([str(self.path)])

    async def _handle_ok(self) -> None:
        """
        Handles the OK button click event.

        Submits the selected file(s) from the grid.
        """
        rows = await ui.run_javascript(
            f"getElement({self.grid.id}).gridOptions.api.getSelectedRows()"
        )
        self.submit([r["path"] for r in rows])


# Example of how to initialize and use LocalFilePicker
if __name__ == "__main__":
    bam_count = defaultdict(lambda: {"counter": 0, "file": {}})
    event_handler = BamEventHandler(bam_count=bam_count)
