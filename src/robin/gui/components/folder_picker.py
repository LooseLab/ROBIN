"""In-dialog folder picker for the filesystem hosting NiceGUI."""

import platform
from pathlib import Path
from typing import List, Optional, Tuple

from nicegui import ui


class local_folder_picker(ui.dialog):
    """Folder picker dialog for selecting a server-side directory."""

    def __init__(
        self,
        directory: str,
        *,
        upper_limit: Optional[str] = ...,
        show_hidden_files: bool = False,
        multiple: bool = False,
    ) -> None:
        """Create a folder picker.

        ``upper_limit=None`` allows navigation across the whole filesystem.
        Otherwise navigation is restricted to the supplied directory, or the
        starting directory when ``upper_limit`` is omitted.
        """
        super().__init__()

        self.path = self._existing_directory(Path(directory).expanduser())
        if upper_limit is None:
            self.upper_limit = None
        else:
            limit = Path(directory if upper_limit == ... else upper_limit).expanduser()
            self.upper_limit = self._existing_directory(limit)
            if not self._is_within(self.path, self.upper_limit):
                self.path = self.upper_limit
        self.show_hidden_files = show_hidden_files
        self.multiple = multiple
        self.selected_paths: List[str] = []

        with self, ui.card().classes(
            "robin-dialog-surface workflow-folder-picker p-4 "
            "w-full max-w-2xl min-w-[18rem]"
        ):
            with ui.row().classes("w-full items-start gap-3 min-w-0"):
                with ui.column().classes("gap-0 flex-1 min-w-0"):
                    ui.label(
                        "Choose folders" if self.multiple else "Choose a folder"
                    ).classes("text-subtitle1 font-medium")
                    ui.label(
                        (
                            "Add folders to the selection, then confirm when ready."
                            if self.multiple
                            else "Open a directory, then select the current folder."
                        )
                    ).classes("workflow-folder-picker-help text-xs")
                ui.button(icon="close", on_click=self.close).props(
                    "flat round dense aria-label=Close"
                )

            self.add_drives_toggle()

            with ui.row().classes("w-full items-center gap-1 mt-2"):
                self.up_button = ui.button(
                    icon="arrow_upward", on_click=self.go_up
                ).props("flat round dense aria-label=Parent-folder")
                self.up_button.tooltip("Parent folder")
                home_button = ui.button(
                    icon="home", on_click=lambda: self.navigate_to(Path.home())
                ).props("flat round dense aria-label=Home-folder")
                home_button.tooltip("Home folder")
                refresh_button = ui.button(
                    icon="refresh", on_click=self.update_list
                ).props("flat round dense aria-label=Refresh")
                refresh_button.tooltip("Refresh")
                self.hidden_toggle = (
                    ui.checkbox(
                        "Hidden folders",
                        value=self.show_hidden_files,
                        on_change=self._toggle_hidden,
                    )
                    .props("dense")
                    .classes("ml-auto text-xs")
                )

            self.path_input = (
                ui.input(label="Current folder", value=str(self.path))
                .props("outlined dense spellcheck=false")
                .classes("w-full workflow-folder-picker-path-input")
            )
            self.path_input.on("keydown.enter", self._navigate_from_input)

            self.status_label = ui.label().classes(
                "workflow-folder-picker-status text-xs"
            )
            self.list_container = ui.column().classes(
                "w-full max-h-[24rem] overflow-y-auto min-w-0 "
                "workflow-folder-picker-list"
            )

            if self.multiple:
                with ui.row().classes("w-full items-center gap-2"):
                    ui.button(
                        "Add current folder",
                        on_click=lambda: self.add_selection(self.path),
                        icon="playlist_add",
                    ).props("flat no-caps").classes("shrink-0")
                    self.selection_count_label = ui.label().classes(
                        "workflow-folder-picker-help text-xs"
                    )
                self.selected_container = ui.column().classes(
                    "w-full gap-1 workflow-folder-picker-selected-list"
                )

            with ui.row().classes(
                "w-full items-center justify-end gap-2 mt-2 flex-wrap"
            ):
                if not self.multiple:
                    self.selection_label = ui.label().classes(
                        "workflow-folder-picker-selection text-xs font-mono "
                        "break-all min-w-0 flex-1"
                    )
                ui.button("Cancel", on_click=self.close).props("flat no-caps outline")
                if self.multiple:
                    self.confirm_button = ui.button(
                        "Use selected folders",
                        on_click=self.submit_selections,
                        icon="check",
                    ).props("color=primary no-caps")
                else:
                    ui.button(
                        "Select this folder",
                        on_click=self.select_current,
                        icon="check",
                    ).props("color=primary no-caps")

            self.update_list()
            if self.multiple:
                self.update_selections()

    @staticmethod
    def _existing_directory(path: Path) -> Path:
        """Return the nearest existing directory, falling back to home."""
        candidate = path
        while not candidate.exists() and candidate != candidate.parent:
            candidate = candidate.parent
        if candidate.exists() and candidate.is_file():
            candidate = candidate.parent
        if not candidate.exists() or not candidate.is_dir():
            return Path.home()
        return candidate

    @staticmethod
    def _is_within(path: Path, parent: Path) -> bool:
        try:
            path.resolve().relative_to(parent.resolve())
            return True
        except (OSError, ValueError):
            return False

    def _can_navigate_to(self, path: Path) -> Tuple[bool, str]:
        try:
            candidate = path.expanduser().resolve()
        except OSError as exc:
            return False, f"Cannot resolve folder: {exc}"
        if not candidate.exists():
            return False, "Folder does not exist."
        if not candidate.is_dir():
            return False, "The selected path is not a folder."
        if self.upper_limit is not None and not self._is_within(
            candidate, self.upper_limit
        ):
            return False, f"Choose a folder inside {self.upper_limit}."
        return True, ""

    def navigate_to(self, path: Path) -> None:
        allowed, message = self._can_navigate_to(path)
        if not allowed:
            self._show_status(message, error=True)
            return
        self.path = path.expanduser().resolve()
        self.update_list()

    def go_up(self) -> None:
        if self.path == self.path.parent:
            return
        if (
            self.upper_limit is not None
            and self.path.resolve() == self.upper_limit.resolve()
        ):
            return
        self.navigate_to(self.path.parent)

    def select_current(self) -> None:
        allowed, message = self._can_navigate_to(self.path)
        if not allowed:
            self._show_status(message, error=True)
            return
        self.submit([str(self.path)])

    def add_selection(self, path: Path) -> None:
        allowed, message = self._can_navigate_to(path)
        if not allowed:
            self._show_status(message, error=True)
            return
        resolved = str(path.expanduser().resolve())
        if resolved not in self.selected_paths:
            self.selected_paths.append(resolved)
        self.update_selections()

    def remove_selection(self, path: str) -> None:
        self.selected_paths = [
            selected for selected in self.selected_paths if selected != path
        ]
        self.update_selections()

    def submit_selections(self) -> None:
        if not self.selected_paths:
            self._show_status("Select at least one folder.", error=True)
            return
        self.submit(list(self.selected_paths))

    def update_selections(self) -> None:
        if not self.multiple:
            return
        count = len(self.selected_paths)
        self.selection_count_label.set_text(
            f"{count} folder{'s' if count != 1 else ''} selected"
        )
        self.confirm_button.set_enabled(count > 0)
        self.selected_container.clear()
        with self.selected_container:
            if not self.selected_paths:
                ui.label("No folders selected yet.").classes(
                    "workflow-folder-picker-muted text-xs italic"
                )
                return
            for path in self.selected_paths:
                with ui.row().classes(
                    "w-full items-center gap-2 px-2 py-1 "
                    "workflow-folder-picker-selected-row"
                ):
                    ui.icon("folder").classes(
                        "workflow-folder-picker-folder-icon shrink-0"
                    )
                    ui.label(path).classes("text-xs font-mono break-all flex-1 min-w-0")
                    ui.button(
                        icon="close",
                        on_click=lambda _, selected=path: self.remove_selection(
                            selected
                        ),
                    ).props("flat round dense aria-label=Remove")

    def _navigate_from_input(self, _=None) -> None:
        value = str(self.path_input.value or "").strip()
        if not value:
            self._show_status("Enter a folder path.", error=True)
            return
        self.navigate_to(Path(value))

    def _toggle_hidden(self, event) -> None:
        self.show_hidden_files = bool(event.value)
        self.update_list()

    def _show_status(self, message: str, *, error: bool = False) -> None:
        self.status_label.set_text(message)
        if error:
            self.status_label.classes(
                add="workflow-folder-picker-status--error",
                remove="workflow-folder-picker-status--muted",
            )
        else:
            self.status_label.classes(
                add="workflow-folder-picker-status--muted",
                remove="workflow-folder-picker-status--error",
            )

    def add_drives_toggle(self) -> None:
        """Add a drive selector on Windows when pywin32 is available."""
        if platform.system() != "Windows":
            return
        try:
            import win32api

            drives = win32api.GetLogicalDriveStrings().split("\000")[:-1]
            if drives:
                self.drives_toggle = ui.toggle(
                    drives, value=drives[0], on_change=self.update_drive
                )
        except ImportError:
            pass

    def update_drive(self) -> None:
        if hasattr(self, "drives_toggle"):
            self.navigate_to(Path(self.drives_toggle.value))

    def _directory_entries(self) -> Tuple[List[Path], int, Optional[str]]:
        directories: List[Path] = []
        file_count = 0
        try:
            for entry in self.path.iterdir():
                if not self.show_hidden_files and entry.name.startswith("."):
                    continue
                try:
                    if entry.is_dir():
                        directories.append(entry)
                    else:
                        file_count += 1
                except OSError:
                    continue
        except (OSError, PermissionError) as exc:
            return [], 0, f"Cannot read this folder: {exc}"
        directories.sort(key=lambda entry: entry.name.casefold())
        return directories, file_count, None

    def update_list(self) -> None:
        """Refresh the folder list and selection summary."""
        self.path_input.value = str(self.path)
        if not self.multiple:
            self.selection_label.set_text(f"Selected: {self.path}")

        at_root = self.path == self.path.parent
        at_limit = (
            self.upper_limit is not None
            and self.path.resolve() == self.upper_limit.resolve()
        )
        self.up_button.set_enabled(not at_root and not at_limit)

        directories, file_count, error = self._directory_entries()
        if error:
            self._show_status(error, error=True)
        elif directories:
            hidden_note = (
                "" if self.show_hidden_files else " Hidden folders are excluded."
            )
            self._show_status(
                f"{len(directories)} folder(s). {file_count} file(s) not shown."
                f"{hidden_note}"
            )
        elif file_count:
            self._show_status(f"No subfolders. {file_count} file(s) not shown.")
        else:
            self._show_status("This folder is empty.")

        self.list_container.clear()
        with self.list_container:
            if error:
                self._empty_state("lock", "This location cannot be opened.")
                return
            if not directories:
                self._empty_state("folder_off", "No subfolders here.")
                return

            for directory in directories:
                with ui.row().classes(
                    "w-full items-center gap-2 px-3 py-2 cursor-pointer "
                    "workflow-folder-picker-row"
                ):
                    ui.icon("folder").classes(
                        "text-lg shrink-0 workflow-folder-picker-folder-icon"
                    )
                    ui.label(directory.name).classes(
                        "font-medium flex-1 min-w-0 break-all cursor-pointer"
                    ).on(
                        "click",
                        lambda _, target=directory: self.navigate_to(target),
                    )
                    if self.multiple:
                        ui.button(
                            icon="add",
                            on_click=lambda _, target=directory: self.add_selection(
                                target
                            ),
                        ).props("flat round dense aria-label=Add-folder").tooltip(
                            "Add this folder"
                        )
                    else:
                        ui.icon("chevron_right").classes(
                            "text-base shrink-0 workflow-folder-picker-muted"
                        )

    @staticmethod
    def _empty_state(icon: str, message: str) -> None:
        with ui.column().classes("w-full items-center justify-center gap-2 py-8"):
            ui.icon(icon).classes("text-2xl workflow-folder-picker-muted")
            ui.label(message).classes("workflow-folder-picker-muted text-sm")
