"""Admin panel: delete or archive completed sample output folders."""

from __future__ import annotations

import asyncio
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from nicegui import ui

from robin.gui import theme
from robin.gui.components.folder_picker import local_folder_picker
from robin.sample_lifecycle import (
    SampleRemovalAssessment,
    archive_sample_data,
    assess_sample_for_removal,
    delete_sample_data,
    result_to_audit_details,
    validate_archive_destination,
)

if TYPE_CHECKING:
    from robin.gui_launcher import GUILauncher

logger = logging.getLogger(__name__)


def _lifecycle_rows(launcher: "GUILauncher") -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    work_dir = (launcher.monitored_directory or "").strip()
    if not work_dir:
        return rows

    with launcher._samples_record_lock:
        records = sorted(
            launcher._samples_master_record.values(),
            key=lambda rec: rec._last_seen_raw,
            reverse=True,
        )

    for record in records:
        assessment = assess_sample_for_removal(
            origin=record.origin,
            active_jobs=record.active_jobs,
            pending_jobs=record.pending_jobs,
        )
        rows.append(
            {
                "sample_id": record.sample_id,
                "test_id": record.test_id or "—",
                "origin": record.origin,
                "last_seen": record.last_seen or "—",
                "active_jobs": record.active_jobs,
                "pending_jobs": record.pending_jobs,
                "status": "Ready" if assessment.removable else assessment.reason,
                "removable": assessment.removable,
            }
        )
    return rows


def _selected_sample_id(table: Any) -> Optional[str]:
    selected = getattr(table, "selected", None) or []
    if not selected:
        return None
    row = selected[0]
    if isinstance(row, dict):
        return str(row.get("sample_id") or "")
    return str(row)


def _require_work_dir(launcher: "GUILauncher") -> Optional[Path]:
    work_dir = (launcher.monitored_directory or "").strip()
    if not work_dir:
        ui.notify("Work directory is not configured.", type="negative")
        return None
    path = Path(work_dir)
    if not path.is_dir():
        ui.notify("Work directory does not exist.", type="negative")
        return None
    return path


def _assessment_for_sample(
    launcher: "GUILauncher", sample_id: str
) -> Optional[SampleRemovalAssessment]:
    with launcher._samples_record_lock:
        record = launcher._samples_master_record.get(sample_id)
    if record is None:
        ui.notify(f"Sample {sample_id!r} is not tracked.", type="warning")
        return None
    assessment = assess_sample_for_removal(
        origin=record.origin,
        active_jobs=record.active_jobs,
        pending_jobs=record.pending_jobs,
    )
    if not assessment.removable:
        ui.notify(assessment.reason, type="warning")
        return None
    return assessment


def build_sample_lifecycle_panel(launcher: "GUILauncher") -> None:
    """Render the sample delete/archive admin tab."""
    state: Dict[str, Any] = {"archive_destination": ""}

    with ui.element("div").classes("classification-insight-card w-full min-w-0"):
        with ui.column().classes("w-full min-w-0 gap-3 p-2 md:p-3"):
            ui.label("Sample management").classes("classification-insight-model")
            ui.label(
                "Delete or archive completed sample output folders. "
                "Live runs and samples with active or pending jobs cannot be removed. "
                "Archives are written as tar.gz files to a folder outside the work directory; "
                "the sample folder is removed after a successful archive."
            ).classes("classification-insight-foot")

            work_dir = (launcher.monitored_directory or "").strip()
            if work_dir:
                ui.label(f"Work directory: {work_dir}").classes(
                    "text-xs font-mono break-all text-slate-600 dark:text-slate-400"
                )
            else:
                ui.label("Work directory is not configured.").classes(
                    "text-sm text-amber-700 dark:text-amber-300"
                )

            with ui.row().classes("w-full items-end gap-2 flex-wrap"):
                archive_path_input = (
                    ui.input(
                        label="Archive destination folder",
                        placeholder="Choose a folder outside the work directory",
                    )
                    .props("outlined dense spellcheck=false")
                    .classes("flex-1 min-w-[14rem]")
                )

                async def pick_archive_folder() -> None:
                    start = (
                        state.get("archive_destination")
                        or work_dir
                        or str(Path.home())
                    )
                    picker = local_folder_picker(start, upper_limit=None)
                    result = await picker
                    if result:
                        selected = result[0]
                        state["archive_destination"] = selected
                        archive_path_input.value = selected

                ui.button(
                    "Browse folders",
                    on_click=pick_archive_folder,
                    icon="folder_open",
                ).props("flat no-caps outline")

            columns = [
                {
                    "name": "sample_id",
                    "label": "Library ID",
                    "field": "sample_id",
                    "align": "left",
                    "sortable": True,
                },
                {
                    "name": "test_id",
                    "label": "Test ID",
                    "field": "test_id",
                    "align": "left",
                    "sortable": True,
                },
                {
                    "name": "origin",
                    "label": "Origin",
                    "field": "origin",
                    "align": "left",
                    "sortable": True,
                },
                {
                    "name": "last_seen",
                    "label": "Last activity",
                    "field": "last_seen",
                    "align": "left",
                    "sortable": True,
                },
                {
                    "name": "status",
                    "label": "Status",
                    "field": "status",
                    "align": "left",
                },
            ]
            _, sample_table = theme.styled_table(
                columns=columns,
                rows=_lifecycle_rows(launcher),
                pagination=15,
                row_key="sample_id",
                selection="single",
            )

            def refresh_table() -> None:
                sample_table.rows = _lifecycle_rows(launcher)
                sample_table.update()

            with ui.row().classes("w-full items-center justify-between gap-2 flex-wrap"):
                ui.button("Refresh", icon="refresh", on_click=refresh_table).props(
                    "flat no-caps outline"
                )
                with ui.row().classes("gap-2"):
                    archive_button = ui.button(
                        "Archive selected",
                        icon="archive",
                        on_click=lambda: None,
                    ).props("color=secondary no-caps")
                    delete_button = ui.button(
                        "Delete selected",
                        icon="delete",
                        on_click=lambda: None,
                    ).props("color=negative no-caps")

            with ui.dialog() as confirm_dialog, ui.card().classes(
                "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md"
            ):
                confirm_title = ui.label("").classes(
                    "text-headline-small text-slate-900 dark:text-slate-50"
                )
                confirm_body = ui.label("").classes(
                    "text-body-medium text-slate-600 dark:text-slate-400"
                )
                confirm_typed = ui.input(
                    label="Type the library ID to confirm",
                    placeholder="",
                ).props("outlined dense")
                confirm_typed.set_visibility(False)

                pending_action: Dict[str, Any] = {"kind": "", "sample_id": ""}

                async def run_confirmed_action() -> None:
                    sample_id = pending_action.get("sample_id") or ""
                    kind = pending_action.get("kind") or ""
                    if not sample_id or not kind:
                        return

                    if confirm_typed.visible:
                        typed = (confirm_typed.value or "").strip()
                        if typed != sample_id:
                            ui.notify(
                                "Confirmation text does not match the library ID.",
                                type="warning",
                            )
                            return

                    work_path = _require_work_dir(launcher)
                    if work_path is None:
                        return

                    if _assessment_for_sample(launcher, sample_id) is None:
                        refresh_table()
                        confirm_dialog.close()
                        return

                    archive_dest: Optional[Path] = None
                    if kind == "archive":
                        dest_raw = (
                            state.get("archive_destination")
                            or (archive_path_input.value or "").strip()
                        )
                        if not dest_raw:
                            ui.notify(
                                "Choose an archive destination folder first.",
                                type="warning",
                            )
                            return
                        try:
                            archive_dest = validate_archive_destination(
                                Path(dest_raw), work_path
                            )
                        except ValueError as exc:
                            ui.notify(str(exc), type="warning")
                            return

                    confirm_dialog.close()

                    with ui.dialog().props("persistent") as progress_dialog:
                        with ui.card().classes("robin-dialog-surface p-4 md:p-5"):
                            with ui.row().classes("items-center gap-3"):
                                ui.spinner(size="lg")
                                progress_label = ui.label("Working…").classes(
                                    "classification-insight-model"
                                )
                    progress_dialog.open()

                    try:
                        if kind == "delete":
                            progress_label.set_text(f"Deleting {sample_id}…")
                            result = await asyncio.to_thread(
                                delete_sample_data, work_path, sample_id
                            )
                            event_type = "sample.deleted"
                        else:
                            assert archive_dest is not None
                            progress_label.set_text(f"Archiving {sample_id}…")
                            result = await asyncio.to_thread(
                                archive_sample_data,
                                work_path,
                                sample_id,
                                archive_dest,
                            )
                            event_type = "sample.archived"

                        launcher._evict_sample_from_gui_state(sample_id)
                        launcher._audit_log(
                            event_type=event_type,
                            user_id=launcher._get_current_user_id(),
                            target_type="sample",
                            target_id=sample_id,
                            details=result_to_audit_details(result),
                        )
                        if kind == "archive":
                            ui.notify(
                                f"Archived {sample_id} to {result.archive_path}",
                                type="positive",
                            )
                        else:
                            ui.notify(f"Deleted sample {sample_id}.", type="positive")
                        refresh_table()
                    except Exception as exc:
                        logger.exception(
                            "Sample %s failed (%s)", sample_id, kind, exc_info=exc
                        )
                        launcher._audit_log(
                            event_type=f"sample.{kind}.failed",
                            result="failure",
                            user_id=launcher._get_current_user_id(),
                            target_type="sample",
                            target_id=sample_id,
                            details={"error": str(exc)},
                            error_code=type(exc).__name__,
                        )
                        ui.notify(str(exc), type="negative")
                    finally:
                        progress_dialog.close()
                        confirm_typed.value = ""
                        confirm_typed.set_visibility(False)

                with ui.row().classes("w-full justify-end gap-2 mt-2"):
                    ui.button("Cancel", on_click=confirm_dialog.close).props(
                        "flat no-caps outline"
                    )
                    ui.button("Confirm", on_click=run_confirmed_action).props(
                        "color=negative no-caps"
                    )

            def _open_confirm(kind: str) -> None:
                sample_id = _selected_sample_id(sample_table)
                if not sample_id:
                    ui.notify("Select a sample first.", type="warning")
                    return
                if _assessment_for_sample(launcher, sample_id) is None:
                    refresh_table()
                    return

                pending_action["kind"] = kind
                pending_action["sample_id"] = sample_id

                if kind == "delete":
                    confirm_title.set_text("Delete sample permanently?")
                    confirm_body.set_text(
                        f"This will permanently delete all output data for "
                        f"{sample_id} under the work directory. This cannot be undone."
                    )
                    confirm_typed.set_visibility(True)
                    confirm_typed.value = ""
                else:
                    dest_raw = (
                        state.get("archive_destination")
                        or (archive_path_input.value or "").strip()
                    )
                    if not dest_raw:
                        ui.notify(
                            "Choose an archive destination folder first.",
                            type="warning",
                        )
                        return
                    work_path = _require_work_dir(launcher)
                    if work_path is None:
                        return
                    try:
                        dest = validate_archive_destination(
                            Path(dest_raw), work_path
                        )
                    except ValueError as exc:
                        ui.notify(str(exc), type="warning")
                        return
                    confirm_title.set_text("Archive sample?")
                    confirm_body.set_text(
                        f"Archive {sample_id} to {dest} as a tar.gz file, then remove "
                        f"the sample folder from the work directory."
                    )
                    confirm_typed.set_visibility(False)

                confirm_dialog.open()

            archive_button.on_click(lambda: _open_confirm("archive"))
            delete_button.on_click(lambda: _open_confirm("delete"))
