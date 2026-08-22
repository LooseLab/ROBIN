"""Per-sample audit history dialog for the GUI."""

from __future__ import annotations

import csv
import io
import json
from typing import TYPE_CHECKING, Any, Dict, List

from nicegui import ui

from robin.gui import theme

if TYPE_CHECKING:
    from robin.gui_launcher import GUILauncher

_AUDIT_COLUMNS = [
    {
        "name": "occurred_at",
        "label": "Time (UTC)",
        "field": "occurred_at",
        "align": "left",
    },
    {"name": "username", "label": "User", "field": "username", "align": "left"},
    {"name": "event_type", "label": "Event", "field": "event_type", "align": "left"},
    {"name": "target", "label": "Target", "field": "target", "align": "left"},
    {"name": "result", "label": "Result", "field": "result", "align": "left"},
    {"name": "ip", "label": "IP", "field": "ip", "align": "left"},
    {"name": "details", "label": "Details", "field": "details", "align": "left"},
]

_EXPORT_FIELDS = [
    "id",
    "occurred_at",
    "user_id",
    "username",
    "event_type",
    "target_type",
    "target_id",
    "result",
    "error_code",
    "ip",
    "user_agent",
    "session_id",
    "request_id",
    "details",
]


def _audit_rows_for_sample(
    launcher: "GUILauncher", sample_id: str, *, limit: int = 500
) -> List[Dict[str, Any]]:
    events = launcher.security_store.query_audit_events(
        sample_id=sample_id,
        limit=limit,
    )
    rows: List[Dict[str, Any]] = []
    for event in events:
        rows.append(
            {
                "occurred_at": event.get("occurred_at", ""),
                "username": event.get("username") or "—",
                "event_type": event.get("event_type", ""),
                "target": f"{event.get('target_type', '')}:{event.get('target_id', '')}".strip(
                    ":"
                ),
                "result": event.get("result", ""),
                "ip": event.get("ip") or "—",
                "details": json.dumps(event.get("details") or {}, ensure_ascii=True),
            }
        )
    return rows


def _export_sample_audit_csv(
    launcher: "GUILauncher", sample_id: str, *, limit: int = 5000
) -> bytes:
    events = launcher.security_store.query_audit_events(
        sample_id=sample_id,
        limit=limit,
    )
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=_EXPORT_FIELDS)
    writer.writeheader()
    for event in events:
        row = dict(event)
        row["details"] = json.dumps(event.get("details") or {}, ensure_ascii=True)
        writer.writerow({k: row.get(k, "") for k in _EXPORT_FIELDS})
    launcher._audit_log(
        event_type="sample.audit.exported",
        user_id=launcher._get_current_user_id(),
        target_type="sample",
        target_id=sample_id,
        details={"rows": len(events), "limit": limit},
    )
    return buf.getvalue().encode("utf-8")


def open_sample_audit_dialog(launcher: "GUILauncher", sample_id: str) -> None:
    """Show audit history for a single sample with CSV export."""
    limit = 500

    with (
        ui.dialog() as dialog,
        ui.card().classes(
            "robin-dialog-surface p-4 md:p-5 w-full max-w-5xl min-w-[18rem]"
        ),
    ):
        ui.label(f"Audit history — {sample_id}").classes(
            "classification-insight-heading text-headline-small q-mb-sm"
        )
        ui.label(
            "Actions recorded for this sample (views, reports, analysis runs, exports)."
        ).classes("classification-insight-foot q-mb-md")

        _, audit_table = theme.styled_table(
            columns=_AUDIT_COLUMNS,
            rows=_audit_rows_for_sample(launcher, sample_id, limit=limit),
            pagination=15,
            class_size="table-sm",
        )

        def _refresh() -> None:
            audit_table.rows = _audit_rows_for_sample(launcher, sample_id, limit=limit)
            audit_table.update()

        def _export() -> None:
            if not launcher._require_export_or_notify():
                return
            payload = _export_sample_audit_csv(launcher, sample_id)
            safe_name = "".join(
                ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in sample_id
            )
            ui.download(payload, f"{safe_name}_audit.csv")

        with ui.row().classes("w-full justify-end gap-2 flex-wrap mt-3"):
            ui.button("Refresh", icon="refresh", on_click=_refresh).props(
                "flat no-caps outline"
            )
            if launcher._current_user_can_export():
                ui.button("Export CSV", icon="download", on_click=_export).props(
                    "flat no-caps outline"
                )
            ui.button("Close", on_click=dialog.close).props("color=primary no-caps")

    launcher._audit_log(
        event_type="sample.audit.viewed",
        user_id=launcher._get_current_user_id(),
        target_type="sample",
        target_id=sample_id,
    )
    dialog.open()
