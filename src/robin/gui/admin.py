"""GUI admin page: user management and audit log viewer."""

from __future__ import annotations

import csv
import io
import json
from typing import TYPE_CHECKING, Any, Callable, Dict, List

from nicegui import ui

from robin.gui import theme
from robin.security import get_consent_version

if TYPE_CHECKING:
    from robin.gui_launcher import GUILauncher


def _user_table_rows(launcher: "GUILauncher") -> List[Dict[str, Any]]:
    store = launcher.security_store
    consent_version = launcher.consent_version
    consent_by_user = {
        row["user_id"]: row for row in store.list_consent_status(consent_version)
    }
    rows: List[Dict[str, Any]] = []
    for user in store.list_users():
        roles = store.get_user_roles(user.id)
        consent = consent_by_user.get(user.id, {})
        rows.append(
            {
                "username": user.username,
                "roles": ", ".join(roles) or "—",
                "active": "yes" if user.is_active else "no",
                "password": "must change" if user.must_change_password else "ok",
                "last_login": user.last_login_at or "never",
                "consent": "accepted" if consent.get("has_consent") else "pending",
                "consent_at": consent.get("agreed_at") or "—",
            }
        )
    return rows


def _audit_table_rows(launcher: "GUILauncher", filters: Dict[str, Any]) -> List[Dict[str, Any]]:
    events = launcher.security_store.query_audit_events(
        username=str(filters.get("username") or ""),
        event_type=str(filters.get("event_type") or ""),
        limit=int(filters.get("limit") or 100),
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


def create_admin_page(launcher: "GUILauncher") -> None:
    """Render the admin users/audit page inside the standard ROBIN shell."""
    consent_version = get_consent_version()
    audit_filters: Dict[str, Any] = {"username": "", "event_type": "", "limit": 100}

    with theme.frame(
        "R.O.B.I.N - Administration",
        smalltitle="Admin",
        batphone=False,
        center=launcher.center,
        setup_notifications=launcher._setup_notification_system,
    ):
        with ui.element("div").classes("w-full min-w-0").props("id=admin-page"):
            with ui.column().classes("w-full max-w-6xl mx-auto gap-3 p-2 md:p-3"):
                with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
                    ui.label("Administration").classes(
                        "classification-insight-heading text-headline-small"
                    )
                    ui.label(
                        f"Manage GUI users and review audit events. "
                        f"Active consent version: {consent_version}."
                    ).classes("classification-insight-foot")

                with ui.tabs().classes("w-full") as tabs:
                    users_tab = ui.tab("users", label="Users")
                    audit_tab = ui.tab("audit", label="Audit log")

                with ui.tab_panels(tabs, value=users_tab).classes("w-full"):
                    with ui.tab_panel(users_tab):
                        _build_users_panel(launcher, consent_version)
                    with ui.tab_panel(audit_tab):
                        _build_audit_panel(launcher, audit_filters)


def _build_users_panel(launcher: "GUILauncher", consent_version: str) -> None:
    with ui.element("div").classes("classification-insight-card w-full min-w-0"):
        with ui.column().classes("w-full min-w-0 gap-3 p-2 md:p-3"):
            with ui.row().classes("w-full items-center justify-between gap-2 flex-wrap"):
                ui.label("User accounts").classes("classification-insight-model")
                with ui.row().classes("gap-2"):
                    refresh_btn = ui.button("Refresh", icon="refresh").props(
                        "flat no-caps outline"
                    )
                    ui.button(
                        "Create user",
                        icon="person_add",
                        on_click=lambda: _open_create_user_dialog(launcher, refresh_users),
                    ).props("color=primary no-caps")

            user_columns = [
                {"name": "username", "label": "Username", "field": "username", "align": "left"},
                {"name": "roles", "label": "Roles", "field": "roles", "align": "left"},
                {"name": "active", "label": "Active", "field": "active", "align": "left"},
                {"name": "password", "label": "Password", "field": "password", "align": "left"},
                {"name": "last_login", "label": "Last login", "field": "last_login", "align": "left"},
                {"name": "consent", "label": "Consent", "field": "consent", "align": "left"},
                {"name": "consent_at", "label": "Consent at", "field": "consent_at", "align": "left"},
                {"name": "actions", "label": "Actions", "field": "actions", "align": "left"},
            ]
            _, user_table = theme.styled_table(
                columns=user_columns,
                rows=_user_table_rows(launcher),
                pagination=15,
            )

            def refresh_users() -> None:
                user_table.rows = _user_table_rows(launcher)
                user_table.update()

            refresh_btn.on_click(refresh_users)

            user_table.add_slot(
                "body-cell-actions",
                """
<q-td key="actions" :props="props">
  <q-btn dense flat no-caps label="Manage" icon="settings"
         @click="$parent.$emit('manage-user', props.row.username)" />
</q-td>
""",
            )

            def _on_manage_user(event: Any) -> None:
                args = getattr(event, "args", None)
                username = args
                if isinstance(args, (list, tuple)) and args:
                    username = args[0]
                if isinstance(args, dict):
                    username = args.get("username") or args.get("value")
                if username:
                    _open_manage_user_dialog(launcher, str(username), refresh_users)

            user_table.on("manage-user", _on_manage_user)


def _build_audit_panel(launcher: "GUILauncher", audit_filters: Dict[str, Any]) -> None:
    with ui.element("div").classes("classification-insight-card w-full min-w-0"):
        with ui.column().classes("w-full min-w-0 gap-3 p-2 md:p-3"):
            ui.label("Audit events").classes("classification-insight-model")

            with ui.row().classes("w-full gap-2 flex-wrap items-end"):
                username_filter = ui.input("Username").classes("min-w-[10rem]").props(
                    "dense outlined clearable"
                )
                event_filter = ui.input("Event type").classes("min-w-[12rem]").props(
                    "dense outlined clearable"
                )
                limit_filter = ui.number(
                    "Limit", value=100, min=1, max=5000, step=1
                ).classes("w-28").props("dense outlined")

            audit_columns = [
                {"name": "occurred_at", "label": "Time (UTC)", "field": "occurred_at", "align": "left"},
                {"name": "username", "label": "User", "field": "username", "align": "left"},
                {"name": "event_type", "label": "Event", "field": "event_type", "align": "left"},
                {"name": "target", "label": "Target", "field": "target", "align": "left"},
                {"name": "result", "label": "Result", "field": "result", "align": "left"},
                {"name": "ip", "label": "IP", "field": "ip", "align": "left"},
                {"name": "details", "label": "Details", "field": "details", "align": "left"},
            ]
            _, audit_table = theme.styled_table(
                columns=audit_columns,
                rows=_audit_table_rows(launcher, audit_filters),
                pagination=25,
                class_size="table-sm",
            )

            def _apply_filters() -> None:
                audit_filters["username"] = str(username_filter.value or "").strip()
                audit_filters["event_type"] = str(event_filter.value or "").strip()
                try:
                    audit_filters["limit"] = int(limit_filter.value or 100)
                except (TypeError, ValueError):
                    audit_filters["limit"] = 100
                audit_table.rows = _audit_table_rows(launcher, audit_filters)
                audit_table.update()

            def _export_csv() -> None:
                _apply_filters()
                events = launcher.security_store.query_audit_events(
                    username=audit_filters["username"],
                    event_type=audit_filters["event_type"],
                    limit=audit_filters["limit"],
                )
                fieldnames = [
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
                buf = io.StringIO()
                writer = csv.DictWriter(buf, fieldnames=fieldnames)
                writer.writeheader()
                for event in events:
                    row = dict(event)
                    row["details"] = json.dumps(event.get("details") or {}, ensure_ascii=True)
                    writer.writerow({k: row.get(k, "") for k in fieldnames})
                launcher._audit_log(
                    event_type="admin.audit.exported",
                    user_id=launcher._get_current_user_id(),
                    target_type="audit",
                    target_id="csv",
                    details={"rows": len(events), **audit_filters},
                )
                ui.download(buf.getvalue().encode("utf-8"), "robin_audit_export.csv")

            with ui.row().classes("w-full gap-2 flex-wrap"):
                ui.button("Apply filters", icon="filter_alt", on_click=_apply_filters).props(
                    "color=primary no-caps"
                )
                ui.button("Export CSV", icon="download", on_click=_export_csv).props(
                    "flat no-caps outline"
                )


def _open_create_user_dialog(
    launcher: "GUILauncher",
    on_created: Callable[[], None] | None = None,
) -> None:
    with ui.dialog() as dialog, ui.card().classes(
        "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md w-full"
    ):
        ui.label("Create user").classes(
            "classification-insight-heading text-headline-small q-mb-sm"
        )
        username_input = ui.input("Username").classes("w-full").props("outlined dense")
        password_input = ui.input("Password").classes("w-full").props(
            "outlined dense type=password"
        )
        confirm_input = ui.input("Confirm password").classes("w-full").props(
            "outlined dense type=password"
        )
        role_select = ui.select(["user", "admin"], value="user", label="Role").classes(
            "w-full"
        ).props("outlined dense")

        def _create() -> None:
            username = str(username_input.value or "").strip()
            password = str(password_input.value or "")
            confirm = str(confirm_input.value or "")
            role = str(role_select.value or "user")
            if not username:
                ui.notify("Username is required", type="negative")
                return
            if not password or password != confirm:
                ui.notify("Passwords must match and cannot be empty", type="negative")
                return
            try:
                user_id = launcher.auth_service.create_user(username, password, role=role)
            except Exception as exc:
                ui.notify(f"Could not create user: {exc}", type="negative")
                return
            launcher._audit_log(
                event_type="admin.user.created",
                user_id=launcher._get_current_user_id(),
                target_type="user",
                target_id=username,
                details={"role": role, "user_id": user_id, "source": "gui"},
            )
            ui.notify(
                f"Created user '{username}'. They must set a new password on first sign-in.",
                type="positive",
            )
            dialog.close()
            if on_created is not None:
                on_created()

        with ui.row().classes("w-full justify-end gap-2 mt-3"):
            ui.button("Cancel", on_click=dialog.close).props("flat no-caps outline")
            ui.button("Create", on_click=_create, icon="check").props("color=primary no-caps")
    dialog.open()


def _open_manage_user_dialog(
    launcher: "GUILauncher",
    username: str,
    on_changed: Callable[[], None],
) -> None:
    store = launcher.security_store
    user = store.get_user_by_username(username)
    if user is None:
        ui.notify(f"User '{username}' not found", type="negative")
        return
    roles = store.get_user_roles(user.id)

    with ui.dialog() as dialog, ui.card().classes(
        "robin-dialog-surface p-4 md:p-5 min-w-[18rem] max-w-md w-full"
    ):
        ui.label(f"Manage {username}").classes(
            "classification-insight-heading text-headline-small q-mb-sm"
        )
        ui.label(f"Roles: {', '.join(roles) or 'none'}").classes("classification-insight-foot")
        ui.label(
            f"Status: {'active' if user.is_active else 'inactive'}"
        ).classes("classification-insight-foot q-mb-md")

        new_password = ui.input("New password (optional)").classes("w-full").props(
            "outlined dense type=password"
        )
        confirm_password = ui.input("Confirm new password").classes("w-full").props(
            "outlined dense type=password"
        )

        def _reset_password() -> None:
            pwd = str(new_password.value or "")
            confirm = str(confirm_password.value or "")
            if not pwd:
                ui.notify("Enter a new password", type="warning")
                return
            if pwd != confirm:
                ui.notify("Passwords do not match", type="negative")
                return
            new_hash = launcher.auth_service.hash_password(pwd)
            if not store.set_user_password_hash(username, new_hash, must_change_password=True):
                ui.notify("Password update failed", type="negative")
                return
            launcher._audit_log(
                event_type="admin.user.password_reset",
                user_id=launcher._get_current_user_id(),
                target_type="user",
                target_id=username,
                details={"source": "gui", "must_change_password": True},
            )
            ui.notify(
                f"Password updated for {username}. They must choose a new password on next sign-in.",
                type="positive",
            )

        def _toggle_active() -> None:
            if user.is_active:
                if store.user_has_role(user.id, "admin") and store.count_active_admins() <= 1:
                    ui.notify("Cannot deactivate the last active admin", type="negative")
                    return
                if not store.set_user_active(username, False):
                    ui.notify("Deactivate failed", type="negative")
                    return
                launcher._audit_log(
                    event_type="admin.user.deactivated",
                    user_id=launcher._get_current_user_id(),
                    target_type="user",
                    target_id=username,
                    details={"source": "gui"},
                )
                ui.notify(f"Deactivated {username}", type="positive")
            else:
                store.set_user_active(username, True)
                launcher._audit_log(
                    event_type="admin.user.activated",
                    user_id=launcher._get_current_user_id(),
                    target_type="user",
                    target_id=username,
                    details={"source": "gui"},
                )
                ui.notify(f"Activated {username}", type="positive")
            dialog.close()
            on_changed()

        def _grant_admin() -> None:
            store.assign_role(user.id, "admin")
            launcher._audit_log(
                event_type="admin.user.role_granted",
                user_id=launcher._get_current_user_id(),
                target_type="user",
                target_id=username,
                details={"role": "admin", "source": "gui"},
            )
            ui.notify("Granted admin role", type="positive")
            dialog.close()
            on_changed()

        def _revoke_admin() -> None:
            if store.user_has_role(user.id, "admin") and store.count_active_admins() <= 1:
                ui.notify("Cannot revoke admin from the last active admin", type="negative")
                return
            if not store.revoke_role(user.id, "admin"):
                ui.notify("User does not have admin role", type="warning")
                return
            launcher._audit_log(
                event_type="admin.user.role_revoked",
                user_id=launcher._get_current_user_id(),
                target_type="user",
                target_id=username,
                details={"role": "admin", "source": "gui"},
            )
            ui.notify("Revoked admin role", type="positive")
            dialog.close()
            on_changed()

        with ui.column().classes("w-full gap-2"):
            ui.button("Update password", on_click=_reset_password, icon="lock").props(
                "flat no-caps outline"
            )
            if user.is_active:
                ui.button("Deactivate user", on_click=_toggle_active, icon="person_off").props(
                    "flat no-caps outline color=negative"
                )
            else:
                ui.button("Activate user", on_click=_toggle_active, icon="person").props(
                    "flat no-caps outline"
                )
            if "admin" in roles:
                ui.button("Revoke admin role", on_click=_revoke_admin, icon="shield").props(
                    "flat no-caps outline"
                )
            else:
                ui.button("Grant admin role", on_click=_grant_admin, icon="shield").props(
                    "flat no-caps outline"
                )

        with ui.row().classes("w-full justify-end gap-2 mt-3"):
            ui.button("Close", on_click=dialog.close).props("flat no-caps outline")
    dialog.open()
