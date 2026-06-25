"""GUI admin page: user management and audit log viewer."""

from __future__ import annotations

import csv
import io
import json
from typing import TYPE_CHECKING, Any, Callable, Dict, List

from nicegui import ui

from robin.gui import theme
from robin.gui.display_config import (
    DISPLAY_GROUP_LABELS,
    DISPLAY_GROUP_ORDER,
    DISPLAY_ROLE_LABELS,
    DISPLAY_ROLES,
    DISPLAY_SECTIONS,
    SampleDisplayConfig,
    config_from_workflow_steps,
    effective_section_map,
)
from robin.security import get_consent_version
from robin.security.user_metadata import CLINICAL_ROLE_KEY, EMAIL_KEY, NOTES_KEY
from robin.security.user_approvals import (
    ADMIN_USER_APPROVALS_UPDATED_EVENT,
    MINKNOW_REMOTE_CONTROL_KEY,
    REPORT_EXPORT_KEY,
    TRAINING_RECEIVED_KEY,
    USER_APPROVAL_FIELDS,
    approval_audit_details,
    default_approvals,
    effective_approvals,
)

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
        effective = effective_approvals(store, user.id)
        rows.append(
            {
                "username": user.username,
                "email": user.metadata.get(EMAIL_KEY) or "—",
                "clinical_role": user.metadata.get(CLINICAL_ROLE_KEY) or "—",
                "roles": ", ".join(roles) or "—",
                "training": "yes" if effective.get(TRAINING_RECEIVED_KEY) else "no",
                "report_export": "yes" if effective.get(REPORT_EXPORT_KEY) else "no",
                "minknow_remote_control": (
                    "yes" if effective.get(MINKNOW_REMOTE_CONTROL_KEY) else "no"
                ),
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
                    display_tab = ui.tab("display", label="Sample display")
                    plotting_tab = ui.tab("plotting", label="Plotting")

                with ui.tab_panels(tabs, value=users_tab).classes("w-full"):
                    with ui.tab_panel(users_tab):
                        _build_users_panel(launcher, consent_version)
                    with ui.tab_panel(audit_tab):
                        _build_audit_panel(launcher, audit_filters)
                    with ui.tab_panel(display_tab):
                        _build_sample_display_panel(launcher)
                    with ui.tab_panel(plotting_tab):
                        _build_plotting_preferences_panel(launcher)


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
                {"name": "email", "label": "Email", "field": "email", "align": "left"},
                {
                    "name": "clinical_role",
                    "label": "Clinical role",
                    "field": "clinical_role",
                    "align": "left",
                },
                {
                    "name": "training",
                    "label": "Training",
                    "field": "training",
                    "align": "left",
                },
                {
                    "name": "report_export",
                    "label": "Report export",
                    "field": "report_export",
                    "align": "left",
                },
                {
                    "name": "minknow_remote_control",
                    "label": "MinKNOW control",
                    "field": "minknow_remote_control",
                    "align": "left",
                },
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


def _build_sample_display_panel(launcher: "GUILauncher") -> None:
    workflow_steps = (
        launcher.workflow_steps if hasattr(launcher, "workflow_steps") else None
    )
    current = (
        launcher.display_config
        if getattr(launcher, "display_config", None) is not None
        else SampleDisplayConfig()
    )
    checkbox_state: Dict[str, Dict[str, Any]] = {
        role: {} for role in DISPLAY_ROLES
    }

    def _initial_visible(section_id: str, role: str) -> bool:
        return effective_section_map(
            workflow_steps,
            current,
            viewer_role=role,
        ).get(section_id, True)

    with ui.element("div").classes("classification-insight-card w-full min-w-0"):
        with ui.column().classes("w-full min-w-0 gap-3 p-2 md:p-3"):
            ui.label("Sample page visibility").classes("classification-insight-model")
            ui.label(
                "Configure what standard users and administrators see on sample pages, "
                "the More details page, and PDF reports. Sections not enabled in the "
                "current workflow cannot be shown."
            ).classes("classification-insight-foot")

            with ui.tabs().classes("w-full") as role_tabs:
                role_tab_items = {
                    role: ui.tab(role, label=DISPLAY_ROLE_LABELS[role])
                    for role in DISPLAY_ROLES
                }

            with ui.tab_panels(role_tabs, value=role_tab_items["user"]).classes("w-full"):
                for role in DISPLAY_ROLES:
                    with ui.tab_panel(role_tab_items[role]):
                        ui.label(
                            "Unchecked sections are hidden for this role."
                        ).classes("classification-insight-foot mb-2")
                        for group in DISPLAY_GROUP_ORDER:
                            ui.label(DISPLAY_GROUP_LABELS[group]).classes(
                                "classification-insight-meta font-medium mt-2"
                            )
                            with ui.column().classes("w-full gap-1 pl-2"):
                                for section in DISPLAY_SECTIONS.values():
                                    if section.group != group:
                                        continue
                                    indent = "pl-4" if section.parent_id else ""
                                    checkbox_state[role][section.id] = ui.checkbox(
                                        section.label,
                                        value=_initial_visible(section.id, role),
                                    ).classes(indent)
                                    if section.parent_id and not _initial_visible(
                                        section.parent_id, role
                                    ):
                                        checkbox_state[role][section.id].disable()

            status_label = ui.label("").classes("classification-insight-meta")

            def _sync_parent_state(role: str) -> None:
                for section_id, section in DISPLAY_SECTIONS.items():
                    if not section.parent_id:
                        continue
                    parent_box = checkbox_state[role].get(section.parent_id)
                    child_box = checkbox_state[role].get(section_id)
                    if parent_box is None or child_box is None:
                        continue
                    if parent_box.value:
                        child_box.enable()
                    else:
                        child_box.set_value(False)
                        child_box.disable()

            for role in DISPLAY_ROLES:
                for section_id, section in DISPLAY_SECTIONS.items():
                    if section.parent_id:
                        parent_box = checkbox_state[role].get(section.parent_id)
                        if parent_box is not None:
                            parent_box.on_value_change(
                                lambda _e, r=role: _sync_parent_state(r)
                            )

            def _save() -> None:
                nonlocal current
                updated = current
                for role in DISPLAY_ROLES:
                    sections = {
                        section_id: bool(box.value)
                        for section_id, box in checkbox_state[role].items()
                    }
                    updated = updated.with_role_updates(role, sections)
                launcher.save_display_config(
                    updated,
                    user_id=launcher._get_current_user_id(),
                )
                current = launcher.display_config
                status_label.set_text(
                    f"Saved at {updated.updated_at or 'now'}"
                    + (f" by {updated.updated_by}" if updated.updated_by else "")
                )
                ui.notify("Sample display settings saved", type="positive")

            def _reset_to_workflow() -> None:
                active_role = role_tabs.value
                if active_role not in role_tab_items:
                    active_role = "user"
                mapped = config_from_workflow_steps(workflow_steps, role=active_role)
                for section_id, box in checkbox_state[active_role].items():
                    visible = mapped.role_sections.get(active_role, {}).get(
                        section_id,
                        DISPLAY_SECTIONS[section_id].default_visible,
                    )
                    box.set_value(visible)
                _sync_parent_state(active_role)
                status_label.set_text(
                    f"Reset {DISPLAY_ROLE_LABELS[active_role]} to workflow defaults "
                    "(not saved yet)"
                )

            with ui.row().classes("w-full gap-2 flex-wrap mt-2"):
                ui.button("Save", icon="save", on_click=_save).props("color=primary no-caps")
                ui.button(
                    "Reset active role to workflow defaults",
                    icon="restart_alt",
                    on_click=_reset_to_workflow,
                ).props("flat no-caps outline")

            if current.updated_at:
                by = f" by {current.updated_by}" if current.updated_by else ""
                status_label.set_text(f"Last saved: {current.updated_at}{by}")


def _build_plotting_preferences_panel(launcher: "GUILauncher") -> None:
    from robin.gui.plotting_preferences import (
        CNV_REPORT_SCALE_LABELS,
        CNV_REPORT_SCALES,
        PlottingPreferencesConfig,
    )

    current = (
        launcher.plotting_preferences
        if getattr(launcher, "plotting_preferences", None) is not None
        else PlottingPreferencesConfig()
    )

    with ui.element("div").classes("classification-insight-card w-full min-w-0"):
        with ui.column().classes("w-full min-w-0 gap-3 p-2 md:p-3"):
            ui.label("Report plotting").classes("classification-insight-model")
            ui.label(
                "Global defaults for figures in generated PDF reports. "
                "CNV plots can use absolute ploidy or log2(ploidy / expected copy number). "
                "Reports from the GUI and reporting CLI "
                "use this setting unless --cnv-normalized-difference is passed on the CLI."
            ).classes("classification-insight-foot mb-2")

            ui.label("Copy number variation (CNV)").classes(
                "classification-insight-meta font-medium mt-2"
            )
            scale_toggle = ui.toggle(
                {
                    scale: CNV_REPORT_SCALE_LABELS[scale]
                    for scale in CNV_REPORT_SCALES
                },
                value=current.cnv_report_scale,
            ).classes("w-full")

            status_label = ui.label("").classes("classification-insight-meta")

            def _save() -> None:
                nonlocal current
                scale = str(scale_toggle.value or current.cnv_report_scale)
                if scale not in CNV_REPORT_SCALES:
                    scale = current.cnv_report_scale
                updated = current.with_updates(cnv_report_scale=scale)
                launcher.save_plotting_preferences(
                    updated,
                    user_id=launcher._get_current_user_id(),
                )
                current = launcher.plotting_preferences
                status_label.set_text(
                    f"Saved at {updated.updated_at or 'now'}"
                    + (f" by {updated.updated_by}" if updated.updated_by else "")
                )
                ui.notify("Plotting preferences saved", type="positive")

            with ui.row().classes("w-full gap-2 flex-wrap mt-2"):
                ui.button("Save", icon="save", on_click=_save).props("color=primary no-caps")

            if current.updated_at:
                by = f" by {current.updated_by}" if current.updated_by else ""
                status_label.set_text(f"Last saved: {current.updated_at}{by}")


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
        email_input = ui.input("Email (optional)").classes("w-full").props(
            "outlined dense"
        )
        clinical_role_input = ui.input("Clinical role (optional)").classes("w-full").props(
            "outlined dense"
        )
        approval_boxes: Dict[str, Any] = {}
        is_admin_role = {"value": str(role_select.value or "user") == "admin"}

        ui.label("Approvals").classes("classification-insight-meta font-medium mt-2")
        for field in USER_APPROVAL_FIELDS:
            approval_boxes[field.key] = ui.checkbox(field.label, value=False).classes(
                "w-full"
            )

        def _on_role_change(e: Any) -> None:
            is_admin_role["value"] = str(getattr(e, "value", e) or "user") == "admin"
            for box in approval_boxes.values():
                if is_admin_role["value"]:
                    box.set_value(True)
                    box.disable()
                else:
                    box.enable()

        role_select.on_value_change(_on_role_change)

        def _create() -> None:
            username = str(username_input.value or "").strip()
            password = str(password_input.value or "")
            confirm = str(confirm_input.value or "")
            role = str(role_select.value or "user")
            metadata = {
                EMAIL_KEY: str(email_input.value or "").strip(),
                CLINICAL_ROLE_KEY: str(clinical_role_input.value or "").strip(),
            }
            approvals = {
                field.key: bool(approval_boxes[field.key].value)
                for field in USER_APPROVAL_FIELDS
            }
            if not username:
                ui.notify("Username is required", type="negative")
                return
            if not password or password != confirm:
                ui.notify("Passwords must match and cannot be empty", type="negative")
                return
            try:
                user_id = launcher.auth_service.create_user(
                    username,
                    password,
                    role=role,
                    metadata=metadata,
                    approvals=approvals if role != "admin" else None,
                )
            except ValueError as exc:
                ui.notify(str(exc), type="negative")
                return
            except Exception as exc:
                ui.notify(f"Could not create user: {exc}", type="negative")
                return
            create_details: Dict[str, Any] = {
                "role": role,
                "user_id": user_id,
                "source": "gui",
                "metadata": metadata,
            }
            if role != "admin":
                create_details["initial_approvals"] = approval_audit_details(
                    default_approvals(),
                    approvals,
                    source="gui",
                )
            launcher._audit_log(
                event_type="admin.user.created",
                user_id=launcher._get_current_user_id(),
                target_type="user",
                target_id=username,
                details=create_details,
            )
            if role != "admin" and approval_audit_details(
                default_approvals(), approvals
            ).get("changes"):
                launcher._audit_log(
                    event_type=ADMIN_USER_APPROVALS_UPDATED_EVENT,
                    user_id=launcher._get_current_user_id(),
                    target_type="user",
                    target_id=username,
                    details=approval_audit_details(
                        default_approvals(),
                        approvals,
                        source="gui",
                        extra={"user_id": user_id, "context": "user_created"},
                    ),
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
    is_admin = store.user_has_role(user.id, "admin")

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

        email_input = ui.input("Email").classes("w-full").props("outlined dense")
        email_input.value = user.metadata.get(EMAIL_KEY, "")
        clinical_role_input = ui.input("Clinical role").classes("w-full").props(
            "outlined dense"
        )
        clinical_role_input.value = user.metadata.get(CLINICAL_ROLE_KEY, "")
        notes_input = ui.textarea("Notes").classes("w-full").props("outlined dense autogrow")
        notes_input.value = user.metadata.get(NOTES_KEY, "")

        ui.label("Approvals").classes("classification-insight-meta font-medium mt-2")
        if is_admin:
            ui.label(
                "Administrators always have all approvals granted."
            ).classes("classification-insight-foot q-mb-sm")
        approval_boxes: Dict[str, Any] = {}
        for field in USER_APPROVAL_FIELDS:
            approval_boxes[field.key] = ui.checkbox(
                field.label,
                value=bool(user.approvals.get(field.key, False)),
            ).classes("w-full")
            if is_admin:
                approval_boxes[field.key].set_value(True)
                approval_boxes[field.key].disable()

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

        def _save_profile() -> None:
            previous_approvals = dict(user.approvals)
            metadata = {
                EMAIL_KEY: str(email_input.value or "").strip(),
                CLINICAL_ROLE_KEY: str(clinical_role_input.value or "").strip(),
                NOTES_KEY: str(notes_input.value or "").strip(),
            }
            approvals = {
                field.key: bool(approval_boxes[field.key].value)
                for field in USER_APPROVAL_FIELDS
            }
            try:
                if not store.update_user_metadata(username, metadata):
                    ui.notify("Profile update failed", type="negative")
                    return
                if not is_admin and not store.update_user_approvals(username, approvals):
                    ui.notify("Approvals update failed", type="negative")
                    return
            except ValueError as exc:
                ui.notify(str(exc), type="negative")
                return
            launcher._audit_log(
                event_type="admin.user.profile_updated",
                user_id=launcher._get_current_user_id(),
                target_type="user",
                target_id=username,
                details={"metadata": metadata, "source": "gui"},
            )
            if not is_admin:
                approval_details = approval_audit_details(
                    previous_approvals,
                    approvals,
                    source="gui",
                    extra={"user_id": user.id},
                )
                if approval_details.get("changes"):
                    launcher._audit_log(
                        event_type=ADMIN_USER_APPROVALS_UPDATED_EVENT,
                        user_id=launcher._get_current_user_id(),
                        target_type="user",
                        target_id=username,
                        details=approval_details,
                    )
            ui.notify(f"Updated profile for {username}", type="positive")
            on_changed()

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
            ui.button("Save profile", on_click=_save_profile, icon="badge").props(
                "color=primary no-caps"
            )
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
