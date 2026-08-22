"""Change-password page for first-login and voluntary password updates."""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from nicegui import app, ui

from robin.gui import theme

if TYPE_CHECKING:
    from robin.gui_launcher import GUILauncher


def create_change_password_page(
    launcher: "GUILauncher",
    *,
    redirect_to: str = "/",
    voluntary: bool = False,
) -> None:
    """Render the password change form."""
    user_id = launcher._get_current_user_id()
    username = launcher._get_current_username() or ""
    forced = bool(
        user_id is not None
        and launcher.security_store.user_must_change_password(user_id)
    )
    if voluntary and forced:
        voluntary = False

    safe_target = (
        redirect_to
        if redirect_to and redirect_to not in ("/login", "/change-password")
        else "/"
    )

    with theme.frame(
        "R.O.B.I.N - Change password",
        smalltitle="Change password",
        batphone=launcher.batman_mode,
        center=launcher.center,
        setup_notifications=launcher._setup_notification_system,
    ):
        with (
            ui.element("div").classes("w-full min-w-0").props("id=change-password-page")
        ):
            with ui.column().classes(
                "w-full max-w-md mx-auto items-center justify-center min-h-[60vh] p-4 gap-3"
            ):
                with ui.element("div").classes(
                    "w-full classification-insight-shell min-w-0"
                ):
                    ui.label("Change password").classes(
                        "classification-insight-heading text-headline-small text-center w-full"
                    )
                    if forced:
                        ui.label(
                            "Your account is using a temporary password. "
                            "Choose a new password before continuing."
                        ).classes("classification-insight-foot q-mb-md")
                    elif username:
                        ui.label(f"Signed in as {username}").classes(
                            "classification-insight-foot q-mb-md"
                        )

                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes("w-full gap-3 p-2 md:p-3"):
                            current_input = None
                            if not forced:
                                current_input = (
                                    ui.input("Current password")
                                    .classes("w-full")
                                    .props(
                                        "outlined dense type=password autocomplete=current-password"
                                    )
                                )
                            new_input = (
                                ui.input("New password")
                                .classes("w-full")
                                .props(
                                    "outlined dense type=password autocomplete=new-password"
                                )
                            )
                            confirm_input = (
                                ui.input("Confirm new password")
                                .classes("w-full")
                                .props(
                                    "outlined dense type=password autocomplete=new-password"
                                )
                            )

                            def _submit() -> None:
                                if user_id is None:
                                    ui.notify(
                                        "Session expired. Please sign in again.",
                                        type="negative",
                                    )
                                    ui.navigate.to("/login")
                                    return
                                new_password = str(new_input.value or "")
                                confirm = str(confirm_input.value or "")
                                if new_password != confirm:
                                    ui.notify(
                                        "New passwords do not match", type="negative"
                                    )
                                    return
                                current_password: Optional[str] = None
                                if current_input is not None:
                                    current_password = str(current_input.value or "")
                                try:
                                    launcher.auth_service.change_password(
                                        int(user_id),
                                        new_password,
                                        current_password=current_password,
                                    )
                                except ValueError as exc:
                                    ui.notify(str(exc), type="negative")
                                    return
                                launcher._audit_log(
                                    event_type="auth.password.changed",
                                    user_id=int(user_id),
                                    target_type="user",
                                    target_id=username,
                                    details={
                                        "forced": forced,
                                        "voluntary": voluntary and not forced,
                                    },
                                )
                                from nicegui import app

                                try:
                                    app.storage.user["must_change_password"] = False
                                except Exception:
                                    pass
                                ui.notify("Password updated", type="positive")
                                launcher._complete_post_login_navigation(
                                    int(user_id),
                                    safe_target,
                                )

                            ui.button(
                                "Update password",
                                on_click=_submit,
                                icon="lock",
                            ).props("color=primary no-caps").classes("w-full")
                            if not forced:
                                ui.button(
                                    "Cancel",
                                    on_click=lambda: ui.navigate.to(safe_target),
                                ).props("flat no-caps outline").classes("w-full")
