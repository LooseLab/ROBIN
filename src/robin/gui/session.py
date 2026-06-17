"""NiceGUI per-browser session helpers for authentication state."""

from __future__ import annotations

from nicegui import app

_AUTH_SESSION_KEYS = (
    "authenticated",
    "user_id",
    "username",
    "roles",
    "must_change_password",
    "_auth_generation",
    "disclaimer_acknowledged",
)


def is_authenticated_session() -> bool:
    """Return True when user storage matches the current server auth generation."""
    if app is None:
        return False
    try:
        gen = app.storage.general.get("_auth_generation")
        return bool(
            app.storage.user.get("authenticated", False)
            and gen is not None
            and app.storage.user.get("_auth_generation") == gen
        )
    except Exception:
        return False


def clear_auth_session_fields() -> None:
    """Remove persisted identity fields while keeping preferences such as dark mode."""
    if app is None:
        return
    try:
        for key in _AUTH_SESSION_KEYS:
            app.storage.user.pop(key, None)
    except Exception:
        pass


def current_session_username() -> str:
    if not is_authenticated_session():
        return ""
    try:
        username = app.storage.user.get("username")
        if username:
            return str(username).strip()
    except Exception:
        pass
    return ""


def current_session_is_admin() -> bool:
    if not is_authenticated_session():
        return False
    try:
        roles = app.storage.user.get("roles") or []
        if isinstance(roles, str):
            roles = [roles]
        return "admin" in set(roles)
    except Exception:
        return False
