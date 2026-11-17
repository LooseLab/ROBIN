from __future__ import annotations

from unittest.mock import MagicMock

import robin.gui.session as session


def _mock_app(*, user: dict, general: dict) -> MagicMock:
    mock_app = MagicMock()

    def user_get(key: str, default=None):
        return user.get(key, default)

    def user_pop(key: str, default=None):
        return user.pop(key, default)

    def general_get(key: str, default=None):
        return general.get(key, default)

    mock_app.storage.user.get = user_get
    mock_app.storage.user.pop = user_pop
    mock_app.storage.general.get = general_get
    return mock_app


def test_stale_username_not_shown_without_valid_session(monkeypatch) -> None:
    monkeypatch.setattr(
        session,
        "app",
        _mock_app(
            user={
                "authenticated": False,
                "username": "admin",
                "roles": ["admin"],
            },
            general={"_auth_generation": "run-1"},
        ),
    )
    assert session.current_session_username() == ""
    assert session.current_session_is_admin() is False


def test_username_shown_only_with_matching_auth_generation(monkeypatch) -> None:
    monkeypatch.setattr(
        session,
        "app",
        _mock_app(
            user={
                "authenticated": True,
                "username": "alice",
                "roles": ["user"],
                "_auth_generation": "run-1",
            },
            general={"_auth_generation": "run-1"},
        ),
    )
    assert session.current_session_username() == "alice"
    assert session.current_session_is_admin() is False


def test_stale_auth_generation_hides_username(monkeypatch) -> None:
    monkeypatch.setattr(
        session,
        "app",
        _mock_app(
            user={
                "authenticated": True,
                "username": "admin",
                "roles": ["admin"],
                "_auth_generation": "old-run",
            },
            general={"_auth_generation": "new-run"},
        ),
    )
    assert session.is_authenticated_session() is False
    assert session.current_session_username() == ""
    assert session.current_session_is_admin() is False


def test_clear_auth_session_fields_preserves_preferences(monkeypatch) -> None:
    user = {
        "authenticated": True,
        "username": "admin",
        "roles": ["admin"],
        "dark_mode": True,
        "_auth_generation": "old-run",
    }
    monkeypatch.setattr(
        session,
        "app",
        _mock_app(user=user, general={"_auth_generation": "new-run"}),
    )
    session.clear_auth_session_fields()
    assert "username" not in user
    assert "roles" not in user
    assert user["dark_mode"] is True
