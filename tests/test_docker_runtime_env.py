"""Non-interactive Docker startup hooks."""

from __future__ import annotations

from pathlib import Path

from robin.cli import _get_user_acknowledgment
from robin.gui_launcher import ensure_default_admin_password_set
from robin.security import AuthService, SecurityStore


def test_research_ack_env_accepts_i_agree(monkeypatch) -> None:
    monkeypatch.setattr("robin.cli.is_development_mode", False)
    monkeypatch.setenv("ROBIN_RESEARCH_ACK", "I agree")

    class _Store:
        def any_active_admin_has_consent(self, _version: str) -> bool:
            return False

    monkeypatch.setattr("robin.security.SecurityStore", lambda: _Store())
    assert _get_user_acknowledgment() is True


def test_research_ack_env_rejects_other_values(monkeypatch) -> None:
    monkeypatch.setattr("robin.cli.is_development_mode", False)
    monkeypatch.setenv("ROBIN_RESEARCH_ACK", "yes")

    class _Store:
        def any_active_admin_has_consent(self, _version: str) -> bool:
            return False

    monkeypatch.setattr("robin.security.SecurityStore", lambda: _Store())
    monkeypatch.setattr("builtins.input", lambda: "nope")
    assert _get_user_acknowledgment() is False


def test_admin_password_from_env(tmp_path: Path, monkeypatch) -> None:
    db_path = tmp_path / "security.db"
    monkeypatch.setenv("ROBIN_ADMIN_PASSWORD", "container-secret")
    store = SecurityStore(db_path)
    auth = AuthService(store)
    assert ensure_default_admin_password_set(auth) is True
    user = store.get_user_by_username("admin")
    assert user is not None
    assert auth.verify_password(user, "container-secret")
