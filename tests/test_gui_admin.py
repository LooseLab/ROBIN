from __future__ import annotations

from pathlib import Path

from robin.gui.admin import _audit_table_rows, _user_table_rows
from robin.security import AuthService, SecurityStore


def test_admin_user_and_audit_row_builders(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    auth = AuthService(store)
    user_id = auth.create_user("admin", "secret", role="admin")
    store.record_consent(user_id, "v1")

    class _Launcher:
        security_store = store
        consent_version = "v1"

    launcher = _Launcher()
    user_rows = _user_table_rows(launcher)  # type: ignore[arg-type]
    assert len(user_rows) == 1
    assert user_rows[0]["username"] == "admin"
    assert user_rows[0]["consent"] == "accepted"

    store.append_audit_event(
        event_type="auth.login.success",
        result="success",
        user_id=user_id,
        target_type="user",
        target_id="admin",
    )
    audit_rows = _audit_table_rows(
        launcher,  # type: ignore[arg-type]
        {"username": "", "event_type": "", "limit": 10},
    )
    assert len(audit_rows) == 1
    assert audit_rows[0]["event_type"] == "auth.login.success"
