from __future__ import annotations

import json
from pathlib import Path

from robin.security.store import SecurityStore


def _store(tmp_path: Path) -> SecurityStore:
    return SecurityStore(db_path=tmp_path / "security.db")


def test_store_user_role_and_consent_roundtrip(tmp_path: Path) -> None:
    store = _store(tmp_path)
    user_id = store.create_user("alice", "hash123")
    store.assign_role(user_id, "admin")

    user = store.get_user_by_username("alice")
    assert user is not None
    assert user.username == "alice"
    assert store.user_has_role(user_id, "admin")
    assert store.count_active_admins() == 1

    assert not store.has_consent(user_id, "v1")
    store.record_consent(user_id, "v1", ip="127.0.0.1")
    assert store.has_consent(user_id, "v1")


def test_store_audit_query_filters(tmp_path: Path) -> None:
    store = _store(tmp_path)
    user_id = store.create_user("bob", "hash123")
    store.assign_role(user_id, "user")

    store.append_audit_event(
        event_type="auth.login.success",
        result="success",
        user_id=user_id,
        target_type="user",
        target_id="bob",
        details={"source": "test"},
    )
    store.append_audit_event(
        event_type="report.exported",
        result="failure",
        user_id=user_id,
        target_type="sample",
        target_id="S1",
        details={"filename": "x.pdf"},
        error_code="file_missing",
    )

    all_events = store.query_audit_events(limit=10)
    assert len(all_events) == 2

    login_events = store.query_audit_events(event_type="auth.login.success", limit=10)
    assert len(login_events) == 1
    assert login_events[0]["event_type"] == "auth.login.success"

    bob_events = store.query_audit_events(username="bob", limit=10)
    assert len(bob_events) == 2


def test_query_audit_events_by_sample_id(tmp_path: Path) -> None:
    store = _store(tmp_path)
    user_id = store.create_user("carol", "hash123")

    store.append_audit_event(
        event_type="sample.viewed",
        result="success",
        user_id=user_id,
        target_type="sample",
        target_id="Sample_103",
        details={},
    )
    store.append_audit_event(
        event_type="report.generated",
        result="success",
        user_id=user_id,
        target_type="sample",
        target_id="Sample_104",
        details={},
    )
    store.append_audit_event(
        event_type="report.generated",
        result="success",
        user_id=user_id,
        target_type="batch",
        target_id="3_samples",
        details={"sample_ids": ["Sample_103", "Sample_105"]},
    )

    rows = store.query_audit_events(sample_id="Sample_103", limit=50)
    assert len(rows) == 2
    event_types = {row["event_type"] for row in rows}
    assert "sample.viewed" in event_types
    assert "report.generated" in event_types
    assert all(
        row["target_id"] == "Sample_103"
        or "Sample_103" in json.dumps(row.get("details") or {})
        for row in rows
    )


def test_store_list_consent_status_pending_and_accepted(tmp_path: Path) -> None:
    store = _store(tmp_path)
    accepted_id = store.create_user("accepted", "hash")
    store.create_user("pending", "hash")
    store.record_consent(accepted_id, "v1")

    rows = {row["username"]: row for row in store.list_consent_status("v1")}
    assert rows["accepted"]["has_consent"] is True
    assert rows["pending"]["has_consent"] is False


def test_any_active_admin_has_consent(tmp_path: Path) -> None:
    store = _store(tmp_path)
    admin_id = store.create_user("admin", "hash")
    store.create_user("user", "hash")
    store.assign_role(admin_id, "admin")

    assert not store.any_active_admin_has_consent("v1")
    store.record_consent(admin_id, "v1")
    assert store.any_active_admin_has_consent("v1")

    store.set_user_active("admin", False)
    assert not store.any_active_admin_has_consent("v1")


def test_must_change_password_flag_and_admin_reset(tmp_path: Path) -> None:
    store = _store(tmp_path)
    user_id = store.create_user("dave", "hash1")
    user = store.get_user_by_id(user_id)
    assert user is not None
    assert user.must_change_password is True

    store.set_user_password_hash("dave", "hash2", must_change_password=False)
    assert not store.user_must_change_password(user_id)

    store.set_user_password_hash("dave", "hash3", must_change_password=True)
    assert store.user_must_change_password(user_id)
