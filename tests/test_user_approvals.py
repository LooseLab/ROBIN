from __future__ import annotations

from pathlib import Path

from robin.security import SecurityStore
from robin.security.user_approvals import (
    MINKNOW_REMOTE_CONTROL_KEY,
    REPORT_EXPORT_KEY,
    TRAINING_RECEIVED_KEY,
    approval_audit_details,
    effective_approvals,
    normalize_approvals,
    parse_approvals_json,
    user_has_approval,
)


def test_normalize_approvals_defaults_false() -> None:
    assert normalize_approvals(None) == {
        TRAINING_RECEIVED_KEY: False,
        REPORT_EXPORT_KEY: False,
        MINKNOW_REMOTE_CONTROL_KEY: False,
    }
    assert normalize_approvals({TRAINING_RECEIVED_KEY: True}) == {
        TRAINING_RECEIVED_KEY: True,
        REPORT_EXPORT_KEY: False,
        MINKNOW_REMOTE_CONTROL_KEY: False,
    }


def test_parse_approvals_json_handles_invalid() -> None:
    assert parse_approvals_json("") == normalize_approvals(None)
    assert parse_approvals_json("not-json") == normalize_approvals(None)


def test_approval_audit_details_records_changes() -> None:
    previous = {TRAINING_RECEIVED_KEY: False, REPORT_EXPORT_KEY: False}
    updated = {TRAINING_RECEIVED_KEY: True, REPORT_EXPORT_KEY: False}
    details = approval_audit_details(previous, updated, source="test")
    assert details["previous"][TRAINING_RECEIVED_KEY] is False
    assert details["updated"][TRAINING_RECEIVED_KEY] is True
    assert details["changes"] == {
        TRAINING_RECEIVED_KEY: {"from": False, "to": True},
    }
    assert approval_audit_details(previous, previous).get("changes") == {}


def test_store_user_approvals_roundtrip(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    user_id = store.create_user(
        "bob",
        "hash",
        approvals={TRAINING_RECEIVED_KEY: True},
    )
    store.assign_role(user_id, "user")

    user = store.get_user_by_username("bob")
    assert user is not None
    assert user.approvals[TRAINING_RECEIVED_KEY] is True
    assert user.approvals[REPORT_EXPORT_KEY] is False

    store.update_user_approvals("bob", {REPORT_EXPORT_KEY: True})
    updated = store.get_user_by_username("bob")
    assert updated is not None
    assert updated.approvals[REPORT_EXPORT_KEY] is True

    public_users = store.list_users()
    assert public_users[0].approvals[REPORT_EXPORT_KEY] is True


def test_admin_has_all_approvals_without_stored_flags(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    user_id = store.create_user("admin1", "hash")
    store.assign_role(user_id, "admin")

    assert user_has_approval(store, user_id, TRAINING_RECEIVED_KEY)
    assert user_has_approval(store, user_id, REPORT_EXPORT_KEY)
    assert user_has_approval(store, user_id, MINKNOW_REMOTE_CONTROL_KEY)
    effective = effective_approvals(store, user_id)
    assert effective[TRAINING_RECEIVED_KEY] is True
    assert effective[REPORT_EXPORT_KEY] is True
    assert effective[MINKNOW_REMOTE_CONTROL_KEY] is True


def test_inactive_user_has_no_approvals(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    user_id = store.create_user(
        "inactive",
        "hash",
        approvals={TRAINING_RECEIVED_KEY: True, REPORT_EXPORT_KEY: True},
    )
    store.assign_role(user_id, "user")
    store.set_user_active("inactive", False)

    assert not user_has_approval(store, user_id, TRAINING_RECEIVED_KEY)
