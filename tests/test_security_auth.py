from __future__ import annotations

from pathlib import Path

from argon2 import PasswordHasher

from robin.security.auth import AuthService
from robin.security.store import SecurityStore


def _store_and_auth(tmp_path: Path) -> tuple[SecurityStore, AuthService]:
    store = SecurityStore(db_path=tmp_path / "security.db")
    return store, AuthService(store)


def test_auth_create_and_verify_login(tmp_path: Path) -> None:
    store, auth = _store_and_auth(tmp_path)
    auth.create_user("carol", "s3cret", role="admin")

    ok_user = auth.verify_login("carol", "s3cret")
    assert ok_user is not None
    assert ok_user.username == "carol"
    assert store.user_has_role(ok_user.id, "admin")

    bad_user = auth.verify_login("carol", "wrong")
    assert bad_user is None


def test_create_user_defaults_must_change_password(tmp_path: Path) -> None:
    store, auth = _store_and_auth(tmp_path)
    user_id = auth.create_user("alice", "temp-pass", role="user")
    user = store.get_user_by_id(user_id)
    assert user is not None
    assert user.must_change_password is True


def test_change_password_forced_skips_current(tmp_path: Path) -> None:
    store, auth = _store_and_auth(tmp_path)
    user_id = auth.create_user("bob", "shared-temp", role="user")
    assert store.user_must_change_password(user_id)

    auth.change_password(user_id, "my-new-secret")
    user = store.get_user_by_id(user_id)
    assert user is not None
    assert not user.must_change_password
    assert auth.verify_login("bob", "my-new-secret") is not None
    assert auth.verify_login("bob", "shared-temp") is None


def test_change_password_voluntary_requires_current(tmp_path: Path) -> None:
    store, auth = _store_and_auth(tmp_path)
    user_id = auth.create_user(
        "carol", "initial", role="user", must_change_password=False
    )
    auth.change_password(user_id, "updated", current_password="initial")
    assert auth.verify_login("carol", "updated") is not None

    try:
        auth.change_password(user_id, "nope", current_password="wrong")
        assert False, "expected ValueError"
    except ValueError as exc:
        assert "incorrect" in str(exc).lower()

    try:
        auth.change_password(user_id, "nope")
        assert False, "expected ValueError"
    except ValueError as exc:
        assert "required" in str(exc).lower()


def test_auth_bootstrap_from_legacy_hash(tmp_path: Path) -> None:
    store, auth = _store_and_auth(tmp_path)
    legacy_path = tmp_path / "gui_password_hash"
    legacy_hash = PasswordHasher().hash("legacy-pass")
    legacy_path.write_text(legacy_hash, encoding="utf-8")

    created = auth.bootstrap_admin_from_legacy_hash(legacy_path)
    assert created
    user = auth.verify_login("admin", "legacy-pass")
    assert user is not None
    assert store.user_has_role(user.id, "admin")
    assert not user.must_change_password


def test_bootstrap_admin_explicit_no_must_change(tmp_path: Path) -> None:
    store, auth = _store_and_auth(tmp_path)
    user_id = auth.create_user(
        "admin", "setup-pass", role="admin", must_change_password=False
    )
    user = store.get_user_by_id(user_id)
    assert user is not None
    assert not user.must_change_password


def test_bootstrap_default_admin_from_password(tmp_path: Path) -> None:
    store, auth = _store_and_auth(tmp_path)

    created = auth.bootstrap_default_admin("setup-pass")

    assert created
    user = auth.verify_login("admin", "setup-pass")
    assert user is not None
    assert store.user_has_role(user.id, "admin")
    assert not user.must_change_password
