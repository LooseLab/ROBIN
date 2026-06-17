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
