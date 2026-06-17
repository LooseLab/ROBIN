from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from robin.cli import main
from robin.security import AuthService, SecurityStore, get_consent_version
from robin.security.constants import CONSENT_VERSION_ENV


@pytest.fixture()
def security_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    db_path = tmp_path / "security.db"
    monkeypatch.setattr(
        "robin.security.store.get_security_db_path",
        lambda: db_path,
    )
    monkeypatch.setattr(
        "robin.security.constants.get_security_db_path",
        lambda: db_path,
    )
    return db_path


def test_get_consent_version_env_override(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(CONSENT_VERSION_ENV, "v2")
    assert get_consent_version() == "v2"


def test_list_consent_status(security_db: Path) -> None:
    store = SecurityStore(db_path=security_db)
    user_id = store.create_user("dana", "hash")
    store.record_consent(user_id, "v1", ip="10.0.0.1")

    rows = store.list_consent_status("v1")
    assert len(rows) == 1
    assert rows[0]["username"] == "dana"
    assert rows[0]["has_consent"] is True
    assert rows[0]["agreed_at"] is not None

    pending = store.list_consent_status("v2")
    assert pending[0]["has_consent"] is False


def test_cli_bootstrap_admin(security_db: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["users", "bootstrap-admin", "--username", "admin"],
        input="s3cret\ns3cret\n",
    )
    assert result.exit_code == 0, result.output

    store = SecurityStore(db_path=security_db)
    auth = AuthService(store)
    user = auth.verify_login("admin", "s3cret")
    assert user is not None
    assert store.user_has_role(user.id, "admin")


def test_cli_bootstrap_admin_refuses_when_users_exist(security_db: Path) -> None:
    store = SecurityStore(db_path=security_db)
    store.create_user("existing", "hash")

    runner = CliRunner()
    result = runner.invoke(main, ["users", "bootstrap-admin"], input="x\nx\n")
    assert result.exit_code != 0
    assert "already exist" in result.output.lower()


def test_cli_consent_status(security_db: Path) -> None:
    store = SecurityStore(db_path=security_db)
    user_id = store.create_user("erin", "hash")
    store.record_consent(user_id, get_consent_version())

    runner = CliRunner()
    result = runner.invoke(main, ["users", "consent-status"])
    assert result.exit_code == 0
    assert "erin" in result.output
    assert "accepted" in result.output


def test_cli_audit_list_after_login_event(security_db: Path) -> None:
    store = SecurityStore(db_path=security_db)
    store.append_audit_event(
        event_type="auth.login.success",
        result="success",
        user_id=None,
        target_type="user",
        target_id="admin",
    )

    runner = CliRunner()
    result = runner.invoke(main, ["audit", "list", "--event", "auth.login.success"])
    assert result.exit_code == 0
    assert "auth.login.success" in result.output
