from __future__ import annotations

from pathlib import Path

import pytest

from robin.security import SecurityStore
from robin.security.user_metadata import (
    CLINICAL_ROLE_KEY,
    EMAIL_KEY,
    normalize_metadata,
    parse_metadata_json,
)


def test_normalize_metadata_strips_and_validates_email() -> None:
    meta = normalize_metadata(
        {
            EMAIL_KEY: "  user@example.com ",
            CLINICAL_ROLE_KEY: "Consultant",
            "unknown": "ignored",
        }
    )
    assert meta == {EMAIL_KEY: "user@example.com", CLINICAL_ROLE_KEY: "Consultant"}

    with pytest.raises(ValueError):
        normalize_metadata({EMAIL_KEY: "not-an-email"})


def test_parse_metadata_json_handles_invalid() -> None:
    assert parse_metadata_json("") == {}
    assert parse_metadata_json("not-json") == {}


def test_store_user_metadata_roundtrip(tmp_path: Path) -> None:
    store = SecurityStore(db_path=tmp_path / "security.db")
    user_id = store.create_user(
        "alice",
        "hash",
        metadata={EMAIL_KEY: "alice@lab.org", CLINICAL_ROLE_KEY: "Scientist"},
    )
    store.assign_role(user_id, "user")

    user = store.get_user_by_username("alice")
    assert user is not None
    assert user.metadata[EMAIL_KEY] == "alice@lab.org"
    assert user.metadata[CLINICAL_ROLE_KEY] == "Scientist"

    store.update_user_metadata(
        "alice",
        {CLINICAL_ROLE_KEY: "Consultant", EMAIL_KEY: ""},
    )
    updated = store.get_user_by_username("alice")
    assert updated is not None
    assert updated.metadata[CLINICAL_ROLE_KEY] == "Consultant"
    assert EMAIL_KEY not in updated.metadata

    public_users = store.list_users()
    assert public_users[0].metadata[CLINICAL_ROLE_KEY] == "Consultant"
