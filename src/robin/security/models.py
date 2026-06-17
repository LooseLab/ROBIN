from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class User:
    id: int
    username: str
    password_hash: str
    is_active: bool
    created_at: str
    last_login_at: Optional[str]


@dataclass(frozen=True)
class UserPublic:
    id: int
    username: str
    is_active: bool
    created_at: str
    last_login_at: Optional[str]
