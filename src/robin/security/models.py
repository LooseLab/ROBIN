from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass(frozen=True)
class User:
    id: int
    username: str
    password_hash: str
    is_active: bool
    created_at: str
    last_login_at: Optional[str]
    must_change_password: bool = False
    metadata: Dict[str, str] = field(default_factory=dict)
    approvals: Dict[str, bool] = field(default_factory=dict)


@dataclass(frozen=True)
class UserPublic:
    id: int
    username: str
    is_active: bool
    created_at: str
    last_login_at: Optional[str]
    must_change_password: bool = False
    metadata: Dict[str, str] = field(default_factory=dict)
    approvals: Dict[str, bool] = field(default_factory=dict)
