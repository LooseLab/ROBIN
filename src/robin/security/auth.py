from __future__ import annotations

from pathlib import Path
from typing import Optional

from argon2 import PasswordHasher
from argon2.exceptions import InvalidHashError, VerifyMismatchError

from .models import User
from .store import SecurityStore


class AuthService:
    def __init__(self, store: SecurityStore):
        self.store = store
        self._hasher = PasswordHasher()

    def hash_password(self, password: str) -> str:
        return self._hasher.hash(password)

    def create_user(
        self,
        username: str,
        password: str,
        *,
        role: str = "user",
        must_change_password: bool = True,
        metadata: Optional[dict] = None,
        approvals: Optional[dict] = None,
    ) -> int:
        user_id = self.store.create_user(
            username=username,
            password_hash=self.hash_password(password),
            must_change_password=must_change_password,
            metadata=metadata,
            approvals=approvals,
        )
        self.store.assign_role(user_id, role)
        return user_id

    def verify_password(self, user: User, password: str) -> bool:
        if not password:
            return False
        try:
            return bool(self._hasher.verify(user.password_hash, password))
        except (VerifyMismatchError, InvalidHashError):
            return False

    def verify_login(self, username: str, password: str) -> Optional[User]:
        user = self.store.get_user_by_username(username)
        if user is None or not user.is_active or not password:
            return None
        if not self.verify_password(user, password):
            return None
        self.store.set_last_login(user.id)
        return self.store.get_user_by_username(username)

    def change_password(
        self,
        user_id: int,
        new_password: str,
        *,
        current_password: Optional[str] = None,
    ) -> None:
        user = self.store.get_user_by_id(user_id)
        if user is None or not user.is_active:
            raise ValueError("User not found or inactive")
        new_password = str(new_password or "")
        if not new_password:
            raise ValueError("Password cannot be empty")

        if user.must_change_password:
            pass
        else:
            if not current_password:
                raise ValueError("Current password is required")
            if not self.verify_password(user, current_password):
                raise ValueError("Current password is incorrect")

        self.store.set_user_password_hash(
            user.username,
            self.hash_password(new_password),
            must_change_password=False,
        )

    def bootstrap_admin_from_legacy_hash(self, legacy_hash_path: Path) -> bool:
        """Create initial admin user from legacy GUI hash file when needed."""
        if self.store.has_users():
            return True
        if not legacy_hash_path.exists():
            return False
        try:
            legacy_hash = legacy_hash_path.read_text(encoding="utf-8").strip()
        except OSError:
            return False
        if not legacy_hash:
            return False
        user_id = self.store.create_user(
            "admin",
            legacy_hash,
            is_active=True,
            must_change_password=False,
        )
        self.store.assign_role(user_id, "admin")
        return True
