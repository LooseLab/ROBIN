"""Filesystem helpers for Docker-created output (root-owned bind mounts)."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_ALPINE_IMAGE = "alpine:latest"


def docker_host_user_spec() -> Optional[str]:
    """Return ``uid:gid`` for the current process, or None on unsupported platforms."""
    if os.name != "posix":
        return None
    try:
        return f"{os.getuid()}:{os.getgid()}"
    except AttributeError:
        return None


def _docker_cli_available() -> bool:
    try:
        result = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _validated_child_name(path: Path) -> tuple[Path, str]:
    target = path.expanduser().resolve()
    if not target.exists():
        raise FileNotFoundError(f"Path not found: {target}")
    name = target.name
    if not name or name in {".", ".."} or "/" in name or "\\" in name:
        raise ValueError(f"Unsafe path name for Docker cleanup: {name!r}")
    return target, name


def tree_has_foreign_ownership(path: Path, *, uid: Optional[int] = None) -> bool:
    """Return True when any entry under ``path`` is not owned by ``uid`` (default: current user)."""
    if uid is None:
        uid = os.getuid()
    root = path.expanduser().resolve()
    if not root.exists():
        return False
    try:
        if root.stat().st_uid != uid:
            return True
    except OSError:
        return True
    try:
        for entry in root.rglob("*"):
            try:
                if entry.stat().st_uid != uid:
                    return True
            except OSError:
                return True
    except OSError:
        return True
    return False


def chown_tree_to_host_user(path: Path, *, timeout_s: int = 600) -> bool:
    """Recursively chown ``path`` to the current user, using Docker when needed.

    Returns True when ownership was normalized (or already correct), False when
    Docker was unavailable and foreign-owned entries remain.
    """
    target, name = _validated_child_name(path)
    host_uid = os.getuid()
    if not tree_has_foreign_ownership(target, uid=host_uid):
        return True

    user_spec = docker_host_user_spec()
    if user_spec is None:
        logger.warning("Cannot chown %s: non-POSIX platform", target)
        return False

    if not _docker_cli_available():
        logger.warning(
            "Cannot chown %s: Docker is not available for ownership fixup", target
        )
        return False

    parent = target.parent
    logger.info(
        "Normalizing ownership of %s to %s via Docker", target, user_spec
    )
    result = subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{parent}:/work",
            _ALPINE_IMAGE,
            "chown",
            "-R",
            user_spec,
            f"/work/{name}",
        ],
        capture_output=True,
        text=True,
        timeout=timeout_s,
        check=False,
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        logger.error("Docker chown failed for %s: %s", target, detail)
        return False
    return not tree_has_foreign_ownership(target, uid=host_uid)


def remove_tree(path: Path, *, timeout_s: int = 600) -> None:
    """Remove a directory tree, falling back to Docker when permission is denied."""
    target, name = _validated_child_name(path)
    parent = target.parent

    def _onerror(func, bad_path, exc_info):
        exc = exc_info[1]
        if isinstance(exc, PermissionError):
            raise exc
        raise exc

    try:
        shutil.rmtree(target, onerror=_onerror)
        return
    except PermissionError:
        logger.warning(
            "Permission denied removing %s with shutil; trying Docker fallback",
            target,
        )

    if not _docker_cli_available():
        raise PermissionError(
            f"Permission denied removing {target} and Docker is not available "
            "to remove root-owned Clair3 (or other Docker) outputs."
        )

    result = subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{parent}:/work",
            _ALPINE_IMAGE,
            "rm",
            "-rf",
            f"/work/{name}",
        ],
        capture_output=True,
        text=True,
        timeout=timeout_s,
        check=False,
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        raise PermissionError(
            f"Failed to remove {target} via Docker: {detail}"
        ) from None
    if target.exists():
        raise OSError(f"Path still exists after Docker removal: {target}")
