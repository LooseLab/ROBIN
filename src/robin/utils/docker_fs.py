"""Filesystem helpers for Docker-created output (root-owned bind mounts)."""

from __future__ import annotations

import errno
import logging
import os
import shutil
import stat
import subprocess
import time
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


def _is_permission_denied(exc: BaseException) -> bool:
    if isinstance(exc, PermissionError):
        return True
    return isinstance(exc, OSError) and exc.errno in (errno.EACCES, errno.EPERM)


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


def _remove_tree_via_docker(
    parent: Path,
    name: str,
    target: Path,
    *,
    timeout_s: int,
) -> None:
    if not _docker_cli_available():
        raise PermissionError(
            f"Permission denied removing {target} and Docker is not available "
            "to remove root-owned Clair3 (or other Docker) outputs."
        )

    logger.info("Removing %s via Docker", target)
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
    logger.info("Normalizing ownership of %s to %s via Docker", target, user_spec)
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


def _make_tree_writable_for_user(path: Path, uid: int) -> None:
    """Best-effort chmod so the current user can unlink entries they own."""
    root = path.expanduser().resolve()
    if not root.exists():
        return
    for dirpath, dirnames, filenames in os.walk(root, topdown=False):
        for name in dirnames + filenames:
            entry = Path(dirpath) / name
            try:
                entry_stat = entry.stat()
            except OSError:
                continue
            if entry_stat.st_uid != uid:
                continue
            mode = entry_stat.st_mode
            desired = mode | stat.S_IWUSR
            if entry.is_dir():
                desired |= stat.S_IXUSR
            if mode != desired:
                try:
                    os.chmod(entry, desired)
                except OSError:
                    pass
    try:
        root_stat = root.stat()
        if root_stat.st_uid == uid:
            desired = root_stat.st_mode | stat.S_IWUSR | stat.S_IXUSR
            if root_stat.st_mode != desired:
                os.chmod(root, desired)
    except OSError:
        pass


def remove_tree(path: Path, *, timeout_s: int = 600, max_attempts: int = 3) -> None:
    """Remove a directory tree, falling back to Docker when permission is denied."""
    target, name = _validated_child_name(path)
    parent = target.parent
    uid = os.getuid()
    last_exc: Optional[BaseException] = None

    if tree_has_foreign_ownership(target):
        logger.warning(
            "Foreign-owned entries detected under %s; using Docker removal",
            target,
        )
        _remove_tree_via_docker(parent, name, target, timeout_s=timeout_s)
        return

    for attempt in range(max_attempts):
        try:
            if attempt:
                time.sleep(0.25 * attempt)
            _make_tree_writable_for_user(target, uid)
            shutil.rmtree(target)
            if not target.exists():
                return
            last_exc = OSError(f"Path still exists after removal attempt: {target}")
        except OSError as exc:
            last_exc = exc
            if not _is_permission_denied(exc):
                raise
            logger.warning(
                "Permission denied removing %s (attempt %s/%s)",
                target,
                attempt + 1,
                max_attempts,
            )

    if last_exc is not None and not _is_permission_denied(last_exc):
        raise OSError(str(last_exc)) from last_exc

    logger.warning(
        "Permission denied removing %s with shutil; trying Docker fallback",
        target,
    )
    _remove_tree_via_docker(parent, name, target, timeout_s=timeout_s)
