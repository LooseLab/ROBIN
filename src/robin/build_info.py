"""Build metadata helpers for ROBIN provenance in reports and exports."""

from __future__ import annotations

import os
import subprocess
from functools import lru_cache
from pathlib import Path

_PKG_ROOT = Path(__file__).resolve().parent


def _find_git_root(start: Path) -> Path | None:
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


@lru_cache(maxsize=1)
def get_git_commit() -> str:
    """Return the ROBIN git commit hash used for this install, if known."""
    env_commit = os.environ.get("ROBIN_GIT_COMMIT", "").strip()
    if env_commit:
        return env_commit

    git_root = _find_git_root(_PKG_ROOT)
    if git_root is None:
        return ""

    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=git_root,
            capture_output=True,
            text=True,
            timeout=2,
            check=False,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return ""
