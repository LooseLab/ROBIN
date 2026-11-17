from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import robin.build_info as build_info


def test_get_git_commit_from_env(monkeypatch) -> None:
    build_info.get_git_commit.cache_clear()
    monkeypatch.setenv("ROBIN_GIT_COMMIT", "deadbeef" * 5)
    monkeypatch.setattr(build_info, "_find_git_root", lambda _start: None)
    assert build_info.get_git_commit() == "deadbeef" * 5


def test_get_git_commit_from_git(monkeypatch, tmp_path: Path) -> None:
    build_info.get_git_commit.cache_clear()
    monkeypatch.delenv("ROBIN_GIT_COMMIT", raising=False)
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / ".git").mkdir()

    def fake_run(cmd, cwd, capture_output, text, timeout, check):
        assert cmd == ["git", "rev-parse", "HEAD"]
        assert cwd == repo
        result = MagicMock()
        result.returncode = 0
        result.stdout = "abc123def456\n"
        return result

    monkeypatch.setattr(build_info.subprocess, "run", fake_run)
    monkeypatch.setattr(build_info, "_find_git_root", lambda _start: repo)

    assert build_info.get_git_commit() == "abc123def456"


def test_get_git_commit_returns_empty_when_unknown(monkeypatch) -> None:
    build_info.get_git_commit.cache_clear()
    monkeypatch.delenv("ROBIN_GIT_COMMIT", raising=False)
    monkeypatch.setattr(build_info, "_find_git_root", lambda _start: None)
    assert build_info.get_git_commit() == ""
