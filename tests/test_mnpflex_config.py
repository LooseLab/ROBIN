"""Tests for MNP-Flex configuration."""

from __future__ import annotations

import pytest

from robin.analysis.mnpflex_config import load_mnpflex_config


def test_load_mnpflex_config_api_from_credentials(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("MNPFLEX_BACKEND", raising=False)
    monkeypatch.delenv("MNPFLEX_DOCKER_IMAGE", raising=False)
    monkeypatch.setenv("MNPFLEX_USERNAME", "user")
    monkeypatch.setenv("MNPFLEX_PASSWORD", "pass")
    config = load_mnpflex_config()
    assert config.backend == "api"
    assert config.validation_error() is None


def test_load_mnpflex_config_docker_requires_image(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("MNPFLEX_BACKEND", "docker")
    monkeypatch.delenv("MNPFLEX_DOCKER_IMAGE", raising=False)
    config = load_mnpflex_config()
    assert config.backend == "docker"
    assert "MNPFLEX_DOCKER_IMAGE" in (config.validation_error() or "")


def test_load_mnpflex_config_docker_with_image(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MNPFLEX_BACKEND", "docker")
    monkeypatch.setenv("MNPFLEX_DOCKER_IMAGE", "mnpflex-synnovis:1.0.0")
    config = load_mnpflex_config()
    assert config.backend == "docker"
    assert config.docker_image == "mnpflex-synnovis:1.0.0"
    assert config.validation_error() is None


def test_load_mnpflex_config_disabled(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MNPFLEX_BACKEND", "disabled")
    monkeypatch.delenv("MNPFLEX_USERNAME", raising=False)
    monkeypatch.delenv("MNPFLEX_PASSWORD", raising=False)
    config = load_mnpflex_config()
    assert config.backend == "disabled"
    assert not config.is_enabled()
