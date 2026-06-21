"""Configuration for MinKNOW monitoring."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional

from robin.minknow.auth import MinKnowAuthConfig


def _env_bool(value: Optional[str], *, default: bool) -> bool:
    if value is None:
        return default
    return value.strip().lower() in ("1", "true", "yes", "on")


@dataclass(frozen=True)
class MinKnowSettings:
    """Settings for monitoring a MinKNOW host from ROBIN."""

    enabled: bool
    host: str
    poll_interval_s: float
    auth: MinKnowAuthConfig
    auto_add_paths: bool = False

    @classmethod
    def from_host(
        cls,
        host: str,
        *,
        enabled: bool = True,
        poll_interval_s: float = 10.0,
        port: Optional[int] = None,
        developer_api_token: Optional[str] = None,
        client_cert_chain_path: Optional[Path] = None,
        client_key_path: Optional[Path] = None,
        ca_cert_path: Optional[Path] = None,
        use_local_token: Optional[bool] = None,
        auto_add_paths: bool = False,
        environ: Optional[Mapping[str, str]] = None,
    ) -> MinKnowSettings:
        auth = MinKnowAuthConfig.from_env(
            host=host,
            port=port,
            developer_api_token=developer_api_token,
            client_cert_chain_path=client_cert_chain_path,
            client_key_path=client_key_path,
            ca_cert_path=ca_cert_path,
            use_local_token=use_local_token,
            environ=environ,
        )
        return cls(
            enabled=enabled,
            host=host,
            poll_interval_s=poll_interval_s,
            auth=auth,
            auto_add_paths=auto_add_paths,
        )

    @classmethod
    def from_environ(
        cls, environ: Optional[Mapping[str, str]] = None
    ) -> MinKnowSettings:
        env = dict(environ or os.environ)
        host = (env.get("MINKNOW_HOST") or "localhost").strip()
        enabled = _env_bool(env.get("MINKNOW_ENABLED"), default=bool(host))
        try:
            poll_interval_s = float(env.get("MINKNOW_POLL_INTERVAL", "10"))
        except ValueError:
            poll_interval_s = 10.0
        if poll_interval_s <= 0:
            poll_interval_s = 10.0

        port = _optional_int(env.get("MINKNOW_API_PORT"))
        auto_add_paths = _env_bool(env.get("MINKNOW_AUTO_WATCH"), default=False)
        return cls.from_host(
            host=host,
            enabled=enabled,
            poll_interval_s=poll_interval_s,
            port=port,
            auto_add_paths=auto_add_paths,
            environ=env,
        )

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> MinKnowSettings:
        """Load settings from a ``[minknow]`` TOML table."""
        host = str(data.get("host") or data.get("minknow_host") or "localhost").strip()
        enabled = data.get("enabled", True)
        if not isinstance(enabled, bool):
            enabled = _env_bool(str(enabled), default=True)

        poll_interval_s = data.get("poll_interval_s", data.get("poll_interval", 10.0))
        try:
            poll_interval_s = float(poll_interval_s)
        except (TypeError, ValueError):
            poll_interval_s = 10.0
        if poll_interval_s <= 0:
            poll_interval_s = 10.0

        port = data.get("port")
        if port is not None:
            port = int(port)

        auto_add_paths = data.get("auto_add_paths", data.get("auto_watch", False))
        if not isinstance(auto_add_paths, bool):
            auto_add_paths = _env_bool(str(auto_add_paths), default=False)

        return cls.from_host(
            host=host,
            enabled=enabled,
            poll_interval_s=poll_interval_s,
            port=port,
            auto_add_paths=auto_add_paths,
            developer_api_token=_optional_str(data.get("developer_api_token")),
            client_cert_chain_path=_optional_path(data.get("client_cert_chain")),
            client_key_path=_optional_path(data.get("client_key")),
            ca_cert_path=_optional_path(data.get("ca_cert")),
            use_local_token=data.get("use_local_token"),
        )


def preset_path_from_environ(
    environ: Optional[Mapping[str, str]] = None,
) -> Optional[Path]:
    """Return preset TOML path from ``MINKNOW_PRESET`` if set."""
    env = dict(environ or os.environ)
    raw = (env.get("MINKNOW_PRESET") or "").strip()
    if not raw:
        return None
    return Path(raw).expanduser()


def workflow_toml_from_environ(
    environ: Optional[Mapping[str, str]] = None,
) -> Optional[Path]:
    """Return workflow TOML path from ``ROBIN_WORKFLOW_TOML`` if set."""
    env = dict(environ or os.environ)
    raw = (env.get("ROBIN_WORKFLOW_TOML") or "").strip()
    if not raw:
        return None
    return Path(raw).expanduser()


def _optional_int(value: Optional[str]) -> Optional[int]:
    if not value:
        return None
    return int(value)


def _optional_str(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _optional_path(value: Any) -> Optional[Path]:
    if value is None:
        return None
    return Path(str(value)).expanduser()
