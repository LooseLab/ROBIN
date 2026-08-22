"""Authentication configuration for remote MinKNOW connections."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, MutableMapping, Optional

LOCAL_HOSTS = frozenset({"localhost", "127.0.0.1", "::1"})


@dataclass(frozen=True)
class MinKnowAuthConfig:
    """Connection and authentication settings for ``minknow_api.manager.Manager``."""

    host: str
    port: Optional[int] = None
    developer_api_token: Optional[str] = None
    client_cert_chain_path: Optional[Path] = None
    client_key_path: Optional[Path] = None
    ca_cert_path: Optional[Path] = None
    use_local_token: Optional[bool] = None

    @classmethod
    def from_env(
        cls,
        host: str,
        *,
        port: Optional[int] = None,
        developer_api_token: Optional[str] = None,
        client_cert_chain_path: Optional[Path] = None,
        client_key_path: Optional[Path] = None,
        ca_cert_path: Optional[Path] = None,
        use_local_token: Optional[bool] = None,
        environ: Optional[Mapping[str, str]] = None,
    ) -> MinKnowAuthConfig:
        env = environ or os.environ
        return cls(
            host=host,
            port=port or _optional_int(env.get("MINKNOW_API_PORT")),
            developer_api_token=(
                developer_api_token
                or env.get("MINKNOW_DEVELOPER_API_TOKEN")
                or env.get("MINKNOW_API_TOKEN")
            ),
            client_cert_chain_path=(
                client_cert_chain_path
                or _optional_path(env.get("MINKNOW_API_CLIENT_CERTIFICATE_CHAIN"))
            ),
            client_key_path=(
                client_key_path or _optional_path(env.get("MINKNOW_API_CLIENT_KEY"))
            ),
            ca_cert_path=(
                ca_cert_path or _optional_path(env.get("MINKNOW_TRUSTED_CA"))
            ),
            use_local_token=_resolve_use_local_token(host, use_local_token, env),
        )

    def manager_kwargs(self, environ: Optional[Mapping[str, str]] = None) -> dict:
        """Build keyword arguments for ``Manager(...)``."""
        merged: MutableMapping[str, str] = dict(environ or os.environ)
        if self.use_local_token is not None:
            merged["MINKNOW_API_USE_LOCAL_TOKEN"] = "1" if self.use_local_token else "0"

        kwargs: dict = {
            "host": self.host,
            "environ": merged,
        }
        if self.port is not None:
            kwargs["port"] = self.port
        if self.developer_api_token:
            kwargs["developer_api_token"] = self.developer_api_token
        if self.client_cert_chain_path is not None:
            kwargs["client_certificate_chain"] = (
                self.client_cert_chain_path.read_bytes()
            )
        if self.client_key_path is not None:
            kwargs["client_private_key"] = self.client_key_path.read_bytes()
        if self.ca_cert_path is not None:
            kwargs["ca_certificate"] = self.ca_cert_path.read_bytes()
        return kwargs


def _optional_path(value: Optional[str]) -> Optional[Path]:
    if not value:
        return None
    return Path(value).expanduser()


def _optional_int(value: Optional[str]) -> Optional[int]:
    if not value:
        return None
    return int(value)


def _resolve_use_local_token(
    host: str,
    explicit: Optional[bool],
    environ: Mapping[str, str],
) -> Optional[bool]:
    if explicit is not None:
        return explicit
    if "MINKNOW_API_USE_LOCAL_TOKEN" in environ:
        return environ["MINKNOW_API_USE_LOCAL_TOKEN"].strip().lower() in (
            "1",
            "true",
            "yes",
            "on",
        )
    if host.strip().lower() not in LOCAL_HOSTS:
        return False
    return None
