"""Shared MinKNOW CLI helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import click

from robin.minknow.auth import MinKnowAuthConfig


def build_auth_config(
    host: str,
    *,
    port: Optional[int] = None,
    api_token: Optional[str] = None,
    client_cert_chain: Optional[Path] = None,
    client_key: Optional[Path] = None,
    ca_cert: Optional[Path] = None,
    use_local_token: Optional[bool] = None,
) -> MinKnowAuthConfig:
    """Build auth config and validate cert options."""
    if (client_cert_chain is None) ^ (client_key is None):
        raise click.BadParameter(
            "--client-cert-chain and --client-key must be provided together."
        )
    return MinKnowAuthConfig.from_env(
        host=host,
        port=port,
        developer_api_token=api_token,
        client_cert_chain_path=client_cert_chain,
        client_key_path=client_key,
        ca_cert_path=ca_cert,
        use_local_token=use_local_token,
    )


def auth_click_options():
    """Return shared Click options for MinKNOW host authentication."""

    def decorator(func):
        options = [
            click.option(
                "--port",
                type=int,
                default=None,
                help="Manager API port (default: 9502, or 9501 with client certificates).",
            ),
            click.option(
                "--api-token",
                default=None,
                help="Developer API token from MinKNOW Host Settings (for remote access).",
            ),
            click.option(
                "--client-cert-chain",
                type=click.Path(exists=True, dir_okay=False, path_type=Path),
                default=None,
                help="PEM client certificate chain for authentication.",
            ),
            click.option(
                "--client-key",
                type=click.Path(exists=True, dir_okay=False, path_type=Path),
                default=None,
                help="PEM private key for the client certificate.",
            ),
            click.option(
                "--ca-cert",
                type=click.Path(exists=True, dir_okay=False, path_type=Path),
                default=None,
                help="Trusted CA certificate (remote MinKNOW installs).",
            ),
            click.option(
                "--use-local-token/--no-local-token",
                default=None,
                help=(
                    "Use MinKNOW local guest token. Defaults to off for remote hosts "
                    "(recommended when connecting by IP)."
                ),
            ),
        ]
        for option in reversed(options):
            func = option(func)
        return func

    return decorator
