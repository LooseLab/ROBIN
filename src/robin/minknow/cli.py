"""CLI commands for MinKNOW integration."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import click

from robin.minknow.auth import MinKnowAuthConfig
from robin.minknow.cli_common import auth_click_options, build_auth_config
from robin.minknow.client import MinKnowConnectionError
from robin.minknow.config import MinKnowSettings
from robin.minknow.models import SequencerStatus
from robin.minknow.monitor import fetch_sequencer_status
from robin.minknow.run import MinKnowStartError, StartRunRequest, start_protocol_run
from robin.minknow.stream_monitor import acquire_stream_monitor
from robin.minknow.toml_config import load_minknow_toml
from robin.minknow.workflow_refs import load_workflow_config_for_refs
from robin.workflow_config import load_workflow_toml
from robin.minknow.watch import process_auto_watch, watch_active_runs


@click.group()
def minknow() -> None:
    """Query and control Oxford Nanopore MinKNOW sequencers."""


@minknow.command("status")
@click.option(
    "--host",
    required=True,
    help="Hostname or IP address of the machine running MinKNOW.",
)
@click.option(
    "--port",
    type=int,
    default=None,
    help="Manager API port (default: 9502, or 9501 with client certificates).",
)
@click.option(
    "--api-token",
    default=None,
    help="Developer API token from MinKNOW Host Settings (for remote access).",
)
@click.option(
    "--client-cert-chain",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="PEM client certificate chain for authentication.",
)
@click.option(
    "--client-key",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="PEM private key for the client certificate.",
)
@click.option(
    "--ca-cert",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Trusted CA certificate (remote MinKNOW installs).",
)
@click.option(
    "--use-local-token/--no-local-token",
    default=None,
    help=(
        "Use MinKNOW local guest token. Defaults to off for remote hosts "
        "(recommended when connecting by IP)."
    ),
)
@click.option("--json", "as_json", is_flag=True, help="Emit machine-readable JSON.")
def status(
    host: str,
    port: Optional[int],
    api_token: Optional[str],
    client_cert_chain: Optional[Path],
    client_key: Optional[Path],
    ca_cert: Optional[Path],
    use_local_token: Optional[bool],
    as_json: bool,
) -> None:
    """Show flow cell positions and active protocol runs on a MinKNOW host."""
    if (client_cert_chain is None) ^ (client_key is None):
        raise click.BadParameter(
            "--client-cert-chain and --client-key must be provided together."
        )

    auth = MinKnowAuthConfig.from_env(
        host=host,
        port=port,
        developer_api_token=api_token,
        client_cert_chain_path=client_cert_chain,
        client_key_path=client_key,
        ca_cert_path=ca_cert,
        use_local_token=use_local_token,
    )

    settings = MinKnowSettings(
        enabled=True,
        host=host,
        poll_interval_s=10.0,
        auth=auth,
    )
    poll_result = fetch_sequencer_status(settings)
    if poll_result.error:
        if "requires the minknow_api package" in poll_result.error:
            click.echo(poll_result.error, err=True)
        else:
            click.echo(f"MinKNOW connection failed: {poll_result.error}", err=True)
        sys.exit(1)
    sequencer_status = poll_result.status
    if sequencer_status is None:
        click.echo("No MinKNOW status returned.", err=True)
        sys.exit(1)

    if as_json:
        click.echo(json.dumps(sequencer_status.to_dict(), indent=2))
        return

    click.echo(format_sequencer_status(sequencer_status))


@minknow.command("watch")
@click.option(
    "--host",
    required=True,
    help="Hostname or IP address of the machine running MinKNOW.",
)
@click.option(
    "--port",
    type=int,
    default=None,
    help="Manager API port (default: 9502, or 9501 with client certificates).",
)
@click.option(
    "--api-token",
    default=None,
    help="Developer API token from MinKNOW Host Settings (for remote access).",
)
@click.option(
    "--client-cert-chain",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="PEM client certificate chain for authentication.",
)
@click.option(
    "--client-key",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="PEM private key for the client certificate.",
)
@click.option(
    "--ca-cert",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Trusted CA certificate (remote MinKNOW installs).",
)
@click.option(
    "--use-local-token/--no-local-token",
    default=None,
    help=(
        "Use MinKNOW local guest token. Defaults to off for remote hosts "
        "(recommended when connecting by IP)."
    ),
)
@click.option(
    "--auto-add-paths",
    is_flag=True,
    help="Keep running and add new protocol runs to the ROBIN watch list automatically.",
)
def watch(
    host: str,
    port: Optional[int],
    api_token: Optional[str],
    client_cert_chain: Optional[Path],
    client_key: Optional[Path],
    ca_cert: Optional[Path],
    use_local_token: Optional[bool],
    auto_add_paths: bool,
) -> None:
    """Add MinKNOW run output directories to the ROBIN workflow watch list."""
    if (client_cert_chain is None) ^ (client_key is None):
        raise click.BadParameter(
            "--client-cert-chain and --client-key must be provided together."
        )

    auth = MinKnowAuthConfig.from_env(
        host=host,
        port=port,
        developer_api_token=api_token,
        client_cert_chain_path=client_cert_chain,
        client_key_path=client_key,
        ca_cert_path=ca_cert,
        use_local_token=use_local_token,
    )
    settings = MinKnowSettings(
        enabled=True,
        host=host,
        poll_interval_s=10.0,
        auth=auth,
        auto_add_paths=auto_add_paths,
    )

    if auto_add_paths:
        _watch_auto_add_paths(settings)
        return

    poll_result = fetch_sequencer_status(settings)
    if poll_result.error:
        click.echo(f"MinKNOW connection failed: {poll_result.error}", err=True)
        sys.exit(1)
    if poll_result.status is None:
        click.echo("No MinKNOW status returned.", err=True)
        sys.exit(1)

    actions = watch_active_runs(poll_result.status)
    if not actions:
        click.echo("No watchable active runs found.")
        return

    exit_code = 0
    for position_name, success, message in actions:
        prefix = "OK" if success else "FAILED"
        click.echo(f"{prefix} {position_name}: {message}")
        if not success:
            exit_code = 1
    sys.exit(exit_code)


@minknow.command("start")
@click.option(
    "--host",
    default=None,
    help="MinKNOW host (overrides [minknow].host in --preset TOML).",
)
@click.option(
    "--position",
    default=None,
    help="Flow cell position name (overrides preset default).",
)
@click.option(
    "--sample-id",
    required=True,
    help="Sample ID for the protocol run (must match ROBIN tracking).",
)
@click.option(
    "--experiment-group",
    default=None,
    help="Protocol group / experiment ID (default: ROBIN_RUN or preset value).",
)
@click.option(
    "--preset",
    "preset_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="TOML file with [minknow.preset], or ROBIN workflow TOML containing [minknow].",
)
@click.option(
    "--workflow-toml",
    "workflow_toml",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="ROBIN workflow TOML for reference and target_panel (defaults to --preset when it is a workflow file).",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Validate preset and print summary without starting a run.",
)
@click.option(
    "--check-paths",
    is_flag=True,
    help="Verify reference/BED paths exist on this machine (ROBIN host). "
    "Preset paths normally live on the MinKNOW host and are not checked by default.",
)
@click.option(
    "--skip-model-check",
    is_flag=True,
    help="Do not verify basecall models against the sequencer.",
)
@auth_click_options()
def start(
    host: Optional[str],
    position: Optional[str],
    sample_id: str,
    experiment_group: Optional[str],
    preset_path: Path,
    workflow_toml: Optional[Path],
    dry_run: bool,
    check_paths: bool,
    skip_model_check: bool,
    port: Optional[int],
    api_token: Optional[str],
    client_cert_chain: Optional[Path],
    client_key: Optional[Path],
    ca_cert: Optional[Path],
    use_local_token: Optional[bool],
) -> None:
    """Start a ROBIN-compliant sequencing run on a MinKNOW host."""
    workflow_config = None
    workflow_path = workflow_toml or preset_path
    try:
        workflow_config = load_workflow_toml(workflow_path)
    except click.BadParameter:
        workflow_config = None

    try:
        workflow_config_loaded = load_minknow_toml(
            preset_path,
            workflow_config=workflow_config,
            prefer_workflow=True,
        )
    except click.BadParameter as exc:
        raise click.ClickException(str(exc)) from exc

    preset = workflow_config_loaded.preset
    if preset is None:
        raise click.ClickException(
            f"No [minknow.preset] section found in {preset_path}"
        )

    resolved_host = (host or workflow_config_loaded.settings.host).strip()
    if not resolved_host:
        raise click.ClickException("--host is required (or set [minknow].host in preset file)")

    resolved_position = (position or preset.position or "").strip()
    if not resolved_position:
        raise click.ClickException(
            "--position is required (or set position in [minknow.preset])"
        )

    auth = build_auth_config(
        resolved_host,
        port=port,
        api_token=api_token,
        client_cert_chain=client_cert_chain,
        client_key=client_key,
        ca_cert=ca_cert,
        use_local_token=use_local_token,
    )

    validation_errors = preset.validate(check_paths=check_paths)
    if validation_errors:
        raise click.ClickException(
            "Preset validation failed:\n  - " + "\n  - ".join(validation_errors)
        )

    click.echo(f"MinKNOW host: {resolved_host}")
    click.echo(f"Position: {resolved_position}")
    click.echo(f"Sample ID: {sample_id}")
    click.echo("Preset:")
    for line in preset.summary_lines():
        click.echo(f"  {line}")

    if dry_run:
        click.echo("Dry run — protocol not started.")
        return

    request = StartRunRequest(
        preset=preset,
        position=resolved_position,
        sample_id=sample_id,
        experiment_group=experiment_group,
    )

    try:
        result = start_protocol_run(
            auth,
            request,
            validate_models=not skip_model_check,
            check_paths=check_paths,
        )
    except MinKnowConnectionError as exc:
        raise click.ClickException(f"MinKNOW connection failed: {exc}") from exc
    except MinKnowStartError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo("Started protocol run:")
    click.echo(f"  run_id: {result.run_id}")
    click.echo(f"  position: {result.position}")
    if result.flow_cell_id:
        click.echo(f"  flow_cell_id: {result.flow_cell_id}")
    click.echo(f"  protocol_id: {result.protocol_id}")
    click.echo(f"  experiment_group: {result.experiment_group}")
    for warning in result.warnings:
        click.echo(f"Warning: {warning}")


@minknow.command("stop")
@click.option(
    "--host",
    required=True,
    help="Hostname or IP address of the machine running MinKNOW.",
)
@click.option(
    "--position",
    required=True,
    help="Flow cell position name (e.g. P2S_000000-A).",
)
@click.option(
    "--run-id",
    "protocol_run_id",
    default=None,
    help="Protocol run ID to stop (defaults to the current run on the position).",
)
@click.option(
    "--wait/--no-wait",
    default=False,
    help="Wait for the protocol to finish after requesting stop.",
)
@auth_click_options()
def stop(
    host: str,
    position: str,
    protocol_run_id: Optional[str],
    wait: bool,
    port: Optional[int],
    api_token: Optional[str],
    client_cert_chain: Optional[Path],
    client_key: Optional[Path],
    ca_cert: Optional[Path],
    use_local_token: Optional[bool],
) -> None:
    """Stop a running MinKNOW protocol on a flow cell position."""
    from robin.minknow.run import MinKnowStopError, StopRunRequest, stop_protocol_run

    auth = build_auth_config(
        host,
        port=port,
        api_token=api_token,
        client_cert_chain=client_cert_chain,
        client_key=client_key,
        ca_cert=ca_cert,
        use_local_token=use_local_token,
    )

    try:
        result = stop_protocol_run(
            auth,
            StopRunRequest(position=position, protocol_run_id=protocol_run_id),
            wait=wait,
        )
    except MinKnowConnectionError as exc:
        raise click.ClickException(f"MinKNOW connection failed: {exc}") from exc
    except MinKnowStopError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Stop requested for {result.position} (run_id={result.protocol_run_id})")
    if result.waited:
        click.echo(f"Protocol finished with state: {result.protocol_state or 'unknown'}")


def _watch_auto_add_paths(settings: MinKnowSettings) -> None:
    import signal
    import time

    monitor, release = acquire_stream_monitor(settings)
    stop = False

    def _handle_signal(_signum, _frame) -> None:
        nonlocal stop
        stop = True

    previous_int = signal.getsignal(signal.SIGINT)
    previous_term = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGINT, _handle_signal)
    signal.signal(signal.SIGTERM, _handle_signal)

    def _on_update(result) -> None:
        if result.status is None:
            if result.error:
                click.echo(f"MinKNOW stream error: {result.error}", err=True)
            return
        for position_name, _success, message in process_auto_watch(result.status):
            click.echo(f"AUTO-WATCH {position_name}: {message}")

    unsubscribe = monitor.subscribe(_on_update)
    click.echo(
        f"Watching MinKNOW at {settings.host} for new runs. Press Ctrl+C to stop."
    )
    try:
        while not stop:
            time.sleep(0.5)
    finally:
        unsubscribe()
        release()
        signal.signal(signal.SIGINT, previous_int)
        signal.signal(signal.SIGTERM, previous_term)


def format_sequencer_status(status: SequencerStatus) -> str:
    """Render sequencer status for terminal output."""
    lines = [
        f"MinKNOW host: {status.host}:{status.port}",
        f"MinKNOW Core: {status.core_version}",
        f"Distribution: {status.distribution_version}",
        f"minknow_api package: {status.minknow_api_version}",
    ]
    if status.version_warning:
        lines.append(f"Warning: {status.version_warning}")

    if not status.positions:
        lines.append("")
        lines.append("No flow cell positions reported.")
        return "\n".join(lines)

    lines.append("")
    lines.append(f"Positions ({len(status.positions)}):")
    for position in status.positions:
        lines.append("")
        lines.extend(_format_position(position))
    return "\n".join(lines)


def _format_position(position) -> list[str]:
    lines = [
        f"  {position.name}",
        f"    state: {position.state}",
        f"    protocol_state: {position.protocol_state}",
        f"    software_running: {position.running}",
    ]
    if position.device_type:
        lines.append(f"    device_type: {position.device_type}")
    if position.flow_cell_id:
        lines.append(f"    flow_cell_id: {position.flow_cell_id}")
    if position.flow_cell_product_code:
        lines.append(f"    flow_cell_product: {position.flow_cell_product_code}")
    if position.sample_id:
        lines.append(f"    sample_id: {position.sample_id}")
    if position.protocol_group_id:
        lines.append(f"    protocol_group_id: {position.protocol_group_id}")
    if position.protocol_run_id:
        lines.append(f"    protocol_run_id: {position.protocol_run_id}")
    if position.protocol_name:
        lines.append(f"    protocol_id: {position.protocol_name}")
    if position.protocol_run_state:
        lines.append(f"    protocol_run_state: {position.protocol_run_state}")
    if position.output_path:
        lines.append(f"    output_path: {position.output_path}")
    if position.output_reads_path:
        lines.append(f"    reads_path: {position.output_reads_path}")
    if position.output_logs_path:
        lines.append(f"    logs_path: {position.output_logs_path}")
    if position.connection_error:
        lines.append(f"    connection_error: {position.connection_error}")
    return lines
