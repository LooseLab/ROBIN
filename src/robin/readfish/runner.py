"""Launch readfish alongside a MinKNOW protocol run."""

from __future__ import annotations

import logging
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import tomli_w

from robin.minknow.auth import MinKnowAuthConfig
from robin.minknow.preset import RobinRunPreset
from robin.readfish.config import ReadfishConfig, resolve_minimap2_index
from robin.readfish.toml_builder import build_readfish_toml_document

LOGGER = logging.getLogger(__name__)

# Brief wait after Popen so immediate crashes are reported instead of silent failure.
_POST_START_CHECK_SECONDS = 1.5


class ReadfishStartError(RuntimeError):
    """Raised when readfish cannot be started."""


@dataclass(frozen=True)
class ReadfishStartResult:
    """Outcome of launching readfish targets."""

    pid: int
    toml_path: str
    log_file: str
    command: tuple[str, ...]
    dorado_address: str
    dorado_config: str
    live_updates_enabled: bool


@dataclass(frozen=True)
class ReadfishPrepareResult:
    """Outcome of writing a readfish experiment TOML without launching."""

    toml_path: str
    live_toml_path: Optional[str]
    targets_bed: str
    minimap2_index: str
    dorado_address: str
    dorado_config: str
    live_updates_registered: bool


def _announce(message: str) -> None:
    """Surface readfish progress on stdout and in logs (CLI / workflow terminal)."""
    text = f"[readfish] {message}"
    print(text, flush=True)
    LOGGER.info("%s", text)


def prepare_readfish_toml(
    *,
    preset: RobinRunPreset,
    config: ReadfishConfig,
    sample_id: str,
    output_dir: Optional[Path] = None,
    work_directory: Optional[Path] = None,
    register_live: bool = True,
    master_bed_path: Optional[Path] = None,
    work_dir: Optional[Path] = None,
    validate: bool = False,
) -> ReadfishPrepareResult:
    """Write a readfish experiment TOML without connecting to MinKNOW.

    Useful for offline testing of TOML generation and master-BED → ``*_live`` updates.
    When ``work_directory`` is set, files go under ``{work_directory}/{sample_id}/``.
    """
    if not preset.readfish_adaptive_sampling_enabled():
        raise ReadfishStartError("Preset is not configured for readfish adaptive sampling")

    targets_bed = preset.effective_read_until_bed_file()
    reference = preset.effective_read_until_reference()
    if not targets_bed or not reference:
        raise ReadfishStartError(
            "readfish requires panel BED and alignment reference "
            "(set alignment_reference / bed_file or workflow reference / target_panel)"
        )

    from robin.readfish.live_updater import (
        resolve_readfish_toml_path,
        resolve_sample_readfish_dir,
    )

    minimap2_index = resolve_minimap2_index(
        alignment_reference=reference,
        explicit_index=config.minimap2_index,
    )
    document = build_readfish_toml_document(
        preset=preset,
        targets_bed=targets_bed,
        minimap2_index=minimap2_index,
        config=config,
    )
    dorado = document["caller_settings"]["dorado"]
    dorado_address = str(dorado["address"])
    dorado_config = str(dorado["config"])

    effective_work = work_directory if work_directory is not None else work_dir
    run_output_dir = resolve_sample_readfish_dir(
        sample_id=sample_id,
        work_directory=effective_work,
        output_dir=output_dir,
    )
    run_output_dir.mkdir(parents=True, exist_ok=True)
    toml_path = resolve_readfish_toml_path(
        sample_id=sample_id,
        work_directory=effective_work,
        output_dir=output_dir,
    )

    with toml_path.open("wb") as handle:
        tomli_w.dump(document, handle)

    _announce(f"Wrote experiment TOML (offline): {toml_path}")
    _announce(f"Dorado address: {dorado_address}")
    _announce(f"Dorado config:  {dorado_config}")
    _announce(f"Targets BED:    {targets_bed}")
    _announce(f"Mapper index:   {minimap2_index}")

    if validate:
        executable = shutil.which(config.readfish_executable)
        if not executable:
            raise ReadfishStartError(
                f"readfish executable {config.readfish_executable!r} not found on PATH"
            )
        _announce(f"Running: {executable} validate {toml_path}")
        _validate_readfish_toml(
            executable=executable,
            toml_path=toml_path,
            prom=config.prom,
        )
        _announce("Validate succeeded")

    live_registered = False
    live_path: Optional[Path] = None
    should_register = register_live and config.live_updates_enabled
    if should_register:
        from robin.readfish.live_updater import (
            ReadfishLiveRegistry,
            ReadfishLiveSession,
            live_toml_path,
            notify_readfish_live_targets,
        )

        ReadfishLiveRegistry.register(
            ReadfishLiveSession(
                sample_id=sample_id,
                base_toml_path=str(toml_path),
                region_name=config.live_region_name,
            )
        )
        live_registered = True
        _announce(
            f"Registered live target updates for sample {sample_id!r} "
            f"(no MinKNOW connection required)"
        )

        if master_bed_path is not None or work_dir is not None:
            live_path = notify_readfish_live_targets(
                sample_id=sample_id,
                master_bed_path=master_bed_path,
                work_dir=work_dir,
            )
        else:
            live_path = live_toml_path(toml_path)
            _announce(
                f"Live TOML will appear at: {live_path} "
                "(after master BED notify / generate_master_bed)"
            )
    elif not config.live_updates_enabled:
        _announce("Live updates disabled in [readfish] (live_updates_enabled=false)")

    return ReadfishPrepareResult(
        toml_path=str(toml_path),
        live_toml_path=str(live_path) if live_path is not None else None,
        targets_bed=str(targets_bed),
        minimap2_index=str(minimap2_index),
        dorado_address=dorado_address,
        dorado_config=dorado_config,
        live_updates_registered=live_registered,
    )


def start_readfish_targets(
    *,
    preset: RobinRunPreset,
    config: ReadfishConfig,
    auth: MinKnowAuthConfig,
    position: str,
    sample_id: str,
    experiment_group: str,
    output_dir: Optional[Path] = None,
    work_directory: Optional[Path] = None,
) -> ReadfishStartResult:
    """Generate readfish TOML and start ``readfish targets`` in the background."""
    if not preset.readfish_adaptive_sampling_enabled():
        raise ReadfishStartError("Preset is not configured for readfish adaptive sampling")

    _announce("Adaptive sampling backend is readfish — preparing launch")

    executable = shutil.which(config.readfish_executable)
    if not executable:
        raise ReadfishStartError(
            f"readfish executable {config.readfish_executable!r} not found on PATH"
        )
    _announce(f"Executable: {executable}")

    targets_bed = preset.effective_read_until_bed_file()
    reference = preset.effective_read_until_reference()
    if not targets_bed or not reference:
        raise ReadfishStartError(
            "readfish requires panel BED and alignment reference on the readfish host"
        )

    from robin.readfish.live_updater import (
        resolve_readfish_toml_path,
        resolve_sample_readfish_dir,
    )

    minimap2_index = resolve_minimap2_index(
        alignment_reference=reference,
        explicit_index=config.minimap2_index,
    )
    document = build_readfish_toml_document(
        preset=preset,
        targets_bed=targets_bed,
        minimap2_index=minimap2_index,
        config=config,
    )
    dorado = document["caller_settings"]["dorado"]
    dorado_address = str(dorado["address"])
    dorado_config = str(dorado["config"])

    run_output_dir = resolve_sample_readfish_dir(
        sample_id=sample_id,
        work_directory=work_directory,
        output_dir=output_dir,
    )
    run_output_dir.mkdir(parents=True, exist_ok=True)
    toml_path = resolve_readfish_toml_path(
        sample_id=sample_id,
        work_directory=work_directory,
        output_dir=output_dir,
    )
    log_file = config.resolve_log_file(sample_id=sample_id, output_dir=run_output_dir)

    with toml_path.open("wb") as handle:
        tomli_w.dump(document, handle)

    _announce(f"Wrote experiment TOML: {toml_path}")
    _announce(f"Dorado address: {dorado_address}")
    _announce(f"Dorado config:  {dorado_config}")
    _announce(f"Targets BED:    {targets_bed}")
    _announce(f"Mapper index:   {minimap2_index}")
    _announce(f"Position:       {position}")
    _announce(f"Experiment:     {experiment_group}")
    _announce(f"Log file:       {log_file}")
    _announce(
        "Live updates:   "
        + ("enabled" if config.live_updates_enabled else "disabled")
    )

    if config.validate_on_start:
        _announce(f"Running: {executable} validate {toml_path}")
        _validate_readfish_toml(
            executable=executable,
            toml_path=toml_path,
            prom=config.prom,
        )
        _announce("Validate succeeded")
    else:
        _announce("Skipping validate (validate_on_start=false)")

    command = _build_readfish_command(
        executable=executable,
        toml_path=toml_path,
        position=position,
        experiment_group=experiment_group,
        log_file=log_file,
        auth=auth,
        prom=config.prom,
    )
    command_text = " ".join(command)
    _announce(f"Launching: {command_text}")

    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
    except OSError as exc:
        raise ReadfishStartError(f"Failed to start readfish: {exc}") from exc

    _announce(f"Spawned process pid={process.pid}")
    _wait_for_startup(process, log_file=log_file)
    _announce(f"Process pid={process.pid} is still running after startup check")

    if config.live_updates_enabled:
        from robin.readfish.live_updater import ReadfishLiveRegistry, ReadfishLiveSession

        ReadfishLiveRegistry.register(
            ReadfishLiveSession(
                sample_id=sample_id,
                base_toml_path=str(toml_path),
                region_name=config.live_region_name,
            )
        )
        _announce(
            f"Registered live target updates for sample {sample_id!r} "
            f"(region={config.live_region_name!r})"
        )

    _announce(
        f"Ready — follow logs with: tail -f {log_file}"
    )
    return ReadfishStartResult(
        pid=process.pid,
        toml_path=str(toml_path),
        log_file=str(log_file),
        command=command,
        dorado_address=dorado_address,
        dorado_config=dorado_config,
        live_updates_enabled=config.live_updates_enabled,
    )


def _wait_for_startup(process: subprocess.Popen, *, log_file: Path) -> None:
    """Fail fast if readfish exits immediately after launch."""
    time.sleep(_POST_START_CHECK_SECONDS)
    exit_code = process.poll()
    if exit_code is None:
        return
    detail = _tail_log(log_file)
    message = (
        f"readfish exited immediately after start "
        f"(pid={process.pid}, exit_code={exit_code})"
    )
    if detail:
        message += f"\n--- {log_file} (tail) ---\n{detail}"
    else:
        message += (
            f"\nLog file empty or missing ({log_file}). "
            "Check Dorado address/config and MinKNOW device name."
        )
    raise ReadfishStartError(message)


def _tail_log(log_file: Path, *, max_chars: int = 4000) -> str:
    try:
        if not log_file.is_file():
            return ""
        text = log_file.read_text(encoding="utf-8", errors="replace").strip()
    except OSError:
        return ""
    if len(text) <= max_chars:
        return text
    return text[-max_chars:]


def _validate_readfish_toml(
    *,
    executable: str,
    toml_path: Path,
    prom: bool,
) -> None:
    command = [executable, "validate", str(toml_path)]
    if prom:
        command.append("--prom")
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exc:
        raise ReadfishStartError(f"Failed to run readfish validate: {exc}") from exc

    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise ReadfishStartError(
            "readfish validate failed"
            + (f": {detail}" if detail else "")
        )


def _build_readfish_command(
    *,
    executable: str,
    toml_path: Path,
    position: str,
    experiment_group: str,
    log_file: Path,
    auth: MinKnowAuthConfig,
    prom: bool,
) -> tuple[str, ...]:
    command: list[str] = [
        executable,
        "targets",
        "--toml",
        str(toml_path),
        "--device",
        position,
        "--experiment-name",
        experiment_group,
        "--log-file",
        str(log_file),
    ]
    if auth.host and auth.host not in {"", "localhost", "127.0.0.1"}:
        command.extend(["--host", auth.host])
    if auth.developer_api_token:
        command.extend(["--api-token", auth.developer_api_token])
    if prom:
        command.append("--prom")
    return tuple(command)


def _safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value) or "run"
