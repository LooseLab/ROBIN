"""Launch readfish alongside a MinKNOW protocol run."""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import tomli_w

from robin.minknow.auth import MinKnowAuthConfig
from robin.minknow.preset import RobinRunPreset
from robin.readfish.config import ReadfishConfig, resolve_minimap2_index
from robin.readfish.toml_builder import build_readfish_toml_document

LOGGER = logging.getLogger(__name__)


class ReadfishStartError(RuntimeError):
    """Raised when readfish cannot be started."""


@dataclass(frozen=True)
class ReadfishStartResult:
    """Outcome of launching readfish targets."""

    pid: int
    toml_path: str
    log_file: str
    command: tuple[str, ...]


def start_readfish_targets(
    *,
    preset: RobinRunPreset,
    config: ReadfishConfig,
    auth: MinKnowAuthConfig,
    position: str,
    sample_id: str,
    experiment_group: str,
    output_dir: Optional[Path] = None,
) -> ReadfishStartResult:
    """Generate readfish TOML and start ``readfish targets`` in the background."""
    if not preset.readfish_adaptive_sampling_enabled():
        raise ReadfishStartError("Preset is not configured for readfish adaptive sampling")

    executable = shutil.which(config.readfish_executable)
    if not executable:
        raise ReadfishStartError(
            f"readfish executable {config.readfish_executable!r} not found on PATH"
        )

    targets_bed = preset.effective_read_until_bed_file()
    reference = preset.effective_read_until_reference()
    if not targets_bed or not reference:
        raise ReadfishStartError(
            "readfish requires panel BED and alignment reference on the readfish host"
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

    run_output_dir = (output_dir or Path.cwd()).expanduser()
    run_output_dir.mkdir(parents=True, exist_ok=True)
    toml_path = run_output_dir / f"readfish_{_safe_name(sample_id)}.toml"
    log_file = config.resolve_log_file(sample_id=sample_id, output_dir=run_output_dir)

    with toml_path.open("wb") as handle:
        tomli_w.dump(document, handle)

    if config.validate_on_start:
        _validate_readfish_toml(
            executable=executable,
            toml_path=toml_path,
            prom=config.prom,
        )

    command = _build_readfish_command(
        executable=executable,
        toml_path=toml_path,
        position=position,
        experiment_group=experiment_group,
        log_file=log_file,
        auth=auth,
        prom=config.prom,
    )

    LOGGER.info("Starting readfish: %s", " ".join(command))
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
    except OSError as exc:
        raise ReadfishStartError(f"Failed to start readfish: {exc}") from exc

    if config.live_updates_enabled:
        from robin.readfish.live_updater import ReadfishLiveRegistry, ReadfishLiveSession

        ReadfishLiveRegistry.register(
            ReadfishLiveSession(
                sample_id=sample_id,
                base_toml_path=str(toml_path),
                region_name=config.live_region_name,
            )
        )

    return ReadfishStartResult(
        pid=process.pid,
        toml_path=str(toml_path),
        log_file=str(log_file),
        command=command,
    )


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
