"""Load and merge TOML configuration for ``robin workflow``."""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any, Mapping, MutableMapping, Optional, Sequence

import click
from click.core import ParameterSource

from robin.minknow.toml_config import MinKnowWorkflowConfig, extract_minknow_config, load_minknow_toml

WORKFLOW_REQUIRED_KEYS = ("path", "workflow", "center", "target_panel")

_PATH_KEYS = frozenset({"path", "work_dir", "reference", "toml"})
_LIST_KEYS = frozenset(
    {"commands", "job_log_level", "deduplicate_jobs", "queue_priority"}
)
_BOOL_KEYS = frozenset(
    {
        "verbose",
        "no_process_existing",
        "no_progress",
        "legacy_analysis_queue",
        "use_ray",
        "use_ray_core",
        "show_priorities",
        "no_watch",
        "with_gui",
        "ray_dashboard",
    }
)
_INT_KEYS = frozenset(
    {
        "analysis_workers",
        "preprocessing_workers",
        "bed_workers",
        "ray_num_cpus",
        "gui_port",
    }
)

WORKFLOW_CONFIG_EXAMPLE = """\
# ROBIN workflow configuration
# Run with: robin workflow --toml my_settings.toml
# CLI flags override values from this file when explicitly provided.

path = "empty_folder"
workflow = "cnv,fusion,target,mgmt,sturgeon,nanodx,pannanodx,random_forest"
center = "NUH"
target_panel = "rCNS2"

work_dir = "../../REF_SAMPLES"
reference = "~/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Optional settings (defaults match CLI)
# log_level = "ERROR"
# analysis_workers = 1
# preset = "standard"
# with_gui = true
# use_ray = true
"""


def load_workflow_toml(path: Path) -> dict[str, Any]:
    """Load workflow settings from a TOML file."""
    if not path.exists():
        raise click.BadParameter(f"TOML config file does not exist: {path}")
    if not path.is_file():
        raise click.BadParameter(f"TOML config path is not a file: {path}")

    try:
        with path.open("rb") as handle:
            raw = tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise click.BadParameter(f"Invalid TOML in {path}: {exc}") from exc

    if not isinstance(raw, dict):
        raise click.BadParameter(f"TOML config must be a table at the top level: {path}")

    return _normalize_config(raw)


def _normalize_config(raw: Mapping[str, Any]) -> dict[str, Any]:
    config: dict[str, Any] = {}
    for key, value in raw.items():
        normalized_key = key.replace("-", "_")
        config[normalized_key] = _normalize_value(normalized_key, value)
    return config


def _normalize_value(key: str, value: Any) -> Any:
    if key == "workflow" and isinstance(value, Sequence) and not isinstance(value, str):
        return ",".join(str(item).strip() for item in value if str(item).strip())

    if key in _LIST_KEYS:
        if value is None:
            return ()
        if isinstance(value, str):
            return (value,)
        if isinstance(value, Sequence):
            return tuple(str(item) for item in value)
        return (str(value),)

    if key in _PATH_KEYS and value is not None:
        return Path(str(value)).expanduser()

    if key in _BOOL_KEYS and not isinstance(value, bool):
        raise click.BadParameter(f"Config key '{key}' must be a boolean")

    if key in _INT_KEYS and value is not None and not isinstance(value, bool):
        try:
            return int(value)
        except (TypeError, ValueError) as exc:
            raise click.BadParameter(f"Config key '{key}' must be an integer") from exc

    return value


def merge_workflow_params(
    ctx: click.Context,
    toml_path: Optional[Path],
    params: MutableMapping[str, Any],
) -> dict[str, Any]:
    """Merge CLI parameters with optional TOML config (CLI wins when explicit)."""
    if toml_path is None:
        _validate_required_params(params)
        return dict(params)

    config = load_workflow_toml(toml_path)
    merged = dict(params)

    for key, config_value in config.items():
        if key == "toml":
            continue
        if _parameter_from_commandline(ctx, key):
            continue
        merged[key] = config_value

    _validate_required_params(merged)
    return merged


def load_minknow_from_workflow_toml(path: Path) -> Optional[MinKnowWorkflowConfig]:
    """Load ``[minknow]`` settings from a workflow TOML file."""
    try:
        config = load_workflow_toml(path)
    except click.BadParameter:
        return None
    minknow_raw = config.get("minknow")
    if not isinstance(minknow_raw, Mapping):
        return None
    return load_minknow_toml(
        path,
        workflow_config=config,
        prefer_workflow=True,
    )


def _parameter_from_commandline(ctx: click.Context, param_name: str) -> bool:
    source = ctx.get_parameter_source(param_name)
    return source in {ParameterSource.COMMANDLINE, ParameterSource.ENVIRONMENT}


def _validate_required_params(params: Mapping[str, Any]) -> None:
    missing = [
        key
        for key in WORKFLOW_REQUIRED_KEYS
        if params.get(key) in (None, "")
    ]
    if missing:
        readable = ", ".join(
            key.replace("_", "-") if key != "path" else "PATH"
            for key in missing
        )
        raise click.BadParameter(
            f"Missing required workflow setting(s): {readable}. "
            "Provide them on the command line or in the TOML file."
        )
