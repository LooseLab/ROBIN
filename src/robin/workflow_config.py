"""Load and merge TOML configuration for ``robin workflow``."""

from __future__ import annotations

import logging
import tomllib
from pathlib import Path
from typing import Any, Mapping, MutableMapping, Optional, Sequence

import click
from click.core import ParameterSource

from robin.minknow.toml_config import MinKnowWorkflowConfig, extract_minknow_config, load_minknow_toml

WORKFLOW_REQUIRED_KEYS = ("path", "workflow", "center", "target_panel")

# Default ruptures KernelCPD penalty for CNV breakpoint detection.
DEFAULT_CNV_PENALTY_VALUE = 10
# Minimum adjacent same-sign bins for a gain/loss region (|CNV| > 0.5).
DEFAULT_CNV_MIN_CONTIGUOUS_BINS = 1

_logger = logging.getLogger("robin.workflow_config")

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

# Optional CNV settings. Defaults: penalty_value=10, min_contiguous_bins=1.
# [cnv]
# penalty_value = 10
# min_contiguous_bins = 3
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
    try:
        return load_minknow_toml(
            path,
            workflow_config=config,
            prefer_workflow=True,
        )
    except click.BadParameter:
        return None


def parse_cnv_penalty_value(
    value: Any,
    *,
    default: int = DEFAULT_CNV_PENALTY_VALUE,
) -> int:
    """Coerce a CNV ruptures penalty to a positive integer, or return ``default``."""
    return _parse_positive_int(
        value,
        default=default,
        setting_name="penalty_value",
    )


def parse_cnv_min_contiguous_bins(
    value: Any,
    *,
    default: int = DEFAULT_CNV_MIN_CONTIGUOUS_BINS,
) -> int:
    """Coerce min contiguous CNV bins to a positive integer, or return ``default``."""
    return _parse_positive_int(
        value,
        default=default,
        setting_name="min_contiguous_bins",
    )


def _parse_positive_int(
    value: Any,
    *,
    default: int,
    setting_name: str,
) -> int:
    if value is None:
        return int(default)
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        _logger.warning(
            "Invalid [cnv] %s %r; using default %s",
            setting_name,
            value,
            default,
        )
        return int(default)
    if parsed <= 0:
        _logger.warning(
            "Invalid [cnv] %s %s (must be > 0); using default %s",
            setting_name,
            parsed,
            default,
        )
        return int(default)
    return parsed


def _load_cnv_section(
    workflow_config: Optional[Mapping[str, Any]] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
) -> Optional[Mapping[str, Any]]:
    config = workflow_config
    if config is None:
        from robin.minknow.config import workflow_toml_from_environ

        toml_path = workflow_toml_from_environ(environ)
        if toml_path is not None and toml_path.is_file():
            try:
                config = load_workflow_toml(toml_path)
            except click.BadParameter as exc:
                _logger.warning(
                    "Could not load workflow TOML for [cnv] settings (%s): %s",
                    toml_path,
                    exc,
                )
                config = None

    if not isinstance(config, Mapping):
        return None
    cnv_section = config.get("cnv")
    if not isinstance(cnv_section, Mapping):
        return None
    return cnv_section


def get_cnv_penalty_value(
    workflow_config: Optional[Mapping[str, Any]] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
    default: int = DEFAULT_CNV_PENALTY_VALUE,
) -> int:
    """
    Resolve CNV breakpoint ``penalty_value`` from workflow config.

    Looks for::

        [cnv]
        penalty_value = 10

    If ``workflow_config`` is omitted, loads the TOML referenced by
    ``ROBIN_WORKFLOW_TOML`` when set. Returns ``default`` (10) when unset
    or invalid.
    """
    cnv_section = _load_cnv_section(workflow_config, environ=environ)
    if cnv_section is None:
        return int(default)
    return parse_cnv_penalty_value(cnv_section.get("penalty_value"), default=default)


def get_cnv_min_contiguous_bins(
    workflow_config: Optional[Mapping[str, Any]] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
    default: int = DEFAULT_CNV_MIN_CONTIGUOUS_BINS,
) -> int:
    """
    Resolve minimum contiguous significant bins for CNV gain/loss regions.

    Looks for::

        [cnv]
        min_contiguous_bins = 3

    Default is 1 (any single significant bin counts). Used when writing
    ``new_file_*.bed`` gain/loss regions from the normalized CNV track.
    """
    cnv_section = _load_cnv_section(workflow_config, environ=environ)
    if cnv_section is None:
        return int(default)
    return parse_cnv_min_contiguous_bins(
        cnv_section.get("min_contiguous_bins"),
        default=default,
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
