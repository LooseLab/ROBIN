"""Load MinKNOW settings and presets from TOML."""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional

import click

from robin.minknow.config import MinKnowSettings
from robin.minknow.preset import RobinRunPreset
from robin.minknow.workflow_refs import (
    apply_workflow_refs_to_preset,
    extract_workflow_ref_keys,
)
from robin.readfish.config import ReadfishConfig


@dataclass(frozen=True)
class MinKnowWorkflowConfig:
    """MinKNOW host settings and optional run preset from workflow TOML."""

    settings: MinKnowSettings
    preset: Optional[RobinRunPreset] = None
    readfish: Optional[ReadfishConfig] = None


def load_minknow_toml(
    path: Path,
    *,
    workflow_config: Optional[Mapping[str, Any]] = None,
    reference: Optional[str] = None,
    target_panel: Optional[str] = None,
    prefer_workflow: bool = True,
) -> MinKnowWorkflowConfig:
    """Load ``[minknow]`` / ``[minknow.preset]`` or a standalone preset file."""
    if not path.exists():
        raise click.BadParameter(f"MinKNOW config file does not exist: {path}")
    if not path.is_file():
        raise click.BadParameter(f"MinKNOW config path is not a file: {path}")

    try:
        with path.open("rb") as handle:
            raw = tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise click.BadParameter(f"Invalid TOML in {path}: {exc}") from exc

    if not isinstance(raw, dict):
        raise click.BadParameter(
            f"TOML config must be a table at the top level: {path}"
        )

    inline_workflow = extract_workflow_ref_keys(raw)
    merged_workflow: dict[str, Any] = dict(inline_workflow)
    if workflow_config:
        merged_workflow.update(workflow_config)

    if "minknow" in raw:
        config = extract_minknow_config(raw)
    else:
        preset = RobinRunPreset.from_mapping(raw)
        host = str(raw.get("host") or "localhost").strip()
        settings = MinKnowSettings.from_host(host)
        config = MinKnowWorkflowConfig(settings=settings, preset=preset)

    readfish = ReadfishConfig.from_mapping(raw.get("readfish"))
    if readfish is not None:
        config = MinKnowWorkflowConfig(
            settings=config.settings,
            preset=config.preset,
            readfish=readfish,
        )

    if config.preset is None:
        return config

    preset = apply_workflow_refs_to_preset(
        config.preset,
        workflow_config=merged_workflow or None,
        reference=reference,
        target_panel=target_panel,
        prefer_workflow=prefer_workflow,
    )
    return MinKnowWorkflowConfig(
        settings=config.settings,
        preset=preset,
        readfish=config.readfish,
    )


def extract_minknow_config(raw: Mapping[str, Any]) -> MinKnowWorkflowConfig:
    """Parse ``[minknow]`` and nested ``[minknow.preset]`` from workflow TOML."""
    minknow_raw = raw.get("minknow")
    if not isinstance(minknow_raw, Mapping):
        raise click.BadParameter("[minknow] must be a TOML table")

    settings = MinKnowSettings.from_mapping(minknow_raw)
    preset = optional_preset_from_mapping(minknow_raw)
    return MinKnowWorkflowConfig(settings=settings, preset=preset)


def resolve_minknow_gui_config(
    workflow_toml: Optional[Path] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
) -> Optional[MinKnowWorkflowConfig]:
    """Return MinKNOW GUI settings when explicitly configured via TOML or env.

    The GUI sequencer page is shown only when this returns a config with
    ``settings.enabled`` true. Sources (in order):

    - ``[minknow]`` in the workflow TOML passed to ``robin workflow --toml``
    - ``[minknow]`` in ``ROBIN_WORKFLOW_TOML``
    - ``MINKNOW_PRESET`` pointing at a preset / workflow TOML file
    - ``MINKNOW_ENABLED=true`` or an explicit ``MINKNOW_HOST`` environment variable
    """
    import os

    from robin.minknow.config import (
        MinKnowSettings,
        _env_bool,
        preset_path_from_environ,
        workflow_toml_from_environ,
    )
    from robin.workflow_config import load_minknow_from_workflow_toml

    env = dict(environ or os.environ)
    workflow_paths: list[Path] = []

    if workflow_toml is not None:
        workflow_paths.append(workflow_toml.expanduser())
    env_workflow = workflow_toml_from_environ(env)
    if env_workflow is not None:
        candidate = env_workflow.expanduser()
        if candidate not in workflow_paths:
            workflow_paths.append(candidate)

    for path in workflow_paths:
        if not path.is_file():
            continue
        config = load_minknow_from_workflow_toml(path)
        if config is not None and config.settings.enabled:
            return config

    preset_path = preset_path_from_environ(env)
    if preset_path is not None and preset_path.is_file():
        try:
            config = load_minknow_toml(preset_path)
        except click.BadParameter:
            config = None
        if config is not None and config.settings.enabled:
            return config

    if _env_bool(env.get("MINKNOW_ENABLED"), default=False):
        settings = MinKnowSettings.from_environ(env)
        if settings.enabled:
            return MinKnowWorkflowConfig(settings=settings)

    if "MINKNOW_HOST" in env and str(env.get("MINKNOW_HOST", "")).strip():
        settings = MinKnowSettings.from_environ(env)
        if settings.enabled:
            return MinKnowWorkflowConfig(settings=settings)

    return None


def minknow_gui_available() -> bool:
    """Return whether the running GUI can expose the MinKNOW page."""
    try:
        from robin.gui.app import get_gui_launcher

        launcher = get_gui_launcher()
    except Exception:
        return False
    return launcher is not None and launcher.minknow_gui_enabled


def minknow_gui_accessible() -> bool:
    """Return whether MinKNOW is available and allowed for the signed-in user."""
    try:
        from robin.gui.app import get_gui_launcher

        launcher = get_gui_launcher()
    except Exception:
        return False
    return launcher is not None and launcher.minknow_gui_accessible()


def optional_preset_from_mapping(data: Mapping[str, Any]) -> Optional[RobinRunPreset]:
    """Return a preset when ``[minknow.preset]`` or preset keys are present."""
    nested = data.get("preset")
    if isinstance(nested, Mapping):
        return RobinRunPreset.from_mapping({"preset": nested})
    if any(key in data for key in _PRESET_KEYS):
        return RobinRunPreset.from_mapping(data)
    return None


_PRESET_KEYS = frozenset(
    {
        "kit",
        "basecall_simplex_model",
        "alignment_reference",
        "reference",
        "bed_file",
        "simulation_bulk_file",
        "simulation",
    }
)
