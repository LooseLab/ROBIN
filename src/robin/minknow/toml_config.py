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


@dataclass(frozen=True)
class MinKnowWorkflowConfig:
    """MinKNOW host settings and optional run preset from workflow TOML."""

    settings: MinKnowSettings
    preset: Optional[RobinRunPreset] = None


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
        raise click.BadParameter(f"TOML config must be a table at the top level: {path}")

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

    if config.preset is None:
        return config

    preset = apply_workflow_refs_to_preset(
        config.preset,
        workflow_config=merged_workflow or None,
        reference=reference,
        target_panel=target_panel,
        prefer_workflow=prefer_workflow,
    )
    return MinKnowWorkflowConfig(settings=config.settings, preset=preset)


def extract_minknow_config(raw: Mapping[str, Any]) -> MinKnowWorkflowConfig:
    """Parse ``[minknow]`` and nested ``[minknow.preset]`` from workflow TOML."""
    minknow_raw = raw.get("minknow")
    if not isinstance(minknow_raw, Mapping):
        raise click.BadParameter("[minknow] must be a TOML table")

    settings = MinKnowSettings.from_mapping(minknow_raw)
    preset = optional_preset_from_mapping(minknow_raw)
    return MinKnowWorkflowConfig(settings=settings, preset=preset)


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
    }
)
