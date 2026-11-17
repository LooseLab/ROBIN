"""BATMAN mode — readfish with live dynamic target updates."""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from robin.minknow.preset import RobinRunPreset
from robin.readfish.config import ReadfishConfig

if TYPE_CHECKING:
    from robin.minknow.toml_config import MinKnowWorkflowConfig


def is_batman_mode(
    *,
    preset: Optional[RobinRunPreset] = None,
    readfish: Optional[ReadfishConfig] = None,
) -> bool:
    """True when readfish adaptive sampling with live TOML updates is enabled."""
    if preset is None or not preset.readfish_adaptive_sampling_enabled():
        return False
    readfish_cfg = readfish if readfish is not None else ReadfishConfig()
    return readfish_cfg.live_updates_enabled


def is_batman_workflow_config(config: Optional[MinKnowWorkflowConfig]) -> bool:
    """Return whether workflow MinKNOW/readfish settings enable BATMAN mode."""
    if config is None:
        return False
    return is_batman_mode(preset=config.preset, readfish=config.readfish)
