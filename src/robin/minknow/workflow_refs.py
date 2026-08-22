"""Resolve MinKNOW reference and panel paths from the active ROBIN workflow."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping, Optional

from robin.minknow.preset import RobinRunPreset
from robin.utils.sequencing_files import resolve_panel_stranded_bed_path


def extract_workflow_ref_keys(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Pull ``reference`` / ``target_panel`` from a workflow TOML table."""
    refs: dict[str, Any] = {}
    reference = raw.get("reference")
    if reference not in (None, ""):
        refs["reference"] = reference
    target_panel = raw.get("target_panel")
    if target_panel not in (None, ""):
        refs["target_panel"] = str(target_panel).strip()
    return refs


def resolve_reference_path(
    *,
    explicit: Any = None,
    workflow_config: Optional[Mapping[str, Any]] = None,
    environ: Optional[Mapping[str, str]] = None,
) -> Optional[str]:
    """Resolve the alignment reference FASTA from workflow config or environment."""
    candidates: list[str] = []

    if explicit not in (None, ""):
        candidates.append(str(Path(explicit).expanduser()))

    if workflow_config:
        reference = workflow_config.get("reference")
        if reference not in (None, ""):
            candidates.append(str(Path(reference).expanduser()))

    env = dict(environ or os.environ)
    env_reference = env.get("robin_REFERENCE")
    if env_reference:
        candidates.append(str(Path(env_reference).expanduser()))

    seen: set[str] = set()
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        return candidate
    return None


def resolve_target_panel(
    *,
    explicit: Optional[str] = None,
    workflow_config: Optional[Mapping[str, Any]] = None,
) -> Optional[str]:
    """Resolve the ROBIN target panel name."""
    if explicit and str(explicit).strip():
        return str(explicit).strip()
    if workflow_config:
        panel = workflow_config.get("target_panel")
        if panel not in (None, ""):
            return str(panel).strip()
    return None


def resolve_panel_bed_file(
    target_panel: Optional[str],
    *,
    reference_path: Optional[str] = None,
) -> Optional[str]:
    """Return the stranded panel BED for MinKNOW / adaptive sampling."""
    if not target_panel:
        return None
    bed_path = resolve_panel_stranded_bed_path(
        target_panel,
        reference_path=reference_path,
    )
    return str(bed_path) if bed_path is not None else None


def workflow_context_from_runner(runner: Any) -> tuple[Optional[str], Optional[str]]:
    """Read reference and target panel from a live workflow runner, if available."""
    if runner is None:
        return None, None

    reference = None
    if hasattr(runner, "get_reference"):
        try:
            reference = runner.get_reference()
        except Exception:
            reference = None
    if reference is None:
        reference = getattr(runner, "reference", None)

    target_panel = None
    if hasattr(runner, "get_target_panel"):
        try:
            target_panel = runner.get_target_panel()
        except Exception:
            target_panel = None
    if target_panel is None:
        target_panel = getattr(runner, "target_panel", None)

    ref_text = (
        str(Path(reference).expanduser()) if reference not in (None, "") else None
    )
    panel_text = str(target_panel).strip() if target_panel not in (None, "") else None
    return ref_text, panel_text


def apply_workflow_refs_to_preset(
    preset: RobinRunPreset,
    *,
    workflow_config: Optional[Mapping[str, Any]] = None,
    reference: Optional[str] = None,
    target_panel: Optional[str] = None,
    prefer_workflow: bool = True,
    environ: Optional[Mapping[str, str]] = None,
) -> RobinRunPreset:
    """Fill preset reference/BED paths from the active ROBIN workflow when possible."""
    resolved_panel = resolve_target_panel(
        explicit=target_panel,
        workflow_config=workflow_config,
    )
    resolved_reference = resolve_reference_path(
        explicit=reference,
        workflow_config=workflow_config,
        environ=environ,
    )
    resolved_bed = resolve_panel_bed_file(
        resolved_panel,
        reference_path=resolved_reference,
    )

    overrides: dict[str, Any] = {}
    if resolved_reference and (prefer_workflow or not preset.alignment_reference):
        overrides["alignment_reference"] = resolved_reference
    if resolved_bed and (prefer_workflow or not preset.bed_file):
        overrides["bed_file"] = resolved_bed

    if not overrides:
        return preset
    return preset.with_overrides(**overrides)


def load_workflow_config_for_refs(path: Path) -> Optional[dict[str, Any]]:
    """Load workflow TOML keys needed for MinKNOW path resolution."""
    if not path.is_file():
        return None
    try:
        from robin.workflow_config import load_workflow_toml

        return load_workflow_toml(path)
    except Exception:
        return None
