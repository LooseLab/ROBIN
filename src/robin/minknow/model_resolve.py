"""Pick MinKNOW basecall simplex models available on the connected host."""

from __future__ import annotations

from typing import Any, Optional

from robin.minknow.preset import RobinRunPreset


def _simplex_version(name: str) -> Optional[str]:
    if "@" not in name:
        return None
    return name.rsplit("@", 1)[-1]


def score_simplex_model(
    name: str,
    *,
    requested: Optional[str] = None,
    requested_version: Optional[str] = None,
) -> int:
    """Higher scores are better for ROBIN (CpG methylation + HAC preferred)."""
    lower = name.lower()
    score = 0

    if "modbases" in lower:
        score += 100
    if "5hmc" in lower and "5mc" in lower:
        score += 50
    if "_cg_" in lower or "5mc_cg" in lower or lower.endswith("_cg"):
        score += 40
    elif "allcontext" in lower:
        score -= 200

    if "hac" in lower:
        score += 30
    elif "sup" in lower:
        score += 20
    elif "fast" in lower:
        score += 10

    version = requested_version or (
        _simplex_version(requested) if requested else None
    )
    if version and name.endswith(f"@{version}"):
        score += 5

    if requested:
        requested_lower = requested.lower()
        for token in ("r10.4.1", "r9.4.1", "400bps", "5khz"):
            if token in requested_lower and token in lower:
                score += 3

    return score


def pick_simplex_model(
    available: list[str],
    requested: str,
) -> tuple[Optional[str], list[str]]:
    """Return the best simplex model name on the host and any warnings."""
    if not available:
        return None, []

    if requested in available:
        return requested, []

    requested_version = _simplex_version(requested)
    resolved = max(
        available,
        key=lambda name: score_simplex_model(
            name,
            requested=requested,
            requested_version=requested_version,
        ),
    )

    warnings = [
        f"Simplex model {requested!r} not available; using {resolved!r} instead."
    ]
    if "modbases" in requested.lower() and "modbases" not in resolved.lower():
        warnings.append(
            "The sequencer does not have a CpG modified-base (5mC/5hmC) model "
            "installed. MGMT, Sturgeon, and other methylation analyses need "
            "modbases models on the MinKNOW host (e.g. via Dorado model updates)."
        )
    return resolved, warnings


def available_simplex_models(configs: Any) -> list[str]:
    names: list[str] = []
    for config in configs:
        for simplex in getattr(config, "simplex_models", ()):
            names.append(simplex.name)
    return sorted(set(names))


def resolve_preset_simplex_model(
    manager: Any,
    preset: RobinRunPreset,
    *,
    product_code: str,
    sample_rate: int,
) -> tuple[RobinRunPreset, list[str], list[str]]:
    """Resolve simplex model against MinKNOW; return (preset, warnings, errors)."""
    from robin.minknow._deps import require_minknow_api

    require_minknow_api()
    from minknow_api.tools import protocols

    warnings: list[str] = []
    try:
        configs = manager.find_basecall_configurations(
            product_code, preset.kit, sample_rate
        )
    except Exception as exc:
        return preset, warnings, [f"Could not query basecall configurations: {exc}"]

    if not configs:
        return (
            preset,
            warnings,
            [
                f"No basecall configurations for kit {preset.kit!r} "
                f"and product code {product_code!r}"
            ],
        )

    available = available_simplex_models(configs)
    try:
        protocols.find_simplex_model(configs, preset.basecall_simplex_model)
        resolved_name = preset.basecall_simplex_model
    except RuntimeError:
        resolved_name, warnings = pick_simplex_model(
            available, preset.basecall_simplex_model
        )
        if resolved_name is None:
            return (
                preset,
                warnings,
                [
                    f"Simplex model {preset.basecall_simplex_model!r} not available "
                    f"and no alternative models were found."
                ],
            )

    updated = preset
    if resolved_name != preset.basecall_simplex_model:
        updated = preset.with_overrides(basecall_simplex_model=resolved_name)

    errors: list[str] = []
    if updated.modified_models:
        try:
            simplex = protocols.find_simplex_model(configs, updated.basecall_simplex_model)
        except RuntimeError:
            simplex = None
        if simplex is not None:
            available_modified = {model.name for model in simplex.modified_models}
            for model_name in updated.modified_models:
                if model_name not in available_modified:
                    errors.append(
                        f"Modified model {model_name!r} not available for "
                        f"{updated.basecall_simplex_model!r}"
                    )

    return updated, warnings, errors
