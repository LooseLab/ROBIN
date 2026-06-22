"""Pick MinKNOW basecall simplex models available on the connected host."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from robin.minknow.preset import RobinRunPreset


def _simplex_version(name: str) -> Optional[str]:
    if "@" not in name:
        return None
    return name.rsplit("@", 1)[-1]


def _modified_model_base(name: str) -> str:
    return name.split("@", 1)[0]


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
    if _requests_cpg_methylation(requested, ()):
        if "modbases" in requested.lower() and "modbases" not in resolved.lower():
            warnings.append(
                "The sequencer does not have a CpG modified-base (5mC/5hmC) model "
                "installed. MGMT, Sturgeon, and other methylation analyses need "
                "modbases models on the MinKNOW host (e.g. via Dorado model updates)."
            )
    return resolved, warnings


def pick_modified_model(
    available: set[str],
    requested: str,
    *,
    simplex_model: Optional[str] = None,
) -> Optional[str]:
    """Match a requested modified model to an installed name (with/without @version).

    MinKNOW often exposes compound names such as
    ``dna_r10.4.1_e8.2_400bps_hac@v5.2.0_5mCG_5hmCG@v2`` rather than ``5mCG_5hmCG``.
    """
    if requested in available:
        return requested

    base = _modified_model_base(requested)
    by_base = sorted(m for m in available if _modified_model_base(m) == base)
    if by_base:
        return _prefer_simplex_prefixed(by_base, simplex_model)

    requested_lower = base.lower()
    version = _simplex_version(requested)
    suffix = f"_{base}" if base else ""
    version_suffix = f"{suffix}@{version}" if version and suffix else ""

    compound_matches = sorted(
        m
        for m in available
        if (
            (suffix and m.endswith(version_suffix))
            or (suffix and m.endswith(suffix))
            or (requested_lower and requested_lower in m.lower())
        )
    )
    if compound_matches:
        if "5mcg" in requested_lower and "5hmc" in requested_lower:
            cpg = [
                m
                for m in compound_matches
                if "5mcg" in m.lower() and "5hmc" in m.lower()
            ]
            if cpg:
                return _prefer_simplex_prefixed(cpg, simplex_model)
        return _prefer_simplex_prefixed(compound_matches, simplex_model)

    if "5mcg" in base.lower() or "5hmc" in base.lower():
        cpg_models = sorted(
            m
            for m in available
            if "5mcg" in m.lower() and "5hmc" in m.lower()
        )
        if cpg_models:
            return _prefer_simplex_prefixed(cpg_models, simplex_model)

    return None


def recommended_cpg_modified_model(
    simplex_model: str,
    available_modified: set[str],
) -> Optional[str]:
    """Return the CpG 5mC/5hmC modified model paired with ``simplex_model``."""
    return pick_modified_model(
        available_modified,
        "5mCG_5hmCG",
        simplex_model=simplex_model,
    )


def _prefer_simplex_prefixed(
    candidates: list[str],
    simplex_model: Optional[str],
) -> str:
    if simplex_model:
        prefixed = [name for name in candidates if name.startswith(simplex_model)]
        if prefixed:
            return prefixed[-1]
    return candidates[-1]


def pick_methylation_simplex(
    available: list[str],
    *,
    reference_simplex: str,
) -> Optional[str]:
    """Return the best installed simplex model with integrated CpG modbases."""
    candidates = [
        name
        for name in available
        if score_simplex_model(name, requested=reference_simplex) >= 140
    ]
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda name: score_simplex_model(name, requested=reference_simplex),
    )


def methylation_capable_simplex_models(available: list[str]) -> list[str]:
    """Simplex models that include CpG 5mC/5hmC calling (integrated modbases)."""
    return sorted(
        name for name in available if score_simplex_model(name) >= 140
    )


def _requests_cpg_methylation(
    simplex_model: str,
    modified_models: tuple[str, ...],
) -> bool:
    if modified_models:
        return True
    lower = simplex_model.lower()
    return "modbases" in lower or "5mc" in lower or "5hmc" in lower


@dataclass(frozen=True)
class SimplexModelInfo:
    name: str
    modified_models: tuple[str, ...]


def collect_simplex_model_info(configs: Any) -> list[SimplexModelInfo]:
    """Return simplex models and attachable modified models from MinKNOW configs."""
    models: list[SimplexModelInfo] = []
    seen: set[str] = set()
    for config in configs:
        for simplex in getattr(config, "simplex_models", ()):
            if simplex.name in seen:
                continue
            seen.add(simplex.name)
            modified = tuple(
                sorted(model.name for model in getattr(simplex, "modified_models", ()))
            )
            models.append(SimplexModelInfo(name=simplex.name, modified_models=modified))
    return sorted(models, key=lambda item: item.name)


def available_simplex_models(configs: Any) -> list[str]:
    return [model.name for model in collect_simplex_model_info(configs)]


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
        resolved_name, pick_warnings = pick_simplex_model(
            available, preset.basecall_simplex_model
        )
        warnings.extend(pick_warnings)
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
            simplex = protocols.find_simplex_model(
                configs, updated.basecall_simplex_model
            )
        except RuntimeError:
            simplex = None

        available_modified: set[str] = set()
        if simplex is not None:
            available_modified = {
                model.name for model in simplex.modified_models
            }

        resolved_modified: list[str] = []
        unresolved: list[str] = []
        for model_name in updated.modified_models:
            picked = pick_modified_model(
                available_modified,
                model_name,
                simplex_model=updated.basecall_simplex_model,
            )
            if picked is not None:
                resolved_modified.append(picked)
                if picked != model_name:
                    warnings.append(
                        f"Modified model {model_name!r} not available; "
                        f"using {picked!r} instead."
                    )
            else:
                unresolved.append(model_name)

        if unresolved:
            methylation_simplex = pick_methylation_simplex(
                available,
                reference_simplex=updated.basecall_simplex_model,
            )
            if methylation_simplex is not None:
                warnings.append(
                    f"Modified models {list(updated.modified_models)!r} are not "
                    f"available for {updated.basecall_simplex_model!r}; "
                    f"using integrated methylation simplex "
                    f"{methylation_simplex!r}."
                )
                updated = updated.with_overrides(
                    basecall_simplex_model=methylation_simplex,
                    modified_models=(),
                )
            else:
                integrated = methylation_capable_simplex_models(available)
                hint = (
                    f"Methylation-capable simplex models on this host: "
                    f"{integrated!r}."
                    if integrated
                    else "No methylation-capable simplex models are installed."
                )
                for model_name in unresolved:
                    errors.append(
                        f"Modified model {model_name!r} not available for "
                        f"{updated.basecall_simplex_model!r}. "
                        f"Attachable modified models: "
                        f"{sorted(available_modified) or 'none'}. "
                        f"{hint} "
                        f"Run `robin minknow models --host <host> --position <pos>` "
                        f"to list installed models."
                    )
        elif tuple(resolved_modified) != updated.modified_models:
            updated = updated.with_overrides(
                modified_models=tuple(resolved_modified)
            )

    return updated, warnings, errors


def query_basecall_models(
    manager: Any,
    *,
    product_code: str,
    kit: str,
    sample_rate: int,
) -> tuple[list[SimplexModelInfo], Optional[str]]:
    """Query MinKNOW for basecall models; return (models, error_message)."""
    try:
        configs = manager.find_basecall_configurations(
            product_code, kit, sample_rate
        )
    except Exception as exc:
        return [], f"Could not query basecall configurations: {exc}"

    if not configs:
        return [], (
            f"No basecall configurations for kit {kit!r} "
            f"and product code {product_code!r}"
        )

    return collect_simplex_model_info(configs), None
