"""Global GUI plotting preferences (admin-controlled, persisted)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

PLOTTING_PREFERENCES_KEY = "plotting_preferences"
PLOTTING_PREFERENCES_SCHEMA_VERSION = 1

CNV_REPORT_SCALE_PLOIDY = "ploidy"
CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE = "normalized_difference"
CNV_REPORT_SCALES = (
    CNV_REPORT_SCALE_PLOIDY,
    CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
)
CNV_REPORT_SCALE_LABELS = {
    CNV_REPORT_SCALE_PLOIDY: "Estimated ploidy (absolute copy number)",
    CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE: (
        "Log2 ratio (ploidy / expected copy number)"
    ),
}
DEFAULT_CNV_REPORT_SCALE = CNV_REPORT_SCALE_PLOIDY
_LEGACY_CNV_SCALE_ALIASES = {
    "log2_ratio": CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
}


@dataclass
class PlottingPreferencesConfig:
    """Admin-controlled defaults for plots in PDF reports."""

    schema_version: int = PLOTTING_PREFERENCES_SCHEMA_VERSION
    cnv_report_scale: str = DEFAULT_CNV_REPORT_SCALE
    updated_at: Optional[str] = None
    updated_by: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {
            "schema_version": self.schema_version,
            "cnv_report_scale": self.cnv_report_scale,
        }
        if self.updated_at:
            out["updated_at"] = self.updated_at
        if self.updated_by:
            out["updated_by"] = self.updated_by
        return out

    @classmethod
    def from_dict(cls, data: Optional[Dict[str, Any]]) -> "PlottingPreferencesConfig":
        if not data:
            return cls()
        scale = str(data.get("cnv_report_scale") or DEFAULT_CNV_REPORT_SCALE)
        scale = _LEGACY_CNV_SCALE_ALIASES.get(scale, scale)
        if scale not in CNV_REPORT_SCALES:
            scale = DEFAULT_CNV_REPORT_SCALE
        return cls(
            schema_version=int(data.get("schema_version") or 1),
            cnv_report_scale=scale,
            updated_at=data.get("updated_at"),
            updated_by=data.get("updated_by"),
        )

    def with_updates(
        self,
        *,
        cnv_report_scale: Optional[str] = None,
        updated_at: Optional[str] = None,
        updated_by: Optional[str] = None,
    ) -> "PlottingPreferencesConfig":
        scale = cnv_report_scale or self.cnv_report_scale
        scale = _LEGACY_CNV_SCALE_ALIASES.get(scale, scale)
        if scale not in CNV_REPORT_SCALES:
            scale = DEFAULT_CNV_REPORT_SCALE
        return PlottingPreferencesConfig(
            schema_version=PLOTTING_PREFERENCES_SCHEMA_VERSION,
            cnv_report_scale=scale,
            updated_at=updated_at or self.updated_at,
            updated_by=updated_by or self.updated_by,
        )


def resolve_cnv_report_scale(scale: Optional[str]) -> str:
    if not scale:
        return DEFAULT_CNV_REPORT_SCALE
    scale = _LEGACY_CNV_SCALE_ALIASES.get(scale, scale)
    if scale in CNV_REPORT_SCALES:
        return scale
    return DEFAULT_CNV_REPORT_SCALE


def cnv_summary_normalized_from_scale(scale: Optional[str]) -> bool:
    """Return True when the configured CNV report scale is log2 ratio mode."""
    return (
        resolve_cnv_report_scale(scale) == CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE
    )


def load_plotting_preferences(store=None) -> PlottingPreferencesConfig:
    """Load persisted admin plotting preferences from the security store."""
    if store is None:
        from robin.security import SecurityStore

        store = SecurityStore()
    raw = store.get_gui_setting(PLOTTING_PREFERENCES_KEY)
    return PlottingPreferencesConfig.from_dict(raw)


def resolve_cnv_summary_normalized(
    explicit: Optional[bool] = None,
    *,
    plotting_preferences: Optional[PlottingPreferencesConfig] = None,
) -> bool:
    """Resolve whether CNV summary/per-chromosome plots use log2 ratio mode.

    An explicit True/False (e.g. CLI ``--cnv-normalized-difference``) overrides
    admin plotting preferences. When ``explicit`` is None, the admin default is used.
    """
    if explicit is not None:
        return bool(explicit)
    prefs = plotting_preferences or load_plotting_preferences()
    return cnv_summary_normalized_from_scale(prefs.cnv_report_scale)


def cnv_report_ylabel(scale: str) -> str:
    if resolve_cnv_report_scale(scale) == CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE:
        return "Log2 ratio (ploidy / expected)"
    return "Estimated copy number / ploidy"


def cnv_report_genome_ylabel(scale: str) -> str:
    if resolve_cnv_report_scale(scale) == CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE:
        return "Log2 ratio (ploidy / expected)"
    return "Estimated ploidy"


def cnv_report_genome_ylabel_mathtext(scale: str) -> str:
    """Matplotlib axis label using the same sans-serif family as other report text."""
    return cnv_report_genome_ylabel(scale)


def cnv_report_plot_caption(scale: str) -> str:
    if resolve_cnv_report_scale(scale) == CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE:
        return (
            "Copy number variation across chromosomes "
            "(log2 ratio of observed ploidy to expected copy number; 0 = normal; "
            "red = gain, blue = loss, grey = within ±0.3). "
            "Panel target points (right axis) mark significantly amplified or lost "
            "genes at their sequencing coverage depth, labelled by gene name; "
            "dashed line = mean coverage."
        )
    return (
        "Copy number variation across chromosomes. "
        "Panel target points (right axis) mark significantly altered genes at "
        "sequencing coverage depth, labelled by gene name; dashed line = mean coverage."
    )
