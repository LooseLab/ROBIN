"""Global GUI plotting preferences (admin-controlled, persisted)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

from robin.reference_contigs import (
    DEFAULT_REFERENCE_CONTIG_SCOPE,
    REFERENCE_CONTIG_SCOPE_LABELS,
    REFERENCE_CONTIG_SCOPES,
    resolve_reference_contig_scope,
)

PLOTTING_PREFERENCES_KEY = "plotting_preferences"
PLOTTING_PREFERENCES_SCHEMA_VERSION = 4

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

# Live CNV GUI control defaults (Administration → Plotting).
CNV_GUI_GENE_COVERAGE_FILTER_ALL = "all"
CNV_GUI_GENE_COVERAGE_FILTER_OUTLIERS = "outliers"
CNV_GUI_GENE_COVERAGE_FILTERS = (
    CNV_GUI_GENE_COVERAGE_FILTER_ALL,
    CNV_GUI_GENE_COVERAGE_FILTER_OUTLIERS,
)
DEFAULT_CNV_GUI_GENE_COVERAGE_FILTER = CNV_GUI_GENE_COVERAGE_FILTER_OUTLIERS

CNV_GUI_COLOR_MODE_CHROMOSOME = "chromosome"
CNV_GUI_COLOR_MODE_VALUE = "value"
CNV_GUI_COLOR_MODES = (
    CNV_GUI_COLOR_MODE_CHROMOSOME,
    CNV_GUI_COLOR_MODE_VALUE,
)
DEFAULT_CNV_GUI_COLOR_MODE = CNV_GUI_COLOR_MODE_CHROMOSOME

DEFAULT_CNV_GUI_SHOW_BREAKPOINTS = True

# Gene name labels on CNV coverage markers (GUI ECharts + report PDF).
DEFAULT_CNV_GENE_LABEL_FONT_SIZE = 12
MIN_CNV_GENE_LABEL_FONT_SIZE = 6
MAX_CNV_GENE_LABEL_FONT_SIZE = 24


def _resolve_gene_label_font_size(value: Any) -> int:
    try:
        size = int(round(float(value)))
    except (TypeError, ValueError):
        return DEFAULT_CNV_GENE_LABEL_FONT_SIZE
    return max(MIN_CNV_GENE_LABEL_FONT_SIZE, min(MAX_CNV_GENE_LABEL_FONT_SIZE, size))


def _resolve_gene_coverage_filter(value: Any) -> str:
    vlow = str(value or "").strip().lower()
    if vlow in (
        CNV_GUI_GENE_COVERAGE_FILTER_ALL,
        "all genes",
    ):
        return CNV_GUI_GENE_COVERAGE_FILTER_ALL
    if vlow in (
        CNV_GUI_GENE_COVERAGE_FILTER_OUTLIERS,
        "outliers only",
        "outliers",
        "≠ average",
        "!= average",
    ):
        return CNV_GUI_GENE_COVERAGE_FILTER_OUTLIERS
    return DEFAULT_CNV_GUI_GENE_COVERAGE_FILTER


def _resolve_color_mode(value: Any) -> str:
    vlow = str(value or "").strip().lower()
    if vlow in (
        CNV_GUI_COLOR_MODE_VALUE,
        "up/down",
        "updown",
        "up_down",
        "up down",
    ):
        return CNV_GUI_COLOR_MODE_VALUE
    if vlow in (CNV_GUI_COLOR_MODE_CHROMOSOME, "chromosomes"):
        return CNV_GUI_COLOR_MODE_CHROMOSOME
    return DEFAULT_CNV_GUI_COLOR_MODE


@dataclass
class PlottingPreferencesConfig:
    """Admin-controlled defaults for plots in the GUI and PDF reports."""

    schema_version: int = PLOTTING_PREFERENCES_SCHEMA_VERSION
    cnv_report_scale: str = DEFAULT_CNV_REPORT_SCALE
    reference_contig_scope: str = DEFAULT_REFERENCE_CONTIG_SCOPE
    cnv_gui_gene_coverage_filter: str = DEFAULT_CNV_GUI_GENE_COVERAGE_FILTER
    cnv_gui_color_mode: str = DEFAULT_CNV_GUI_COLOR_MODE
    cnv_gui_show_breakpoints: bool = DEFAULT_CNV_GUI_SHOW_BREAKPOINTS
    cnv_gene_label_font_size: int = DEFAULT_CNV_GENE_LABEL_FONT_SIZE
    updated_at: Optional[str] = None
    updated_by: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {
            "schema_version": self.schema_version,
            "cnv_report_scale": self.cnv_report_scale,
            "reference_contig_scope": self.reference_contig_scope,
            "cnv_gui_gene_coverage_filter": self.cnv_gui_gene_coverage_filter,
            "cnv_gui_color_mode": self.cnv_gui_color_mode,
            "cnv_gui_show_breakpoints": bool(self.cnv_gui_show_breakpoints),
            "cnv_gene_label_font_size": int(self.cnv_gene_label_font_size),
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
        contig_scope = resolve_reference_contig_scope(
            data.get("reference_contig_scope")
        )
        gene_filter = _resolve_gene_coverage_filter(
            data.get(
                "cnv_gui_gene_coverage_filter",
                DEFAULT_CNV_GUI_GENE_COVERAGE_FILTER,
            )
        )
        color_mode = _resolve_color_mode(
            data.get("cnv_gui_color_mode", DEFAULT_CNV_GUI_COLOR_MODE)
        )
        show_bp_raw = data.get(
            "cnv_gui_show_breakpoints", DEFAULT_CNV_GUI_SHOW_BREAKPOINTS
        )
        if isinstance(show_bp_raw, bool):
            show_bp = show_bp_raw
        else:
            show_bp = str(show_bp_raw).strip().lower() not in (
                "0",
                "false",
                "no",
                "off",
                "hide",
            )
        gene_font = _resolve_gene_label_font_size(
            data.get("cnv_gene_label_font_size", DEFAULT_CNV_GENE_LABEL_FONT_SIZE)
        )
        schema_version = int(data.get("schema_version") or 1)
        if schema_version < PLOTTING_PREFERENCES_SCHEMA_VERSION:
            schema_version = PLOTTING_PREFERENCES_SCHEMA_VERSION
        return cls(
            schema_version=schema_version,
            cnv_report_scale=scale,
            reference_contig_scope=contig_scope,
            cnv_gui_gene_coverage_filter=gene_filter,
            cnv_gui_color_mode=color_mode,
            cnv_gui_show_breakpoints=show_bp,
            cnv_gene_label_font_size=gene_font,
            updated_at=data.get("updated_at"),
            updated_by=data.get("updated_by"),
        )

    def with_updates(
        self,
        *,
        cnv_report_scale: Optional[str] = None,
        reference_contig_scope: Optional[str] = None,
        cnv_gui_gene_coverage_filter: Optional[str] = None,
        cnv_gui_color_mode: Optional[str] = None,
        cnv_gui_show_breakpoints: Optional[bool] = None,
        cnv_gene_label_font_size: Optional[int] = None,
        updated_at: Optional[str] = None,
        updated_by: Optional[str] = None,
    ) -> "PlottingPreferencesConfig":
        scale = cnv_report_scale or self.cnv_report_scale
        scale = _LEGACY_CNV_SCALE_ALIASES.get(scale, scale)
        if scale not in CNV_REPORT_SCALES:
            scale = DEFAULT_CNV_REPORT_SCALE
        contig_scope = resolve_reference_contig_scope(
            reference_contig_scope or self.reference_contig_scope
        )
        gene_filter = _resolve_gene_coverage_filter(
            cnv_gui_gene_coverage_filter
            if cnv_gui_gene_coverage_filter is not None
            else self.cnv_gui_gene_coverage_filter
        )
        color_mode = _resolve_color_mode(
            cnv_gui_color_mode
            if cnv_gui_color_mode is not None
            else self.cnv_gui_color_mode
        )
        show_bp = (
            self.cnv_gui_show_breakpoints
            if cnv_gui_show_breakpoints is None
            else bool(cnv_gui_show_breakpoints)
        )
        gene_font = _resolve_gene_label_font_size(
            self.cnv_gene_label_font_size
            if cnv_gene_label_font_size is None
            else cnv_gene_label_font_size
        )
        return PlottingPreferencesConfig(
            schema_version=PLOTTING_PREFERENCES_SCHEMA_VERSION,
            cnv_report_scale=scale,
            reference_contig_scope=contig_scope,
            cnv_gui_gene_coverage_filter=gene_filter,
            cnv_gui_color_mode=color_mode,
            cnv_gui_show_breakpoints=show_bp,
            cnv_gene_label_font_size=gene_font,
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
    return resolve_cnv_report_scale(scale) == CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE


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


def resolve_plotting_reference_contig_scope(
    plotting_preferences: Optional[PlottingPreferencesConfig] = None,
) -> str:
    """Resolve the contig scope used for coverage and CNV plots."""
    prefs = plotting_preferences or load_plotting_preferences()
    return resolve_reference_contig_scope(prefs.reference_contig_scope)


def resolve_cnv_gui_gene_coverage_filter(
    plotting_preferences: Optional[PlottingPreferencesConfig] = None,
) -> str:
    """Default Coverage genes filter for the live CNV GUI."""
    prefs = plotting_preferences or load_plotting_preferences()
    return _resolve_gene_coverage_filter(prefs.cnv_gui_gene_coverage_filter)


def resolve_cnv_gui_color_mode(
    plotting_preferences: Optional[PlottingPreferencesConfig] = None,
) -> str:
    """Default Color-by mode for the live CNV GUI."""
    prefs = plotting_preferences or load_plotting_preferences()
    return _resolve_color_mode(prefs.cnv_gui_color_mode)


def resolve_cnv_gui_show_breakpoints(
    plotting_preferences: Optional[PlottingPreferencesConfig] = None,
) -> bool:
    """Default Breakpoints visibility for single-chromosome CNV GUI view."""
    prefs = plotting_preferences or load_plotting_preferences()
    return bool(prefs.cnv_gui_show_breakpoints)


def resolve_cnv_gui_y_scale(
    plotting_preferences: Optional[PlottingPreferencesConfig] = None,
) -> str:
    """Default Y-axis mode for the live CNV GUI (``linear`` or ``log``)."""
    if resolve_cnv_summary_normalized(None, plotting_preferences=plotting_preferences):
        return "log"
    return "linear"


def resolve_cnv_gene_label_font_size(
    plotting_preferences: Optional[PlottingPreferencesConfig] = None,
) -> int:
    """Font size for gene name labels on CNV coverage markers (GUI + report)."""
    prefs = plotting_preferences or load_plotting_preferences()
    return _resolve_gene_label_font_size(prefs.cnv_gene_label_font_size)


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
