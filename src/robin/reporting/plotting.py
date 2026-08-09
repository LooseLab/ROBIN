"""
plotting.py

This module contains functions for creating plots used in the PDF report.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpath_effects
import io
import textwrap
import matplotlib.font_manager as fm
import os
from robin.gui import fonts

import natsort

from matplotlib import gridspec

import logging
from typing import List, Optional, Sequence, Tuple, Dict, Any

logger = logging.getLogger(__name__)

from robin.reference_contigs import is_visible_contig

# Define consistent color scheme and style
MODERN_COLORS = {
    "primary": "#2C3E50",  # Dark blue-grey (matching report text)
    "secondary": "#E2E8F0",  # Light grey (matching table grid)
    "background": "#F8FAFC",  # Light background (matching table alternate rows)
    "accent": "#3498DB",  # Blue accent
    "grid": "#E2E8F0",  # Grid color
}

CNV_COLORS = {
    "points": "#5B7C9D",
    "median": "#1E3A5F",
    "diploid": "#3D5A80",
    "reference": "#C5CDD6",
    "gain_fill": "#D4EDDA",
    "loss_fill": "#F8D7DA",
    "gain_edge": "#3D7A4A",
    "loss_edge": "#A94442",
    "gene": "#5C6BC0",
    "plot_gain": "#DC2626",
    "plot_loss": "#2563EB",
    "plot_neutral": "#9CA3AF",
}

# Scatter styling for dense CNV tracks.
CNV_POINT_MARKER = "."
CNV_POINT_ALPHA = 1.0
CNV_POINT_ALPHA_NEUTRAL = CNV_POINT_ALPHA
CNV_POINT_ALPHA_CALLED = CNV_POINT_ALPHA
CNV_POINT_ALPHA_DEFAULT = CNV_POINT_ALPHA

_CNV_GENE_LABEL_PATH_EFFECTS = [
    mpath_effects.withStroke(linewidth=2.6, foreground="white", alpha=0.95),
]
_CNV_GENE_LABEL_BBOX = {
    "boxstyle": "round,pad=0.18",
    "facecolor": "white",
    "edgecolor": "none",
    "alpha": 0.70,
}

CNV_FONT = {
    "title": 11,
    "subtitle": 8,
    "axis": 9.5,
    "tick": 7.5,
    "annotation": 6.5,
}

CNV_TEXT = {
    "primary": MODERN_COLORS["primary"],
    "muted": "#5A6B7D",
}

CNV_CHROMOSOME_FIG_WIDTH = 7.5
CNV_CHROMOSOME_SCATTER_SIZE = 6
CNV_CHROMOSOME_PLOTS_PER_PAGE = 4
CNV_CHROMOSOME_PLOT_SPACER_PT = 6
# ReportLab's default Frame uses 6pt padding on each side.
CNV_REPORT_FRAME_PADDING_PT = 12
# Genome-wide summary plot sits on a dedicated A4 landscape page.
CNV_GENOME_LANDSCAPE_FIG_WIDTH = 16.0
CNV_GENOME_LANDSCAPE_FIG_HEIGHT = 6.0
CNV_GENOME_LANDSCAPE_CAPTION_RESERVE_PT = 42.0

_CNV_FONT_REGULAR: Optional[fm.FontProperties] = None
_CNV_FONT_BOLD: Optional[fm.FontProperties] = None
_CNV_FONTS_READY = False


def _setup_cnv_fonts() -> None:
    """Register report plot fonts with a safe fallback when bundled TTFs are absent."""
    global _CNV_FONT_REGULAR, _CNV_FONT_BOLD, _CNV_FONTS_READY
    if _CNV_FONTS_READY:
        return

    fonts_dir = os.path.dirname(os.path.abspath(fonts.__file__))
    regular_path = os.path.join(fonts_dir, "fira-sans-v16-latin-regular.ttf")
    bold_path = os.path.join(fonts_dir, "fira-sans-v16-latin-700.ttf")
    family = "DejaVu Sans"

    if os.path.isfile(regular_path):
        try:
            fm.fontManager.addfont(regular_path)
            _CNV_FONT_REGULAR = fm.FontProperties(fname=regular_path)
            family = _CNV_FONT_REGULAR.get_name()
        except Exception:
            _CNV_FONT_REGULAR = fm.FontProperties(family=family)
    else:
        _CNV_FONT_REGULAR = fm.FontProperties(family=family)

    if os.path.isfile(bold_path):
        try:
            fm.fontManager.addfont(bold_path)
            _CNV_FONT_BOLD = fm.FontProperties(fname=bold_path)
        except Exception:
            _CNV_FONT_BOLD = fm.FontProperties(family=family, weight="bold")
    else:
        _CNV_FONT_BOLD = fm.FontProperties(family=family, weight="bold")

    plt.rcParams["font.family"] = family
    # Keep math labels on the same sans family as axis text (no Computer Modern mix).
    plt.rcParams["mathtext.fontset"] = "dejavusans"
    _CNV_FONTS_READY = True


def set_modern_style():
    """Set consistent modern style for all plots"""
    _setup_cnv_fonts()

    plt.style.use("seaborn-v0_8-whitegrid")
    sns.set_theme(style="whitegrid", font=plt.rcParams["font.family"])

    plt.rcParams.update(
        {
            # Figure settings
            "figure.facecolor": "white",
            "figure.dpi": 300,
            # Axes settings
            "axes.facecolor": "white",
            "axes.edgecolor": MODERN_COLORS["primary"],
            "axes.labelcolor": MODERN_COLORS["primary"],
            "axes.titlecolor": MODERN_COLORS["primary"],
            "axes.grid": True,
            "axes.labelsize": CNV_FONT["axis"],
            "axes.titlesize": CNV_FONT["title"],
            # Grid settings
            "grid.color": MODERN_COLORS["grid"],
            "grid.linestyle": "--",
            "grid.linewidth": 0.5,
            "grid.alpha": 0.5,
            # Tick settings
            "xtick.color": MODERN_COLORS["primary"],
            "ytick.color": MODERN_COLORS["primary"],
            "xtick.labelsize": CNV_FONT["tick"],
            "ytick.labelsize": CNV_FONT["tick"],
            # Legend settings
            "legend.frameon": True,
            "legend.facecolor": "white",
            "legend.edgecolor": MODERN_COLORS["grid"],
            "legend.fontsize": CNV_FONT["tick"],
            # Line settings
            "lines.linewidth": 1.5,
            "lines.markersize": 6,
        }
    )


def _apply_cnv_chromosome_axes(
    ax,
    x_max_mb: float,
    y_max: float,
    *,
    y_min: float = 0.0,
    xlabel: str,
    ylabel: str,
    show_xlabel: bool = True,
) -> None:
    """Pin axes at the origin with no leading x padding."""
    _setup_cnv_fonts()
    ax.set_xlim(0, x_max_mb)
    ax.set_ylim(y_min, y_max)
    ax.margins(x=0, y=0)
    ax.autoscale(enable=False)
    if show_xlabel:
        ax.set_xlabel(
            xlabel,
            fontsize=CNV_FONT["axis"],
            color=CNV_TEXT["primary"],
            labelpad=10,
            fontproperties=_CNV_FONT_REGULAR,
        )
    else:
        ax.set_xlabel("")
    ax.set_ylabel(
        ylabel,
        fontsize=CNV_FONT["axis"],
        color=CNV_TEXT["primary"],
        labelpad=10,
        fontproperties=_CNV_FONT_REGULAR,
    )
    ax.spines["left"].set_position(("data", 0))
    if y_min < 0:
        # Keep the Mb axis at the bottom of the panel; log2 reference stays at y=0 inside.
        ax.spines["bottom"].set_position(("axes", 0.0))
    else:
        ax.spines["bottom"].set_position(("data", 0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.grid(True, axis="y", color=MODERN_COLORS["grid"], linestyle="--", linewidth=0.4, alpha=0.55)
    ax.grid(False, axis="x")
    ax.tick_params(colors=CNV_TEXT["primary"], labelsize=CNV_FONT["tick"])


def _chromosome_display_name(contig: str) -> str:
    return contig.replace("chr", "")


def _plot_cnv_track(
    ax,
    cnv_df: pd.DataFrame,
    x_max_mb: float,
    *,
    point_size: float = 4,
) -> None:
    """Render CNV bin values as scatter with a rolling median trace."""
    ax.scatter(
        cnv_df["position_mb"],
        cnv_df["ploidy"],
        s=point_size,
        c=CNV_COLORS["points"],
        marker=CNV_POINT_MARKER,
        alpha=CNV_POINT_ALPHA_DEFAULT,
        linewidths=0,
        edgecolors="none",
        rasterized=True,
        zorder=2,
    )
    window = max(5, len(cnv_df) // 35)
    if len(cnv_df) >= window:
        median = cnv_df["ploidy"].rolling(window, center=True, min_periods=1).median()
        ax.plot(
            cnv_df["position_mb"],
            median,
            color=CNV_COLORS["median"],
            linewidth=1.4,
            alpha=0.95,
            solid_capstyle="round",
            zorder=4,
        )


PANEL_GENE_CNV_OUTLIER_SD = 3


def _is_gene_cnv_outlier(cnv_val: float, mean_cnv: float, std_cnv: float) -> bool:
    """Return True when a panel target copy number deviates >3 SD from the chromosome mean."""
    if std_cnv < 1e-6:
        return abs(cnv_val - mean_cnv) > 0.5
    return abs(cnv_val - mean_cnv) > PANEL_GENE_CNV_OUTLIER_SD * std_cnv


SIGNIFICANT_CNV_REGION_TYPES = {"GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"}


def _should_label_panel_gene(
    point: Dict[str, Any],
    mean_cnv: float,
    std_cnv: float,
) -> bool:
    """Label panel genes more than 3 SD from the chromosome mean."""
    return _is_gene_cnv_outlier(point["cnv_val"], mean_cnv, std_cnv)


def _is_highlighted_panel_gene(
    point: Dict[str, Any],
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
    off_scale_mode: bool,
    y_max: float,
    *,
    y_min: float = 0.0,
) -> tuple[bool, bool]:
    """Return (highlight, off_scale) for a panel target."""
    cnv_val = point["cnv_val"]
    if y_min < 0:
        off_scale = off_scale_mode and abs(cnv_val) > y_max * 0.97
    else:
        off_scale = off_scale_mode and cnv_val > y_max * 0.97
    highlight = off_scale or _should_label_panel_gene(point, mean_cnv, std_cnv)
    return highlight, off_scale


def _add_cnv_regions_on_plot(ax, regions: List[Dict[str, Any]], y_max: float) -> None:
    """Shade called regions on the main track with edge guides and top brackets."""
    for region in regions:
        start_mb = region["start_pos"] / 1_000_000
        end_mb = region["end_pos"] / 1_000_000
        is_gain = region["type"] in ("GAIN", "HIGH_GAIN")
        fill_color = CNV_COLORS["gain_fill"] if is_gain else CNV_COLORS["loss_fill"]
        edge_color = CNV_COLORS["gain_edge"] if is_gain else CNV_COLORS["loss_edge"]

        ax.axvspan(start_mb, end_mb, color=fill_color, alpha=0.55, zorder=0, linewidth=0)
        ax.axvline(start_mb, color=edge_color, linestyle="--", linewidth=0.9, alpha=0.75, zorder=1)
        ax.axvline(end_mb, color=edge_color, linestyle="--", linewidth=0.9, alpha=0.75, zorder=1)
        _draw_region_bracket(
            ax, start_mb, end_mb, y_max * 0.992, edge_color, height_frac=y_max * 0.028,
        )


PANEL_LABEL_MAX_CHARS = 10
LOLLIPOP_LABEL_FONT_SIZE = 6


def _truncate_panel_label(label: str, max_chars: int = PANEL_LABEL_MAX_CHARS) -> str:
    label = str(label).strip()
    if len(label) <= max_chars:
        return label
    return f"{label[: max_chars - 1]}…"


def _wrap_chromosome_status_text(status_text: str, *, width: int = 88) -> str:
    parts = [part.strip() for part in status_text.split(";") if part.strip()]
    if not parts:
        return status_text
    lines: List[str] = []
    current = ""
    for part in parts:
        candidate = part if not current else f"{current}; {part}"
        if len(candidate) <= width:
            current = candidate
            continue
        if current:
            lines.append(current)
        current = textwrap.fill(part, width=width)
    if current:
        lines.append(current)
    return "\n".join(lines[:3])


def _chromosome_figure_margins(
    *,
    has_status: bool,
    compact: bool = False,
) -> Dict[str, float]:
    if compact:
        return {
            "top": 0.70 if has_status else 0.80,
            "bottom": 0.20,
            "right": 0.97,
            "left": 0.13,
        }
    return {
        "top": 0.76 if has_status else 0.86,
        "bottom": 0.18,
        "right": 0.97,
        "left": 0.13,
    }


def _normalise_coverage_to_cnv_axis(
    coverage: float,
    *,
    mean_cov: float,
    scale_mean_cnv: float,
    use_log: bool,
) -> Optional[float]:
    """Map target depth onto the CNV y-axis for visual comparison.

    Linear/ploidy: ``scale_mean_cnv * (cov / mean_cov)`` so mean coverage sits
    at the CNV mean. Log2: ``log2(cov / mean_cov)`` so mean coverage sits at 0.
    """
    if not np.isfinite(coverage) or coverage <= 0:
        return None
    if not np.isfinite(mean_cov) or mean_cov <= 0:
        return None
    ratio = float(coverage) / float(mean_cov)
    if use_log:
        return float(np.log2(ratio))
    if not np.isfinite(scale_mean_cnv):
        return None
    return float(scale_mean_cnv) * ratio


def _panel_label_matches_configured(label: str, configured_genes: Sequence[str]) -> bool:
    """True when a panel target label matches a configured ``[cnv].genes`` symbol."""
    key = str(label).strip().casefold()
    if not key:
        return False
    for gene in configured_genes:
        want = str(gene).strip().casefold()
        if not want:
            continue
        if key == want:
            return True
        if key.startswith(f"{want}_") or key.startswith(f"{want}-") or key.startswith(
            f"{want} "
        ):
            return True
    return False


def _layout_panel_coverage_point_labels(
    coverage_points: List[Dict[str, Any]],
    y_min: float,
    y_max: float,
    x_max: float,
    *,
    x_key: str,
    min_x_spacing: float,
) -> Dict[tuple[str, float], float]:
    """Place gene labels above gains / below losses, staggering overlaps on the CNV axis."""
    layouts: Dict[tuple[str, float], float] = {}
    occupied: List[tuple[float, float]] = []
    x_spacing = max(x_max * 0.020, min_x_spacing)
    y_span = max(y_max - y_min, 1e-6)
    y_step = y_span * 0.055
    pad = y_span * 0.025

    for point in sorted(coverage_points, key=lambda item: item[x_key]):
        x_pos = float(point[x_key])
        y_head = float(point["y_norm"])
        above = point.get("direction") != "loss"
        label_y = y_head + pad if above else y_head - pad
        attempts = 0
        while any(
            abs(x_pos - ox) < x_spacing and abs(label_y - oy) < y_step
            for ox, oy in occupied
        ):
            label_y += y_step if above else -y_step
            attempts += 1
            if attempts > 8:
                break
        label_y = min(max(label_y, y_min + y_span * 0.02), y_max - y_span * 0.02)
        occupied.append((x_pos, label_y))
        layouts[(point["label"], x_pos)] = label_y
    return layouts


def _annotate_significant_panel_points(
    panel_points: List[Dict[str, Any]],
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
    *,
    use_log2: bool,
) -> List[Dict[str, Any]]:
    """Keep panel targets that are >3 SD outliers and have finite coverage."""
    significant_points: List[Dict[str, Any]] = []
    for point in panel_points:
        if not _should_label_panel_gene(point, mean_cnv, std_cnv):
            continue
        coverage_val = point.get("coverage_val")
        if coverage_val is None or not np.isfinite(coverage_val):
            continue
        significant_points.append(
            {
                **point,
                "direction": _gene_cnv_direction(point, regions, use_log2=use_log2),
            }
        )
    return significant_points


def _annotate_configured_panel_points(
    panel_points: List[Dict[str, Any]],
    regions: List[Dict[str, Any]],
    *,
    use_log2: bool,
    configured_genes: Sequence[str],
) -> List[Dict[str, Any]]:
    """Keep configured ``[cnv].genes`` targets that have finite coverage."""
    annotated: List[Dict[str, Any]] = []
    for point in panel_points:
        if not _panel_label_matches_configured(point["label"], configured_genes):
            continue
        coverage_val = point.get("coverage_val")
        if coverage_val is None or not np.isfinite(coverage_val):
            continue
        annotated.append(
            {
                **point,
                "direction": _gene_cnv_direction(point, regions, use_log2=use_log2),
            }
        )
    return annotated


def _select_panel_coverage_points(
    panel_points: List[Dict[str, Any]],
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
    *,
    use_log2: bool,
    configured_genes: Sequence[str] = (),
) -> List[Dict[str, Any]]:
    """Prefer configured gene list; otherwise keep >3 SD CNV outliers."""
    if configured_genes:
        return _annotate_configured_panel_points(
            panel_points,
            regions,
            use_log2=use_log2,
            configured_genes=configured_genes,
        )
    return _annotate_significant_panel_points(
        panel_points,
        mean_cnv,
        std_cnv,
        regions,
        use_log2=use_log2,
    )


def _attach_normalised_coverage(
    points: List[Dict[str, Any]],
    *,
    mean_cov: Optional[float],
    scale_mean_cnv: float,
    use_log2: bool,
) -> List[Dict[str, Any]]:
    """Add ``y_norm`` / ``baseline_y`` so coverage can share the CNV axis."""
    if not points:
        return []
    if mean_cov is None or not np.isfinite(mean_cov) or mean_cov <= 0:
        coverage_vals = np.asarray(
            [float(p["coverage_val"]) for p in points if p.get("coverage_val") is not None],
            dtype=float,
        )
        coverage_vals = coverage_vals[np.isfinite(coverage_vals) & (coverage_vals > 0)]
        if len(coverage_vals) == 0:
            return []
        mean_cov = float(np.mean(coverage_vals))

    baseline_y = 0.0 if use_log2 else float(scale_mean_cnv)
    normalised: List[Dict[str, Any]] = []
    for point in points:
        y_norm = _normalise_coverage_to_cnv_axis(
            float(point["coverage_val"]),
            mean_cov=float(mean_cov),
            scale_mean_cnv=float(scale_mean_cnv),
            use_log=use_log2,
        )
        if y_norm is None:
            continue
        normalised.append(
            {
                **point,
                "y_norm": float(y_norm),
                "baseline_y": float(baseline_y),
            }
        )
    return normalised


def _expand_ylim_for_coverage_points(
    y_min: float,
    y_max: float,
    coverage_points: Sequence[Dict[str, Any]],
) -> Tuple[float, float]:
    """Widen CNV axis limits so normalised coverage markers stay in view."""
    if not coverage_points:
        return y_min, y_max
    ys = [float(p["y_norm"]) for p in coverage_points if np.isfinite(p.get("y_norm", np.nan))]
    if not ys:
        return y_min, y_max
    pad = max((y_max - y_min) * 0.08, 0.15)
    return min(y_min, min(ys) - pad), max(y_max, max(ys) + pad)


def _add_panel_coverage_points(
    ax_cnv,
    coverage_points: List[Dict[str, Any]],
    x_max: float,
    *,
    x_key: str,
    min_x_spacing: float,
    y_min: float,
    y_max: float,
    label_font_size: Optional[int] = None,
) -> bool:
    """Plot panel genes as coverage normalised onto the shared CNV axis."""
    if not coverage_points:
        return False

    font_size = (
        int(label_font_size)
        if label_font_size is not None
        else LOLLIPOP_LABEL_FONT_SIZE
    )

    for point in coverage_points:
        x_pos = float(point[x_key])
        y_head = float(point["y_norm"])
        y_base = float(point.get("baseline_y", 0.0))
        color = (
            CNV_COLORS["plot_gain"]
            if point.get("direction") == "gain"
            else CNV_COLORS["plot_loss"]
        )
        ax_cnv.vlines(
            x_pos,
            min(y_base, y_head),
            max(y_base, y_head),
            colors=color,
            linewidths=0.8,
            alpha=1.0,
            zorder=5,
            clip_on=True,
        )

    for direction, color in (
        ("gain", CNV_COLORS["plot_gain"]),
        ("loss", CNV_COLORS["plot_loss"]),
    ):
        subset = [p for p in coverage_points if p.get("direction") == direction]
        if not subset:
            continue
        ax_cnv.scatter(
            [float(point[x_key]) for point in subset],
            [float(point["y_norm"]) for point in subset],
            s=18,
            color=color,
            zorder=6,
            edgecolors="white",
            linewidths=0.35,
            alpha=0.9,
            clip_on=True,
        )

    head_label_y = _layout_panel_coverage_point_labels(
        coverage_points,
        y_min,
        y_max,
        x_max,
        x_key=x_key,
        min_x_spacing=min_x_spacing,
    )
    for point in coverage_points:
        x_pos = float(point[x_key])
        color = (
            CNV_COLORS["plot_gain"]
            if point.get("direction") == "gain"
            else CNV_COLORS["plot_loss"]
        )
        above = point.get("direction") != "loss"
        ax_cnv.text(
            x_pos,
            head_label_y[(point["label"], x_pos)],
            _truncate_panel_label(point["label"]),
            ha="center",
            va="bottom" if above else "top",
            fontsize=font_size,
            color=color,
            fontweight="bold",
            zorder=8,
            clip_on=False,
            bbox=_CNV_GENE_LABEL_BBOX,
            path_effects=_CNV_GENE_LABEL_PATH_EFFECTS,
        )
    return True


def _gene_cnv_direction(
    point: Dict[str, Any],
    regions: List[Dict[str, Any]],
    *,
    use_log2: bool,
) -> str:
    """Classify a significant panel gene as gain-like or loss-like for point colour."""
    mid_bp = float(point["mid_mb"]) * 1_000_000
    for region in regions:
        if region.get("type") not in SIGNIFICANT_CNV_REGION_TYPES:
            continue
        start_bp = float(region["start_pos"])
        end_bp = float(region["end_pos"])
        if start_bp <= mid_bp <= end_bp:
            if region["type"] in ("GAIN", "HIGH_GAIN"):
                return "gain"
            if region["type"] in ("LOSS", "DEEP_LOSS"):
                return "loss"
    cnv_val = float(point["cnv_val"])
    if use_log2:
        return "gain" if cnv_val >= 0.0 else "loss"
    return "gain" if cnv_val >= 0.0 else "loss"


def _mean_target_coverage(
    target_coverage_df: Optional[pd.DataFrame],
    contig: Optional[str] = None,
) -> Optional[float]:
    """Mean sequencing coverage across panel targets (optionally one chromosome)."""
    if target_coverage_df is None or target_coverage_df.empty:
        return None
    if "coverage" not in target_coverage_df.columns:
        return None
    subset = target_coverage_df
    if contig is not None and "chrom" in target_coverage_df.columns:
        chrom_subset = target_coverage_df[target_coverage_df["chrom"] == contig]
        if not chrom_subset.empty:
            subset = chrom_subset
    vals = np.asarray(subset["coverage"], dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) == 0:
        return None
    return float(np.mean(vals))


def _add_genome_panel_coverage_points(
    ax_cnv,
    panel_points: List[Dict[str, Any]],
    x_max_bp: float,
    *,
    y_min: float,
    y_max: float,
    label_font_size: Optional[int] = None,
) -> bool:
    """Plot panel gene coverage markers on the genome-wide summary plot."""
    return _add_panel_coverage_points(
        ax_cnv,
        panel_points,
        x_max_bp,
        x_key="position_bp",
        min_x_spacing=1_800_000.0,
        y_min=y_min,
        y_max=y_max,
        label_font_size=label_font_size,
    )


def _add_chromosome_panel_coverage_overlay(
    ax_cnv,
    panel_points: List[Dict[str, Any]],
    x_max_mb: float,
    *,
    y_min: float,
    y_max: float,
    label_font_size: Optional[int] = None,
) -> bool:
    """Add coverage-normalised panel markers on the shared chromosome CNV axis."""
    if not panel_points:
        return False
    return _add_panel_coverage_points(
        ax_cnv,
        panel_points,
        float(ax_cnv.get_xlim()[1]) if ax_cnv.get_xlim()[1] > 0 else x_max_mb,
        x_key="mid_mb",
        min_x_spacing=1.8,
        y_min=y_min,
        y_max=y_max,
        label_font_size=label_font_size,
    )


def _collect_chromosome_significant_panel_points(
    panel_genes_df: Optional[pd.DataFrame],
    contig: str,
    values_array: np.ndarray,
    analysis_bin_width: int,
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
    target_coverage_df: Optional[pd.DataFrame],
    *,
    use_log2: bool,
    configured_genes: Sequence[str] = (),
) -> List[Dict[str, Any]]:
    """Collect panel coverage markers for one chromosome (configured genes or outliers)."""
    panel_points = _collect_panel_gene_points(
        panel_genes_df,
        contig,
        values_array,
        analysis_bin_width,
        use_max_abs=use_log2,
        target_coverage_df=target_coverage_df,
    )
    selected = _select_panel_coverage_points(
        panel_points,
        mean_cnv,
        std_cnv,
        regions,
        use_log2=use_log2,
        configured_genes=configured_genes,
    )
    return _attach_normalised_coverage(
        selected,
        # Use whole-panel mean coverage so gene markers match the genome summary.
        mean_cov=_mean_target_coverage(target_coverage_df),
        scale_mean_cnv=mean_cnv,
        use_log2=use_log2,
    )


def _collect_genome_significant_panel_points(
    panel_genes_df: Optional[pd.DataFrame],
    cnv_source: Dict[str, np.ndarray],
    ordered_contigs: List[str],
    chrom_start_offsets: Dict[str, float],
    analysis_bin_width: int,
    significant_regions: Optional[Dict[str, List[Dict[str, Any]]]],
    target_coverage_df: Optional[pd.DataFrame],
    *,
    use_log2: bool,
    scale_mean_cnv: float,
    configured_genes: Sequence[str] = (),
) -> List[Dict[str, Any]]:
    """Collect panel coverage markers with genome-wide bp positions."""
    if panel_genes_df is None or panel_genes_df.empty:
        return []

    significant_regions = significant_regions or {}
    genome_points: List[Dict[str, Any]] = []
    mean_cov = _mean_target_coverage(target_coverage_df)

    for contig in ordered_contigs:
        if contig not in cnv_source:
            continue
        values_array = np.asarray(cnv_source[contig], dtype=float)
        finite_values = values_array[np.isfinite(values_array)]
        if len(finite_values) == 0:
            continue

        mean_cnv = float(np.mean(finite_values))
        std_cnv = float(np.std(finite_values))
        regions = significant_regions.get(contig, [])
        panel_points = _collect_panel_gene_points(
            panel_genes_df,
            contig,
            values_array,
            analysis_bin_width,
            use_max_abs=use_log2,
            target_coverage_df=target_coverage_df,
        )

        chrom_offset = float(chrom_start_offsets.get(contig, 0.0))
        for point in _select_panel_coverage_points(
            panel_points,
            mean_cnv,
            std_cnv,
            regions,
            use_log2=use_log2,
            configured_genes=configured_genes,
        ):
            genome_points.append(
                {
                    **point,
                    "position_bp": chrom_offset + float(point["mid_mb"]) * 1_000_000,
                }
            )

    return _attach_normalised_coverage(
        genome_points,
        mean_cov=mean_cov,
        scale_mean_cnv=scale_mean_cnv,
        use_log2=use_log2,
    )


def _add_cnv_log2_reference_line(ax, x_max_mb: float) -> None:
    """Add a reference line at zero log2 ratio (no copy-number change)."""
    ax.axhline(
        0.0,
        color=CNV_COLORS["diploid"],
        linestyle="--",
        linewidth=1.3,
        alpha=0.95,
        zorder=1,
    )
    ax.text(
        x_max_mb * 0.008,
        0.0,
        "0",
        fontsize=CNV_FONT["annotation"],
        color=CNV_COLORS["diploid"],
        va="center",
        ha="left",
        zorder=1,
    )


def _add_cnv_ploidy_reference_lines(ax, y_max: float, x_max_mb: float) -> None:
    """Add dashed reference lines at integer ploidy values."""
    for ploidy in (1, 2, 3):
        if ploidy > y_max:
            continue
        is_diploid = ploidy == 2
        ax.axhline(
            ploidy,
            color=CNV_COLORS["diploid"] if is_diploid else CNV_COLORS["reference"],
            linestyle="--",
            linewidth=1.3 if is_diploid else 0.7,
            alpha=0.95 if is_diploid else 0.55,
            zorder=1,
        )
        ax.text(
            x_max_mb * 0.008,
            ploidy,
            str(ploidy),
            fontsize=CNV_FONT["annotation"],
            color=CNV_COLORS["diploid"] if is_diploid else CNV_COLORS["reference"],
            va="center",
            ha="left",
            zorder=1,
        )


def _draw_region_bracket(
    ax,
    start_mb: float,
    end_mb: float,
    y_level: float,
    color: str,
    height_frac: float = 0.04,
) -> None:
    """Draw a horizontal genomic bracket."""
    y_top = y_level
    y_drop = y_level - height_frac
    ax.plot(
        [start_mb, start_mb, end_mb, end_mb],
        [y_drop, y_top, y_top, y_drop],
        color=color,
        linewidth=1.1,
        solid_capstyle="butt",
        zorder=5,
        clip_on=False,
    )


def _apply_cnv_axes_style(ax, *, xlabel: str, ylabel: str, title: Optional[str] = None) -> None:
    """Apply consistent seaborn-inspired styling to a CNV axes."""
    _setup_cnv_fonts()
    ax.set_facecolor("white")
    ax.set_xlabel(
        xlabel,
        fontsize=CNV_FONT["axis"],
        color=CNV_TEXT["primary"],
        labelpad=10,
        fontproperties=_CNV_FONT_REGULAR,
    )
    ax.set_ylabel(
        ylabel,
        fontsize=CNV_FONT["axis"],
        color=CNV_TEXT["primary"],
        labelpad=10,
        fontproperties=_CNV_FONT_REGULAR,
    )
    if title:
        ax.set_title(
            title,
            fontsize=CNV_FONT["title"],
            color=CNV_TEXT["primary"],
            pad=10,
            fontproperties=_CNV_FONT_BOLD,
        )
    sns.despine(ax=ax, top=True, right=True)
    ax.grid(True, axis="y", color=MODERN_COLORS["grid"], linestyle="--", linewidth=0.4, alpha=0.55)
    ax.grid(False, axis="x")
    ax.tick_params(colors=CNV_TEXT["primary"], labelsize=CNV_FONT["tick"])


def _apply_cnv_genome_overview_axes(
    ax,
    *,
    xlabel: str,
    ylabel: str,
    title: Optional[str] = None,
    x_max_bp: float,
) -> None:
    """Genome-wide CNV panel: y-axis at x=0, no bottom axis line, no genomic tick labels."""
    _setup_cnv_fonts()
    ax.set_facecolor("white")
    ax.set_xlim(0, x_max_bp)
    ax.margins(x=0)
    ax.set_xlabel(
        xlabel,
        fontsize=CNV_FONT["axis"],
        color=CNV_TEXT["primary"],
        labelpad=18,
        fontproperties=_CNV_FONT_REGULAR,
    )
    ax.set_ylabel(
        ylabel,
        fontsize=CNV_FONT["axis"],
        color=CNV_TEXT["primary"],
        labelpad=10,
        fontproperties=_CNV_FONT_REGULAR,
    )
    if title:
        ax.set_title(
            title,
            fontsize=CNV_FONT["title"],
            color=CNV_TEXT["primary"],
            pad=10,
            fontproperties=_CNV_FONT_BOLD,
        )
    ax.spines["left"].set_position(("data", 0))
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("none")
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax.yaxis.set_ticks_position("left")
    ax.grid(True, axis="y", color=MODERN_COLORS["grid"], linestyle="--", linewidth=0.4, alpha=0.55)
    ax.grid(False, axis="x")
    ax.tick_params(axis="y", colors=CNV_TEXT["primary"], labelsize=CNV_FONT["tick"])


def _chromosome_cnv_dataframe(positions_mb, values) -> pd.DataFrame:
    return pd.DataFrame(
        {"position_mb": positions_mb, "ploidy": pd.Series(values, dtype=float)}
    )


def _add_cnv_reference_lines(ax, mean_cnv: float, std_cnv: float, y_min: float, y_max: float) -> None:
    """Genome-wide reference guides (mean plus optional spread)."""
    ax.axhline(y=mean_cnv, color=CNV_COLORS["reference"], linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)
    for offset in (std_cnv, 2 * std_cnv):
        if y_min < mean_cnv + offset <= y_max:
            ax.axhline(y=mean_cnv + offset, color=CNV_COLORS["reference"], linestyle=":", linewidth=0.6, alpha=0.4, zorder=1)
        if y_min <= mean_cnv - offset < y_max:
            ax.axhline(y=mean_cnv - offset, color=CNV_COLORS["reference"], linestyle=":", linewidth=0.6, alpha=0.4, zorder=1)


def _log2_linear_axis_limits(
    log2_values: np.ndarray,
    *,
    percentile: float = 97.5,
    min_span: float = 2.0,
    max_span: float = 4.0,
    pad: float = 0.08,
) -> Tuple[float, float]:
    """Data-driven symmetric linear limits for log2 ratio plots."""
    vals = np.asarray(log2_values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return -min_span, min_span
    span = float(np.percentile(np.abs(vals), percentile)) + pad
    span = max(span, min_span)
    span = min(span, max_span)
    return -span, span


def _cnv_plot_point_state(
    value: float,
    chromosome: str,
    sex_estimate: str,
) -> str:
    """Classify a log2 CNV value as gain, loss, or neutral using calling thresholds."""
    if not np.isfinite(value):
        return "neutral"
    from robin.classification_config import get_cnv_thresholds

    gain_thr, loss_thr = get_cnv_thresholds(chromosome, sex_estimate)
    if value > gain_thr:
        return "gain"
    if value < loss_thr:
        return "loss"
    return "neutral"


def _scatter_cnv_chromosome_points(
    ax,
    cnv_df: pd.DataFrame,
    *,
    color_by_state: bool,
    point_size: float = CNV_CHROMOSOME_SCATTER_SIZE,
) -> None:
    """Scatter per-chromosome CNV points, optionally coloured by threshold state."""
    if color_by_state and "state" in cnv_df.columns:
        for state, color, zorder in (
            ("neutral", CNV_COLORS["plot_neutral"], 1),
            ("loss", CNV_COLORS["plot_loss"], 2),
            ("gain", CNV_COLORS["plot_gain"], 2),
        ):
            subset = cnv_df[cnv_df["state"] == state]
            if subset.empty:
                continue
            ax.scatter(
                subset["position_mb"],
                subset["ploidy"],
                c=color,
                s=point_size,
                marker=CNV_POINT_MARKER,
                alpha=(
                    CNV_POINT_ALPHA_NEUTRAL
                    if state == "neutral"
                    else CNV_POINT_ALPHA_CALLED
                ),
                linewidth=0,
                edgecolors="none",
                rasterized=True,
                zorder=zorder,
            )
        return

    _plot_cnv_track(
        ax,
        cnv_df,
        float(cnv_df["position_mb"].max()),
        point_size=point_size,
    )


def _scatter_cnv_genome_points(ax, df: pd.DataFrame, *, color_by_state: bool) -> None:
    """Scatter genome-wide CNV points, optionally coloured by threshold state."""
    if color_by_state and "state" in df.columns:
        for state, color, zorder in (
            ("neutral", CNV_COLORS["plot_neutral"], 1),
            ("loss", CNV_COLORS["plot_loss"], 2),
            ("gain", CNV_COLORS["plot_gain"], 2),
        ):
            subset = df[df["state"] == state]
            if subset.empty:
                continue
            ax.scatter(
                subset["position_bp"],
                subset["ploidy"],
                c=color,
                s=4,
                marker=CNV_POINT_MARKER,
                alpha=(
                    CNV_POINT_ALPHA_NEUTRAL
                    if state == "neutral"
                    else CNV_POINT_ALPHA_CALLED
                ),
                linewidth=0,
                edgecolors="none",
                rasterized=True,
                zorder=zorder,
            )
        return

    palette = sns.color_palette("muted", n_colors=max(df["contig"].nunique(), 3))
    contig_palette = dict(zip(sorted(df["contig"].unique()), palette))
    for contig, color in contig_palette.items():
        subset = df[df["contig"] == contig]
        ax.scatter(
            subset["position_bp"],
            subset["ploidy"],
            c=[color],
            s=4,
            marker=CNV_POINT_MARKER,
            alpha=CNV_POINT_ALPHA_DEFAULT,
            linewidth=0,
            edgecolors="none",
            rasterized=True,
        )


def target_distribution_plot(df):
    """
    Creates a target distribution plot.

    Args:
        df (pd.DataFrame): DataFrame containing the target distribution data.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    set_modern_style()

    df["chrom"] = pd.Categorical(
        df["chrom"], categories=natsort.natsorted(df["chrom"].unique()), ordered=True
    )
    df = df.sort_values("chrom")

    # Generate the plot
    plt.figure(figsize=(16, 8))
    boxplot = sns.boxplot(
        x="chrom",
        y="coverage",
        data=df,
        color=MODERN_COLORS["accent"],
        flierprops={"marker": "o", "markerfacecolor": MODERN_COLORS["primary"]},
    )

    plt.title(
        "Distribution of Target Coverage on Each Chromosome",
        fontsize=12,
        color=MODERN_COLORS["primary"],
        pad=20,
    )
    plt.xlabel("Chromosome", color=MODERN_COLORS["primary"])
    plt.ylabel("Coverage", color=MODERN_COLORS["primary"])
    plt.xticks(rotation=45)

    # Identify and annotate outliers
    def annotate_outliers(df, boxplot):
        # Calculate quartiles and IQR
        for chrom in df["chrom"].unique():
            chrom_data = df[df["chrom"] == chrom]
            Q1 = chrom_data["coverage"].quantile(0.25)
            Q3 = chrom_data["coverage"].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR

            # Find outliers
            outliers = chrom_data[
                (chrom_data["coverage"] < lower_bound)
                | (chrom_data["coverage"] > upper_bound)
            ]

            for idx in outliers.index:
                outlier = outliers.loc[idx]
                boxplot.annotate(
                    outlier["name"],
                    xy=(df.loc[idx, "chrom"], df.loc[idx, "coverage"]),
                    xytext=(10, 10),  # Offset the text more from the point
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=12,
                    color="red",
                )  # Increase font size

    annotate_outliers(df, boxplot)

    plt.tight_layout()
    # Save the plot as a JPG file
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300)
    buf.seek(0)
    return buf


def _create_empty_cnv_buffer():
    """Create a minimal valid JPEG buffer for empty CNV plots."""
    try:
        plt.figure(figsize=(16, 4))
        plt.text(0.5, 0.5, "No CNV data available", 
                ha='center', va='center', transform=plt.gca().transAxes,
                fontsize=14, color='gray')
        plt.title("Copy Number Changes")
        plt.axis('off')
        
        buf = io.BytesIO()
        fig = plt.gcf()
        plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        
        # Validate buffer contains data
        if buf.getvalue():
            return buf
        else:
            raise ValueError("Empty buffer")
    except Exception:
        # If even this fails, create a minimal valid JPEG programmatically
        try:
            # Create a minimal 1x1 white JPEG
            from PIL import Image as PILImage
            img = PILImage.new('RGB', (1, 1), color='white')
            buf = io.BytesIO()
            img.save(buf, format='JPEG')
            buf.seek(0)
            return buf
        except Exception:
            # Last resort: return a minimal valid JPEG binary directly
            # This is a valid 1x1 white JPEG
            jpeg_bytes = bytes([
                0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10, 0x4A, 0x46, 0x49, 0x46, 0x00, 0x01,
                0x01, 0x01, 0x00, 0x48, 0x00, 0x48, 0x00, 0x00, 0xFF, 0xDB, 0x00, 0x43,
                0x00, 0x08, 0x06, 0x06, 0x07, 0x06, 0x05, 0x08, 0x07, 0x07, 0x07, 0x09,
                0x09, 0x08, 0x0A, 0x0C, 0x14, 0x0D, 0x0C, 0x0B, 0x0B, 0x0C, 0x19, 0x12,
                0x13, 0x0F, 0x14, 0x1D, 0x1A, 0x1F, 0x1E, 0x1D, 0x1A, 0x1C, 0x1C, 0x20,
                0x24, 0x2E, 0x27, 0x20, 0x22, 0x2C, 0x23, 0x1C, 0x1C, 0x28, 0x37, 0x29,
                0x2C, 0x30, 0x31, 0x34, 0x34, 0x34, 0x1F, 0x27, 0x39, 0x3D, 0x38, 0x32,
                0x3C, 0x2E, 0x33, 0x34, 0x32, 0xFF, 0xC0, 0x00, 0x0B, 0x08, 0x00, 0x01,
                0x00, 0x01, 0x01, 0x01, 0x11, 0x00, 0xFF, 0xC4, 0x00, 0x14, 0x00, 0x01,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x08, 0xFF, 0xC4, 0x00, 0x14, 0x10, 0x01, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0xFF, 0xDA, 0x00, 0x08, 0x01, 0x01, 0x00, 0x00, 0x3F, 0x00,
                0xD2, 0xCF, 0x20, 0xFF, 0xD9
            ])
            buf = io.BytesIO(jpeg_bytes)
            buf.seek(0)
            return buf


def create_CNV_plot(
    result,
    cnv_dict,
    normalized_cnv=None,
    *,
    use_normalized_difference: bool = False,
    plot_bin_width: Optional[int] = None,
    sex_estimate: str = "Unknown",
    panel_genes_df: Optional[pd.DataFrame] = None,
    target_coverage_df: Optional[pd.DataFrame] = None,
    significant_regions: Optional[Dict[str, List[Dict[str, Any]]]] = None,
    reference_contig_scope: Optional[str] = None,
    configured_genes: Sequence[str] = (),
):
    """
    Creates a CNV plot.

    Args:
        result (Result): CNV result object.
        cnv_dict (dict): Dictionary containing CNV data.
        normalized_cnv (dict, optional): Per-chromosome log2(ploidy / expected copy number)
            values derived from the absolute CNV track.
        use_normalized_difference (bool): When True, plot log2(ploidy / expected) instead of
            absolute ploidy for the genome-wide summary chart.
        plot_bin_width (int, optional): Display bin width in bp for genome-wide plot.
            Defaults to 1 Mb. Values below the analysis bin width are ignored.
        sex_estimate: Sex estimate for per-chromosome log2 calling thresholds.
        configured_genes: Optional ``[cnv].genes`` list. When non-empty, those panel
            targets are marked; otherwise >3 SD CNV outliers are used.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    try:
        from robin.analysis.cnv_analysis import (
            CNV_REPORT_GENOME_PLOT_BIN_WIDTH,
            downsample_cnv_for_plot,
            resolve_cnv_plot_bin_width,
        )
        from robin.gui.plotting_preferences import (
            CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE,
            CNV_REPORT_SCALE_PLOIDY,
            cnv_report_genome_ylabel_mathtext,
            resolve_cnv_gene_label_font_size,
        )

        set_modern_style()

        gene_label_font_size = resolve_cnv_gene_label_font_size()
        cnv_source = result.cnv if hasattr(result, "cnv") else None
        plot_normalized = use_normalized_difference and normalized_cnv
        if use_normalized_difference and not normalized_cnv:
            logger.warning(
                "Log2 CNV summary requested but no relative data available; "
                "falling back to absolute ploidy"
            )
        if plot_normalized:
            cnv_source = normalized_cnv
        scale = (
            CNV_REPORT_SCALE_NORMALIZED_DIFFERENCE
            if plot_normalized
            else CNV_REPORT_SCALE_PLOIDY
        )

        # Check if result has CNV data
        if not cnv_source:
            logger.warning("No CNV data available for plotting")
            return _create_empty_cnv_buffer()

        # Prepare data for plotting
        analysis_bin_width = int(cnv_dict["bin_width"])
        display_bin_width = resolve_cnv_plot_bin_width(
            analysis_bin_width,
            plot_bin_width or CNV_REPORT_GENOME_PLOT_BIN_WIDTH,
        )

        plot_rows = []
        offset_bp = 0.0
        contig_centers = {}
        contig_boundaries = []
        chrom_start_offsets: Dict[str, float] = {}
        log2_values_for_limits: List[float] = []
        ordered_contigs = [
            contig
            for contig in natsort.natsorted(cnv_source.keys())
            if is_visible_contig(contig, reference_contig_scope)
        ]

        for contig in ordered_contigs:
            values = np.asarray(cnv_source[contig], dtype=float)
            chrom_start_offsets[contig] = offset_bp
            chrom_span_bp = len(values) * analysis_bin_width
            x_local, plot_values = downsample_cnv_for_plot(
                values, analysis_bin_width, display_bin_width
            )
            x_global = offset_bp + x_local
            for position_bp, y_value in zip(x_global, plot_values):
                y_value = float(y_value)
                if plot_normalized and not np.isfinite(y_value):
                    continue
                if plot_normalized:
                    log2_values_for_limits.append(y_value)
                row = {
                    "contig": contig,
                    "position_bp": float(position_bp),
                    "ploidy": y_value,
                }
                if plot_normalized:
                    row["state"] = _cnv_plot_point_state(y_value, contig, sex_estimate)
                plot_rows.append(row)
            contig_centers[contig] = offset_bp + (chrom_span_bp / 2)
            offset_bp += chrom_span_bp
            contig_boundaries.append(offset_bp)

        if not plot_rows:
            logger.warning("No plot data available for CNV plot")
            return _create_empty_cnv_buffer()

        df = pd.DataFrame(plot_rows)
        mean_value = float(df["ploidy"].mean())
        std_value = float(df["ploidy"].std())
        if plot_normalized:
            y_min, y_max = _log2_linear_axis_limits(np.asarray(log2_values_for_limits))
        else:
            y_min = max(0.0, float(df["ploidy"].min()) - 0.25)
            y_max = max(mean_value + (4 * std_value), mean_value * 1.35, 2.5)

        width = CNV_GENOME_LANDSCAPE_FIG_WIDTH
        fig, ax = plt.subplots(
            figsize=(width, CNV_GENOME_LANDSCAPE_FIG_HEIGHT)
        )
        genome_panel_points = _collect_genome_significant_panel_points(
            panel_genes_df,
            cnv_source,
            ordered_contigs,
            chrom_start_offsets,
            analysis_bin_width,
            significant_regions,
            target_coverage_df,
            use_log2=plot_normalized,
            scale_mean_cnv=mean_value,
            configured_genes=configured_genes,
        )
        has_lollipops = bool(genome_panel_points)
        if has_lollipops:
            y_min, y_max = _expand_ylim_for_coverage_points(
                y_min, y_max, genome_panel_points
            )

        _scatter_cnv_genome_points(ax, df, color_by_state=plot_normalized)

        for boundary in contig_boundaries[:-1]:
            ax.axvline(
                boundary,
                color=MODERN_COLORS["grid"],
                linewidth=0.6,
                linestyle="-",
                alpha=0.8,
                zorder=0,
            )

        if plot_normalized:
            ax.axhline(
                y=0.0,
                xmin=0.0,
                xmax=1.0,
                color=CNV_TEXT["primary"],
                linestyle="-",
                linewidth=0.8,
                alpha=0.35,
                zorder=1,
                clip_on=True,
            )
        else:
            _add_cnv_reference_lines(ax, mean_value, std_value, y_min, y_max)

        label_y = y_min + (y_max - y_min) * 0.03
        for contig, center_bp in contig_centers.items():
            ax.text(
                center_bp,
                label_y,
                _chromosome_display_name(contig),
                fontsize=CNV_FONT["tick"],
                ha="center",
                va="bottom",
                rotation=0,
                color=CNV_TEXT["primary"],
                fontproperties=_CNV_FONT_REGULAR,
                clip_on=False,
            )

        ax.set_ylim(y_min, y_max)
        _apply_cnv_genome_overview_axes(
            ax,
            xlabel="Chromosome",
            ylabel=cnv_report_genome_ylabel_mathtext(scale),
            title="Copy number variation across chromosomes",
            x_max_bp=offset_bp,
        )
        if has_lollipops:
            _add_genome_panel_coverage_points(
                ax,
                genome_panel_points,
                offset_bp,
                y_min=y_min,
                y_max=y_max,
                label_font_size=gene_label_font_size,
            )

        buf = io.BytesIO()
        fig.savefig(buf, format="jpg", dpi=300, bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
        buf.seek(0)

        if not buf.getvalue():
            logger.warning("Empty buffer for CNV plot")
            return _create_empty_cnv_buffer()

        buf_data = buf.getvalue()
        if len(buf_data) < 2 or buf_data[:2] != b'\xff\xd8':
            logger.warning("Invalid JPEG data for CNV plot")
            return _create_empty_cnv_buffer()

        buf.seek(0)
        return buf
    except Exception as e:
        logger.error(f"Error creating CNV plot: {str(e)}")
        plt.close('all')
        return _create_empty_cnv_buffer()


REPORTABLE_CHROMOSOMES = ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]
PANEL_Y_EXPANSION_FACTOR = 6.0


def cnv_chromosome_fig_height_for_page(
    page_height_inch: float,
    *,
    plots_per_page: int = CNV_CHROMOSOME_PLOTS_PER_PAGE,
    spacer_pt: float = CNV_CHROMOSOME_PLOT_SPACER_PT,
    frame_padding_pt: float = CNV_REPORT_FRAME_PADDING_PT,
) -> float:
    """Matplotlib figure height so ``plots_per_page`` chromosome plots fit one PDF page.

    ``page_height_inch`` should be ``SimpleDocTemplate.height`` in inches. ReportLab
    frames reserve ``frame_padding_pt`` (default 12) of vertical space for padding.
    """
    usable_inch = page_height_inch - frame_padding_pt / 72.0
    spacer_inch = max(plots_per_page - 1, 0) * spacer_pt / 72.0
    return (usable_inch - spacer_inch) / plots_per_page


def cnv_chromosome_fig_width_for_page(
    page_width_inch: float,
    *,
    frame_padding_pt: float = CNV_REPORT_FRAME_PADDING_PT,
) -> float:
    """Matplotlib figure width that fits the printable ReportLab frame."""
    return page_width_inch - frame_padding_pt / 72.0


def cnv_genome_landscape_image_size_pt(
    *,
    left_margin_pt: float = 72.0,
    right_margin_pt: float = 72.0,
    top_margin_pt: float = 1.35 * 72.0,
    bottom_margin_pt: float = 72.0,
    frame_padding_pt: float = CNV_REPORT_FRAME_PADDING_PT,
    caption_reserve_pt: float = CNV_GENOME_LANDSCAPE_CAPTION_RESERVE_PT,
) -> tuple[float, float]:
    """Return (width, height) in points for the genome-wide CNV plot on A4 landscape."""
    from reportlab.lib.pagesizes import A4, landscape as rl_landscape

    page_w, page_h = rl_landscape(A4)
    frame_w = page_w - left_margin_pt - right_margin_pt - frame_padding_pt
    frame_h = (
        page_h
        - top_margin_pt
        - bottom_margin_pt
        - frame_padding_pt
        - caption_reserve_pt
    )
    # Prefer a wide landscape aspect (~2.4:1) while filling available width.
    width = max(frame_w, 1.0)
    target_height = width / 2.4
    height = min(max(target_height, 1.0), max(frame_h, 1.0))
    return width, height


def _panel_target_label(gene_row) -> str:
    """Return the display label for a target-panel region."""
    for key in ("gene", "name", "target"):
        if key in gene_row.index and pd.notna(gene_row[key]):
            raw = str(gene_row[key]).strip()
            if raw and raw.lower() != "nan":
                return raw.split(",")[0].strip()
    return ""


def _coverage_for_panel_target(
    gene_row: pd.Series,
    contig: str,
    target_coverage_df: Optional[pd.DataFrame],
) -> Optional[float]:
    """Resolve target coverage for a panel BED interval."""
    if target_coverage_df is None or target_coverage_df.empty:
        return None
    chrom_matches = target_coverage_df[target_coverage_df["chrom"] == contig]
    if chrom_matches.empty:
        return None

    start_bp = float(gene_row["start_pos"])
    end_bp = float(gene_row["end_pos"])
    overlapping = chrom_matches[
        (chrom_matches["endpos"] >= start_bp) & (chrom_matches["startpos"] <= end_bp)
    ]
    if not overlapping.empty:
        return float(overlapping["coverage"].max())

    label = _panel_target_label(gene_row)
    if not label:
        return None
    for _, row in chrom_matches.iterrows():
        names = [name.strip() for name in str(row["name"]).split(",") if name.strip()]
        if label in names:
            return float(row["coverage"])
    return None


def _collect_panel_gene_points(
    panel_genes_df: Optional[pd.DataFrame],
    contig: str,
    values,
    bin_width: int,
    *,
    use_max_abs: bool = False,
    target_coverage_df: Optional[pd.DataFrame] = None,
) -> List[Dict[str, Any]]:
    """Collect panel target positions and peak CNV across each target region."""
    if panel_genes_df is None or panel_genes_df.empty:
        return []

    values_array = np.array(values)
    genes = panel_genes_df[panel_genes_df["chrom"] == contig].sort_values("start_pos")
    points: List[Dict[str, Any]] = []

    for _, gene_row in genes.iterrows():
        label_text = _panel_target_label(gene_row)
        if not label_text:
            continue

        start_bp = float(gene_row["start_pos"])
        end_bp = float(gene_row["end_pos"])
        mid_bp = (start_bp + end_bp) / 2.0
        mid_mb = mid_bp / 1_000_000
        start_bin = max(0, int(start_bp // bin_width))
        end_bin = min(len(values_array) - 1, int(end_bp // bin_width))
        if end_bin < start_bin:
            cnv_val = 0.0
        else:
            region_vals = values_array[start_bin : end_bin + 1]
            if len(region_vals):
                if use_max_abs:
                    cnv_val = float(region_vals[np.nanargmax(np.abs(region_vals))])
                else:
                    cnv_val = float(np.nanmax(region_vals))
            else:
                cnv_val = 0.0

        coverage_val = _coverage_for_panel_target(gene_row, contig, target_coverage_df)

        points.append(
            {
                "label": label_text,
                "mid_mb": mid_mb,
                "cnv_val": cnv_val,
                "coverage_val": coverage_val,
            }
        )

    merged: Dict[str, Dict[str, Any]] = {}
    for point in points:
        existing = merged.get(point["label"])
        if existing is None:
            merged[point["label"]] = point
            continue
        if point["cnv_val"] > existing["cnv_val"]:
            existing["cnv_val"] = point["cnv_val"]
        if point.get("coverage_val") is not None:
            prev = existing.get("coverage_val")
            existing["coverage_val"] = max(prev or 0.0, point["coverage_val"])

    return list(merged.values())


def _compute_cnv_y_limits(
    values_array: np.ndarray,
    panel_points: List[Dict[str, Any]],
) -> tuple[float, float, bool, float, float]:
    """
    Choose y-axis limits that preserve bulk CNV detail while surfacing
    amplified panel targets when feasible.
    """
    mean_cnv = float(np.mean(values_array))
    std_cnv = float(np.std(values_array))
    baseline_max = max(mean_cnv + (2 * std_cnv), mean_cnv * 1.4, 2.5)

    if not panel_points:
        return 0.0, baseline_max, False, mean_cnv, std_cnv

    panel_max = max(point["cnv_val"] for point in panel_points)
    if panel_max <= baseline_max * 1.05:
        return 0.0, baseline_max, False, mean_cnv, std_cnv

    if panel_max <= baseline_max * PANEL_Y_EXPANSION_FACTOR:
        return 0.0, panel_max * 1.08, False, mean_cnv, std_cnv

    return 0.0, baseline_max, True, mean_cnv, std_cnv


def _compute_log2_y_limits(
    values_array: np.ndarray,
    panel_points: List[Dict[str, Any]],
) -> tuple[float, float, bool, float, float]:
    """Symmetric log2 ratio limits for per-chromosome difference plots."""
    mean_cnv = float(np.mean(values_array))
    std_cnv = float(np.std(values_array))
    y_min, y_max = _log2_linear_axis_limits(values_array)
    off_scale_mode = False
    if panel_points:
        panel_extreme = max(abs(point["cnv_val"]) for point in panel_points)
        if panel_extreme > y_max * 0.97:
            off_scale_mode = True
    return y_min, y_max, off_scale_mode, mean_cnv, std_cnv


def create_CNV_plot_per_chromosome(
    result,
    cnv_dict,
    significant_regions=None,
    chromosomes: Optional[List[str]] = None,
    panel_genes_df: Optional[pd.DataFrame] = None,
    chromosome_status: Optional[Dict[str, str]] = None,
    normalized_cnv: Optional[Dict[str, np.ndarray]] = None,
    target_coverage_df: Optional[pd.DataFrame] = None,
    *,
    use_log2_ratio: bool = False,
    plot_bin_width: Optional[int] = None,
    sex_estimate: str = "Unknown",
    fig_height: Optional[float] = None,
    fig_width: Optional[float] = None,
    reference_contig_scope: Optional[str] = None,
    configured_genes: Sequence[str] = (),
):
    """Creates CNV plots per chromosome.

    Args:
        result (Result): CNV result object.
        cnv_dict (dict): Dictionary containing CNV data.
        significant_regions (dict): Dictionary mapping chromosomes to lists of significant regions.
        chromosomes (list, optional): Ordered chromosome names to plot.
        panel_genes_df (pd.DataFrame, optional): Target panel genes.
        chromosome_status (dict, optional): Status text keyed by chromosome.
        normalized_cnv (dict, optional): Per-chromosome log2(ploidy / expected) values.
        sex_estimate (str): Sample sex estimate for log2 threshold colouring.
        use_log2_ratio (bool): When True, plot log2(ploidy / expected) instead of absolute ploidy.
        plot_bin_width (int, optional): Display bin width in bp. Defaults to the
            sample analysis bin width (finest resolution available).
        fig_height (float, optional): Matplotlib figure height in inches. Defaults
            to a height that fits four plots per PDF page.
        fig_width (float, optional): Matplotlib figure width in inches. Defaults
            to ``CNV_CHROMOSOME_FIG_WIDTH``.
        configured_genes: Optional ``[cnv].genes`` list. When non-empty, those panel
            targets are marked; otherwise >3 SD CNV outliers are used.

    Returns:
        List[Tuple[str, io.BytesIO]]: List of tuples containing chromosome names and plot buffers.
    """
    plots = []
    chromosome_status = chromosome_status or {}
    try:
        from robin.analysis.cnv_analysis import (
            downsample_cnv_chromosome_track,
            resolve_cnv_plot_bin_width,
        )
        from robin.gui.plotting_preferences import (
            cnv_report_genome_ylabel_mathtext,
            resolve_cnv_gene_label_font_size,
        )

        set_modern_style()

        plot_log2 = use_log2_ratio and normalized_cnv
        if use_log2_ratio and not normalized_cnv:
            logger.warning(
                "Log2 per-chromosome CNV plots requested but no relative data available; "
                "falling back to absolute ploidy"
            )

        cnv_source = normalized_cnv if plot_log2 else (result.cnv if hasattr(result, "cnv") else None)
        if not cnv_source:
            logger.warning("No CNV data available for per-chromosome plotting")
            return plots

        gene_label_font_size = resolve_cnv_gene_label_font_size()

        if chromosomes is None:
            chromosomes = [
                contig
                for contig in natsort.natsorted(cnv_source.keys())
                if is_visible_contig(contig, reference_contig_scope)
            ]

        scale = "normalized_difference" if plot_log2 else "ploidy"
        ylabel = cnv_report_genome_ylabel_mathtext(scale)
        analysis_bin_width = int(cnv_dict["bin_width"])
        report_plot_bin_width = resolve_cnv_plot_bin_width(
            analysis_bin_width,
            plot_bin_width,
        )
        compact_layout = True
        total_height = fig_height or cnv_chromosome_fig_height_for_page(9.34)
        total_width = fig_width or CNV_CHROMOSOME_FIG_WIDTH

        for contig in chromosomes:
            if contig not in cnv_source:
                continue
            values = cnv_source[contig]

            values_array = np.asarray(values, dtype=float)
            finite_values = values_array[np.isfinite(values_array)]
            if len(finite_values) == 0:
                continue

            positions_mb, plot_values, x_max_mb = downsample_cnv_chromosome_track(
                values_array,
                analysis_bin_width,
                report_plot_bin_width,
            )
            plot_mask = np.isfinite(plot_values)
            positions_mb = positions_mb[plot_mask]
            plot_values = plot_values[plot_mask]
            if len(plot_values) == 0:
                continue
            if plot_log2:
                y_min, y_max, _, mean_cnv, std_cnv = _compute_log2_y_limits(
                    finite_values, [],
                )
            else:
                y_min, y_max, _, mean_cnv, std_cnv = _compute_cnv_y_limits(
                    finite_values, [],
                )
                y_min = 0.0
            regions = (significant_regions or {}).get(contig, [])
            significant_panel_points = _collect_chromosome_significant_panel_points(
                panel_genes_df,
                contig,
                values_array,
                analysis_bin_width,
                mean_cnv,
                std_cnv,
                regions,
                target_coverage_df,
                use_log2=plot_log2,
                configured_genes=configured_genes,
            )
            if significant_panel_points:
                y_min, y_max = _expand_ylim_for_coverage_points(
                    y_min, y_max, significant_panel_points
                )
            chrom_name = _chromosome_display_name(contig)
            status_text = chromosome_status.get(contig, "No significant CNV change")

            has_status = bool(
                status_text and status_text != "No significant CNV change"
            )
            fig = plt.figure(figsize=(total_width, total_height))
            ax = fig.add_subplot(111)

            fig.subplots_adjust(
                **_chromosome_figure_margins(
                    has_status=has_status,
                    compact=compact_layout,
                )
            )
            fig.suptitle(
                f"Chromosome {chrom_name}",
                fontsize=CNV_FONT["title"],
                color=CNV_TEXT["primary"],
                y=0.99 if has_status else 0.97,
                fontproperties=_CNV_FONT_BOLD,
            )
            if has_status:
                wrapped_status = _wrap_chromosome_status_text(status_text)
                fig.text(
                    0.5,
                    0.88 if compact_layout else 0.90,
                    wrapped_status,
                    ha="center",
                    va="top",
                    fontsize=CNV_FONT["subtitle"],
                    color=CNV_TEXT["muted"],
                    linespacing=1.2,
                    fontproperties=_CNV_FONT_REGULAR,
                )

            cnv_df = _chromosome_cnv_dataframe(positions_mb, plot_values)
            if plot_log2:
                cnv_df["state"] = [
                    _cnv_plot_point_state(value, contig, sex_estimate)
                    for value in cnv_df["ploidy"]
                ]
            if plot_log2:
                _add_cnv_log2_reference_line(ax, x_max_mb)
            else:
                _add_cnv_ploidy_reference_lines(ax, y_max, x_max_mb)
            _scatter_cnv_chromosome_points(
                ax, cnv_df, color_by_state=plot_log2,
            )
            _apply_cnv_chromosome_axes(
                ax,
                x_max_mb,
                y_max,
                y_min=y_min,
                xlabel="Position (Mb)",
                ylabel=ylabel,
                show_xlabel=True,
            )
            if significant_panel_points:
                logger.debug(
                    "Adding coverage overlay with %d panel targets on %s",
                    len(significant_panel_points),
                    contig,
                )
                _add_chromosome_panel_coverage_overlay(
                    ax,
                    significant_panel_points,
                    x_max_mb,
                    y_min=y_min,
                    y_max=y_max,
                    label_font_size=gene_label_font_size,
                )

            try:
                pad_inches = 0.10 if compact_layout else 0.12
                if significant_panel_points:
                    pad_inches = 0.12 if compact_layout else 0.15
                buf = io.BytesIO()
                fig.savefig(
                    buf,
                    format="jpg",
                    dpi=300,
                    bbox_inches="tight",
                    pad_inches=pad_inches,
                )
                plt.close(fig)
                buf.seek(0)

                if not buf.getvalue():
                    logger.warning(f"Empty buffer for chromosome {contig} CNV plot")
                    continue

                buf_data = buf.getvalue()
                if len(buf_data) < 2 or buf_data[:2] != b'\xff\xd8':
                    logger.warning(f"Invalid JPEG data for chromosome {contig} CNV plot")
                    continue

                buf.seek(0)
                plots.append((contig, buf))
            except Exception as e:
                logger.error(f"Error creating CNV plot for chromosome {contig}: {str(e)}")
                plt.close('all')
                continue

    except Exception as e:
        logger.error(f"Error in create_CNV_plot_per_chromosome: {str(e)}")
        plt.close('all')

    return plots


def classification_plot(df, title, threshold):
    """
    Creates a classification plot.

    Args:
        df (pd.DataFrame): DataFrame containing the classification data.
        title (str): Title of the plot.
        threshold (float): Threshold value for filtering classifications.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    set_modern_style()

    df["timestamp"] = pd.to_datetime(df["timestamp"], unit="ms", utc=True)

    # Reshape the data to long format
    df_melted = df.melt(id_vars=["timestamp"], var_name="Condition", value_name="Value")
    meta_conditions = {
        "number_probes",
        "covered_cpgs",
        "temperature",
        "diagnostic",
        "probes",
    }
    df_melted = df_melted[
        ~df_melted["Condition"].astype(str).str.strip().str.lower().isin(meta_conditions)
    ]

    # Filter conditions that cross the threshold
    top_conditions = df_melted.groupby("Condition")["Value"].max().nlargest(10).index
    df_filtered = df_melted[df_melted["Condition"].isin(top_conditions)]

    conditions_above_threshold = df_filtered[df_filtered["Value"] > threshold][
        "Condition"
    ].unique()
    df_filtered = df_filtered[df_filtered["Condition"].isin(conditions_above_threshold)]

    # Create figure with adjusted size and margins
    fig = plt.figure(figsize=(10, 6))

    # Only create legend if we have data to plot
    if not df_filtered.empty:
        sns.lineplot(
            data=df_filtered, x="timestamp", y="Value", hue="Condition", palette="Set2"
        )

        # Move the legend below the plot with adjusted position
        plt.legend(
            title="Condition", bbox_to_anchor=(0.5, -0.3), loc="upper center", ncol=3
        )
    else:
        # If no data, create an empty plot
        plt.plot([])
        plt.text(
            0.5,
            0.5,
            "No classification data above threshold",
            horizontalalignment="center",
            verticalalignment="center",
        )

    plt.title(f"{title} Classifications over Time")
    plt.xlabel("Timestamp")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    # Format the x-axis with custom date format
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%Y-%m-%d %H:%M"))

    # Adjust layout with explicit margins
    plt.subplots_adjust(bottom=0.25, left=0.1, right=0.9, top=0.9)

    # Save the plot as a JPG file with reduced DPI
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    plt.close(fig)  # Close the figure to free memory
    buf.seek(0)
    return buf


def coverage_plot(df):
    """
    Creates a coverage plot.

    Args:
        df (pd.DataFrame): DataFrame containing coverage data.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    set_modern_style()

    # df = df[df["#rname"] != "chrM"].copy()
    df = df[
        df["#rname"].isin(["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"])
    ].copy()

    # Sort chromosomes naturally
    df["#rname"] = pd.Categorical(
        df["#rname"], categories=natsort.natsorted(df["#rname"].unique()), ordered=True
    )
    df = df.sort_values("#rname")

    # Create figure with adjusted size
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace=0.4)

    # Font size settings
    title_fontsize = 8
    label_fontsize = 8
    tick_fontsize = 6

    # Plot number of reads per chromosome
    ax0 = plt.subplot(gs[0])
    sns.barplot(
        x="#rname", y="numreads", data=df, ax=ax0, color=MODERN_COLORS["accent"]
    )
    ax0.set_title("Number of Reads per Chromosome", fontsize=title_fontsize)
    ax0.set_xlabel("", fontsize=label_fontsize)
    ax0.set_ylabel("Number of Reads", fontsize=label_fontsize)
    ax0.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax0.tick_params(axis="y", labelsize=tick_fontsize)

    # Plot number of bases per chromosome
    ax1 = plt.subplot(gs[1])
    sns.barplot(
        x="#rname", y="covbases", data=df, ax=ax1, color=MODERN_COLORS["accent"]
    )
    ax1.set_title("Number of Bases per Chromosome", fontsize=title_fontsize)
    ax1.set_xlabel("", fontsize=label_fontsize)
    ax1.set_ylabel("Number of Bases", fontsize=label_fontsize)
    ax1.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax1.tick_params(axis="y", labelsize=tick_fontsize)

    # Plot mean depth per chromosome
    ax2 = plt.subplot(gs[2])
    sns.barplot(
        x="#rname", y="meandepth", data=df, ax=ax2, color=MODERN_COLORS["accent"]
    )
    ax2.set_title("Mean Depth per Chromosome", fontsize=title_fontsize)
    ax2.set_xlabel("Chromosome", fontsize=label_fontsize)
    ax2.set_ylabel("Mean Depth", fontsize=label_fontsize)
    ax2.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax2.tick_params(axis="y", labelsize=tick_fontsize)

    # Adjust layout with explicit margins
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, hspace=0.5)

    # Save the plot as a JPG file with reduced DPI
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    plt.close(fig)  # Close the figure to free memory
    buf.seek(0)
    return buf


def plot_classification_timeline(
    df: pd.DataFrame, classification_level: str = "class", title: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot classification changes over time.

    Args:
        df: DataFrame containing classification data
        classification_level: Level of classification to plot ('superfamily', 'family', 'class', 'subclass')
        title: Optional title for the plot

    Returns:
        Tuple of (figure, axes) objects
    """
    if df.empty:
        logger.warning("No data available for plotting")
        return None, None

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot classification scores
    ax.plot(
        df["timestamp"],
        df[f"{classification_level}_score"],
        marker="o",
        linestyle="-",
        label="Score",
    )

    # Add labels for each point
    for x, y, label in zip(
        df["timestamp"],
        df[f"{classification_level}_score"],
        df[f"{classification_level}_label"],
    ):
        ax.annotate(label, (x, y), xytext=(5, 5), textcoords="offset points")

    # Customize plot
    ax.set_xlabel("Time")
    ax.set_ylabel("Classification Score")
    ax.set_title(title or f"{classification_level.title()} Classification Over Time")
    ax.grid(True, alpha=0.3)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Adjust layout
    plt.tight_layout()

    return fig, ax


def plot_mgmt_timeline(
    df: pd.DataFrame, title: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot MGMT methylation changes over time.

    Args:
        df: DataFrame containing MGMT data
        title: Optional title for the plot

    Returns:
        Tuple of (figure, axes) objects
    """
    if df.empty:
        logger.warning("No data available for plotting")
        return None, None

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot methylation percentage
    ax.plot(
        df["timestamp"],
        df["mgmt_methylation"],
        marker="o",
        linestyle="-",
        label="Methylation %",
    )

    # Add status labels
    for x, y, status in zip(df["timestamp"], df["mgmt_methylation"], df["mgmt_status"]):
        ax.annotate(status, (x, y), xytext=(5, 5), textcoords="offset points")

    # Customize plot
    ax.set_xlabel("Time")
    ax.set_ylabel("Methylation Percentage")
    ax.set_title(title or "MGMT Methylation Over Time")
    ax.grid(True, alpha=0.3)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Adjust layout
    plt.tight_layout()

    return fig, ax


def save_plot(fig: plt.Figure, output_path: str, dpi: int = 300):
    """Save a plot to a file.

    Args:
        fig: Figure object to save
        output_path: Path to save the plot
        dpi: DPI for the output image
    """
    try:
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Save plot
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    except Exception as e:
        logger.error(f"Error saving plot to {output_path}: {str(e)}")


def plot_to_bytes(fig: plt.Figure, format: str = "png", dpi: int = 300) -> bytes:
    """Convert a plot to bytes.

    Args:
        fig: Figure object to convert
        format: Output format ('png', 'jpg', etc.)
        dpi: DPI for the output image

    Returns:
        Bytes containing the plot image
    """
    try:
        buf = io.BytesIO()
        fig.savefig(buf, format=format, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()

    except Exception as e:
        logger.error(f"Error converting plot to bytes: {str(e)}")
        return None
