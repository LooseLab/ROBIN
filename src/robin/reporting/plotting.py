"""
plotting.py

This module contains functions for creating plots used in the PDF report.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import textwrap
import matplotlib.font_manager as fm
import os
from robin.gui import fonts

import natsort

from matplotlib import gridspec

import logging
from typing import List, Optional, Tuple, Dict, Any

logger = logging.getLogger(__name__)

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


def _plot_cnv_track(ax, cnv_df: pd.DataFrame, x_max_mb: float) -> None:
    """Render CNV bin values as scatter with a rolling median trace."""
    ax.scatter(
        cnv_df["position_mb"],
        cnv_df["ploidy"],
        s=4,
        c=CNV_COLORS["points"],
        alpha=0.35,
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


def _is_gene_cnv_outlier(cnv_val: float, mean_cnv: float, std_cnv: float) -> bool:
    """Return True when a panel target copy number deviates >2 SD from the chromosome mean."""
    if std_cnv < 1e-6:
        return abs(cnv_val - mean_cnv) > 0.5
    return abs(cnv_val - mean_cnv) > 2 * std_cnv


SIGNIFICANT_CNV_REGION_TYPES = {"GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"}


def _gene_overlaps_cnv_region(mid_bp: float, regions: List[Dict[str, Any]]) -> bool:
    """Return True when a gene midpoint falls inside a called gain/loss region."""
    for region in regions:
        if region.get("type") not in SIGNIFICANT_CNV_REGION_TYPES:
            continue
        start_bp = float(region["start_pos"])
        end_bp = float(region["end_pos"])
        if start_bp <= mid_bp <= end_bp:
            return True
    return False


def _should_label_panel_gene(
    point: Dict[str, Any],
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
) -> bool:
    """Label outliers and panel genes inside called CNV regions."""
    if _is_gene_cnv_outlier(point["cnv_val"], mean_cnv, std_cnv):
        return True
    return _gene_overlaps_cnv_region(point["mid_mb"] * 1_000_000, regions)


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
    highlight = off_scale or _should_label_panel_gene(point, mean_cnv, std_cnv, regions)
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
    has_coverage_axis: bool,
) -> Dict[str, float]:
    return {
        "top": 0.76 if has_status else 0.86,
        "bottom": 0.18,
        "right": 0.88 if has_coverage_axis else 0.97,
        "left": 0.14 if has_coverage_axis else 0.13,
    }


def _layout_coverage_head_labels(
    coverage_points: List[Dict[str, Any]],
    cov_ylim: float,
    x_max_mb: float,
) -> Dict[tuple[str, float], float]:
    """Place lollipop labels above each marker, staggering overlaps in coverage space."""
    layouts: Dict[tuple[str, float], float] = {}
    occupied: List[tuple[float, float]] = []
    x_spacing = max(x_max_mb * 0.020, 1.8)
    y_step = cov_ylim * 0.055

    for point in sorted(coverage_points, key=lambda item: item["mid_mb"]):
        mid_mb = float(point["mid_mb"])
        y_top = min(float(point["coverage_val"]), cov_ylim * 0.90)
        label_y = y_top + cov_ylim * 0.025
        attempts = 0
        while any(
            abs(mid_mb - ox) < x_spacing and abs(label_y - oy) < y_step
            for ox, oy in occupied
        ):
            label_y += y_step
            attempts += 1
            if attempts > 8:
                break
        label_y = min(label_y, cov_ylim * 0.97)
        occupied.append((mid_mb, label_y))
        layouts[(point["label"], mid_mb)] = label_y
    return layouts


def _add_panel_gene_lollipops(
    fig: plt.Figure,
    ax_cnv,
    panel_points: List[Dict[str, Any]],
    x_max_mb: float,
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
) -> None:
    """Draw panel lollipops on a right-hand coverage axis (0 to max target coverage)."""
    coverage_points = [
        point
        for point in panel_points
        if point.get("coverage_val") is not None and np.isfinite(point["coverage_val"])
    ]
    if not coverage_points:
        return

    ax_cov = ax_cnv.twinx()
    cov_max = max(float(point["coverage_val"]) for point in coverage_points)
    cov_ylim = max(cov_max * 1.22, cov_max + 1.0)
    ax_cov.set_ylim(0.0, cov_ylim)
    ax_cov.set_ylabel(
        "Coverage (x)",
        fontsize=CNV_FONT["axis"],
        color=CNV_TEXT["primary"],
        labelpad=8,
        fontproperties=_CNV_FONT_REGULAR,
    )
    ax_cov.tick_params(
        colors=CNV_TEXT["primary"],
        labelsize=CNV_FONT["tick"],
        pad=2,
    )
    ax_cov.spines["right"].set_color(CNV_TEXT["primary"])
    ax_cov.spines["top"].set_visible(False)
    ax_cov.grid(False)

    head_label_y = _layout_coverage_head_labels(coverage_points, cov_ylim, x_max_mb)

    for point in coverage_points:
        mid_mb = point["mid_mb"]
        coverage_val = float(point["coverage_val"])
        highlight = _should_label_panel_gene(point, mean_cnv, std_cnv, regions)
        color = AMPLIFIED_COLOR if highlight else CNV_COLORS["gene"]
        stem_alpha = 0.9 if highlight else 0.45
        stem_width = 1.0 if highlight else 0.6
        y_top = min(coverage_val, cov_ylim * 0.90)

        ax_cov.plot(
            [mid_mb, mid_mb],
            [0.0, y_top],
            color=color,
            linewidth=stem_width,
            alpha=stem_alpha,
            zorder=5,
            solid_capstyle="round",
            clip_on=True,
        )
        ax_cov.scatter(
            [mid_mb],
            [y_top],
            s=24 if highlight else 12,
            color=color,
            zorder=6,
            edgecolors="white",
            linewidths=0.35,
            alpha=0.95 if highlight else 0.7,
        )

        short_label = _truncate_panel_label(point["label"])
        label_y = head_label_y[(point["label"], mid_mb)]
        ax_cov.text(
            mid_mb,
            label_y,
            short_label,
            ha="center",
            va="bottom",
            fontsize=LOLLIPOP_LABEL_FONT_SIZE,
            color=color,
            fontweight="bold" if highlight else "normal",
            zorder=8,
            clip_on=False,
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
            Defaults to 500 kb. Values below the analysis bin width are ignored.

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
        )

        set_modern_style()

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
        log2_values_for_limits: List[float] = []
        reportable = ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]
        ordered_contigs = [
            contig for contig in natsort.natsorted(cnv_source.keys()) if contig in reportable
        ]

        for contig in ordered_contigs:
            values = np.asarray(cnv_source[contig], dtype=float)
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
                plot_rows.append(
                    {
                        "contig": contig,
                        "position_bp": float(position_bp),
                        "ploidy": y_value,
                    }
                )
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

        palette = sns.color_palette("muted", n_colors=max(len(ordered_contigs), 3))
        contig_palette = dict(zip(ordered_contigs, palette))

        width = 16
        fig, ax = plt.subplots(figsize=(width, width / 4))
        sns.scatterplot(
            data=df,
            x="position_bp",
            y="ploidy",
            hue="contig",
            palette=contig_palette,
            ax=ax,
            legend=False,
            s=4,
            alpha=0.28,
            linewidth=0,
            edgecolor=None,
            rasterized=True,
        )

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
                color=CNV_COLORS["reference"],
                linestyle="--",
                linewidth=0.9,
                alpha=0.7,
                zorder=1,
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
            )

        ax.set_ylim(y_min, y_max)
        _apply_cnv_axes_style(
            ax,
            xlabel="Genomic position (bp)",
            ylabel=cnv_report_genome_ylabel_mathtext(scale),
            title="Copy number variation across chromosomes",
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
CNV_CHROMOSOME_FIG_WIDTH = 7.5
PANEL_Y_EXPANSION_FACTOR = 6.0
AMPLIFIED_COLOR = "#C62828"


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
        use_log2_ratio (bool): When True, plot log2(ploidy / expected) instead of absolute ploidy.
        plot_bin_width (int, optional): Display bin width in bp. Defaults to 500 kb
            for report plots (GUI default is the analysis bin width).

    Returns:
        List[Tuple[str, io.BytesIO]]: List of tuples containing chromosome names and plot buffers.
    """
    plots = []
    chromosome_status = chromosome_status or {}
    try:
        from robin.analysis.cnv_analysis import (
            CNV_REPORT_GENOME_PLOT_BIN_WIDTH,
            downsample_cnv_chromosome_track,
        )
        from robin.gui.plotting_preferences import cnv_report_genome_ylabel_mathtext

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

        if chromosomes is None:
            chromosomes = [
                contig
                for contig in natsort.natsorted(cnv_source.keys())
                if contig in REPORTABLE_CHROMOSOMES
            ]

        scale = "normalized_difference" if plot_log2 else "ploidy"
        ylabel = cnv_report_genome_ylabel_mathtext(scale)
        analysis_bin_width = int(cnv_dict["bin_width"])
        report_plot_bin_width = (
            plot_bin_width
            if plot_bin_width is not None
            else CNV_REPORT_GENOME_PLOT_BIN_WIDTH
        )

        for contig in chromosomes:
            if contig not in cnv_source:
                continue
            values = cnv_source[contig]
            if contig not in REPORTABLE_CHROMOSOMES:
                continue

            values_array = np.asarray(values, dtype=float)
            finite_values = values_array[np.isfinite(values_array)]
            if len(finite_values) == 0:
                continue

            panel_points = _collect_panel_gene_points(
                panel_genes_df,
                contig,
                values_array,
                analysis_bin_width,
                use_max_abs=plot_log2,
                target_coverage_df=target_coverage_df,
            )

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
            chrom_name = _chromosome_display_name(contig)
            status_text = chromosome_status.get(contig, "No significant CNV change")

            has_coverage_lollipops = any(
                point.get("coverage_val") is not None and np.isfinite(point["coverage_val"])
                for point in panel_points
            )
            has_status = bool(
                status_text and status_text != "No significant CNV change"
            )
            total_height = CNV_CHROMOSOME_FIG_WIDTH * 0.30
            fig = plt.figure(figsize=(CNV_CHROMOSOME_FIG_WIDTH, total_height))
            ax = fig.add_subplot(111)

            fig.subplots_adjust(
                **_chromosome_figure_margins(
                    has_status=has_status,
                    has_coverage_axis=has_coverage_lollipops,
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
                    0.90,
                    wrapped_status,
                    ha="center",
                    va="top",
                    fontsize=CNV_FONT["subtitle"],
                    color=CNV_TEXT["muted"],
                    linespacing=1.2,
                    fontproperties=_CNV_FONT_REGULAR,
                )

            cnv_df = _chromosome_cnv_dataframe(positions_mb, plot_values)
            if regions:
                _add_cnv_regions_on_plot(ax, regions, y_max)
            if plot_log2:
                _add_cnv_log2_reference_line(ax, x_max_mb)
            else:
                _add_cnv_ploidy_reference_lines(ax, y_max, x_max_mb)
            _plot_cnv_track(ax, cnv_df, x_max_mb)
            _apply_cnv_chromosome_axes(
                ax,
                x_max_mb,
                y_max,
                y_min=y_min,
                xlabel="Position (Mb)",
                ylabel=ylabel,
                show_xlabel=True,
            )
            _add_panel_gene_lollipops(
                fig, ax, panel_points, x_max_mb, mean_cnv, std_cnv, regions,
            )

            try:
                buf = io.BytesIO()
                fig.savefig(buf, format="jpg", dpi=300, bbox_inches="tight", pad_inches=0.12)
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
    df_melted = df_melted[df_melted["Condition"].ne("number_probes")]

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
